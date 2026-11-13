"""Replay saved readouts and audit results; never refit models or establish external validity."""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from regjepa.evaluation import paired_bootstrap, ranking_records, summarize_records
from regjepa.experiment import (average_repetitions, candidate_subset, component_score_diagnostics,
                               grid_edge_diagnostics, penalty_candidates, penalty_vector)
from regjepa.feature_io import load_feature_bank
from regjepa.linear import LinearScorer, LinearTransform, linear_objective
from regjepa.protocol import file_sha256, load_split_plan, write_json


def compare_values(actual: object, expected: object) -> None:
    """Compare JSON structures, allowing only tiny differences in numeric reductions."""
    if isinstance(expected, dict):
        assert isinstance(actual, dict) and actual.keys() == expected.keys()
        for key in expected:
            compare_values(actual[key], expected[key])
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected)
        for left, right in zip(actual, expected):
            compare_values(left, right)
    elif isinstance(expected, float):
        np.testing.assert_allclose(actual, expected, rtol=1e-11, atol=1e-11)
    else:
        assert actual == expected, (actual, expected)


def check_convergence(diagnostic: dict, optimizer: dict) -> None:
    assert diagnostic["converged"] is True
    assert 0 <= diagnostic["iterations"] <= optimizer["max_iterations"]
    assert np.isfinite(diagnostic["objective"])
    assert 0 <= diagnostic["gradient_inf_norm"] <= optimizer["accepted_gradient"]


def validate_fit(values: np.ndarray, queries: list[dict], block: str, outer: dict,
                 training: dict, config: dict, output: Path) -> tuple[dict, dict]:
    tag = f"{block}_fold{outer['fold']}_budget{training['budget']}_seed{training['subset_seed']}"
    result = json.loads((output / f"{tag}.json").read_text(encoding="utf-8"))
    assert result["block"] == block and result["outer_fold"] == outer["fold"]
    assert result["budget"] == training["budget"] and result["subset_seed"] == training["subset_seed"]
    assert result["actual_training_labels"] == training["actual_labels"]
    choices = penalty_candidates(block, config)
    assert len(result["search"]) == len(choices)
    diagnostics = []
    for choice, search in zip(choices, result["search"]):
        assert search["penalties"] == choice
        assert search["validation_molecules"] == training["actual_labels"]
        assert np.isfinite(search["validation_nll"])
        assert len(search["fits"]) == len(training["inner_folds"])
        diagnostics.extend(search["fits"])
    inner_fits = len(diagnostics)
    diagnostics.append(result["final_fit"])
    for diagnostic in diagnostics:
        check_convergence(diagnostic, config["optimizer"])
    selected = min(result["search"], key=lambda search: search["validation_nll"])["penalties"]
    assert result["selected_penalties"] == selected
    assert result["regularization_grid_edges"] == grid_edge_diagnostics(selected, config)
    model_path = output / f"{tag}.npz"
    assert file_sha256(model_path) == result["model_sha256"]
    with np.load(model_path, allow_pickle=False) as archive:
        assert set(archive.files) == {"keep", "scale", "weights"}
        keep, scale, weights = (archive[key].copy() for key in ["keep", "scale", "weights"])
    assert keep.dtype == np.bool_ and np.isfinite(scale).all() and np.all(scale > 0)
    assert weights.shape == (int(keep.sum()),) and np.isfinite(weights).all()
    assert result["final_fit"]["active_coefficients"] == len(weights)
    x_train, _, train_table = candidate_subset(values, queries, training["train_compounds"])
    x_test, test_queries, test_table = candidate_subset(values, queries, outer["test_compounds"])
    transform = LinearTransform.fit(x_train, train_table, standardize=block != "morgan")
    np.testing.assert_array_equal(keep, transform.keep)
    np.testing.assert_array_equal(scale, transform.scale)
    scorer = LinearScorer(LinearTransform(keep, scale, block != "morgan"), weights, result["final_fit"])
    objective, gradient = linear_objective(weights, scorer.transform.apply(x_train, train_table),
                                           train_table, penalty_vector(selected, values.shape[1])[keep])
    gradient_norm = float(np.max(np.abs(gradient))) if len(gradient) else 0.0
    np.testing.assert_allclose(objective, result["final_fit"]["objective"], rtol=1e-11, atol=1e-11)
    np.testing.assert_allclose(gradient_norm, result["final_fit"]["gradient_inf_norm"], rtol=1e-9, atol=1e-11)
    assert gradient_norm <= config["optimizer"]["accepted_gradient"] + 1e-11
    scores = scorer.predict(x_test, test_table)
    assert len(result["predictions"]) == len(test_queries)
    for prediction, query in zip(result["predictions"], test_queries):
        for key in ["query_id", "compound_id", "structure_key", "reported_positive"]:
            assert prediction[key] == query[key]
    saved_scores = np.array([prediction["score"] for prediction in result["predictions"]])
    np.testing.assert_array_equal(scores, saved_scores)
    compare_values(result["molecular_metrics"], ranking_records(scores, test_table))
    if block == "jepa_plus_qm":
        compare_values(result["component_scores"], component_score_diagnostics(scorer, x_test, test_table))
    else:
        assert result["component_scores"] is None
    evidence = {
        "tag": tag, "model_sha256": result["model_sha256"], "inner_fits": inner_fits,
        "test_molecules": len(test_table.rows), "test_rows": len(test_queries),
        "max_saved_score_replay_error": float(np.max(np.abs(scores - saved_scores))),
        "objective_replay_error": abs(objective - result["final_fit"]["objective"]),
        "gradient_replay_error": abs(gradient_norm - result["final_fit"]["gradient_inf_norm"]),
        "max_iterations_in_task": max(item["iterations"] for item in diagnostics),
        "max_gradient_in_task": max(item["gradient_inf_norm"] for item in diagnostics),
        "fits_over_original_1000_cap": sum(item["iterations"] > 1000 for item in diagnostics),
    }
    return result, evidence


def validate_results(split_path: Path, output: Path) -> dict:
    plan = load_split_plan(ROOT, split_path)
    config = plan["config"]
    assert json.loads((output / "split_manifest.json").read_text()) == plan
    summary = json.loads((output / "summary.json").read_text())
    assert summary["proposal_id"] == plan["proposal_id"]
    assert summary["split_sha256"] == file_sha256(split_path)
    accumulated = defaultdict(list)
    evidence, fit_diagnostics = [], []
    for block in config["blocks"]:
        bank = load_feature_bank(ROOT, cohort=config["cohort"], blocks=[block])
        for outer in plan["outer_folds"]:
            for training in outer["training_sets"]:
                result, checked = validate_fit(bank.X[block], bank.queries, block, outer, training, config, output)
                evidence.append(checked)
                accumulated[(block, str(training["budget"]))].extend(result["molecular_metrics"])
                fit_diagnostics.append({key: result[key] for key in ["block", "outer_fold", "budget", "subset_seed",
                                                                     "selected_penalties", "regularization_grid_edges", "component_scores"]})
        print(f"Validated all saved {block} readouts", flush=True)
    tags = {item["tag"] for item in evidence}
    assert len(tags) == summary["completed_final_fits"]
    assert {path.stem for path in output.glob("*.npz")} == tags
    assert {path.stem for path in output.glob("*.json")} == tags | {"split_manifest", "summary"}
    compare_values(summary["fit_diagnostics"], fit_diagnostics)
    expected_ids = set(plan["compound_to_group"])
    summaries, molecular = {}, {}
    for (block, budget), records in accumulated.items():
        averaged = average_repetitions(records, expected_ids)
        repetitions = 1 if budget == "full" else len(config["subset_seeds"])
        assert all(item["repetitions"] == repetitions for item in averaged)
        summaries[f"{block}:budget{budget}"] = {"metrics": summarize_records(averaged), "molecular_metrics": averaged}
        molecular[(block, budget)] = {item["compound_id"]: item["top1"] for item in averaged}
    compare_values(summary["results"], summaries)
    pairs = [("atom_jepa", "random_backbone"), ("atom_jepa", "morgan"), ("atom_jepa", "electronic"),
             ("jepa_plus_qm", "atom_jepa"), ("jepa_plus_qm", "electronic")]
    comparisons = []
    for budget in [str(item) for item in config["label_budgets"]] + ["full"]:
        for left, right in pairs:
            if (left, budget) in molecular and (right, budget) in molecular:
                comparison = paired_bootstrap(molecular[(left, budget)], molecular[(right, budget)], plan["compound_to_group"],
                                              seed=config["bootstrap_seed"], repetitions=config["bootstrap_repetitions"])
                comparisons.append({"left": left, "right": right, "budget": budget, **comparison})
    compare_values(summary["paired_top1_comparisons"], comparisons)
    planned_inner = sum(len(outer["training_sets"]) for outer in plan["outer_folds"]) * config["inner_folds"] * sum(len(penalty_candidates(block, config)) for block in config["blocks"])
    assert sum(item["inner_fits"] for item in evidence) == planned_inner
    return {
        "schema_version": 1, "proposal_id": plan["proposal_id"], "status": "passed",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "summary_sha256": file_sha256(output / "summary.json"), "frozen_artifact_sha256": file_sha256(split_path),
        "completed_final_fits": len(evidence), "completed_inner_fits": planned_inner,
        "molecules_per_setting": len(expected_ids), "leakage_groups": len(set(plan["compound_to_group"].values())),
        "result_settings": len(summaries), "paired_comparisons_recomputed": len(comparisons),
        "all_model_hashes_valid": True, "all_train_only_transforms_reproduced": True,
        "all_saved_query_scores_reproduced_exactly": True, "all_molecular_metrics_and_summaries_reproduced": True,
        "all_component_diagnostics_reproduced": True, "all_fits_converged": True,
        "largest_iteration_count": max(item["max_iterations_in_task"] for item in evidence),
        "largest_gradient_inf_norm": max(item["max_gradient_in_task"] for item in evidence),
        "fits_over_original_1000_cap": sum(item["fits_over_original_1000_cap"] for item in evidence),
        "max_score_replay_error": max(item["max_saved_score_replay_error"] for item in evidence),
        "max_objective_replay_error": max(item["objective_replay_error"] for item in evidence),
        "max_gradient_replay_error": max(item["gradient_replay_error"] for item in evidence),
        "models_refitted": 0, "scope": "Same-runtime artifact replay, not an independent scientific replication",
        "fit_checks": evidence,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="Existing completed run directory")
    parser.add_argument("--report", type=Path, required=True, help="New JSON audit report outside the run directory")
    args = parser.parse_args()
    output = args.output.resolve()
    report_path = args.report.resolve()
    if output == report_path or output in report_path.parents:
        raise ValueError("Keep validation reports outside the immutable experiment directory")
    if report_path.exists() or report_path.with_suffix(".md").exists():
        raise FileExistsError("Validation evidence is never overwritten")
    report = validate_results(args.splits.resolve(), output)
    write_json(report_path, report, exclusive=True)
    lines = ["# Saved-model replay validation", "", f"Proposal: `{report['proposal_id']}`. Status: **passed**.", "",
             f"Verified {report['completed_inner_fits']} inner fits and {report['completed_final_fits']} final fits; all converged.",
             f"Every saved model hash, training-only transform and query-aligned score reproduced; maximum score error {report['max_score_replay_error']}.",
             f"All {report['result_settings']} molecular summaries and {report['paired_comparisons_recomputed']} paired bootstrap comparisons reproduced.",
             f"Each setting covers {report['molecules_per_setting']} molecules in {report['leakage_groups']} leakage components.",
             f"Largest iteration count: {report['largest_iteration_count']}; fits exceeding the original 1,000-iteration cap: {report['fits_over_original_1000_cap']}.",
             f"Largest gradient infinity norm: {report['largest_gradient_inf_norm']:.12g}.", "",
             "No models were refitted. Numerical checks use the same reviewed implementation and runtime; they are an artifact audit, not an independent replication.",
             "The JSON sibling includes per-fit evidence and hashes.", "", "## Not done", "",
             "Independent scientific review, pretraining-overlap audit and stronger external/series validation remain outstanding.", ""]
    report_path.with_suffix(".md").write_text("\n".join(lines), encoding="utf-8", newline="\n")
    print(json.dumps({key: value for key, value in report.items() if key != "fit_checks"}, indent=2))


if __name__ == "__main__":
    main()
