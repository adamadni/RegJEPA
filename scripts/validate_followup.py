"""Replay saved follow-up results; never refit or claim an independent scientific replication."""
import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from regjepa.evaluation import ranking_records
from regjepa.experiment import candidate_subset
from regjepa.feature_io import load_feature_bank
from regjepa.linear import linear_objective
from regjepa.protocol import file_sha256, write_json
from regjepa_followup.experiment import choices, summarize, vector, workload
from regjepa_followup.features import materialize
from regjepa_followup.models import Scorer, Transform, nonlinear_objective
from regjepa_followup.protocol import load


def check_diagnostics(diagnostic: dict, head: str, config: dict) -> int:
    expected_seeds = config["mlp_restart_seeds"] if head == "mlp" else [None]
    attempts = diagnostic["attempts"]
    assert [a["seed"] for a in attempts] == expected_seeds
    for a in attempts:
        assert a["converged"] and np.isfinite(a["objective"])
        assert 0 <= a["iterations"] <= config["max_iterations"]
        assert 0 <= a["gradient_inf_norm"] <= config["accepted_gradient"]
    assert diagnostic["selected_restart"] == min(range(len(attempts)), key=lambda i: attempts[i]["objective"])
    return len(attempts)


def validate(root: Path, split: Path, output: Path) -> dict:
    plan = load(root, split)
    config = plan["config"]
    assert json.loads((output / "split_manifest.json").read_text()) == plan
    result = json.loads((output / "summary.json").read_text())
    assert result["frozen_sha256"] == file_sha256(split)
    assert result["proposal_id"] == plan["proposal_id"]
    bank = load_feature_bank(root, cohort=config["cohort"], blocks=["electronic", "atom_jepa", "random_backbone", "morgan"])
    accumulated = defaultdict(list)
    tags, fit_diagnostics = [], []
    inner_attempts = final_attempts = 0
    for arm in config["arms"]:
        values, groups, morgan = materialize(bank, arm)
        for outer in plan["outer_folds"]:
            for training in outer["training_sets"]:
                if training["budget"] not in arm["budgets"]:
                    continue
                tag = f"{arm['id']}_fold{outer['fold']}_budget{training['budget']}_seed{training['subset_seed']}"
                record = json.loads((output / f"{tag}.json").read_text())
                for key, expected in {"arm": arm["id"], "head": arm["head"], "outer_fold": outer["fold"],
                                      "budget": training["budget"], "subset_seed": training["subset_seed"],
                                      "actual_training_labels": training["actual_labels"]}.items():
                    assert record[key] == expected
                assert len(record["search"]) == len(choices(arm, config))
                for choice, search in zip(choices(arm, config), record["search"]):
                    assert search["penalties"] == choice and np.isfinite(search["validation_nll"])
                    assert search["validation_molecules"] == training["actual_labels"]
                    assert len(search["fits"]) == len(training["inner_folds"])
                    inner_attempts += sum(check_diagnostics(d, arm["head"], config["optimizer"]) for d in search["fits"])
                selected = min(record["search"], key=lambda s: s["validation_nll"])["penalties"]
                assert record["selected_penalties"] == selected
                final_attempts += check_diagnostics(record["final_fit"], arm["head"], config["optimizer"])
                assert record["model_sha256"] == file_sha256(output / f"{tag}.npz")
                x_train, _, train_table = candidate_subset(values, bank.queries, training["train_compounds"])
                x_test, queries, test_table = candidate_subset(values, bank.queries, outer["test_compounds"])
                expected = Transform.fit(x_train, train_table, arm["head"], morgan)
                with np.load(output / f"{tag}.npz", allow_pickle=False) as saved:
                    assert set(saved.files) == {"parameters", "keep", "mean", "scale", "hidden"}
                    for name in ["keep", "mean", "scale"]:
                        np.testing.assert_array_equal(saved[name], getattr(expected, name))
                    hidden = int(saved["hidden"])
                    assert hidden == (config["optimizer"]["mlp_hidden"] if arm["head"] == "mlp" else 0)
                    model = Scorer(expected, saved["parameters"].copy(), hidden, record["final_fit"])
                assert model.parameters.shape == ((int(expected.keep.sum()) * hidden + 2 * hidden) if hidden else int(expected.keep.sum()),)
                data = expected.apply(x_train, train_table)
                penalties = vector(selected, groups)
                objective = (nonlinear_objective(model.parameters, data, train_table, penalties[0], hidden)
                             if hidden else linear_objective(model.parameters, data, train_table, penalties[expected.keep]))
                final = record["final_fit"]["attempts"][record["final_fit"]["selected_restart"]]
                gradient = float(np.max(np.abs(objective[1]))) if len(objective[1]) else 0.0
                np.testing.assert_allclose([objective[0], gradient], [final["objective"], final["gradient_inf_norm"]], atol=1e-10, rtol=1e-9)
                assert gradient <= config["optimizer"]["accepted_gradient"] + 1e-10
                assert record["final_fit"]["parameter_count"] == len(model.parameters)
                assert record["final_fit"]["active_features"] == int(expected.keep.sum())
                expected_edges = {k: {"value": v, "at_minimum": v == min(config["l2_grid"]), "at_maximum": v == max(config["l2_grid"])} for k, v in selected.items()}
                assert record["grid_edges"] == expected_edges
                scores = model.predict(x_test, test_table)
                assert len(record["predictions"]) == len(queries)
                for p, q in zip(record["predictions"], queries):
                    assert all(p[k] == q[k] for k in ["query_id", "compound_id", "structure_key", "reported_positive"])
                np.testing.assert_array_equal(scores, [p["score"] for p in record["predictions"]])
                assert ranking_records(scores, test_table) == record["molecular_metrics"]
                accumulated[(arm["id"], str(training["budget"]))].extend(record["molecular_metrics"])
                tags.append(tag)
                fit_diagnostics.append({"tag": tag, "arm": arm["id"], "budget": training["budget"],
                                        "grid_edges": record["grid_edges"], "final_fit": record["final_fit"]})
        print(f"Validated {arm['id']}", flush=True)
    expected = workload(plan)
    assert len(tags) == result["completed_final_models"] == expected["final_models"]
    assert result["workload"] == expected
    assert inner_attempts == expected["inner_optimizer_attempts"] and final_attempts == expected["final_optimizer_attempts"]
    assert {p.stem for p in output.glob("*.npz")} == set(tags)
    assert {p.stem for p in output.glob("*.json")} == set(tags) | {"split_manifest", "summary"}
    assert result["fit_diagnostics"] == fit_diagnostics
    recalculated = summarize(accumulated, plan)
    assert result["results"] == recalculated["results"]
    assert result["paired_comparisons"] == recalculated["paired_comparisons"]
    assert result["descriptive_diagnostics"] == recalculated["descriptive_diagnostics"]
    return {"proposal_id": plan["proposal_id"], "status": "passed", "final_models_replayed": len(tags),
            "inner_optimizer_diagnostics_checked": inner_attempts, "final_optimizer_diagnostics_checked": final_attempts,
            "selected_final_gradients_replayed": len(tags), "models_refitted": 0,
            "summary_sha256": file_sha256(output / "summary.json"),
            "scope": "Same-runtime artifact replay; nonselected and inner gradients checked from saved diagnostics only."}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    output, report = args.output.resolve(), args.report.resolve()
    if output in report.parents or report.exists() or report.with_suffix(".md").exists():
        raise ValueError("Use new audit paths outside the immutable run directory")
    result = validate(ROOT, args.splits.resolve(), output)
    write_json(report, result, exclusive=True)
    text = (f"# Follow-up artifact validation\n\nProposal `{result['proposal_id']}`: **passed**.\n\n"
            f"Replayed {result['final_models_replayed']} selected final readouts and their gradients without refitting. "
            f"Checked diagnostics for {result['inner_optimizer_diagnostics_checked']} inner optimizer attempts and "
            f"{result['final_optimizer_diagnostics_checked']} final attempts. All molecular summaries and paired intervals reproduced.\n\n"
            "Inner and nonselected restart gradients are checked from saved diagnostics only. This is same-runtime replay, not independent scientific replication.\n\n"
            "## Not done\n\nIndependent scientific replication and external confirmation.\n")
    report.with_suffix(".md").write_text(text, encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
