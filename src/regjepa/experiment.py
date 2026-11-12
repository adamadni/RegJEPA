"""Nested linear probes; evaluation does not establish a JEPA-objective or mechanism claim."""
from __future__ import annotations

from collections import defaultdict
from itertools import product
from pathlib import Path

import numpy as np

from .evaluation import CandidateTable, choice_loss, paired_bootstrap, ranking_records, summarize_records
from .feature_io import load_feature_bank
from .linear import LinearScorer, fit_linear, linear_scores
from .protocol import file_sha256, load_split_plan, write_json


def penalty_candidates(block: str, config: dict) -> list[dict]:
    grid = sorted(config["l2_grid"], reverse=True)
    if block == "jepa_plus_qm" and config["combined_penalties"] == "separate_jepa_qm":
        return [{"jepa": a, "qm": b} for a, b in product(grid, grid)]
    return [{"shared": value} for value in grid]


def penalty_vector(choice: dict, width: int) -> np.ndarray:
    if "shared" in choice:
        return np.full(width, choice["shared"], dtype=np.float64)
    if width != 774:
        raise ValueError("Separate JEPA/QM penalties require the literal 768+6 feature layout")
    return np.concatenate([np.full(768, choice["jepa"]), np.full(6, choice["qm"])])


def grid_edge_diagnostics(choice: dict, config: dict) -> dict:
    low, high = min(config["l2_grid"]), max(config["l2_grid"])
    return {name: {"value": value, "at_minimum": value == low, "at_maximum": value == high}
            for name, value in choice.items()}


def component_score_diagnostics(model: LinearScorer, values: np.ndarray, table: CandidateTable) -> dict:
    """Relative score amplitudes, not unique information attribution or physical rates."""
    if values.shape[1] != 774:
        raise ValueError("JEPA/QM score decomposition requires the original 768+6 columns")
    features = model.transform.apply(values, table)
    original_columns = np.flatnonzero(model.transform.keep)
    jepa_columns = original_columns < 768
    jepa = linear_scores(features[:, jepa_columns], model.weights[jepa_columns])
    qm = linear_scores(features[:, ~jepa_columns], model.weights[~jepa_columns])
    mean_product = lambda a, b: float(np.mean([np.mean(a[indices] * b[indices]) for indices in table.rows]))
    jepa_rms = float(np.sqrt(mean_product(jepa, jepa)))
    qm_rms = float(np.sqrt(mean_product(qm, qm)))
    total = jepa_rms + qm_rms
    return {"test_molecules": len(table.rows), "jepa_rms": jepa_rms, "qm_rms": qm_rms,
            "total_score_rms": float(np.sqrt(mean_product(jepa + qm, jepa + qm))),
            "jepa_qm_mean_product": mean_product(jepa, qm),
            "qm_relative_rms": qm_rms / total if total > 0 else None,
            "interpretation": "equal-molecule score RMS; correlated components prevent unique information attribution"}


def candidate_subset(values: np.ndarray, queries: list[dict], compounds: list[str]) -> tuple[np.ndarray, list[dict], CandidateTable]:
    selected = set(compounds)
    if len(selected) != len(compounds):
        raise ValueError("Repeated training or validation record ID")
    indices = [i for i, q in enumerate(queries) if q["compound_id"] in selected]
    subset = [queries[i] for i in indices]
    if {q["compound_id"] for q in subset} != selected:
        raise ValueError("Requested record missing from feature queries")
    table = CandidateTable.from_queries(subset)
    return values[indices], subset, table


def fit_nested_split(values: np.ndarray, queries: list[dict], outer: dict, training: dict,
                     block: str, config: dict) -> tuple[dict, object]:
    """Tune only inside the selected label budget, then predict the untouched outer set."""
    selected_ids = set(training["train_compounds"])
    if not selected_ids <= set(outer["train_compounds"]) or selected_ids & set(outer["test_compounds"]):
        raise ValueError("Training budget escapes outer training set")
    validation_ids = [name for fold in training["inner_folds"] for name in fold["test_compounds"]]
    if len(validation_ids) != len(selected_ids) or set(validation_ids) != selected_ids:
        raise ValueError("Inner validation labels are missing or repeated")
    candidates = penalty_candidates(block, config)
    searches = []
    for choice in candidates:
        weighted_loss, molecules = 0.0, 0
        diagnostics = []
        for inner in training["inner_folds"]:
            allowed = set(training["train_compounds"])
            if (set(inner["train_compounds"]) | set(inner["test_compounds"]) != allowed or
                    set(inner["train_compounds"]) & set(inner["test_compounds"])):
                raise ValueError("Inner selection uses labels outside its budget or overlaps validation")
            x_train, _, train_table = candidate_subset(values, queries, inner["train_compounds"])
            x_val, _, val_table = candidate_subset(values, queries, inner["test_compounds"])
            scorer = fit_linear(x_train, train_table, penalty_vector(choice, values.shape[1]),
                                standardize=block != "morgan", **config["optimizer"])
            loss, _ = choice_loss(scorer.predict(x_val, val_table), val_table)
            weighted_loss += loss * len(val_table.rows)
            molecules += len(val_table.rows)
            diagnostics.append(scorer.diagnostics)
        if molecules != training["actual_labels"]:
            raise ValueError("Inner validation labels are missing or repeated")
        searches.append({"penalties": choice, "validation_nll": weighted_loss / molecules,
                         "validation_molecules": molecules, "fits": diagnostics})
    # Candidate order prefers stronger shrinkage on an exact validation-loss tie.
    best = min(range(len(searches)), key=lambda i: searches[i]["validation_nll"])
    selected = searches[best]["penalties"]
    x_train, _, train_table = candidate_subset(values, queries, training["train_compounds"])
    x_test, test_queries, test_table = candidate_subset(values, queries, outer["test_compounds"])
    scorer = fit_linear(x_train, train_table, penalty_vector(selected, values.shape[1]),
                        standardize=block != "morgan", **config["optimizer"])
    scores = scorer.predict(x_test, test_table)
    per_molecule = ranking_records(scores, test_table)
    predictions = [{"query_id": q["query_id"], "compound_id": q["compound_id"],
                    "structure_key": q["structure_key"], "reported_positive": q["reported_positive"],
                    "score": float(score)} for q, score in zip(test_queries, scores)]
    return {"block": block, "outer_fold": outer["fold"], "budget": training["budget"],
            "subset_seed": training["subset_seed"], "actual_training_labels": len(train_table.rows),
            "selected_penalties": selected, "search": searches, "final_fit": scorer.diagnostics,
            "regularization_grid_edges": grid_edge_diagnostics(selected, config),
            "component_scores": component_score_diagnostics(scorer, x_test, test_table) if block == "jepa_plus_qm" else None,
            "molecular_metrics": per_molecule, "predictions": predictions}, scorer


def average_repetitions(records: list[dict], expected_ids: set[str]) -> list[dict]:
    by_compound: dict[str, list[dict]] = defaultdict(list)
    for record in records:
        by_compound[record["compound_id"]].append(record)
    if set(by_compound) != expected_ids:
        raise ValueError("Out-of-fold molecular coverage differs between methods")
    if len({len(items) for items in by_compound.values()}) != 1:
        raise ValueError("Molecules have inconsistent repetition counts")
    metrics = ["top1", "top2", "reciprocal_rank", "nll", "uniform_top1", "uniform_nll", "information_gain_nats"]
    result = []
    for compound in sorted(by_compound):
        items = by_compound[compound]
        first = items[0]
        if any((r["structure_key"], r["n_candidates"]) != (first["structure_key"], first["n_candidates"]) for r in items):
            raise ValueError("Candidate sets changed between repetitions")
        result.append({"compound_id": compound, "structure_key": first["structure_key"],
                       "n_candidates": first["n_candidates"], "repetitions": len(items),
                       **{name: float(np.mean([r[name] for r in items])) for name in metrics}})
    return result


def run_experiment(root: Path, split_path: Path, output: Path) -> dict:
    plan = load_split_plan(root, split_path)
    config = plan["config"]
    if output.exists():
        raise FileExistsError("Use a new output directory; evaluation results are never overwritten")
    output.mkdir(parents=True)
    write_json(output / "split_manifest.json", plan, exclusive=True)
    accumulated: dict[tuple[str, str], list[dict]] = defaultdict(list)
    completed = []
    fit_diagnostics = []
    for block in config["blocks"]:
        bank = load_feature_bank(root, cohort=config["cohort"], blocks=[block])
        for outer in plan["outer_folds"]:
            for training in outer["training_sets"]:
                result, model = fit_nested_split(bank.X[block], bank.queries, outer, training, block, config)
                tag = f"{block}_fold{outer['fold']}_budget{training['budget']}_seed{training['subset_seed']}"
                np.savez_compressed(output / f"{tag}.npz", weights=model.weights,
                                    keep=model.transform.keep, scale=model.transform.scale)
                result["model_sha256"] = file_sha256(output / f"{tag}.npz")
                write_json(output / f"{tag}.json", result, exclusive=True)
                accumulated[(block, str(training["budget"]))].extend(result["molecular_metrics"])
                fit_diagnostics.append({key: result[key] for key in ["block", "outer_fold", "budget", "subset_seed",
                                                                    "selected_penalties", "regularization_grid_edges", "component_scores"]})
                completed.append(tag)
                print(f"Completed {tag}; active coefficients={model.diagnostics['active_coefficients']}", flush=True)
    expected = set(plan["compound_to_group"])
    summaries, molecular = {}, {}
    for (block, budget), records in accumulated.items():
        averaged = average_repetitions(records, expected)
        key = f"{block}:budget{budget}"
        molecular[(block, budget)] = {r["compound_id"]: r["top1"] for r in averaged}
        summaries[key] = {"metrics": summarize_records(averaged), "molecular_metrics": averaged}
    comparisons = []
    pairs = [("atom_jepa", "random_backbone"), ("atom_jepa", "morgan"), ("atom_jepa", "electronic"),
             ("jepa_plus_qm", "atom_jepa"), ("jepa_plus_qm", "electronic")]
    for budget in [str(b) for b in config["label_budgets"]] + ["full"]:
        for left, right in pairs:
            if (left, budget) not in molecular or (right, budget) not in molecular:
                continue
            comparison = paired_bootstrap(molecular[(left, budget)], molecular[(right, budget)], plan["compound_to_group"],
                                          seed=config["bootstrap_seed"], repetitions=config["bootstrap_repetitions"])
            comparisons.append({"left": left, "right": right, "budget": budget, **comparison})
    summary = {"proposal_id": plan["proposal_id"], "split_sha256": file_sha256(split_path),
               "completed_final_fits": len(completed), "results": summaries, "paired_top1_comparisons": comparisons,
               "fit_diagnostics": fit_diagnostics,
               "interpretation": "Frozen-feature linear readouts; source-annotation prediction, not product fractions or mechanism inference.",
               "pretraining_overlap_status": "Uni-Mol corpus identified; benchmark membership has not been audited."}
    write_json(output / "summary.json", summary, exclusive=True)
    lines = ["# Linear frozen-feature evaluation", "", f"Split proposal: `{plan['proposal_id']}`.", "",
             "| Block / label budget | Molecules | Top-1 | Top-2 | Top-2 (k>2) | MRR | NLL | Uniform NLL | Information gain |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for name, value in summaries.items():
        m = value["metrics"]
        nontrivial = f"{m['top2_nontrivial']:.4f}" if m["top2_nontrivial"] is not None else "n/a"
        lines.append(f"| {name} | {m['molecules']} | {m['top1']:.4f} | {m['top2']:.4f} | {nontrivial} | {m['reciprocal_rank']:.4f} | {m['nll']:.4f} | {m['uniform_nll']:.4f} | {m['information_gain_nats']:.4f} |")
    lines.extend(["", "## Paired top-1 differences", "", "| Left / right | Budget | Difference | Conditional 95% interval |",
                  "|---|---|---:|---|"])
    for comparison in comparisons:
        low, high = comparison["ci95"]
        lines.append(f"| {comparison['left']} / {comparison['right']} | {comparison['budget']} | {comparison['mean_delta']:.4f} | [{low:.4f}, {high:.4f}] |")
    edge_counts: dict[tuple[str, str], list[bool]] = defaultdict(list)
    for fit in fit_diagnostics:
        edge_counts[(fit["block"], str(fit["budget"]))].append(
            any(value["at_minimum"] or value["at_maximum"] for value in fit["regularization_grid_edges"].values()))
    lines.extend(["", "## Regularization-grid diagnostics", "", "| Block / budget | Final fits selecting a grid edge |", "|---|---:|"])
    for (block, budget), flags in edge_counts.items():
        lines.append(f"| {block} / {budget} | {sum(flags)} / {len(flags)} |")
    lines.extend(["", "Edge selections are diagnostic; the grid is not automatically changed using outer-test performance.",
                  "", "## JEPA + QM score components", "",
                  "Relative RMS measures score amplitude, not independent information or product fractions; covariance is retained in JSON.", "",
                  "| Fold / budget / subset seed | JEPA RMS | QM RMS | QM relative RMS |", "|---|---:|---:|---:|"])
    for fit in fit_diagnostics:
        diagnostic = fit["component_scores"]
        if diagnostic is not None:
            share = f"{diagnostic['qm_relative_rms']:.4f}" if diagnostic["qm_relative_rms"] is not None else "n/a"
            lines.append(f"| {fit['outer_fold']} / {fit['budget']} / {fit['subset_seed']} | {diagnostic['jepa_rms']:.4f} | {diagnostic['qm_rms']:.4f} | {share} |")
    lines.extend(["", "Exact ties receive expected metrics. Repetition metrics are averaged per molecule; logits are not ensembled.",
                  "Paired confidence intervals resample leakage components and are conditional on the fitted models.",
                  "Intervals are computed for the primary top-1 comparison; secondary metrics are descriptive.",
                  "All source records are retained. Bibliographic uncertainty and pretraining overlap limit generalization claims.",
                  "", "Full diagnostics, molecular predictions and paired intervals are in summary.json and per-fit JSON files.", ""])
    (output / "summary.md").write_text("\n".join(lines), encoding="utf-8", newline="\n")
    return summary
