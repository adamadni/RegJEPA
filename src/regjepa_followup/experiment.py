"""Nested diagnostic readouts; results cannot establish causal descriptor importance or mechanisms."""
from __future__ import annotations

from collections import defaultdict
from itertools import product
from pathlib import Path

import numpy as np

from regjepa.evaluation import choice_loss, paired_bootstrap, ranking_records, summarize_records
from regjepa.experiment import average_repetitions, candidate_subset
from regjepa.feature_io import load_feature_bank
from regjepa.protocol import file_sha256, write_json

from .features import materialize
from .models import Scorer, fit
from .protocol import load


def choices(arm: dict, config: dict) -> list[dict]:
    grid = sorted(config["l2_grid"], reverse=True)
    if arm["penalty_mode"] == "shared":
        return [{"shared": penalty} for penalty in grid]
    return [dict(zip(["representation", "qm"], values)) for values in product(grid, repeat=2)]


def vector(choice: dict, groups: list[str]) -> np.ndarray:
    if "shared" in choice:
        return np.full(len(groups), choice["shared"])
    return np.array([choice[group] for group in groups], dtype=np.float64)


def workload(plan: dict) -> dict:
    inner_models = final_models = inner_attempts = final_attempts = 0
    arms = []
    for arm in plan["config"]["arms"]:
        settings = [t for f in plan["outer_folds"] for t in f["training_sets"] if t["budget"] in arm["budgets"]]
        candidates = len(choices(arm, plan["config"]))
        inners = sum(len(t["inner_folds"]) for t in settings) * candidates
        restarts = len(plan["config"]["optimizer"]["mlp_restart_seeds"]) if arm["head"] == "mlp" else 1
        inner_models += inners; final_models += len(settings)
        inner_attempts += inners * restarts; final_attempts += len(settings) * restarts
        arms.append({"id": arm["id"], "head": arm["head"], "settings": len(settings),
                     "width": sum(len(p["columns"]) for p in arm["parts"]),
                     "penalty_candidates": candidates, "inner_optimizer_attempts": inners * restarts,
                     "final_optimizer_attempts": len(settings) * restarts})
    return {"arms": arms, "inner_models": inner_models, "final_models": final_models,
            "inner_optimizer_attempts": inner_attempts, "final_optimizer_attempts": final_attempts}


def fit_setting(values: np.ndarray, queries: list[dict], groups: list[str], morgan: np.ndarray,
                arm: dict, outer: dict, training: dict, config: dict) -> tuple[dict, Scorer]:
    selected_ids = set(training["train_compounds"])
    if not selected_ids <= set(outer["train_compounds"]) or selected_ids & set(outer["test_compounds"]):
        raise ValueError("Label budget escapes or overlaps the outer split")
    validation = [c for f in training["inner_folds"] for c in f["test_compounds"]]
    if len(validation) != len(selected_ids) or set(validation) != selected_ids:
        raise ValueError("Inner held-out labels are missing or duplicated")
    for inner in training["inner_folds"]:
        a, b = set(inner["train_compounds"]), set(inner["test_compounds"])
        if a | b != selected_ids or a & b:
            raise ValueError("Inner split escapes or overlaps its label budget")
    searches = []
    for choice in choices(arm, config):
        losses, count, diagnostics = 0.0, 0, []
        for inner in training["inner_folds"]:
            x_train, _, train_table = candidate_subset(values, queries, inner["train_compounds"])
            x_val, _, val_table = candidate_subset(values, queries, inner["test_compounds"])
            scorer = fit(x_train, train_table, morgan, arm["head"], vector(choice, groups), config["optimizer"])
            loss, _ = choice_loss(scorer.predict(x_val, val_table), val_table)
            losses += loss * len(val_table.rows)
            count += len(val_table.rows)
            diagnostics.append(scorer.diagnostics)
        if count != training["actual_labels"]:
            raise ValueError("Inner coverage differs from the exact label budget")
        searches.append({"penalties": choice, "validation_nll": losses / count,
                         "validation_molecules": count, "fits": diagnostics})
    best = min(range(len(searches)), key=lambda i: searches[i]["validation_nll"])
    selected = searches[best]["penalties"]
    x_train, _, train_table = candidate_subset(values, queries, training["train_compounds"])
    x_test, test_queries, test_table = candidate_subset(values, queries, outer["test_compounds"])
    scorer = fit(x_train, train_table, morgan, arm["head"], vector(selected, groups), config["optimizer"])
    scores = scorer.predict(x_test, test_table)
    low, high = min(config["l2_grid"]), max(config["l2_grid"])
    result = {"arm": arm["id"], "head": arm["head"], "outer_fold": outer["fold"],
              "budget": training["budget"], "subset_seed": training["subset_seed"],
              "actual_training_labels": len(train_table.rows), "selected_penalties": selected,
              "search": searches, "final_fit": scorer.diagnostics,
              "grid_edges": {k: {"value": v, "at_minimum": v == low, "at_maximum": v == high} for k, v in selected.items()},
              "active_columns_by_group": {g: sum(keep and group == g for keep, group in zip(scorer.transform.keep.tolist(), groups)) for g in sorted(set(groups))},
              "molecular_metrics": ranking_records(scores, test_table),
              "predictions": [{"query_id": q["query_id"], "compound_id": q["compound_id"],
                               "structure_key": q["structure_key"], "reported_positive": q["reported_positive"],
                               "score": float(s)} for q, s in zip(test_queries, scores)]}
    return result, scorer


def summarize(accumulated: dict, plan: dict) -> dict:
    expected = set(plan["compound_to_group"])
    results, descriptive_diagnostics = {}, {}
    for (arm, budget), records in accumulated.items():
        averaged = average_repetitions(records, expected)
        # Coverage is additionally checked against the exact configured subset seeds.
        expected_repetitions = 1 if budget == "full" else len({t["subset_seed"] for t in plan["outer_folds"][0]["training_sets"] if str(t["budget"]) == budget})
        if any(r["repetitions"] != expected_repetitions for r in averaged):
            raise ValueError("Unexpected molecular repetition count")
        results[f"{arm}:budget{budget}"] = {"metrics": summarize_records(averaged), "molecular_metrics": averaged}
        descriptive_diagnostics[f"{arm}:budget{budget}"] = {
            "negative_information_gain_molecules": sum(r["information_gain_nats"] < 0 for r in averaged),
            "evaluated_molecule_fits": len(records),
            "reported_outcome_tied_molecule_fits": sum(r["positive_tie_size"] > 1 for r in records),
            "candidate_count_breakdown": {str(k): summarize_records([r for r in averaged if r["n_candidates"] == k])
                                          for k in sorted({r["n_candidates"] for r in averaged})}}
    comparisons = []
    for left, right in plan["config"]["comparisons"]:
        for budget in ["50", "100", "full"]:
            keys = [f"{arm}:budget{budget}" for arm in [left, right]]
            if not all(key in results for key in keys):
                continue
            for metric in plan["config"]["metrics_with_paired_intervals"]:
                a, b = [{r["compound_id"]: r[metric] for r in results[key]["molecular_metrics"]} for key in keys]
                interval = paired_bootstrap(a, b, plan["compound_to_group"], seed=plan["config"]["bootstrap_seed"],
                                            repetitions=plan["config"]["bootstrap_repetitions"])
                comparisons.append({"left": left, "right": right, "budget": budget, "metric": metric,
                                    "direction": "positive favors left" if metric == "top1" else "negative favors left", **interval})
    return {"results": results, "paired_comparisons": comparisons, "descriptive_diagnostics": descriptive_diagnostics}


def run(root: Path, split_path: Path, output: Path) -> dict:
    plan = load(root, split_path)
    if output.exists():
        raise FileExistsError("Use a new output directory; follow-up runs are never overwritten or resumed")
    output.mkdir(parents=True)
    write_json(output / "split_manifest.json", plan, exclusive=True)
    bank = load_feature_bank(root, cohort=plan["config"]["cohort"], blocks=["electronic", "atom_jepa", "random_backbone", "morgan"])
    accumulated = defaultdict(list)
    completed = []
    diagnostics = []
    for arm in plan["config"]["arms"]:
        values, groups, morgan = materialize(bank, arm)
        for outer in plan["outer_folds"]:
            for training in outer["training_sets"]:
                if training["budget"] not in arm["budgets"]:
                    continue
                tag = f"{arm['id']}_fold{outer['fold']}_budget{training['budget']}_seed{training['subset_seed']}"
                try:
                    result, scorer = fit_setting(values, bank.queries, groups, morgan, arm, outer, training, plan["config"])
                except Exception as error:
                    write_json(output / "failure.json", {"task": tag, "error": repr(error), "complete": False}, exclusive=True)
                    raise
                model_path = output / f"{tag}.npz"
                np.savez_compressed(model_path, parameters=scorer.parameters, keep=scorer.transform.keep,
                                    mean=scorer.transform.mean, scale=scorer.transform.scale,
                                    hidden=np.array(scorer.hidden))
                result["model_sha256"] = file_sha256(model_path)
                write_json(output / f"{tag}.json", result, exclusive=True)
                accumulated[(arm["id"], str(training["budget"]))].extend(result["molecular_metrics"])
                diagnostics.append({"tag": tag, "arm": arm["id"], "budget": training["budget"],
                                    "grid_edges": result["grid_edges"], "final_fit": result["final_fit"]})
                completed.append(tag)
                print(f"Completed {tag}", flush=True)
    expected = workload(plan)
    if len(completed) != expected["final_models"]:
        raise ValueError("Incomplete frozen arm schedule")
    summary = {"proposal_id": plan["proposal_id"], "frozen_sha256": file_sha256(split_path),
               "workload": expected, "completed_final_models": len(completed),
               "fit_diagnostics": diagnostics, **summarize(accumulated, plan),
               "interpretation": "Adaptive exploratory follow-up; unadjusted conditional intervals, no causal descriptor or JEPA-objective attribution."}
    write_json(output / "summary.json", summary, exclusive=True)
    lines = ["# Exploratory follow-up results", "", f"Proposal: `{plan['proposal_id']}`.", "",
             "| Arm / budget | Molecules | Top-1 | NLL | Information gain | Negative gain molecules | Reported outcome ties / molecule-fits |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for key, result in summary["results"].items():
        m = result["metrics"]
        d = summary["descriptive_diagnostics"][key]
        lines.append(f"| {key} | {m['molecules']} | {m['top1']:.4f} | {m['nll']:.4f} | {m['information_gain_nats']:.4f} | {d['negative_information_gain_molecules']} | {d['reported_outcome_tied_molecule_fits']} / {d['evaluated_molecule_fits']} |")
    lines += ["", "## Paired conditional differences", "", "| Left / right | Budget | Metric | Difference | Unadjusted 95% interval |", "|---|---|---|---:|---|"]
    for c in summary["paired_comparisons"]:
        lines.append(f"| {c['left']} / {c['right']} | {c['budget']} | {c['metric']} | {c['mean_delta']:.4f} | [{c['ci95'][0]:.4f}, {c['ci95'][1]:.4f}] |")
    lines += ["", "Positive top-1 differences and negative NLL differences favor the left method. Intervals resample leakage components and are conditional on fitted models. Many comparisons are exploratory, unadjusted and evaluated on a previously inspected dataset.",
              "", "## Per-component grid endpoints", "", "| Arm / budget | Component | Minimum selections | Maximum selections | Final models |", "|---|---|---:|---:|---:|"]
    edges = defaultdict(list)
    for d in diagnostics:
        for component, flags in d["grid_edges"].items():
            edges[(d["arm"], str(d["budget"]), component)].append(flags)
    for (arm, budget, component), flags in edges.items():
        lines.append(f"| {arm} / {budget} | {component} | {sum(f['at_minimum'] for f in flags)} | {sum(f['at_maximum'] for f in flags)} | {len(flags)} |")
    lines += ["", "Endpoint selections are warnings, not proof that an outside-grid optimum or performance gain exists. The grid is not expanded using outer performance.",
              "", "Full candidate-count breakdowns, top-2/MRR, restart diagnostics, molecular predictions and feature masks are in JSON/NPZ artifacts.",
              "", "## Not done", "", "Independent review and external confirmation; no pretraining-objective or mechanistic conclusion follows.", ""]
    (output / "summary.md").write_text("\n".join(lines), encoding="utf-8", newline="\n")
    return summary
