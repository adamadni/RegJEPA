import json

import numpy as np
import pytest

from regjepa.experiment import run_experiment
from regjepa.feature_io import FeatureBank
from regjepa.splits import build_split_plan


def test_full_runner_preserves_held_out_coverage_and_paired_artifacts_on_synthetic_data(monkeypatch, tmp_path):
    queries, values = [], []
    names = [f"m{i:02}" for i in range(12)]
    for i, compound in enumerate(names):
        for site in range(2):
            positive = site == i % 2
            queries.append({"query_id": f"{compound}:p{site}", "compound_id": compound,
                            "structure_key": f"s{i}", "outcome_group_index": site,
                            "reported_positive": positive, "n_reported_outcomes": 1, "n_candidate_outcomes": 2})
            values.append([1. if positive else -1., float(i), 7.])
    values = np.array(values)
    original = values.copy()
    config = {"cohort": "synthetic_only", "blocks": ["atom_jepa", "random_backbone"],
              "outer_folds": 3, "inner_folds": 2, "split_seed": 42, "subset_seeds": [17, 29],
              "label_budgets": [4], "l2_grid": [.1, 1.], "combined_penalties": "shared",
              "optimizer": {"max_iterations": 1000}, "bootstrap_seed": 42, "bootstrap_repetitions": 30}
    groups = {compound: f"g{i}" for i, compound in enumerate(names)}
    plan = {"status": "frozen", "proposal_id": "synthetic", "config": config, "compound_to_group": groups,
            "outer_folds": build_split_plan(names, groups, config)}
    split_path = tmp_path / "synthetic-splits.json"
    split_path.write_text(json.dumps(plan))
    monkeypatch.setattr("regjepa.experiment.load_split_plan", lambda root, path: plan)
    monkeypatch.setattr("regjepa.experiment.load_feature_bank", lambda root, cohort, blocks:
                        FeatureBank({blocks[0]: values}, queries, {}))
    output = tmp_path / "synthetic-run"
    summary = run_experiment(tmp_path, split_path, output)
    assert summary["completed_final_fits"] == 18
    assert len(list(output.glob("*.npz"))) == 18
    assert (output / "summary.md").is_file()
    assert json.loads((output / "summary.json").read_text())["proposal_id"] == "synthetic"
    for setting in summary["results"].values():
        assert setting["metrics"]["molecules"] == 12
        assert {r["compound_id"] for r in setting["molecular_metrics"]} == set(names)
    assert all(r["repetitions"] == 2 for r in summary["results"]["atom_jepa:budget4"]["molecular_metrics"])
    assert all(r["repetitions"] == 1 for r in summary["results"]["atom_jepa:budgetfull"]["molecular_metrics"])
    assert len(summary["paired_top1_comparisons"]) == 2
    np.testing.assert_array_equal(values, original)
    with pytest.raises(FileExistsError):
        run_experiment(tmp_path, split_path, output)
