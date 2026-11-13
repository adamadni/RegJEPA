import copy
import json
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from regjepa.evaluation import CandidateTable, choice_loss, ranking_records
from regjepa.feature_io import FeatureBank
from regjepa.splits import build_split_plan
from regjepa_followup.experiment import fit_setting, run, vector, workload
from regjepa_followup.features import materialize, validate_arm
from regjepa_followup.models import Scorer, Transform, fit, nonlinear_objective, nonlinear_scores
from regjepa_followup.protocol import freeze, load, validate_config


def chemistry_context(molecules=12):
    queries, values = [], []
    for i in range(molecules):
        context = 1.0 if i % 2 else -1.0
        for site, local in enumerate([1.0, -1.0]):
            queries.append({"query_id": f"m{i:02}:o{site}", "compound_id": f"m{i:02}",
                            "structure_key": f"s{i}", "outcome_group_index": site,
                            "n_candidate_outcomes": 2, "n_reported_outcomes": 1,
                            "reported_positive": bool(local * context > 0)})
            values.append([local, context, 7.0])
    return queries, np.array(values)


def optimizer():
    return {"max_iterations": 5000, "gradient_tolerance": 1e-7, "accepted_gradient": 1e-5,
            "mlp_hidden": 16, "mlp_restart_seeds": [11, 23]}


def test_nonlinear_gradient_matches_finite_differences():
    q, x = chemistry_context(4)
    table = CandidateTable.from_queries(q)
    rng = np.random.default_rng(8)
    hidden = 3
    p = rng.normal(0, .2, x.shape[1] * hidden + 2 * hidden)
    loss, derivative = nonlinear_objective(p, x, table, .03, hidden)
    numerical = np.zeros_like(p)
    for index in range(len(p)):
        delta = np.zeros_like(p); delta[index] = 1e-6
        numerical[index] = (nonlinear_objective(p + delta, x, table, .03, hidden)[0] -
                            nonlinear_objective(p - delta, x, table, .03, hidden)[0]) / 2e-6
    assert np.isfinite(loss)
    np.testing.assert_allclose(derivative, numerical, atol=1e-8, rtol=1e-6)


def test_global_transform_preserves_molecular_context_and_is_train_only():
    q, x = chemistry_context(4)
    table = CandidateTable.from_queries(q)
    transform = Transform.fit(x, table, "mlp", np.zeros(3, dtype=bool))
    np.testing.assert_array_equal(transform.keep, [True, True, False])
    for rows in table.rows:
        assert np.ptp(transform.apply(x, table)[rows, 1]) == 0
    before = (transform.mean.copy(), transform.scale.copy())
    transform.apply(x * 100 + 42, table)
    np.testing.assert_array_equal(transform.mean, before[0])
    np.testing.assert_array_equal(transform.scale, before[1])
    altered = [dict(r, reported_positive=not r["reported_positive"]) for r in q]
    other = Transform.fit(x, CandidateTable.from_queries(altered), "mlp", np.zeros(3, dtype=bool))
    np.testing.assert_array_equal(other.mean, transform.mean)
    np.testing.assert_array_equal(other.scale, transform.scale)


def test_scaling_weights_molecules_not_rows_and_morgan_stays_unscaled():
    q, _ = chemistry_context(2)
    q[-1]["n_candidate_outcomes"] = q[-2]["n_candidate_outcomes"] = 4
    q.extend([dict(q[-1], query_id=f"m01:o{i}", outcome_group_index=i, reported_positive=False) for i in [2, 3]])
    x = np.array([[-1., -1.], [1., 1.], [-3., -3.], [3., 3.], [-3., -3.], [3., 3.]])
    transform = Transform.fit(x, CandidateTable.from_queries(q), "linear", np.array([False, True]))
    np.testing.assert_allclose(transform.scale, [np.sqrt(5.), 1.])


def test_nonlinear_head_can_learn_context_interaction_that_linear_cannot():
    q, x = chemistry_context()
    table = CandidateTable.from_queries(q)
    mask = np.zeros(x.shape[1], dtype=bool)
    nonlinear = fit(x, table, mask, "mlp", np.full(3, .01), optimizer())
    linear = fit(x, table, mask, "linear", np.full(3, .01), optimizer())
    assert choice_loss(nonlinear.predict(x, table), table)[0] < .3
    assert choice_loss(linear.predict(x, table), table)[0] == pytest.approx(np.log(2), abs=1e-7)
    attempts = nonlinear.diagnostics["attempts"]
    assert len(attempts) == 2 and all(a["converged"] for a in attempts)
    assert attempts[nonlinear.diagnostics["selected_restart"]]["objective"] == min(a["objective"] for a in attempts)


def test_fixed_nonlinear_scores_preserve_identical_rows_and_permutation():
    rng = np.random.default_rng(3)
    x = rng.normal(size=(12, 8)); x[8] = x[0]; x[11] = x[0]
    p = rng.normal(size=8 * 4 + 8)
    scores = nonlinear_scores(x, p, 4)[0]
    assert scores[0] == scores[8] == scores[11]
    order = rng.permutation(len(x))
    np.testing.assert_array_equal(nonlinear_scores(x[order], p, 4)[0], scores[order])


def test_feature_registry_rejects_metadata_and_respects_literal_qm_columns():
    q, x = chemistry_context(2)
    arm = {"id": "test", "head": "linear", "budgets": ["full"],
           "parts": [{"block": "electronic", "columns": [1, 3], "penalty_group": "qm"}]}
    bank = FeatureBank({"electronic": np.arange(len(q) * 6).reshape(len(q), 6)}, q, {})
    values, groups, morgan = materialize(bank, arm)
    np.testing.assert_array_equal(values, bank.X["electronic"][:, [1, 3]])
    assert groups == ["qm", "qm"] and not morgan.any()
    for block in ["reported_positive", "n_reported_outcomes", "product", "yield"]:
        bad = copy.deepcopy(arm); bad["parts"][0]["block"] = block
        with pytest.raises(ValueError):
            validate_arm(bad)
    bad = copy.deepcopy(arm); bad["parts"][0]["columns"] = [1, 1]
    with pytest.raises(ValueError):
        validate_arm(bad)


def test_full_configuration_preserves_descriptor_masks_and_matched_controls():
    config = json.loads((Path(__file__).parents[1] / "configs/followup.json").read_text())
    validate_config(config)
    arms = {a["id"]: a for a in config["arms"]}
    assert len(arms) == 37
    assert arms["qm_fminus_linear"]["parts"][0]["columns"] == [3]
    assert arms["qm_fukui_pair_linear"]["parts"][0]["columns"] == [2, 3]
    assert arms["qm_no_fukui_linear"]["parts"][0]["columns"] == [0, 1, 5]
    assert arms["jepa_scalar_linear"]["parts"][0]["columns"] == list(range(256))
    assert arms["random_scalar_linear"]["parts"][0]["columns"] == list(range(256))
    assert arms["random_qm_linear"]["penalty_mode"] == "separate"
    assert ["jepa_qm_linear", "random_qm_linear"] in config["comparisons"]
    assert ["jepa_qm_mlp", "qm_all_mlp"] in config["comparisons"]
    np.testing.assert_array_equal(vector({"representation": 10., "qm": .001}, ["representation", "qm", "qm"]), [10., .001, .001])


def test_frozen_binding_rejects_status_edits_and_changed_proposals(monkeypatch, tmp_path):
    proposal = {"schema_version": 1, "status": "proposal", "proposal_id": "synthetic", "config": {}}
    monkeypatch.setattr("regjepa_followup.protocol.prepare", lambda root, config: copy.deepcopy(proposal))
    proposal_path = tmp_path / "proposal.json"
    proposal_path.write_text(json.dumps(proposal))
    forged = dict(proposal, status="frozen")
    frozen_path = tmp_path / "frozen.json"
    frozen_path.write_text(json.dumps(forged))
    with pytest.raises(ValueError):
        load(tmp_path, frozen_path)
    frozen_path.unlink()
    freeze(tmp_path, proposal_path)
    assert load(tmp_path, frozen_path)["status"] == "frozen"
    proposal_path.write_text(json.dumps(dict(proposal, config={"changed": True})))
    with pytest.raises(ValueError):
        load(tmp_path, frozen_path)


def test_runner_and_saved_models_cover_both_metrics_without_refitting(monkeypatch, tmp_path):
    q, x = chemistry_context(12)
    config = {"cohort": "synthetic", "outer_folds": 3, "inner_folds": 2, "split_seed": 42,
              "subset_seeds": [17], "label_budgets": [4], "l2_grid": [.1, 1.], "optimizer": optimizer(),
              "bootstrap_seed": 42, "bootstrap_repetitions": 30, "metrics_with_paired_intervals": ["top1", "nll"]}
    part = [{"block": "electronic", "columns": [0, 1, 2], "penalty_group": "qm"}]
    config["arms"] = [{"id": head, "head": head, "parts": part, "budgets": ["full"], "penalty_mode": "shared"} for head in ["linear", "mlp"]]
    config["comparisons"] = [["mlp", "linear"]]
    names = sorted({r["compound_id"] for r in q})
    groups = {c: c for c in names}
    plan = {"proposal_id": "synthetic", "status": "frozen", "config": config,
            "outer_folds": build_split_plan(names, groups, config), "compound_to_group": groups}
    split = tmp_path / "frozen.json"; split.write_text(json.dumps(plan))
    monkeypatch.setattr("regjepa_followup.experiment.load", lambda root, path: plan)
    monkeypatch.setattr("regjepa_followup.experiment.load_feature_bank", lambda *args, **kwargs: FeatureBank({"electronic": x}, q, {}))
    out = tmp_path / "run"
    result = run(tmp_path, split, out)
    assert result["completed_final_models"] == 6
    assert {c["metric"] for c in result["paired_comparisons"]} == {"top1", "nll"}
    assert workload(plan)["final_optimizer_attempts"] == 9
    spec = importlib.util.spec_from_file_location("followup_validator", Path(__file__).parents[1] / "scripts/validate_followup.py")
    validator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(validator)
    monkeypatch.setattr(validator, "load", lambda root, path: plan)
    monkeypatch.setattr(validator, "load_feature_bank", lambda *args, **kwargs: FeatureBank({"electronic": x}, q, {}))
    audit = validator.validate(tmp_path, split, out)
    assert audit["final_models_replayed"] == 6 and audit["models_refitted"] == 0
    for path in out.glob("*_fold*.json"):
        fitted = json.loads(path.read_text())
        outer = next(o for o in plan["outer_folds"] if o["fold"] == fitted["outer_fold"])
        row_ids = [i for i, r in enumerate(q) if r["compound_id"] in outer["test_compounds"]]
        table = CandidateTable.from_queries([q[i] for i in row_ids])
        with np.load(path.with_suffix(".npz"), allow_pickle=False) as saved:
            model = Scorer(Transform(fitted["head"], saved["keep"], saved["mean"], saved["scale"]), saved["parameters"], int(saved["hidden"]), fitted["final_fit"])
            scores = model.predict(x[row_ids], table)
        np.testing.assert_array_equal(scores, [p["score"] for p in fitted["predictions"]])
        assert ranking_records(scores, table) == fitted["molecular_metrics"]
    with pytest.raises(FileExistsError):
        run(tmp_path, split, out)
    bad = copy.deepcopy(plan["outer_folds"][0]["training_sets"][-1])
    bad["inner_folds"][0]["test_compounds"].append("outside")
    with pytest.raises(ValueError, match="Inner held-out"):
        fit_setting(x, q, ["qm"] * 3, np.zeros(3, dtype=bool), config["arms"][0], plan["outer_folds"][0], bad, config)
