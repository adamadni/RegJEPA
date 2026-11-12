import copy

import numpy as np
import pytest

from regjepa.evaluation import CandidateTable, choice_loss
from regjepa.experiment import fit_nested_split, penalty_candidates, penalty_vector
from regjepa.linear import LinearTransform, centered_features, fit_linear, linear_objective


def candidates(molecules=8):
    q, values = [], []
    for molecule in range(molecules):
        for site in range(2):
            positive = site == molecule % 2
            q.append({"query_id": f"m{molecule}:p{site}", "compound_id": f"m{molecule}",
                      "structure_key": f"s{molecule}", "outcome_group_index": site,
                      "n_candidate_outcomes": 2, "n_reported_outcomes": 1, "reported_positive": positive})
            values.append([1.0 if positive else -1.0, float(molecule), 7.0])
    return q, np.array(values)


def test_analytic_gradient_matches_finite_differences_with_variable_candidate_sets():
    q, _ = candidates(2)
    q[-1]["n_candidate_outcomes"] = q[-2]["n_candidate_outcomes"] = 3
    q.append(dict(q[-1], query_id="m1:p2", outcome_group_index=2, reported_positive=False))
    table = CandidateTable.from_queries(q)
    values = np.random.default_rng(7).normal(size=(5, 3))
    weights = np.array([.2, -.1, .3])
    penalties = np.array([.1, .2, .5])
    _, analytic = linear_objective(weights, values, table, penalties)
    numerical = []
    for i in range(3):
        step = np.zeros(3); step[i] = 1e-6
        numerical.append((linear_objective(weights+step, values, table, penalties)[0] -
                          linear_objective(weights-step, values, table, penalties)[0])/2e-6)
    np.testing.assert_allclose(analytic, numerical, atol=1e-8)


def test_training_transform_discards_only_unusable_linear_directions_and_is_not_refitted():
    q, values = candidates()
    table = CandidateTable.from_queries(q)
    transform = LinearTransform.fit(values, table, standardize=True)
    np.testing.assert_array_equal(transform.keep, [True, False, False])
    before = transform.scale.copy()
    shifted = values + np.repeat(np.arange(8)*1000, 2)[:, None]
    np.testing.assert_allclose(transform.apply(values, table), transform.apply(shifted, table))
    transform.apply(values*100, table)
    np.testing.assert_array_equal(before, transform.scale)


def test_transform_statistics_weight_molecules_equally_not_candidate_rows():
    q, _ = candidates(2)
    q[-1]["n_candidate_outcomes"] = q[-2]["n_candidate_outcomes"] = 4
    q.extend([dict(q[-1], query_id=f"m1:p{i}", outcome_group_index=i, reported_positive=False) for i in [2, 3]])
    values = np.array([[-1.], [1.], [-3.], [3.], [-3.], [3.]])
    transform = LinearTransform.fit(values, CandidateTable.from_queries(q), standardize=True)
    assert transform.scale[0] == pytest.approx(np.sqrt((1+9)/2))


def test_linear_fit_learns_shared_feature_relation_and_candidate_order_is_irrelevant():
    q, values = candidates()
    table = CandidateTable.from_queries(q)
    model = fit_linear(values, table, .1, standardize=True)
    scores = model.predict(values, table)
    assert all(scores[indices][positive] > scores[indices][1-positive] for indices, positive in zip(table.rows, table.positives))
    assert model.diagnostics["converged"]
    assert choice_loss(scores, table)[0] < np.log(2)
    order = np.random.default_rng(1).permutation(len(q))
    reordered = CandidateTable.from_queries([q[i] for i in order])
    other = fit_linear(values[order], reordered, .1, standardize=True)
    np.testing.assert_allclose(other.weights, model.weights, atol=1e-10)
    np.testing.assert_allclose(other.predict(values[order], reordered), scores[order])


def test_no_variable_training_features_produces_uniform_scores_without_division_by_zero():
    q, values = candidates()
    table = CandidateTable.from_queries(q)
    model = fit_linear(np.ones_like(values), table, .1, standardize=True)
    assert model.weights.size == 0
    np.testing.assert_array_equal(model.predict(np.ones_like(values), table), np.zeros(len(q)))


def test_optimizer_failure_is_not_silently_reported_as_a_valid_fit():
    q, values = candidates()
    with pytest.raises(RuntimeError, match="did not converge"):
        fit_linear(values, CandidateTable.from_queries(q), .01, standardize=True, max_iterations=1)


def test_concat_penalties_are_applied_to_the_correct_components():
    vector = penalty_vector({"jepa": 10., "qm": .01}, 774)
    np.testing.assert_array_equal(vector[:768], np.full(768, 10.))
    np.testing.assert_array_equal(vector[768:], np.full(6, .01))
    config = {"l2_grid": [.1, 1.], "combined_penalties": "separate_jepa_qm"}
    assert len(penalty_candidates("jepa_plus_qm", config)) == 4
    assert len(penalty_candidates("atom_jepa", config)) == 2


def test_nested_selection_and_fitted_weights_do_not_see_outer_test_labels_or_features():
    q, values = candidates(10)
    train = [f"m{i}" for i in range(8)]
    outer = {"fold": 0, "train_compounds": train, "test_compounds": ["m8", "m9"]}
    training = {"budget": "full", "subset_seed": None, "actual_labels": 8, "train_compounds": train,
                "inner_folds": [{"train_compounds": train[:4], "test_compounds": train[4:]},
                                {"train_compounds": train[4:], "test_compounds": train[:4]}]}
    config = {"l2_grid": [.01, 1.], "combined_penalties": "shared", "optimizer": {"max_iterations": 1000}}
    first, model = fit_nested_split(values, q, outer, training, "electronic", config)
    changed = copy.deepcopy(q)
    for row in changed:
        if row["compound_id"] in outer["test_compounds"]:
            row["reported_positive"] = not row["reported_positive"]
    modified = values.copy(); modified[-4:] *= 500
    second, other = fit_nested_split(modified, changed, outer, training, "electronic", config)
    assert first["selected_penalties"] == second["selected_penalties"]
    assert first["search"] == second["search"]
    np.testing.assert_array_equal(model.weights, other.weights)
    np.testing.assert_array_equal(model.transform.scale, other.transform.scale)
    assert {r["compound_id"] for r in first["predictions"]} == {"m8", "m9"}
