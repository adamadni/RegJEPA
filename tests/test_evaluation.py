import copy
import json
from pathlib import Path

import numpy as np
import pytest

from regjepa.evaluation import CandidateTable, choice_loss, paired_bootstrap, ranking_records, summarize_records
from regjepa.feature_io import load_feature_bank
from regjepa.linear import LinearScorer, LinearTransform, centered_features

ROOT = Path(__file__).resolve().parents[1]


def queries(counts=(2, 3), positives=(0, 1)):
    result = []
    for molecule, (count, positive) in enumerate(zip(counts, positives)):
        for index in range(count):
            result.append({"query_id": f"m{molecule}:p{index}", "compound_id": f"m{molecule}",
                           "structure_key": f"s{molecule}", "outcome_group_index": index,
                           "n_candidate_outcomes": count, "n_reported_outcomes": 1,
                           "reported_positive": index == positive})
    return result


def test_exact_top_tie_has_expected_reciprocal_rank_not_reciprocal_average_rank():
    table = CandidateTable.from_queries(queries((2,), (0,)))
    got = ranking_records(np.array([1., 1.]), table)[0]
    assert got["top1"] == 0.5
    assert got["top2"] == 1.0
    assert got["reciprocal_rank"] == 0.75


def test_middle_tie_has_fractional_top2_and_expected_reciprocal_rank():
    table = CandidateTable.from_queries(queries((3,), (1,)))
    got = ranking_records(np.array([3., 1., 1.]), table)[0]
    assert got["top1"] == 0.0
    assert got["top2"] == 0.5
    assert got["reciprocal_rank"] == pytest.approx((1/2 + 1/3)/2)


def test_loss_is_stable_at_large_logits_and_unchanged_by_molecule_offsets():
    table = CandidateTable.from_queries(queries())
    scores = np.array([1000., 999., -1000., -999., -1001.])
    first = choice_loss(scores, table)
    shifted = scores + np.array([700., 700., -500., -500., -500.])
    second = choice_loss(shifted, table)
    assert first[0] == pytest.approx(second[0])
    np.testing.assert_allclose(first[1], second[1])
    for indices in table.rows:
        assert first[1][indices].sum() == pytest.approx(0, abs=1e-15)


def test_each_molecule_contributes_one_loss_independent_of_candidate_count():
    table = CandidateTable.from_queries(queries((2, 5), (0, 0)))
    probabilities = np.array([0.8, 0.2, 0.6, 0.1, 0.1, 0.1, 0.1])
    loss, _ = choice_loss(np.log(probabilities), table)
    assert loss == pytest.approx((-np.log(.8)-np.log(.6))/2)


def test_metrics_and_loss_are_invariant_to_candidate_row_order():
    q = queries()
    scores = np.array([0., 0., 1., 1., 0.])
    order = np.array([4, 1, 3, 0, 2])
    first = CandidateTable.from_queries(q)
    second = CandidateTable.from_queries([q[i] for i in order])
    assert ranking_records(scores, first) == ranking_records(scores[order], second)
    assert choice_loss(scores, first)[0] == pytest.approx(choice_loss(scores[order], second)[0])


@pytest.mark.parametrize("corruption", ["missing_candidate", "multiple_positive", "no_positive", "duplicate_query", "no_choice"])
def test_invalid_or_partial_candidate_sets_are_rejected(corruption):
    q = queries()
    if corruption == "missing_candidate":
        q.pop()
    elif corruption == "multiple_positive":
        q[1]["reported_positive"] = True
    elif corruption == "no_positive":
        q[0]["reported_positive"] = False
    elif corruption == "duplicate_query":
        q[1]["query_id"] = q[0]["query_id"]
    else:
        q = queries((1,), (0,))
    with pytest.raises(ValueError):
        CandidateTable.from_queries(q)


def test_constant_real_bank_predictions_recover_analytic_baseline_and_not_source_order():
    bank = load_feature_bank(ROOT, cohort="conservative_bromination", blocks=["morgan"])
    table = CandidateTable.from_queries(bank.queries)
    got = summarize_records(ranking_records(np.zeros(table.n_rows), table))
    assert got["top1"] == pytest.approx(sum(1/len(rows) for rows in table.rows)/299)
    assert got["top1"] == pytest.approx(0.28314752879970273)
    assert sum(positive == 0 for positive in table.positives) == 92
    assert got["top1"] != pytest.approx(92/299)
    duplicate_molecules = []
    for compound, indices in zip(table.compound_ids, table.rows):
        values = bank.X["morgan"][indices]
        if len(np.unique(values, axis=0)) < len(values):
            duplicate_molecules.append(compound)
    assert "comp40" in duplicate_molecules
    assert len(duplicate_molecules) == 6


def test_identical_real_morgan_rows_remain_exact_score_ties_with_this_platform_blas():
    bank = load_feature_bank(ROOT, cohort="conservative_bromination", blocks=["morgan"])
    table = CandidateTable.from_queries(bank.queries)
    values = centered_features(bank.X["morgan"], table)
    pairs = [(a, b) for indices in table.rows for j, a in enumerate(indices) for b in indices[j+1:]
             if np.array_equal(values[a], values[b])]
    assert pairs
    # Fixed random coefficients test numerical equality, not predictive efficacy;
    # no coefficient is fitted to these real-data labels.
    for seed in [7, 19, 43]:
        rng = np.random.default_rng(seed)
        coefficients = rng.normal(size=values.shape[1])
        transform = LinearTransform(np.ones(values.shape[1], dtype=bool), np.ones(values.shape[1]), False)
        scorer = LinearScorer(transform, coefficients, {})
        for order in [np.arange(len(values)), rng.permutation(len(values))]:
            permuted = CandidateTable.from_queries([bank.queries[i] for i in order])
            scores = scorer.predict(bank.X["morgan"][order], permuted)
            inverse = np.argsort(order)
            assert all(scores[inverse[a]] == scores[inverse[b]] for a, b in pairs)


def test_bootstrap_resamples_whole_leakage_components_and_preserves_pairing():
    left = {"a": .5, "b": .8, "c": .4}
    right = {key: value-.2 for key, value in left.items()}
    result = paired_bootstrap(left, right, {"a": "g1", "b": "g1", "c": "g2"}, repetitions=100)
    assert result["molecules"] == 3
    assert result["resampled_leakage_groups"] == 2
    np.testing.assert_allclose(result["ci95"], [.2, .2])
    with pytest.raises(ValueError):
        paired_bootstrap(left, {"a": .1}, {"a": "g1"}, repetitions=10)
