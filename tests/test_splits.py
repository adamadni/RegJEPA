import copy

import pytest

from regjepa.splits import build_split_plan, grouped_folds, nested_budget_sets, validate_split_plan


def identity():
    names = [f"m{i:02}" for i in range(32)]
    return names, {name: f"g{i//2}" for i, name in enumerate(names)}


def test_outer_and_inner_folds_keep_complete_groups_and_nested_budgets():
    names, groups = identity()
    config = {"outer_folds": 4, "inner_folds": 3, "split_seed": 42, "subset_seeds": [17, 29], "label_budgets": [8, 16]}
    plan = build_split_plan(names, groups, config)
    validate_split_plan(plan, names, groups, config)
    assert plan == build_split_plan(list(reversed(names)), groups, config)
    assert sorted(x for f in plan for x in f["test_compounds"]) == names
    for fold in plan:
        assert len(fold["training_sets"]) == 5
        for training in fold["training_sets"]:
            for inner in training["inner_folds"]:
                assert not set(inner["test_compounds"]) & set(fold["test_compounds"])


def test_label_budget_never_splits_a_component_to_reach_requested_count():
    names, groups = identity()
    budgets = nested_budget_sets(names, groups, [5, 9], seed=7)
    assert [x["actual_labels"] for x in budgets] == [4, 8]
    assert set(budgets[0]["train_compounds"]) <= set(budgets[1]["train_compounds"])


def test_impossible_grouped_fold_request_is_not_silently_relaxed():
    with pytest.raises(ValueError, match="Cannot make"):
        grouped_folds(["a", "b"], {"a": "same", "b": "same"}, folds=2, seed=42)


def test_corrupt_split_that_separates_same_structure_is_rejected():
    names, groups = identity()
    config = {"outer_folds": 4, "inner_folds": 2, "split_seed": 42, "subset_seeds": [17], "label_budgets": [8]}
    plan = build_split_plan(names, groups, config)
    bad = copy.deepcopy(plan)
    held = bad[0]["test_compounds"].pop()
    bad[0]["train_compounds"].append(held)
    with pytest.raises(ValueError, match="Leakage group"):
        validate_split_plan(bad, names, groups, config)
