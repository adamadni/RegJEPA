"""Label-blind grouped folds and nested label budgets; cohorts are not test splits."""
from __future__ import annotations

from collections import defaultdict

import numpy as np


def grouped_folds(compounds: list[str], leakage_groups: dict[str, str], *, folds: int, seed: int) -> list[dict]:
    members: dict[str, list[str]] = defaultdict(list)
    if len(compounds) != len(set(compounds)) or not set(compounds) <= set(leakage_groups):
        raise ValueError("Unique compounds with known leakage groups are required")
    for compound in sorted(compounds):
        members[leakage_groups[compound]].append(compound)
    if folds < 2 or len(members) < folds:
        raise ValueError(f"Cannot make {folds} grouped folds from {len(members)} leakage groups")
    rng = np.random.default_rng(seed)
    tie_breaks = dict(zip(sorted(members), rng.random(len(members))))
    ordered = sorted(members, key=lambda g: (-len(members[g]), tie_breaks[g], g))
    bins: list[list[str]] = [[] for _ in range(folds)]
    for group in ordered:
        choices = [i for i in range(folds) if len(bins[i]) == min(map(len, bins))]
        target = int(rng.choice(choices))
        bins[target].extend(members[group])
    result = []
    universe = set(compounds)
    for i, held_out in enumerate(bins):
        test = sorted(held_out)
        train = sorted(universe - set(test))
        if {leakage_groups[x] for x in train} & {leakage_groups[x] for x in test}:
            raise ValueError("Leakage group crosses fold boundary")
        result.append({"fold": i, "train_compounds": train, "test_compounds": test})
    return result


def nested_budget_sets(compounds: list[str], leakage_groups: dict[str, str], budgets: list[int],
                       *, seed: int) -> list[dict]:
    """Sample whole components; actual counts may be below the requested budget."""
    if any(b <= 0 for b in budgets) or len(set(budgets)) != len(budgets):
        raise ValueError("Budgets must be distinct positive integers")
    members: dict[str, list[str]] = defaultdict(list)
    for compound in sorted(compounds):
        members[leakage_groups[compound]].append(compound)
    ordered = list(np.random.default_rng(seed).permutation(sorted(members)))
    selected: set[str] = set()
    result = []
    for budget in sorted(budgets):
        if budget >= len(compounds):
            raise ValueError("Numeric label budgets must be smaller than the full training fold")
        for group in ordered:
            if group not in selected and sum(len(members[g]) for g in selected) + len(members[group]) <= budget:
                selected.add(group)
        ids = sorted(compound for group in selected for compound in members[group])
        if not ids:
            raise ValueError(f"No whole leakage group fits label budget {budget}")
        result.append({"budget": budget, "actual_labels": len(ids), "train_compounds": ids})
    return result


def build_split_plan(compounds: list[str], leakage_groups: dict[str, str], config: dict) -> list[dict]:
    outer = grouped_folds(compounds, leakage_groups, folds=config["outer_folds"], seed=config["split_seed"])
    for fold in outer:
        training_sets = []
        for subset_seed in config["subset_seeds"]:
            subsets = nested_budget_sets(fold["train_compounds"], leakage_groups, config["label_budgets"], seed=subset_seed)
            for subset in subsets:
                subset["subset_seed"] = subset_seed
                training_sets.append(subset)
        training_sets.append({"budget": "full", "actual_labels": len(fold["train_compounds"]),
                              "subset_seed": None, "train_compounds": fold["train_compounds"]})
        for training in training_sets:
            training["inner_folds"] = grouped_folds(training["train_compounds"], leakage_groups,
                                                   folds=config["inner_folds"], seed=config["split_seed"] + fold["fold"] + 1)
        fold["training_sets"] = training_sets
    validate_split_plan(outer, compounds, leakage_groups, config)
    return outer


def validate_split_plan(plan: list[dict], compounds: list[str], leakage_groups: dict[str, str], config: dict) -> None:
    universe = set(compounds)
    if len(universe) != len(compounds) or not universe <= set(leakage_groups):
        raise ValueError("Expected compounds must be unique and have known leakage groups")
    _check_fold_ids(plan, config["outer_folds"])
    expected_settings = {(budget, seed) for budget in config["label_budgets"] for seed in config["subset_seeds"]}
    expected_settings.add(("full", None))
    tests = []
    for fold in plan:
        train, test = _unique_ids(fold["train_compounds"]), _unique_ids(fold["test_compounds"])
        _check_partition(train, test, universe, leakage_groups)
        tests.extend(test)
        training_sets = fold.get("training_sets")
        if not isinstance(training_sets, list):
            raise ValueError("Configured training sets are missing")
        settings = [(item["budget"], item["subset_seed"]) for item in training_sets]
        if len(settings) != len(expected_settings) or set(settings) != expected_settings:
            raise ValueError("Configured training sets are missing, duplicated or unexpected")
        by_seed: dict[int, list[dict]] = defaultdict(list)
        for training in training_sets:
            selected = _unique_ids(training["train_compounds"])
            if not selected <= train or len(selected) != training["actual_labels"]:
                raise ValueError("Training set escapes outer training fold or has incorrect label count")
            if training["budget"] == "full":
                if selected != train:
                    raise ValueError("Full budget omits training records")
            else:
                if len(selected) > training["budget"]:
                    raise ValueError("Training label budget exceeded")
                by_seed[training["subset_seed"]].append(training)
            selected_groups = {leakage_groups[x] for x in selected}
            if any(leakage_groups[x] in selected_groups for x in train - selected):
                raise ValueError("A training subset splits a leakage group")
            inner_tests = []
            _check_fold_ids(training["inner_folds"], config["inner_folds"])
            for inner in training["inner_folds"]:
                a, b = _unique_ids(inner["train_compounds"]), _unique_ids(inner["test_compounds"])
                _check_partition(a, b, selected, leakage_groups)
                inner_tests.extend(b)
            if len(inner_tests) != len(selected) or set(inner_tests) != selected:
                raise ValueError("Inner validation coverage is incomplete or repeated")
        for subsets in by_seed.values():
            ordered = sorted(subsets, key=lambda x: x["budget"])
            if any(not set(a["train_compounds"]) <= set(b["train_compounds"]) for a, b in zip(ordered, ordered[1:])):
                raise ValueError("Learning-curve subsets are not nested")
    if len(tests) != len(universe) or set(tests) != universe:
        raise ValueError("Outer held-out coverage is incomplete or repeated")


def _check_fold_ids(folds: list[dict], expected: int) -> None:
    if len(folds) != expected or {fold["fold"] for fold in folds} != set(range(expected)):
        raise ValueError("Configured fold counts or fold identifiers do not match")


def _unique_ids(values: list[str]) -> set[str]:
    result = set(values)
    if len(result) != len(values):
        raise ValueError("A fold or training set contains duplicate record IDs")
    return result


def _check_partition(train: set[str], test: set[str], universe: set[str], groups: dict[str, str]) -> None:
    if not train or not test or train & test or train | test != universe:
        raise ValueError("Invalid train/test partition")
    if {groups[x] for x in train} & {groups[x] for x in test}:
        raise ValueError("Leakage group crosses a fold boundary")
