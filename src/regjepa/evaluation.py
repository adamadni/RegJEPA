"""Rank distinct outcome rows; scores do not represent measured product fractions."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CandidateTable:
    """Candidate sets use record IDs; split independence always uses leakage groups."""

    compound_ids: tuple[str, ...]
    structure_keys: tuple[str, ...]
    rows: tuple[np.ndarray, ...]
    positives: tuple[int, ...]
    n_rows: int

    @classmethod
    def from_queries(cls, queries: list[dict]) -> CandidateTable:
        if not queries:
            raise ValueError("An empty ranking dataset is not valid")
        if len({q["query_id"] for q in queries}) != len(queries):
            raise ValueError("Duplicate query IDs")
        by_record: dict[str, list[int]] = defaultdict(list)
        for i, query in enumerate(queries):
            if not isinstance(query["reported_positive"], (bool, np.bool_)):
                raise ValueError("Recorded-site labels must be boolean, not truthy text or numbers")
            by_record[query["compound_id"]].append(i)
        ids, structures, groups, positives = [], [], [], []
        for compound in sorted(by_record):
            indices = by_record[compound]
            records = [queries[i] for i in indices]
            keys = {r["structure_key"] for r in records}
            declared = {r["n_candidate_outcomes"] for r in records}
            observed = [i for i, r in enumerate(records) if r["reported_positive"]]
            if len(keys) != 1 or declared != {len(indices)}:
                raise ValueError(f"Incomplete or inconsistent candidate set: {compound}")
            if len(indices) < 2 or len(observed) != 1 or any(r["n_reported_outcomes"] != 1 for r in records):
                raise ValueError(f"Single-outcome ranking requires >=2 candidates and one positive: {compound}")
            if len({r["outcome_group_index"] for r in records}) != len(indices):
                raise ValueError(f"Repeated outcome group: {compound}")
            ids.append(compound)
            structures.append(next(iter(keys)))
            groups.append(np.array(indices, dtype=np.int64))
            positives.append(observed[0])
        return cls(tuple(ids), tuple(structures), tuple(groups), tuple(positives), len(queries))

    def check_scores(self, scores: np.ndarray) -> np.ndarray:
        values = np.asarray(scores, dtype=np.float64)
        if values.shape != (self.n_rows,) or not np.isfinite(values).all():
            raise ValueError("Scores must be finite and aligned to every candidate row")
        return values


def choice_loss(scores: np.ndarray, table: CandidateTable) -> tuple[float, np.ndarray]:
    """Mean conditional cross-entropy and gradient, equally weighted per molecule."""
    values = table.check_scores(scores)
    derivative = np.zeros_like(values)
    losses = []
    for indices, positive in zip(table.rows, table.positives):
        logits = values[indices]
        shifted = logits - logits.max()
        exp = np.exp(shifted)
        probabilities = exp / exp.sum()
        losses.append(float(np.log(exp.sum()) - shifted[positive]))
        probabilities[positive] -= 1.0
        derivative[indices] = probabilities / len(table.rows)
    return float(np.mean(losses)), derivative


def ranking_records(scores: np.ndarray, table: CandidateTable) -> list[dict]:
    """Expected metrics under uniform ordering of exact ties; no atom-score pooling."""
    values = table.check_scores(scores)
    result = []
    for compound, structure, indices, positive in zip(table.compound_ids, table.structure_keys, table.rows, table.positives):
        logits = values[indices]
        target = logits[positive]
        above = int((logits > target).sum())
        tied = int((logits == target).sum())
        shifted = logits - logits.max()
        nll = float(np.log(np.exp(shifted).sum()) - shifted[positive])
        ranks = above + np.arange(1, tied + 1, dtype=np.float64)
        count = len(indices)
        result.append({
            "compound_id": compound, "structure_key": structure, "n_candidates": count,
            "top1": 1.0 / tied if above == 0 else 0.0,
            "top2": float(np.clip(2 - above, 0, tied) / tied),
            "reciprocal_rank": float(np.mean(1.0 / ranks)),
            "nll": nll, "uniform_top1": 1.0 / count,
            "uniform_nll": float(np.log(count)), "information_gain_nats": float(np.log(count) - nll),
            "positive_tie_size": tied, "n_above_positive": above,
        })
    return result


def summarize_records(records: list[dict]) -> dict:
    if not records or len({r["compound_id"] for r in records}) != len(records):
        raise ValueError("Summaries require one complete prediction per source record")
    names = ["top1", "top2", "reciprocal_rank", "nll", "uniform_top1", "uniform_nll", "information_gain_nats"]
    summary = {name: float(np.mean([r[name] for r in records])) for name in names}
    nontrivial = [r for r in records if r["n_candidates"] > 2]
    summary.update({"molecules": len(records), "structures": len({r["structure_key"] for r in records}),
                    "top2_nontrivial_molecules": len(nontrivial),
                    "top2_nontrivial": float(np.mean([r["top2"] for r in nontrivial])) if nontrivial else None})
    return summary


def paired_bootstrap(left: dict[str, float], right: dict[str, float], leakage_groups: dict[str, str],
                     *, seed: int = 42, repetitions: int = 10000) -> dict:
    """Paired cluster bootstrap of molecular metrics, conditional on fitted models."""
    if not left or set(left) != set(right) or not set(left) <= set(leakage_groups):
        raise ValueError("Paired comparisons require identical molecules and known leakage groups")
    if repetitions < 2:
        raise ValueError("At least two bootstrap repetitions are required")
    grouped: dict[str, list[float]] = defaultdict(list)
    for compound in sorted(left):
        delta = float(left[compound]) - float(right[compound])
        if not np.isfinite(delta):
            raise ValueError("Nonfinite paired metric")
        grouped[leakage_groups[compound]].append(delta)
    sums = np.array([sum(grouped[key]) for key in sorted(grouped)])
    counts = np.array([len(grouped[key]) for key in sorted(grouped)])
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(repetitions):
        sample = rng.integers(len(sums), size=len(sums))
        draws.append(float(sums[sample].sum() / counts[sample].sum()))
    low, high = np.quantile(draws, [0.025, 0.975])
    return {"mean_delta": float(sums.sum() / counts.sum()), "ci95": [float(low), float(high)],
            "molecules": len(left), "resampled_leakage_groups": len(sums),
            "repetitions": repetitions, "seed": seed,
            "scope": "paired molecular mean; cluster-resampled; conditional on fitted out-of-fold models"}
