"""Load aligned distinct-position feature queries without fitting or transforming them."""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class FeatureBank:
    X: dict[str, np.ndarray]
    queries: list[dict]
    manifest: dict

    @property
    def observed_membership(self) -> np.ndarray:
        """Recorded-site membership, not measured selectivity or proof of negative outcomes."""
        return np.array([row["reported_positive"] for row in self.queries], dtype=bool)

    @property
    def compound_ids(self) -> np.ndarray:
        return np.array([row["compound_id"] for row in self.queries])


def load_feature_bank(root: Path | str, *, cohort: str, blocks: list[str] | None = None) -> FeatureBank:
    root = Path(root)
    latest = json.loads((root / "data/features/latest.json").read_text(encoding="utf-8"))
    manifest_path = root / latest["manifest"]
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if cohort not in manifest["cohorts"]:
        raise ValueError(f"Unknown cohort: {cohort}")
    with (manifest_path.parent / "queries.csv").open(encoding="utf-8", newline="") as handle:
        queries = list(csv.DictReader(handle))
    for row in queries:
        for field in ["source_atom_indices", "canonical_atom_indices", "cohorts"]:
            row[field] = json.loads(row[field])
        for field in ["row_id", "outcome_group_index", "n_reported_outcomes", "n_candidate_outcomes"]:
            row[field] = int(row[field])
        row["reported_positive"] = row["reported_positive"] == "True"
    mask = np.array([cohort in row["cohorts"] for row in queries])
    ids = np.array([row["query_id"] for row in queries])
    matrices = {}
    for name in blocks or list(manifest["blocks"]):
        if name not in manifest["blocks"]:
            raise ValueError(f"Unknown feature block: {name}")
        with np.load(manifest_path.parent / f"{name}.npz", allow_pickle=False) as arrays:
            if not np.array_equal(arrays["query_ids"], ids):
                raise ValueError(f"Feature/query alignment mismatch: {name}")
            if not arrays["present"][mask].all() or not np.isfinite(arrays["X"][mask]).all():
                raise ValueError(f"Incomplete features for {name}; select a complete cohort or review missing records. Individual candidate positions must not be silently dropped.")
            matrices[name] = arrays["X"][mask].copy()
    return FeatureBank(matrices, [row for row, selected in zip(queries, mask) if selected], manifest)
