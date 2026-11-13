"""Explicit feature selections; metadata and observed outcomes never become inputs."""
from __future__ import annotations

import numpy as np
from regjepa.feature_io import FeatureBank

WIDTHS = {"electronic": 6, "atom_jepa": 768, "random_backbone": 768, "morgan": 6144}
QM_NAMES = ["Mulliken", "CM5", "fplus", "fminus", "fzero", "Wiberg"]


def validate_arm(arm: dict) -> None:
    if arm["head"] not in {"linear", "mlp"} or arm["budgets"] not in [["full"], [50, 100, "full"]]:
        raise ValueError("Unsupported head or budget schedule")
    if not arm["parts"]:
        raise ValueError("An arm needs explicit feature columns")
    occupied = set()
    for part in arm["parts"]:
        block, columns = part["block"], part["columns"]
        if block not in WIDTHS or not columns or len(columns) != len(set(columns)):
            raise ValueError("Unknown block or repeated/empty columns")
        if any(type(c) is not int or not 0 <= c < WIDTHS[block] for c in columns):
            raise ValueError("Feature column out of range")
        if part["penalty_group"] not in {"qm", "representation"}:
            raise ValueError("Unsupported penalty component")
        for column in columns:
            if (block, column) in occupied:
                raise ValueError("Repeated feature across parts")
            occupied.add((block, column))


def materialize(bank: FeatureBank, arm: dict) -> tuple[np.ndarray, list[str], np.ndarray]:
    validate_arm(arm)
    arrays, groups, morgan = [], [], []
    for part in arm["parts"]:
        block, columns = part["block"], part["columns"]
        values = bank.X[block][:, columns]
        if not np.isfinite(values).all():
            raise ValueError("Nonfinite feature input")
        arrays.append(values)
        groups.extend([part["penalty_group"]] * len(columns))
        morgan.extend([block == "morgan"] * len(columns))
    return np.concatenate(arrays, axis=1).astype(np.float64), groups, np.array(morgan, dtype=bool)
