"""Bind exploratory follow-up plans; this does not create a new unseen test population."""
from __future__ import annotations

import json
import re
from pathlib import Path

from regjepa.grouping import stable_digest
from regjepa.protocol import file_sha256, load_split_plan, write_json

from .features import validate_arm


def validate_config(config: dict) -> None:
    if config["schema_version"] != 1 or config["cohort"] != "conservative_bromination":
        raise ValueError("Only the current cohort is supported")
    expected = {"dense_scaling": "train_only_equal_molecule_global_rms",
                "morgan_linear_scaling": "unscaled", "tie_policy": "uniform_expected_exact_ties",
                "scope": "adaptive_exploratory_followup_on_previously_evaluated_molecules"}
    if any(config[k] != v for k, v in expected.items()):
        raise ValueError("Unsupported protocol/scaling/claim policy")
    if config["l2_grid"] != [1e-6, 1e-5, 1e-4, .001, .01, .1, 1., 10., 100., 1000.]:
        raise ValueError("This protocol fixes the widened ten-point grid")
    if config["metrics_with_paired_intervals"] != ["top1", "nll"]:
        raise ValueError("Both paired metrics are required")
    if config["bootstrap_repetitions"] != 10000 or config["bootstrap_seed"] != 42:
        raise ValueError("This proposal fixes the paired bootstrap settings")
    if config["optimizer"] != {"max_iterations": 20000, "gradient_tolerance": 1e-7,
                                "accepted_gradient": 1e-5, "mlp_hidden": 16,
                                "mlp_restart_seeds": [11, 23]}:
        raise ValueError("Unknown optimizer/readout specification")
    ids = [arm["id"] for arm in config["arms"]]
    if len(ids) != len(set(ids)) or not ids:
        raise ValueError("Repeated/missing arms")
    for arm in config["arms"]:
        validate_arm(arm)
        if not re.fullmatch("[a-z][a-z0-9_]*", arm["id"]):
            raise ValueError("Unsafe arm identifier")
        if arm["penalty_mode"] not in {"shared", "separate"} or (arm["head"] == "mlp" and arm["penalty_mode"] != "shared"):
            raise ValueError("Unsupported penalty mode")
        if arm["penalty_mode"] == "separate" and {p["penalty_group"] for p in arm["parts"]} != {"qm", "representation"}:
            raise ValueError("Separate penalties require both representation and QM columns")
    pairs = [tuple(pair) for pair in config["comparisons"]]
    if len(pairs) != len(set(pairs)) or any(a not in ids or b not in ids or a == b for a, b in pairs):
        raise ValueError("Invalid paired contrasts")


def prepare(root: Path, config: dict) -> dict:
    validate_config(config)
    baseline_path = root / config["baseline_frozen"]
    if config["baseline_frozen"] != "data/splits/af05d8348e60ad38b6e4/frozen.json":
        raise ValueError("The reviewed follow-up preserves this exact baseline split schedule")
    baseline = load_split_plan(root, baseline_path)
    paths = sorted((root / "src/regjepa_followup").glob("*.py"))
    paths += [root / f"scripts/{name}" for name in ["prepare_followup.py", "run_followup.py", "validate_followup.py"]]
    paths += [root / "requirements-followup.txt", baseline_path]
    plan = {"schema_version": 1, "status": "proposal", "config": config,
            "baseline_proposal_id": baseline["proposal_id"],
            "baseline_frozen_sha256": file_sha256(baseline_path),
            "feature_run_key": baseline["feature_run_key"],
            "compound_to_group": baseline["compound_to_group"], "outer_folds": baseline["outer_folds"],
            "input_sha256": {p.relative_to(root).as_posix(): file_sha256(p) for p in paths},
            "baseline_versions": baseline["versions"],
            "claim": "Exploratory analysis on an already inspected benchmark; not fresh confirmatory or OOD evidence."}
    plan["proposal_id"] = stable_digest(plan)[:20]
    return plan


def load(root: Path, path: Path, *, frozen: bool = True) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("status") not in {"proposal", "frozen"} or (frozen and value["status"] != "frozen"):
        raise ValueError("An explicitly frozen follow-up is required for fitting")
    rebuilt = prepare(root, value["config"])
    if value["status"] == "proposal":
        if value != rebuilt:
            raise ValueError("Proposal cannot be exactly reproduced")
    else:
        proposal_path = path.with_name("proposal.json")
        if not proposal_path.is_file() or json.loads(proposal_path.read_text()) != rebuilt:
            raise ValueError("Missing/changed bound proposal")
        expected = dict(rebuilt, status="frozen", frozen_proposal_sha256=file_sha256(proposal_path))
        if value != expected:
            raise ValueError("Frozen follow-up binding cannot be reproduced")
    return value


def freeze(root: Path, path: Path) -> Path:
    plan = load(root, path, frozen=False)
    if plan["status"] != "proposal":
        raise ValueError("Freeze an explicit proposal, not an existing frozen artifact")
    target = path.with_name("frozen.json")
    write_json(target, dict(plan, status="frozen", frozen_proposal_sha256=file_sha256(path)), exclusive=True)
    return target
