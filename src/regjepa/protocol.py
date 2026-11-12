"""Reviewable split proposals and input locks; preparation never trains a model."""
from __future__ import annotations

import hashlib
import json
import platform
from collections import Counter
from importlib.metadata import version
from pathlib import Path

import numpy as np
from rdkit import rdBase

from .evaluation import CandidateTable, ranking_records, summarize_records
from .feature_io import load_feature_bank
from .grouping import build_grouping, stable_digest
from .splits import build_split_plan, validate_split_plan


BLOCK_WIDTHS = {"atom_jepa": 768, "random_backbone": 768, "morgan": 6144,
                "electronic": 6, "jepa_plus_qm": 774}


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def locked_input_sha256(root: Path, relative: str) -> str:
    data = (root / relative).read_bytes()
    if relative.startswith("data/processed/"):
        # Derived tables were generated with native Windows line endings, while
        # Git stores them as LF. Lock the same text across fresh checkouts without
        # rewriting those files. Raw snapshots and feature exports stay exact.
        data = data.replace(b"\r\n", b"\n")
    return hashlib.sha256(data).hexdigest()


def write_json(path: Path, value: object, *, exclusive: bool = False) -> None:
    text = json.dumps(value, indent=2, allow_nan=False) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    if exclusive and path.exists():
        if path.read_text(encoding="utf-8") != text:
            raise FileExistsError(f"Refusing to overwrite an immutable artifact: {path}")
        return
    with path.open("x" if exclusive else "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def validate_config(config: dict) -> None:
    if config["schema_version"] != 1 or config["cohort"] != "conservative_bromination":
        raise ValueError("This protocol version supports only the conservative single-outcome pilot")
    if not config["blocks"] or len(config["blocks"]) != len(set(config["blocks"])) or not set(config["blocks"]) <= set(BLOCK_WIDTHS):
        raise ValueError("Unknown or repeated feature blocks; metadata must never be model inputs")
    if config["regime"] != "source_series_structure":
        raise ValueError("This freezable protocol requires combined identity/source/series/structure constraints")
    if not config["include_full_budget"]:
        raise ValueError("This protocol includes the full training budget")
    if not config["subset_seeds"] or len(set(config["subset_seeds"])) != len(config["subset_seeds"]):
        raise ValueError("Distinct subset seeds are required")
    if not config["l2_grid"] or len(set(config["l2_grid"])) != len(config["l2_grid"]) or any(not np.isfinite(v) or v <= 0 for v in config["l2_grid"]):
        raise ValueError("L2 grid must be finite and positive")
    if config["combined_penalties"] not in {"shared", "separate_jepa_qm"}:
        raise ValueError("Unknown concatenated-feature penalty policy")
    expected = {"tie_policy": "uniform_expected_metrics_for_exact_ties",
                "dense_scaling": "train_only_equal_molecule_within_rms", "morgan_scaling": "unscaled"}
    if any(config[key] != value for key, value in expected.items()):
        raise ValueError("Unsupported metric or scaling policy")


def input_locks(root: Path, manifest_path: Path, blocks: list[str]) -> dict[str, str]:
    paths = [manifest_path, manifest_path.parent / "queries.csv", root / "data/processed/molecules.jsonl",
             root / "data/raw/manifest.json", root / "data/raw/db_regioselectivity/compounds_doi.csv",
             root / "data/raw/db_regioselectivity/bib/library.bib", root / "requirements-evaluation.txt",
             root / "requirements.txt"]
    paths.extend(manifest_path.parent / f"{block}.npz" for block in blocks)
    paths.extend(sorted((root / "src/regjepa").glob("*.py")))
    paths.extend(root / "scripts" / name for name in ["prepare_evaluation.py", "run_linear_evaluation.py"])
    relatives = [path.relative_to(root).as_posix() for path in paths]
    return {relative: locked_input_sha256(root, relative) for relative in relatives}


def verify_inputs(root: Path, plan: dict) -> None:
    for relative, digest in plan["input_sha256"].items():
        path = root / relative
        if not path.is_file() or locked_input_sha256(root, relative) != digest:
            raise ValueError(f"Locked input changed or disappeared: {relative}")
    latest = json.loads((root / "data/features/latest.json").read_text(encoding="utf-8"))
    if latest["manifest"] != plan["feature_manifest"]:
        raise ValueError("Latest feature bank differs from the locked split proposal")
    current = {"python": platform.python_version(), "numpy": np.__version__,
               "rdkit": rdBase.rdkitVersion, "scipy": version("scipy")}
    if current != plan["versions"]:
        raise ValueError("Runtime versions differ from the locked proposal; prepare a new proposal")


def prepare_protocol(root: Path, config: dict) -> tuple[dict, dict, dict]:
    validate_config(config)
    bank = load_feature_bank(root, cohort=config["cohort"], blocks=config["blocks"])
    table = CandidateTable.from_queries(bank.queries)
    latest = json.loads((root / "data/features/latest.json").read_text(encoding="utf-8"))
    manifest_path = root / latest["manifest"]
    for block in config["blocks"]:
        if bank.X[block].shape != (table.n_rows, BLOCK_WIDTHS[block]):
            raise ValueError(f"Unexpected feature dimensions for {block}")
        if file_sha256(manifest_path.parent / f"{block}.npz") != bank.manifest["blocks"][block]["sha256"]:
            raise ValueError(f"Feature block checksum mismatch: {block}")
    if file_sha256(manifest_path.parent / "queries.csv") != bank.manifest["queries_sha256"]:
        raise ValueError("Query-table checksum mismatch")
    records = [json.loads(line) for line in (root / "data/processed/molecules.jsonl").read_text(encoding="utf-8").splitlines()]
    source = {record["compound_id"]: record for record in records}
    for query in bank.queries:
        record = source[query["compound_id"]]
        if (query["structure_key"] != stable_digest(record["canonical_isomeric_smiles"])[:20] or
                query["n_candidate_outcomes"] != record["n_candidate_outcomes"] or
                query["reported_positive"] != (query["outcome_group_index"] in record["observed_outcome_group_indices"])):
            raise ValueError("Processed source identity/annotation does not match the frozen feature queries")
    if file_sha256(root / "data/raw/manifest.json") != bank.manifest["source_manifest_sha256"]:
        raise ValueError("Raw source manifest differs from feature-extraction provenance")
    grouping = build_grouping(records, set(table.compound_ids),
                              (root / "data/raw/db_regioselectivity/bib/library.bib").read_text(encoding="utf-8"),
                              config["fingerprint"], regime=config["regime"])
    folds = build_split_plan(list(table.compound_ids), grouping["compound_to_group"], config)
    plan = {"schema_version": 1, "status": "proposal", "config": config,
            "feature_manifest": latest["manifest"], "feature_run_key": bank.manifest["run_key"],
            "compound_to_group": grouping["compound_to_group"], "outer_folds": folds,
            "grouping_digest": stable_digest(grouping),
            "input_sha256": input_locks(root, manifest_path, config["blocks"]),
            "input_hash_policy": "exact bytes; data/processed/** only canonicalizes CRLF to LF as Git does",
            "versions": {"python": platform.python_version(), "numpy": np.__version__,
                         "rdkit": rdBase.rdkitVersion, "scipy": version("scipy")},
            "source_limitations": grouping["limitations"],
            "fit_status": "No real-data scorers trained; split preparation and analytical baselines only."}
    plan["proposal_id"] = stable_digest(plan)[:20]
    baseline = summarize_records(ranking_records(np.zeros(table.n_rows), table))
    duplicates = []
    if "morgan" in bank.X:
        for compound, indices in zip(table.compound_ids, table.rows):
            features = bank.X["morgan"][indices]
            pairs = [[i, j] for i in range(len(indices)) for j in range(i + 1, len(indices))
                     if np.array_equal(features[i], features[j])]
            if pairs:
                duplicates.append({"compound_id": compound, "local_row_pairs": pairs})
    sizes = [len(group["selected_members"]) for group in grouping["groups"]]
    summary = {"proposal_id": plan["proposal_id"], "split_status": "proposal_not_frozen",
               "cohort": config["cohort"], "rows": table.n_rows, "molecules": len(table.rows),
               "structures": len(set(table.structure_keys)),
               "candidate_count_distribution": dict(sorted(Counter(map(len, table.rows)).items())),
               "analytical_uniform_metrics": baseline,
               "constant_first_row_top1": float(np.mean([positive == 0 for positive in table.positives])),
               "identical_morgan_candidate_vectors": duplicates,
               "leakage_groups": len(sizes), "largest_leakage_group": max(sizes),
               "leakage_group_size_distribution": dict(sorted(Counter(sizes).items())),
               "maximum_between_group_tanimoto": grouping["maximum_between_group_selected_tanimoto"],
               "unresolved_bibliographic_records": sorted(k for k, v in grouping["source_entities"].items()
                                                          if v["bibliographic_identity_unresolved"]),
               "outer_counts": [{"fold": f["fold"], "train": len(f["train_compounds"]), "test": len(f["test_compounds"])} for f in folds],
               "actual_learning_curve_labels": sorted({(str(t["budget"]), t["actual_labels"]) for f in folds for t in f["training_sets"]}),
               "geometry_fallback_records": sorted({q["compound_id"] for q in bank.queries if q["geometry_method"] != "MMFF94"}),
               "source_limitations": grouping["limitations"], "full_experiment_run": False}
    return plan, grouping, summary


def load_split_plan(root: Path, path: Path, *, require_frozen: bool = True) -> dict:
    plan = json.loads(path.read_text(encoding="utf-8"))
    if require_frozen and plan.get("status") != "frozen":
        raise ValueError("Experiment requires an explicitly frozen split artifact; proposals cannot be used for fitting")
    validate_config(plan["config"])
    verify_inputs(root, plan)
    validate_split_plan(plan["outer_folds"], sorted(plan["compound_to_group"]), plan["compound_to_group"])
    proposal = {k: v for k, v in plan.items() if k not in {"proposal_id", "frozen_proposal_sha256"}}
    proposal["status"] = "proposal"
    if stable_digest(proposal)[:20] != plan["proposal_id"]:
        raise ValueError("Split proposal content does not match its identifier")
    return plan


def freeze_protocol(root: Path, proposal_path: Path) -> Path:
    plan = load_split_plan(root, proposal_path, require_frozen=False)
    if plan["status"] != "proposal":
        raise ValueError("Only a reviewed proposal can be explicitly frozen")
    grouping_path = proposal_path.with_name("grouping.json")
    grouping = json.loads(grouping_path.read_text(encoding="utf-8"))
    if stable_digest(grouping) != plan["grouping_digest"] or grouping["compound_to_group"] != plan["compound_to_group"]:
        raise ValueError("Grouping audit differs from split proposal")
    # Regenerate from locked metadata/graphs, not labels or model performance.
    rebuilt, rebuilt_grouping, _ = prepare_protocol(root, plan["config"])
    if rebuilt != plan or rebuilt_grouping != grouping:
        raise ValueError("Split proposal cannot be reproduced from locked inputs")
    frozen = dict(plan, status="frozen", frozen_proposal_sha256=file_sha256(proposal_path))
    output = proposal_path.with_name("frozen.json")
    write_json(output, frozen, exclusive=True)
    return output
