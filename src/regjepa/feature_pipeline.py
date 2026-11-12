"""Resumable feature banks; one row per symmetrically distinct position."""
from __future__ import annotations

import csv
import hashlib
import json
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from rdkit import rdBase

from .representations import canonical_group_indices, digest_json, generate_geometry, group_features, morgan_atom_features

STAGES = ["conservative_bromination", "conservative_halogenation", "broad_single_outcome", "all_source_records"]


def cohort_memberships(record: dict) -> list[str]:
    result = []
    if record["eligible_conservative_bromination"]:
        result.append(STAGES[0])
    if (record["eligible_broad_single_outcome"] and
            record["reaction_class"] in {"bromination_only_evidence", "chlorination_only_evidence", "iodination_only_evidence"} and
            not (set(record["review_flags"]) - {"not_bromination_only_evidence"})):
        result.append(STAGES[1])
    if record["eligible_broad_single_outcome"]:
        result.append(STAGES[2])
    result.append(STAGES[3])
    return result


def prepare(root: Path, config: dict, cohort: str = "all_source_records", limit: int | None = None) -> dict:
    records = [json.loads(line) for line in (root / "data/processed/molecules.jsonl").read_text(encoding="utf-8").splitlines()]
    selected = [r for r in records if cohort in cohort_memberships(r)]
    selected.sort(key=lambda r: (min(STAGES.index(s) for s in cohort_memberships(r)), r["source_line"]))
    structures = {}
    for record in selected:
        key = digest_json(record["canonical_isomeric_smiles"])[:20]
        record["structure_key"] = key
        if key not in structures:
            structures[key] = {"key": key, "smiles": record["canonical_isomeric_smiles"]}
    if limit is not None:
        structures = dict(list(structures.items())[:limit])
        selected = [r for r in selected if r["structure_key"] in structures]
    provenance = json.loads((root / "data/feature_tools.json").read_text(encoding="utf-8"))
    signature = {"config": config, "rdkit": rdBase.rdkitVersion, "tools": provenance}
    run_key = digest_json(signature)[:20]
    cache = root / ".cache/feature_runs" / run_key
    cache.mkdir(parents=True, exist_ok=True)
    return {"root": root, "config": config, "cohort": cohort, "records": selected,
            "structures": list(structures.values()), "tools": provenance,
            "run_key": run_key, "cache": cache}


def save_arrays(path: Path, **arrays) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".part.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(path)


def cache_path(state: dict, block: str, key: str) -> Path:
    return state["cache"] / block / f"{key}.npz"


def _geometry_job(path: str, smiles: str, config: dict, xtb_executable: str | None = None) -> dict:
    target = Path(path)
    if target.exists() and target.with_suffix(".json").exists():
        return json.loads(target.with_suffix(".json").read_text(encoding="utf-8"))
    arrays, metadata = generate_geometry(smiles, config)
    if metadata["forcefield"] == "GFN2_xTB_pending":
        from .quantum import relax_geometry
        if not xtb_executable:
            raise ValueError("No xTB executable for the requested geometry fallback")
        arrays, metadata = relax_geometry(Path(xtb_executable), arrays, metadata, target.parent.parent / "geometry_relax" / target.stem)
    save_arrays(target, **arrays)
    target.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return metadata


def geometry_phase(state: dict, workers: int = 2) -> None:
    started = time.perf_counter()
    failures = {}
    statuses = {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        xtb_executable = str(state["root"] / state["tools"]["xtb"]["executable"])
        futures = {pool.submit(_geometry_job, str(cache_path(state, "geometry", s["key"])), s["smiles"], state["config"]["conformers"], xtb_executable): s for s in state["structures"]}
        for number, future in enumerate(as_completed(futures), 1):
            structure = futures[future]
            try:
                statuses[structure["key"]] = future.result()
            except Exception as error:
                failures[structure["key"]] = repr(error)
                print(f"Geometry failed {structure['key']}: {error}", flush=True)
            if number % 25 == 0 or number == len(futures):
                print(f"Geometry {number}/{len(futures)}; failures={len(failures)}; elapsed={time.perf_counter()-started:.1f}s", flush=True)
    (state["cache"] / "geometry_status.json").write_text(json.dumps({"statuses": statuses, "failures": failures}, indent=2) + "\n", encoding="utf-8")


def morgan_phase(state: dict) -> None:
    for structure in state["structures"]:
        path = cache_path(state, "morgan", structure["key"])
        if not path.exists():
            values = morgan_atom_features(structure["smiles"], state["config"]["morgan"])
            save_arrays(path, values=values)
    print(f"Rooted Morgan features ready for {len(state['structures'])} unique structures", flush=True)


def electronic_phase(state: dict, workers: int = 2) -> None:
    from .quantum import electronic_structure
    executable = state["root"] / state["tools"]["xtb"]["executable"]
    failures = {}
    started = time.perf_counter()
    def job(structure):
        target = cache_path(state, "electronic", structure["key"])
        if target.exists():
            return
        geometry = cache_path(state, "geometry", structure["key"])
        if not geometry.exists():
            raise FileNotFoundError("No valid conformers")
        with np.load(geometry, allow_pickle=False) as arrays:
            metadata = json.loads(geometry.with_suffix(".json").read_text(encoding="utf-8"))
            values = electronic_structure(executable, arrays, metadata, state["config"]["electronic"],
                                          min(state["config"]["feature_conformers"], len(arrays["positions"])), state["cache"] / "xtb_jobs" / structure["key"])
        save_arrays(target, values=values)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(job, s): s for s in state["structures"]}
        for number, future in enumerate(as_completed(futures), 1):
            structure = futures[future]
            try:
                future.result()
            except Exception as error:
                failures[structure["key"]] = repr(error)
                print(f"Electronic features failed {structure['key']}: {error}", flush=True)
            if number % 25 == 0 or number == len(futures):
                print(f"Electronic {number}/{len(futures)}; failures={len(failures)}; elapsed={time.perf_counter()-started:.1f}s", flush=True)
    (state["cache"] / "electronic_status.json").write_text(json.dumps({"failures": failures}, indent=2) + "\n", encoding="utf-8")


def jepa_phase(state: dict, *, random_init: bool = False, device: str = "cpu") -> None:
    import gc
    import torch
    from atom_jepa import AtomJEPA
    from atom_jepa.checkpoint import load_pretrained_encoder
    torch.set_num_threads(state["config"]["atom_jepa"]["threads"])
    block = "random_backbone" if random_init else "atom_jepa"
    pending = [s for s in state["structures"] if not cache_path(state, block, s["key"]).exists()
               and cache_path(state, "geometry", s["key"]).exists()]
    if not pending:
        print(f"{block}: no uncached structures with available geometry", flush=True)
        return
    directory = str(state["root"] / state["config"]["atom_jepa"]["checkpoint_directory"])
    torch.manual_seed(state["config"]["atom_jepa"]["random_seed"])
    if random_init:
        encoder, cfg, _ = load_pretrained_encoder(directory, device, load_weights=False)
        model = AtomJEPA(encoder, cfg).to(device)
        del _
    else:
        model = AtomJEPA.from_pretrained(directory, device=device, precision="fp32", cuequivariance=False, compile=False)
    model.eval().requires_grad_(False)
    assert all(not p.requires_grad for p in model.parameters())
    sha = hashlib.sha256()
    for name, parameter in model.named_parameters():
        sha.update(name.encode())
        sha.update(parameter.detach().cpu().contiguous().numpy().tobytes())
    from dataclasses import asdict
    model_info = {"block": block, "architecture": asdict(model.config), "parameter_sha256": sha.hexdigest(),
                  "parameter_count": sum(p.numel() for p in model.parameters()), "trainable_parameter_count": 0,
                  "seed": state["config"]["atom_jepa"]["random_seed"] if random_init else None,
                  "device": device, "precision": "fp32", "initialization": "upstream architecture constructor" if random_init else "published molecular checkpoint"}
    (state["cache"] / f"{block}_model.json").write_text(json.dumps(model_info, indent=2) + "\n", encoding="utf-8")
    failures = {}
    started = time.perf_counter()
    for number, structure in enumerate(pending, 1):
        path = cache_path(state, "geometry", structure["key"])
        try:
            with np.load(path, allow_pickle=False) as geometry:
                n = int(geometry["n_original_atoms"])
                samples = [(geometry["atomic_numbers"], position) for position in geometry["positions"][:state["config"]["feature_conformers"]]]
                outputs = model.embed(samples, per_atom=True, degrees=state["config"]["atom_jepa"]["degrees"],
                                      invariant=state["config"]["atom_jepa"]["invariant"],
                                      batch_size=state["config"]["atom_jepa"]["batch_size"])
                values = np.stack([output[:n].numpy().reshape(n, -1) for output in outputs]).mean(axis=0)
            if not np.isfinite(values).all():
                raise ValueError("Nonfinite encoder output")
            save_arrays(cache_path(state, block, structure["key"]), values=values)
        except Exception as error:
            failures[structure["key"]] = repr(error)
            print(f"{block} failed {structure['key']}: {error}", flush=True)
        if number % 10 == 0 or number == len(pending):
            print(f"{block} {number}/{len(pending)}; failures={len(failures)}; elapsed={time.perf_counter()-started:.1f}s", flush=True)
    (state["cache"] / f"{block}_status.json").write_text(json.dumps({"failures": failures}, indent=2) + "\n", encoding="utf-8")
    del model
    gc.collect()


def assemble(state: dict) -> dict:
    from .quantum import ELECTRONIC_NAMES
    output = state["root"] / "data/features" / state["run_key"]
    output.mkdir(parents=True, exist_ok=True)
    queries = []
    for record in state["records"]:
        geometry_path = cache_path(state, "geometry", record["structure_key"])
        geometry_info = json.loads(geometry_path.with_suffix(".json").read_text(encoding="utf-8")) if geometry_path.exists() else {}
        if geometry_path.exists():
            with np.load(geometry_path, allow_pickle=False) as geometry_arrays:
                selected_statuses = geometry_arrays["forcefield_status"][:state["config"]["feature_conformers"]]
                selected_converged = bool((selected_statuses == 0).all())
        else:
            selected_converged = False
        for group_index, (source_group, canonical_group) in enumerate(zip(record["candidate_outcome_groups"], canonical_group_indices(record))):
            queries.append({"row_id": len(queries), "query_id": f"{record['compound_id']}:outcome{group_index}",
                            "compound_id": record["compound_id"], "structure_key": record["structure_key"],
                            "outcome_group_index": group_index, "source_atom_indices": source_group,
                            "canonical_atom_indices": canonical_group,
                            "structure_position_key": record["structure_key"] + ":" + ",".join(map(str, canonical_group)),
                            "reported_positive": group_index in record["observed_outcome_group_indices"],
                            "n_reported_outcomes": record["n_observed_outcomes"],
                            "n_candidate_outcomes": record["n_candidate_outcomes"],
                            "cohorts": cohort_memberships(record),
                            "geometry_method": geometry_info.get("forcefield", "unavailable"),
                            "geometry_conformers_available": geometry_info.get("num_conformers", 0),
                            "feature_conformers_used": min(state["config"]["feature_conformers"], geometry_info.get("num_conformers", 0)),
                            "selected_geometry_converged": selected_converged})
    with (output / "queries.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(queries[0]), lineterminator="\n")
        writer.writeheader()
        for row in queries:
            writer.writerow({k: json.dumps(v) if isinstance(v, list) else v for k, v in row.items()})
    architecture = json.loads((state["root"] / state["config"]["atom_jepa"]["checkpoint_directory"] / "config.json").read_text())["eqv3_cfg"]
    degree_count = architecture["lmax"]+1 if state["config"]["atom_jepa"]["degrees"] == "all" else 1
    latent_width = architecture["num_channels"] * degree_count
    dimensions = {"morgan": len(state["config"]["morgan"]["radii"])*state["config"]["morgan"]["bits_per_radius"],
                  "electronic": len(ELECTRONIC_NAMES), "atom_jepa": latent_width, "random_backbone": latent_width}
    matrices = {}
    for block, dimension in dimensions.items():
        values = np.full((len(queries), dimension), np.nan, dtype=np.float32)
        for record in state["records"]:
            path = cache_path(state, block, record["structure_key"])
            if not path.exists():
                continue
            with np.load(path, allow_pickle=False) as cached:
                grouped = group_features(cached["values"], canonical_group_indices(record))
            if grouped.shape[1] != dimension:
                raise ValueError(f"Unexpected feature dimension: {block} {grouped.shape}")
            indices = [row["row_id"] for row in queries if row["compound_id"] == record["compound_id"]]
            values[indices] = grouped
        matrices[block] = values
    matrices["jepa_plus_qm"] = np.concatenate([matrices["atom_jepa"], matrices["electronic"]], axis=1)
    blocks = {}
    for name, values in matrices.items():
        present = np.isfinite(values).all(axis=1)
        path = output / f"{name}.npz"
        save_arrays(path, X=values, present=present, query_ids=np.array([q["query_id"] for q in queries]))
        blocks[name] = {"shape": list(values.shape), "available_queries": int(present.sum()),
                        "available_records": len({q["compound_id"] for q, ok in zip(queries, present) if ok}),
                        "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    matched = np.logical_and.reduce([np.isfinite(values).all(axis=1) for values in matrices.values()])
    cohort_info = {}
    for cohort in STAGES:
        mask = np.array([cohort in row["cohorts"] for row in queries])
        cohort_info[cohort] = {"records": len({row["compound_id"] for row in queries if cohort in row["cohorts"]}),
                               "queries": int(mask.sum()), "matched_queries": int((mask & matched).sum()),
                               "matched_records": len({row["compound_id"] for row, ok in zip(queries, mask & matched) if ok})}
    model_info = {name: json.loads(path.read_text(encoding="utf-8")) for name in ["atom_jepa", "random_backbone"]
                  if (path := state["cache"] / f"{name}_model.json").exists()}
    if len(model_info) == 2:
        if model_info["atom_jepa"]["architecture"] != model_info["random_backbone"]["architecture"]:
            raise ValueError("Pretrained and random encoders have different architecture configurations")
        if model_info["atom_jepa"]["parameter_sha256"] == model_info["random_backbone"]["parameter_sha256"]:
            raise ValueError("Random encoder unexpectedly matches pretrained parameters")
    # Portable, pickle-free geometry archive: cache paths and external tools are not
    # required to inspect or reuse the exact coordinates supplied to the encoders.
    geometry_info = []
    atom_numbers, positions, atom_offsets, conformer_offsets, structure_indices = [], [], [0], [0], []
    conformer_statuses, conformer_energies, original_counts = [], [], []
    for structure in state["structures"]:
        path = cache_path(state, "geometry", structure["key"])
        if not path.exists():
            continue
        metadata = json.loads(path.with_suffix(".json").read_text(encoding="utf-8"))
        with np.load(path, allow_pickle=False) as arrays:
            atom_numbers.append(arrays["atomic_numbers"])
            atom_offsets.append(atom_offsets[-1] + len(arrays["atomic_numbers"]))
            original_counts.append(int(arrays["n_original_atoms"]))
            for coords, status, energy in zip(arrays["positions"], arrays["forcefield_status"], arrays["forcefield_energy_kcal_mol"]):
                positions.append(coords)
                conformer_offsets.append(conformer_offsets[-1] + len(coords))
                structure_indices.append(len(geometry_info))
                conformer_statuses.append(status)
                conformer_energies.append(energy)
        geometry_info.append(dict(structure_key=structure["key"], **metadata))
    if geometry_info:
        save_arrays(output / "conformers.npz", structure_keys=np.array([g["structure_key"] for g in geometry_info]),
                    atomic_numbers=np.concatenate(atom_numbers), atom_offsets=np.array(atom_offsets),
                    n_original_atoms=np.array(original_counts), positions=np.concatenate(positions),
                    conformer_atom_offsets=np.array(conformer_offsets), conformer_structure_index=np.array(structure_indices),
                    forcefield_status=np.array(conformer_statuses), forcefield_energy_kcal_mol=np.array(conformer_energies))
        (output / "geometry_metadata.jsonl").write_text("".join(json.dumps(g)+"\n" for g in geometry_info), encoding="utf-8", newline="\n")
    manifest = {"schema_version": 1, "run_key": state["run_key"], "config": state["config"], "tools": state["tools"],
                "source_records": len(state["records"]), "unique_structures": len(state["structures"]),
                "query_rows": len(queries), "cohorts": cohort_info, "blocks": blocks, "models": model_info,
                "queries_sha256": hashlib.sha256((output / "queries.csv").read_bytes()).hexdigest(),
                "geometry": {"unique_structures": len(geometry_info), "conformers": len(conformer_statuses),
                             "archive": "conformers.npz", "archive_sha256": hashlib.sha256((output / "conformers.npz").read_bytes()).hexdigest() if geometry_info else None,
                             "metadata_sha256": hashlib.sha256((output / "geometry_metadata.jsonl").read_bytes()).hexdigest() if geometry_info else None},
                "electronic_feature_names": ELECTRONIC_NAMES,
                "query_policy": "one row per exact substitution outcome; arithmetic mean of member-atom features, no multiplicity weights",
                "labels": "reported_positive is an observed-site membership annotation, not a measured product fraction; multiple outcomes remain multi-positive",
                "scaling": "none; all fitted scalers must be trained within downstream training folds",
                "split_status": "cohorts are cumulative data stages; train/validation/test indices are not inferred from feature extraction",
                "source_manifest_sha256": hashlib.sha256((state["root"] / "data/raw/manifest.json").read_bytes()).hexdigest(),
                "dependency_lock_sha256": hashlib.sha256((state["root"] / "requirements-features.lock.txt").read_bytes()).hexdigest(),
                "implementation_sha256": {p.relative_to(state["root"]).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in (state["root"] / "src/regjepa").glob("*.py")}}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n")
    (state["root"] / "data/features/latest.json").write_text(json.dumps({"run_key": state["run_key"], "manifest": (output / "manifest.json").relative_to(state["root"]).as_posix()}, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"output": str(output), "cohorts": cohort_info, "blocks": blocks}, indent=2), flush=True)
    return manifest
