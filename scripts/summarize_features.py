"""Verify feature alignment, frozen-model controls and geometry coverage; write a report."""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    latest = json.loads((ROOT / "data/features/latest.json").read_text())
    manifest_path = ROOT / latest["manifest"]
    manifest = json.loads(manifest_path.read_text())
    directory = manifest_path.parent
    assert hashlib.sha256((directory / "queries.csv").read_bytes()).hexdigest() == manifest["queries_sha256"]
    assert hashlib.sha256((directory / "geometry_metadata.jsonl").read_bytes()).hexdigest() == manifest["geometry"]["metadata_sha256"]
    with (directory / "queries.csv").open(encoding="utf-8", newline="") as handle:
        queries = list(csv.DictReader(handle))
    ids = np.array([q["query_id"] for q in queries])
    assert len(ids) == len(set(ids)) == manifest["query_rows"]
    metrics = {}
    for name, spec in manifest["blocks"].items():
        path = directory / f"{name}.npz"
        assert hashlib.sha256(path.read_bytes()).hexdigest() == spec["sha256"]
        with np.load(path, allow_pickle=False) as arrays:
            assert list(arrays["X"].shape) == spec["shape"]
            np.testing.assert_array_equal(arrays["query_ids"], ids)
            assert arrays["present"].dtype == np.dtype(bool)
            assert arrays["present"].shape == (len(ids),)
            np.testing.assert_array_equal(arrays["present"], np.isfinite(arrays["X"]).all(axis=1))
            assert int(arrays["present"].sum()) == spec["available_queries"]
            assert np.isfinite(arrays["X"][arrays["present"]]).all()
            metrics[name] = {"rows": len(ids), "columns": arrays["X"].shape[1],
                             "available": int(arrays["present"].sum()),
                             "all_rows_available": bool(arrays["present"].all())}
    with np.load(directory / "atom_jepa.npz", allow_pickle=False) as jepa, np.load(directory / "electronic.npz", allow_pickle=False) as qm, np.load(directory / "jepa_plus_qm.npz", allow_pickle=False) as combined:
        np.testing.assert_array_equal(combined["X"], np.concatenate([jepa["X"], qm["X"]], axis=1))
    a, b = manifest["models"]["atom_jepa"], manifest["models"]["random_backbone"]
    assert a["architecture"] == b["architecture"]
    assert a["parameter_count"] == b["parameter_count"]
    assert a["parameter_sha256"] != b["parameter_sha256"]
    assert a["trainable_parameter_count"] == b["trainable_parameter_count"] == 0
    geometries = [json.loads(line) for line in (directory / "geometry_metadata.jsonl").read_text().splitlines()]
    assert hashlib.sha256((directory / "conformers.npz").read_bytes()).hexdigest() == manifest["geometry"]["archive_sha256"]
    with np.load(directory / "conformers.npz", allow_pickle=False) as coords:
        assert np.isfinite(coords["positions"]).all()
        assert len(coords["structure_keys"]) == len(geometries)
        assert coords["atom_offsets"][-1] == len(coords["atomic_numbers"])
        assert coords["conformer_atom_offsets"][-1] == len(coords["positions"])
        assert len(coords["conformer_structure_index"]) == manifest["geometry"]["conformers"]
        for i, g in enumerate(geometries):
            assert sum(coords["conformer_structure_index"] == i) == g["num_conformers"]
    methods = dict(Counter(g["forcefield"] for g in geometries))
    counts = dict(Counter(g["num_conformers"] for g in geometries))
    summary = {"run_key": manifest["run_key"], "blocks": metrics, "cohorts": manifest["cohorts"],
               "geometry_methods": methods, "conformers_per_structure": counts,
               "geometry_conformers_total": manifest["geometry"]["conformers"],
               "source_records": manifest["source_records"], "unique_structures": manifest["unique_structures"],
               "selected_geometry_nonconverged_records": len({q["compound_id"] for q in queries if q["selected_geometry_converged"] == "False"}),
               "frozen_parameter_count": a["parameter_count"],
               "all_checks_passed": True}
    reports = ROOT / "reports"
    (reports / "feature_extraction_summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    block_rows = "\n".join(f"| {name} | {m['columns']} | {m['available']} / {m['rows']} |" for name, m in metrics.items())
    cohort_rows = "\n".join(f"| {name} | {info['records']} | {info['queries']} | {info['matched_queries']} |" for name, info in manifest["cohorts"].items())
    method_rows = "\n".join(f"| {name} | {count} |" for name, count in methods.items())
    text = f'''# Frozen feature extraction: first shared-conformer bank

Run: `{manifest['run_key']}`. All validation checks passed. No site-prediction head or downstream scaling has been fitted.

## Completed artifacts

Features cover {manifest['source_records']} source records / {manifest['unique_structures']} unique substrates. The {manifest['query_rows']} rows are symmetrically distinct candidate outcomes per source record. Atom indices, source annotations and class memberships are in `queries.csv`. Repeated source records share the same canonical-structure computations.

| Block | Columns | Available distinct-position rows |
|---|---:|---:|
{block_rows}

Atom-JEPA features concatenate the final-layer l=0 scalars and per-channel norms of l=1 and l=2, 256 channels each. The pretrained and random models have exactly the same configuration and {a['parameter_count']:,} frozen parameters each. The random backbone uses upstream initialization with seed 42 and loads no pretrained parameters; its parameter checksum differs from the pretrained checksum.

Rooted Morgan fingerprints use radii 1, 2 and 3 in separate 2048-bit blocks with chirality enabled. Electronic descriptors are native GFN1-xTB Mulliken charge, CM5 charge, Fukui f(+), f(-), f(0), and total Wiberg bond order. They are single-point calculations at the shared geometry with GBSA methanol, electronic temperature 300 K, and explicit molecular charge/spin. Atomic element/index correspondence, molecular charge sums and Fukui sum rules were checked. This is a standardized descriptor environment, not a reconstruction of the literature reaction conditions. JEPA + QM is the literal 768+6 concatenation.

## Approved cumulative data stages

| Stage | Source records | Distinct-position queries | Queries shared by every block |
|---|---:|---:|---:|
{cohort_rows}

The second stage adds 28 conservative chlorination and 47 conservative iodination records to the initial 299 brominations. Later stages retain mixed/ambiguous and multi-outcome cases with their original flags. No-choice or multi-positive records in the final stage must not be silently treated as ordinary single-label ranking examples. These stages are cohorts; actual supervised train/validation/test assignments remain a downstream step.

## Geometry protocol and exceptions

ETKDGv3, explicit hydrogens, fixed seed 42, no RMS pruning, enforced specified chirality, and MMFF94 relaxation up to 500 iterations follow the Atom-JEPA procedure. Random coordinates are used only if normal embedding fails. Source atom order is mapped to a canonical molecular graph before adding H, and original-to-canonical correspondence is kept for outcome queries.

The portable `conformers.npz` archive contains **{manifest['geometry']['conformers']:,} conformers** for every unique substrate; `geometry_metadata.jsonl` records methods, convergence and any rejected fallback conformers.

| Relaxation method | Unique structures |
|---|---:|
{method_rows}

For 13 molecules without MMFF parameters, a fully parameterized UFF fallback is explicit. Three SF5 substrates use xTB GFN-FF relaxation; fallback conformers that break an input-bond distance check are discarded without using outcome labels. They retain 6, 8 and 9 valid conformers. One selenium substrate uses GFN2-xTB relaxation. These are recorded departures from the paper's main MMFF94 procedure and can be excluded in a method-sensitivity analysis.

MMFF nonconvergence is recorded, and last coordinates are retained as in the upstream package. {summary['selected_geometry_nonconverged_records']} source record(s) use a nonconverged first geometry; the metadata allow excluding them consistently across every feature block.

## Compute scope and comparison policy

The completed bank uses the same **first retained conformer** for both encoders and electronic descriptors. This CPU-only first pass is not a ten-conformer ensemble. All generated conformers are saved, so `feature_conformers` can later be raised to 10 to produce matched ensemble means; molecules with fewer valid fallback conformers use all available ones. The configuration change produces a separate run key.

Each symmetry group is queried once. Member-atom descriptors are averaged into its single vector; group size is not a model input, weight, probability multiplier or extra training example. No outcome-product structure, positive-site label or measured yield is an input to the extractors. All feature arrays are unscaled; fit scalers only within training folds. Features currently describe the substrate, not the electrophile/condition identity.

## Files and rerun

Feature directory: `data/features/{manifest['run_key']}/`.

- `queries.csv`: aligned outcome-query identifiers, canonical/source atom correspondence, annotations and geometry QC.
- `atom_jepa.npz`, `random_backbone.npz`, `morgan.npz`, `electronic.npz`, `jepa_plus_qm.npz`: `X`, `present`, and `query_ids`, readable with `numpy.load(..., allow_pickle=False)`.
- `conformers.npz`, `geometry_metadata.jsonl`: portable coordinates and metadata.
- `manifest.json`: config, checkpoint revision/checksums, model initialization checksums, block hashes, row coverage and source-code fingerprints.

Install `requirements-features.txt`, run `scripts/setup_feature_tools.py`, then `scripts/extract_features.py --phase all`. The script resumes cached work. Run `scripts/summarize_features.py` to verify the exported bank. Large model weights, native binaries and per-calculation scratch logs remain in ignored `.cache/`; the compact feature/geometry exports and provenance are retained in the repository.

Sources: [Atom-JEPA paper, Appendix H](https://arxiv.org/html/2610.08400v1), [reference conformer implementation](https://github.com/khelverskovp/atom-jepa/blob/main/atom_jepa/conformers.py), [xTB properties](https://xtb-docs.readthedocs.io/en/latest/properties.html), [xTB Fukui calculations](https://xtb-docs.readthedocs.io/en/latest/sp.html).
'''
    (reports / "feature_extraction.md").write_text(text, encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
