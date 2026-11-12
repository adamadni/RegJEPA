# Frozen feature extraction: first shared-conformer bank

Run: `d76c42313dcf5e17ac96`. All validation checks passed. No site-prediction head or downstream scaling has been fitted.

## Completed artifacts

Features cover 535 source records / 534 unique substrates. The 2290 rows are symmetrically distinct candidate outcomes per source record. Atom indices, source annotations and class memberships are in `queries.csv`. Repeated source records share the same canonical-structure computations.

| Block | Columns | Available distinct-position rows |
|---|---:|---:|
| morgan | 6144 | 2290 / 2290 |
| electronic | 6 | 2290 / 2290 |
| atom_jepa | 768 | 2290 / 2290 |
| random_backbone | 768 | 2290 / 2290 |
| jepa_plus_qm | 774 | 2290 / 2290 |

Atom-JEPA features concatenate the final-layer l=0 scalars and per-channel norms of l=1 and l=2, 256 channels each. The pretrained and random models have exactly the same configuration and 39,035,904 frozen parameters each. The random backbone uses upstream initialization with seed 42 and loads no pretrained parameters; its parameter checksum differs from the pretrained checksum.

Rooted Morgan fingerprints use radii 1, 2 and 3 in separate 2048-bit blocks with chirality enabled. Electronic descriptors are native GFN1-xTB Mulliken charge, CM5 charge, Fukui f(+), f(-), f(0), and total Wiberg bond order. They are single-point calculations at the shared geometry with GBSA methanol, electronic temperature 300 K, and explicit molecular charge/spin. Atomic element/index correspondence, molecular charge sums and Fukui sum rules were checked. This is a standardized descriptor environment, not a reconstruction of the literature reaction conditions. JEPA + QM is the literal 768+6 concatenation.

## Approved cumulative data stages

| Stage | Source records | Distinct-position queries | Queries shared by every block |
|---|---:|---:|---:|
| conservative_bromination | 299 | 1264 | 1264 |
| conservative_halogenation | 374 | 1597 | 1597 |
| broad_single_outcome | 495 | 2133 | 2133 |
| all_source_records | 535 | 2290 | 2290 |

The second stage adds 28 conservative chlorination and 47 conservative iodination records to the initial 299 brominations. Later stages retain mixed/ambiguous and multi-outcome cases with their original flags. No-choice or multi-positive records in the final stage must not be silently treated as ordinary single-label ranking examples. These stages are cohorts; actual supervised train/validation/test assignments remain a downstream step.

## Geometry protocol and exceptions

ETKDGv3, explicit hydrogens, fixed seed 42, no RMS pruning, enforced specified chirality, and MMFF94 relaxation up to 500 iterations follow the Atom-JEPA procedure. Random coordinates are used only if normal embedding fails. Source atom order is mapped to a canonical molecular graph before adding H, and original-to-canonical correspondence is kept for outcome queries.

The portable `conformers.npz` archive contains **5,333 conformers** for every unique substrate; `geometry_metadata.jsonl` records methods, convergence and any rejected fallback conformers.

| Relaxation method | Unique structures |
|---|---:|
| MMFF94 | 517 |
| UFF_fallback | 13 |
| GFNFF_xTB_fallback | 3 |
| GFN2_xTB_fallback | 1 |

For 13 molecules without MMFF parameters, a fully parameterized UFF fallback is explicit. Three SF5 substrates use xTB GFN-FF relaxation; fallback conformers that break an input-bond distance check are discarded without using outcome labels. They retain 6, 8 and 9 valid conformers. One selenium substrate uses GFN2-xTB relaxation. These are recorded departures from the paper's main MMFF94 procedure and can be excluded in a method-sensitivity analysis.

MMFF nonconvergence is recorded, and last coordinates are retained as in the upstream package. 1 source record(s) use a nonconverged first geometry; the metadata allow excluding them consistently across every feature block.

## Compute scope and comparison policy

The completed bank uses the same **first retained conformer** for both encoders and electronic descriptors. This CPU-only first pass is not a ten-conformer ensemble. All generated conformers are saved, so `feature_conformers` can later be raised to 10 to produce matched ensemble means; molecules with fewer valid fallback conformers use all available ones. The configuration change produces a separate run key.

Each symmetry group is queried once. Member-atom descriptors are averaged into its single vector; group size is not a model input, weight, probability multiplier or extra training example. No outcome-product structure, positive-site label or measured yield is an input to the extractors. All feature arrays are unscaled; fit scalers only within training folds. Features currently describe the substrate, not the electrophile/condition identity.

## Files and rerun

Feature directory: `data/features/d76c42313dcf5e17ac96/`.

- `queries.csv`: aligned outcome-query identifiers, canonical/source atom correspondence, annotations and geometry QC.
- `atom_jepa.npz`, `random_backbone.npz`, `morgan.npz`, `electronic.npz`, `jepa_plus_qm.npz`: `X`, `present`, and `query_ids`, readable with `numpy.load(..., allow_pickle=False)`.
- `conformers.npz`, `geometry_metadata.jsonl`: portable coordinates and metadata.
- `manifest.json`: config, checkpoint revision/checksums, model initialization checksums, block hashes, row coverage and source-code fingerprints.

Install `requirements-features.txt`, run `scripts/setup_feature_tools.py`, then `scripts/extract_features.py --phase all`. The script resumes cached work. Run `scripts/summarize_features.py` to verify the exported bank. Large model weights, native binaries and per-calculation scratch logs remain in ignored `.cache/`; the compact feature/geometry exports and provenance are retained in the repository.

Sources: [Atom-JEPA paper, Appendix H](https://arxiv.org/html/2610.08400v1), [reference conformer implementation](https://github.com/khelverskovp/atom-jepa/blob/main/atom_jepa/conformers.py), [xTB properties](https://xtb-docs.readthedocs.io/en/latest/properties.html), [xTB Fukui calculations](https://xtb-docs.readthedocs.io/en/latest/sp.html).
