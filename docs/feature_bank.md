# Feature-bank dictionary

`data/features/latest.json` points to the exported run manifest. Each run directory contains the same outcome-query order in all five blocks. No scaler, prediction head, or supervised split is fitted during extraction.

## Queries and labels

`queries.csv` has one row per source record and symmetrically distinct aromatic C-H substitution outcome. `query_id` identifies that row; `structure_position_key` identifies the canonical substrate and its candidate atom group across duplicate source records. `source_atom_indices` follow the authoritative zero-based source SMILES. `canonical_atom_indices` follow the canonical SMILES in `geometry_metadata.jsonl`, before adding hydrogens.

`reported_positive` means that the source annotated at least one member atom of the outcome group. It is recorded-site membership, not a measured selectivity probability. Unannotated outcomes are not experimentally verified failures. `n_reported_outcomes` and `n_candidate_outcomes` preserve multi-positive and no-choice cases for appropriate downstream treatment. These labels are metadata and never enter the feature extractor.

`cohorts` is a JSON list of cumulative stages. It does not assign a train or test fold. Split all source records, positions and conformers of a substrate together using `structure_key`; keep source/publication groups together where the evaluation design requires it.

## Feature arrays

Every block is a pickle-free NumPy `.npz` archive containing `X` (float32, query rows by feature columns), `present` (one boolean per row), and `query_ids` (strings). Missing features, if a future run is partial, remain NaN with `present=False`. Never remove only an unavailable candidate from an otherwise retained reaction. The loader rejects incomplete selected cohorts.

| Block | Column order |
|---|---|
| `atom_jepa` | 256 final-layer l=0 scalar channels, 256 l=1 per-channel norms, 256 l=2 per-channel norms |
| `random_backbone` | Same extraction and architecture as `atom_jepa`, with upstream random initialization and frozen weights |
| `morgan` | Three separate 2048-bit atom-centered Morgan blocks: radius 1, radius 2, radius 3, with chirality |
| `electronic` | GFN1-xTB Mulliken charge, CM5 charge, f(+), f(-), f(0), total Wiberg bond order |
| `jepa_plus_qm` | Literal `atom_jepa` columns followed by `electronic` columns |

The electronic charges are in elementary charge units; Fukui indices and Wiberg orders are dimensionless. Descriptor values come from xTB's printed atomic tables, so their native print precision applies. GBSA methanol and the default 300 K electronic temperature provide a standardized descriptor environment. They do not reconstruct literature conditions. Current blocks describe the substrate only, without reagent or reaction-condition inputs.

For each position group, arithmetic averaging gives one vector over equivalent member atoms. Group size is not an added feature, training weight or probability multiplier. The first bank uses the same first retained conformer for the two encoders and the electronic descriptors. Averaging across additional conformers is a separate, configurable experiment.

Use the raw features with scalers fitted on training folds only. The large numerical and dimensional differences between blocks make global pre-split scaling or an unregularized high-dimensional head unsuitable comparisons.

## Portable geometry

`conformers.npz` stores all generated valid conformers, including explicit hydrogens, without Python objects:

- `structure_keys[s]` matches one line of `geometry_metadata.jsonl`.
- `atomic_numbers[atom_offsets[s]:atom_offsets[s+1]]` is the atom-number vector for substrate `s`.
- The first `n_original_atoms[s]` atoms follow its canonical SMILES; subsequently added hydrogens follow them.
- `conformer_structure_index[c]` gives the substrate index for conformer `c`.
- `positions[conformer_atom_offsets[c]:conformer_atom_offsets[c+1]]` gives its `[atoms, 3]` coordinates in Angstrom.
- `forcefield_status[c]` is zero for converged relaxation; nonzero retained MMFF/UFF statuses are explicitly flagged.
- `forcefield_energy_kcal_mol[c]` records the method's energy in kcal/mol. Energies from different methods are not comparable and are not input features.

Within each substrate, conformer order is embedding order after any rejected fallback conformers. The first bank uses the first retained one; it does not select conformers using labels or minimum energy. Geometry method, retained/rejected fallback indices, unspecified stereocenter counts and convergence summaries are in the metadata. `selected_geometry_converged` in the query table allows consistently filtering all blocks together.

## Provenance and comparison limits

The manifest records config, pinned checkpoint revision and file hashes, pretrained/random parameter checksums, frozen parameter counts, block hashes, source-manifest hash and implementation fingerprints. Full dependency versions are in `requirements-features.lock.txt`. Weights and executable binaries remain in the ignored cache; the compact features and exact coordinates are committed.

The random backbone controls architecture and extraction, not the pretraining corpus or training compute. Comparing it with the pretrained encoder does not isolate the JEPA objective. Audit pretraining overlap before making generalization claims. Feature extraction itself does not demonstrate predictive benefit; that requires matched downstream heads on grouped held-out data.
