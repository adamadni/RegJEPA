# RegJEPA initial EAS data audit

Audit scope: pinned public sources; RDKit 2026.09.2. This is a structural and annotation audit, not independent re-extraction of every literature experiment.

## Decision

Suitable for a frozen-feature feasibility experiment. Insufficient alone for training a JEPA from scratch, establishing superiority across reaction classes, or modeling quantitative product ratios. Use the conservative bromination cohort initially and keep the other records visible as review and sensitivity cohorts.

## Acquired data and integrity

- 535 source rows, 535 IDs, 534 distinct isomeric structures.
- All 535 substrate SMILES parse and every observed index is an aromatic carbon bearing H under the zero-based source convention.
- 0 structural/label validation failures.
- 69 rows have an empty reference field upstream. They are retained and flagged, never forward-filled.
- The SDF contains 535 records: 534 match the source graph/stereochemistry and atom order. comp485 has additional alkene stereochemistry in the SDF. The SMILES source remains authoritative.
- All SDF coordinates are 2D; these are not suitable 3D Atom-JEPA inputs without conformer generation.
- Raw files, upstream licenses, commit IDs, Git blob hashes and SHA256 digests are preserved in `data/raw/` and `data/sources.lock.json`.

## Reaction evidence from the original free-text notes

Although the RegioSQM20 paper describes a 535-molecule bromination benchmark, the original annotation notes contain other electrophiles and multiple experiments. The following counts are heuristic text classifications; reaction conditions and product-to-experiment assignments have not been verified from the cited sources.

| Note-derived category | Molecule records |
|---|---:|
| bromination_only_evidence | 390 |
| chlorination_only_evidence | 41 |
| iodination_only_evidence | 65 |
| multiple_EAS_classes_evidenced | 30 |
| nitration_only_evidence | 1 |
| unresolved_reagent | 8 |

The 30 mixed-class records belong to one category each; the table does not count them repeatedly in the single-class categories. Ambiguous strings such as `12`, `C12`, `ICI` and `Nal` are flagged rather than silently corrected. Recognized reagent evidence can coexist with ambiguous typography. This source has no enolization, 1,2/1,4 addition or cycloaddition benchmark.

## Molecular coverage

| Graph-derived family | Records |
|---|---:|
| carbocyclic_aromatic_only | 84 |
| heteroaromatic_only | 262 |
| mixed_carbo_and_heteroaromatic | 189 |

Family refers to the aromatic rings, not to whether the rest of the molecule contains heteroatoms. The mixed class has both carbocyclic aromatic and heteroaromatic rings. Ring signatures below are RDKit ring-basis descriptors rather than named chemical series; fused systems can contribute more than one signature.

| Ring containing an annotated site (size:element counts) | Records containing that signature |
|---|---:|
| 5:C3N2 | 96 |
| 5:C3N1O1 | 8 |
| 5:C3N1S1 | 14 |
| 6:C5N1 | 96 |
| 6:C4N2 | 78 |
| 5:C4N1 | 60 |
| 5:C2N3 | 6 |
| 6:C6 | 115 |
| 5:C4S1 | 7 |
| 5:C4O1 | 57 |

Heavy atoms: {'min': 5, 'median': 15, 'max': 44}; molecular weight: {'min': 68.079, 'median': 217.276, 'max': 661.505} Da. There are 19 records with assigned tetrahedral stereocenters and 31 with unassigned potential centers. This does not count all forms of stereochemistry or establish stereoselectivity labels.

## Symmetry findings

- 2683 eligible aromatic C-H atoms reduce to 2290 distinct single-substitution outcomes across records.
- 166 molecules contain at least one symmetry-equivalent pair of candidate atoms.
- 16 molecules have at least one equivalent counterpart of an observed site omitted from the source labels.
- 64 records list multiple atom labels, but 26 of those collapse to one distinct marker outcome. The remaining 38 records have multiple distinct annotated outcomes and need experiment/mixture/multiple-substitution review.
- 2 molecules have no nontrivial site choice: comp485, comp588. Exclude them from the primary ranking metric and report them separately.
- Graph canonical symmetry ranks and explicit marker-product partitions agree for every record. Removing specified stereochemistry changes no within-molecule partition in this dataset; retain stereo-aware handling for later datasets.
- comp13 and comp598 are the same substrate and the same marked substitution outcome under different halogenation annotations. Keep them in one split group; do not treat them as independent chemistry.

### Recommended policy

Define equivalence by identical hypothetical single-Br marker products on the input graph, preserving its tautomer, isotope labels and specified stereochemistry. No 3D distance tolerance is used. Green source labels and blue equivalent counterparts are visible in `molecule_browser.html` and `symmetry_examples.png`.

Retain every raw atom label. Add equivalence groups and expanded positive sets in derived files. Average site features or logits within each group, then normalize/rank across distinct groups. This avoids giving a group with two symmetric atoms an automatic counting advantage. A degeneracy/rate correction is a different physical modeling question and requires quantitative kinetic data. Keep the molecule as one evaluation unit.

For a genuinely multisite observation, do not assume equal probabilities or infer that all labels describe a mixture in one experiment. Those rows remain a separate review cohort. The product marker is a symmetry test, not a claim that every row is experimentally brominated.

## Proposed starting cohorts (no train/test split frozen yet)

| Cohort | Records | Unique structures | Records with equivalent candidates |
|---|---:|---:|---:|
| All source records | 535 | 534 | 166 |
| Broad single-outcome | 495 | 494 | 155 |
| Conservative bromination | 299 | 299 | 111 |

Broad single-outcome requires valid structure/site labels, at least two distinct candidate outcomes, and exactly one annotated outcome after symmetry reduction; it can include mixed reaction evidence, missing references and repeated-experiment notes, so it is a sensitivity cohort. Conservative bromination additionally requires only bromination evidence, a nonempty reference, no multiple-yield, explicit mixture, repeated-substitution, ambiguous-reagent, SDF-mismatch or symmetry-disagreement flag. These are screening rules, not a chemical revalidation of all experiments. Further manual review can reduce or enlarge the cohort.

Uniform random choice among distinct candidate outcomes yields a mean top-1 hit rate of 28.3% on the conservative cohort (computed from candidate counts; not a trained baseline).

## Is this large/high-quality enough?

The data are traceable, literature-derived, and reusable with an upstream MIT license, and the source atom indices pass structural checks. That is enough to start a small frozen-encoder comparison with local fingerprints and QM descriptors. It is not a uniform-condition, independently replicated, mixture-resolved set of 535 experiments. Metadata curation is needed before interpreting learned condition effects or product ratios.

The conservative cohort has 299 molecules; its candidate atoms are not independent additional training examples. At a hypothetical 80% accuracy over 299 independent predictions, an illustrative 95% Wilson interval is 75.1%–84.1%. A single small held-out test is less precise, and structural/source correlation reduces effective sample size further. Use grouped outer cross-validation, inner-only model selection, low-label learning curves and paired molecule/group resampling. A gain of a few cases is exploratory, not a broad result.

Keep identical structures, tautomer variants and records from the same experiment/source together. Canonical tautomer IDs are supplied as grouping aids, not as substituted feature inputs. Resolve multi-reference strings and group publication/patent families before freezing splits. The large benzene core makes a simplistic scaffold-only split inadequate; evaluate fingerprint-cluster splits and report family coverage. Audit pretrained corpus overlap and label exposure before benchmarking published checkpoints.

## Files and reproduction

- `data/processed/molecules.csv`: one row per source molecule, raw annotations, descriptor summary, groups, flags and cohort inclusion.
- `data/processed/molecules.jsonl`: the same records with arrays and booleans preserved.
- `data/processed/candidate_sites.csv`: one row per candidate atom, with group index, source-positive and expanded-positive flags, and hypothetical marker product.
- `data/processed/summary.json`: machine-readable counts.
- `molecule_browser.html`: offline searchable structure/annotation browser covering all 535 records.
- `symmetry_examples.png`: selected examples with source atom indices.
- `docs/data_dictionary.md` and `docs/benchmark_decisions.md`: semantics, limitations, and pending decisions.

Run `python scripts/download_data.py --verify-only`, `python scripts/validate_benchmark.py`, `python scripts/build_reports.py`, and `python -m pytest -q` in the configured environment. No encoder training, 3D conformer generation or experimental relabeling has been performed.

Sources: [original EAS paper](https://doi.org/10.1039/C7SC04156J), [RegioSQM20](https://doi.org/10.1186/s13321-021-00490-7), [original data](https://github.com/jensengroup/db-regioselectivity), [RegioSQM20 inputs](https://github.com/jensengroup/SI_RegioSQM20).
