# Benchmark data dictionary

## Units and provenance

One molecule record corresponds to one row/ID in the original source, not necessarily one independent experiment. Several free-text annotations pool reagents, references or products. There are 535 rows but 534 distinct isomeric substrate graphs. Candidate atoms and hypothetical products are nested within molecule records; they are not independent observations.

Raw files are copied without newline, whitespace, index or chemical-standardization changes. The files named `compounds_smiles.csv` and `compounds_yield.csv` are actually whitespace-delimited records; `compounds_doi.csv` is comma-delimited. Parsing uses Unicode whitespace, so the nonbreaking space in comp4 is handled without changing its raw file. Trailing commas in label lists are accepted. No product ratios are inferred from isolated yields.

## Molecule table / JSONL

- `compound_id`, `source_line`: stable upstream ID and one-based line position.
- `source_smiles`, `source_label_text`, `reference_raw`, `yield_reagent_note_raw`: original source text, with only outer record delimiters removed.
- `observed_atom_indices`: **zero-based** indices in `Chem.MolFromSmiles(source_smiles)` before any standardization or hydrogen addition.
- `atom_mapped_smiles`: equivalent mapped structure; map number = original zero-based index + 1. Re-parsing a canonical SMILES changes atom order: recover correspondence from map numbers, never reuse the old indices blindly.
- `canonical_isomeric_smiles`: identity key for deduplication, not a replacement for source indexing.
- `constitution_key`: identity after removal of specified stereochemistry, preserving isotope labels.
- `tautomer_key`: RDKit canonical tautomer identity, used only to help prevent split leakage. It is not used to change input structures or positive-site assignments.
- `same_structure_ids`, `same_constitution_ids`, `same_canonical_tautomer_ids`: record IDs sharing the corresponding identity.
- `suggested_identity_group`: stable digest of the canonical tautomer key; a grouping aid, not a finalized split.
- `murcko_scaffold`: graph scaffold descriptor. Many arenes share this scaffold; it is not sufficient as the only split criterion.
- `formula`, `molecular_weight`, `n_heavy_atoms`, `formal_charge`, `n_fragments`: graph-derived molecular summary; weight is in Da.
- `n_assigned_chiral_centers`, `n_unassigned_chiral_centers`: tetrahedral-center checks; do not enumerate all stereochemical phenomena.
- `molecule_family`: aromatic-ring coverage: carbocyclic only, heteroaromatic only, or mixed. Nonaromatic heteroatoms do not determine this class.
- `aromatic_ring_signatures`, `observed_ring_signatures`: ring-basis sizes and atom-element counts, e.g. `6:C5N1`. These are descriptive graph signatures, not complete named ring-system classifications.
- `candidate_atom_indices`: all aromatic carbons bearing H in the supplied graph.
- `candidate_outcome_groups`: lists of candidate indices whose hypothetical single-Br substitutions give the same canonical isomeric product. Group indices are zero-based; original atom indices are retained.
- `constitution_only_outcome_groups`: diagnostic partition after removing specified stereochemistry. This does not govern the default scoring.
- `observed_outcome_group_indices`: candidate outcome groups touched by the source labels.
- `expanded_equivalent_observed_atoms`: union of all members of observed groups. This derived set avoids penalizing an unlabeled equivalent counterpart; it does not overwrite source labels.
- `n_candidate_atoms`, `n_candidate_outcomes`, `n_observed_atoms`, `n_observed_outcomes`: counts before/after grouping.
- `uniform_outcome_top1_chance`: number of annotated distinct outcomes divided by candidate distinct outcomes; an analytical random-choice reference, not a measured model score.
- `has_equivalent_candidate_atoms`, `has_equivalent_unannotated_positive`: symmetry-audit flags.
- `symmetry_rank_product_agree`: independent partition cross-check against RDKit canonical atom symmetry ranks.
- `stereo_changes_outcome_partition`: whether removal of specified stereo alters this molecule's partition.
- `sdf_same_structure_by_row`, `sdf_same_atom_order_by_row`: checks against the supplementary SDF at the same row position. SDF records have no IDs, so these are verified row matches, not an independent identifier join.
- `sdf_has_3d_coordinates`: SDF conformer dimensionality; all supplied conformers are 2D.
- `validation_errors`: structural/index failures; blank references are metadata review flags instead.
- `review_flags`: record-level reasons requiring care, possibly overlapping.

## Reaction annotations

- `reaction_class_evidence`: reaction-class names inferred from recognizable reagent strings in the raw notes.
- `reaction_class`: exclusive evidence bucket. `multiple_EAS_classes_evidenced` means several recognized classes occur in one record. `unresolved_reagent` means no supported assignment; this is not silently called bromination.
- `ambiguous_reagent_typography`: strings such as `12`, `C12`, `ICI`, `Nal` whose intended identity requires checking the original source.
- `multiple_yield_values`: more than one percent sign; a cautious pooled-experiment/mixture warning, not proof of distinct conditions.
- `explicit_repeated_substitution_note`, `explicit_mixture_note`: heuristic free-text review flags. False negatives/positives are possible; they are screening aids.
- `yield_values_percent`: numerical percentages copied from the note without assigning them to products, conditions or ee. Do not normalize these into selectivity ratios.

## Cohorts

- `eligible_broad_single_outcome`: valid structure/site labels; at least two distinct candidate outcomes; one source-annotated outcome after symmetry grouping. May contain provenance or reaction-condition ambiguities.
- `eligible_conservative_bromination`: the broad criterion plus only bromination evidence, nonempty source reference, and no review flags. This is a conservative pilot cohort, not a guarantee of fully verified monobromination data.

## Candidate atom table

`candidate_sites.csv` contains the source ID, original atom index, outcome-group index and size, raw-positive and symmetry-expanded-positive flags, plus `hypothetical_Br_marker_product`. The marker product is a synthetic bookkeeping structure used to compare outcomes. It must not be passed off as an experimentally observed product or used as evidence for a reaction mechanism.

