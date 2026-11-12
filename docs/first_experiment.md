# First experiment: frozen JEPA features for EAS site ranking

## Question

Do fixed Atom-JEPA molecular atom features improve prediction of reported EAS sites relative to local structural and electronic descriptors, especially with limited labels and unfamiliar substrates?

This first study tests the utility of the published representation. It does not isolate the JEPA loss from the effects of backbone size, pretraining corpus or compute.

## Data gate

Follow the approved cumulative cohorts: conservative bromination (299 records), conservative halogenation including Cl/I (374), broad single-outcome (495), then all source records (535). Audit duplicate/source grouping and pretrained-data overlap. The latter two stages retain the source ambiguity flags, and the final stage contains multi-positive/no-choice examples requiring separate metrics. All conformers and duplicate records of one substrate must remain in one split.

Fresh 3D conformers and the first shared-conformer feature bank have been extracted from authoritative SMILES with source-to-canonical atom correspondence. See `reports/feature_extraction.md` for exact coverage, protocol exceptions and frozen model provenance. Actual supervised train/validation/test assignments remain downstream of feature extraction.

## Compared representations

- Published structure-only Atom-JEPA `molecules` checkpoint, frozen; start with invariant per-atom features.
- Randomly initialized frozen version of the same backbone, same features/head.
- Atom-centered Morgan fingerprints and a simple scorer.
- Electronic atomic descriptors (e.g. xTB CM5 charges or established predicted QM features) with a suitable simple scorer.
- JEPA + QM feature concatenation to test complementarity.

Choose a linear site scorer and one small nonlinear head with hyperparameters selected only on inner training/validation groups. Each symmetrically distinct position already has one feature vector; rank those queries directly. Do not provide the known product to the input feature extractor. Fit all scalers only within the training fold; the exported features are unscaled.

## Evaluation

Use grouped outer cross-validation and inner-only model selection. The implemented linear pilot uses identical nested 50- and 100-label training subsets plus the full available training fold, with three subset seeds; 200 is omitted because full training contains only 239-240 records. Label budgets include inner-selection labels. Report actual cohort/fold counts and compare paired molecular predictions. The linear optimizer is deterministic and is not given artificial initialization repetitions.

See `docs/evaluation_protocol.md` and `reports/evaluation_preflight.md` for the reviewable source/series graph, explicit whole-molecule similarity constraints and split proposals. No production splits are frozen and no real-data scorers have been trained yet. Nonlinear heads remain a later comparison.

Report outcome-level top-1, top-2 where meaningful, reciprocal rank, subgroup results, and training-size learning curves. Use molecule/group-level uncertainty, not atom-level confidence intervals. The conservative cohort's analytical uniform-outcome baseline is 28.3%, not 50%.

## Follow-up gate

Reproducible low-label and structurally separated-fold gains motivate matched pretraining on the same unlabeled corpus/backbone with JEPA versus reconstruction/denoising. Confirm gains on an independent dataset before extending the claim to other reaction families. A small improvement on this pilot alone remains exploratory.

No new mechanism or selectivity annotations are required for the first frozen-feature study. Any manually reviewed corrections must be stored as overrides with original evidence and never applied to the raw snapshot.

