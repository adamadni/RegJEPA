# First experiment: frozen JEPA features for EAS site ranking

## Question

Do fixed Atom-JEPA molecular atom features improve prediction of reported EAS sites relative to local structural and electronic descriptors, especially with limited labels and unfamiliar substrates?

This first study tests the utility of the published representation. It does not isolate the JEPA loss from the effects of backbone size, pretraining corpus or compute.

## Data gate

Use the proposed conservative bromination cohort after inspecting flags and molecule depictions. Audit duplicate/source grouping and pretrained-data overlap. Keep the broad single-outcome cohort as a separately disclosed sensitivity analysis. Generate fresh 3D conformers from authoritative SMILES while tracking atom-map correspondence; the supplied SDF is 2D. All conformers of one substrate must remain in one split.

## Compared representations

- Published structure-only Atom-JEPA `molecules` checkpoint, frozen; start with invariant per-atom features.
- Randomly initialized frozen version of the same backbone, same features/head.
- Atom-centered Morgan fingerprints and a simple scorer.
- Electronic atomic descriptors (e.g. xTB CM5 charges or established predicted QM features) with a suitable simple scorer.
- Optional JEPA + QM feature concatenation to test complementarity.

Choose a linear site scorer and one small nonlinear head with hyperparameters selected only on inner training/validation groups. Average site scores over equivalent candidate groups and rank the resulting outcomes. Do not provide the known product to the input feature extractor.

## Evaluation

Use grouped outer cross-validation and inner-only model selection. Evaluate identical 50-, 100-, and 200-label training subsets where fold size allows, plus the full available training fold. Report the actual cohort/fold counts; do not promise training sizes that a group split cannot provide. Repeat training seeds and compare paired molecular predictions.

Report outcome-level top-1, top-2 where meaningful, reciprocal rank, subgroup results, and training-size learning curves. Use molecule/group-level uncertainty, not atom-level confidence intervals. The conservative cohort's analytical uniform-outcome baseline is 28.3%, not 50%.

## Follow-up gate

Reproducible low-label and structurally separated-fold gains motivate matched pretraining on the same unlabeled corpus/backbone with JEPA versus reconstruction/denoising. Confirm gains on an independent dataset before extending the claim to other reaction families. A small improvement on this pilot alone remains exploratory.

No new mechanism or selectivity annotations are required for the first frozen-feature study. Any manually reviewed corrections must be stored as overrides with original evidence and never applied to the raw snapshot.

