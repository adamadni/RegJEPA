# Evaluation preflight: grouping and split proposal

No real-data scorer has been fitted. Splits are proposals, not frozen experiment assignments.

Proposal: `83cde96ea4d814325d8c`. File: `data/splits/83cde96ea4d814325d8c/proposal.json`.

Cohort: 299 records / 299 structures / 1264 outcome queries.
Uniform expected top-1: 0.2831475288; first-row tie artifact: 0.3076923077.

Identity + explicit source/series + whole-molecule similarity components: **250**; largest selected component: **5** records.
Maximum Tanimoto between different components: **0.566667**, below the fixed 0.6 boundary.

All 535 source records participate in linkage before cohort filtering. Pointer tokens and named targets are joined explicitly, including mixed DOI/patent + pointer records. Suffix aliases and international patent families are not guessed.

Unverified bibliographic identities remain for **38** selected records. They are listed in the JSON sibling; known explicit links are constrained, but this is not a verified publication-family holdout.

| Outer fold | Training records | Test records |
|---|---:|---:|
| 0 | 240 | 59 |
| 1 | 239 | 60 |
| 2 | 239 | 60 |
| 3 | 239 | 60 |
| 4 | 239 | 60 |

Candidate budgets are nested within each subset seed and contain whole leakage components. All tuning labels are inside the stated budget; full training is fitted once per fold.
Actual label counts: `[('100', 100), ('50', 50), ('full', 239), ('full', 240)]`.

Tie handling: exact ties receive uniform expected top-1/top-2 and expected reciprocal rank. No raw-atom helper is used on bank rows.
Linear scaling: subtract each molecule's candidate-feature mean, then scale dense columns using equal-molecule, training-only within-molecule RMS. Morgan columns remain unscaled. Columns without training variation are dropped. No intercept is fitted.
JEPA + QM has independently tuned block penalties (36 combinations); each other block has 6 candidates. This larger search is reported and contained inside inner folds.
Planned workload: 6300 inner fits plus 175 final fits. These have not been run.

The primary cohort retains all 299 records. Seven geometry-fallback cases are flagged, not removed to improve a metric. No geometry sensitivity model is fitted in this preflight.

Read docs/evaluation_protocol.md and reports/claude_review_prompt.txt before freezing or running the experiment.
