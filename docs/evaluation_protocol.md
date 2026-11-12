# Linear frozen-feature evaluation protocol

This pilot tests whether fixed representations support prediction of source-reported EAS outcomes. It does not measure product ratios, infer mechanisms, isolate the JEPA objective, or establish performance after encoder fine-tuning. All current feature-bank and raw-data bytes remain unchanged.

## Scope and dependencies

`configs/evaluation.json` specifies the conservative bromination cohort, all five feature blocks, five outer folds, three inner folds, nested 50/100-label budgets with subset seeds 17/29/43, and one full-training setting per outer fold. No nonlinear head is included yet.

`requirements-evaluation.txt` explicitly pins NumPy 2.5.4 and SciPy 1.18.1 alongside the existing audit dependencies. SciPy supplies L-BFGS-B; no scikit-learn, pandas, plotting library or training framework is imported. The feature extra already used these NumPy/SciPy versions, so no new package installation was needed in the working environment.

## Grouping before split assignment

The grouping graph uses all 535 source records before selecting the 299-record cohort. This preserves links through records outside the cohort. It joins:

1. Identical canonical structures and canonical tautomer identities.
2. Explicit, case-normalized DOI/citation tokens and numbered patent identifiers. Patent kind codes are normalized, but different national applications are not assumed to belong to the same family.
3. Shared `comp...` series tokens and their explicitly named record targets, including mixed patent/DOI + series strings. Number ranges are not guessed; named series membership is preserved as supplied. Suffix aliases such as `comp222c` are not silently changed to `comp222`.
4. Every pair of whole-molecule radius-2, 2048-bit, chiral Morgan fingerprints with Tanimoto >=0.6. Connected components are used rather than centroid clustering, so even transitive links stay together and no above-threshold neighbour can cross folds.

These constraints are combined by graph union, not applied independently and then broken to balance folds. Group IDs are deterministic hashes of the complete source-record membership, including out-of-cohort members. The grouping audit retains edge reasons, explicit identifiers, unresolved citation keys and complete component memberships. It reads no observed positions or yields to decide grouping.

The raw metadata do not recover every publication or international patent family. Pointer-only and opaque-key records remain explicitly flagged; no source is repaired by guessing. The resulting evaluation is a holdout of the known source/series and structural components, not a proven independent-publication-family benchmark. Unverified source identities are a review item before freezing and a limitation to retain in any eventual report.

Largest-first assignment balances molecular counts across whole components, with fixed seeded tie-breaking. The same grouping is used for inner folds and whole-component learning-curve subsets. Requested fold counts are not silently reduced. A budget is an upper bound on source records, not candidate rows; actual counts are reported if atomic groups prevent an exact budget. Every supervised label used for inner selection is inside that budget.

## Linear readout and scaling

Each bank row is already one distinct candidate outcome. The shared score is `w.T @ x`, with no intercept; the encoder and descriptors remain fixed. Scores are normalized by softmax within each source record. Training minimizes mean molecular negative log likelihood plus `0.5 * sum(lambda_j * w_j**2)`. Unreported sites compete in this objective but are not claimed to be experimentally confirmed negative products.

For this linear-only protocol, subtract each molecule's mean candidate feature vector. This leaves its softmax probabilities unchanged for any fixed coefficients. Dense features then use per-column RMS of these centered values, fitted only on training molecules with equal total weight per molecule. Morgan features stay unscaled. Directions with training RMS <=1e-12 are omitted; validation/test sets cannot reactivate them or refit the scale. Prediction-time centering uses only the new molecule's candidate features, never its observed label. A future nonlinear head must reconsider centering because molecular context may matter through interactions.

The six L2 strengths are 1e-4, 1e-3, 1e-2, 0.1, 1 and 10, selected separately for each feature block using molecule-averaged inner-validation NLL. For JEPA + QM, independently search the Cartesian product of JEPA and QM strengths (36 candidates). The larger search budget is explicit and remains inside inner folds. A single penalty would still be a valid model; separate penalties make the tested complementarity less dependent on a shared shrinkage choice.

L-BFGS-B starts at zero with analytic gradients, at most 1000 iterations, gradient tolerance 1e-7, relative objective tolerance 1e-14 and accepted final gradient infinity norm <=1e-5. Failed convergence aborts a fit; it is never counted as a valid result. No early stopping uses outer-test data. An exact inner-NLL tie prefers stronger shrinkage in a fixed candidate order.

## Exact ties and evaluation units

The harness ranks outcome rows directly; it never calls `chemistry.mean_group_scores` or `chemistry.top1_hit` on bank rows. Those legacy helpers take atom-indexed scores and are outside this evaluation path.

If a reported outcome has `a` rows strictly above it and is tied with `t` rows including itself, top-k credit is `clip(k-a, 0, t)/t`. Reciprocal rank is the mean of `1/(a+1)` through `1/(a+t)`, not the reciprocal of their mean rank. Ties are exact equality in float64 logits; no implicit NumPy relative tolerance is used. Reordering candidate rows must not change any metric.

Report molecule-average top-1, top-2, top-2 restricted to candidate sets larger than two, MRR, raw NLL, uniform NLL `log(k)`, and information gain `log(k)-NLL`. Training loss is not divided by `log(k)`, which would change molecular weights. Repeated subset metrics are averaged per molecule before aggregate reporting; logits are not ensembled.

Paired 95% intervals use 10,000 bootstrap resamples of the full leakage components, preserving all member molecules and method pairing. The metric remains a molecular mean, and candidate rows are never resampled independently. These intervals are conditional on the fitted out-of-fold models, not a complete estimate of uncertainty from repeating model selection. Comparisons are exploratory and are not presented as adjusted confirmatory significance tests.

All 299 records stay in the primary result. Seven geometry-fallback records are flagged for a possible predefined secondary sensitivity analysis using matched cases across every block. No sensitivity models are implemented or run in this milestone.

## Preparation, review, explicit freeze, then run

```text
python scripts/prepare_evaluation.py
python scripts/run_linear_evaluation.py --splits data/splits/PROPOSAL_ID/proposal.json
```

Both commands above perform preparation/validation only. The first saves immutable `proposal.json` and `grouping.json`, plus `reports/evaluation_preflight.md` and its JSON sibling. Module/script source, raw-reference, processed-record, feature-block and manifest hashes are locked into the proposal. Only derived files under `data/processed/` canonicalize CRLF to LF for their input hashes, matching Git's text normalization without rewriting the files. Raw snapshots and feature exports always use exact bytes. Actual counts and source limitations are in the preflight. Review with `reports/claude_review_prompt.txt` before running the full experiment.

After that review, explicitly freeze the specific proposal and run:

```text
python scripts/prepare_evaluation.py --freeze data/splits/PROPOSAL_ID/proposal.json
python scripts/run_linear_evaluation.py --splits data/splits/PROPOSAL_ID/frozen.json --run
```

Freezing validates hashes, graph constraints and fold/subset coverage, and regenerates the plan from locked metadata to check reproducibility. It writes a separate immutable frozen artifact; it does not overwrite the proposal. Fitting rejects proposals, changed inputs and changed latest-bank pointers. A later protocol or implementation change requires preparing a new proposal. The runner also requires `--run`; merely inspecting an artifact never fits a model.

The full run writes per-fit coefficients, masks/scales, convergence diagnostics, hyperparameter searches and query-aligned held-out scores. Markdown/JSON summaries and paired intervals are stored under `reports/evaluation_runs/PROPOSAL_ID/`. Existing output directories are not overwritten. Per-fit artifacts remain inspectable after an interrupted run, but automatic resume is not implemented.

## Outstanding scientific checks

The Atom-JEPA paper identifies Uni-Mol molecular pretraining; benchmark membership and label exposure have not been audited. The random backbone uses only the existing encoder seed 42. No conclusion isolates objective choice from corpus, architecture and pretraining compute. Additional review of bibliographic identities, independent confirmation data and later multi-positive/no-choice tasks remains necessary before broader claims.

Implementation references: [SciPy L-BFGS-B](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html), [Atom-JEPA paper](https://arxiv.org/html/2610.08400v1). Structural grouping here deliberately uses all-pair connected components rather than [Butina centroid clusters](https://www.rdkit.org/docs/source/rdkit.ML.Cluster.Butina.html).
