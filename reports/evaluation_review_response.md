# Response to the pre-experiment Claude review

Reviewed baseline: `2524f802de4b93b82cb00bafd9c9770a516f62fc`.
Corrected proposal: `09ebadea135ec50dc873` (`data/splits/09ebadea135ec50dc873/proposal.json`). Production splits remain unfrozen; no real-data scorer or full experiment has run.

## Confirmed and corrected

- **Freeze bypass:** every frozen load now checks the exact sibling proposal hash, requires its grouping audit, and regenerates the canonical graph and full schedule. Configured outer/inner folds, all budget/seed settings and the full-training setting are mandatory. Regression tests reject status-only edits, rehashed dissolved groups with a forged matching audit, omitted settings and corrupt bindings. This enforces a reproducible artifact contract, not proof that a human approved it.
- **DOI semicolon:** comp504 now retains its complete DOI and no longer creates the spurious `d` citation key. Unresolved-record count changes from 38 to 37; raw annotations are unchanged.
- **Windows exact ties:** fixed-coefficient tests reproduced BLAS differences up to 4.44e-16 for identical Morgan rows. The objective, predictions and component scores now use a fixed per-row float64 reduction; exact feature ties stay tied without a score tolerance or merging candidate rows.
- **Reporting:** final fits flag grid boundaries and record JEPA/QM RMS contributions and their mean product. Relative RMS is not unique information attribution. Markdown includes every descriptive metric and primary top-1 paired intervals. The grid and statistical selection procedure are unchanged.
- **Fresh-clone tests:** README and reviewer instructions explicitly create the temporary directory parent before invoking pytest.

## Recomputed checks

- **90 tests pass**, including all reported freeze probes, synthetic fits, gradients, candidate-order invariance and platform-specific identical-input ties.
- **299 source records, 1,264 rows, 250 leakage groups** remain.
- Group assignments and the **complete outer/inner/budget/subset schedule are identical** to the previous proposal, as is the fitting configuration. The proposal ID changes because locked implementation files and the grouping audit change.
- All tracked raw-data and feature-bank files are byte-identical to the reviewed baseline.
- Median nearest outer-training Tanimoto is **0.344828**; **27/299** test molecules have a training neighbour at >=0.5. This does not establish an OOD benchmark.

## Scientific limitations retained

comp480 and comp481 remain in outer folds 1 and 0. Their papers are consecutive parts of a broader Furopyridines series; the [publisher's XXVIII page](https://onlinelibrary.wiley.com/doi/10.1002/jhet.5570350602/full) identifies the preceding part. The current algorithm holds out explicit source/series tokens and similarity components, not all independently recovered publication families or broader research series. Adding curated series edges would be a new documented grouping protocol, rather than a silent change during this correction.

A common numerical L2 grid does not imply matched effective shrinkage across representations. Grid-edge flags expose bounded-search concerns; finite penalties also do not guarantee that a fitted concatenation reproduces a single-block fit exactly. No grid extension is made after inspecting outer-test performance.

Real-data convergence and automatic resume remain untested/unimplemented. Paired intervals remain limited to the primary top-1 metric and are conditional on fitted models. Pointer ranges are not guessed. Pretraining overlap, one random encoder seed, and incomplete minor-product annotations remain limitations.

## Not done

No production freeze, real-data scorer fit, full experiment, independent-series reconstruction or pretraining-overlap audit. The updated `reports/claude_review_prompt.txt` is ready for a focused re-review.
