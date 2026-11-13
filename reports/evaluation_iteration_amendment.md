# Approved iteration-budget amendment

Original approved proposal: `09ebadea135ec50dc873`. It stopped at the 1,000-iteration cap during inner JEPA+QM fitting after saving 146 of 175 final fits.

Approved replacement: **`af05d8348e60ad38b6e4`**. Status: **user-approved and frozen**. The exact approval is "Approve cap 5,000 and rerun". Its proposal-hash binding is in `data/splits/af05d8348e60ad38b6e4/approval.json`; current execution state is in `reports/evaluation_amended_execution_status.json`.

The only configuration change is `optimizer.max_iterations`: **1,000 → 5,000**. Both the 1e-7 optimizer target and the successful-termination/accepted-gradient checks remain unchanged. Optimizer type, initial coefficients, objective, feature blocks, transformations, hyperparameter grids, seeds, data, leakage groups, outer/inner folds and nested label budgets are identical. The larger cap permits additional iterations; it does not force already converged fits to continue.

The execution record attributes the change to the convergence failure rather than outer-test performance. This is an execution account, not an independently verifiable access history: the original run had already written held-out metrics for 146 fits. The completed reporting audit verifies that all 146 original model archives and fit JSON files are byte-identical in the amended run. This establishes unchanged outputs, not what was viewed or how the cap was chosen. See `evaluation_reporting_audit.md` and `.json`. The amended run executed the same 6,300 inner fits and 175 final fits in a new directory, without resuming or overwriting the original attempt.

The original-cap replay reproduced the failure exactly: inner fold 1, 160 training molecules, JEPA penalty 0.0001 and QM penalty 10. Its diagnostics are in `evaluation_convergence_diagnosis.md` and `.json`. The artifact audit in `evaluation_failed_attempt_integrity.md` and `.json` independently checked the cap-only difference, identical groups/folds/locks and the saved model hashes. These checks do not prove the larger cap will converge.

The subsequently approved fresh run **completed all 6,300 inner fits and 175 final fits** and passed saved-model replay validation. Nine fits needed more than 1,000 iterations; the largest used 1,186. The original outputs remain unchanged and separate. Current results: `evaluation_results.md` and `.json`; validation: `evaluation_result_validation.md` and `.json`.

Claude review prompt: `reports/claude_iteration_amendment_review_prompt.txt`.

Files: `configs/evaluation_iteration_extension.json`, `data/splits/af05d8348e60ad38b6e4/proposal.json` and its sibling `grouping.json`.

## Not done

The supplied results review has been received and its reporting findings incorporated; focused re-review and unchecked independent-review gates remain outstanding. No convergence threshold was relaxed and no failed fit was counted as valid.
