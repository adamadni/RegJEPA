# Approved iteration-budget amendment

Original approved proposal: `09ebadea135ec50dc873`. It stopped at the 1,000-iteration cap during inner JEPA+QM fitting after saving 146 of 175 final fits.

Approved replacement: **`af05d8348e60ad38b6e4`**. Status: **user-approved and frozen**. The exact approval is "Approve cap 5,000 and rerun". Its proposal-hash binding is in `data/splits/af05d8348e60ad38b6e4/approval.json`; current execution state is in `reports/evaluation_amended_execution_status.json`.

The only configuration change is `optimizer.max_iterations`: **1,000 → 5,000**. Both the 1e-7 optimizer target and the successful-termination/accepted-gradient checks remain unchanged. Optimizer type, initial coefficients, objective, feature blocks, transformations, hyperparameter grids, seeds, data, leakage groups, outer/inner folds and nested label budgets are identical. The larger cap permits additional iterations; it does not force already converged fits to continue.

The change is driven by the convergence failure, without inspecting outer-test performance. A successful full run would execute the same 6,300 inner fits and 175 final fits in a new output directory. The existing runner has no resume path, so it will rerun all five blocks and preserve the original incomplete attempt separately.

The original-cap replay reproduced the failure exactly: inner fold 1, 160 training molecules, JEPA penalty 0.0001 and QM penalty 10. Its diagnostics are in `evaluation_convergence_diagnosis.md` and `.json`. The artifact audit in `evaluation_failed_attempt_integrity.md` and `.json` independently checked the cap-only difference, identical groups/folds/locks and the saved model hashes. These checks do not prove the larger cap will converge.

Claude review prompt: `reports/claude_iteration_amendment_review_prompt.txt`.

Files: `configs/evaluation_iteration_extension.json`, `data/splits/af05d8348e60ad38b6e4/proposal.json` and its sibling `grouping.json`.

## Not done

The fresh rerun and its verification are pending. No convergence threshold has been relaxed and no failed fit has been counted as valid.
