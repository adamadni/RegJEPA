# Approved amended evaluation execution status

Proposal **`af05d8348e60ad38b6e4`** is user-approved, frozen and **complete with artifact validation passed**. Frozen approval/runtime artifacts were committed at `f938849` before fitting.

- Completed inner fits: **6,300 / 6,300**.
- Completed final fits: **175 / 175**; **35 for each of the five blocks**.
- All fits converged under the unchanged acceptance rule. Largest iteration count: **1186**; **9** fits exceeded the original cap.
- Every saved model hash, training-only transform and query score replayed, and all summaries and paired intervals reproduced without refitting. Maximum score, objective and gradient replay errors were **0**.
- Original incomplete attempt `09ebadea135ec50dc873` remains separate; no old fits were reused.
- Fitting start (UTC): `2026-11-13T04:30:03.917871+00:00`.
- Fitting end (UTC): `2026-11-13T04:41:45.932490+00:00`; elapsed **702.0 seconds**.
- Validation finished (UTC): `2026-11-13T04:42:39.657665+00:00`.

Results: [evaluation_results.md](evaluation_results.md) and its JSON sibling. Numerical audit: [evaluation_result_validation.md](evaluation_result_validation.md) and its JSON sibling. Independent review prompt: [claude_results_review_prompt.txt](claude_results_review_prompt.txt).

The only amendment was a larger iteration cap, 1,000 to 5,000. Convergence tolerances, source/data locks, groups, folds, label budgets, hyperparameter grids and model code remained unchanged. The same matching runtime is recorded in `evaluation_amended_runtime.json`.

## Not done

Independent Claude results review and the scientific checks listed in the results report. No new cohort, nonlinear model or posthoc tuning has been run.
