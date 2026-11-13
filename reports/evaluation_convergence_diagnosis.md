# Original-cap convergence diagnosis

Proposal `09ebadea135ec50dc873` failed during inner fitting for `jepa_plus_qm_fold0_budgetfull_seedNone`. A diagnostic replay used the same frozen protocol, runtime, zero initialization, transformations, penalties and **1,000-iteration** limit. It did not save replacement models or compare held-out predictions.

The replay reproduced the exact original failure:

| Quantity | Observed value |
|---|---:|
| Inner fold | 1 |
| Training molecules | 160 |
| JEPA L2 penalty | 0.0001 |
| QM L2 penalty | 10 |
| Active coefficients | 774 |
| Iterations | 1,000 |
| Regularized training objective | 0.0007126250584747123 |
| Gradient infinity norm | 3.269827779833867e-7 |

L-BFGS-B reported `STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT`. Its gradient target is 1e-7 and relative-objective tolerance is 1e-14. The additional acceptance guard is gradient at most 1e-5, together with successful termination and finite coefficients/objective. Although the gradient guard passed, optimizer success did not; the fit correctly raised an exception.

The replay completed 91 inner fits and failed on fit 92. This is diagnostic work under the original cap, separate from the original run's 3,168 inner fits recorded in completed tasks. Its full fit-level diagnostics are in `evaluation_convergence_diagnosis.json`. The execution account states that cap selection used convergence diagnostics rather than outer-test performance; the saved artifacts do not independently establish what was viewed. The later reporting audit verifies unchanged outputs for all 146 already-completed fits.

## Not done

A larger iteration cap has not been tested or approved. Reproducing the original failure does not establish that 5,000 iterations will suffice.
