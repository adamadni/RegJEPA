# Approved evaluation execution status

Approved proposal: `09ebadea135ec50dc873`; approved source commit: `977f1206632a846638ca43606cc6a7375ddbeadf`.

State: **stopped by the reviewed convergence guard; full experiment incomplete**.

- Completed final fits: **146 / 175**.
- Completed inner fits recorded in finished tasks: **3,168 / 6,300** (additional inner fits in the failed task are not yet counted).
- All 35 tasks each for Atom-JEPA, random backbone, Morgan and electronic descriptors are saved; six JEPA+QM tasks are saved.
- Failed task: `jepa_plus_qm_fold0_budgetfull_seedNone`, during inner fitting.
- Stop: **1,000 iterations reached**, with gradient infinity norm **3.26982777983e-07**.
- L-BFGS-B uses gradient target **1e-7** and relative-objective tolerance **1e-14**. The reviewed acceptance rule requires successful optimizer termination, finite objective/coefficients, and gradient at most **1e-5**. This fit failed because the optimizer reported its iteration limit; passing the gradient guard alone is insufficient. No result was silently accepted or replaced.

An original-cap diagnostic replay reproduced the exact failure at inner fold **1**, with **160** training molecules and penalties **JEPA = 0.0001, QM = 10**. It completed 91 inner fits before failing on fit 92. Replay fits are not replacements for the original attempt and are not mixed into its saved outputs. See `evaluation_convergence_diagnosis.json`.

An artifact audit verified all 146 saved model hashes and convergence diagnostics for all 3,168 recorded inner fits plus 146 final fits. The execution account attributes the cap amendment to convergence diagnostics rather than outer-test performance; it does not provide an independently auditable access history. A later reporting audit confirms that all 146 original model archives and fit JSON files are byte-identical in the amended run. See `evaluation_failed_attempt_integrity.json` and `evaluation_reporting_audit.json`.

Frozen artifacts, approval and matching Linux runtime are committed at `75cb053`. The approved code and data remain unchanged. Partial artifacts are preserved under `reports/evaluation_runs/09ebadea135ec50dc873/`; they are not presented as a completed five-model study.

## Not done

Full comparison, final paired intervals, and successful completion of all JEPA+QM settings. An iteration-budget amendment requires a new proposal and explicit user approval.
