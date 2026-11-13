# Approved amended evaluation execution status

Proposal **`af05d8348e60ad38b6e4`** is user-approved, frozen, and ready for a fresh all-five-block run. The exact approval is "Approve cap 5,000 and rerun"; its binding is in `data/splits/af05d8348e60ad38b6e4/approval.json`.

The sole change from the original reviewed protocol is the iteration cap: **1,000 to 5,000**. Convergence criteria, input locks, split assignments and model settings are unchanged. Frozen loading and deterministic reconstruction passed in the same exact-version Linux runtime. The reviewed source remains at `977f120` with subsequent artifact-only commits.

Planned workload: **6,300 inner fits and 175 final fits**. The original `09ebadea135ec50dc873` attempt remains separately preserved with 146 final fits; none will be reused. This amended run has not started yet.

Runtime: `evaluation_amended_runtime.json`. Machine-readable execution state: `evaluation_amended_execution_status.json`.

## Not done

Fresh fitting, convergence verification, full result comparison and results review.
