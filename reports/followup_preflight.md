# Exploratory follow-up preflight

Proposal: `1ddb9c393580ec2e232f`; **not frozen, no real-data fits**.

Observed: 299 records, 299 structures, 1264 outcome rows, 250 groups.
The complete original grouped fold/subset schedule and existing bank are preserved. This is an adaptive follow-up on already evaluated molecules, not a fresh confirmatory test.

Planned workload: 75,300 inner models and 935 final models; including MLP restarts, **84,750 inner and 1,250 final optimizer attempts**.

| Arm | Head | Columns | Outer settings | Penalty choices | Inner optimizer attempts |
|---|---|---:|---:|---:|---:|
| qm_all_linear | linear | 6 | 35 | 10 | 1050 |
| qm_mulliken_linear | linear | 1 | 5 | 10 | 150 |
| qm_cm5_linear | linear | 1 | 5 | 10 | 150 |
| qm_fplus_linear | linear | 1 | 5 | 10 | 150 |
| qm_fminus_linear | linear | 1 | 5 | 10 | 150 |
| qm_fzero_linear | linear | 1 | 5 | 10 | 150 |
| qm_wiberg_linear | linear | 1 | 5 | 10 | 150 |
| qm_charges_linear | linear | 2 | 35 | 10 | 1050 |
| qm_fukui_linear | linear | 3 | 35 | 10 | 1050 |
| qm_fukui_pair_linear | linear | 2 | 35 | 10 | 1050 |
| qm_no_charges_linear | linear | 4 | 35 | 10 | 1050 |
| qm_no_fukui_linear | linear | 3 | 35 | 10 | 1050 |
| qm_without_mulliken_linear | linear | 5 | 5 | 10 | 150 |
| qm_without_cm5_linear | linear | 5 | 5 | 10 | 150 |
| qm_without_fplus_linear | linear | 5 | 5 | 10 | 150 |
| qm_without_fminus_linear | linear | 5 | 5 | 10 | 150 |
| qm_without_fzero_linear | linear | 5 | 5 | 10 | 150 |
| qm_without_wiberg_linear | linear | 5 | 5 | 10 | 150 |
| jepa_linear | linear | 768 | 35 | 10 | 1050 |
| jepa_qm_linear | linear | 774 | 35 | 100 | 10500 |
| random_linear | linear | 768 | 35 | 10 | 1050 |
| random_qm_linear | linear | 774 | 35 | 100 | 10500 |
| morgan_linear | linear | 6144 | 35 | 10 | 1050 |
| morgan_qm_linear | linear | 6150 | 35 | 100 | 10500 |
| jepa_scalar_linear | linear | 256 | 35 | 10 | 1050 |
| jepa_scalar_qm_linear | linear | 262 | 35 | 100 | 10500 |
| random_scalar_linear | linear | 256 | 35 | 10 | 1050 |
| random_scalar_qm_linear | linear | 262 | 35 | 100 | 10500 |
| qm_all_mlp | mlp | 6 | 35 | 10 | 2100 |
| jepa_mlp | mlp | 768 | 35 | 10 | 2100 |
| random_mlp | mlp | 768 | 35 | 10 | 2100 |
| jepa_qm_mlp | mlp | 774 | 35 | 10 | 2100 |
| random_qm_mlp | mlp | 774 | 35 | 10 | 2100 |
| jepa_scalar_mlp | mlp | 256 | 35 | 10 | 2100 |
| random_scalar_mlp | mlp | 256 | 35 | 10 | 2100 |
| jepa_scalar_qm_mlp | mlp | 262 | 35 | 10 | 2100 |
| random_scalar_qm_mlp | mlp | 262 | 35 | 10 | 2100 |

Maximum native-precision residual f0 - (f+ + f-)/2: 0.00050000846. Fukui channels are not independent; leave-one-out scores are not unique or causal importance.

Read docs/followup_protocol.md and reports/claude_followup_review_prompt.txt before freezing.

## Not done

Independent protocol/code review, explicit freeze, real-data execution and all-layer extraction.
