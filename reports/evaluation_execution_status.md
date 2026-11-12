# Approved evaluation execution status

Proposal `09ebadea135ec50dc873` was explicitly approved by the user for freezing and the full linear experiment.

- Reviewed source commit: `977f1206632a846638ca43606cc6a7375ddbeadf`.
- Frozen artifact: `data/splits/09ebadea135ec50dc873/frozen.json`.
- Frozen artifact SHA256: `85291f5e62470fe73a084e5851733b54972af3e3b12e62857c78daeeb91314b5`.
- Approval: `data/splits/09ebadea135ec50dc873/approval.json`.
- State: **frozen and ready to launch; no real-data fits completed yet**.
- Planned workload: 6,300 inner fits plus 175 final fits, using the reviewed configuration.

Windows Application Control blocked the existing Windows RDKit DLLs. The existing WSL2 Linux environment was given an isolated cache-local Python 3.12.14 environment with the exact pinned NumPy 2.5.4, RDKit 2026.9.2 and SciPy 1.18.1 packages. It rebuilt the exact approved proposal and passed all 90 tests. No Windows security settings, reviewed implementation, raw data or feature-bank files were changed.

Runtime details and validation are in `reports/evaluation_runtime.json`; the machine-readable execution state is in `reports/evaluation_execution_status.json`.

## Not done

The full experiment has not yet started. Automatic resume, pretraining-overlap recovery and independent publication-family reconstruction are outside this run.
