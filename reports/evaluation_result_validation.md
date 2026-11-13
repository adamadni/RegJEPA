# Saved-model replay validation

Proposal: `af05d8348e60ad38b6e4`. Status: **passed**.

Verified 6300 inner fits and 175 final fits; all converged.
Every saved model hash, training-only transform and query-aligned score reproduced; maximum score error 0.0.
All 15 molecular summaries and 15 paired bootstrap comparisons reproduced.
Each setting covers 299 molecules in 250 leakage components.
Largest iteration count: 1186; fits exceeding the original 1,000-iteration cap: 9.
Largest gradient infinity norm: 7.82815737198e-07.

No models were refitted. Numerical checks use the same reviewed implementation and runtime; they are an artifact audit, not an independent replication.
The JSON sibling includes per-fit evidence and hashes.

## Not done

Independent scientific review, pretraining-overlap audit and stronger external/series validation remain outstanding.
