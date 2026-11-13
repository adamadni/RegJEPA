# Reporting-review response and independent recomputation

Proposal: `af05d8348e60ad38b6e4`. Baseline: `b9c94f48003cf8f45e68324922e4833770132cfb`. Computed UTC: `2026-11-13T09:23:11.328626+00:00`.

The supplied review correctly identified missing endpoint directions and descriptive reporting. All reported counts, subgroup examples, full-budget win/loss counts and 146 original-fit comparisons reproduced. No discrepancy with the supplied numerical findings was found. Frozen loading and deterministic grouping/split reconstruction also passed in the pinned WSL runtime.

## Changes

- Added per-component minimum/maximum selections, and qualified the 50-label random-control and JEPA+QM complementarity interpretations.
- Added Morgan ties with molecule-fit versus unique-molecule denominators and directly checked their causes using saved masks/features/weights.
- Added per-molecule negative-information-gain counts after repetition averaging.
- Added candidate-count breakdowns at every budget, with exact-k detail in JSON and pooled 7+ reporting in Markdown.
- Replaced an unverifiable inference about information access with an attributed execution account and independently checkable unchanged-output evidence.

## Morgan tie causes

Of 141/70/5 tied molecule-fits at 50/100/full, 137/67/4 involve distinct original fingerprints becoming identical under the training mask; all differing dropped columns were zero throughout training. Another 3/3/1 involve identical full fingerprints. One remaining 50-label event, comp414 at fold 2 / seed 17, ties despite differing retained vectors with two differing nonzero-weight columns. No scorer or tie rule was changed. The four other blocks have no reported-positive ties.

The mask tests within-molecule training variation, not merely whether a bit was ever nonzero. Unseen bits explain the observed mask-collapsed cases, but are not the full definition of the mask. The audit JSON includes event identifiers and pair classifications.

These are 6,144-column atom-centered fingerprints, distinct from the whole-molecule fingerprints used for leakage grouping. For training-zero columns, the likelihood gradient is zero and the positive L2 term is minimized at zero coefficient. Retaining those columns would therefore not by itself fix these ties under the same linear objective.

## Full-budget paired outcomes

| Left / right | Left alone correct | Right alone correct | Both correct | Both incorrect |
|---|---:|---:|---:|---:|
| Atom-JEPA / Random backbone | 34 | 31 | 192 | 42 |
| JEPA + QM / Electronic (xTB) | 12 | 8 | 248 | 31 |

## Interpretation qualifications

An endpoint selection is a warning about the bounded search, not proof that extending it improves held-out performance. The possible bias toward Atom-JEPA at 50 labels is unquantified. The original pretrained/random and combined/QM top-1 contrasts remain unresolved; no equivalence or no-benefit claim is made.

The review's posthoc paired NLL intervals were not recomputed or promoted into the approved results. Bootstrap tail shares were not adopted as multiplicity-adjusted hypothesis tests; the original paired top-1 intervals remain exploratory and unadjusted. Broader NLL inference, grids and random-backbone seeds belong in a new prospective protocol, not a modification of this frozen run.

All 146 original archives and fit JSON files are byte-identical to their amended-run counterparts. That establishes unchanged outputs, not whether metrics were viewed. Historical outer_performance_inspected flags are execution statements, not access-log evidence.

## Verification scope

Checked all 175 saved model hashes; loaded and rebuilt the frozen protocol in the pinned runtime; replayed all 35 saved Morgan models exactly to diagnose ties; recomputed descriptive tables; compared 146 original/amended model/score/penalty artifacts. Models refitted: **0**. New bootstrap intervals: **0**. The previous all-model replay audit remains unchanged. All numerical details and evidence are in the JSON sibling.

## Not done

No wider-grid experiment, additional encoder seed, new cohort, pretraining-membership audit, new inferential NLL analysis or independent complete review signoff. A focused Claude re-review prompt accompanies these reporting changes.
