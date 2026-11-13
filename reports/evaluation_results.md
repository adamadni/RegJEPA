# Frozen-feature linear evaluation results

Approved proposal **`af05d8348e60ad38b6e4`** completed **6,300 inner fits and 175 final fits**. All fits passed the unchanged convergence rule. The iteration-cap amendment was approved before this fresh run; the original incomplete attempt remains separate.

The comparison covers **299 conservative bromination records / 299 structures / 1,264 distinct-outcome rows**, with **250** leakage components. Five outer folds test 59 or 60 molecules each; full training uses 239 or 240 labels. The two smaller budgets use exactly 50 or 100 labels, including all labels available to inner selection. Three subset repetitions are averaged per molecule at 50/100; full training is fitted once per outer fold.

Electronic descriptors are the strongest single block under this linear readout. Atom-JEPA outperforms atom-centered Morgan at all three budgets, while electronic descriptors outperform Atom-JEPA; the corresponding exploratory paired top-1 intervals exclude zero throughout.

The pilot **does not resolve a top-1 advantage from pretraining over this random backbone, or from adding JEPA to QM**. All six conditional intervals include zero. The 50-label Atom-JEPA-versus-random contrast is additionally sensitive to the L2 search boundary: random selected the strongest available penalty in 7 of 15 final fits. Its apparent advantage must not be treated as directional evidence of pretraining benefit. At full training, Atom-JEPA minus random is **+1.00 percentage point [−4.73, +6.86]**; JEPA+QM minus electronic is **+1.34 points [−2.37, +5.12]**. These are inconclusive comparisons, not evidence of equivalence or of no possible benefit. JEPA+QM selected the strongest tested JEPA penalty in 4/15, 2/15 and 1/5 fits at 50/100/full, respectively; that tendency toward stronger JEPA shrinkage belongs alongside any complementarity interpretation.

JEPA+QM has lower descriptive NLL than electronic descriptors at each budget: **0.5405 vs 0.5581** (50), **0.4231 vs 0.5126** (100), and **0.4017 vs 0.4551** (full). This suggests a possible improvement in probabilistic prediction worth reviewing, but no paired NLL interval was part of this protocol. Its softmax values describe source-annotation choices, not measured product fractions. The comparisons assess these feature blocks and this regularized linear head; they do not establish how all possible Morgan, QM or neural models would perform.

## Primary top-1 result

Values are molecule-average expected top-1 percentages, with fractional credit for exact ties. Each cell evaluates the same 299 molecules; repetitions are not extra test examples.

| Feature block | 50 labels | 100 labels | Full (239–240 labels) |
|---|---:|---:|---:|
| Atom-JEPA | 64.55% | 71.35% | 75.59% |
| Random backbone | 61.43% | 67.34% | 74.58% |
| Atom-centered Morgan | 44.06% | 51.34% | 63.55% |
| Electronic (xTB) | 82.61% | 85.28% | 85.62% |
| JEPA + QM | 82.94% | 85.95% | 86.96% |
| Uniform over distinct outcomes | 28.31% | 28.31% | 28.31% |

## Paired top-1 differences

Differences and conditional 95% intervals are in percentage points (left minus right). Intervals use 10,000 paired resamples of leakage components. The 15 comparisons are exploratory and unadjusted; these intervals exclude uncertainty from repeating the entire training/model-selection procedure.

| Left minus right | Label budget | Difference | Conditional 95% interval |
|---|---|---:|---|
| Atom-JEPA minus Random backbone | 50 | +3.12 | [-1.95, +8.31] |
| Atom-JEPA minus Atom-centered Morgan | 50 | +20.48 | [+14.94, +25.90] |
| Atom-JEPA minus Electronic (xTB) | 50 | -18.06 | [-23.47, -12.67] |
| JEPA + QM minus Atom-JEPA | 50 | +18.39 | [+13.61, +23.17] |
| JEPA + QM minus Electronic (xTB) | 50 | +0.33 | [-2.06, +2.51] |
| Atom-JEPA minus Random backbone | 100 | +4.01 | [-1.00, +8.86] |
| Atom-JEPA minus Atom-centered Morgan | 100 | +20.01 | [+14.35, +25.45] |
| Atom-JEPA minus Electronic (xTB) | 100 | -13.94 | [-19.11, -8.95] |
| JEPA + QM minus Atom-JEPA | 100 | +14.60 | [+10.50, +18.89] |
| JEPA + QM minus Electronic (xTB) | 100 | +0.67 | [-2.11, +3.49] |
| Atom-JEPA minus Random backbone | full | +1.00 | [-4.73, +6.86] |
| Atom-JEPA minus Atom-centered Morgan | full | +12.04 | [+5.61, +18.64] |
| Atom-JEPA minus Electronic (xTB) | full | -10.03 | [-16.11, -3.93] |
| JEPA + QM minus Atom-JEPA | full | +11.37 | [+6.21, +16.72] |
| JEPA + QM minus Electronic (xTB) | full | +1.34 | [-2.37, +5.12] |

## Secondary descriptive metrics

Top-2 restricted to more than two candidates avoids the 55 trivial two-candidate cases (244 molecules remain). NLL and information gain are in nats per molecule; uniform NLL is 1.3525221466. Negative per-molecule information gain means worse NLL than uniform; such cases occur in every block even though every aggregate gain is positive. No secondary-metric intervals were computed.

| Feature block | Budget | Top-2 | Top-2 (k>2) | MRR | NLL | Information gain |
|---|---|---:|---:|---:|---:|---:|
| Atom-JEPA | 50 | 0.8606 | 0.8292 | 0.7955 | 0.8755 | 0.4770 |
| Atom-JEPA | 100 | 0.9064 | 0.8852 | 0.8382 | 0.7343 | 0.6182 |
| Atom-JEPA | full | 0.9298 | 0.9139 | 0.8649 | 0.6167 | 0.7358 |
| Random backbone | 50 | 0.8439 | 0.8087 | 0.7749 | 1.0289 | 0.3236 |
| Random backbone | 100 | 0.9030 | 0.8811 | 0.8171 | 0.8550 | 0.4975 |
| Random backbone | full | 0.9465 | 0.9344 | 0.8627 | 0.6651 | 0.6875 |
| Atom-centered Morgan | 50 | 0.6983 | 0.6302 | 0.6537 | 1.2474 | 0.1051 |
| Atom-centered Morgan | 100 | 0.7897 | 0.7423 | 0.7120 | 1.1259 | 0.2266 |
| Atom-centered Morgan | full | 0.9047 | 0.8832 | 0.7967 | 0.9624 | 0.3902 |
| Electronic (xTB) | 50 | 0.9331 | 0.9180 | 0.8991 | 0.5581 | 0.7944 |
| Electronic (xTB) | 100 | 0.9554 | 0.9454 | 0.9168 | 0.5126 | 0.8400 |
| Electronic (xTB) | full | 0.9699 | 0.9631 | 0.9211 | 0.4551 | 0.8974 |
| JEPA + QM | 50 | 0.9521 | 0.9413 | 0.9052 | 0.5405 | 0.8120 |
| JEPA + QM | 100 | 0.9688 | 0.9617 | 0.9236 | 0.4231 | 0.9295 |
| JEPA + QM | full | 0.9666 | 0.9590 | 0.9279 | 0.4017 | 0.9508 |

## Molecules scoring worse than uniform

Counts use the per-molecule information gain after averaging the three subset repetitions at 50/100; equivalently, a molecule's mean NLL exceeds its own `log(n_candidates)`. Full training has one repetition. Each denominator is 299; these are not candidate-row or molecule-fit counts. The JSON sibling retains zero/positive counts and negative-case identifiers.

| Feature block | 50 labels | 100 labels | Full |
|---|---:|---:|---:|
| Atom-JEPA | 61 | 51 | 44 |
| Random backbone | 84 | 75 | 46 |
| Atom-centered Morgan | 109 | 81 | 64 |
| Electronic (xTB) | 30 | 25 | 22 |
| JEPA + QM | 30 | 25 | 21 |

## Exact ties involving the reported Morgan outcome

A molecule-fit is one held-out molecule under one training-subset repetition. Low-label budgets contain 897 molecule-fits (299 molecules times three repetitions); full training contains 299. Fractional tie credit is unchanged. The table separates tied evaluations from distinct molecules affected.

Here, a full fingerprint means all 6,144 columns of the atom-centered Morgan feature block, not the whole-molecule fingerprint used for grouping.

| Budget | Tied molecule-fits | Distinct molecules affected | Dropped unseen columns explain tie* | Identical full fingerprints* | Distinct retained vectors* |
|---|---:|---:|---:|---:|---:|
| 50 | 141 / 897 | 79 / 299 | 137 | 3 | 1 |
| 100 | 70 / 897 | 47 / 299 | 67 | 3 | 0 |
| full | 5 / 299 | 5 / 299 | 4 | 1 | 0 |

*The last three columns partition the tied molecule-fits. In the mask category, at least one tied alternative differs in the full fingerprint, but every tied alternative equals the reported outcome after the saved training mask. All differing dropped columns in these cases were zero throughout training. The mask rule itself is broader: it drops columns without within-molecule training variation, including potentially seen-but-constant columns. No such seen-column explanation was needed in these observed cases.

Simply retaining the unseen columns would not resolve these ties under this linear objective: columns that are zero throughout training have no likelihood gradient, and positive L2 regularization gives them zero optimal coefficients. The mask exposes a training-coverage limitation rather than demonstrating an avoidable scoring defect.

The one remaining 50-label event is comp414 under fold 2 / subset seed 17: two retained columns differ and have nonzero coefficients, yet the saved linear scores tie exactly. It is classified as a projection tie, not a masking or full-fingerprint collision. These diagnostics support a vocabulary-coverage explanation for most low-label Morgan ties, not all ties or all of its performance gap. The other four blocks have no ties involving the reported outcome in the saved molecule-fits.

## Descriptive breakdown by candidate count

Values are molecular top-1 percentages, using the same per-molecule repetition averages as the main result. The 7+ bin pools 28 seven-candidate, four eight-candidate, four nine-candidate and one twelve-candidate molecule. This pooling was added after review for readable descriptive reporting; no subgroup intervals, tests or model choices were added. Exact-k counts and all descriptive metrics are in the JSON sibling.

| Budget | Distinct outcomes | Molecules | Atom-JEPA | Random | Morgan | Electronic | JEPA + QM |
|---|---|---:|---:|---:|---:|---:|---:|
| 50 | 2 | 55 | 79.39% | 80.00% | 70.00% | 93.94% | 92.73% |
| 50 | 3 | 67 | 55.72% | 47.26% | 35.57% | 85.57% | 83.08% |
| 50 | 4 | 62 | 66.13% | 69.89% | 40.95% | 79.03% | 76.88% |
| 50 | 5 | 39 | 63.25% | 61.54% | 42.81% | 79.49% | 76.07% |
| 50 | 6 | 39 | 62.39% | 55.56% | 40.74% | 72.65% | 79.49% |
| 50 | 7+ | 37 | 59.46% | 51.35% | 30.93% | 80.18% | 89.19% |
| 100 | 2 | 55 | 83.64% | 83.03% | 73.94% | 94.55% | 92.73% |
| 100 | 3 | 67 | 60.20% | 50.75% | 41.79% | 84.58% | 84.58% |
| 100 | 4 | 62 | 73.66% | 73.12% | 49.19% | 86.02% | 83.87% |
| 100 | 5 | 39 | 71.79% | 64.10% | 48.86% | 81.20% | 78.63% |
| 100 | 6 | 39 | 70.09% | 65.81% | 47.72% | 73.50% | 81.20% |
| 100 | 7+ | 37 | 70.27% | 69.37% | 45.05% | 88.29% | 94.59% |
| full | 2 | 55 | 81.82% | 89.09% | 69.09% | 94.55% | 94.55% |
| full | 3 | 67 | 65.67% | 59.70% | 52.24% | 83.58% | 82.09% |
| full | 4 | 62 | 79.03% | 77.42% | 65.32% | 87.10% | 87.10% |
| full | 5 | 39 | 82.05% | 74.36% | 60.26% | 82.05% | 79.49% |
| full | 6 | 39 | 69.23% | 69.23% | 58.97% | 74.36% | 82.05% |
| full | 7+ | 37 | 78.38% | 81.08% | 81.08% | 89.19% | 97.30% |

At full training, random exceeds Atom-JEPA on the 55 two-candidate molecules (89.09% versus 81.82%). JEPA+QM's net advantage over electronic descriptors is three additional hits among the 37 molecules with seven or more candidates (97.30% versus 89.19%), and one additional net hit among the other 262. This localizes the observed point difference; these small, post-review subgroups do not establish generalization or a mechanism.

## Numerical and model-selection diagnostics

Every saved model hash, train-only transform and held-out score replayed successfully. All molecular summaries and paired intervals reproduced without fitting any model. Maximum score-replay error: **0.0**. Largest iteration count: **1186**; **9** fits exceeded 1,000 iterations. Largest saved gradient infinity norm: **7.82815737198e-07**.

L-BFGS-B uses gradient target 1e-7 and relative-objective tolerance 1e-14. Acceptance additionally requires successful termination, finite coefficients/objective and gradient at most 1e-5. A gradient above 1e-7 is valid only with successful termination by the reviewed rule.

The L2 grid ranges from **0.0001** (minimum, weaker shrinkage) to **10** (maximum, stronger shrinkage). Counts below are final fits selecting each endpoint, separated by penalty component. No grid was widened or model refitted in response to this review.

| Feature block | Budget | Penalty component | Minimum selected | Maximum selected | Final fits |
|---|---|---|---:|---:|---:|
| Atom-JEPA | 50 | shared | 0 | 0 | 15 |
| Atom-JEPA | 100 | shared | 0 | 0 | 15 |
| Atom-JEPA | full | shared | 0 | 0 | 5 |
| Random backbone | 50 | shared | 0 | 7 | 15 |
| Random backbone | 100 | shared | 0 | 0 | 15 |
| Random backbone | full | shared | 0 | 0 | 5 |
| Atom-centered Morgan | 50 | shared | 0 | 2 | 15 |
| Atom-centered Morgan | 100 | shared | 0 | 0 | 15 |
| Atom-centered Morgan | full | shared | 0 | 0 | 5 |
| Electronic (xTB) | 50 | shared | 1 | 0 | 15 |
| Electronic (xTB) | 100 | shared | 0 | 0 | 15 |
| Electronic (xTB) | full | shared | 2 | 0 | 5 |
| JEPA + QM | 50 | jepa | 0 | 4 | 15 |
| JEPA + QM | 50 | qm | 2 | 0 | 15 |
| JEPA + QM | 100 | jepa | 0 | 2 | 15 |
| JEPA + QM | 100 | qm | 0 | 0 | 15 |
| JEPA + QM | full | jepa | 0 | 1 | 5 |
| JEPA + QM | full | qm | 0 | 0 | 5 |

Random's maximum-penalty selections in 7/15 low-label fits flag a potentially underregularized control and a comparison that may favor Atom-JEPA. A selected endpoint does **not** prove that an outside-grid optimum exists, quantify the bias, or show that more shrinkage would improve held-out top-1. The positive 50-label point difference is not evidence of pretraining benefit. Morgan also selects maximum shrinkage in 2/15 such fits.

JEPA+QM repeatedly selects maximum JEPA shrinkage (4/15, 2/15, 1/5 at 50/100/full), while its QM component selects the minimum in 2/15 low-label fits. Electronic alone selects the minimum in 1/15 low-label and 2/5 full-training fits. These directions matter when interpreting complementarity and physical-descriptor baselines. Finite L2 penalties do not remove an entire block, and boundary selections do not establish an optimal zero contribution. Any grid extension requires a new prospective proposal.

JEPA+QM uses separate penalties and a 36-candidate inner search; single blocks use 6. Its QM-relative RMS diagnostics are score amplitudes, not independent information shares or product fractions. Correlated components can reinforce or cancel.

| Budget | Minimum QM relative RMS | Median | Maximum |
|---|---:|---:|---:|
| 50 | 0.5162 | 0.7853 | 0.9087 |
| 100 | 0.4675 | 0.7322 | 0.8852 |
| full | 0.5150 | 0.6730 | 0.9035 |

## What the iteration-amendment artifacts establish

The execution record states that the larger cap was chosen from convergence diagnostics rather than held-out performance. That statement is an execution account, not an independently verifiable access history: the original attempt had already written held-out metrics for 146 fits.

All **146** saved NPZ archives and per-fit JSON files from the original attempt are byte-identical in the amended run, including weights, masks, scales, scores and selected penalties. The amendment changed none of those outputs and none of the four complete single-block comparisons. This does not establish what anyone viewed or prove the decision process was independent of available metrics. See `evaluation_reporting_audit.md` and `.json` for the per-fit comparison.

## Interpretation limits

- Exploratory, unadjusted paired top-1 intervals are conditional on the fitted models; other metrics are descriptive.
- Bibliographic/patent families are incomplete; comp480/comp481 span folds. This is not proven independent-series or OOD evaluation.
- Pretraining overlap is unaudited and the random backbone has only one initialization seed.
- Features describe substrates without reaction conditions; unreported sites are not confirmed negative products.
- The bank uses one shared conformer and retains all seven geometry-fallback records.
- The numerical L2 grid can imply different effective shrinkage across blocks; edge selections are diagnostic, not corrected posthoc.
- The study does not isolate the JEPA pretraining objective or infer mechanisms, stereoselectivity or product fractions.

Frozen-feature gains would demonstrate accessible predictive information under this readout and split. They would not identify which pretraining objective caused that information or establish broader chemistry/mechanistic generalization.

Source artifacts: `evaluation_runs/af05d8348e60ad38b6e4/summary.md` and `.json`, query-aligned per-fit predictions and coefficients, `evaluation_result_validation.md` and `.json`, and the JSON sibling of this report. Independent review prompt: `claude_results_review_prompt.txt`.

## Not done

The supplied Claude review's reporting findings have been independently checked and incorporated. Re-review of these amendments and the reviewer's unchecked gates remain outstanding, as do pretraining-membership auditing, additional random-backbone seeds and stronger external/series validation. A new protocol must specify paired NLL intervals, wider L2 searches and encoder-seed replication before the next cohort. Claude's posthoc NLL intervals and adjustment claims have not been adopted as approved experimental results. No next cohort, nonlinear head or posthoc tuning has been run.
