# Frozen-feature linear evaluation results

Approved proposal **`af05d8348e60ad38b6e4`** completed **6,300 inner fits and 175 final fits**. All fits passed the unchanged convergence rule. The iteration-cap amendment was approved before this fresh run; the original incomplete attempt remains separate.

The comparison covers **299 conservative bromination records / 299 structures / 1,264 distinct-outcome rows**, with **250** leakage components. Five outer folds test 59 or 60 molecules each; full training uses 239 or 240 labels. The two smaller budgets use exactly 50 or 100 labels, including all labels available to inner selection. Three subset repetitions are averaged per molecule at 50/100; full training is fitted once per outer fold.

Electronic descriptors are the strongest single block under this linear readout. Atom-JEPA outperforms atom-centered Morgan at all three budgets, while electronic descriptors outperform Atom-JEPA; the corresponding exploratory paired top-1 intervals exclude zero throughout.

The pilot **does not resolve a top-1 advantage from pretraining over this random backbone, or from adding JEPA to QM**. Both contrasts have positive point estimates at every budget, but all six conditional intervals include zero. At full training, Atom-JEPA minus random is **+1.00 percentage point [−4.73, +6.86]**; JEPA+QM minus electronic is **+1.34 points [−2.37, +5.12]**. These are inconclusive comparisons, not evidence of equivalence or of no possible benefit.

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

Top-2 restricted to more than two candidates avoids the 55 trivial two-candidate cases (244 molecules remain). NLL and information gain are in nats per molecule; uniform NLL is 1.3525221466. Negative information gain, if present, means worse NLL than uniform. No secondary-metric intervals were computed.

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

## Numerical and model-selection diagnostics

Every saved model hash, train-only transform and held-out score replayed successfully. All molecular summaries and paired intervals reproduced without fitting any model. Maximum score-replay error: **0.0**. Largest iteration count: **1186**; **9** fits exceeded 1,000 iterations. Largest saved gradient infinity norm: **7.82815737198e-07**.

L-BFGS-B uses gradient target 1e-7 and relative-objective tolerance 1e-14. Acceptance additionally requires successful termination, finite coefficients/objective and gradient at most 1e-5. A gradient above 1e-7 is valid only with successful termination by the reviewed rule.

The following final fits selected at least one endpoint of the unchanged L2 grid. Endpoints were flagged, not expanded using outer-test performance. The JSON sibling retains per-component minimum/maximum counts.

| Feature block | Budget | Fits selecting a grid edge |
|---|---|---:|
| Atom-JEPA | 50 | 0 / 15 |
| Atom-JEPA | 100 | 0 / 15 |
| Atom-JEPA | full | 0 / 5 |
| Random backbone | 50 | 7 / 15 |
| Random backbone | 100 | 0 / 15 |
| Random backbone | full | 0 / 5 |
| Atom-centered Morgan | 50 | 2 / 15 |
| Atom-centered Morgan | 100 | 0 / 15 |
| Atom-centered Morgan | full | 0 / 5 |
| Electronic (xTB) | 50 | 1 / 15 |
| Electronic (xTB) | 100 | 0 / 15 |
| Electronic (xTB) | full | 2 / 5 |
| JEPA + QM | 50 | 6 / 15 |
| JEPA + QM | 100 | 2 / 15 |
| JEPA + QM | full | 1 / 5 |

JEPA+QM uses separate penalties and a 36-candidate inner search; single blocks use 6. Its QM-relative RMS diagnostics are score amplitudes, not independent information shares or product fractions. Correlated components can reinforce or cancel.

| Budget | Minimum QM relative RMS | Median | Maximum |
|---|---:|---:|---:|
| 50 | 0.5162 | 0.7853 | 0.9087 |
| 100 | 0.4675 | 0.7322 | 0.8852 |
| full | 0.5150 | 0.6730 | 0.9035 |

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

Independent Claude results review, pretraining-membership audit, additional random-backbone seeds and stronger external/series validation. No next cohort, nonlinear head or posthoc tuning has been run.
