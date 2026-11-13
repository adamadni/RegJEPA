# Exploratory descriptor, hybrid and readout follow-up

This protocol implements the user's request to identify useful QM descriptors, run the missing hybrid controls and compare JEPA readouts. It does not perform reactivity pretraining, identify a mechanism, estimate product fractions or establish causal descriptor importance.

Review proposal **`1ddb9c393580ec2e232f`**, generated from `configs/followup.json`. It is **not frozen**. No real-data model has been fitted under it. Exact arm widths and fit counts are in `reports/followup_preflight.md` and `.json`.

## Dataset, grouping and interpretation

Reuse the 299-record conservative bromination cohort and feature run `d76c42313dcf5e17ac96`. Every arm preserves the complete candidate set of each molecule. Bank rows are already symmetry-averaged outcomes; no second atom pooling or symmetry-multiplicity weighting is introduced.

The entire original five-outer/three-inner grouped schedule, 50/100-label subsets, subset seeds 17/29/43 and full-training sets are copied exactly from `af05d8348e60ad38b6e4`. All training, scaling and inner selection remain inside the selected label budget. Full training uses 239 or 240 records per outer fold. Existing source/series uncertainties and comp480/comp481 remain limitations. Source and similarity components, not compound IDs alone, determine separation.

This is an **adaptive exploratory follow-up on previously evaluated molecules**. Freezing a new analysis does not turn the reused outer folds into a fresh confirmatory test. New models are never selected using outer outcomes; comparisons among fixed arms can generate hypotheses for independent data. No stronger source-family, unfamiliar-substrate or out-of-distribution claim is made.

## QM panel: 18 linear arms

The literal electronic columns are Mulliken charge, CM5 charge, f(+), f(-), f(0) and total Wiberg bond order. Evaluate:

- All six descriptors.
- Each of the six single descriptors.
- Charges together; the three Fukui columns together; the two-column f(+)/f(-) pair.
- All six minus each individual descriptor.
- All six minus the charge group; all six minus the Fukui group.

The six singles and six individual leave-one-out arms run at full training only. The full panel, charge/Fukui groups, independent Fukui pair and group omissions retain the 50/100/full learning curves. This defines the analysis in advance rather than running every possible subset and reporting its winner.

The native bank satisfies f(0) approximately equal to (f(+) + f(-))/2, with maximum observed rounding residual 0.0005000085. A near-zero leave-one-out effect for a Fukui column cannot establish irrelevance: the others may reconstruct it. Charges can also be correlated. Group ablations and individual-descriptor sufficiency are complementary diagnostics. Redundant columns can change effective regularization as well as information content, so none of these contrasts estimates unique physical causation.

## Structural and hybrid controls: ten linear arms

Standalone blocks: the published frozen Atom-JEPA, frozen random backbone seed 42, atom-centered Morgan, scalar-only Atom-JEPA and scalar-only random backbone. Pair each with the full six-column QM block. This supplies QM+random and QM+Morgan controls in addition to QM+Atom-JEPA.

The current invariant feature vector has 768 columns: 256 final-layer l=0 scalars, then 256 channel norms each for l=1 and l=2. Scalar-only uses the first 256 columns of the existing bank, exactly matching the pinned upstream API's flattening order. The random and pretrained definitions are identical. Both encoders and the shared first conformer remain fixed.

These controls ask whether any hybrid gain is specific to pretraining rather than generic geometry, architecture, extra features or structural fingerprints. No new random-backbone initialization has been extracted in this phase; head restart seeds are not encoder seeds. Pretraining-specific conclusions still require later encoder-seed replication and overlap auditing.

## Nonlinear readouts: nine additional arms

Use one 16-unit tanh hidden layer and a shared scalar site score, without an output bias. The encoder remains frozen. Evaluate QM-only; pretrained/random 768-column blocks and their QM hybrids; pretrained/random 256-column scalar blocks and their QM hybrids.

Each fit uses initializations 11 and 23 and retains the lowest regularized **training** objective, within the same hyperparameter choice. Both restarts must converge. Inner validation selects the penalty; neither restart choice nor stopping sees the outer data. Nonlinear optimization is local and two restarts do not establish a global optimum.

Every hidden weight, hidden bias and output weight receives the shared positive L2 penalty. Linear hybrids instead use independent representation/QM penalties. Thus hybrid head contrasts compare complete readout procedures, including their regularization parameterization, rather than isolating nonlinear activation alone. QM-only and standalone representation head comparisons have the same shared-penalty policy, though the parameter norms and model capacities differ.

## Preprocessing shared across dense readouts

Fit per-column means and RMS deviations on training data only, weighting every molecule equally and then every candidate within its molecule equally. Use these global training statistics for both linear and nonlinear dense inputs. Drop columns with RMS at most 1e-12. Morgan columns in linear arms remain unscaled; their fitted mean can be subtracted without changing linear within-molecule choices.

Do not subtract each molecule's input feature mean before a nonlinear head. A feature constant across one molecule's candidates can still modulate site preferences through an interaction. The prior pilot's within-molecule RMS is therefore replaced by common global RMS for dense follow-up arms, making preprocessing consistent between heads while preserving molecular context.

This change alters effective regularization for linear models. Every reference model is rerun under the follow-up definition; old and new fit metrics must not be pooled or treated as identical procedures. It is a named methodological change, not a silent repair of the original study.

## Optimization, tuning and workload

The widened L2 grid is 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100 and 1,000. Single-block linear and all MLP arms search ten shared penalties. Linear hybrids search all 100 independent representation/QM pairs. Exact inner-NLL ties prefer stronger penalties in the fixed descending candidate order.

Minimize equal-molecule listwise NLL plus the specified L2 term. L-BFGS-B starts at zero for linear models and seeded fan-in initialization for MLPs. The wider grid can be harder to optimize, so the proposed ceiling is 20,000 iterations, with max function evaluations equal to 50 times that ceiling, max line-search steps 50 and maxcor 10. Gradient target 1e-7, ftol 1e-14, successful termination, finite parameters/objective and accepted gradient at most 1e-5 retain the prior acceptance semantics. Failure aborts; no failed or nonconverged restart is ignored.

The schedule contains **75,300 inner model fits and 935 final model fits**. Accounting for two MLP restarts, it requires **84,750 inner optimizer attempts and 1,250 final optimizer attempts**. This is substantially larger than the original pilot. No runtime or cloud-cost estimate has been established. Review this workload before approving execution.

## Reporting and uncertainty

Predeclare paired top-1 and NLL reporting, plus descriptive top-2, nontrivial top-2, MRR, information gain, negative-gain molecule counts, positive-outcome ties and candidate-count breakdowns. Rank outcome rows directly with the original exact-tie expectation rule. Fixed per-row float64 reductions preserve identical-input ties for linear and nonlinear scoring; no approximate tie threshold is used.

At 50/100 labels, average metrics over the three subset repetitions per molecule. Do not ensemble logits or treat repetitions as independent molecules. Full training has one repetition. All matched arms/budgets evaluate the same 299 records.

The configuration names 49 paired contrasts, yielding 123 pair/budget cells and **246 conditional intervals** over top-1 and NLL. Use 10,000 paired bootstrap draws of the original leakage components, seed 42. Positive top-1 differences favor the left arm; negative NLL differences favor it. Intervals are exploratory and unadjusted, conditional on the fitted models, not a family-wise confirmatory test or a complete estimate of retraining uncertainty. Report all configured comparisons, not only favorable ones. No new NLL inference is retroactively added to the original frozen experiment.

Report minimum versus maximum penalty selections per component, including denominator counts. A boundary selection is a warning, not proof that an outside-grid value improves performance. Do not adapt the grid from outer outcomes. No component-amplitude diagnostic is interpreted as a unique information share; nonlinear scores do not have the prior linear additive decomposition.

## Isolation, dependencies and execution

The new code lives under `src/regjepa_followup`, leaving `src/regjepa` unchanged so its source-glob freeze bindings remain reproducible. Dependencies are `requirements-followup.txt`, which includes the already pinned evaluation requirements. No sklearn, pandas, PyTorch or new runtime dependency is imported by these readouts. Use the existing exact-version Linux evaluation runtime; source-feature generation has separate dependencies.

Preparation and inspection fit no real-data model:

```text
python -m pytest -q tests/test_followup.py --basetemp .cache/pytest-followup-review
python scripts/prepare_followup.py
python scripts/run_followup.py --splits data/followup_splits/1ddb9c393580ec2e232f/proposal.json
```

After independent review and approval of the exact proposal, explicitly freeze and execute:

```text
python scripts/prepare_followup.py --freeze data/followup_splits/1ddb9c393580ec2e232f/proposal.json
python scripts/run_followup.py --splits data/followup_splits/1ddb9c393580ec2e232f/frozen.json --run
python scripts/validate_followup.py --splits data/followup_splits/1ddb9c393580ec2e232f/frozen.json --output reports/followup_runs/1ddb9c393580ec2e232f --report reports/followup_validation.json
```

Frozen loading reconstructs the original plan, current follow-up source locks and exact sibling-proposal binding. It establishes reproducibility, not human approval. Record the approved proposal ID and source HEAD separately. Existing run directories are never overwritten and there is no automatic resume. Per-fit artifacts survive failure, and the failing task is recorded. The validator replays selected final models, transforms, scores and gradients and reproduces summaries/intervals; inner and nonselected restart gradients are checked from saved diagnostics only.

## Deferred extraction and later research

All-layer or intermediate-layer features cannot be recovered from the current last-layer arrays. The pinned upstream API can extract them, but Windows currently blocks its xxhash dependency through Application Control. No security settings were changed. A separate compatible feature runtime and a new, provenance-bound export are needed before that experiment. The original bank and latest pointer must remain intact. The present readout tests concern scalar versus current invariant features and head families, not which encoder layer is best.

Reactivity-focused pretraining, additional encoder seeds, conformer ensembles, new cohorts and pretraining-membership auditing remain separate work. First use these controls to identify useful targets and distinguish electronic information, generic structural complements and adaptation of the existing embedding.

## Not done

Independent Claude review, explicit approval/freeze and real-data execution. The implemented synthetic checks are validation of the harness, not scientific results. This proposal has not been approved or run.
