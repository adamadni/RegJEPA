# Follow-up preparation validation

Proposal `1ddb9c393580ec2e232f` is prepared, **not frozen or run**.

- **9 synthetic tests passed**, 0 failed, in 6.49 seconds in the pinned Linux runtime.
- Tested analytic nonlinear gradients, training-only molecule-balanced scaling, retention of context, nonlinear context interaction, exact input ties/permutations, literal descriptor columns, metadata rejection, nested-label guards, freeze-binding failure cases, saved-model/gradient replay, both paired metrics and immutable outputs.
- Inspection fitted **0** models. The original frozen plan also loaded and rebuilt exactly, fitting **0** models.
- All original grouped outer/inner folds and subset assignments are identical. The old source package, raw/processed data, feature bank/latest pointer, original splits and completed results have no Git diff from `b466fd2`.
- Planned reporting: **49** named pairs, **123** pair/budget cells, **246** exploratory conditional intervals over top-1 and NLL.
- Workload: **75,300** inner models / **935** final models; **84,750** inner and **1,250** final optimizer attempts including MLP restarts.
- No new runtime dependency was installed or imported. New requirements file explicitly includes the existing pinned evaluation dependencies.

The synthetic runner/validator checks use generated examples, not benchmark labels. These tests do not establish efficacy, real-data convergence of every new arm, or completion of an independent Claude review.

## Not done

Independent review, exact proposal approval and freeze, real-data execution/validation, all-layer extraction, additional encoder seeds and reactivity pretraining.
