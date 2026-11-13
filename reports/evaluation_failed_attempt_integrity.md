# Failed-attempt artifact integrity

The original frozen run `09ebadea135ec50dc873` contains **146** complete JSON/model pairs. Every saved model's SHA256 agrees with its JSON record. All **3,168** recorded inner fits and **146** final fits satisfy the original convergence rule. The largest saved iteration count is **940**; the largest saved gradient infinity norm is **4.436536773287124e-7**. No completed-study Markdown or JSON summary exists.

| Block | Saved final fits | Planned final fits |
|---|---:|---:|
| Atom-JEPA | 35 | 35 |
| Random backbone | 35 | 35 |
| Morgan | 35 | 35 |
| Electronic | 35 | 35 |
| JEPA + QM | 6 | 35 |
| Total | 146 | 175 |

The proposed replacement `af05d8348e60ad38b6e4` differs in exactly one configuration leaf: `optimizer.max_iterations`, from **1,000** to **5,000**. All other proposal fields, apart from its ID and embedded configuration, compare equal. This includes feature/data/source input locks, software versions, leakage groups, all outer/inner folds and every training subset. The grouping audit files are byte-identical. The replacement is not frozen.

The audit checked hashes and optimization diagnostics without inspecting predictive metrics. The machine-readable sibling `evaluation_failed_attempt_integrity.json` contains the per-model hash checks. Reviewed source, scripts, original configuration, raw data and feature bank also have no Git diff from `977f120`.

## Not done

This audit does not evaluate predictive performance or complete the interrupted study. The cap-only amendment awaits explicit approval.
