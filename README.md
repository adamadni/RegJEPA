# RegJEPA

Investigating whether frozen JEPA-derived molecular features improve prediction of reported EAS regioselectivity with limited reaction labels.

The current milestone is a **completed frozen-feature linear pilot on 299 conservative bromination records**. We compare a published structure-only Atom-JEPA encoder, a randomly initialized frozen copy, atom-centered Morgan fingerprints, electronic descriptors, and JEPA + QM. All 6,300 inner fits and 175 final fits converged; every saved model and score passed replay validation.

Electronic descriptors are the strongest single block. Adding JEPA gives small positive top-1 point estimates, but its gain over QM and pretrained JEPA's gain over the random backbone remain unresolved by the exploratory paired intervals. Read the [results and limits](reports/evaluation_results.md), [numerical audit](reports/evaluation_result_validation.md), and [Claude results-review prompt](reports/claude_results_review_prompt.txt).

The [reporting-review response](reports/evaluation_reporting_audit.md) adds penalty-edge direction, Morgan tie causes, per-molecule negative information gain and candidate-count breakdowns without refitting. At 50 labels, random selects the strongest tested penalty in 7 of 15 fits; that bounded search limits interpretation of the small pretrained-versus-random point difference.

The next [descriptor/hybrid/readout protocol](docs/followup_protocol.md) is implemented and prepared as proposal `1ddb9c393580ec2e232f`, **not frozen or run**. It covers QM ablations, QM+random/Morgan controls, scalar-only features and small nonlinear readouts. See the [exact workload](reports/followup_preflight.md), [preparation checks](reports/followup_preparation_validation.md), and [Claude review prompt](reports/claude_followup_review_prompt.txt) before authorizing the larger experiment.

See the [feature extraction report](reports/feature_extraction.md) and `data/features/latest.json` for the aligned feature bank, geometry archive, checkpoint provenance and coverage.

## Initial data decision

- 535 EAS source records, representing 534 distinct isomeric substrates.
- Every source structure parses and every source site label identifies an aromatic C-H under zero-based indexing.
- The original reagent/yield notes contain bromination, chlorination, iodination, nitration, and pooled experiments. The data must not be treated as 535 uniform-condition brominations.
- 166 records have equivalent candidate sites; 16 have equivalent observed-site counterparts omitted from their raw labels.
- A conservative, rule-filtered bromination cohort contains 299 records. It is suitable for an exploratory frozen-feature benchmark, subject to the documented limitations.

Read the [quantitative audit](reports/data_audit.md), [decision record](docs/benchmark_decisions.md), and [data dictionary](docs/data_dictionary.md). Download and open [the offline molecule browser](reports/molecule_browser.html) in a browser to inspect all structures and annotations. GitHub's file page does not execute the browser.

![Selected symmetry and annotation cases](reports/symmetry_examples.png)

Green marks source-annotated sites; blue marks additional equivalent sites. Atom indices follow the original zero-based SMILES atom order. Br-marker products in the tables are hypothetical symmetry probes, not asserted observed products.

## Reproduce

Python 3.11+; the initial audit used Python 3.12.14 and RDKit 2026.09.2. Commands below assume a virtual environment is active; on Windows its interpreter can also be invoked directly as `.venv\Scripts\python.exe`.

```text
python -m venv .venv
python -m pip install -r requirements.txt
python scripts/download_data.py
python scripts/download_data.py --verify-only
python scripts/validate_benchmark.py
python scripts/build_reports.py
python -m pytest -q
```

The data downloader needs network access when raw data are absent. It verifies pinned Git blob IDs and records SHA256 checksums. Upstream data and license files are committed verbatim. Generated tables retain the original labels alongside derived symmetry information; they do not overwrite upstream annotations.

`requirements.lock.txt` additionally pins the full dependency set used for this audit; install it instead of `requirements.txt` for that exact environment. `.gitattributes` disables line-ending conversion for raw snapshots and hashed feature exports so their checksums remain valid on Windows and Linux.

## Layout

```text
data/sources.lock.json          Upstream commits, file paths, and blob hashes
data/raw/                      Immutable source snapshots and manifest
data/processed/molecules.csv   Molecule-level audit and cohort flags
data/processed/molecules.jsonl Typed molecule records
data/processed/candidate_sites.csv  Candidate atom and outcome-group table
data/processed/summary.json    Aggregate counts
data/features/                Aligned frozen feature blocks and portable 3D geometry
data/splits/                  Immutable proposals, approvals and frozen split bindings
src/regjepa/                   Ingestion, chemistry, features and linear evaluation
scripts/                       Download, audit, and report commands
tests/                         Chemistry and source-integration checks
reports/                       Data audits, molecule browser and evaluation results
docs/                          Definitions, pending decisions, and experiment plan
```

## Frozen features

```text
python -m pip install -r requirements-features.txt
python scripts/setup_feature_tools.py
python scripts/extract_features.py --phase all --workers 2
python scripts/summarize_features.py
```

The setup step downloads the public molecular checkpoint at a pinned revision and, on Windows, the official checksum-verified xTB binary. Large weights, binaries and scratch calculations are cached in ignored `.cache/`. Dependencies are fully recorded in `requirements-features.lock.txt`.

The completed first bank uses one shared conformer for all 3D feature blocks. Up to 10 relaxed conformers are generated and saved per substrate, enabling a separate ensemble experiment. MMFF94 is the main relaxation method; UFF or xTB fallbacks and convergence/geometry rejections are explicitly recorded for unsupported molecules.

Each source record contributes one row per symmetrically distinct outcome. Member-atom features are averaged into one vector. Group multiplicity is not a probability weight or an additional example. The five blocks are Atom-JEPA (768 columns), random frozen backbone (768), atom-centered Morgan radii 1/2/3 (6144), xTB electronic descriptors (6), and JEPA + QM concatenation (774). Arrays are unscaled; downstream scaling belongs inside training folds.

Load an aligned cohort through `regjepa.feature_io` after installing the project (`python -m pip install -e .`), or with `src` on `PYTHONPATH`:

```python
from pathlib import Path
from regjepa.feature_io import load_feature_bank

bank = load_feature_bank(Path.cwd(), cohort="conservative_bromination")
X = bank.X["atom_jepa"]  # one row per distinct-position query
recorded_sites = bank.observed_membership
structure_keys = [row["structure_key"] for row in bank.queries]
```

The loader verifies row alignment and rejects incomplete candidate sets. Use `structure_key` from `bank.queries` to keep duplicate substrates together when assigning supervised folds. See the [feature-bank dictionary](docs/feature_bank.md) for array and geometry conventions.

## Approved progression and next experiment

The cumulative cohorts are conservative bromination (299 records), conservative halogenation including Cl/I (374), broad single-outcome (495), and all source records (535). All cohorts use the same feature definitions and query identifiers. Multi-outcome and no-choice records remain flagged in the final stage.

Compare these frozen features through simple site-scoring heads on identical grouped splits. Each molecule contributes one evaluation unit. See [the first-experiment plan](docs/first_experiment.md).

## Linear evaluation harness

The reviewed linear scorer and nested evaluation harness have completed proposal `af05d8348e60ad38b6e4`. Read the [evaluation protocol](docs/evaluation_protocol.md), [completed results](reports/evaluation_results.md), and [current results-review prompt](reports/claude_results_review_prompt.txt). The [preflight](reports/evaluation_preflight.md) and original review prompt are historical preparation records.

```text
python -m pip install -r requirements-evaluation.txt
python scripts/prepare_evaluation.py
```

Preparation writes an immutable, label-blind split **proposal** and a grouping audit under `data/splits/`, plus Markdown/JSON preflight reports. All 535 records participate in linkage before filtering the 299-record cohort. Exact identity, tautomer, explicit source/series links and whole-molecule Morgan Tanimoto >=0.6 neighbours stay together. This is a conservative metadata constraint, not verified international patent-family recovery.

After reviewing the proposal, use the explicit freeze and run commands in the protocol. The runner inspects artifacts without fitting unless `--run` is supplied, and fitting rejects unfrozen or changed inputs. Legacy raw-atom scoring helpers in `chemistry.py` are not used on feature-bank rows; the harness has independent exact-tie-aware metrics.

The original 1,000-iteration attempt remains separately preserved after a convergence failure. The user-approved [cap-only amendment](reports/evaluation_iteration_amendment.md) raised the cap to 5,000; the fresh run completed without relaxing convergence criteria. See the [split index](data/splits/README.md), exact approval alongside the amended frozen artifact, and [runtime record](reports/evaluation_amended_runtime.json). This split holds out known source/series/similarity components; it is not a verified independent-publication-family or out-of-distribution evaluation.

To audit the committed results without fitting, use the exact locked runtime and a new report path:

```text
python scripts/validate_evaluation_results.py --splits data/splits/af05d8348e60ad38b6e4/frozen.json --output reports/evaluation_runs/af05d8348e60ad38b6e4 --report .cache/replay-audit.json
```

A fresh model rerun must specify a new `--output` directory: the committed experiment outputs are never overwritten. Validation reuses the implementation and checks saved artifacts; it is not an independent scientific replication.

For tests with a workspace-local temporary directory, create its parent first (also works in a fresh clone):

```text
python -c "from pathlib import Path; Path('.cache').mkdir(exist_ok=True)"
python -m pytest -q --basetemp .cache/pytest-evaluation
```

## Source attribution

The vendored benchmark is from [jensengroup/db-regioselectivity](https://github.com/jensengroup/db-regioselectivity), associated with [Kromann et al., Chemical Science, 2018](https://doi.org/10.1039/C7SC04156J). Additional input structures are from [jensengroup/SI_RegioSQM20](https://github.com/jensengroup/SI_RegioSQM20), associated with [Ree et al., Journal of Cheminformatics, 2021](https://doi.org/10.1186/s13321-021-00490-7). Their MIT licenses and upstream copyright notices are preserved with the corresponding raw data. Repository code licensing has not yet been selected.

