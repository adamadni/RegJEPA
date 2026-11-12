# RegJEPA

Investigating whether JEPA-derived molecular features improve regioselectivity prediction, especially with limited reaction labels and unfamiliar substrates.

The current milestone is **data acquisition, formatting, structural validation, and symmetry auditing**. No JEPA or regioselectivity predictor has been trained here yet.

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

The downloader is the only workflow step that needs network access when data are absent. It verifies pinned Git blob IDs and records SHA256 checksums. Upstream data and license files are committed verbatim. Generated tables retain the original labels alongside derived symmetry information; they do not overwrite upstream annotations.

## Layout

```text
data/sources.lock.json          Upstream commits, file paths, and blob hashes
data/raw/                      Immutable source snapshots and manifest
data/processed/molecules.csv   Molecule-level audit and cohort flags
data/processed/molecules.jsonl Typed molecule records
data/processed/candidate_sites.csv  Candidate atom and outcome-group table
data/processed/summary.json    Aggregate counts
src/regjepa/                   Reusable ingestion, validation, and symmetry code
scripts/                       Download, audit, and report commands
tests/                         Chemistry and source-integration checks
reports/                       Audit, offline molecule browser, and examples
docs/                          Definitions, pending decisions, and experiment plan
```

## Next experiment

Compare fixed Atom-JEPA atom features, local fingerprints, QM descriptors, and random-backbone features through simple site-scoring heads on identical grouped splits. Primary ranking should operate on distinct substitution outcomes, with each molecule contributing one evaluation unit. See [the first-experiment plan](docs/first_experiment.md).

## Source attribution

The vendored benchmark is from [jensengroup/db-regioselectivity](https://github.com/jensengroup/db-regioselectivity), associated with [Kromann et al., Chemical Science, 2018](https://doi.org/10.1039/C7SC04156J). Additional input structures are from [jensengroup/SI_RegioSQM20](https://github.com/jensengroup/SI_RegioSQM20), associated with [Ree et al., Journal of Cheminformatics, 2021](https://doi.org/10.1186/s13321-021-00490-7). Their MIT licenses and upstream copyright notices are preserved with the corresponding raw data. Repository code licensing has not yet been selected.

