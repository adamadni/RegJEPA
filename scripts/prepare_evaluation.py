"""Prepare grouping/split review artifacts, or explicitly freeze a reviewed proposal; never fit models."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from regjepa.experiment import penalty_candidates
from regjepa.protocol import freeze_protocol, prepare_protocol, write_json


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/evaluation.json")
    parser.add_argument("--freeze", type=Path, help="Explicitly freeze this existing, reviewed proposal")
    args = parser.parse_args()
    if args.freeze:
        print(freeze_protocol(ROOT, args.freeze.resolve()))
        return
    config = json.loads(args.config.read_text(encoding="utf-8"))
    plan, grouping, summary = prepare_protocol(ROOT, config)
    output = ROOT / "data/splits" / plan["proposal_id"]
    write_json(output / "proposal.json", plan, exclusive=True)
    write_json(output / "grouping.json", grouping, exclusive=True)
    summary["proposal_path"] = (output / "proposal.json").relative_to(ROOT).as_posix()
    training_sets = sum(len(f["training_sets"]) for f in plan["outer_folds"])
    candidates = {b: len(penalty_candidates(b, config)) for b in config["blocks"]}
    summary["planned_inner_fits"] = training_sets * config["inner_folds"] * sum(candidates.values())
    summary["planned_final_fits"] = training_sets * len(config["blocks"])
    summary["penalty_candidates_per_block"] = candidates
    write_json(ROOT / "reports/evaluation_preflight_summary.json", summary)
    lines = ["# Evaluation preflight: grouping and split proposal", "",
             "No real-data scorer has been fitted. Splits are proposals, not frozen experiment assignments.", "",
             f"Proposal: `{plan['proposal_id']}`. File: `{summary['proposal_path']}`.", "",
             f"Cohort: {summary['molecules']} records / {summary['structures']} structures / {summary['rows']} outcome queries.",
             f"Uniform expected top-1: {summary['analytical_uniform_metrics']['top1']:.10f}; first-row tie artifact: {summary['constant_first_row_top1']:.10f}.", "",
             f"Identity + explicit source/series + whole-molecule similarity components: **{summary['leakage_groups']}**; largest selected component: **{summary['largest_leakage_group']}** records.",
             f"Maximum Tanimoto between different components: **{summary['maximum_between_group_tanimoto']:.6f}**, below the fixed {config['fingerprint']['similarity_threshold']} boundary.", "",
             "All 535 source records participate in linkage before cohort filtering. Pointer tokens and named targets are joined explicitly, including mixed DOI/patent + pointer records. Suffix aliases and international patent families are not guessed.", "",
             f"Unverified bibliographic identities remain for **{len(summary['unresolved_bibliographic_records'])}** selected records. They are listed in the JSON sibling; known explicit links are constrained, but this is not a verified publication-family holdout.", "",
             "| Outer fold | Training records | Test records |", "|---|---:|---:|"]
    lines.extend(f"| {f['fold']} | {f['train']} | {f['test']} |" for f in summary["outer_counts"])
    lines.extend(["", "Candidate budgets are nested within each subset seed and contain whole leakage components. All tuning labels are inside the stated budget; full training is fitted once per fold.",
                  f"Actual label counts: `{summary['actual_learning_curve_labels']}`.", "",
                  "Tie handling: exact ties receive uniform expected top-1/top-2 and expected reciprocal rank. No raw-atom helper is used on bank rows.",
                  "Linear scaling: subtract each molecule's candidate-feature mean, then scale dense columns using equal-molecule, training-only within-molecule RMS. Morgan columns remain unscaled. Columns without training variation are dropped. No intercept is fitted.",
                  "JEPA + QM has independently tuned block penalties (36 combinations); each other block has 6 candidates. This larger search is reported and contained inside inner folds.",
                  f"Planned workload: {summary['planned_inner_fits']} inner fits plus {summary['planned_final_fits']} final fits. These have not been run.", "",
                  "The primary cohort retains all 299 records. Seven geometry-fallback cases are flagged, not removed to improve a metric. No geometry sensitivity model is fitted in this preflight.", "",
                  "Read docs/evaluation_protocol.md and reports/claude_review_prompt.txt before freezing or running the experiment.", ""])
    (ROOT / "reports/evaluation_preflight.md").write_text("\n".join(lines), encoding="utf-8", newline="\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
