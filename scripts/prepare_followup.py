"""Prepare follow-up review artifacts or explicitly freeze a reviewed proposal; never train."""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from regjepa.feature_io import load_feature_bank
from regjepa.protocol import write_json
from regjepa_followup.experiment import workload
from regjepa_followup.protocol import freeze, prepare


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/followup.json")
    parser.add_argument("--freeze", type=Path)
    args = parser.parse_args()
    if args.freeze:
        print(freeze(ROOT, args.freeze.resolve()))
        return
    config = json.loads(args.config.read_text())
    plan = prepare(ROOT, config)
    bank = load_feature_bank(ROOT, cohort=config["cohort"], blocks=["electronic"])
    values = bank.X["electronic"]
    import numpy as np
    residual = values[:, 4] - .5 * (values[:, 2] + values[:, 3])
    report = {"proposal_id": plan["proposal_id"], "status": "proposal_not_frozen_no_real_data_models_fitted",
              "source_records": len(set(bank.compound_ids)), "structures": len({q['structure_key'] for q in bank.queries}),
              "outcome_rows": len(bank.queries), "leakage_groups": len(set(plan["compound_to_group"].values())),
              "identical_original_fold_and_subset_schedule": True, "feature_bank_unchanged": True,
              "fzero_minus_half_fplus_fminus_max_abs": float(np.max(np.abs(residual))),
              "workload": workload(plan)}
    directory = ROOT / "data/followup_splits" / plan["proposal_id"]
    write_json(directory / "proposal.json", plan, exclusive=True)
    write_json(ROOT / "reports/followup_preflight.json", report)
    lines = ["# Exploratory follow-up preflight", "", f"Proposal: `{plan['proposal_id']}`; **not frozen, no real-data fits**.", "",
             f"Observed: {report['source_records']} records, {report['structures']} structures, {report['outcome_rows']} outcome rows, {report['leakage_groups']} groups.",
             "The complete original grouped fold/subset schedule and existing bank are preserved. This is an adaptive follow-up on already evaluated molecules, not a fresh confirmatory test.", "",
             f"Planned workload: {report['workload']['inner_models']:,} inner models and {report['workload']['final_models']:,} final models; including MLP restarts, **{report['workload']['inner_optimizer_attempts']:,} inner and {report['workload']['final_optimizer_attempts']:,} final optimizer attempts**.", "",
             "| Arm | Head | Columns | Outer settings | Penalty choices | Inner optimizer attempts |", "|---|---|---:|---:|---:|---:|"]
    for a in report["workload"]["arms"]:
        lines.append(f"| {a['id']} | {a['head']} | {a['width']} | {a['settings']} | {a['penalty_candidates']} | {a['inner_optimizer_attempts']} |")
    lines += ["", f"Maximum native-precision residual f0 - (f+ + f-)/2: {report['fzero_minus_half_fplus_fminus_max_abs']:.8g}. Fukui channels are not independent; leave-one-out scores are not unique or causal importance.",
              "", "Read docs/followup_protocol.md and reports/claude_followup_review_prompt.txt before freezing.",
              "", "## Not done", "", "Independent protocol/code review, explicit freeze, real-data execution and all-layer extraction.", ""]
    (ROOT / "reports/followup_preflight.md").write_text("\n".join(lines), encoding="utf-8", newline="\n")
    print(json.dumps({k: v for k, v in report.items() if k != "workload"}, indent=2))


if __name__ == "__main__":
    main()
