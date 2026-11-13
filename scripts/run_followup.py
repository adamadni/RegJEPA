"""Inspect a follow-up proposal by default; fitting requires --run and a frozen artifact."""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from regjepa_followup.experiment import run, workload
from regjepa_followup.protocol import load


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    plan = load(ROOT, args.splits.resolve(), frozen=args.run)
    if not args.run:
        print(json.dumps({"proposal_id": plan["proposal_id"], "status": plan["status"],
                          "models_fitted": 0, "workload": workload(plan)}, indent=2))
        return
    output = args.output or ROOT / "reports/followup_runs" / plan["proposal_id"]
    result = run(ROOT, args.splits.resolve(), output.resolve())
    print(json.dumps({"proposal_id": plan["proposal_id"], "completed_final_models": result["completed_final_models"]}))


if __name__ == "__main__":
    main()
