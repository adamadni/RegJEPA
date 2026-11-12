"""Inspect a split artifact by default; --run requires frozen splits and fits only linear readouts."""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from regjepa.experiment import run_experiment
from regjepa.protocol import load_split_plan


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--run", action="store_true", help="Fit the full linear experiment; requires frozen splits")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    path = args.splits.resolve()
    plan = load_split_plan(ROOT, path, require_frozen=args.run)
    if not args.run:
        print(json.dumps({"proposal_id": plan["proposal_id"], "status": plan["status"],
                          "outer_folds": len(plan["outer_folds"]), "models_fitted": 0}, indent=2))
        return
    output = args.output or ROOT / "reports/evaluation_runs" / plan["proposal_id"]
    summary = run_experiment(ROOT, path, output.resolve())
    print(json.dumps({"completed_final_fits": summary["completed_final_fits"], "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
