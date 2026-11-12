"""Extract label-free descriptor blocks and assemble aligned distinct-position queries."""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from regjepa.feature_pipeline import STAGES, prepare, geometry_phase, morgan_phase, electronic_phase, jepa_phase, assemble


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=["geometry", "morgan", "electronic", "jepa", "random", "assemble", "all"], default="all")
    parser.add_argument("--cohort", choices=STAGES, default="all_source_records")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--device", default="cpu")
    args = parser.parse_args()
    config = json.loads((ROOT / "configs/features.json").read_text(encoding="utf-8"))
    state = prepare(ROOT, config, args.cohort, args.limit)
    print(f"Run {state['run_key']}: {len(state['records'])} records, {len(state['structures'])} structures, phase={args.phase}", flush=True)
    if args.phase in ["geometry", "all"]:
        geometry_phase(state, args.workers)
    if args.phase in ["morgan", "all"]:
        morgan_phase(state)
    if args.phase in ["electronic", "all"]:
        electronic_phase(state, args.workers)
    if args.phase in ["jepa", "all"]:
        jepa_phase(state, device=args.device)
    if args.phase in ["random", "all"]:
        jepa_phase(state, random_init=True, device=args.device)
    if args.phase in ["assemble", "all"]:
        assemble(state)


if __name__ == "__main__":
    main()
