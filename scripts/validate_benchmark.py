"""Run the deterministic molecular audit; raw source files are never rewritten."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from regjepa.audit import run

if __name__ == "__main__":
    print(json.dumps(run(ROOT), indent=2))

