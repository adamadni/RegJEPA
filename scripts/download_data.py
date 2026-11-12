"""Download immutable upstream files, verifying Git blob hashes and SHA256 checksums."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]


def git_blob_sha(data: bytes) -> str:
    return hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    lock = json.loads((ROOT / "data/sources.lock.json").read_text(encoding="utf-8"))
    manifest_path = ROOT / "data/raw/manifest.json"
    previous = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else {}
    entries = []
    for source in lock["sources"]:
        for relative, expected_blob in source["files"].items():
            path = ROOT / source["destination"] / relative
            url = f"https://raw.githubusercontent.com/{source['repository']}/{source['commit']}/{relative}"
            if not path.exists():
                if args.verify_only:
                    raise FileNotFoundError(path)
                request = Request(url, headers={"User-Agent": "RegJEPA-data-audit/0.1"})
                with urlopen(request, timeout=60) as response:
                    data = response.read()
                if git_blob_sha(data) != expected_blob:
                    raise ValueError(f"Upstream Git blob checksum mismatch: {url}")
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(data)
            data = path.read_bytes()
            if git_blob_sha(data) != expected_blob:
                raise ValueError(f"Local Git blob checksum mismatch: {path}")
            entries.append({"path": path.relative_to(ROOT).as_posix(), "url": url,
                            "git_blob_sha1": expected_blob, "sha256": hashlib.sha256(data).hexdigest(),
                            "bytes": len(data), "citation": source["citation"]})
            print(f"Verified {path.relative_to(ROOT)} ({len(data):,} bytes)")
    manifest = {"schema_version": 1,
                "retrieved_at_utc": previous.get("retrieved_at_utc", datetime.now(timezone.utc).isoformat()),
                "files": entries}
    if args.verify_only:
        if previous.get("files") != entries:
            raise ValueError("Manifest does not match files; rerun download_data.py")
    else:
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
