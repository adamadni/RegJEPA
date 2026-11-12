"""Fetch pinned public Atom-JEPA weights and the official Windows xTB binary."""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
import zipfile
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]


def get(url: str, path: Path) -> str:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        request = Request(url, headers={"User-Agent": "RegJEPA/0.2"})
        with urlopen(request, timeout=180) as response, path.with_suffix(path.suffix + ".part").open("wb") as handle:
            shutil.copyfileobj(response, handle)
        path.with_suffix(path.suffix + ".part").replace(path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    print(f"Ready {path.relative_to(ROOT)} ({path.stat().st_size:,} bytes)", flush=True)
    return digest


def main() -> None:
    provenance_path = ROOT / "data/feature_tools.json"
    old = json.loads(provenance_path.read_text(encoding="utf-8")) if provenance_path.exists() else {}
    if "atom_jepa" in old:
        revision = old["atom_jepa"]["revision"]
    else:
        with urlopen(Request("https://huggingface.co/api/models/atom-jepa/atom-jepa", headers={"User-Agent": "RegJEPA/0.2"}), timeout=60) as response:
            revision = json.load(response)["sha"]
    model_files = {}
    for name in ["config.json", "model.safetensors"]:
        url = f"https://huggingface.co/atom-jepa/atom-jepa/resolve/{revision}/molecules/{name}"
        path = ROOT / ".cache/models/atom-jepa/molecules" / name
        digest = get(url, path)
        previous = old.get("atom_jepa", {}).get("files", {}).get(name)
        if previous and previous["sha256"] != digest:
            raise ValueError(f"Cached model checksum mismatch: {name}")
        model_files[name] = {"url": url, "sha256": digest, "bytes": path.stat().st_size}
    data = {"atom_jepa": {"repository": "atom-jepa/atom-jepa", "revision": revision,
                          "checkpoint": "molecules", "package_version": "0.4.0", "files": model_files}}
    if sys.platform == "win32":
        base = "https://github.com/grimme-lab/xtb/releases/download/v6.7.1/xtb-6.7.1pre-windows-x86_64.zip"
        archive = ROOT / ".cache/tools/xtb-6.7.1pre-windows-x86_64.zip"
        get(base + ".sha256", archive.with_suffix(".zip.sha256"))
        expected = archive.with_suffix(".zip.sha256").read_text().split()[0]
        digest = get(base, archive)
        if digest != expected:
            raise ValueError("xTB release archive failed its publisher SHA256 check")
        target = ROOT / ".cache/tools/xtb"
        target.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(archive) as z:
            for member in z.infolist():
                destination = (target / member.filename).resolve()
                if not destination.is_relative_to(target.resolve()):
                    raise ValueError("Archive path escapes the tool directory")
                if not destination.exists():
                    z.extract(member, target)
        executables = list(target.rglob("xtb.exe"))
        if len(executables) != 1:
            raise ValueError(f"Expected one xtb.exe, found {executables}")
        data["xtb"] = {"release": "v6.7.1", "binary_version": "6.7.1pre Windows",
                       "archive_url": base, "archive_sha256": digest,
                       "executable": executables[0].relative_to(ROOT).as_posix()}
        print(f"xTB executable: {executables[0]}", flush=True)
    elif "xtb" in old:
        data["xtb"] = old["xtb"]
    provenance_path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
