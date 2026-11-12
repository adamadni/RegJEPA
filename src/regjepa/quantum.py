"""Verified xTB GFN1 atomic electronic descriptors on the shared conformers."""
from __future__ import annotations

import json
import os
import re
import subprocess
import shutil
from pathlib import Path

import numpy as np

from .representations import write_xyz

ELECTRONIC_NAMES = ["gfn1_mulliken_charge_e", "gfn1_cm5_charge_e", "gfn1_fukui_fplus",
                    "gfn1_fukui_fminus", "gfn1_fukui_fradical", "gfn1_total_wiberg_order"]
NUMBER = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?"


def _number(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def relax_geometry(executable: Path, arrays: dict, metadata: dict, job_root: Path) -> tuple[dict, dict]:
    """Explicit, recorded GFN2 geometry fallback for unsupported MMFF94/UFF elements."""
    from rdkit import Chem
    outputs, energies, retained_indices, rejected = [], [], [], {}
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    mol = Chem.AddHs(Chem.MolFromSmiles(metadata["canonical_smiles"]))
    sf5_preopt = any(a.GetAtomicNum() == 16 and sum(n.GetAtomicNum() == 9 for n in a.GetNeighbors()) >= 5 for a in mol.GetAtoms())
    for index, positions in enumerate(arrays["positions"]):
        job = job_root / f"conf_{index:03d}"
        job.mkdir(parents=True, exist_ok=True)
        write_xyz(job / "geometry.xyz", arrays["atomic_numbers"], positions)
        if sf5_preopt:
            precommand = [str(executable.resolve()), "geometry.xyz", "--gfnff", "--opt", "normal", "--cycles", "500",
                          "--chrg", str(metadata["formal_charge"])]
            pre = subprocess.run(precommand, cwd=job, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                 text=True, errors="replace", timeout=180,
                                 creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
            (job / "gfnff_preoptimization.txt").write_text(pre.stdout, encoding="utf-8")
            if pre.returncode != 0 or "GEOMETRY OPTIMIZATION CONVERGED" not in pre.stdout:
                rejected[index] = "GFN-FF optimization did not converge"
                continue
            shutil.copyfile(job / "xtbopt.xyz", job / "geometry.xyz")
        command = [str(executable.resolve()), "geometry.xyz", "--gfn", "2", "--opt", "normal", "--cycles", "500",
                   "--chrg", str(metadata["formal_charge"]), "--uhf", str(metadata["unpaired_electrons"]), "--iterations", "500"]
        # GFN2 can change SF5 connectivity during optimization. Retain the validated
        # GFN-FF geometry for these FF-unsupported substrates instead of dissociating them.
        result = pre if sf5_preopt else subprocess.run(
            command, cwd=job, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, errors="replace", timeout=180,
            creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        (job / "stdout.txt").write_text(result.stdout, encoding="utf-8")
        if result.returncode != 0 or "GEOMETRY OPTIMIZATION CONVERGED" not in result.stdout:
            rejected[index] = "xTB geometry optimization did not converge"
            continue
        lines = (job / "xtbopt.xyz").read_text().splitlines()
        symbols = [line.split()[0] for line in lines[2:] if line.strip()]
        expected = [Chem.GetPeriodicTable().GetElementSymbol(int(z)) for z in arrays["atomic_numbers"]]
        if int(lines[0]) != len(expected) or symbols != expected:
            raise ValueError("Fallback optimization changed atom order")
        coords = np.array([[float(value) for value in line.split()[1:4]] for line in lines[2:] if line.strip()])
        if not np.isfinite(coords).all():
            rejected[index] = "Nonfinite optimized coordinates"
            continue
        periodic = Chem.GetPeriodicTable()
        bad_bond = False
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            radius = periodic.GetRcovalent(int(arrays["atomic_numbers"][i])) + periodic.GetRcovalent(int(arrays["atomic_numbers"][j]))
            if not 0.55 * radius < np.linalg.norm(coords[i]-coords[j]) < 1.7 * radius:
                rejected[index] = f"Implausible optimized bond {i}-{j}: {np.linalg.norm(coords[i]-coords[j]):.3f} Angstrom"
                bad_bond = True
                break
        if bad_bond:
            continue
        outputs.append(coords)
        energies.append(float(re.findall(r"TOTAL ENERGY\s+([-\d.]+)", result.stdout)[-1]) * 627.509474)
        retained_indices.append(index)
    if not outputs:
        raise ValueError(f"No graph-preserving converged fallback conformer; rejected={rejected}")
    arrays["positions"] = np.stack(outputs)
    arrays["forcefield_status"] = np.zeros(len(outputs), dtype=np.int16)
    arrays["forcefield_energy_kcal_mol"] = np.array(energies)
    metadata.update({"forcefield": "GFNFF_xTB_fallback" if sf5_preopt else "GFN2_xTB_fallback", "nonconverged_conformers": 0,
                     "fallback_geometry_solvent": "gas phase",
                     "fallback_geometry_method": "GFN-FF native --opt normal" if sf5_preopt else "GFN2-xTB native --opt normal",
                     "requested_conformer_count": len(metadata.get("fallback_retained_input_indices", [])) or metadata["num_conformers"],
                     "num_conformers": len(outputs), "fallback_retained_input_indices": retained_indices,
                     "fallback_rejected_conformers": rejected})
    return arrays, metadata


def parse_xtb_output(text: str, atomic_numbers: np.ndarray, expected_charge: int) -> np.ndarray:
    """Parse named property tables; never use a zero-filled fallback for missing QM."""
    from rdkit import Chem
    if "normal termination of xtb" not in text:
        raise ValueError("xTB did not report normal termination")
    if re.search(r"SCC (?:is )?not converged|SCC convergence failed", text, re.I):
        raise ValueError("xTB SCC convergence failure")
    n = len(atomic_numbers)
    expected = [Chem.GetPeriodicTable().GetElementSymbol(int(z)).lower() for z in atomic_numbers]
    if text.count("Mulliken/CM5 charges") != 1 or text.count("Fukui functions:") != 1:
        raise ValueError("Ambiguous or missing property-table headers")
    charges = text.split("Mulliken/CM5 charges", 1)[1].split("Wiberg/Mayer", 1)[0]
    charge_rows = re.findall(rf"^\s*(\d+)\s*([A-Za-z]+)\s+({NUMBER})\s+({NUMBER})\s+{NUMBER}\s+{NUMBER}\s+{NUMBER}", charges, re.M)
    fukui = text.split("Fukui functions:", 1)[1].split("Property Printout", 1)[0]
    fukui_rows = re.findall(rf"^\s*(\d+)\s*([A-Za-z]+)\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})", fukui, re.M)
    wbo_text = text.split("Wiberg/Mayer (AO) data.", 1)[1]
    wbo_rows = re.findall(rf"^\s*(\d+)\s+(\d+)\s+([A-Za-z]+)\s+({NUMBER})\s+--", wbo_text, re.M)
    if any(len(rows) != n for rows in [charge_rows, fukui_rows, wbo_rows]):
        raise ValueError(f"Atomic table sizes differ: charges={len(charge_rows)}, Fukui={len(fukui_rows)}, WBO={len(wbo_rows)}, expected={n}")
    values = np.empty((n, len(ELECTRONIC_NAMES)), dtype=np.float64)
    for index, (charge, fukui, wbo) in enumerate(zip(charge_rows, fukui_rows, wbo_rows)):
        if int(charge[0]) != index+1 or int(fukui[0]) != index+1 or int(wbo[0]) != index+1:
            raise ValueError("xTB atomic indices do not match XYZ order")
        if charge[1].lower() != expected[index] or fukui[1].lower() != expected[index] or wbo[2].lower() != expected[index]:
            raise ValueError("xTB atomic elements do not match the input")
        if int(wbo[1]) != int(atomic_numbers[index]):
            raise ValueError("xTB atomic number mismatch")
        values[index] = [_number(charge[2]), _number(charge[3]), _number(fukui[2]),
                         _number(fukui[3]), _number(fukui[4]), _number(wbo[3])]
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite xTB descriptors")
    if not np.allclose(values[:, :2].sum(axis=0), expected_charge, atol=2e-3):
        raise ValueError(f"Atomic charges do not sum to molecular charge {expected_charge}")
    if not np.allclose(values[:, 2:5].sum(axis=0), 1.0, atol=max(0.01, 0.0006*n)):
        raise ValueError("Condensed Fukui functions fail the electron-count sum rule")
    return values


def electronic_structure(executable: Path, geometry, metadata: dict, config: dict, count: int, job_root: Path) -> np.ndarray:
    outputs = []
    numbers = geometry["atomic_numbers"]
    n = int(geometry["n_original_atoms"])
    for index, positions in enumerate(geometry["positions"][:count]):
        job = job_root / f"conf_{index:03d}"
        job.mkdir(parents=True, exist_ok=True)
        write_xyz(job / "geometry.xyz", numbers, positions)
        command = [str(executable.resolve()), "geometry.xyz", "--gfn", "1", "--vfukui",
                   "--chrg", str(metadata["formal_charge"]), "--uhf", str(metadata["unpaired_electrons"]),
                   "--acc", str(config["accuracy"]), "--iterations", str(config["max_iterations"])]
        if config["solvent"]:
            command += ["--gbsa", config["solvent"]]
        command_path = job / "command.json"
        command_path.write_text(json.dumps({"command": command, "electronic_temperature_K": 300,
                                           "geometry": "single point on the shared RDKit relaxed conformer"}, indent=2)+"\n", encoding="utf-8")
        log = job / "stdout.txt"
        if log.exists():
            text = log.read_text(encoding="utf-8")
            try:
                outputs.append(parse_xtb_output(text, numbers, metadata["formal_charge"])[:n])
                continue
            except ValueError:
                pass
        env = dict(os.environ, OMP_NUM_THREADS=str(config["threads_per_process"]),
                   OPENBLAS_NUM_THREADS=str(config["threads_per_process"]), MKL_NUM_THREADS=str(config["threads_per_process"]))
        result = subprocess.run(command, cwd=job, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, errors="replace", timeout=config["timeout_seconds"],
                                creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        log.write_text(result.stdout, encoding="utf-8")
        if result.returncode != 0:
            raise RuntimeError(f"xTB exited with {result.returncode}; see {log}")
        outputs.append(parse_xtb_output(result.stdout, numbers, metadata["formal_charge"])[:n])
    if len(outputs) != count:
        raise ValueError("Electronic features do not cover all requested conformers")
    return np.stack(outputs).mean(axis=0).astype(np.float32)
