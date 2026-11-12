"""Label-free molecular geometry and descriptor operations for distinct outcome queries."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdFingerprintGenerator


def digest_json(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def canonical_structure(record: dict) -> tuple[Chem.Mol, dict[int, int]]:
    """Map original source atoms to canonical-graph atoms before adding any hydrogens."""
    original = Chem.MolFromSmiles(record["source_smiles"])
    canonical = Chem.MolFromSmiles(record["canonical_isomeric_smiles"])
    if original is None or canonical is None:
        raise ValueError("Invalid structure")
    canonical_to_original = original.GetSubstructMatch(canonical, useChirality=True)
    if len(canonical_to_original) != original.GetNumAtoms():
        raise ValueError(f"No complete stereo-preserving atom correspondence for {record['compound_id']}")
    return canonical, {source: index for index, source in enumerate(canonical_to_original)}


def canonical_group_indices(record: dict) -> list[list[int]]:
    _, mapping = canonical_structure(record)
    return [sorted(mapping[i] for i in group) for group in record["candidate_outcome_groups"]]


def group_features(atom_features: np.ndarray, groups: list[list[int]]) -> np.ndarray:
    """One equally weighted feature vector per distinct outcome, without multiplicity."""
    values = np.asarray(atom_features)
    if values.ndim != 2 or not np.isfinite(values).all():
        raise ValueError("Atom features must be a finite [atoms, features] matrix")
    if any(not g or min(g) < 0 or max(g) >= len(values) for g in groups):
        raise ValueError("Invalid outcome-group atom indices")
    return np.stack([values[group].mean(axis=0) for group in groups]).astype(np.float32)


def generate_geometry(smiles: str, config: dict) -> tuple[dict[str, np.ndarray], dict]:
    """Explicit-H ETKDGv3, fixed seed, no RMS pruning, MMFF94 up to 500 iterations.

    Partial embedding is retried with random starting coordinates as in the paper.
    UFF, when enabled, is a recorded deviation for molecules without MMFF94 parameters.
    Nonconverged relaxed conformers are retained and flagged, matching upstream behavior.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("A single sanitized molecular component is required")
    n_original = mol.GetNumAtoms()
    original_numbers = [a.GetAtomicNum() for a in mol.GetAtoms()]
    mol = Chem.AddHs(mol)
    assert [a.GetAtomicNum() for a in mol.GetAtoms()][:n_original] == original_numbers
    params = AllChem.ETKDGv3()
    params.randomSeed = config["seed"]
    params.pruneRmsThresh = -1.0
    params.enforceChirality = True
    params.numThreads = 1
    params.timeout = 30
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=config["num_conformers"], params=params))
    fallback = len(ids) != config["num_conformers"]
    if fallback:
        params.useRandomCoords = True
        ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=config["num_conformers"], params=params))
    if len(ids) != config["num_conformers"]:
        raise RuntimeError(f"Embedded only {len(ids)} of {config['num_conformers']} requested conformers")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        forcefield = "MMFF94"
        status_energy = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1, maxIters=config["mmff_max_iterations"], mmffVariant="MMFF94")
    elif config.get("allow_uff_fallback") and AllChem.UFFHasAllMoleculeParams(mol):
        forcefield = "UFF_fallback"
        status_energy = AllChem.UFFOptimizeMoleculeConfs(mol, numThreads=1, maxIters=config["mmff_max_iterations"])
    elif config.get("allow_xtb_fallback"):
        forcefield = "GFN2_xTB_pending"
        status_energy = [(-1, float("nan")) for _ in ids]
    else:
        raise ValueError("MMFF94 parameters unavailable and no permitted fully parameterized UFF fallback")
    coords = np.stack([mol.GetConformer(i).GetPositions() for i in ids])
    if not np.isfinite(coords).all():
        raise ValueError("Nonfinite conformer coordinates")
    original = Chem.RemoveHs(Chem.Mol(mol))
    if Chem.MolToSmiles(original, isomericSmiles=True) != Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True):
        raise ValueError("Hydrogen addition/removal changed the specified graph or stereochemistry")
    arrays = {"atomic_numbers": np.array([a.GetAtomicNum() for a in mol.GetAtoms()], dtype=np.int64),
              "positions": coords, "n_original_atoms": np.array(n_original),
              "forcefield_status": np.array([s for s, _ in status_energy], dtype=np.int16),
              "forcefield_energy_kcal_mol": np.array([e for _, e in status_energy])}
    metadata = {"canonical_smiles": smiles, "rdkit_version": rdBase.rdkitVersion,
                "config": config, "forcefield": forcefield, "random_coordinate_fallback": fallback,
                "num_conformers": len(ids), "n_original_atoms": n_original,
                "n_atoms_with_H": mol.GetNumAtoms(), "formal_charge": Chem.GetFormalCharge(mol),
                "unpaired_electrons": sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()),
                "nonconverged_conformers": sum(s != 0 for s, _ in status_energy),
                "n_unassigned_chiral_centers": sum(label == "?" for _, label in Chem.FindMolChiralCenters(mol, includeUnassigned=True))}
    return arrays, metadata


def morgan_atom_features(smiles: str, config: dict) -> np.ndarray:
    """Rooted circular bit fingerprints, with radii kept in separate feature blocks."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")
    width = config["bits_per_radius"]
    values = np.zeros((mol.GetNumAtoms(), width * len(config["radii"])), dtype=np.float32)
    for j, radius in enumerate(config["radii"]):
        generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=width,
                                                              includeChirality=config["include_chirality"])
        for index in range(mol.GetNumAtoms()):
            fp = generator.GetFingerprintAsNumPy(mol, fromAtoms=[index])
            values[index, j*width:(j+1)*width] = fp
    return values


def write_xyz(path: Path, atomic_numbers: np.ndarray, positions: np.ndarray) -> None:
    periodic = Chem.GetPeriodicTable()
    lines = [str(len(atomic_numbers)), "RegJEPA canonical atom order, Angstrom"]
    lines.extend(f"{periodic.GetElementSymbol(int(z))} {x:.12f} {y:.12f} {zcoord:.12f}"
                 for z, (x, y, zcoord) in zip(atomic_numbers, positions))
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
