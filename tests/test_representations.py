import copy
import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from regjepa.feature_pipeline import STAGES, cohort_memberships
from regjepa.quantum import parse_xtb_output
from regjepa.representations import canonical_group_indices, canonical_structure, generate_geometry, group_features, morgan_atom_features

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def records():
    return {r["compound_id"]: r for r in map(json.loads, (ROOT / "data/processed/molecules.jsonl").read_text().splitlines())}


def test_staged_cohorts_are_nested(records):
    sets = [{key for key, record in records.items() if stage in cohort_memberships(record)} for stage in STAGES]
    assert all(a <= b for a, b in zip(sets, sets[1:]))
    assert len(sets[0]) == 299
    assert len(sets[1]) > len(sets[0])
    assert len(sets[2]) == 495
    assert len(sets[3]) == 535


def test_equivalent_atom_features_yield_only_one_query():
    atoms = np.array([[1., 3.], [3., 1.], [8., 9.]])
    got = group_features(atoms, [[0, 1], [2]])
    np.testing.assert_array_equal(got, [[2., 2.], [8., 9.]])
    assert len(got) == 2


def test_feature_extraction_does_not_read_outcome_labels(records):
    source = records["comp18"]
    changed = copy.deepcopy(source)
    changed["observed_atom_indices"] = [1]
    changed["observed_outcome_group_indices"] = [0]
    assert canonical_group_indices(source) == canonical_group_indices(changed)
    config = {"radii": [1, 2, 3], "bits_per_radius": 128, "include_chirality": True}
    np.testing.assert_array_equal(morgan_atom_features(source["canonical_isomeric_smiles"], config),
                                  morgan_atom_features(changed["canonical_isomeric_smiles"], config))


def test_duplicate_with_different_source_atom_order_gives_identical_queries(records):
    config = {"radii": [1, 2, 3], "bits_per_radius": 128, "include_chirality": True}
    fingerprints = morgan_atom_features(records["comp13"]["canonical_isomeric_smiles"], config)
    a = group_features(fingerprints, canonical_group_indices(records["comp13"]))
    b = group_features(fingerprints, canonical_group_indices(records["comp598"]))
    assert sorted(a.tobytes() for a in a) == sorted(b.tobytes() for b in b)


def test_rooted_fingerprints_distinguish_anisole_positions():
    config = {"radii": [1, 2, 3], "bits_per_radius": 2048, "include_chirality": True}
    got = morgan_atom_features("COc1ccccc1", config)
    assert got.shape == (8, 6144)
    np.testing.assert_array_equal(got[3], got[7])
    np.testing.assert_array_equal(got[4], got[6])
    assert not np.array_equal(got[3], got[5])


def test_fixed_seed_conformers_are_repeatable_and_keep_source_atoms():
    config = {"num_conformers": 2, "seed": 42, "mmff_max_iterations": 500, "allow_uff_fallback": False}
    a, metadata = generate_geometry("N[C@H](C)c1ccccc1", config)
    b, _ = generate_geometry("N[C@H](C)c1ccccc1", config)
    np.testing.assert_array_equal(a["atomic_numbers"], b["atomic_numbers"])
    np.testing.assert_allclose(a["positions"], b["positions"], atol=1e-10)
    assert metadata["forcefield"] == "MMFF94"
    assert metadata["num_conformers"] == 2
    assert a["positions"].shape[0] == 2
    expected = [atom.GetAtomicNum() for atom in Chem.MolFromSmiles("N[C@H](C)c1ccccc1").GetAtoms()]
    assert list(a["atomic_numbers"][:len(expected)]) == expected


def test_real_xtb_properties_parse_and_obey_charge_sum():
    text = (ROOT / "tests/fixtures/xtb_gfn1_anisole_properties.txt").read_text()
    mol = Chem.AddHs(Chem.MolFromSmiles("COc1ccccc1"))
    numbers = np.array([a.GetAtomicNum() for a in mol.GetAtoms()])
    got = parse_xtb_output(text, numbers, expected_charge=0)
    assert got.shape == (16, 6)
    np.testing.assert_allclose(got[0, :2], [0.16048, 0.03921])
    np.testing.assert_allclose(got[:, :2].sum(axis=0), 0, atol=1e-4)
    np.testing.assert_allclose(got[:, 2:5].sum(axis=0), 1, atol=0.003)


@pytest.mark.parametrize("corruption", ["no_termination", "wrong_elements", "missing_charge_row", "nonneutral_charges"])
def test_corrupted_qm_output_is_rejected(corruption):
    text = (ROOT / "tests/fixtures/xtb_gfn1_anisole_properties.txt").read_text()
    mol = Chem.AddHs(Chem.MolFromSmiles("COc1ccccc1"))
    numbers = np.array([a.GetAtomicNum() for a in mol.GetAtoms()])
    if corruption == "no_termination":
        text = text.replace("normal termination of xtb", "failed")
    elif corruption == "wrong_elements":
        numbers[0] = 7
    elif corruption == "missing_charge_row":
        text = text.replace("     1C     0.16048  0.03921   0.975  2.864  0.000", "")
    else:
        text = text.replace("0.16048", "1.16048")
    with pytest.raises(ValueError):
        parse_xtb_output(text, numbers, expected_charge=0)
