import csv
import json

import numpy as np
import pytest

from regjepa.feature_io import load_feature_bank


@pytest.fixture
def bank_root(tmp_path):
    directory = tmp_path / "data/features/test"
    directory.mkdir(parents=True)
    (directory.parent / "latest.json").write_text(json.dumps({"manifest": "data/features/test/manifest.json"}))
    (directory / "manifest.json").write_text(json.dumps({"cohorts": {"first": {}, "all": {}}, "blocks": {"a": {}, "b": {}}}))
    rows = []
    for i in range(3):
        rows.append({"row_id": i, "query_id": f"q{i}", "compound_id": "mol1" if i<2 else "mol2",
                     "outcome_group_index": i if i<2 else 0, "source_atom_indices": "[0,1]" if i==0 else "[2]",
                     "canonical_atom_indices": "[0,1]" if i==0 else "[2]", "n_reported_outcomes": 1,
                     "n_candidate_outcomes": 2 if i<2 else 1, "reported_positive": i==0,
                     "cohorts": json.dumps(["first", "all"] if i<2 else ["all"])})
    with (directory / "queries.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    for name in ["a", "b"]:
        np.savez(directory / f"{name}.npz", X=np.arange(6).reshape(3,2).astype(float), present=np.ones(3,dtype=bool), query_ids=np.array([f"q{i}" for i in range(3)]))
    return tmp_path


def test_loader_selects_aligned_cohort_and_preserves_queries(bank_root):
    bank = load_feature_bank(bank_root, cohort="first")
    assert bank.X["a"].shape == bank.X["b"].shape == (2,2)
    assert bank.queries[0]["source_atom_indices"] == [0,1]
    np.testing.assert_array_equal(bank.observed_membership, [True,False])
    np.testing.assert_array_equal(bank.compound_ids, ["mol1","mol1"])


def test_loader_rejects_permuted_feature_rows(bank_root):
    path = bank_root / "data/features/test/b.npz"
    np.savez(path, X=np.zeros((3,2)), present=np.ones(3,dtype=bool), query_ids=np.array(["q1","q0","q2"]))
    with pytest.raises(ValueError,match="alignment"):
        load_feature_bank(bank_root, cohort="first")


def test_loader_never_drops_an_unavailable_candidate_silently(bank_root):
    path = bank_root / "data/features/test/a.npz"
    np.savez(path, X=np.array([[1.,2.],[np.nan,np.nan],[3.,4.]]), present=np.array([True,False,True]), query_ids=np.array(["q0","q1","q2"]))
    with pytest.raises(ValueError,match="Incomplete"):
        load_feature_bank(bank_root, cohort="first")


def test_loader_rejects_unknown_cohort(bank_root):
    with pytest.raises(ValueError,match="Unknown cohort"):
        load_feature_bank(bank_root, cohort="typo")
