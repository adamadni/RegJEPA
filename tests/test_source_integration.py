from pathlib import Path

import pytest
from rdkit import Chem

from regjepa.audit import audit_records
from regjepa.chemistry import substitution_key


@pytest.fixture(scope="module")
def audited():
    records, summary = audit_records(Path(__file__).resolve().parents[1])
    return {r['compound_id']: r for r in records}, summary


def test_source_atom_convention_and_full_join(audited):
    rows, summary = audited
    assert len(rows) == 535
    assert summary['validation_error_records'] == 0
    assert summary['metadata']['extra_note_ids'] == []
    assert summary['metadata']['extra_reference_ids'] == []
    assert rows['comp14']['observed_atom_indices'] == [0]


def test_missing_symmetry_counterpart_is_preserved_and_expanded(audited):
    rows, _ = audited
    assert rows['comp18']['observed_atom_indices'] == [2]
    assert rows['comp18']['expanded_equivalent_observed_atoms'] == [2, 8]


def test_renumbered_duplicate_has_same_observed_outcome(audited):
    rows, _ = audited
    a, b = rows['comp13'], rows['comp598']
    assert a['canonical_isomeric_smiles'] == b['canonical_isomeric_smiles']
    assert substitution_key(Chem.MolFromSmiles(a['source_smiles']), a['observed_atom_indices'][0]) == substitution_key(
        Chem.MolFromSmiles(b['source_smiles']), b['observed_atom_indices'][0])
    assert a['suggested_identity_group'] == b['suggested_identity_group']


def test_missing_reference_is_metadata_review_not_bad_chemistry(audited):
    rows, _ = audited
    assert rows['comp155']['validation_errors'] == []
    assert 'missing_reference' in rows['comp155']['review_flags']
    assert not rows['comp155']['eligible_conservative_bromination']


def test_sdf_stereo_discrepancy_is_not_silently_used(audited):
    rows, _ = audited
    assert not rows['comp485']['sdf_same_structure_by_row']
    assert 'sdf_row_structure_differs' in rows['comp485']['review_flags']
    assert not rows['comp485']['eligible_broad_single_outcome']


def test_canonical_mapped_smiles_retains_source_index_mapping(audited):
    rows, _ = audited
    for row in rows.values():
        mol = Chem.MolFromSmiles(row['atom_mapped_smiles'])
        atom_by_map = {a.GetAtomMapNum()-1: a for a in mol.GetAtoms()}
        assert set(atom_by_map) == set(range(mol.GetNumAtoms()))
        for idx in row['observed_atom_indices']:
            atom = atom_by_map[idx]
            assert atom.GetAtomicNum() == 6 and atom.GetIsAromatic() and atom.GetTotalNumHs() > 0
