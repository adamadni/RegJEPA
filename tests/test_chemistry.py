import pytest
from rdkit import Chem
from regjepa.chemistry import candidate_sites, outcome_groups, substitution_key, mean_group_scores, top1_hit
from regjepa.audit import reaction_evidence


def test_benzene_has_one_distinct_outcome():
    groups = outcome_groups(Chem.MolFromSmiles("c1ccccc1"))
    assert groups == [[0, 1, 2, 3, 4, 5]]


def test_anisole_equivalent_ortho_and_meta_sites():
    mol = Chem.MolFromSmiles("COc1ccccc1")
    assert sorted(map(len, outcome_groups(mol))) == [1, 2, 2]


def test_explicit_aromatic_hydrogen_handled():
    mol = Chem.MolFromSmiles("[cH]1ccccc1")
    assert len(outcome_groups(mol)) == 1
    assert Chem.MolFromSmiles(substitution_key(mol, 0)) is not None


def test_not_all_positions_in_heteroaromatic_are_equivalent():
    mol = Chem.MolFromSmiles("n1ccccc1")
    assert sorted(map(len, outcome_groups(mol))) == [1, 2, 2]
    assert 0 not in candidate_sites(mol)


def test_isotopic_label_breaks_benzene_symmetry():
    mol = Chem.MolFromSmiles("[13cH]1ccccc1")
    assert len(outcome_groups(mol)) == 4


def test_atom_map_numbers_do_not_break_equivalence():
    mol = Chem.MolFromSmiles("[cH:17]1[cH:3]cccc1")
    assert len(outcome_groups(mol)) == 1


def test_renumbering_preserves_product_identity_partition():
    mol = Chem.MolFromSmiles("COc1cccnc1")
    order = list(reversed(range(mol.GetNumAtoms())))
    renumbered = Chem.RenumberAtoms(mol, order)
    original_keys = {i: substitution_key(mol, i) for i in candidate_sites(mol)}
    for new, old in enumerate(order):
        if old in original_keys:
            assert substitution_key(renumbered, new) == original_keys[old]


def test_invalid_candidate_rejected():
    with pytest.raises(ValueError):
        substitution_key(Chem.MolFromSmiles("COc1ccccc1"), 0)


def test_product_keys_preserve_specified_stereochemistry():
    right = Chem.MolFromSmiles("N[C@H](C)c1ccccc1")
    left = Chem.MolFromSmiles("N[C@@H](C)c1ccccc1")
    right_keys = {substitution_key(right, i) for i in candidate_sites(right)}
    left_keys = {substitution_key(left, i) for i in candidate_sites(left)}
    assert right_keys.isdisjoint(left_keys)
    assert {substitution_key(right, i, keep_stereo=False) for i in candidate_sites(right)} == {
        substitution_key(left, i, keep_stereo=False) for i in candidate_sites(left)}


def test_marker_generation_does_not_change_input():
    mol = Chem.MolFromSmiles("COc1ccccc1")
    before = Chem.MolToSmiles(mol)
    outcome_groups(mol)
    assert Chem.MolToSmiles(mol) == before


def test_group_scoring_has_no_symmetry_size_bonus():
    groups = [[0, 1], [2]]
    assert mean_group_scores({0: 0.4, 1: 0.4, 2: 0.5}, groups) == [0.4, 0.5]
    assert top1_hit({0: 0.4, 1: 0.4, 2: 0.5}, groups, [2])
    assert top1_hit({0: 0.7, 1: 0.7, 2: 0.5}, groups, [1])


@pytest.mark.parametrize("text,expected", [
    ("100% (NBS)", "bromination_only_evidence"),
    ("97% (NIS, TFA)", "iodination_only_evidence"),
    ("74% (SO2Cl2)", "chlorination_only_evidence"),
    ("98% (HNO3)", "nitration_only_evidence"),
    ("80% (NBS), 63% (NIS)", "multiple_EAS_classes_evidenced"),
    ("12, 94%", "unresolved_reagent")])
def test_reaction_class_evidence(text, expected):
    assert reaction_evidence(text)["reaction_class"] == expected


def test_no_invention_of_ratios_or_reagents():
    result = reaction_evidence("12, 94%")
    assert result["ambiguous_reagent_typography"]
    assert result["reaction_class_evidence"] == []
    assert reaction_evidence("2xBr 96% (NBS)")["explicit_repeated_substitution_note"]

