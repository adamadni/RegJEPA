import copy

from regjepa.grouping import bibliography_entities, build_grouping, reference_entities

FP = {"radius": 2, "bits": 2048, "include_chirality": True, "similarity_threshold": .6}


def records():
    smiles = ["c1ccccc1", "COc1ccccc1", "c1ccncc1", "c1ccoc1"]
    return [{"compound_id": f"comp{i+1}", "canonical_isomeric_smiles": s, "tautomer_key": s,
             "reference_raw": f"10.1234/source{i}"} for i, s in enumerate(smiles)]


def test_explicit_identifier_normalization_preserves_meaning_and_joins_kind_codes():
    bib = bibliography_entities('@article{Example,\n doi={10.1234/AbC},\n}\n')
    info = reference_entities("https://doi.org/10.1234/ABC WO2011055115A1 Example", bib)
    assert "doi:10.1234/abc" in info["entities"]
    assert "patent:wo2011055115" in info["entities"]
    assert not info["unresolved_citation_keys"]
    assert reference_entities("WO2011055115A3", bib)["entities"] == ["patent:wo2011055115"]


def test_mixed_pointer_and_patent_links_through_out_of_cohort_record():
    r = records()
    r[0]["reference_raw"] = "comp1-2"
    r[1]["reference_raw"] = "WO2007039297A1 comp1-2"
    r[2]["reference_raw"] = "WO2007039297A3"
    grouped = build_grouping(r, {"comp1", "comp3", "comp4"}, "", FP, regime="source_series")
    mapping = grouped["compound_to_group"]
    assert mapping["comp1"] == mapping["comp3"]
    assert mapping["comp1"] != mapping["comp4"]
    assert grouped["source_entities"]["comp1"]["bibliographic_identity_unresolved"]


def test_suffix_aliases_are_not_guessed_to_be_unsuffixed_record_ids():
    r = records()
    r[0]["reference_raw"] = "comp2a"
    grouped = build_grouping(r, {x["compound_id"] for x in r}, "", FP, regime="source_series")
    assert grouped["compound_to_group"]["comp1"] != grouped["compound_to_group"]["comp2"]


def test_structural_linkage_is_transitive_and_label_blind(monkeypatch):
    r = records()[:3]
    smiles = [row["canonical_isomeric_smiles"] for row in r]
    monkeypatch.setattr("regjepa.grouping.molecular_similarities", lambda records, config: (smiles, [[], [.7], [.2, .7]]))
    grouped = build_grouping(r, {"comp1", "comp3"}, "", FP, regime="source_series_structure")
    assert grouped["compound_to_group"]["comp1"] == grouped["compound_to_group"]["comp3"]
    changed = copy.deepcopy(r)
    for record in changed:
        record["observed_atom_indices"] = [999]
        record["yield_reagent_note_raw"] = "unrelated changed label text"
    assert build_grouping(changed, {"comp1", "comp3"}, "", FP, regime="source_series_structure") == grouped


def test_source_record_order_does_not_change_group_assignments():
    r = records()
    first = build_grouping(r, {x["compound_id"] for x in r}, "", FP, regime="source_series_structure")
    second = build_grouping(list(reversed(r)), {x["compound_id"] for x in r}, "", FP, regime="source_series_structure")
    assert first == second
