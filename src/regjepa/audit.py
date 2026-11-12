"""Validate immutable EAS benchmark records and make derived, traceable tables."""
from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

from .chemistry import candidate_sites, outcome_groups, substitution_key


def reaction_evidence(note: str) -> dict:
    """Conservative text classification, retaining ambiguity rather than inventing conditions."""
    evidence = []
    patterns = {
        "bromination": r"\bNBS\b|\bBr2\b|Br2Me2|\bNBSD\b",
        "chlorination": r"\bNCS\b|SO2Cl2|\bCl2\b",
        "iodination": r"\bNIS\b|\bI2\b|\bICl\b",
        "nitration": r"\bHNO3\b",
    }
    for name, pattern in patterns.items():
        if re.search(pattern, note, flags=re.I):
            evidence.append(name)
    ambiguous_token = bool(re.search(r"\b(?:12|C12|ICI|Nal)\b", note))
    if len(evidence) > 1:
        category = "multiple_EAS_classes_evidenced"
    elif evidence:
        category = evidence[0] + "_only_evidence"
    else:
        category = "unresolved_reagent"
    repeated = bool(re.search(r"[234]x\s*(?:Br|Cl|I)|di(?:Br|Cl|I)|tri(?:Br|Cl|I)|2,3-diCl|4\+5-pos|4\+6-di", note, re.I))
    mixture = bool(re.search(r"regioisomer|two mono|\b1:1\b|\b96:4\b|ortho\s*\+\s*para|\d+%\s*(?:&|/|\+)", note, re.I))
    return {"reaction_class_evidence": evidence, "reaction_class": category,
            "ambiguous_reagent_typography": ambiguous_token,
            "multiple_yield_values": note.count("%") > 1,
            "explicit_repeated_substitution_note": repeated,
            "explicit_mixture_note": mixture,
            "yield_values_percent": [float(x) for x in re.findall(r"(\d+(?:\.\d+)?)\s*%", note)]}


def ring_summary(mol: Chem.Mol, sites: list[int]) -> dict:
    aromatic_rings = [tuple(r) for r in mol.GetRingInfo().AtomRings()
                      if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    signatures = []
    reactive_signatures = set()
    benzene = hetero = 0
    for ring in aromatic_rings:
        counts = Counter(mol.GetAtomWithIdx(i).GetSymbol() for i in ring)
        signature = f"{len(ring)}:" + "".join(f"{k}{v}" for k, v in sorted(counts.items()))
        signatures.append(signature)
        if set(ring).intersection(sites):
            reactive_signatures.add(signature)
        if len(counts) == 1 and "C" in counts:
            benzene += 1
        else:
            hetero += 1
    family = ("mixed_carbo_and_heteroaromatic" if benzene and hetero else
              "heteroaromatic_only" if hetero else "carbocyclic_aromatic_only")
    return {"molecule_family": family, "aromatic_ring_signatures": signatures,
            "observed_ring_signatures": sorted(reactive_signatures),
            "n_aromatic_rings": len(aromatic_rings), "n_heteroaromatic_rings": hetero}


def read_sources(root: Path) -> tuple[list[dict], dict]:
    directory = root / "data/raw/db_regioselectivity"
    with (directory / "compounds_doi.csv").open(encoding="utf-8", newline="") as handle:
        refs = defaultdict(list)
        for row in list(csv.reader(handle))[1:]:
            if row:
                refs[row[0].strip()].append(",".join(row[1:]).strip())
    notes = defaultdict(list)
    for line in (directory / "compounds_yield.csv").read_text(encoding="utf-8").splitlines():
        key, text = line.split(maxsplit=1)
        notes[key].append(text)
    records = []
    for line_number, line in enumerate((directory / "compounds_smiles.csv").read_text(encoding="utf-8").splitlines(), 1):
        key, smiles, labels = line.split()
        records.append({"compound_id": key, "source_line": line_number, "source_smiles": smiles,
                        "source_label_text": labels,
                        "observed_atom_indices": [int(i) for i in labels.split(",") if i],
                        "reference_raw": "; ".join(refs[key]), "yield_reagent_note_raw": "; ".join(notes[key]),
                        "source_has_non_ascii_whitespace": any(ord(c) > 127 and c.isspace() for c in line)})
    metadata = {"doi_rows": sum(map(len, refs.values())), "yield_note_rows": sum(map(len, notes.values())),
                "extra_reference_ids": sorted(set(refs) - {r["compound_id"] for r in records}),
                "extra_note_ids": sorted(set(notes) - {r["compound_id"] for r in records})}
    return records, metadata


def audit_records(root: Path) -> tuple[list[dict], dict]:
    records, metadata = read_sources(root)
    tautomerizer = rdMolStandardize.TautomerEnumerator()
    tautomerizer.SetMaxTautomers(100)
    sdf = list(Chem.SDMolSupplier(str(root / "data/raw/regiosqm20/input_mols.sdf")))
    canonical_ids = defaultdict(list)
    stereo_free_ids = defaultdict(list)
    tautomer_ids = defaultdict(list)
    all_sites = []
    for position, record in enumerate(records):
        mol = Chem.MolFromSmiles(record["source_smiles"])
        record["parse_ok"] = mol is not None
        if mol is None:
            record["validation_errors"] = ["smiles_parse_failed"]
            continue
        sites = candidate_sites(mol)
        observed = record["observed_atom_indices"]
        errors = []
        if len(set(observed)) != len(observed):
            errors.append("duplicate_observed_atom_index")
        if not observed:
            errors.append("missing_observed_atom_index")
        if not set(observed).issubset(sites):
            errors.append("observed_atom_not_aromatic_CH")
        groups = outcome_groups(mol)
        no_stereo_groups = outcome_groups(mol, keep_stereo=False)
        group_of = {atom: j for j, group in enumerate(groups) for atom in group}
        observed_groups = sorted({group_of[i] for i in observed if i in group_of})
        expanded = sorted({i for j in observed_groups for i in groups[j]})
        canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
        constitution = Chem.Mol(mol)
        Chem.RemoveStereochemistry(constitution)
        constitution_key = Chem.MolToSmiles(constitution, isomericSmiles=True)
        tautomer = tautomerizer.Canonicalize(mol)
        tautomer_key = Chem.MolToSmiles(tautomer, isomericSmiles=True)
        canonical_ids[canonical].append(record["compound_id"])
        stereo_free_ids[constitution_key].append(record["compound_id"])
        tautomer_ids[tautomer_key].append(record["compound_id"])
        sdf_mol = sdf[position] if position < len(sdf) else None
        sdf_matches = sdf_mol is not None and Chem.MolToSmiles(sdf_mol, isomericSmiles=True) == canonical
        sdf_order_matches = sdf_matches and sdf_mol.GetNumAtoms() == mol.GetNumAtoms() and all(
            (sdf_mol.GetAtomWithIdx(i).GetAtomicNum(), sdf_mol.GetAtomWithIdx(i).GetIsAromatic()) ==
            (mol.GetAtomWithIdx(i).GetAtomicNum(), mol.GetAtomWithIdx(i).GetIsAromatic())
            for i in range(mol.GetNumAtoms())) and sorted(
                (min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), str(b.GetBondType())) for b in sdf_mol.GetBonds()) == sorted(
                (min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), str(b.GetBondType())) for b in mol.GetBonds())
        classes = reaction_evidence(record["yield_reagent_note_raw"])
        review = []
        if not record["reference_raw"]:
            review.append("missing_reference")
        if not record["yield_reagent_note_raw"]:
            review.append("missing_note")
        if len(observed_groups) > 1:
            review.append("multiple_distinct_annotated_outcomes")
        for field in ["ambiguous_reagent_typography", "multiple_yield_values", "explicit_repeated_substitution_note", "explicit_mixture_note"]:
            if classes[field]:
                review.append(field)
        if classes["reaction_class"] != "bromination_only_evidence":
            review.append("not_bromination_only_evidence")
        if not sdf_matches:
            review.append("sdf_row_structure_differs")
        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=True, includeIsotopes=True))
        rank_groups = defaultdict(list)
        for i in sites:
            rank_groups[ranks[i]].append(i)
        rank_agrees = sorted(sorted(g) for g in rank_groups.values()) == sorted(sorted(g) for g in groups)
        if not rank_agrees:
            review.append("symmetry_rank_product_disagreement")
        record.update(classes)
        record.update(ring_summary(mol, observed))
        mapped = Chem.Mol(mol)
        for atom in mapped.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)
        record.update({"validation_errors": errors, "review_flags": review,
                       "atom_mapped_smiles": Chem.MolToSmiles(mapped, isomericSmiles=True),
                       "canonical_isomeric_smiles": canonical, "constitution_key": constitution_key,
                       "tautomer_key": tautomer_key, "murcko_scaffold": MurckoScaffold.MurckoScaffoldSmiles(mol=mol),
                       "formula": rdMolDescriptors.CalcMolFormula(mol), "molecular_weight": round(Descriptors.MolWt(mol), 3),
                       "n_heavy_atoms": mol.GetNumHeavyAtoms(), "formal_charge": Chem.GetFormalCharge(mol),
                       "n_fragments": len(Chem.GetMolFrags(mol)),
                       "n_assigned_chiral_centers": sum(label != "?" for _, label in Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
                       "n_unassigned_chiral_centers": sum(label == "?" for _, label in Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
                       "candidate_atom_indices": sites, "candidate_outcome_groups": groups,
                       "constitution_only_outcome_groups": no_stereo_groups,
                       "observed_outcome_group_indices": observed_groups,
                       "expanded_equivalent_observed_atoms": expanded,
                       "n_candidate_atoms": len(sites), "n_candidate_outcomes": len(groups),
                       "n_observed_atoms": len(observed), "n_observed_outcomes": len(observed_groups),
                       "has_equivalent_candidate_atoms": any(len(g) > 1 for g in groups),
                       "has_equivalent_unannotated_positive": bool(set(expanded) - set(observed)),
                       "symmetry_rank_product_agree": rank_agrees,
                       "stereo_changes_outcome_partition": groups != no_stereo_groups,
                       "sdf_same_structure_by_row": sdf_matches, "sdf_same_atom_order_by_row": sdf_order_matches,
                       "sdf_has_3d_coordinates": bool(sdf_mol and sdf_mol.GetNumConformers() and sdf_mol.GetConformer().Is3D()),
                       "eligible_broad_single_outcome": not errors and len(groups) >= 2 and len(observed_groups) == 1,
                       "eligible_conservative_bromination": not errors and len(groups) >= 2 and len(observed_groups) == 1 and not review,
                       "uniform_outcome_top1_chance": round(len(observed_groups) / len(groups), 6) if groups else None})
        for j, group in enumerate(groups):
            for atom_index in group:
                all_sites.append({"compound_id": record["compound_id"], "atom_index": atom_index,
                                  "outcome_group_index": j, "group_size": len(group),
                                  "source_observed": atom_index in observed,
                                  "symmetry_expanded_observed": atom_index in expanded,
                                  "hypothetical_Br_marker_product": substitution_key(mol, atom_index)})
    for record in records:
        if not record["parse_ok"]:
            continue
        record["same_structure_ids"] = canonical_ids[record["canonical_isomeric_smiles"]]
        record["same_constitution_ids"] = stereo_free_ids[record["constitution_key"]]
        record["same_canonical_tautomer_ids"] = tautomer_ids[record["tautomer_key"]]
        record["suggested_identity_group"] = hashlib.sha256(record["tautomer_key"].encode()).hexdigest()[:16]
    summary = {"schema_version": 1, "rdkit_version": rdBase.rdkitVersion, "records": len(records),
               "unique_ids": len({r["compound_id"] for r in records}), "metadata": metadata,
               "unique_isomeric_structures": len(canonical_ids),
               "exact_duplicate_groups": [v for v in canonical_ids.values() if len(v) > 1],
               "tautomer_identity_groups_with_multiple_records": [v for v in tautomer_ids.values() if len(v) > 1],
               "sdf_records": len(sdf), "candidate_atom_records": len(all_sites)}
    for key in ["parse_ok", "has_equivalent_candidate_atoms", "has_equivalent_unannotated_positive",
                "stereo_changes_outcome_partition", "sdf_same_structure_by_row", "sdf_same_atom_order_by_row",
                "sdf_has_3d_coordinates", "eligible_broad_single_outcome", "eligible_conservative_bromination"]:
        summary[key] = sum(bool(r.get(key)) for r in records)
    summary["validation_error_records"] = sum(bool(r["validation_errors"]) for r in records)
    for key in ["reaction_class", "molecule_family", "n_candidate_atoms", "n_candidate_outcomes", "n_observed_atoms", "n_observed_outcomes"]:
        summary[key + "_counts"] = dict(sorted(Counter(r.get(key) for r in records).items(), key=lambda x: str(x[0])))
    summary["review_flag_counts"] = dict(Counter(flag for r in records for flag in r.get("review_flags", [])))
    summary["observed_ring_signature_counts"] = dict(Counter(sig for r in records for sig in r.get("observed_ring_signatures", [])))
    summary["reference_type_counts"] = dict(Counter("missing" if not r["reference_raw"] else "patent" if re.match(r"(?:WO|US|EP|JP|CN|GB)", r["reference_raw"]) else "doi" if r["reference_raw"].startswith("10.") else "bibliographic_or_other" for r in records))
    summary["unique_nonempty_reference_strings"] = len({r["reference_raw"] for r in records if r["reference_raw"]})
    for key in ["n_heavy_atoms", "molecular_weight"]:
        vals = sorted(r[key] for r in records if key in r)
        summary[key + "_range"] = {"min": vals[0], "median": vals[len(vals)//2], "max": vals[-1]}
    for label, selection in [("all", records), ("broad_single_outcome", [r for r in records if r.get("eligible_broad_single_outcome")]),
                             ("conservative_bromination", [r for r in records if r.get("eligible_conservative_bromination")])]:
        summary[label + "_uniform_outcome_top1_chance"] = sum(r["uniform_outcome_top1_chance"] for r in selection) / len(selection) if selection else None
        summary[label + "_cohort"] = {
            "records": len(selection),
            "unique_structures": len({r["canonical_isomeric_smiles"] for r in selection}),
            "has_equivalent_candidate_atoms": sum(r["has_equivalent_candidate_atoms"] for r in selection),
            "has_equivalent_unannotated_positive": sum(r["has_equivalent_unannotated_positive"] for r in selection),
            "molecule_family_counts": dict(Counter(r["molecule_family"] for r in selection)),
            "reaction_class_counts": dict(Counter(r["reaction_class"] for r in selection))}
    summary["multiple_atom_labels_same_outcome"] = sum(r["n_observed_atoms"] > 1 and r["n_observed_outcomes"] == 1 for r in records)
    summary["no_competing_outcome_ids"] = [r["compound_id"] for r in records if r["n_candidate_outcomes"] < 2]
    summary["total_distinct_candidate_outcomes"] = sum(r["n_candidate_outcomes"] for r in records)
    summary["assigned_stereocenter_records"] = sum(r["n_assigned_chiral_centers"] > 0 for r in records)
    summary["unassigned_stereocenter_records"] = sum(r["n_unassigned_chiral_centers"] > 0 for r in records)
    summary["non_ascii_whitespace_records"] = sum(r["source_has_non_ascii_whitespace"] for r in records)
    summary["sites"] = all_sites
    return records, summary


def write_csv(path: Path, records: list[dict]) -> None:
    if not records:
        return
    fields = list(dict.fromkeys(k for r in records for k in r))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for record in records:
            writer.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v for k, v in record.items()})


def run(root: Path) -> dict:
    records, summary = audit_records(root)
    processed = root / "data/processed"
    processed.mkdir(parents=True, exist_ok=True)
    write_csv(processed / "molecules.csv", records)
    sites = summary.pop("sites")
    write_csv(processed / "candidate_sites.csv", sites)
    with (processed / "molecules.jsonl").open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=False) + "\n")
    (processed / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    cohorts = {"schema_version": 1, "split_status": "No train/validation/test assignments frozen",
               "broad_single_outcome": [r["compound_id"] for r in records if r.get("eligible_broad_single_outcome")],
               "conservative_bromination": [r["compound_id"] for r in records if r.get("eligible_conservative_bromination")],
               "no_competing_outcomes": [r["compound_id"] for r in records if r.get("n_candidate_outcomes", 0) < 2],
               "multiple_distinct_observed_outcomes": [r["compound_id"] for r in records if r.get("n_observed_outcomes", 0) > 1]}
    (processed / "cohorts.json").write_text(json.dumps(cohorts, indent=2) + "\n", encoding="utf-8")
    return summary

