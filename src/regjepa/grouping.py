"""Conservative provenance/structure groups; pointer tokens are not recovered citations."""
from __future__ import annotations

import hashlib
import json
import re
from collections import defaultdict

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


def stable_digest(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


class Components:
    def __init__(self, names: list[str]):
        self.parent = {name: name for name in names}

    def find(self, name: str) -> str:
        self.parent.setdefault(name, name)
        root = name
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[name] != name:
            previous = self.parent[name]
            self.parent[name] = root
            name = previous
        return root

    def join(self, left: str, right: str) -> None:
        a, b = self.find(left), self.find(right)
        self.parent[max(a, b)] = min(a, b)


DOI = re.compile(r"10\.\d{4,9}/[^\s{}\"']+", re.I)
PATENT = re.compile(r"\b(?:WO|US|EP|CN|JP|KR|DE|FR|GB|IN)\d[A-Za-z0-9]*\b", re.I)
POINTER = re.compile(r"\bcomp\d+[a-z]*(?:-\d+[a-z]*)*\b", re.I)


def doi_entity(token: str) -> str:
    """Preserve internal suffix punctuation, including legacy DOI semicolons."""
    return "doi:" + token.rstrip(".,;").casefold()


def bibliography_entities(text: str) -> dict[str, list[str]]:
    entries = {}
    for match in re.finditer(r"^@\w+\{\s*([^,\s]+)\s*,(.*?)(?=^@\w+\{|\Z)", text, re.M | re.S):
        key, body = match.groups()
        entries[key.casefold()] = sorted({doi_entity(token) for token in DOI.findall(body)})
    return entries


def reference_entities(raw: str, bibliography: dict[str, list[str]]) -> dict:
    """Normalize explicit identifiers only; never infer an international patent family."""
    remaining = re.sub(r"https?://(?:dx\.)?doi\.org/", "", raw, flags=re.I)
    entities, pointers, targets, opaque = set(), [], set(), []
    for match in DOI.findall(remaining):
        entities.add(doi_entity(match))
    remaining = DOI.sub(" ", remaining)
    for match in PATENT.findall(remaining):
        # Kind codes (A/A1/A3/B2...) distinguish publications of the same numbered
        # application. Do not guess links between different national applications.
        base = re.sub(r"[A-Z]\d?$", "", match.upper())
        entities.add("patent:" + base.casefold())
    remaining = PATENT.sub(" ", remaining)
    for match in POINTER.findall(remaining):
        token = match.casefold()
        entities.add("series:" + token)
        pointers.append(token)
        pieces = token.removeprefix("comp").split("-")
        targets.update("comp" + piece for piece in pieces)
    remaining = POINTER.sub(" ", remaining)
    for token in re.findall(r"[A-Za-z][A-Za-z0-9_-]*", remaining):
        key = token.casefold()
        entities.add("citekey:" + key)
        if key in bibliography:
            entities.update(bibliography[key])
        else:
            opaque.append(key)
    return {"entities": sorted(entities), "pointer_tokens": sorted(set(pointers)),
            "pointer_targets": sorted(targets), "unresolved_citation_keys": sorted(set(opaque)),
            "missing_reference": not raw.strip(),
            "bibliographic_identity_unresolved": bool(pointers or opaque or not raw.strip())}


def molecular_similarities(records: list[dict], config: dict) -> tuple[list[str], list[list[float]]]:
    smiles = sorted({r["canonical_isomeric_smiles"] for r in records})
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=config["radius"], fpSize=config["bits"],
                                                          includeChirality=config["include_chirality"])
    fingerprints = []
    for value in smiles:
        mol = Chem.MolFromSmiles(value)
        if mol is None:
            raise ValueError("Grouping input contains invalid SMILES")
        fingerprints.append(generator.GetFingerprint(mol))
    # The lower triangle is small for this benchmark and retained only in memory.
    similarities = [list(DataStructs.BulkTanimotoSimilarity(fp, fingerprints[:i]))
                    for i, fp in enumerate(fingerprints)]
    return smiles, similarities


def build_grouping(records: list[dict], selected_ids: set[str], bibliography: str,
                   fingerprint_config: dict, *, regime: str) -> dict:
    """Build globally, then filter; excluded records can still link selected series."""
    if regime not in {"identity", "source_series", "source_series_structure"}:
        raise ValueError("Unknown grouping regime")
    names = [r["compound_id"] for r in records]
    if len(names) != len(set(names)) or not selected_ids <= set(names):
        raise ValueError("Duplicate or unknown source record ID")
    graph = Components(names)
    bib = bibliography_entities(bibliography)
    sources = {}
    edges = []
    for record in sorted(records, key=lambda r: r["compound_id"]):
        compound = record["compound_id"]
        identity = ["structure:" + record["canonical_isomeric_smiles"], "tautomer:" + record["tautomer_key"]]
        info = reference_entities(record["reference_raw"], bib)
        sources[compound] = info
        entities = identity + (info["entities"] if regime != "identity" else [])
        for entity in entities:
            graph.join(compound, entity)
            edges.append({"record": compound, "entity": entity})
        if regime != "identity":
            for target in info["pointer_targets"]:
                # Suffix aliases remain virtual nodes: comp222c is not guessed to
                # mean comp222. Explicit series-token matches still join records.
                graph.join(compound, target)
                edges.append({"record": compound, "entity": target, "reason": "explicit_pointer_target"})
    smiles, similarities = molecular_similarities(records, fingerprint_config)
    owners = {}
    for record in sorted(records, key=lambda r: r["compound_id"]):
        owners.setdefault(record["canonical_isomeric_smiles"], record["compound_id"])
    threshold = fingerprint_config["similarity_threshold"]
    if not 0 < threshold <= 1:
        raise ValueError("Invalid structural similarity threshold")
    structural_edges = []
    if regime == "source_series_structure":
        for i, values in enumerate(similarities):
            for j, similarity in enumerate(values):
                if similarity >= threshold:
                    left, right = owners[smiles[i]], owners[smiles[j]]
                    graph.join(left, right)
                    structural_edges.append({"left": left, "right": right, "tanimoto": similarity})
    all_members: dict[str, list[str]] = defaultdict(list)
    for name in sorted(names):
        all_members[graph.find(name)].append(name)
    groups, mapping = [], {}
    for members in sorted(all_members.values(), key=lambda x: x[0]):
        selected = sorted(selected_ids.intersection(members))
        if not selected:
            continue
        key = stable_digest(members)[:20]
        groups.append({"group_id": key, "selected_members": selected, "all_source_members": members})
        mapping.update({member: key for member in selected})
    cross_max = 0.0
    selected_smiles = {r["canonical_isomeric_smiles"]: r["compound_id"] for r in records if r["compound_id"] in selected_ids}
    for i, values in enumerate(similarities):
        left = selected_smiles.get(smiles[i])
        if not left:
            continue
        for j, similarity in enumerate(values):
            right = selected_smiles.get(smiles[j])
            if right and mapping[left] != mapping[right]:
                cross_max = max(cross_max, similarity)
    if regime == "source_series_structure" and cross_max >= threshold:
        raise ValueError("A structural neighbour crosses leakage components")
    return {"regime": regime, "universe_records": len(records), "selected_records": len(selected_ids),
            "groups": groups, "compound_to_group": mapping,
            "source_entities": {name: sources[name] for name in sorted(selected_ids)},
            "entity_edges": edges, "structural_edges": structural_edges,
            "fingerprint": fingerprint_config,
            "structural_policy": "connected components of all Tanimoto >= threshold edges, joined with identity and source/series edges",
            "maximum_between_group_selected_tanimoto": cross_max,
            "limitations": ["Pointer tokens preserve explicit series links but do not recover missing publications.",
                            "National patent identifiers are not independently verified international patent families.",
                            "Opaque citation keys are grouped exactly; their bibliographic identities remain unverified."]}
