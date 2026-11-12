"""Chemistry operations that preserve the source atom order and specified stereochemistry."""
from __future__ import annotations

from collections import defaultdict
from rdkit import Chem


def candidate_sites(mol: Chem.Mol) -> list[int]:
    """Aromatic carbon sites bearing H; this is an explicit first-EAS task definition."""
    return [a.GetIdx() for a in mol.GetAtoms()
            if a.GetAtomicNum() == 6 and a.GetIsAromatic() and a.GetTotalNumHs() > 0]


def substitution_key(mol: Chem.Mol, site: int, *, keep_stereo: bool = True) -> str:
    """Canonical identity after replacing one aromatic C-H by a Br marker.

    The marker is a symmetry probe, NOT an asserted observed product or a reaction
    prediction. Single substitution preserves the input tautomer and isotopes.
    """
    if site not in candidate_sites(mol):
        raise ValueError(f"Atom {site} is not an eligible aromatic C-H")
    edited = Chem.RWMol(mol)
    atom = edited.GetAtomWithIdx(site)
    if atom.GetNumExplicitHs():
        atom.SetNumExplicitHs(atom.GetNumExplicitHs() - 1)
    marker = edited.AddAtom(Chem.Atom(35))
    edited.AddBond(site, marker, Chem.BondType.SINGLE)
    product = edited.GetMol()
    for a in product.GetAtoms():
        a.SetAtomMapNum(0)
    Chem.SanitizeMol(product)
    if not keep_stereo:
        Chem.RemoveStereochemistry(product)
    return Chem.MolToSmiles(product, canonical=True, isomericSmiles=True)


def outcome_groups(mol: Chem.Mol, *, keep_stereo: bool = True) -> list[list[int]]:
    """Group sites yielding identical marker products, not geometric near-symmetry."""
    groups: dict[str, list[int]] = defaultdict(list)
    for site in candidate_sites(mol):
        groups[substitution_key(mol, site, keep_stereo=keep_stereo)].append(site)
    return sorted((sorted(g) for g in groups.values()), key=lambda g: g[0])


def mean_group_scores(atom_scores: dict[int, float], groups: list[list[int]]) -> list[float]:
    """Average site logits so an outcome receives no automatic multiplicity bonus."""
    return [sum(atom_scores[i] for i in group) / len(group) for group in groups]


def top1_hit(atom_scores: dict[int, float], groups: list[list[int]], observed: list[int]) -> bool:
    """Outcome-level top-1 membership; does not claim a product-ratio prediction."""
    scores = mean_group_scores(atom_scores, groups)
    best = max(range(len(groups)), key=lambda i: scores[i])
    return bool(set(groups[best]).intersection(observed))

