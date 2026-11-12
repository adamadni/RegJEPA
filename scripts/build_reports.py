"""Build an offline molecule browser, molecular examples, and a quantitative audit."""
from __future__ import annotations

import csv
import html
import json
import math
from collections import Counter
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "reports"
GREEN = (0.55, 0.89, 0.65)
BLUE = (0.55, 0.75, 1.0)


def depiction(record: dict) -> str:
    mol = Chem.MolFromSmiles(record["source_smiles"])
    colors = {i: BLUE for i in record["expanded_equivalent_observed_atoms"]}
    colors.update({i: GREEN for i in record["observed_atom_indices"]})
    drawer = rdMolDraw2D.MolDraw2DSVG(360, 230)
    drawer.drawOptions().addAtomIndices = True
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, highlightAtoms=list(colors), highlightAtomColors=colors)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg[svg.index("<svg"):]


def wilson(n: int, p: float = 0.8) -> tuple[float, float]:
    z = 1.959963984540054
    denom = 1 + z*z/n
    center = (p+z*z/(2*n))/denom
    half = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/denom
    return center-half, center+half


def main() -> None:
    REPORTS.mkdir(exist_ok=True)
    records = [json.loads(line) for line in (ROOT / "data/processed/molecules.jsonl").read_text(encoding="utf-8").splitlines()]
    summary = json.loads((ROOT / "data/processed/summary.json").read_text(encoding="utf-8"))
    cards = []
    for r in records:
        escaped = lambda value: html.escape(str(value))
        cards.append(f'''<article data-family="{r['molecule_family']}" data-class="{r['reaction_class']}"
 data-sym="{str(r['has_equivalent_candidate_atoms']).lower()}" data-primary="{str(r['eligible_conservative_bromination']).lower()}"
 data-search="{escaped(' '.join([r['compound_id'], r['source_smiles'], r['reference_raw'], r['yield_reagent_note_raw']]))}">
<h2>{r['compound_id']} <small>{r['formula']}</small></h2>{depiction(r)}
<p><b>{r['reaction_class'].replace('_',' ')}</b><br>{r['molecule_family'].replace('_',' ')}</p>
<p>Source atoms: <b>{r['observed_atom_indices']}</b>; equivalent positives: <b>{r['expanded_equivalent_observed_atoms']}</b><br>
{r['n_candidate_atoms']} candidate atoms → {r['n_candidate_outcomes']} distinct marker outcomes</p>
<p>Candidate atom groups: <code>{r['candidate_outcome_groups']}</code></p>
<p class="note">{escaped(r['yield_reagent_note_raw'])}</p>
<p>Reference: {escaped(r['reference_raw']) or '<b>missing upstream</b>'}</p>
<p class="flags">{'Conservative bromination cohort' if r['eligible_conservative_bromination'] else escaped(', '.join(r['review_flags']) or 'No competing site')}</p>
<details><summary>Structure and source identity</summary><code>{escaped(r['source_smiles'])}</code>
<p>Source line {r['source_line']}; duplicate group {r['same_structure_ids']}; canonical tautomer group {r['same_canonical_tautomer_ids']}.</p></details></article>''')
    options = lambda values: ''.join(f'<option value="{html.escape(x)}">{html.escape(x.replace("_", " "))}</option>' for x in sorted(values))
    browser = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>RegJEPA · EAS molecule audit</title><style>
body{font:15px/1.5 system-ui,sans-serif;background:#f3f5f7;color:#182531;margin:0}header{padding:24px 5%;background:#183d4f;color:white}
header h1{margin:0}main{padding:20px 3%}.controls{position:sticky;top:0;background:white;padding:14px;display:flex;gap:10px;flex-wrap:wrap;border-bottom:1px solid #ccc;z-index:1}
input,select{padding:8px;max-width:100%;border:1px solid #aab8c2;border-radius:5px}#cards{display:grid;grid-template-columns:repeat(auto-fill,minmax(365px,1fr));gap:16px}
article{background:white;border:1px solid #d8e0e5;border-radius:10px;padding:16px;overflow-wrap:anywhere}article[hidden]{display:none}h2{font-size:19px;margin:0}small{font-weight:400;color:#456}
svg{max-width:100%;height:auto}p{margin:9px 0}.note{padding:8px;background:#eef3f7}.flags{color:#725414;font-size:13px}code{font-size:12px}a{color:#176185}#count{margin:16px 0}.legend{font-size:14px}
</style><header><h1>RegJEPA · 535-record EAS audit</h1><p>Graph and annotation audit · original source atom numbering is zero-based</p>
<p class="legend">Green = source-annotated site. Blue = equivalent site missing from the source labels. Equivalence means the same hypothetical single-Br marker product on the specified graph.</p>
<p>Reaction classes below are text-derived evidence, not verified condition-resolved experiments. The source contains multiple halogenations and pooled experiments.</p></header>
<div class="controls"><input id="search" placeholder="Search ID, SMILES, reference, reagent"><select id="family"><option value="">All molecule families</option>FAMILYOPTIONS</select>
<select id="class"><option value="">All reaction evidence</option>CLASSOPTIONS</select><select id="sym"><option value="">All symmetry cases</option><option value="true">Equivalent candidate sites</option><option value="false">No equivalent candidate sites</option></select>
<select id="primary"><option value="">All cohorts</option><option value="true">Conservative bromination cohort</option><option value="false">Other / review</option></select></div>
<main><p id="count"></p><section id="cards">CARDS</section></main><script>
const fields=['search','family','class','sym','primary'];function filter(){let n=0;for(const card of document.querySelectorAll('article')){const query=document.getElementById('search').value.toLowerCase();let show=card.dataset.search.toLowerCase().includes(query);for(const key of fields.slice(1)){const v=document.getElementById(key).value;if(v&&card.dataset[key]!==v)show=false;}card.hidden=!show;if(show)n++;}document.getElementById('count').textContent=`${n} of 535 records shown`;}
for(const key of fields)document.getElementById(key).addEventListener('input',filter);filter();</script></html>'''
    browser = browser.replace('FAMILYOPTIONS', options({r['molecule_family'] for r in records})).replace('CLASSOPTIONS', options({r['reaction_class'] for r in records})).replace('CARDS', '\n'.join(cards))
    (REPORTS / "molecule_browser.html").write_text(browser, encoding="utf-8")
    chosen = ['comp18','comp21','comp488','comp492','comp274','comp158','comp13','comp598','comp485']
    chosen_records = [next(r for r in records if r['compound_id'] == key) for key in chosen]
    options = rdMolDraw2D.MolDrawOptions()
    options.addAtomIndices = True
    options.legendFontSize = 15
    image = Draw.MolsToGridImage([Chem.MolFromSmiles(r['source_smiles']) for r in chosen_records], molsPerRow=3, subImgSize=(420, 330),
                                legends=[f"{r['compound_id']}: labels {r['observed_atom_indices']}\n{r['n_candidate_atoms']} sites / {r['n_candidate_outcomes']} outcomes" for r in chosen_records],
                                highlightAtomLists=[r['expanded_equivalent_observed_atoms'] for r in chosen_records],
                                highlightAtomColors=[{i: GREEN if i in r['observed_atom_indices'] else BLUE for i in r['expanded_equivalent_observed_atoms']} for r in chosen_records],
                                drawOptions=options)
    image.save(str(REPORTS / "symmetry_examples.png"))
    table = lambda counts: '\n'.join(f'| {key} | {value} |' for key,value in counts.items())
    rows = []
    for label, predicate in [('All source records', lambda r: True), ('Broad single-outcome', lambda r: r['eligible_broad_single_outcome']), ('Conservative bromination', lambda r: r['eligible_conservative_bromination'])]:
        selected = [r for r in records if predicate(r)]
        rows.append(f"| {label} | {len(selected)} | {len({r['canonical_isomeric_smiles'] for r in selected})} | {sum(r['has_equivalent_candidate_atoms'] for r in selected)} |")
    lo, hi = wilson(summary['eligible_conservative_bromination'])
    text = f'''# RegJEPA initial EAS data audit

Audit scope: pinned public sources; RDKit {summary['rdkit_version']}. This is a structural and annotation audit, not independent re-extraction of every literature experiment.

## Decision

Suitable for a frozen-feature feasibility experiment. Insufficient alone for training a JEPA from scratch, establishing superiority across reaction classes, or modeling quantitative product ratios. Use the conservative bromination cohort initially and keep the other records visible as review and sensitivity cohorts.

## Acquired data and integrity

- {summary['records']} source rows, {summary['unique_ids']} IDs, {summary['unique_isomeric_structures']} distinct isomeric structures.
- All {summary['parse_ok']} substrate SMILES parse and every observed index is an aromatic carbon bearing H under the zero-based source convention.
- {summary['validation_error_records']} structural/label validation failures.
- {summary['review_flag_counts'].get('missing_reference',0)} rows have an empty reference field upstream. They are retained and flagged, never forward-filled.
- The SDF contains {summary['sdf_records']} records: {summary['sdf_same_structure_by_row']} match the source graph/stereochemistry and atom order. comp485 has additional alkene stereochemistry in the SDF. The SMILES source remains authoritative.
- All SDF coordinates are 2D; these are not suitable 3D Atom-JEPA inputs without conformer generation.
- Raw files, upstream licenses, commit IDs, Git blob hashes and SHA256 digests are preserved in `data/raw/` and `data/sources.lock.json`.

## Reaction evidence from the original free-text notes

Although the RegioSQM20 paper describes a 535-molecule bromination benchmark, the original annotation notes contain other electrophiles and multiple experiments. The following counts are heuristic text classifications; reaction conditions and product-to-experiment assignments have not been verified from the cited sources.

| Note-derived category | Molecule records |
|---|---:|
{table(summary['reaction_class_counts'])}

The {summary['reaction_class_counts']['multiple_EAS_classes_evidenced']} mixed-class records belong to one category each; the table does not count them repeatedly in the single-class categories. Ambiguous strings such as `12`, `C12`, `ICI` and `Nal` are flagged rather than silently corrected. Recognized reagent evidence can coexist with ambiguous typography. This source has no enolization, 1,2/1,4 addition or cycloaddition benchmark.

## Molecular coverage

| Graph-derived family | Records |
|---|---:|
{table(summary['molecule_family_counts'])}

Family refers to the aromatic rings, not to whether the rest of the molecule contains heteroatoms. The mixed class has both carbocyclic aromatic and heteroaromatic rings. Ring signatures below are RDKit ring-basis descriptors rather than named chemical series; fused systems can contribute more than one signature.

| Ring containing an annotated site (size:element counts) | Records containing that signature |
|---|---:|
{table(summary['observed_ring_signature_counts'])}

Heavy atoms: {summary['n_heavy_atoms_range']}; molecular weight: {summary['molecular_weight_range']} Da. There are {summary['assigned_stereocenter_records']} records with assigned tetrahedral stereocenters and {summary['unassigned_stereocenter_records']} with unassigned potential centers. This does not count all forms of stereochemistry or establish stereoselectivity labels.

## Symmetry findings

- {summary['candidate_atom_records']} eligible aromatic C-H atoms reduce to {summary['total_distinct_candidate_outcomes']} distinct single-substitution outcomes across records.
- {summary['has_equivalent_candidate_atoms']} molecules contain at least one symmetry-equivalent pair of candidate atoms.
- {summary['has_equivalent_unannotated_positive']} molecules have at least one equivalent counterpart of an observed site omitted from the source labels.
- {sum(r['n_observed_atoms']>1 for r in records)} records list multiple atom labels, but {summary['multiple_atom_labels_same_outcome']} of those collapse to one distinct marker outcome. The remaining {sum(r['n_observed_outcomes']>1 for r in records)} records have multiple distinct annotated outcomes and need experiment/mixture/multiple-substitution review.
- {len(summary['no_competing_outcome_ids'])} molecules have no nontrivial site choice: {', '.join(summary['no_competing_outcome_ids'])}. Exclude them from the primary ranking metric and report them separately.
- Graph canonical symmetry ranks and explicit marker-product partitions agree for every record. Removing specified stereochemistry changes no within-molecule partition in this dataset; retain stereo-aware handling for later datasets.
- comp13 and comp598 are the same substrate and the same marked substitution outcome under different halogenation annotations. Keep them in one split group; do not treat them as independent chemistry.

### Recommended policy

Define equivalence by identical hypothetical single-Br marker products on the input graph, preserving its tautomer, isotope labels and specified stereochemistry. No 3D distance tolerance is used. Green source labels and blue equivalent counterparts are visible in `molecule_browser.html` and `symmetry_examples.png`.

Retain every raw atom label. Add equivalence groups and expanded positive sets in derived files. Average site features or logits within each group, then normalize/rank across distinct groups. This avoids giving a group with two symmetric atoms an automatic counting advantage. A degeneracy/rate correction is a different physical modeling question and requires quantitative kinetic data. Keep the molecule as one evaluation unit.

For a genuinely multisite observation, do not assume equal probabilities or infer that all labels describe a mixture in one experiment. Those rows remain a separate review cohort. The product marker is a symmetry test, not a claim that every row is experimentally brominated.

## Proposed starting cohorts (no train/test split frozen yet)

| Cohort | Records | Unique structures | Records with equivalent candidates |
|---|---:|---:|---:|
{chr(10).join(rows)}

Broad single-outcome requires valid structure/site labels, at least two distinct candidate outcomes, and exactly one annotated outcome after symmetry reduction; it can include mixed reaction evidence, missing references and repeated-experiment notes, so it is a sensitivity cohort. Conservative bromination additionally requires only bromination evidence, a nonempty reference, no multiple-yield, explicit mixture, repeated-substitution, ambiguous-reagent, SDF-mismatch or symmetry-disagreement flag. These are screening rules, not a chemical revalidation of all experiments. Further manual review can reduce or enlarge the cohort.

Uniform random choice among distinct candidate outcomes yields a mean top-1 hit rate of {summary['conservative_bromination_uniform_outcome_top1_chance']:.1%} on the conservative cohort (computed from candidate counts; not a trained baseline).

## Is this large/high-quality enough?

The data are traceable, literature-derived, and reusable with an upstream MIT license, and the source atom indices pass structural checks. That is enough to start a small frozen-encoder comparison with local fingerprints and QM descriptors. It is not a uniform-condition, independently replicated, mixture-resolved set of 535 experiments. Metadata curation is needed before interpreting learned condition effects or product ratios.

The conservative cohort has {summary['eligible_conservative_bromination']} molecules; its candidate atoms are not independent additional training examples. At a hypothetical 80% accuracy over {summary['eligible_conservative_bromination']} independent predictions, an illustrative 95% Wilson interval is {lo:.1%}–{hi:.1%}. A single small held-out test is less precise, and structural/source correlation reduces effective sample size further. Use grouped outer cross-validation, inner-only model selection, low-label learning curves and paired molecule/group resampling. A gain of a few cases is exploratory, not a broad result.

Keep identical structures, tautomer variants and records from the same experiment/source together. Canonical tautomer IDs are supplied as grouping aids, not as substituted feature inputs. Resolve multi-reference strings and group publication/patent families before freezing splits. The large benzene core makes a simplistic scaffold-only split inadequate; evaluate fingerprint-cluster splits and report family coverage. Audit pretrained corpus overlap and label exposure before benchmarking published checkpoints.

## Files and reproduction

- `data/processed/molecules.csv`: one row per source molecule, raw annotations, descriptor summary, groups, flags and cohort inclusion.
- `data/processed/molecules.jsonl`: the same records with arrays and booleans preserved.
- `data/processed/candidate_sites.csv`: one row per candidate atom, with group index, source-positive and expanded-positive flags, and hypothetical marker product.
- `data/processed/summary.json`: machine-readable counts.
- `molecule_browser.html`: offline searchable structure/annotation browser covering all 535 records.
- `symmetry_examples.png`: selected examples with source atom indices.
- `docs/data_dictionary.md` and `docs/benchmark_decisions.md`: semantics, limitations, and pending decisions.

Run `python scripts/download_data.py --verify-only`, `python scripts/validate_benchmark.py`, `python scripts/build_reports.py`, and `python -m pytest -q` in the configured environment. No encoder training, 3D conformer generation or experimental relabeling has been performed.

Sources: [original EAS paper](https://doi.org/10.1039/C7SC04156J), [RegioSQM20](https://doi.org/10.1186/s13321-021-00490-7), [original data](https://github.com/jensengroup/db-regioselectivity), [RegioSQM20 inputs](https://github.com/jensengroup/SI_RegioSQM20).
'''
    (REPORTS / "data_audit.md").write_text(text, encoding="utf-8")
    with (REPORTS / "reaction_by_molecule_family.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        families = sorted(summary['molecule_family_counts'])
        writer.writerow(['reaction_class']+families+['total'])
        for category in summary['reaction_class_counts']:
            selected = [r for r in records if r['reaction_class']==category]
            counts = Counter(r['molecule_family'] for r in selected)
            writer.writerow([category]+[counts[f] for f in families]+[len(selected)])
    print(f"Built audit, searchable browser, cross-tabulation, and symmetry examples in {REPORTS}")


if __name__ == "__main__":
    main()
