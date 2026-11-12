# Initial benchmark decisions

## Adopted for the data audit

1. Preserve all upstream rows, atom indices and raw metadata. Never repair a reference or reagent by guessing.
2. Use exact single-substitution product identity as the primary definition of equivalent candidate positions. Preserve specified stereochemistry, isotopes and the supplied tautomer.
3. Preserve raw positive indices and add expanded positive sets in separate fields.
4. Rank distinct outcome groups. Average atom features/logits within a group to avoid an automatic symmetry multiplicity advantage.
5. Keep molecules as evaluation units; group duplicate substrates before splitting.
6. Exclude no-choice and multiple-distinct-outcome records from the initial single-outcome ranking task. Keep them available for separate analyses.
7. Use the rule-filtered 299-record conservative bromination cohort as a proposed pilot, with the 495-record broad single-outcome cohort as a sensitivity analysis.
8. Use SMILES indexing as authoritative. Do not substitute the SDF structure for comp485 silently; its alkene stereochemistry differs. Supplementary coordinates are 2D throughout.

## Specific cases for review

| Case | Observation | Proposed handling |
|---|---|---|
| comp18, comp21 | One label misses an exactly equivalent counterpart | Accept either atom as the same outcome; retain the original label |
| comp488 | Nitrobenzene has five candidate C-H atoms but three distinct outcomes | Rank the three outcomes, not five atom positions |
| comp492 | Source labels two distinct outcome groups; notes say `2xBr` | Review sequential/multiple substitution; exclude from first single-outcome model |
| comp274 | Two annotated para positions in equivalent phenyl groups; notes say repeated iodination | Symmetry alone does not turn a multiple-substitution experiment into a monofunctionalization example |
| comp158 | Two genuinely distinct annotated outcomes and paired yields | Retain for literature review; do not invent an equal-probability target |
| comp13 / comp598 | Same substrate and same marker outcome; NBS vs NCS notes | One substrate group across all splits, separate experimental records |
| comp485 | Source SMILES lacks alkene stereo specified in SDF; only one eligible site | Flag the source discrepancy and keep out of the nontrivial ranking metric |
| 69 empty-reference rows | The upstream table leaves the field blank | Preserve/flag; do not forward-fill from the previous citation |
| Ambiguous reagent text | Examples `12`, `C12`, `ICI`, `Nal` | Do not autocorrect to I2/Cl2/ICl/NaI without source review |

## Still requiring an executive decision

- Is the first scientific claim about traditional substrate-driven EAS susceptibility, or condition-specific bromination? The current source is stronger for the former; full conditions are not normalized.
- Should the 69 missing-reference rows be restored to a training cohort after literature recovery? This would enlarge the conservative cohort, but reference recovery must precede that decision.
- Should symmetry-expanded positives be pooled before or after the site head? Both should agree for a linear head; averaging logits is the specified baseline.
- What independent second dataset will be used to confirm a result? The current dataset is too small and narrow for the final broad claim.
- How will source/patent-family groups and fingerprint clusters be balanced for grouped cross-validation? No final test split is committed yet to avoid freezing it before these issues are resolved.

The current policy does not merge chemically similar sites, resonance-related sites in different tautomers, near-symmetric 3D conformations, or stereochemically distinct products. Atropisomerism and other stereo not encoded in the source cannot be recovered by this audit.

