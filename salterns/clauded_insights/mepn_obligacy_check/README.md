# Are the three MePn-obligate models a modelling artefact? (checked 2026-11-04/05)

**Question.** In the original saltern metabolic models (fliu's 24 MePn-consumer models, `cliff_mags/mac/salterns`), three models cannot grow without methylphosphonate (MePn): Roseovarius.9, Roseovarius.11 and Albimonas.1. Between them they produce 72 % of the simulated methane. Is this obligacy real biology or a gapfilling/template artefact?

**Answer: an artefact.** It was confirmed by a model analysis (`analysis/`) and an independent biology check (`biology_check/`), which also reproduced every growth value and the cell-27 methane values with GLPK.

- None of the three has a biological reason to need MePn.
- After a single-reaction repair, each grows without MePn and gives the same ratio-bound methane as every other model.
- Model-vs-measured methane, on the consistent areal basis of `../modelling_consistent_basis_r232.py`, falls from **0.564 to 0.190** of measured at r = 0.07 (MePn/Pi), and from **0.999 to 0.721** at r = 0.33.
- At r = 0.07 the three models' share of simulated methane falls from 72 % to 18 %. Spearman ρ goes from −0.20 to −0.09; neither is significant.
- "About half of the measured methane" is overstated: it is roughly one fifth at r = 0.07 and about 0.7 at r = 0.33.

**In the manuscript review (2026-11-08).** The corrected numbers are written into `../../Metabolic Manuscript - R232 review (tracked changes).docx` as tracked edits dated 2026-11-08 (see the review notes' 2026-11-08 section). They come from the `MePn_*_20_repaired` sets of `../modelling_consistent_basis_r232.py`, which asserts them against `analysis/step5_consistent_basis.json`. Figure 7 (`../figure7_r232.py`) shows the repaired models filled and the as-built models open. An independent re-check on 2026-11-08, with the two repair reactions written by hand instead of taken from the template, reproduced the obligacy, the rescue and the per-biomass methane.

## Why each model needs MePn

| Model (MAG) | Without MePn | Cause in the model | Gene evidence (MAG) |
|---|---|---|---|
| Roseovarius.9 (R2_B_D2 metabat.55; ~99 % complete) | glycine blocked (also Thr and Ile, which the model makes from glycine) | **One-carbon surplus with no outlet.** GlyA releases CH2-THF, but the model has no PurU, Fhs is forward-only in the template, glycine cleavage runs one way, and there is no homoserine kinase. MePn breakdown through fliu's hand-built C-P lyase (rxs43120) makes 5'-deoxyadenosine, which the template drain SK_cpd03091 removes. Each MePn therefore forces a purine to be rebuilt, which works as a one-carbon sink. Closing that drain makes the model unable to grow with or without MePn. | Has fhs (CDS.706), FDH, GcvPTH, MtaP and Apt. No homoserine kinase in any of 16 Rhodobacteraceae MAGs, nor in the KEGG genomes of R. pomeroyi, D. shibae, R. sphaeroides and Jannaschia: a lineage-wide pathway hole, not an auxotrophy. |
| Albimonas.1 (R2_restored_C_black metabat.8; ~99 %) | as Roseovarius.9 | Same mechanism as Roseovarius.9. | **Its C-P lyase genes are broken:** PhnJ, PhnI, PhnH and PhnK are frameshifted, there is no PhnCDE transporter, and the locus sits on a low-coverage contig. It may not use MePn at all. |
| Roseovarius.11 (R2_C_D2 concoct.45; **81.5 % complete, 19.9 % contamination**, below MIMAG medium) | PRPP blocked (so ATP, NAD, CoA, nucleotides, His, Trp) | **No PRPP synthetase.** The MePn route (PhnGHIL → PhnM → C-P lyase → PhnP → PhnN) is its only PRPP source. It works only with AMP nucleosidase running backwards, a hydrolysis that ModelSEED lists as forward-only but the template sets reversible. Forcing it forward-only leaves the model unable to grow at all. | Prs absent in this MAG only, plus the 79.5 %-complete Phaeodactylibacter.1; present in 22 of 24 MePn MAGs. Most likely missed by binning. |

**Why gapfilling did not catch it.** fliu gapfilled on a Pi-free medium where MePn was the only P source (`gapfill.ipynb` cells 92–93). There, forced MePn uptake supplied the one-carbon sink and the PRPP for free, so these gaps were never filled. In Roseovarius.9 and .11 no gapfilled reaction carries the MePn-dependent route; in Albimonas.1 one does (rxn01643).

## Repairs

Each rescuing repair gives exactly the same methane: 0.043209 per unit biomass at r = 0.07 and 0.163880 at r = 0.33, which is r/(1+r) × 0.660485 with MePn/Pi = r. None creates free ATP.

- **Roseovarius.9, Albimonas.1:**
  - homoserine kinase (rxn01300; the minimal repair, 1 of 51 single additions that rescue);
  - or the gene-supported options: Fhs opened toward formate (rxn00690), or glycine cleavage reversed.
  - In *R. pomeroyi*, Fhs releases formate (Lidbury et al. 2015).
- **Roseovarius.11:** PRPP synthetase (rxn00770), which is also MSGapfill's choice.
- **Across all 29 metabolic repair variants** (single additions, the annotation-supported reversals, and MSGapfill sets), methane is identical.
- **Artificial "export" repairs** (excreting a cofactor or nucleotide) can leave MePn as the cheapest one-carbon sink. They raise the total to at most 0.305 at r = 0.07 and 0.804 at r = 0.33.

**Consequence for the modelling claim.** After repair, all 20 simulated models give the same methane per unit biomass, fixed by the template biomass P content. The model-vs-measured comparison then tests only the MePn-consumer biomass estimates and the assumed r, not the genomes.

## Corrected totals (consistent areal basis; `analysis/step5_consistent_basis.json`)

| | r = 0.07 before | r = 0.07 after | r = 0.33 before | r = 0.33 after |
|---|---|---|---|---|
| Model / measured, all 6 cores | 0.564 | **0.190** | 0.999 | **0.721** |
| Per-core range | 0.178–1.412 | 0.057–0.460 | 0.305–2.451 | 0.215–1.744 |
| Spearman ρ (p) | −0.20 (0.70) | −0.086 (0.87) | −0.086 (0.87) | −0.086 (0.87) |
| Share of the 3 models | 72.3 % | 18.0 % | 40.8 % | 18.0 % |

Other results:
- **Relative-abundance-only sensitivity:** 0.540 → 0.186 (r = 0.07) and 0.967 → 0.704 (r = 0.33).
- **Excluding Albimonas.1** (broken C-P lyase): 0.186 and 0.705.

## Caveats

- **Repair universe:** repairs were searched in the local ModelSEED GramNeg v6 template with fliu's cell-72 edits re-applied. fliu's own KBase template (GramNegModelTemplateV6_test, ws 12218) is not available locally.
- **Gapfilling not redone:** the drafts were not re-gapfilled from scratch on a Pi-containing, MePn-free medium. That would be the cleanest test.
- **Albimonas.1 locus:** its frameshifts were not checked against reads.
- **Thiolinea.1:** its stored methane is not unique at the optimum. Its share is about 0.01 %.
- **Rebuilt model:** the scouting rebuild of Roseovarius.9 grows without MePn, but it cannot use MePn at all. It lacks PhnP, because the RAST name does not match the template role.

## Files

- `analysis/`: the model analysis.
  - Scripts `step0`–`step6` with their logs.
  - `summary.json`, with the per-model results, corrected totals and caveats.
  - `summary_per_model.csv`.
  - `step5_consistent_basis.json`, the corrected totals.
  - The scripts ran from the session scratchpad. They read the models and data in this repository and wrote nothing into it.
- `biology_check/`: the independent biology verifier's `notes.md`, logs `b1`–`b14` and scripts: gene evidence (RAST, KOfam/TIGRFAM HMMs, six-frame searches), KEGG reference genomes, model perturbations, the hand one-carbon balance, and the Albimonas.1 frameshift map.
- Not kept here (large; regenerable with the scripts): the repaired model JSONs (~95 MB) and `step4c_repair_and_rerun.json`.
