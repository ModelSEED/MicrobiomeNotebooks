# MePn obligacy of Roseovarius.9, Roseovarius.11 and Albimonas.1: biology and annotation check

2026-11-05. Read-only on the salterns repo. GLPK only, at most 8 threads, no KBase calls.
VB = /tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns/a24fbb2e-9d33-5e05-aeea-a07d777ccde7/scratchpad/mepn_verify/verify_bio
Scripts are in VB/scripts/b*.py. Each writes VB/b*.log (and .json).

## Verdict

The obligacy is a modelling artefact in all three models. There is no biological basis for any of these
organisms needing MePn for growth. The stored high MePn fluxes (0.4697, 0.5534 and 0.4697 mmol per unit
biomass) should be treated as artefacts in the manuscript.

Every biology-motivated repair I tested gives the same ratio-bound methane as the other agent's repairs:
0.043209 per unit biomass at r = 0.07 and 0.16388 at r = 0.33. The corrected totals (0.190 and 0.721) therefore
stand.

I disagree with the analysis on mechanism and gene support, not on the numbers. Points 3, 4 and 6 below give
the details.

## 1. Roseovarius.11: the PRPP block comes from an incomplete MAG and a hydrolase run backwards

- **MAG quality:** 81.45% complete, 19.88% contamination (guilds_r232.csv row Roseovarius.11). The contamination
  is above the MIMAG medium-quality limit of 10%.
- **Prs is missing from this MAG only:**
  - RAST has no ribose-phosphate pyrophosphokinase (b1_rast_grep.log line 3, "Prs" column = 0).
  - KOfam K00948 finds no hit (b2_hmm_kofam.log line 233).
  - A six-frame search of all 852 contigs finds only a 0.40x-threshold match to the shared PRTase domain, HMM
    205-311 (b4_sixframe.log line 55).
  - Prs is present in 22 of the 24 MePn MAGs. It is missing only in Roseovarius.11 and Phaeodactylibacter.1
    (79.5% complete), the two least complete MAGs (b2 lines 230-256).
  - Prs is essential and near-universal. At 81% completeness, the chance of missing any one gene is about 19%.
- **How the model makes PRPP:** it uses PhnI (rxn26456: ATP + MePn -> adenine + RPnTP) followed by PhnM, PhnJ,
  PhnP and PhnN, which gives PRPP. The adenine goes back onto ribose through rxn00131 (AMP nucleosidase) running
  backwards.
  - rxn00131 is a hydrolysis. The ModelSEED biochemistry database lists it as forward-only (reversibility '>').
    The local GramNeg V6 template sets it to '=' (b13_msdb_and_template.log lines 2 and 22).
  - With rxn00131 limited to the forward direction, Roseovarius.11 cannot grow at all, even with MePn: growth
    a/b/c = 0/0/0 (b6_bio_perturbations.log line 20).
  - So 100% of its MePn flux depends on a thermodynamically implausible reversal. That flux is 41.4% of all
    stored simulated methane at r = 0.07 (step5_consistent_basis.json, per_model_total: 1657.8 of 4000.7).
- **Known biochemistry does not rescue it:**
  - The C-P lyase pathway does end in PRPP (PhnN). But its ribose comes from ATP, so it is not a de novo source.
  - In E. coli, Δprs strains need NAD/NMN plus guanosine, uridine, His and Trp. Suppressors that use the PhnN route
    (R5P -> R1P -> R1,5P2 -> PRPP) regain only NAD prototrophy and "still require guanosine, uridine, histidine,
    and tryptophan" (Hove-Jensen et al. 2003, J Bacteriol 185:2793; PMC154390, fetched).
  - Roseovarius.11 also lacks phosphopentomutase (K01839: none, b2 line 233), so even that route is unavailable.
- **Repair:** add PRPP synthetase. That is a completeness fill, not a RAST-supported one. Growth becomes
  4.68/4.68/4.67 and methane becomes ratio-bound (b6 line 24).

## 2. Roseovarius.9 and Albimonas.1: the model has no way to dispose of excess one-carbon units

- **MAG quality:** both MAGs are near-complete with low contamination (98.87/4.01 and 99.0/2.4). Missing genes are
  not a completeness effect here.
- **No recognisable homoserine kinase, in these two MAGs or anywhere in the clade:**
  - RAST finds none (b1 lines 2 and 4).
  - None of the three KOfam HSK KOs (K00872, K02204, K02203) has a hit at threshold (b2 lines 3-31). The best
    sub-threshold hits are an APH-family "putative phosphotransferase" at 0.54x and 0.77x (b2 decisive section).
  - Six-frame search of the contigs also finds nothing at threshold (b4).
  - TIGRFAM thrB, thrB_alt and the "proposed homoserine kinase" TIGR02535 find nothing in any of the
    16 Rhodobacteraceae MAGs (b9_tigr_hsk.log).
  - KEGG assigns no gene to any HSK KO in four cultured Rhodobacteraceae: Ruegeria pomeroyi DSS-3, Dinoroseobacter
    shibae, Rhodobacter sphaeroides 2.4.1 and Jannaschia CCS1. Each of them does have Hom (K00003) and ThrC (K01733)
    (b10_kegg_reference_genomes.log). The controls work: E. coli thrB b0003 and P. aeruginosa thrB2 PA5495 are
    found.
  - This is a lineage-wide pathway hole in organisms that make threonine. The enzyme is unidentified or
    non-orthologous; it is not a real auxotrophy. Filling rxn01300 is a standard hole fill.
- **PurU is absent across the clade:** K01433 is found only in Saliniramus.1 among the 24 MAGs (b2 line 174ff),
  and R. pomeroyi has none in KEGG (b10).
  - In R. pomeroyi, Fhs "is predicted to convert formyl-H4F to formate". It is needed to oxidise the methyl groups
    from choline and glycine betaine, and an fhs null mutant "failed to grow on GBT" (Lidbury et al. 2015, Environ
    Microbiol; PMC4744692, fetched).
  - Fhs is present in both MAGs (Roseovarius.9 CDS.706, KOfam 1.05x; Albimonas.1 CDS.2814, 1.07x; b2 decisive
    section). In both models, rxn00690 is forward-only, (0,1000) (b5), as in the template.
- **Gene-supported single repairs exist.** The analysis says there are none.
  - P2, Fhs reversible: growth without MePn is 5.1456 for Roseovarius.9 (b6 line 3) and 5.1024 for Albimonas.1
    (b6 line 12).
  - P3, glycine cleavage (rxn06377 + rxn06600) reversible, with GcvPTH genes present: 5.1456 (line 4) and 5.1220
    (line 13). Caveat: the glycine-synthase direction is used natively mainly by anaerobes.
  - P5, HSK hole fill: 5.1251 and 5.1012 (lines 7 and 15).
  - All of them give ratio-bound methane: 0.043209 at r = 0.07 and 0.16388 at r = 0.33, with MePn/Pi = r.
  - step4a missed the Fhs and GCS reversals because it only reversed reactions whose template GapfillDirection
    allows it (analysis/step4a_single_repair_scan.py docstring lines 8-9). The template sets rxn00690 to '>' and
    rxn06377 to '>' (b13 line 23; b5 template section).
- **The 5'-deoxyadenosine sink makes the obligacy.** Closing SK_cpd03091_c0 stops all growth, with or without MePn,
  in both models: 0/0/0 (b6 lines 2 and 11).
  - This sink is a template drain (template_neg.drains includes cpd03091_c, (0,1000)). It is added after gapfilling by
    `_add_dm` (../fliu_notebooks_dump.txt lines 2821-2828, 2448-2474). Gapfilling used M_GENOME_SCALE, which has
    MePn and no Pi exchange (lines 2389-2416).
  - No model reaction salvages 5'-dA (b5). The template does have a 5'-dA nucleosidase, rxn13649, typed
    "gapfilling" with no roles.
  - In real cells, 5'-dA from radical-SAM enzymes is cleaved and the adenine is salvaged (Beaudoin et al. 2018 Nat
    Commun). Roseovarius.9 carries MtaP (CDS.1406, K00772 1.25x) three genes from Apt (CDS.1409, K00759 3.2x).
  - In the model, MePn therefore acts as a one-carbon sink only through the loss of purines.
- **The magnitude matches a one-carbon balance computed independently:**
  - Biomass stoichiometry only, no FBA (b11_c1_balance.log). With Thr and Ile made from glycine, the one-carbon
    surplus is +0.489 per unit biomass. The FBA MePn flux is 0.4697, about 4% lower.
  - With Thr made from aspartate, the surplus drops to +0.05.
- **No known biochemistry links MePn to glycine or one-carbon balance:**
  - PhnY/PhnZ gives glycine + Pi, but from 2-aminoethylphosphonate, not MePn (McSorley et al. 2012 JACS).
  - Glyphosate C-P lyase gives sarcosine, which would add one-carbon units.
  - MePn C-P lyase gives CH4 + PRcP + Met + 5'-dA.
  - C-P lyase belongs to the Pho regulon (phn derepression in Pho-constitutive mutants; Hove-Jensen 2003). So a
    MePn uptake of 2.5x Pi on a Pi-replete medium is also physiologically unexpected.
- **The model misses the ethylmalonyl-CoA pathway the genomes encode.** All three MAGs have the complete pathway
  (K14446, K14447, K14448, K14449, K14470, K08691 and K14451 all pass; b2 line 135ff). It is the
  Rhodobacteraceae route to glyoxylate without isocitrate lyase.
  - The template has the roles (e.g. "Crotonyl-CoA carboxylase/reductase, ethylmalonyl-CoA producing") but no
    reactions attached to them, so the pathway never enters the models (b13 lines 25-33).

## 3. Albimonas.1 may not degrade MePn at all (new finding)

- **The C-P lyase genes on NODE_10693 are broken by frameshifts:**
  - Each HMM matches as colinear pieces in different frames (b7_frameshift_map.log lines 1-28).
  - PhnH is in 2 frames, PhnI is split by 3 frameshifts, and PhnJ is in 5 pieces with at least 3 frameshifts and
    one in-frame stop. PhnK is also broken.
  - Controls: the Roseovarius.9 and Roseovarius.11 loci give single-frame PhnI and PhnJ at 524-605 bits (b7 lines
    29-60).
- **There is no intact copy anywhere in the MAG.** No KOfam K06163, K06164 or K06165 hit at threshold in the whole
  proteome (b2 line 271). The best six-frame hits are PhnJ 0.54x, PhnI 0.90x and PhnH 0.85x, all on the same
  contig (b4 lines 127-134). There is no PhnCDE transporter at threshold (b4).
- **The contig is a low-coverage outlier:** 3.50x against a MAG median of 5.04x. 75% GC. Three of the eight
  junctions sit at G-homopolymers (b8_locus_context.log).
- Assembly error or real pseudogene: I cannot tell which without reads. Either way, the genome as assembled does
  not support MePn use.
- **Effect on the numbers:** dropping Albimonas.1 after repair changes the corrected total from 0.190 to 0.186 at
  r = 0.07, and from 0.721 to 0.705 at r = 0.33 (step5_consistent_basis.json, per_model_total).
- Its MePn-consumer flag in guilds_r232.csv rests on RAST names attached to these fragments.

## 4. Annotation details where I disagree with the analysis

- **"GlyA EC 2.1.2.1 x4 (CDS.1374/1553/1665/2886)" in Roseovarius.9:** all four are 60-68 aa fragments at contig
  edges, covering HMM 125-172 (x2) and 477-544 (x2) of 623. The MAG has no full-length glyA, even in a six-frame
  search (b12_shmt_r9.log). SHMT is surely present in vivo; the GPR rests on fragments.
- **PhnN in Albimonas.1:** CDS.3804 ("ATP-binding protein PhnN") passes KOfam K05774 at 1.95x (b3). So rxn04651 is
  gene-supported even though its GPR is empty.

## 5. Reproduction (GLPK, cell-5 medium checked against `git show 22f4629`)

Growth a/b/c, with methane per unit biomass in cell 27:

| Model | Growth a / b / c | Methane per unit biomass | MePn/Pi |
|---|---|---|---|
| Roseovarius.9 | 4.6998 / 0 / 4.6981 | 0.469737 | 2.463 |
| Albimonas.1 | 4.6629 / 0 / 4.6544 | 0.469737 | 2.463 |
| Roseovarius.11 | 4.6734 / 0 / 4.6696 | 0.553445 | 5.17 |

Source: b6 lines 1, 10 and 18. These match the analysis.

## 6. Energy-loop check of the repairs

I ran every repair on the cell-5 medium without lactate and without MePn, with bio1 fixed at 0, and maximised
ATPM. The maximum is 0 for P0, P2 (Fhs reversible), P3 (GCS reversible), P5 (HSK) and P6 (Prs) in all three
models (b14_loop_check.log). None of these repairs creates free ATP.

## 7. Limits

- **Literature evidence:** Hove-Jensen 2003 and Lidbury 2015 were read from fetched PMC pages. Kamat 2013 (PhnJ is
  a radical-SAM enzyme that produces 5'-dA and Met), Beaudoin 2018, McSorley 2012 and GCS reversibility come from
  abstracts or search snippets only.
- **KEGG reference genomes:** "no gene assigned" is a database statement. Its weight comes from these being
  curated genomes of organisms that make threonine.
- **Albimonas.1 phn locus:** not checked against reads. The 2026 reassembly is the place to resolve it.
- **Not done:**
  - I did not re-gapfill from the drafts.
  - I did not test whether adding the missing EMC core reactions alone would repair the models.
  - I did not compute rho for the variant without Albimonas.1.
