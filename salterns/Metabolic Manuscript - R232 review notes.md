# R232 review of the saltern manuscript (2026-10-21, revised 2026-10-26)

**File:** `Metabolic Manuscript - R232 review (tracked changes).docx`. It is a copy of `Metabolic Manuscript.docx` (2 July 2026) in which every edit is a Word tracked change by "Claude (R232 review)". "Reject all" restores the original text, tables and figures.
**Accompanying file:** `Table S5 - MAG taxonomy, quality and guilds (GTDB R232).xlsx`, one row per MAG (310).
**Rebuild:** `clauded_insights/make_r232_review_docx.py`. Its inputs are `correlational/r232/` and `data/`, plus `clauded_insights/r232_review_edits_2026-10-26.json` (39 tracked edits and comments, each with its rationale and source) and `clauded_insights/figure7_r232.png`.

> **The draft may be out of date.** The "Metabolic Manuscript" Google Doc in the "JGI Salterns Modeling" Drive folder was modified on 2026-08-27, after this .docx. Edits made there since July will need to be reconciled with these changes.

> **Author fixes (2026-10-21).** The author reports fixing the three draft errors (Table 3, P69, P35) in their own copy, which is not the local .docx or the Drive Google Doc. The tracked edits at P35 and P69 now duplicate those fixes. The R232 Table 3 replacement is a separate update. A hand fix of the old Table 3 also needs new Unrestored values, because that column pooled the R2A cores (see the Table 3 row below).

## 2026-10-26 revision

What this revision changed, and where the evidence is:

- **Modelling text and numbers (P18, P36–P40, P69–P77, P95, P99, Figure 7 caption; comment on P100).** The headline finding of this revision: **on a dimensionally consistent basis, the simulated methane does not track the measured emissions among the six unrestored cores.**
  - **Why the old fits looked good.** The stored simulations fixed each model's biomass to relative abundance ÷ DNA, averaged over depths. Two depths have almost no DNA (R2_B_D1 and R2_C_D2, 0.16 µg g⁻¹), so every simulated value carries a 1/DNA term. All earlier fits came from that term:
    - the draft's R² = 0.978 / 0.968 / 0.862, with both axes divided by DNA (0.966 matches no stored fit);
    - an intermediate "areal" correction from July (`figure7_corrected.py`, ρ = 0.60), where multiplying by the core's DNA did not cancel the term.

    The factor VW DNA × mean(1/DNA) alone gives ρ = 0.60 against measured methane.
  - **Consistent basis.** Each MAG's flux per unit of simulated abundance is rescaled to its biomass in the core: relative abundance × DNA × layer thickness, summed over D1 and D2.
  - **MePn models:**
    - ρ = −0.20, exact p = 0.71, n = 6, R² = 0.08 at ratio 0.07, with 56 % of measured emissions in total (18–141 % by core);
    - ρ = −0.09, R² = 0.06 at 0.33, with 100 % in total (30–245 %);
    - so the models can supply methane of the measured magnitude but not the pattern among cores.
  - **Betaine models** (TMA mostly at the uptake bound, so converted as before): ρ = −0.09, R² = 0.10.
  - **Other numbers:**
    - CH₄ : MePn is exactly 1 : 1, not 1.4 : 1.
    - Figure S7 (DNA vs CH₄) is R² = 0.22, ρ = 0.52, p = 0.08, so "significantly" became "weakly and not significantly".
    - The "unconstrained" R² = 0.970 run is not stored. It matches the three MAGs (*Roseovarius*.9, .11 and *Albimonas*.1) for which the ratio constraint does not bind; they supply 72 % of the simulated methane.
  - **Sources:** `clauded_insights/modelling_consistent_basis_r232.py` → `.json`. `modelling_numbers_r232.json` traces the draft's values but its "areal" basis is superseded.
- **Equations P70, P72, P75 (Office Math).** Each is replaced as a tracked deletion of the old equation plus a tracked insertion of the new fit:
  - P70: −0.20x + 907;
  - P72: −0.32x + 1563;
  - P75: −0.71x + 3192.

  The units drop the 1/(µg DNA g⁻¹ soil) term. Word shows these as tracked changes; LibreOffice ignores revisions on equations and shows both versions.
- **Betaine models.** These are restricted to the nine *grdH* carriers. A re-simulation reproduces the stored fluxes exactly, because no model fix was supported by the annotations. The fit is ρ = −0.09, p = 0.87, R² = 0.10. Other findings:
  - The betaine reductase reaction (rxn17220_c0) was added by hand to all 15 models with no gene association.
  - 33 of 54 fluxes sit at the 0.1 uptake bound, so modelled betaine methane tracks DNA density (R² = 0.97), not abundance.
  - Only *Halanaerobium*.2 is forced to zero TMA. The old review comment that "the *Halanaerobium* MAGs make no TMA" was wrong for *Halanaerobium*.1 and .3.
  - Evidence: `correlational/betaine_grdH_r232/`.
- **Figure 7 replaced** by `clauded_insights/figure7_r232.png` (`figure7_r232.py`). It is a single panel on the consistent areal basis with a 1:1 line, MePn (20 models) and betaine (9 *grdH* models). The old panel b is superseded by Figure 9, Figure S8 and Table 5.
- **Other corrections from the verification pass:**
  - Methods cited Figure S9 for gap-filling; it now cites Figure S10.
  - A Discussion sentence attributing the community-wide correlation to "higher metabolic activity" is deleted; the revised Results attribute it to the site gradient.
  - Table 4's last column is now "In the curated modeling set but lacking the marker", noting which of those MAGs were simulated.
  - The new co-occurrence wording is tightened (filters, edge counts, maxT).
  - The deleted and new Table 3 are separated by a paragraph so that Word does not merge them.
  - A comment on the Conclusions suggests hedging "primary source … precise mechanism".
- **Halomonas isolates (P47, P67, Figure 5 caption).** GTDB-Tk 2.8.0 / R232 keeps both isolates in *Halomonas*.
  - I1 is *Halomonas* sp. D, the C–P lyase carrier. It is deposited as GCA_032334135.1 and is its own species cluster, *Halomonas* sp032334135.
  - I4 is *H. saccharevitans*: 98.5 % ANI to the type strain.
  - R232 places 11 of the 23 Figure 5 references labelled *Halomonas* in other genera. Three of the six MPn-positive references remain *Halomonas*.
  - Evidence: `MAG_QC/gtdbtk_isolates_r232/` (`run_gtdbtk_isolates.sh`, README).
- **Co-occurrence paragraph and Figure 8 caption (R232 insertions) qualified** by `correlational/r232/prefilter_sensitivity/within_REPORT.md`:
  - The network's size and modularity fall within Graphical Lasso fits to permuted residuals; only the excess of edges between related taxa does not.
  - The key taxa have no edges because the lasso acts on unstandardized residuals (degree tracks residual variance, ρ = 0.90).
  - Relaxing the prevalence filter gives no edge to a guild taxon that lacked one; new edges mostly join rarely detected taxa.
  - At most two pairs reach the maxT threshold, depending on the number of permutations.
- **Smaller corrections:** "Saligentibacter" → *Salegentibacter* (P67), and in the Figure S11 caption "red" → "orange".

## What the new analysis changed

- **Quality (CheckM):** the CheckM re-run reproduced the original values exactly. 275 of 310 MAGs meet MIMAG medium quality; 35 exceed 10 % contamination, including *Roseovarius*.11 and .12.
- **Taxonomy (GTDB-Tk 2.8.0, R232):**
  - 280 bacterial and 30 archaeal MAGs.
  - 290 MAGs resolved to genus (originally 201) and 176 to species (originally 19).
  - The genera central to the MPn argument are unchanged. One *Halomonas* MAG is now *Vreelandella*.
- **Gene-based guilds:**
  - MPn consumers (*phnJ*): 22, of the 27 in the modeling set.
  - Betaine reducers (*grdH*): 9, of the 15 in the modeling set.
  - Methanogens (*mcrA*): 1, *Methanolobus* sp020055655.
  - *cutC* is not annotated in any MAG.
  - The before/after membership is in `correlational/r232/guild_membership_2020_vs_r232.csv` / `.md`.
- **Guild abundance vs methane:**
  - On the areal basis, every guild correlates across the 12 cores, but none beats random MAG sets.
  - On the DNA-normalized basis, only MPn consumers stand out (ρ = 0.92, above 99 % of random sets).
  - Nothing is significant within the six unrestored cores.
- **Methanogen:** *Methanolobus*.1 occurs only in the unrestored cores, and within them it does not track methane (ρ = −0.54, p = 0.27).
- **Co-occurrence:** the archaea now enter the network (169 taxa, 13 archaeal). No pair passes FDR. The network has 6 Louvain modules (Q = 0.747), including one haloarchaeal module. Its size and modularity are within the permutation null, and the key taxa have no edges, for the variance reason above.

## Tracked text changes (P-numbers refer to the original paragraphs)

| Where | Change |
|---|---|
| Abstract (P18) | Modelling clause: "reproduced the rank order, though not the magnitude" |
| Methods, after Metagenomics (P30) | **New subsection** "MAG quality, taxonomy, and functional guilds": CheckM, GTDB-Tk R232, iterativeID naming, gene rules for the guilds |
| Methods, Metabolic Modeling (P35) | "dereplicated medium and high quality MAGs" → "310 dereplicated MAGs (275 meet MIMAG medium quality; Table S5)" |
| Methods (P36) | Phosphate ratio = a minimum MePn : Pi uptake ratio; 20 of 27 MePn consumers simulated (why); betaine models = the nine *grdH* carriers; 0.1 = maximum uptake (mmol g~DW~⁻¹ h⁻¹, unit to confirm) |
| Methods (P38, P40) | Regression and Spearman; conversion with each core's DNA density (areal basis) |
| Methods, isolates (P47) | GTDB-Tk v2.8.0 / R232 reclassification sentence; BioProject and assembly accessions |
| Methods, Analyses (P52) | Replaced: correlation design (12 cores and 6 unrestored, two bases, 4,000-set random comparison, BH-FDR) |
| Methods, Analyses (after P53) | **New paragraph:** co-occurrence methods (CLR, within-site residuals, permutation FDR, graphical lasso, Louvain) |
| Results, after Metagenomics (P57) | **New subsection** "MAG quality, taxonomy, and functional guilds" |
| Results, culturing (P67) | R232 placement of both isolates and of the six MPn-positive references; *Salegentibacter* spelling |
| Results, modeling (P69–P77) | Areal-basis statistics replace every old R² (see above); equations P70/P72/P75 swapped; P76 rewritten for the nine *grdH* betaine models; P77 stoichiometry 1 : 1 |
| Results (before P78) | **New heading** "Guild abundances and methane flux" |
| Results (P78) | Replaced: the old R² values (0.974/0.966/0.852/0.798) give way to the R232 guild correlations and random-set comparison |
| Results (after P78) | **New paragraphs:** methanogen vs methane; new heading "Co-occurrence" with the network results and their permutation/variance caveats |
| Discussion (P83) | Inserted sentence: the only *mcrA* MAG is *Methanolobus* sp020055655; it is restricted to the unrestored salterns and does not track methane |
| Discussion (P86) | Inserted sentence: the *Nitrosopumilus* MAGs occur only in the reference wetland, and RAST annotated no *mpnS* in any MAG |
| Discussion (P95) | First sentence replaced (MPn specificity only on the DNA-normalized basis); modelling sentences replaced (rank order, 1 : 1, betaine models not correlated; models translate abundance into methane) |
| Conclusions (P99) | "Secondly…" and "Thirdly…" sentences replaced to match |
| Data Availability (P102) | Added the repository location of the MAG QC and correlation analyses |

## Figures and tables

| Item | Status |
|---|---|
| Table 3 | **Replaced.** Recomputed with R232 families for all 310 MAGs (archaea included; † marks archaeal families) and correct unbinned shares (52 / 88 / 89 %). **Cause of the old table's errors:** it came from `correlational/abundances.ipynb` cell 1, which had three faults:<br>• `'R1' in col or 'R2' in col` also matches the six R2A cores, so "unrestored" pooled R1 + R2 + R2A.<br>• The 30 archaea had no taxonomy entry and were filled as `unbinned_reads`, so they were dropped with the unbinned reads; "68 %" is unbinned plus archaea over those 18 cores.<br>• The notebook outputs unrestored / reference / restored, which was copied under the headers Unrestored / Restored / Reference.<br>The notebook is now fixed (2026-10-21) and prints the corrected 2020-family values. They are also in `Table 3 - corrected values (2020 GTDB families).csv`. |
| Table 4 | **New.** Functional guilds: marker gene, members, and what was dropped from the modeling sets. |
| Table 5 | **New.** Guild–methane Spearman correlations on both bases, 12 cores and 6 unrestored, with the random-set fractions. |
| Figure 5 | Caption: R232 genus placement of the reference genomes labelled *Halomonas* (comment on the one inferred genome). |
| Figure 7 | **Replaced** (2026-10-26) by a single panel on the areal basis (`clauded_insights/figure7_r232.png`); caption rewritten. |
| Figure 8 | **Replaced.** Both old panels and the caption are deleted. Panel b (the Louvain view) was re-rendered with R232 labels and archaea and is now the single panel; panel a (the plain graphical-lasso layout of the same edges) was not re-rendered. |
| Figure 9 | **New.** Guild abundance vs methane, areal basis (`correlational/r232/guild_methane_raw_r232.png`). |
| Figure S8 | **Replaced** with the DNA-normalized version of Figure 9; the stray caption fragment (old P151) is deleted. |
| Figure S11 | Not re-rendered. Caption colour fixed; comments explain the misdescribed asterisk rule and what regeneration needs (`correlational/r232/FigureS11_regeneration.md`). |
| Table S5 | **New caption.** Content is in the Excel file. |
| Table S6 | **New.** Per-core methane flux, *Methanolobus*.1 abundance, and total archaeal abundance. |

## Comments (18, by "Claude (R232 review)")

1. **Start of document:** summary of this review and revision.
2. **Abstract:** the modelling clause is now qualified; consider qualifying any guild-correlation claim too.
3. **Methods P47:** KBase narrative 139158 is not readable by the author's account and has no public read entry, so it does not appear to be public.
4. **Methods P47:** reference 68 is the GTDB-Tk v1 paper; add the GTDB-Tk v2 paper.
5. **Results P67:** "*Halomonas* sp. D" is isolate I1 (NCBI "Halomonas sp. I1"); use one name.
6. **Results P67:** the C–P lyase *Halomonas* is supported by the isolate, not by the MAGs (neither Halomonadaceae MAG carries *phn* genes).
7. **Results P69, P71:** what the replaced equations were and why the old R² (and the intermediate ρ = 0.60) were inflated by the 1/DNA term.
8. **Results P74:** the replaced betaine equation and the betaine-model caveats (hand-added reaction, O₂ dependence, non-unique TMA optima, contamination).
9. **Discussion P83:** reconcile *Ca.* Methanosalis SBSPR1A (ref. 58) with GTDB's *Methanolobus* sp020055655.
10. **Discussion P93:** the 2021 *Halomonas* species count predates the split into *Halomonas* and *Vreelandella*.
11. **Table 3:** the three confirmed causes of the old table's errors.
12. **Figure 5 caption:** the R232 counts assume the *H. campaniensis* tip is strain SK03.
13. **Figure 8 caption:** the figure is now a single panel; the old panel a (plain glasso layout) was not re-rendered.
14. **Figure S11 caption:** what regeneration needs, and the misdescribed asterisk rule.
15. **Table S3:** not in the document; consider adding the isolates' genome-based names and accessions.
16. **Conclusions:** consider hedging "primary source … precise mechanism" in light of Figure 7 and Table 5.

(Some topics carry two comments, which makes 18 in total.) Resolved and removed in this revision: the isolate-reclassification request, the modelling-R² comment, the betaine-model comment, and the Figure 7 panel-b comment.

## Not changed

- **FBA simulations of the MPn models** were not re-run. They use the depth-mean DNA-normalised abundances of commit 22f4629, not the current DNA-weighted merge; since each model's methane is proportional to its abundance, a re-run would only rescale them.
- **Figure S11:** not re-rendered (see comment 14 and `correlational/r232/FigureS11_regeneration.md`).
- **Figure 5:** tip labels are the original KBase names (the tree object in narrative 139158 could not be read).
- **Citation formatting:** citations in inserted text, e.g. (51), (53), (67), (68), are plain text, not reference-manager links. Relink them if you use one.
