# R232 review of the saltern manuscript (2026-10-21)

**File:** `Metabolic Manuscript - R232 review (tracked changes).docx`. It is a copy of `Metabolic Manuscript.docx` (2 July 2026) in which every edit is a Word tracked change by "Claude (R232 review)". "Reject all" restores the original text, tables and figures.
**Accompanying file:** `Table S5 - MAG taxonomy, quality and guilds (GTDB R232).xlsx`, one row per MAG (310).
**Rebuild:** `clauded_insights/make_r232_review_docx.py`. Its inputs are `correlational/r232/` and `data/`.

> **The draft may be out of date.** The "Metabolic Manuscript" Google Doc in the "JGI Salterns Modeling" Drive folder was modified on 2026-08-27, after this .docx. Edits made there since July will need to be reconciled with these changes.

> **Author fixes (2026-10-21).** The author reports fixing the three draft errors (Table 3, P69, P35) in their own copy, which is not the local .docx or the Drive Google Doc. The tracked edits at P35 and P69, and comment 6, now duplicate those fixes. The R232 Table 3 replacement is a separate update. A hand fix of the old Table 3 also needs new Unrestored values, because that column pooled the R2A cores (see the Table 3 row below).

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
- **Guild abundance vs methane:**
  - On the areal basis, every guild correlates across the 12 cores, but none beats random MAG sets.
  - On the DNA-normalized basis, only MPn consumers stand out (ρ = 0.92, above 99 % of random sets).
  - Nothing is significant within the six unrestored cores.
- **Methanogen:** *Methanolobus*.1 occurs only in the unrestored cores, and within them it does not track methane (ρ = −0.54, p = 0.27).
- **Co-occurrence:** the archaea now enter the network (169 taxa, 13 archaeal). No pair passes FDR. The network has 6 Louvain modules (Q = 0.747), including one haloarchaeal module, and the key taxa have no edges.

## Tracked text changes (P-numbers refer to the original paragraphs)

| Where | Change |
|---|---|
| Methods, after Metagenomics (P30) | **New subsection** "MAG quality, taxonomy, and functional guilds": CheckM, GTDB-Tk R232, iterativeID naming, gene rules for the guilds |
| Methods, Metabolic Modeling (P35) | "dereplicated medium and high quality MAGs" → "310 dereplicated MAGs (275 meet MIMAG medium quality; Table S5)" |
| Methods, Analyses (P52) | Replaced: correlation design (12 cores and 6 unrestored, two bases, 4,000-set random comparison, BH-FDR) |
| Methods, Analyses (after P53) | **New paragraph:** co-occurrence methods (CLR, within-site residuals, permutation FDR, graphical lasso, Louvain) |
| Results, after Metagenomics (P57) | **New subsection** "MAG quality, taxonomy, and functional guilds" |
| Results, modeling (P69) | "(Figure S8)" → "(Figure 7a)": S8 is no longer the modeling figure |
| Results (before P78) | **New heading** "Guild abundances and methane flux" |
| Results (P78) | Replaced: the old R² values (0.974/0.966/0.852/0.798) give way to the R232 guild correlations and random-set comparison |
| Results (after P78) | **New paragraphs:** methanogen vs methane; new heading "Co-occurrence" with the network results |
| Discussion (P83) | Inserted sentence: the only *mcrA* MAG is *Methanolobus* sp020055655; it is restricted to the unrestored salterns and does not track methane |
| Discussion (P86) | Inserted sentence: the *Nitrosopumilus* MAGs occur only in the reference wetland, and RAST annotated no *mpnS* in any MAG |
| Discussion (P95) | First sentence replaced: MPn specificity holds only on the DNA-normalized basis |
| Conclusions (P99) | "Secondly…" sentence replaced to match |
| Data Availability (P102) | Added the repository location of the MAG QC and correlation analyses |

## Figures and tables

| Item | Status |
|---|---|
| Table 3 | **Replaced.** Recomputed with R232 families for all 310 MAGs (archaea included; † marks archaeal families) and correct unbinned shares (52 / 88 / 89 %). **Cause of the old table's errors (confirmed after the docx was built):** it came from `correlational/abundances.ipynb` cell 1. There, `'R1' in col or 'R2' in col` also matches the six R2A cores, so "unrestored" pooled R1 + R2 + R2A. The 30 archaea had no taxonomy entry and were filled as `unbinned_reads`, so they were dropped with the unbinned reads; "68 %" is unbinned plus archaea over those 18 cores. The notebook outputs unrestored / reference / restored, which was copied under the headers Unrestored / Restored / Reference. Rebuilt under these three faults, the old values match to within 5 × 10⁻⁷. The notebook is now fixed and prints the corrected 2020-family values. They are also in `Table 3 - corrected values (2020 GTDB families).csv` (`clauded_insights/table3_2020_corrected.py`) alongside the old values. Haloferacaceae (5.9 %) and Haloarculaceae (2.5 %) enter the top 10. "Undetermined" (all MAGs without a family) drops from 8.2 % to 6.8 %, and Flavobacteriaceae drops from 2.6 % to 1.9 %. |
| Table 4 | **New.** Functional guilds: marker gene, members, and what was dropped from the modeling sets. |
| Table 5 | **New.** Guild–methane Spearman correlations on both bases, 12 cores and 6 unrestored, with the random-set fractions. |
| Figure 8 | **Replaced.** Both old panels and the caption are deleted. Panel b (the Louvain view) was re-rendered with R232 labels and archaea and is now the single panel; panel a (the plain graphical-lasso layout of the same edges) was not re-rendered. |
| Figure 9 | **New.** Guild abundance vs methane, areal basis (`correlational/r232/guild_methane_raw_r232.png`). |
| Figure S8 | **Replaced** with the DNA-normalized version of Figure 9; the stray caption fragment (old P151) is deleted. |
| Table S5 | **New caption.** Content is in the Excel file. |
| Table S6 | **New.** Per-core methane flux, *Methanolobus*.1 abundance, and total archaeal abundance. |

## Comments added (flagged, not changed)

1. **Start of document:** summary of this review.
2. **Abstract:** consider qualifying the modeling/correlation claim.
3. **Methods, isolates (GTDB-Tk v1.7.0):** reclassify the two *Halomonas* isolates with R232, since GTDB has split *Halomonas*.
4. **Results P69:** the modeling R² values (0.966 / 0.968 / 0.970 / 0.862; Discussion 0.978) still don't match Figure 7. They are corrected in the July file `clauded_insights/Metabolic Manuscript (tracked changes).docx` and were not changed here.
5. **Results P76:** the *Halanaerobium* MAGs carry *grdH* but make no TMA in the models, and 6 of the 15 betaine models lack *grdH*. Consider re-running the betaine models on the 9 *grdH* carriers.
6. **Table 3:** the 2020 column swap, the omitted archaea, and the 68 % vs 52 % unbinned discrepancy. The comment says the swap "appears"; its cause has since been confirmed (see the Table 3 row).
7. **Figure 7 caption:** panel b is superseded by Figure 9, Figure S8 and Table 5.
8. **Discussion P83:** reconcile *Ca.* Methanosalis SBSPR1A (ref. 58) with GTDB's *Methanolobus* sp020055655.
9. **Figure 8 caption:** the figure is now a single panel; the old panel a (plain glasso layout) was not re-rendered.
10. **Figure S11 caption:** not regenerated with R232 labels; the R232 matrix is `correlational/r232/C_within_r232.csv`.

## Not changed

- **Modeling text and numbers, and Figure 7:** the FBA modeling was not re-run (see comments 4 and 7).
- **Figure S11:** not re-rendered.
- **Abstract:** comment only.
- **Citation formatting:** citations in inserted text, e.g. (51), (53), (67), (68), are plain text, not reference-manager links. Relink them if you use one.
