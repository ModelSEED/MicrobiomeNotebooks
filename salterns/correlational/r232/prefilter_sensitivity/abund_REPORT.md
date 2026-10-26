# Abundance and prevalence prefilters: what they remove, and what happens without them

Scope: the filters that sit **outside** the within-condition pipeline (CLR residuals, permutation p, BH, Graphical Lasso). That pipeline's own permutation and FDR stage is covered by the within-condition sensitivity analysis. Here it is touched only where Figure S11's prevalence filter feeds it.

Everything below comes from `../prefilter_sensitivity_abundance.py`, which writes the `abund_*` files in this folder. Inferences are marked *(inference)*.

**Reproduction checks built into the script.** Each one is an assert that passes:
- The 2020 heatmap (2020 IDs, `defined`, `reduced`, at least 3 of 4 sites) is rebuilt with 94 taxa and 4,371 tested pairs. Of those, 1,006 have p = 0 and BH rejects 1,006, the same as the 2020 `pvals.txt`.
- The r232 heatmap has 109 taxa and 1,285 BH rejections (`heatmap_r232.log`).
- The S11 / within-condition matrix has 140 taxa with 2020 IDs (`C_within.csv`) and 169 with r232 IDs (`C_within_r232.csv`).
- Per-MAG BH on the raw basis gives 254 of 310.
- The guild ρ values match `guild_methane_correlations_r232.csv`.

---

## Short answers

1. **How the prefilters limit inclusion.** Two filters do almost all of the excluding.
   - **(a) The 2020 "defined iterativeID" filter** removes all 30 archaea. That includes the only methanogen (*Methanolobus*.1) and 3 MAGs of the 2020 MePn set itself (*Natrinema*.1, *Haloquadratum*.4, *Halorubrum*.3).
   - **(b) The 4-site heatmap filter** (site mean > 1e-5 at 3 or more of the 4 sites) removes, by construction, every MAG found only in the two unrestored salterns, where the methane is. There are 52 such MAGs. It also removes MAGs whose detection outside those salterns is marginal. Together this excludes *Methanolobus*.1, 8 of the 9 *grdH* betaine reducers (all three *Halanaerobium*, plus M55B127 and UBA8670, which reach 1.2–2.3 % of reads), and *Vreelandella*.1.
   - **(c) The S11 filter** (at least 6 of 12 cores) is much milder. With r232 labels it keeps all 22 *phnJ* MePn consumers, 8 of the 9 *grdH* reducers and *Methanolobus*.1. However, 41 unrestored-only MAGs pass at exactly the threshold (all 6 unrestored cores).
   - **(d) The guild-vs-CH₄ and per-MAG-vs-CH₄ correlations** have no prevalence filter at all.
   - *Roseovarius* (×12), *Albimonas*.1 and *Marivita* (×3) are detected in 12 of 12 cores and pass every filter.

2. **How essential the filters are.**
   - **(a) is pure bookkeeping.** It plays no statistical role.
   - **(b) cannot be rescued by any filter.**
     - With n = 4 sites, the smallest exact two-sided Spearman p is 1/12. So no pair can reach p < 0.05, filtered or not.
     - Every "significant" box in that heatmap is a |ρ| = 1 pair that scipy's t-approximation scores as p = 0.
   - **(c) needs some detection floor, but the "6 of 12" cut is a choice.**
     - CLR-plus-Pearson needs a floor: sub-detection cells carry 40–45 % of the CLR variance, and censored zeros produce identical CLR vectors that are perfectly correlated.
     - Downstream, the pair rule (co-occurrence ≥ 6) already enforces prevalence on every tested pair.
   - **(d) needs no filter.** Removing organisms costs far more than the BH power the filter buys. On these data, across families of 72 to 310 tests, the largest rejected p on the raw basis stays between 0.039 and 0.042.

3. **Effect of removing them.**
   - **Heatmap:**
     - |ρ| = 1 pairs rise from 1,285 to 3,930. In 2,680 of these, at least one taxon would have been filtered out.
     - Every one of those pairs would still be "significant" under scipy's p, and none would be under the exact p.
   - **S11 matrix:**
     - It grows from 169 to 216 taxa (13 to 18 archaeal), and the share of cells ≤ 1e-5 goes from 24 % to 36 %.
     - There are still no ties and no constant vectors, because the "zeros" are not zeros (see §3.2).
     - The filtered taxa's r_within values shift (mean |Δr| = 0.07; 517 pairs cross |r| = 0.7).
   - **Guild results:** they barely move, because the guild sums are dominated by their abundant members. What does change is the random-subset specificity null, which depends on which MAGs it may draw. When the null is drawn only from widespread MAGs, the 2020 MePn set no longer stands out (12.5 % of random sets match it), while the *phnJ* set still does (3.6 %) (§3.3).

---

## 1. The filters, where they are, and what they keep

Data: `data/Cliff_310MAG_relabund.tsv`, 310 MAGs plus `unbinned_reads`; unbinned reads are 42–94 % of each sample. Detection means relative abundance > 1e-5 of all reads. Abundances are DNA-weighted D1/D2 merges into 12 cores. Site values are means of 3 cores.

What ZERO = 1e-5 censors: 38.4 % of the 310 × 24 sample cells are ≤ 1e-5, yet only 22 cells are exactly 0. After the depth merge, 37.1 % of the 310 × 12 MAG × core cells are ≤ 1e-5 and only 1 is exactly 0. At the taxon level (216 taxa × 12 cores) the figure is 36.0 %, with no exact zeros.

| Filter | Rule | Code | MAGs kept / removed | Taxa kept | Archaeal MAGs kept | MePn *phnJ* (22) | MePn 2020 (27) | Betaine *grdH* (9) | Betaine 2020 (15) | *Methanolobus* | TMA candidates (4) | Binned abundance kept |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| (a) 2020 defined ID | drop MAGs with no 2020 iterativeID | `analysis_2_marimo.py:1348-1350` (also :1795-1796, :2000-2001) | 280 / 30 | 179 of 179 (2020 labels) | 0 / 30 | 20 | 24 | 9 | 15 | 0 | 4 | 93.9 % |
| (b) 2020 heatmap = (a) + ≥3/4 sites | site mean > 1e-5 at ≥ 3 sites, per MAG, before summing into taxa | `analysis_2_marimo.py:1344-1345` | 145 / 165 | 94 | 0 | 19 | 22 | 1 | 1 | 0 | 2 | 66.9 % |
| (b) r232 heatmap | same | `rerun_correlational_r232.py:283` | 164 / 146 | 109 of 216 | 19 | 20 | 24 | 1 | 1 | 0 | 2 | 72.5 % |
| (c) 2020 S11 = (a) + ≥6/12 cores | > 1e-5 in ≥ 6 of 12 cores, per MAG | `analysis_2_marimo.py:1996` (`_wc_run`, `MIN_COOCC=6`) | 222 / 88 | 140 | 0 | 20 | 23 | 8 | 14 | 0 | 4 | 91.7 % |
| (c) r232 within-condition | same | `rerun_correlational_r232.py:386` | 245 / 65 | 169 | 23 | 22 | 26 | 8 | 14 | 1 | 4 | 97.5 % |
| (c-old) ≥3/12 cores | older within-condition cell | `analysis_2_marimo.py:1790-1791`; `analysis(2).ipynb` cell 55 | 278 / 32 | 193 | 28 | 22 | 27 | 9 | 15 | 1 | 4 | 99.96 % |
| (d) guild sums, per-MAG vs CH₄ | none; per-MAG needs only a non-constant vector | `rerun_correlational_r232.py:164`, `:192` | 310 / 0 | 216 | 30 | 22 | 27 | 9 | 15 | 1 | 4 | 100 % |

Source: `abund_filter_summary.csv`; the per-MAG detail is in `abund_mag_prevalence.csv`. "Binned abundance kept" is the share of the summed MAG relative abundance, averaged over the 12 cores.

**What the removed MAGs look like:**
- **(b) removes 146 MAGs:**
  - 52 detected only in R1/R2 cores. This is every unrestored-only MAG in the catalogue; it is structural, because two sites can never make three.
  - 50 detected only in SF2/R2A cores.
  - 25 never above 1e-5 in any core.
  - 19 detected at both site types but above 1e-5 in the site mean at only 2 sites.
- **(c) removes 65 MAGs:**
  - 25 never detected;
  - 15 confined to one site;
  - 23 at two sites;
  - 2 at three sites.
  - 11 of the 65 are unrestored-only.
- **The 25 never-detected MAGs were all binned from water-column (H₂O) metagenomes.** They are *Aquiluna*, *Planktomarina*, *Actinomarina*, HIMB11 and similar. *(Inference: planktonic organisms that are not part of the sediment community.)* Removing them is uncontroversial.

**Partial-sum taxa** (`abund_partial_taxa.csv`). Both (b) and (c) are applied per MAG before MAGs are summed into a taxon. A taxon can therefore survive as a partial sum of its MAGs:
- (b) produces 8 such taxa; for example, *Psychroflexus* keeps 2 of its 4 MAGs and only 7 % of its abundance.
- (c) produces 5, including BM714, which loses its *grdH* carrier BM714.2, and *Haloquadratum*, which loses *Haloquadratum*.1.
- No guild-coloured taxon lost all of its gene-carrying MAGs, so the guild label colouring in the heatmaps is not misleading on this count.

## 2. Which key organisms each filter removes, and why

`abund_key_organism_prevalence.png` shows all 76 named MAGs, with detected cores per site and pass/fail for each filter. `abund_key_organisms.csv` holds the numbers.

**(a) 2020 defined ID** removes every archaeon:
- *Natrinema*.1 and *Haloquadratum*.4 (both *phnJ* carriers);
- *Halorubrum*.3 (2020 MePn set, accessory *phn* genes only);
- *Methanolobus*.1 (*mcrA*);
- 26 other haloarchaea and others.
- The 2020 heatmap and S11 therefore could never show 3 of the 27 MAGs of the 2020 MePn set.

**(b) ≥ 3 of 4 sites**, r232 labels:

| Organism | Role | Cores detected R1/R2/SF2/R2A | Max rel. abund. (12 cores) | Why it fails |
|---|---|---|---|---|
| *Haloquadratum*.4 | MePn *phnJ*, archaeon | 3/3/0/0 | 0.31 % | 2 sites |
| *Phaeodactylibacter*.1 | MePn *phnJ* | 3/3/0/0 | 0.35 % | 2 sites |
| *Thiolinea*.1 | MePn 2020 set only | 2/3/0/0 | 0.007 % | 2 sites |
| *Halanaerobium*.1, .2, .3 | betaine *grdH* | 3/3/0/0 each | 0.39 %, 0.35 %, 0.24 % | 2 sites |
| BM714.2, BM714.4 | betaine *grdH* | 3/2/0/0; 3/3/0/0 | 0.22 %, 0.35 % | 2 sites |
| M55B127.1, M55B127.2, UBA8670.2 | betaine *grdH* | 3/3/0/1 each | 1.8 %, 1.2 %, 2.3 % | the one R2A core leaves the R2A site mean ≤ 1e-5 |
| BM714.1, BM714.3, CAZFJG01.1, CAZFJG01.2, DATFIM01.1 | betaine 2020 set only (CAZFJG01.2 is also a TMA candidate) | 3/3/0/0 | 0.30–0.80 % | 2 sites |
| UBA8670.1 | betaine 2020 set only | 3/3/0/1 | 0.24 % | R2A site mean ≤ 1e-5 |
| *Methanolobus*.1 | methanogen *mcrA* | 3/3/0/0 | 0.30 % | 2 sites |
| UBA11579.1 | TMA candidate | 0/0/3/3 | 0.28 % | 2 sites (restored and reference only) |
| *Vreelandella*.1 (2020: *Halomonas*.2) | named | 3/3/0/0 | 0.04 % | 2 sites |
| 11 archaea | context | various | various | 2 or fewer sites |

The 11 archaea are *Haloquadratum*.1–.4, *Methanolobus*.1, B1-Br10-U2g21.1, *Nanosalina*.1, *Nanohalovita*.1, *Kariarchaeum*.1 and *Nitrosopumilus*.1, .2.

Kept by (b):
- *Natrinema*.1 (3 sites) and *Halorubrum*.3;
- DATFJG01.1, the only *grdH* MAG kept;
- J153.1 and T1Sed10-7.1 (TMA candidates);
- *Halomonas*.1;
- all *Roseovarius*, *Albimonas* and *Marivita* MAGs.

The betaine reducers that (b) removes are not rare. UBA8670.2 reaches 2.3 % and M55B127.1 1.8 % of reads in its highest core. **The filter selects for being widespread, not for being abundant.**

**(c) ≥ 6 of 12 cores** removes:
- *Thiolinea*.1 (5 cores; 2020 MePn set);
- BM714.2 (*grdH*; cores 3/2/0/0, so 5); the BM714 taxon survives through BM714.1, .3 and .4;
- 7 archaea: B1-Br10-U2g21.1 (5 cores), *Haloquadratum*.1 (3), *Nanosalina*.1 (3), *Nanohalovita*.1 (2), *Kariarchaeum*.1 (1), *Nitrosopumilus*.1 and .2 (5 each; SF2/R2A only).

It keeps, **at exactly the threshold**:
- *Methanolobus*.1;
- the three *Halanaerobium*;
- *Haloquadratum*.4 and *Phaeodactylibacter*.1;
- BM714.1, .3 and .4;
- CAZFJG01.1 and .2, DATFIM01.1;
- *Vreelandella*.1.

All of these are present in the 6 unrestored cores and nowhere else, so missing one core would remove them. UBA11579.1 passes at 6 with the SF2 and R2A cores.

With 2020 IDs, (c) is applied after (a), so the 2020 S11 also lacks every archaeon.

**(d)** removes nothing.

## 3. Without the prefilters

### 3.1 The 4-site Spearman heatmap (filter b)

With n = 4 there are only 24 rank orders. Enumerating them (`abund_spearman_n4_exact_null.csv`) gives:
- ρ takes exactly 11 values: 0, ±0.2, ±0.4, ±0.6, ±0.8 and ±1. The orders per value are 1, 3, 1, 4, 2, 2, 2, 4, 1, 3, 1 from −1 to +1.
- The exact two-sided p at |ρ| = 1 is **2/24 = 1/12 = 0.083**. At |ρ| = 0.8 it is 8/24 = 0.33.
- The one-sided p at ρ = +1 is 1/24 = 0.042.
- `scipy.stats.spearmanr` (scipy 1.17.1) returns **p = 0.0** at |ρ| = 1 and p = 0.2 at |ρ| = 0.8. These were verified by calling it directly.

Variants (`abund_heatmap_variants.csv`, `abund_heatmap_rho_values.csv`). The co-present pairs are the BH family the code uses: pairs present at a common site.

| Variant | Taxa (archaeal) | All pairs | Co-present pairs m | \|ρ\| = 1 pairs (m/12 expected if independent) | BH rejections, scipy p | BH rejections, exact p |
|---|---|---|---|---|---|---|
| 2020 heatmap reproduced | 94 (0) | 4,371 | 4,371 | 1,006 (364) | 1,006 | 0 |
| 2020 `reduced` without `defined` | 107 (13 pseudo-taxa) | 5,671 | 5,671 | 1,204 | 1,204 | 0 |
| r232 heatmap, as published in `r232/` | 109 (11) | 5,886 | 5,886 | 1,285 (490) | 1,285 | 0 |
| r232, **no prevalence filter** | 216 (18) | 23,220 | 17,383 | 3,930 (1,449) | 3,930 | 0 |
| r232, no filter, values ≤ 1e-5 set to 0 | 216 (18); 19 constant, 107 with ties | 19,306 finite | 17,162 | 2,025 | 2,025 | 0 |
| 2020 IDs (defined), no filter | 179 (0) | 15,931 | 11,884 | 2,758 | 2,758 | 0 |

**What this shows:**
- **The BH verdict does not depend on m or on the filter.**
  - With scipy p, every |ρ| = 1 pair has p = 0, which lies below any BH line. So the rejections equal the |ρ| = 1 count, in every variant.
  - With exact p, the smallest attainable p is 1/12 > 0.05, so BH rejects nothing at any m. Even at k = m the BH line is 0.05. At k = the number of |ρ| = 1 pairs it is 0.006–0.012.
  - A one-sided exact test would also fail: 1/24 = 0.042 against a BH line of 0.0064 at k = 748 positive ±1 pairs, m = 5,886.
- **There are more |ρ| = 1 pairs than independence predicts.** For example, 1,285 against 490. *(Inference: many taxa share the same unrestored-versus-restored/reference ordering of the four sites. A |ρ| = 1 pair therefore marks a shared response to the site gradient, not an interaction.)*
- **Removing the filter adds noise-level taxa to the boxes.**
  - 2,680 of the 3,930 |ρ| = 1 pairs involve a taxon the filter would have removed, and 618 involve two.
  - *Methanolobus* gets 25 "significant" partners out of 150 tested, *Halanaerobium* 39 of 150 and *Haloquadratum* 36 of 189.
  - Only 22 of the 7,440 sample cells are exactly 0, so the sub-1e-5 site means are almost all non-zero. *(Inference: these are mapping-noise-level values, so without the filter Spearman ranks noise.)*
  - Censoring at 1e-5 creates ties and 19 all-zero columns. ρ then takes 30 distinct values, such as ±0.949 and ±0.738, and 19 taxa have undefined ρ.
- **The 2020 `reduced` variant without `defined`** (behind the 2020 file `abundance_correlatons_one_triangle_reduced.png`) did include the 19 archaea that pass the filter. But it summed them under bin-name stems: 12 pseudo-taxa such as `..._R2_A_D2_MG_DASTool_bins_concoct_out`, which equals *Halonotius*.1 + *Halosegnis*.1 + *Halorubrum*.2. It also added an `unbinned` column.

### 3.2 The S11 / within-condition matrix (filter c)

Figure S11 is `within_condition_triangle.png`:
- It is drawn from `wc_C_within` (`analysis_2_marimo.py:2288`, saved at :2384).
- `_wc_run` builds that matrix; its prevalence filter is at :1996.
- The caption's "below 1e-5 in more than half of the 12 samples removed" equals keeping taxa present in at least 6 of 12 cores, which matches the code (140 taxa with 2020 IDs, the same as `C_within.csv`).
- The older cell (:1790-1791, ≥ 3 of 12, pseudocount +0.5) would give 159 taxa with 2020 IDs.
- The caption's blank asterisk threshold is, in the code, a Graphical-Lasso edge with |partial r| > 0.05. That is `PARTIAL_FOR_ASTERISKS` at :2338, with the asterisk colour `#90EE90` set at :2225.

`abund_s11_matrix_properties.csv`, r232 labels (2020 labels behave the same):

| | ≥ 6 of 12 (S11) | ≥ 3 of 12 (old) | No filter |
|---|---|---|---|
| Taxa (archaeal) | 169 (13) | 193 (16) | 216 (18) |
| Cells ≤ 1e-5 | 24.0 % | 28.9 % | 36.0 % |
| Exact zeros after the depth merge (taxon level) | 0 | 0 | 0 |
| Taxa with tied values / constant / NaN r_within (uncensored) | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |
| Pseudocount (half the minimum non-zero, as `_wc_run`) | 3.7e-8 | 1.4e-9 | 1.3e-9 |
| CLR range (ln units) | 13.5 | 15.3 | 16.9 |
| Share of ΣCLR² from cells ≤ 1e-5 | 39.9 % | 43.2 % | 45.2 % |
| Pairs: all / co-occurring in ≥ 6 cores | 14,196 / 11,573 | 18,528 / 11,860 | 23,220 / 11,860 |
| Censored at 1e-5: all-zero taxa → pairs with r_within = ±1 exactly | 0 → 0 | 0 → 0 | 15 → 105 (= C(15,2)) |

- **Without the filter, no "zero problem" appears unless one is created.**
  - The raw matrix has no exact zeros, ties or constant vectors, so CLR and Pearson stay defined.
  - The problem only appears if values below detection are treated as absent. Then every all-absent taxon gets the same CLR vector (log pseudocount minus the sample's log geometric mean). That yields 105 exact ±1 within-condition correlations that are pure artefact.
- **The prefilter changes the CLR reference and the pseudocount.**
  - Recomputing r_within for the 169 S11 taxa with CLR taken over all 216 taxa gives r = 0.978 between the two versions, with mean |Δr| = 0.073 and max |Δr| = 0.57.
  - 517 of 14,196 pairs cross the |r| = 0.7 effect threshold (1,779 pairs above it filtered, 1,852 unfiltered).
- **Sub-detection cells dominate the CLR geometry either way.** They carry about 40 % of ΣCLR² even with the filter.
- **The co-occurrence ≥ 6 pair rule makes the MAG-level prefilter almost redundant for the tests.**
  - The rule is at `rerun_correlational_r232.py:434-435` and `analysis_2_marimo.py:2058`.
  - Dropping the prefilter adds 9,024 pairs, but only 287 of them (+2.5 %) become testable.
  - *(Inference: the prefilter's independent influence is on the CLR reference, the pseudocount and the Graphical-Lasso input, which uses every taxon with no pair rule. The within-condition analysis quantifies that.)*
- **The older cell's +0.5 pseudocount shrinks the CLR range to 0.15.** Every value lies within 0.15 ln units, so CLR becomes nearly linear in abundance. `_wc_run`'s own comment says this.

### 3.3 Guild and per-MAG correlations with CH₄ (filter d)

**There is no prevalence filter in these correlations:**
- `guild_correlations` sums every member (`rerun_correlational_r232.py:164`).
- `ZERO` only feeds `cores_detected` (:169, :174) and the open "below detection" markers (:217).
- `per_mag_correlations` drops only constant vectors (:192), and none of the 310 MAGs is constant.
- The 2020 cells (`analysis_2_marimo.py:321-327`, per MAG on 24 DNA-normalised samples, no p-values; :350-357 guild sums) and `clauded_insights/reassessment_figure.py` also use every member.
- CH₄ is the same in the D1 and D2 rows of each core (a core-level flux), so merging CH₄ changes nothing.

**Robustness: what if members were filtered or censored?** From `abund_guild_sum_sensitivity.csv`, Spearman ρ over the 12 cores:

| Guild | As run (raw / DNA-norm.) | (a) members only | (b) members only | (c) members only | Censored ≤ 1e-5 → 0 |
|---|---|---|---|---|---|
| MePn *phnJ* | 0.755 / 0.916 | 0.727 / 0.881 (20) | 0.755 / 0.895 (20) | 0.755 / 0.916 (22) | 0.755 / 0.916 |
| Betaine *grdH* | 0.678 / 0.769 | same (9) | 0.650 / 0.804 (1 MAG = 9 % of guild abundance) | 0.678 / 0.776 (8) | 0.666 / 0.723 |
| *Methanolobus*.1 | 0.776 / 0.804 | absent | absent | 0.776 / 0.804 | 0.735 / 0.750 |
| All archaea | 0.769 / 0.636 | absent | 0.853 / 0.699 (19) | 0.832 / 0.699 (23) | 0.769 / 0.636 |

Guild sums are dominated by their abundant members. Restricting membership changes ρ by ≤ 0.04 for MePn, unless the filter removes the guild entirely (archaea and the methanogen under (a); betaine under (b)).

**The specificity null depends on the pool it draws from.** The table gives the same guild sums with random same-size sets drawn only from MAGs that pass a filter; values are the share of random sets reaching the guild's ρ (DNA-normalised basis):

| Guild | All 310 (as run) | 280 bacteria (a) | 164 MAGs (b) | 245 MAGs (c) |
|---|---|---|---|---|
| MePn *phnJ* | 0.9 % (`README`: 0.7 %) | 0.3 % | 3.6 % | 1.2 % |
| MePn 2020 set | 2.7 % (`README`: 2.5 %) | 0.8 % | 12.5 % | 4.2 % |
| Betaine *grdH* | 19 % (`README`: 20 %) | 15 % | 39 % | 22 % |

*(Inference: the MePn sets are mostly MAGs found in all 12 cores. A null drawn from equally widespread MAGs is the fairer comparison.)* Under it:
- the *phnJ* set stays below 5 % (3.6 %);
- the 2020 set does not (12.5 %).

**BH and the number of tests m**, per-MAG vs CH₄ over 12 cores (`abund_bh_vs_m.csv`):

| Family (prefilter) | m | Raw: discoveries (largest rejected p) | DNA-norm.: discoveries | *grdH* tested | Methanogen tested |
|---|---|---|---|---|---|
| All 310 (as run) | 310 | 254 (0.039) | 227 | 9 | 1 |
| (a) 2020-defined | 280 | 227 (0.039) | 201 | 9 | 0 |
| (b) ≥ 3/4 sites | 164 | 131 (0.039) | 119 | 1 | 0 |
| (c) ≥ 6/12 cores | 245 | 211 (0.042) | 193 | 8 | 1 |
| ≥ 3/12 cores | 278 | 237 (0.042) | 218 | 9 | 1 |
| Guild members only | 72 | 69 (0.039) | 67 | 9 | 1 |

- **Shrinking m barely moves BH here.** On the raw basis, the largest rejected p stays between 0.039 and 0.042, and the BH line at the last rejection between 0.040 and 0.048. The reason is that most p-values are tiny: the restoration gradient makes 262 of 310 MAGs nominally significant (raw basis).
- **The real effect of a filter is to leave organisms out of the test.**
- For contrast, a synthetic family of 5 tests at p = 1e-4 survives BH up to m = 1,000 and fails at m = 5,886, where the BH line at k = 5 is 4.2e-5. Shrinking m matters when the signal is weak relative to m, which is not the case here.

### 3.4 The DNA-weighted depth merge (d)

From `abund_depth_merge_sensitivity.csv`.

**How the merge weights the depths:**
- `merge_depths_by_DNA` (`analysis_2_marimo.py:65-66`; `rerun_correlational_r232.py:86-89`) weights D1 and D2 by DNA concentration (µg/g) only.
- The weight on D1 ranges from 0.054 (R2_B) to 0.955 (R2_C). So in R2_B the "core" is essentially D2, and in R2_C essentially D1.
- The same file defines a 5:10 thickness weighting (`_volume_weighted_sample_dna`, :53) but does not use it for abundances. P27 gives 15 cm cores and P59 gives D1 = 0–5 cm. *(D2 = 5–15 cm is an inference.)*
- On the DNA-normalised basis the merge reduces to (rel_D1 + rel_D2)/(DNA_D1 + DNA_D2).

**What changes under other merges:**
- 110 of 3,720 MAG × core cells have one depth above 1e-5 but a merged value ≤ 1e-5.
- MAGs passing (b) range from 158 to 176 across merges (plain mean, 5:10, DNA × 5:10, maximum of the two depths); MAGs passing (c) range from 243 to 255.
- Key organisms that flip:
  - Merging by the maximum of the two depths would let *Methanolobus*.1, UBA8670.2 and M55B127.2 pass (b).
  - B1-Br10-U2g21.1 drops out of (c) under the thickness merges.
- Guild ρ (raw, 12 cores) across merges:

  | Guild | Range of ρ |
  |---|---|
  | MePn *phnJ* | 0.734–0.769 |
  | Betaine *grdH* | 0.678–0.755 |
  | Betaine 2020 set | 0.664–0.762 |
  | TMA candidates | 0.13–0.25 |
  | *Methanolobus* | 0.776–0.797 |

  So the merge is a modest choice for the guild results: ρ moves by ≤ 0.08, or 0.11 for the 4 TMA candidates. At the 1e-5 boundary it matters more.

## 4. How essential is each filter?

| Filter | Statistical role given the downstream methods | Cost in coverage | Verdict |
|---|---|---|---|
| (a) defined 2020 ID | None. It is bookkeeping, because archaea had no 2020 label, and no test needs it. | All 30 archaea, including the methanogen and 3 of the 2020 MePn MAGs | **Arbitrary and harmful.** It is obsolete now that r232 labels cover all 310 MAGs. |
| (b) ≥ 3 of 4 sites (heatmap) | None that works. With n = 4 no pair can be significant under an exact test, so no filter makes the boxes valid. A detection floor only decides whether noise gets ranked. | All 52 unrestored-only MAGs; *Methanolobus*; 8 of 9 *grdH* reducers; 2 *phnJ* MAGs; 2 TMA candidates; 11 archaea | **Arbitrary.** It biases the figure towards cosmopolitan taxa and against the methane-site specialists. The significance boxes are a t-approximation artefact either way. |
| (c) ≥ 6 of 12 cores (S11) | Partly necessary. A detection floor is needed for CLR with a pseudocount: sub-detection cells carry about 40 % of the CLR variance, and censored zeros create identical vectors and exact ±1 r. The pair rule (co-occurrence ≥ 6) already enforces prevalence on tested pairs. | 1 *grdH* MAG (BM714.2), *Thiolinea*.1, 7 archaea. 41 unrestored-only MAGs sit exactly on the cut. | **The floor is justified; the cut-off of 6 is a choice.** Report it as a sensitivity parameter, along with how r_within shifts with the CLR reference. |
| (d) ZERO in guild and per-MAG correlations | None needed. Guild sums are robust (censoring changes ρ by ≤ 0.06). BH gains almost nothing from a smaller m on these data. | None as run | **No filter, which is correct.** The limiting factor is the between-site confound, not filtering. The random-subset null should be prevalence-matched. |
| DNA-weighted depth merge | A weighting choice. It gives up to 95 % weight to one depth. | 110 cells pushed below detection; 3 key MAGs on the (b) boundary | **A defensible but undocumented choice.** It shifts guild ρ by ≤ 0.08 (0.11 for the TMA candidates). |

## 5. Implications for the manuscript (suggestions)

- **Figure 8a caption.** It says nodes are organisms "with presence (≥ 1e-5) in ≥ 3 experimental sites". The network (`within_condition_glasso_network`, `analysis_2_marimo.py:2605-2856`) is built from `_wc_run`, whose filter is ≥ 6 of 12 cores (:1996). The caption describes the heatmap filter (b), not the one used. This is established from the code path only; the docx image was not checked here.
- **Figure S11 caption.**
  - "Pairs whose correlations exceed [blank]" should describe the asterisks as Graphical-Lasso edges with |partial r| > 0.05 (:2338).
  - The 2020 version also excluded all archaea through (a). This should be stated, or the figure regenerated with r232 labels (169 taxa, 13 archaeal).
- **The 4-site "all v all" heatmap.** If it is used anywhere, drop the p < 0.05 boxes, or state that they mark |ρ| = 1 at n = 4, whose exact p is 1/12. Also state that the ≥ 3 of 4 sites filter removes every organism confined to the unrestored salterns.
- **Methods text should state that guild sums include every member** (no prevalence filter), and that 1e-5 is only a display and detection label there.

## Files

| File | Content |
|---|---|
| `../prefilter_sensitivity_abundance.py` | the analysis; run it from `correlational/r232/` with `~/Documents/py_venv/bin/python prefilter_sensitivity_abundance.py` (about 30 s). It imports the inputs from `rerun_correlational_r232.py` and reads `guilds_r232.csv`, `mag_methane_correlations_r232.csv` and `guild_methane_correlations_r232.csv`. |
| `abund_mag_prevalence.csv` | all 310 MAGs: detected samples/cores/sites, site means, pass/fail per filter |
| `abund_key_organisms.csv` | the 76 named MAGs: guild members, *Methanolobus*, TMA candidates, *Roseovarius*/*Albimonas*/*Marivita*, *Halomonas*/*Vreelandella*, *Halanaerobium*, archaea |
| `abund_filter_summary.csv` | §1 table, plus counts of removed MAGs by site pattern |
| `abund_partial_taxa.csv` | taxa that survive (b)/(c) as partial sums; guild-colouring check |
| `abund_spearman_n4_exact_null.csv` | the n = 4 exact null against scipy's p |
| `abund_heatmap_variants.csv`, `abund_heatmap_rho_values.csv` | §3.1 |
| `abund_s11_matrix_properties.csv` | §3.2, including the CLR-reference comparison (`CLRref_*` columns) |
| `abund_bh_vs_m.csv` | §3.3, BH vs m (per-MAG families and the synthetic illustration) |
| `abund_guild_sum_sensitivity.csv` | §3.3, guild ρ by member filter, censoring and null pool (20,000 random sets per row, seed 0) |
| `abund_depth_merge_sensitivity.csv` | §3.4 |
| `abund_key_organism_prevalence.png`, `abund_n4_spearman_bh.png` | the two figures |

Caveats:
- The random-subset null fractions are Monte-Carlo estimates from 20,000 draws, giving a standard error of 0.35 points or less. They differ from `guild_methane_correlations_r232.csv` (4,000 draws) by up to about 2 percentage points. For example, 5 repeats of 20,000 draws for the *phnJ* set (DNA-normalised basis) gave 0.76–0.91 %.
- Detection at 1e-5 is a fraction of all reads, and binned MAGs cover only 6–58 % of reads per sample. *(Inference: the same threshold is therefore relatively stricter for MAGs at the restored and reference sites, where the binned fraction is lowest.)*
