# Within-condition co-occurrence: what the prefilters do (GTDB r232)

This note covers the within-condition pipeline: 4 conditions × 3 cores, followed by Graphical Lasso and Louvain.
- **Code:** `../prefilter_sensitivity_within.py`, a copy of `within_condition()` from `../rerun_correlational_r232.py` with the filters made into parameters. That file is neither imported nor modified.
- **Outputs:** the `within_*` files in this folder.
- **Baseline check:** the baseline reproduces the committed r232 files exactly (`within_baseline_verification.json`):
  - the 169 taxa;
  - all 11,573 tested pairs (max |Δ| of r_within, perm_p and q ≤ 1.1e-16);
  - the 326 Graphical Lasso edges (max |Δpartial r| 1e-16);
  - the Louvain modules (82 nodes, 186 edges, 6 modules, Q = 0.747).
- **Not approximated:** nothing was subsampled. Every setting uses the full 1,296-permutation null and 500 maxT shuffles.

## Short answer

1. **In this pipeline, the prevalence prefilter (MIN_COOCC = 6) removes almost none of the key organisms.**
   - All of these are among the 169 taxa:
     - every *phnJ* MePn-consumer taxon, including the archaeal *Natrinema* and *Haloquadratum*;
     - every *grdH* betaine-reducer taxon and all four TMA-route candidates;
     - *Methanolobus*, *Roseovarius*, *Albimonas*, *Marivita*, *Halomonas*, *Vreelandella* and *Halanaerobium*.
   - Lost at the taxon level:
     - *Thiolinea* (2020 MePn set only, present in 5 cores);
     - five archaeal taxa: *Nitrosopumilus*, B1-Br10-U2g21, *Nanosalina*, *Nanohalovita* and *Kariarchaeum*.
   - Lost at the MAG level only: BM714.2, a *grdH* carrier. Its taxon stays in through BM714.1/.3/.4, and BM714.4 also carries *grdH*.
   - This differs from the 4-site heatmap. Its "≥ 3 of 4 sites" filter drops *Methanolobus*, 4 of the 5 betaine-reducer taxa and 2 of the 10 MePn taxa (`../README.md`).
2. **Many key MAGs pass only just.** 69 MAGs are detected in exactly 6 cores: 41 only in the unrestored R1 + R2, and 25 only in R2A + SF2.
   - These include *Methanolobus*.1, *Halanaerobium*.1–3, BM714.1/.3/.4, CAZFJG01.1/.2, DATFIM01.1, *Haloquadratum*.2–4, *Phaeodactylibacter*.1 and *Vreelandella*.1 (unrestored only), plus UBA11579.1 (R2A + SF2 only).
   - A threshold of 7 would remove all of these MAGs, and with them the *Methanolobus*, *Halanaerobium*, BM714, CAZFJG01, DATFIM01, *Haloquadratum*, *Phaeodactylibacter*, *Vreelandella* and UBA11579 taxa.
   - For 66 of the 69, the value 6 means "present in all cores of one restoration block and nowhere else".
3. **Key taxa are missing from the network for downstream reasons, not because of the prefilter.**
   - **BH-FDR cannot pass any pair at any MIN_COOCC.** The smallest possible permutation p is 2/1297. BH therefore needs at least m/32.4 of the m tested pairs tied at that floor:
     - 357 needed at baseline, 26 observed;
     - 717 needed with no filter, 63 observed.
   - **Graphical Lasso degree tracks each taxon's residual variance (Spearman +0.90 at baseline).** The fit runs on unstandardised residuals, so taxa with low variance are shrunk out first.
     - *Roseovarius* sits at the 8th percentile of residual variance, *Marivita* the 14th, *Halanaerobium* the 19th, *Albimonas* the 22nd, *Methanolobus* the 34th and *Halomonas* the 35th. All have 0 edges.
     - The hub, DATFJG01 (36 edges), has the largest variance.
     - On standardised residuals, every key taxon gets 3–11 edges with |partial r| > 0.05.
4. **Removing the prefilters gives no guild organism (MePn, betaine, TMA, methanogen, *Halomonas*/*Vreelandella*) more edges.** DATFJG01 (23 → 9–11 edges above 0.05), M55B127 (2 → 1) and UBA11579 (3 → 1) lose most of their baseline partners and pick up a few new ones (DATFJG01 4 at MIN 4 and 6 at MIN 0; M55B127 and UBA11579 1 each), all of them newly admitted rare taxa. The taxa that gain edges are almost all newly admitted rare taxa, mostly bacteria (JAKAGZ01, CAZFAK01, DATFNI01), plus the nanohaloarchaea once MIN ≤ 3 admits them. The only exceptions are single-edge gains by two non-guild baseline taxa (M0040, through JAKAGZ01; DATFMH01 at MIN 4). Run with MIN_COOCC = 0, 1, 2, 3 and 4:
   - **Edges fall instead of rising:** 326 → 161–190, because the cross-validated α jumps from 0.51 to 0.76.
   - **The network is rewired:** only 26–30 of the 186 baseline edges above 0.05 survive (Jaccard 0.09–0.11), and module agreement on shared nodes is ARI 0.44–0.48.
   - **New edges go to newly admitted rare taxa:** 73–87 of the 105–116 edges touch taxa found in fewer than 6 cores. By degree above 0.05 these are chiefly the bacteria JAKAGZ01 (40–42), CAZFAK01 (18–20) and DATFNI01 (10–11), then the nanohaloarchaea *Nanosalina* (8–10, MIN ≤ 3) and *Nanohalovita* (5, MIN ≤ 2). Of the 18 taxa that gain edges at MIN 0, 15 are bacteria and 3 archaea.
   - *Haloquadratum* and the haloarchaeal module lose edges.
5. **The baseline network is not distinguishable from noise in size or modularity.**
   - Graphical Lasso was refitted on 20 within-condition-shuffled residual matrices (coupling destroyed, design kept).
   - **Edges:** the null fits give 167–534 (median 264), against 326 observed; 4 of the 20 have at least as many.
   - **Modularity:** Louvain on the null fits gives Q 0.60–0.80 (median 0.773), against 0.747 observed; 16 of the 20 are at least as modular.
   - **What is above the null is which taxa the edges join:**
     - 26.9 % of the edges above 0.05 join taxa of the same phylum, against 9.4–14.9 % in the 20 shuffles (further independent shuffles reach about 16 %, still well below);
     - there are 18 archaea–archaea edges (1–8 in the shuffles);
     - 73 % of edges are positive, against 54–69 % in the shuffles. This margin is thin (the largest of 20 shuffles is 4 points below), so the sign excess is suggestive only.
   - **The same-phylum excess is not mostly haloarchaeal.** Of the 50 same-phylum edges, 18 are Halobacteriota–Halobacteriota. The other 32 lie within bacterial phyla: Bacteroidota 15, Desulfobacterota 8, Pseudomonadota 5, Actinomycetota 3 and Bacillota 1. Among bacteria–bacteria edges alone, 19.8 % (32 of 162) are same-phylum, which is also above all of the further shuffles checked (see Caveats).
   - The same pattern holds with no prefilter (MIN 0; there, 16 of the 25 same-phylum edges are Bacteroidota and 5 are Halobacteriota) and with standardised residuals (table below).
6. **The other two switches:**
   - **EFFECT gate (|r_within| ≥ 0.7):** only a label. It does not feed Graphical Lasso or Louvain. 1,589 pairs pass it, but 784 are expected under the null (mean of the maxT shuffles; 95th percentile 823). With 0 FDR passes, turning it off changes nothing.
   - **2020 "defined IDs only" restriction (drops all 30 archaea):** at MIN_COOCC = 6 the network shrinks from 326 to 70 edges (α 0.51 → 0.82). Louvain is left with 25 nodes in 4 modules (Q 0.533). Adding or removing 13 archaeal taxa changes the whole network, not just the archaeal part.

**How essential is MIN_COOCC?**
- **It does not create power.** FDR is resolution-limited at every threshold, and the maxT 95th percentile is 0.992–0.995 everywhere.
- **What it does is limit garbage-in.** At baseline, sub-detection entries carry 10.8 % of the residual sum of squares; without the filter they carry 20.7 %. The rare taxa it admits are the most variable ones, and they set the Graphical Lasso penalty scale.
- **"Without prefilters" therefore does not recover the key organisms.** It swaps one noise-level network for another, centred on taxa seen in 2–5 cores.
- **The binding constraint is the design:** 3 replicates per condition leave 8 residual degrees of freedom for 169–216 taxa.

## Step by step (baseline, MIN_COOCC = 6)

`rerun:` = `../rerun_correlational_r232.py`, `marimo:` = `../../analysis_2_marimo.py` (2020), `louvain:` = `../../louvain_cooccurrence_network.py`.

| # | Step | Rule (file:line) | What survives | Key organisms affected |
|---|---|---|---|---|
| 1 | DNA-weighted depth merge; unbinned reads dropped | rerun:86-89, 253; marimo:1970-1982, 1997-2001 | 310 MAGs × 12 cores | Unbinned reads (44–93 % of each core) are not a CLR part, so the CLR runs over binned taxa only |
| 2 | Detection level ZERO = 1e-5 | rerun:49; marimo:1959 | Defines "present" only | 37.0 % of the 310 × 12 MAG entries are in (0, 1e-5]; only 0.03 % are exactly 0. These sub-detection values stay in the CLR as values (e.g. *Methanolobus* CLR −4.6 to −1.1 in R2A/SF2, where it is "absent") |
| 3 | MAG prevalence ≥ MIN_COOCC cores | rerun:386; marimo:1996 | 245 / 310 MAGs (23 / 30 archaeal); 169 / 216 taxa (13 / 18 archaeal) | Lost: *Thiolinea*.1 (2020 MePn), BM714.2 (*grdH*), *Haloquadratum*.1, and 5 archaeal taxa. 16 taxa (25 MAGs) are never above 1e-5 in any core. 69 MAGs sit exactly at 6 (point 2 above) |
| 4 | 2020 defined-ID restriction | marimo:1997-2001 (the r232 port keeps all 310, rerun:253) | If on: 222 MAGs, 156 taxa, 0 archaea | All 30 archaea, including *Methanolobus*, *Natrinema* and *Haloquadratum* |
| 5 | MAGs summed per taxon | rerun:254-255, 387; marimo:2002-2003 | 245 MAGs → 169 taxa | *Roseovarius* = 12 MAGs (11 *phnJ*); *Marivita* = 3 (2 *phnJ*); *Haloquadratum* = .2/.3/.4 (only .4 has *phnJ*). Summed sub-detection values cross 1e-5: the *Haloquadratum* taxon is "detected" in 9 cores and *Halanaerobium* in 7, though no single MAG is in more than 6 |
| 6 | Pseudocount + CLR | rerun:389-390; marimo:2013-2016 | Pseudocount 3.705e-8 (half the smallest value of the filtered matrix) | There are no exact zeros left at taxon level, so the pseudocount acts as a soft floor under the sub-detection values rather than replacing zeros |
| 7 | Within-site residuals | rerun:394-397; marimo:2020-2023 | 12 rows − 4 means → 8 residual degrees of freedom; rank ≤ 8 for 169 taxa | For unrestored-only taxa, half of the residual vector is sub-detection noise. *Methanolobus* residuals reach 0.67 in the R2A/SF2 cores, where it is "absent", against 0.89 in R1/R2 |
| 8 | Pair co-occurrence ≥ MIN_COOCC | rerun:435; marimo:2058 | 11,573 / 14,196 pairs | *Methanolobus* is tested against 123 of 168 partners, *Halanaerobium* 130, UBA11579 96, *Roseovarius* 168 |
| 9 | Exact 1,296-permutation null + BH-FDR 0.05 over all tested pairs | rerun:403-448; marimo:2041-2077 | 0 pairs (min q = 0.40) | All. Key-taxon minimum p is 2/1297–0.015 (`within_key_taxa_by_setting.csv`) |
| 10 | EFFECT gate \|r_within\| ≥ 0.7 | rerun:449; marimo:2082 | 1,589 pairs labelled (null mean 784) | Label only; nothing downstream reads it |
| 11 | maxT (500 within-condition shuffles, seed 42) | rerun:453-462; marimo:2095-2108 | 95th percentile 0.9930 | 1 pair is above it: DAHVFK01–Multivoradales, r = +0.9936 (its own exact permutation p = 4/1297). The runner-up, *Humimonas*–JAFDBW01 (r = −0.99292), is only 0.00003 below, so "one pair" reflects the 500 draws at seed 42. With 10,000 draws the 95th percentile is about 0.9927 and both pairs pass, at maxT-adjusted p = 0.032 and 0.047 (CAZWDL01–DT-110 next at 0.072). Neither pair involves a key organism |
| 12 | GraphicalLassoCV(alphas=10, max_iter=200) on unstandardised residuals | rerun:463-469; marimo:2133-2142 | α = 0.506; 326 edges | 0 edges: *Methanolobus*, *Roseovarius*, *Albimonas*, *Marivita*, *Halomonas*, *Vreelandella*, *Halanaerobium*, BM714, CAZFJG01, DATFIM01, T1Sed10-7, *Desulfotignum*, JAVJPJ01, *Roseibacterium*, CAZRZN01, DAINOF01, *Gaopeijia*. *Natrinema*: 1 edge (Halohasta +0.003). *Saliniramus*: 1 edge (M55B111 +0.045) |
| 13 | \|partial r\| > 0.05, then Louvain on positive edges (seed 42) | louvain:42, 95-125 | 186 edges, 82 nodes, 6 modules + 8 singletons, Q = 0.747 | See the note below |

**Who is in the step-13 network:**
- DATFJG01 (*grdH*; module_2; 23 edges above 0.05; partners such as *Christiangramia* +0.28 and *Gracilimonas* +0.28);
- UBA8670 (*grdH*; module_2, through DATFJG01 +0.056);
- M55B127 (*grdH*; module_3; CAZUKV01 +0.22);
- *Phaeodactylibacter* (*phnJ*; module_5; DT-110 +0.20);
- J153 (TMA candidate; module_1; JADGBM01 +0.06);
- *Haloquadratum* (*phnJ* via .4; module_4 with 8 other haloarchaea; *Halosegnis* +0.19, *Natronomonas* +0.15);
- UBA11579 (TMA candidate; a singleton, because all its edges are negative: *Humimonas* −0.12).

Per-MAG prevalence for every guild MAG, by site and with the threshold each passes, is in `within_key_MAG_prevalence.csv`.

## The sensitivity grid

Columns:
- **Pseudo:** the CLR pseudocount.
- **≤ ZERO:** the fraction of matrix entries at or below 1e-5.
- **SS ≤ ZERO:** the share of the residual sum of squares carried by those entries.
- **Floor obs / need:** pairs observed at the p floor against the number BH needs.
- **Effect obs (null):** pairs passing |r_within| ≥ 0.7, against the mean from the null shuffles.
- **GL:** Graphical Lasso.
- **J(base):** Jaccard overlap of edges with |partial r| > 0.05 against the baseline.
- **ARI:** Louvain agreement with the baseline on shared nodes.

There are no exact zeros at the taxon level in any setting. BH passes 0 pairs, and the gate-on / gate-off pass counts are 0, in all 14 settings; the same holds with the pair filter kept at 6. Full columns are in `within_settings_summary.csv`.

| Setting | MAGs | Taxa (arch.) | Pseudo | ≤ ZERO | SS ≤ ZERO | Pairs tested | Floor obs / need | Effect obs (null) | maxT 95 | GL α | GL edges (> 0.05) | Louvain nodes / modules / Q | J(base) | ARI |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **base, MIN 6** | 245 | 169 (13) | 3.7e-8 | 0.240 | 0.108 | 11,573 | 26 / 357 | 1,589 (784) | 0.993 | 0.506 | 326 (186) | 82 / 6 / 0.747 | 1 | 1 |
| MIN 4 | 269 | 186 (15) | 1.4e-9 | 0.273 | 0.165 | 14,533 | 35 / 449 | 2,133 (1,059) | 0.994 | 0.760 | 161 (105) | 68 / 6 / 0.649 | 0.111 | 0.48 |
| MIN 3 | 278 | 193 (16) | 1.4e-9 | 0.289 | 0.169 | 16,096 | 44 / 497 | 2,462 (1,239) | 0.993 | 0.762 | 175 (110) | 71 / 6 / 0.658 | 0.096 | 0.47 |
| MIN 2 | 280 | 195 (17) | 1.4e-9 | 0.294 | 0.176 | 16,922 | 45 / 522 | 2,532 (1,304) | 0.994 | 0.762 | 188 (116) | 73 / 6 / 0.683 | 0.098 | 0.47 |
| MIN 1 ("no prevalence filter") | 285 | 200 (18) | 1.4e-9 | 0.310 | 0.180 | 18,179 | 51 / 561 | 2,659 (1,390) | 0.993 | 0.763 | 190 (115) | 73 / 7 / 0.684 | 0.095 | 0.44 |
| MIN 0 (all 310 MAGs; also pairs that never co-occur) | 310 | 216 (18) | 1.3e-9 | 0.360 | 0.207 | 23,220 | 63 / 717 | 2,966 (1,573) | 0.994 | 0.764 | 183 (115) | 74 / 6 / 0.682 | 0.095 | 0.47 |
| MIN 6, 2020 IDs only | 222 | 156 (0) | 3.7e-8 | 0.243 | 0.108 | 9,809 | 22 / 303 | 1,447 (734) | 0.993 | 0.821 | 70 (36) | 25 / 4 / 0.533 | 0.150 | 0.77 |
| MIN 4, 2020 IDs | 243 | 171 (0) | 1.4e-9 | 0.275 | 0.169 | 12,227 | 29 / 378 | 1,989 (1,002) | 0.992 | 0.761 | 143 (86) | 58 / 5 / 0.634 | 0.079 | 0.44 |
| MIN 3, 2020 IDs | 250 | 177 (0) | 1.4e-9 | 0.290 | 0.165 | 13,455 | 31 / 415 | 2,295 (1,175) | 0.992 | 0.763 | 144 (93) | 63 / 5 / 0.624 | 0.073 | 0.42 |
| MIN 2, 2020 IDs | 251 | 178 (0) | 1.4e-9 | 0.293 | 0.168 | 14,044 | 33 / 434 | 2,364 (1,237) | 0.992 | 0.763 | 151 (94) | 63 / 5 / 0.629 | 0.073 | 0.42 |
| MIN 1, 2020 IDs | 255 | 182 (0) | 1.4e-9 | 0.307 | 0.171 | 15,044 | 30 / 464 | 2,464 (1,322) | 0.994 | 0.763 | 149 (97) | 66 / 5 / 0.628 | 0.072 | 0.41 |
| MIN 0, 2020 IDs | 280 | 198 (0) | 1.3e-9 | 0.362 | 0.200 | 19,503 | 55 / 602 | 2,732 (1,487) | 0.993 | 0.765 | 149 (100) | 66 / 5 / 0.621 | 0.075 | 0.41 |
| MIN 1, pseudocount pinned at 3.7e-8 | 285 | 200 (18) | 3.7e-8 | 0.310 | 0.160 | 18,179 | 48 / 561 | 2,654 (1,396) | 0.993 | 0.758 | 171 (111) | 71 / 6 / 0.677 | 0.112 | 0.48 |
| MIN 0, pseudocount pinned at 3.7e-8 | 310 | 216 (18) | 3.7e-8 | 0.360 | 0.180 | 23,220 | 56 / 717 | 2,952 (1,564) | 0.995 | 0.760 | 169 (109) | 71 / 6 / 0.672 | 0.105 | 0.43 |

**Does any key taxon gain Graphical Lasso edges when the filter is relaxed?** No guild taxon gains in number. DATFJG01, M55B127 and UBA11579 trade most of their baseline partners for a few newly admitted rare taxa. The taxa that gain are the newly admitted ones in the last two rows, mostly the bacteria JAKAGZ01, CAZFAK01 and DATFNI01 (`within_key_taxa_by_setting.csv` and `within_glasso_edges_by_setting.csv`; edges with |partial r| > 0.05):

| Taxon (role) | base 6 | 4 | 3 | 2 | 1 | 0 | 6, 2020 IDs | 0, 2020 IDs |
|---|---|---|---|---|---|---|---|---|
| *Methanolobus* (*mcrA*) | 0 | 0 | 0 | 0 | 0 | 0 | not in matrix | not in matrix |
| *Natrinema* (*phnJ*, archaeal) | 0 | 0 | 0 | 0 | 0 | 0 | not in matrix | not in matrix |
| *Haloquadratum* (*phnJ* in .4, archaeal) | 2 | 1 | 0 | 0 | 0 | 0 | not in matrix | not in matrix |
| *Roseovarius*, *Albimonas*, *Marivita*, *Desulfotignum*, JAVJPJ01, *Roseibacterium*, *Saliniramus* (*phnJ*) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| *Phaeodactylibacter* (*phnJ*) | 4 | 3 | 3 | 4 | 4 | 4 | 0 | 2 |
| *Thiolinea* (2020 MePn) | not in matrix | 0 | 0 | 0 | 0 | 0 | not in matrix | 0 |
| *Halanaerobium*, BM714 (*grdH*) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| DATFJG01 (*grdH*) | 23 | 9 | 9 | 10 | 10 | 11 | 12 | 11 |
| M55B127 (*grdH*) | 2 | 1 | 1 | 1 | 1 | 1 | 0 | 1 |
| UBA8670 (*grdH*) | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| J153 (TMA candidate) | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| UBA11579 (TMA candidate) | 3 | 1 | 1 | 1 | 1 | 1 | 0 | 1 |
| CAZFJG01, T1Sed10-7 (TMA candidates) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| *Halomonas*, *Vreelandella* | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Haloarchaeal module (*Halonotius*, *Halosegnis*, *Natronomonas*, *Halohasta*, *Haloarcula*, *Halorubrum*, *Haloplanus*) | 4–7 each | 0–5 | 0–5 | 0–6 | 0–5 | 0–6 | not in matrix | not in matrix |
| Newly admitted bacteria: JAKAGZ01 (4 cores) / CAZFAK01 (5) / DATFNI01 (4) | not in matrix | 40 / 20 / 10 | 41 / 19 / 11 | 41 / 18 / 11 | 41 / 18 / 11 | 42 / 18 / 11 | not in matrix | 44 / 17 / 11 |
| Newly admitted archaea: *Nanosalina* (3 cores) / *Nanohalovita* (2) / B1-Br10-U2g21 (5) | not in matrix | – / – / 1 | 8 / – / 2 | 10 / 5 / 1 | 10 / 5 / 1 | 9 / 5 / 1 | not in matrix | not in matrix |

**Partners of the newly admitted taxa, at MIN 0:**
- **Nanohaloarchaea:** the positive partners above 0.05 are haloarchaea (*Natronomonas*, *Halonotius*, *Halorubrum*, *Halosegnis*), the bacterium CAZFAK01, and each other. *Nanosalina* also has JAKAGZ01 (+0.06, a bacterium) and B1-Br10-U2g21 (+0.09, Nanohalarchaeota rather than a haloarchaeon). For example, *Nanosalina*–*Nanohalovita* is +0.21.
- **Caution:** both are present in only 2–3 cores (R2), so in the other 9–10 cores their residuals are log-ratios of sub-detection values.
- **The pinned-pseudocount runs show part of this is CLR scaling:** with the baseline pseudocount, *Nanosalina* has 6 Graphical Lasso edges (any size) at MIN 0 instead of 13, and the *Nanosalina*–*Nanohalovita* edge falls to +0.08. So part of these edges comes from the smaller pseudocount (1.3e-9 against 3.7e-8), not from more data.

## Why the statistics cannot be rescued by relaxing (or keeping) MIN_COOCC

- **Permutation resolution and BH.**
  - With 3 replicates per condition there are only 6⁴ = 1,296 within-condition permutations. The identity permutation always counts, so the smallest p is 2/1297 = 1.5e-3.
  - BH rejects the k smallest p only if p_(k) ≤ 0.05·k/m. That needs at least m·(2/1297)/0.05 = m/32.4 pairs tied at the floor (357 at baseline, 717 at MIN 0).
  - The observed counts at the first five p steps are 26 / 23 / 29 / 25 / 25 at baseline, against m/1296 = 8.9 expected per step under a uniform null.
  - Relaxing the filter raises m roughly as V² and the requirement grows with it, so the gap does not close. At MIN 0: 63 observed against 717 needed.
  - The roughly 3× excess at small p (and the about 2× excess of |r| ≥ 0.7 over the shuffle null) shows there is within-condition co-variation, or shared structure such as the common CLR denominator. No single pair can be resolved with this design. (Inference: the data run here cannot separate those two sources.)
- **maxT.**
  - With 8 residual degrees of freedom, the null maximum |r| over 12,000–23,000 pairs is 0.992–0.995 in every setting.
  - At 500 draws, only DAHVFK01–Multivoradales clears it, and only in the 6-core and some 2020-ID settings (`within_maxT_exceedances.csv`).
  - That count sits at the resolution of the 500 draws. At baseline, 10,000 draws put the 95th percentile at about 0.9927, and *Humimonas*–JAFDBW01 also passes, just (maxT-adjusted p = 0.047, against 0.032 for DAHVFK01–Multivoradales). Neither pair involves a guild taxon. At MIN 0 and MIN 4, 5,000 draws still pass no pair (smallest adjusted p ≈ 0.09 and ≈ 0.06).
- **Graphical Lasso with p ≫ n.**
  - R has rank ≤ 8 for 156–216 taxa. The default 5-fold cross-validation is unshuffled, so in two of the five folds a whole condition (R1 or R2) is held out.
  - The fit uses the raw covariance of the residuals, so the L1 penalty is not scale-invariant:
    - degree tracks residual variance (Spearman +0.90 at baseline, +0.80–0.82 in the relaxed settings; `within_glasso_degree_vs_variance_summary.csv`);
    - the median variance is 0.20 for taxa with 0 edges and 0.63 for taxa with ≥ 1 edge;
    - rare taxa admitted by relaxing are among the most variable (top residual variances at MIN 0: JAKAGZ01 2.72, CAZFAK01 2.59, *Nanosalina* 1.76);
    - they raise the α grid (maximum 1.41 → 2.11–2.13) and the chosen α (0.51 → 0.76), which shrinks the established edges away.
  - Summing several MAGs into one taxon (*Roseovarius*: 12) averages out variance, which works against the most-binned guilds. (Inference from the variance ranks above.)
- **Graphical Lasso against its own null** (`within_glasso_null_extended.csv`, `within_glasso_null_standardised.csv`):

| Fit | GL edges | > 0.05 | Louvain modules / Q | Same-phylum share (> 0.05) | Archaea–archaea edges | Positive share |
|---|---|---|---|---|---|---|
| Baseline (MIN 6), observed | 326 | 186 | 6 / 0.747 | 0.269 | 18 | 0.73 |
| Baseline, 20 shuffles: min – median – max | 167 – 264 – 534 | 108 – 165 – 308 | 6–11 / 0.600 – 0.773 – 0.800 | 0.094 – 0.120 – 0.149 | 1 – 4 – 8 | 0.54 – 0.60 – 0.69 |
| MIN 0 (no filter), observed | 183 | 115 | 6 / 0.682 | 0.217 | 14 | 0.74 |
| MIN 0, 20 shuffles | 124 – 206 – 330 | 91 – 130 – 197 | 7–11 / 0.648 – 0.714 – 0.776 | 0.081 – 0.116 – 0.147 | 0 – 3 – 6 | 0.47 – 0.56 – 0.64 |
| Baseline on standardised residuals, observed | 1,312 | 604 | 10 / 0.757 | 0.285 | 23 | 0.66 |
| Standardised, 10 shuffles | 1,038 – 1,057 – 1,085 | 567 – 585 – 608 | 10–14 / 0.747 – 0.764 – 0.786 | 0.154 – 0.168 – 0.186 | 2 – 5.5 – 8 | 0.51 – 0.53 – 0.56 |

How to read this:
- **Size and modularity are noise-level.** In the unstandardised fits, the edge count, the edges above 0.05 and Louvain Q all fall inside the shuffle null, both at baseline and with no filter. Neither the 326-edge network nor its 6 modules shows more structure than shuffled data.
- **Edge identity is not.** In every fit, the same-phylum share and the archaea–archaea count are above every shuffle in the table (1.5–1.8× and 2.3–2.9× the largest shuffle value). Those largest values are sample maxima from 10–20 draws, not bounds; further independent baseline shuffles exceed the same-phylum one but stay far below 0.269. The positive share is also above every shuffle in the table, but at baseline only by 4 points (0.73 against 0.69), so the sign excess is weak. So the residuals do carry co-variation among related taxa, among the haloarchaea (18 baseline edges) and within bacterial phyla such as Bacteroidota (15). But how many edges appear, and which ones, is set by the penalty (α and the variance scale), not resolved by the data.
- **The standardised fit** has more edges than its null (1,312 against at most 1,085), but the same number above 0.05 and the same Q. Its excess is again phylogenetic.

## What "without prefilters" would mean for the manuscript

- **Taking away MIN_COOCC alone does not bring the MePn consumers, betaine reducers or *Methanolobus* into the network.** It dilutes the network towards rare taxa seen in 2–5 cores, and in the r232-label settings it removes the one or two marginal pairs that clear maxT.
- **Taking away the 2020 "defined IDs" restriction** (keeping archaea, as the r232 port already does) is what produces the haloarchaeal module and the 326-edge network. That network rests on an α choice that moves with which taxa are included.
- **Defensible use of this analysis.** Neither the filtered nor the unfiltered within-condition network supports pairwise interaction claims for the guild organisms:
  - no pair passes FDR in any setting;
  - one or two marginal pairs pass maxT, depending on the number of shuffles, and neither involves a guild taxon;
  - the network's size and modularity match the shuffle null.
- **The only feature above the null is that edges join related taxa** (same phylum; archaea to archaea). That fits a shared response to within-site heterogeneity as well as it fits interaction. (Inference.)
- **If a network figure is kept,** standardising the residuals before Graphical Lasso removes the variance artefact. Every key taxon then gets edges, for example:
  - *Roseovarius*–*Marivita* +0.13 and *Roseovarius*–JAVJPJ01 +0.15, all Rhodobacteraceae *phnJ* carriers;
  - *Methanolobus*–*Chromatocurvus* +0.21 and *Methanolobus*–*Halomonas* +0.08.

  The standardised network has more weak edges than its shuffle null, but no more edges above 0.05 and no higher modularity. Individual edges cannot be tested with this design, so these partners are leads, not evidence of interaction.

## Caveats

- **Glasso null size:** 20 shuffles each for the baseline and MIN 0, and 10 for the standardised baseline. These are enough to place the observed values inside the null range, not to give precise tail probabilities. The quoted min–max ranges are those of these samples, not bounds, and more shuffles will widen them. Twelve further independent baseline shuffles (an ad hoc check with other seeds, not stored here) gave same-phylum shares up to about 16 %, above the 14.9 % quoted but far below the observed 26.9 %, and a bacteria-only same-phylum share of at most about 17 %. Their Q (0.70–0.79) stayed around the observed 0.747. Their edge counts (189–316) all fell just short of 326, so 4 of the 32 baseline shuffles in all reach the observed size, which is still not a significant excess.
- **maxT percentiles** come from 500 random draws, as in the 2020 code. At that resolution, the particular draw decides whether the second baseline pair passes. The 10,000-draw (baseline) and 5,000-draw (MIN 0, MIN 4) checks quoted above were run ad hoc with the script's shuffle loop and other seeds; they are not in the script outputs.
- **Graphical Lasso numerics:** GraphicalLassoCV emitted overflow / invalid-value RuntimeWarnings inside some CV folds (5–746 per fit), as in `../network_r232.log`. Every final precision matrix was finite, and there were no convergence warnings.
- **Not varied:** ZERO itself, the EFFECT value, the |partial r| > 0.05 cut, the CV folds, and including unbinned reads as a CLR part.
- **2020-ID settings** keep r232 labels on the 280 bacterial MAGs, so they are not the 140-taxon 2020 reproduction. The README documents that one separately.
- **Runtime:** about 15 min for the 14-setting grid on this machine, and about 1–2.5 min per null Graphical Lasso fit.

## Files (all `within_*` in this folder)

- `within_settings_summary.csv`: one row per setting (all columns above and more).
- `within_key_taxa_by_setting.csv`: per setting × key taxon. Columns:
  - MAGs included and per-MAG prevalence;
  - pairs tested, minimum p and q, effect passes;
  - all Graphical Lasso partners with sign and partial r;
  - degree above 0.05, Louvain module and co-members.
- `within_key_MAG_prevalence.csv`: every guild MAG (plus *Halomonas*/*Vreelandella*, *Roseovarius*, *Marivita*, *Albimonas* and *Halanaerobium*), with cores detected per site and the thresholds it passes.
- `within_taxa_membership.csv`: all 216 taxa, with roles, prevalence, and Graphical Lasso degree and module per setting.
- `within_glasso_edges_by_setting.csv`, `within_louvain_modules_by_setting.json`: the networks.
- `within_glasso_null_shuffles.csv` (3 shuffles, from the main run), `within_glasso_null_extended.csv` (20 shuffles) and `within_glasso_observed_coherence.csv`: the null comparison.
- `within_glasso_degree_vs_variance*.csv`, `within_glasso_standardised*.csv` and `within_glasso_null_standardised.csv`: the variance check and the standardised refit.
- `within_maxT_exceedances.csv`, `within_baseline_verification.json`.
- `within_run*.log`: console output.

Run, from `correlational/r232/`:
1. `~/Documents/py_venv/bin/python prefilter_sensitivity_within.py`
2. then `--null 20`
3. then `--variance`
4. then `--null-std 10`
