# Correlational studies rerun with GTDB r232 labels (2026-10-21)

Everything here is produced by `rerun_correlational_r232.py`. No 2020 output in `correlational/` was changed.
Labels come from `../iterativeIDs_r232.json` (GTDB-Tk 2.8.0, r232), which covers all 310 dereplicated MAGs. The 30 archaea
are therefore included; the 2020 analyses silently dropped them because they had no iterativeID.
Not rerun or re-rendered:
- **FBA simulations** (`fits*.json` / `saltern_fluxes*.json`): these are modelling, not correlation.
- **2020 renderings of the co-occurrence results:** `within_condition_glasso_network.png` (the source of manuscript Figure 8), `within_condition_triangle.png`, and the plotly/pyvis `cooccurrence_network_*.html` files. Their r232 input edges are in `glasso_direct_edges_r232.csv`.

## Guilds (`guilds_r232.csv`, one row per MAG, with the RAST evidence strings)

| Guild | Rule | MAGs | vs. 2020 set |
|---|---|---|---|
| MePn consumer | *phnJ* (C–P lyase, EC 4.7.1.1) | 22 (2 archaeal: *Natrinema*.1, *Haloquadratum*.4) | 2020 had 27; 5 had only accessory *phn* genes: CAZRZN01.1, DAINOF01.1, *Gaopeijia*.1, *Halorubrum*.3, *Thiolinea*.1 |
| Betaine reducer | *grdH* (betaine reductase B) | 9 (*Halanaerobium* ×3, BM714 ×2, M55B127 ×2, DATFJG01, UBA8670) | 2020 had 15 model-curated MAGs; 6 lack *grdH*: BM714.1, BM714.3, CAZFJG01.1, CAZFJG01.2, DATFIM01.1, UBA8670.1 |
| TMA-route candidate | *cntA*, or genes "clustered with choline TMA-lyase" | 4 (weak evidence; no MAG has CutC itself annotated) | — |
| Methanogen | *mcrA* | 1 (*Methanolobus*.1, the only Methanosarcinales) | — |
| Context | all archaea | 30 | — |

## Guild abundance vs CH₄ (`guild_methane_correlations_r232.csv`, `guild_methane_{raw,dnanorm}_r232.png`)

Spearman ρ on both bases:
- **Raw:** DNA-weighted depth merge versus areal flux.
- **DNA-normalised:** abundance/DNA versus CH₄/DNA, as in the submitted Figure 7.

"Null" is the share of random same-size MAG sets whose 12-core ρ is at least as high.

| Guild | Raw: 12 cores | Raw: 6 unrestored | Raw: null | DNA-norm: 12 cores | DNA-norm: 6 unrestored | DNA-norm: null |
|---|---|---|---|---|---|---|
| MePn consumers (*phnJ*, 22) | +0.76 (p = 0.005) | +0.49 (p = 0.33) | 32 % | +0.92 (p < 0.001) | +0.77 (p = 0.07) | **0.7 %** |
| MePn 2020 set (27) | +0.73 (p = 0.007) | +0.49 | 45 % | +0.87 | +0.60 | 2.5 % |
| Betaine reducers (*grdH*, 9) | +0.68 (p = 0.015) | −0.26 | 44 % | +0.77 | +0.20 | 20 % |
| Betaine 2020 set (15) | +0.70 (p = 0.011) | −0.09 | 46 % | +0.75 | +0.03 | 21 % |
| *Methanolobus*.1 (methanogen) | +0.78 (p = 0.003) | −0.54 (p = 0.27) | 24 % | +0.80 | −0.43 | 38 % |
| All archaea (30) | +0.77 (p = 0.003) | −0.31 | 32 % | +0.64 | +0.37 | 45 % |

How to read this:
- **On the raw basis, the 12-core correlations reflect the restoration gradient, not guild specificity.** No guild beats the random-subset null, and 254 of 310 individual MAGs have q < 0.05 (`mag_methane_correlations_r232.csv`).
- **On the DNA-normalised basis, the MePn consumers are the exception.**
  - Only 0.7 % of random 22-MAG sets match the *phnJ* set; for the 2020 set it is 2.5 %. Betaine, the methanogen and archaea stay at 20–45 %.
  - The null uses the same basis, so the shared 1/DNA divisor is applied to the random sets too.
  - Still, whether that basis is the right one is a judgement for the paper. Specificity appears only there, and the stricter *phnJ* definition sharpens it on both bases.
- **Within the six unrestored cores nothing is significant on either basis.** The best is MePn on the DNA-normalised basis, p = 0.07.
- ***Methanolobus*.1 is the only methanogen.** It is detected (relative abundance > 1e-5) only in the six unrestored cores, at up to 0.30 % of reads (`methanogen_archaea_by_core_r232.csv`). Within those cores the correlation with methane is negative but not significant (ρ = −0.54, p = 0.27, n = 6).
- **The DNA-normalised basis gives higher ρ** (e.g. MePn 0.92), because methane and abundance share the 1/DNA divisor.
- **Reference cores (R2A) have net CH₄ uptake** (negative flux), so the figures use a symmetric-log x-axis.

## Co-occurrence

- **`abundance_correlatons_one_triangle_reduced_defined_r232.png`:** the 2020 "all v all" heatmap, recomputed on 109 taxa across the 4 sites.
  - Its p < 0.05 boxes only mark perfect ±1 rank agreements. With n = 4, scipy's t-approximation assigns them p = 0, although they occur by chance about once in 12. The 2020 figure has the same caveat.
  - The filter keeps only taxa present at ≥ 3 of the 4 sites, which drops most guild taxa:
    - the methanogen *Methanolobus*;
    - 4 of the 5 betaine-reducer taxa (*Halanaerobium*, BM714, M55B127, UBA8670);
    - 2 of the 10 MePn taxa (*Haloquadratum*, *Phaeodactylibacter*).
- **Within-condition pipeline (4 conditions × 3 cores, 2020 parameters):** `all_candidate_pairs_r232.csv`, `direct_interactions_r232.csv` (empty), `glasso_direct_edges_r232.csv`, `C_*_r232.csv` and `marginal_vs_within_r232.png`.
  - The analysis covers 169 taxa: 13 archaeal and 156 bacterial, against 140 bacterial in 2020. Only part of the difference is the archaea. The rest is r232 resolution: MAGs that 2020 lumped at family or order level now split into genera, and 93 bacterial taxon names are new.
  - The port is verified: with the 2020 IDs, this code reproduces the 2020 `C_within.csv` exactly (140 taxa, max difference 1e-16).
  - 0 of 11,573 pairs pass FDR. The maxT 95th-percentile |r| is 0.993 and one pair just exceeds it (DAHVFK01–Multivoradales, r = 0.9936). Apart from that pair, within-condition coupling is not distinguishable from noise; 2020 was the same.
  - The permutation p-values are an exact vectorisation of the 2020 loop. They differ from it only in floating-point tie handling, at most one permutation step (1/1297).
- **`within_condition_louvain_network_r232.png` and `within_condition_louvain_modules_r232.json`:** Graphical-Lasso edges (326), drawn by `../louvain_cooccurrence_network.py`.
  - The network has 82 nodes and 186 edges (|partial r| > 0.05), split into 6 modules plus 8 singletons, Q = 0.747. 2020 had 5 + 12, Q = 0.701.
  - All 9 archaeal taxa in the network form one haloarchaeal module (module_4).
  - *Methanolobus*, *Roseovarius*, *Albimonas*, *Marivita*, *Halanaerobium* and *Halomonas*/*Vreelandella* have no edge, so they do not appear in the network.

Run from this folder: `~/Documents/py_venv/bin/python rerun_correlational_r232.py [--parts guilds,corr,heatmap,network]`

## Follow-up analyses (2026-10-26)

Each has its own script in this folder. None changes the outputs above.

- **Guild membership before vs after** (`guild_membership_2020_vs_r232.csv` / `.md`, `make_guild_membership_comparison.py`):
  - one row per MAG in the 2020 or R232 MePn and betaine sets, plus the methanogen and the TMA-route candidates;
  - for each: the change and its gene evidence, quality, abundance by site type, and whether it was modelled in 2020;
  - the reconciliation of the 27 curated MePn consumers with the 20 simulated.
- **Prefilter sensitivity** (`prefilter_sensitivity/abund_REPORT.md`, `prefilter_sensitivity/within_REPORT.md`; scripts `prefilter_sensitivity_abundance.py`, `prefilter_sensitivity_within.py`). For each filter: which key organisms it removes, what happens without it, and how essential it is given the downstream FDR, permutation and Graphical Lasso steps.
  - **4-site heatmap filter:** the ≥ 3-of-4-sites filter drops *Methanolobus* and 8 of 9 *grdH* MAGs.
  - **4-site heatmap statistics:** with n = 4, no pair can be significant under exact testing (smallest two-sided p = 1/12). The heatmap's p < 0.05 boxes are an artefact of scipy's t-approximation.
  - **Within-condition 6-of-12-cores filter:** keeps nearly every key taxon. Removing it gives no edge to a guild taxon that lacked one; the new edges mostly join rarely detected taxa.
  - **Why the key taxa have no edges:** Graphical Lasso on unstandardised residuals, where edge degree tracks residual variance (ρ = 0.90).
  - **Network vs a shuffle null:** size and modularity fall within the null; only the excess of edges between related taxa does not.
- **Archaeal MePn consumers and TMA-route candidates** (`archaeal_mepn_and_tma_candidates.md`, `_genes.csv`, `make_archaeal_mepn_tma_report.py`):
  - **Archaeal MePn consumers:** gene-level look at *Natrinema*.1 and *Haloquadratum*.4, which have complete C–P lyase loci on the annotations; *Halorubrum*.3 has none.
  - **TMA-route candidates:** the four candidates show no TMA-forming enzyme on the annotations.
- **Figure S11 regeneration** (`FigureS11_regeneration.md`): what the figure is, the code that drew it, what an R232 version needs, and a prototype renderer.
