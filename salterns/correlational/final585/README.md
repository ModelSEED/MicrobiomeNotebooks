# Correlational analyses on the final re-assembled MAG set (585 species)

On 2026-11-20 the analyses of `../r232/` were rerun on the final MAG set, with the same statistics and parameters. The community changes:
- **Old:** the 310 original dereplicated MAGs.
- **New:** the 585 species representatives of the re-assembled set. These come from 935 second-pass (Arm C) MAGs plus 520 kept MAGs of the first-pass site co-assemblies, clustered at skani 95 % ANI / 50 % AF (`../../reassembly_2026/comparison_rollout/`).

The manuscript review (`../../clauded_insights/make_r232_review_docx.py`) reads its Tables 3–6, Figures 8, 9, S8 and S11, and the abundance, guild and co-occurrence text from this folder.

## Inputs

| What | Source |
|---|---|
| Abundance | `../../reassembly_2026/abundance_final/final585_relabund_24sediment.tsv`. CoverM 0.7.0 genome mode, minimap2 2.30, ≥ 95 % identity, ≥ 75 % aligned, mapping the BBDuk-trimmed read pairs; read count / total reads; 24 sediment samples plus an `unbinned_reads` row. |
| Guilds | `../../reassembly_2026/comparison_rollout/final_guilds_species.tsv` (KOfam 2026-10-26, DRAM rule, Sec read-through). A species is in a guild when **any** of its MAGs carries the marker; the representative's own call is kept as `<guild>_rep`. |
| Taxonomy | GTDB-Tk 2.8.0 / R232 of the representatives (`mags_three_way.tsv`). |
| CH4, DNA, metadata | `../../data/` (unchanged). |

## Scripts (run in this order, from this folder, with `~/Documents/py_venv/bin/python`)

1. **`make_ids_final585.py`** writes `iterativeIDs_final585.{json,csv}`, `iterativeID_levels_final585.json` and `phylogeny_final585.json`.
   - Each label is the deepest named rank down to genus plus `.n`.
   - **Inherited labels:** a species containing an original MAG keeps that MAG's R232 label when the taxon matches (260 of 260). With several originals in one cluster, it keeps the lowest number: BM714.3 absorbs BM714.4, Roseovarius.10 absorbs .11, DSBN01.1 absorbs .2/.3/.6, Halanaerobium.2 absorbs .3, JAAXQD01.1 absorbs .2.
   - **New labels:** other species are numbered after the highest number already used in their taxon, by decreasing mean abundance. Lost originals' numbers are never reused.
2. **`rerun_correlational_final585.py [--parts guilds,corr,heatmap,network,table3]`** writes the `*_final585.*` outputs:
   - guild table;
   - guild–CH4 correlations with the 4,000-set random null;
   - per-species correlations (BH);
   - Figures 9 and S8 (`guild_methane_{raw,dnanorm}_final585.png`);
   - the 4-site heatmap;
   - the within-site permutation and Graphical Lasso pipeline, with the Louvain figure (`within_condition_louvain_network_final585.png`, via `../louvain_cooccurrence_network.py`);
   - `table3_final585.json`;
   - `summary_final585.json`.
3. **`figureS11_final585.py`** renders Figure S11 (`figureS11_final585.png`, `.order.json`). It is the R232 renderer with a canvas scaled to the taxon count.
4. **`network_sensitivity_final585.py --null 20 --workers 10`** (run on TF, `~/tf_runs/salterns_network_sens_final585/`, venv `~/venvs/salterns_correlational`) writes `network_sensitivity_final585/`:
   - a 20-shuffle Graphical Lasso null;
   - residual variance vs degree, and a standardised refit;
   - prefilter relaxation (MIN_COOCC 1 and 0).

   The baseline refit on TF reproduced the local network exactly (194 edges, same α).

## Results (2026-11-20)

**Guilds** (species / MAGs carrying the marker):

| Guild | Species | MAGs |
|---|---:|---:|
| MePn (*phnJ*) | 39 | 109 |
| Betaine (*grdH*/*grdI*) | 18 | 53 |
| Methanogen (*Methanolobus* sp020055655) | 1 | 4 |
| TMA candidates (*cutC* or *cntA* only) | 36 | 83 |
| Archaeal species | 57 | |

There are no complete TMA pairs.

**Guild vs CH4, 12 cores, areal:**

| Guild | ρ | p | Random sets ≥ ρ |
|---|---:|---:|---:|
| MePn | 0.75 | 0.005 | 27 % |
| Betaine | 0.64 | 0.024 | 46 % |
| Methanogen | 0.73 | 0.007 | 37 % |
| Archaea | 0.76 | 0.005 | 27 % |

- 480 of the 585 species are individually correlated (q < 0.05).
- **DNA-normalised:** MePn ρ = 0.78 (9 % of random sets ≥ ρ) and betaine ρ = 0.76 (14 %).
- **This differs from the 310-MAG analysis,** where MePn had ρ = 0.92 on the DNA-normalised basis, beating 99 % of random sets. That specificity claim does not hold on the final set.
- **Six unrestored cores:** only the DNA-normalised MePn sum is nominally correlated (ρ = 0.83, p = 0.042, n = 6; null 19 %).

**Co-occurrence:**
- 287 species present in ≥ 6 cores form 201 taxa (23 archaeal). No pair passes FDR, and one pair is ≥ the maxT 95th percentile (0.985).
- Glasso: 194 edges; 117 have |partial r| > 0.05, among 62 taxa.
- Louvain: 7 modules, Q = 0.741. Eight of the 10 archaeal taxa in the network are in one module (16 taxa, 8 of them bacterial).
- **Null (20 shuffles):** 236–617 edges, so the observed network is sparser than every shuffle. Q is 0.69–0.76, so modularity is not exceptional. Same-phylum edges are 21 % observed against at most 13 % in the null.
- **Guild taxa:** 8 of the 32 in the panel have edges. MePn: *Haloquadratum*, *Haloplanus*, *Halonotius*, *Halorubrum*. Betaine: DATFJG01, QPJE01, BM714, *Halonatronomonas*.
- *Roseovarius*, *Marivita*, *Methanolobus*, *Halanaerobium* and M55B127 have no edges. Residual variance tracks degree (ρ = 0.80), and the standardised refit gives these taxa 3–10 edges.
- Relaxing the prefilter adds no edge to an in-panel guild taxon. Of the 284 new edges, 282 involve below-filter taxa.

**4-site heatmap:** 58 taxa against 109 for the originals, probably because the stricter CoverM mapping leaves fewer species detected at ≥ 3 of 4 sites. The heatmap is not in the manuscript.

**Unbinned reads** (sample means): SF2 89 %, R2A 85 %, R2 60 %, R1 51 %.
