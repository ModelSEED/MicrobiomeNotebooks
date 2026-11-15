# Weisenhorn wetland MAG census: results (2026-11-13)

This covers handoff Phases 3–5 as revised in `HANDOFF_REVIEW.md`:
- **Population:** the 12,859 MAGs of KBase narrative 272312, from 702 of the 770 metagenome assemblies.
- **Quality:** CheckM2 1.1.0 for every MAG.
- **Taxonomy:** GTDB-Tk 2.8.0 / GTDB R232.
- **Clustering:** skani 0.3.1, ANI ≥ 95% and max alignment fraction ≥ 50% for species, ≥ 99% for strains.

Outputs are in `~/gtdbtk_kbase/census_final/`, with a copy in `census/outputs/` in this repo:
- `mag_census.tsv`: one row per MAG, with all three quality estimates, R232 fields, cluster ids, novel-DB / r220 hits and R232-adjacency labels;
- `clusters.tsv`, `strain_edges.tsv`, `qc_known_novel_edges.tsv`;
- `headline.json`;
- `rep_worklist.tsv`, `cluster_reps.tsv` and `rep_summary.json` (Phase 5).

Poplar compute is in `/scratch1/afreiburger/wetland_census`.

## Phase 4: quality (CheckM2)
| | MAGs |
|---|---|
| High quality (> 90% complete, < 5% contaminated) | 4,198 |
| Medium quality or better (≥ 50 / < 10) | 11,070 |
| Below medium quality by CheckM2 (all passed DAS Tool's ≥ 50/<10) | 1,789 |
| Model used | General (gradient boost) 10,676; specific (neural network) 2,183 |

**CheckM2 rescues CPR.** For the 471 Patescibacteriota MAGs, 78.6% pass > 70/<10 under CheckM2, against 17.4% under
KBase's CheckM v1, with median completeness 89.3 vs 57.6. The DPANN phyla (Nanobdellota, Aenigmatarchaeota,
Micrarchaeota, Nanohalarchaeota, SpSt-1190, Altiarchaeota; 147 MAGs) go from 42.2% to 68.0%, driven by Nanobdellota
(67 MAGs, 17.9% → 64.2%). All other MAGs are nearly unchanged (68.9% → 70.6%). Spearman correlation of v1 vs v2
completeness is 0.88 over all 12,859 MAGs. CheckM2 is therefore
the authoritative quality source; CheckM v1 and DAS Tool SCG are kept as columns.

## Phase 3: phylogenetic census
| | MAGs | Species clusters |
|---|---|---|
| Known GTDB R232 species | 6,497 (50.5%) | 1,893 (31.2%) |
| New species-level clusters (de novo, 95% ANI / 50% AF) | 6,362 | 4,171 |
| **Total** | **12,859** | **6,064** (6,889 strain clusters at 99%) |
| Clusters with ≥ 1 MAG at > 70/<10 | | 4,160 (1,589 known = 38.2%) |

- **No GTDB reference at all:** for 4,481 MAGs (34.8%), GTDB-Tk reports no R232 reference genome of any kind (no
  closest genome, placement reference, related reference or AF warning). This is a result about these communities, not a failure.
- **Novelty depth of the new clusters,** taken from the first unnamed R232 rank:
  - new species in named genera: 3,602;
  - new genera: 495;
  - new families: 47;
  - new orders: 17;
  - new classes: 4;
  - new phyla: 6.

  The deep calls are mostly single low-quality MAGs. Of the 27 clusters at order level or deeper, only 6 have a
  member at > 70/<10, so claim deep novelty only for those.
- **By habitat**, as the share of MAGs in a known species: permafrost 59.3%, peatlands 56.3%, wetlands 42.4%.
  - The wetland metagenomes are the least characterized.
  - Peatlands contribute 6,736 MAGs across 2,822 clusters, and wetlands 5,475 across 2,864.
- **Cluster sizes:** 3,856 of the 6,064 clusters are singletons and 154 have ≥ 10 MAGs. The largest are *Methanoflorens
  stordalenmirensis* (61 MAGs), Fen-455 sp003140595 (60) and Bog-375 sp003142185 (58).
- **Novel-genome references, used as labels:**
  - `novel_all_representatives` (78,132 genomes: 74,005 SPIRE/MGnify MAGs, about 95%, plus 4,127 AllTheBacteria/BV-BRC isolates): 2,707 MAGs have a species-level hit, and it contains 506 of the 4,171 new clusters, so about 12% are already known from other metagenome studies.
  - The handoff's `novel_hq_v2` (isolates): only 4 MAGs hit it, and it contains none of the new clusters.
- **The handoff's GTDB r220 route vs R232:** r220 puts 3,141 MAGs in a species, against 6,497 under R232. 3,089 are in
  both, 3,408 are known only to R232 and 52 only to r220. An r220-based census would have reported about 61% (1,156 of 1,893)
  of the known species as novel, and 52.5% of the known-species MAGs.
- **Boundary cases, which are kept visible:**
  - **MAGs next to a known species:** 120 new-cluster MAGs touch a known species through another MAG. GTDB-Tk explicitly rejected 109 of them against the species representative: 69 below the ANI radius, and 40 at or above it but with max AF < 0.5. For the other 11, GTDB-Tk reports no value to that representative. Of the 40 low-AF MAGs, 10 are < 50% complete, 10 are ≥ 50% complete with ≥ 10% contamination, and 20 are ≥ 50% / < 10%; why they fail the AF test is not established.
  - **Low-AF-only clusters:** 35 new clusters consist only of such low-AF MAGs, so the "new clusters" count excluding them is 4,136.
  - **Known↔new edges:** `qc_known_novel_edges.tsv` lists the 504 species-level edges between known and new clusters; 128 of them are also strain-level (≥ 99% ANI).

## Phase 5: representative worklist for the functional census
| | |
|---|---|
| Worklist | **5,356 MAGs** from 4,996 clusters: 84 named phyla plus 2 representatives with no phylum assignment (4,919 bacterial, 437 archaeal) |
| Per cluster | 1 rep in 4,756 clusters; 2–10 in 240 (one per 99% strain, max 10; only 1 cluster hits the cap) |
| Relaxed reps | 902 clusters had no eligible MAG (> 70/<10 without a GTDB-Tk multiple-hit warning; 901 have no MAG at > 70/<10 at all, and 1 has only a flagged one), so they get one rep at ≥ 50/<10, flagged `relaxed` |
| Census-only clusters | 1,068: 1,057 with nothing ≥ 50/<10 without a GTDB-Tk multiple-hit warning, and 11 whose MAGs at ≥ 50/<10 are all ≥ 99% ANI to a rep already chosen (7 of them have MAGs at > 70/<10; the other 4 have only floor-level MAGs) |
| MIMAG high quality by completeness/contamination | 2,161 of the worklist (rRNA/tRNA not yet called; Phase 2 build provides them) |
| Excluded for GTDB-Tk "markers with multiple hits" (contamination/chimera flag) | 405 MAGs |
| Habitats | peatlands 2,521, wetlands 2,475, permafrost 360 |

Rules are in `select_reps.py`:
- **Score:** completeness − 5 × contamination + 0.5·log10 N50.
- **Strains and sites:** one MAG per strain, preferring a new sampling site (coordinates merged within 2 km) as a tie-break.
- **No duplicates:** no two representatives are linked at ≥ 99% ANI with max alignment fraction ≥ 50% (the strain
  rule, `strain_edges.tsv`), in any cluster. 35 representative pairs, each in a different species cluster, reach
  ≥ 99% ANI only over max AF < 50%.

The checks pass: every MAG is in exactly one cluster, strains are nested in species, and no representative pair is ≥ 99% ANI at max AF ≥ 50%.

## Cost of the next phases (estimates)
- **Phase 2, KBDL genome builds for the 5,356 representatives:** roughly 1–2 days of queue with `bakta`, or
  roughly 4–7 h with `inprocess`, at the service's limit of 4 concurrent builds, which is shared with other users.
  These are rough: they are scaled by MAG count (× 0.42) from extrapolated, unmeasured whole-set estimates
  (2.6–4.9 days for `bakta`; 8–17 h for `inprocess`, depending on whether the per-genome skani search is counted).
  The pilot build should measure the real rate.
- **Phase 6, annotation of the worklist:** KBase DRAM at ≥ 25 min per MAG, run one after another, would take at least
  93 node-days (5,356 × 25 min). That needs a batch route such as DRAM2 on h100 or KBDLGenomeAnnotation.
