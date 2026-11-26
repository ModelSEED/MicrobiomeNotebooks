# DRAM2 results: 1,206 high-quality, multi-sample wetland MAG representatives (run 4, 2026-11-25)

DRAM2 v2.0.0-beta17 was run natively on poplar (setup and scripts in `../`).
- **Inputs:** one representative per species cluster, CheckM2 >90% complete / <5% contaminated, cluster seen in ≥2
  samples, from the run-4 census of 13,164 MAGs (`census/outputs/dram2_hq_multisample_run4.tsv`).
  - The same rule on the first census (12,859 MAGs) gave the original 1,184.
  - Run 4 replaced 12 of them with higher-quality MAGs of the same species and added 22 newly qualifying species:
    34 additions, run as chunk 6.
- **Gene calling:** DRAM2's own Prodigal.
- **Databases:** KOfam, dbCAN, MEROPS, VOG, CAMPER, CANT-HYD, methyl, sulfur, FeGenie.
- **Runs:** 6 chunks (1,218 MAGs; 12 of them superseded), all verified. The merged table was cut to the 1,206
  (`../bin/subset_merged.py`), then one distill/product run went over exactly those (`runs/all_mags_v2`;
  `../bin/finalize1206.sh`).

Files here (the 1,206 set):
- `DISTILL/distillate.xlsx`: DRAM2 distill (topics default, camper; ecosystems eng_sys, ag).
- `PRODUCT/product.tsv`: DRAM2 product (module/function completeness per MAG).
- **Coding and non-coding genes.** DRAM2 writes them separately: coding genes in the raw annotations, non-coding
  genes in its tRNA/rRNA scans.
  - `noncoding/`: DRAM2's non-coding calls, namely `combined_trna_scan.tsv.gz` (tRNAscan-SE, per tRNA),
    `combined_rrna_scan.tsv` (barrnap 5S/16S/23S), and the per-MAG count matrices `collected_{t,r}rnas.tsv`.
  - `dram2_per_mag_genes.tsv`: per MAG counts of CDS, CDS with any DB hit, with a KO, with a threshold-passing KO,
    tRNA, distinct tRNA amino acids, and 5S/16S/23S rRNA.
  - **Totals (1,206 MAGs):**
    - 4,007,770 CDS, of which 2,997,696 have a DB hit;
    - 54,428 tRNA, present in every MAG (median 19 amino acids per MAG);
    - 943 5S, 757 16S and 530 23S rRNA.
  - DRAM2 calls no other non-coding RNA (no tmRNA, Rfam ncRNA or CRISPR); Bakta (`../../bakta/results`) does.
  - Two local-only tables:
    - `raw-annotations.all.tsv.gz` (207 MB): DRAM2's per-CDS annotations;
    - `dram2_genes_coding_noncoding.tsv.gz`: one row per gene with gene_type CDS/tRNA/rRNA, location and product,
      built by `../merge_coding_noncoding.py`.
- `v1_1184/`: the same files for the original 1,184-MAG set (2026-11-17/18), kept for traceability.
- **Residual outputs** (local only, `residual/`; listed in `RESIDUAL_MANIFEST.tsv` with sha256): everything else
  DRAM2 wrote for all 6 chunks and both distill runs. That covers:
  - Prodigal gene calls (faa/fna/gff);
  - per-database HMM and MMseqs2 hits, and the raw KOfam hits (for exact re-thresholding);
  - per-chunk raw annotations, QUAST and the processed inputs;
  - run logs and configs.

Full outputs on poplar, `/scratch1/afreiburger/wetland_dram2/`:
- `merged/raw-annotations.all.tsv` (all 1,218 MAGs run) and `merged_set1206/` (the set);
- `runs/all_mags_v2/out/{DISTILL,PRODUCT}`, plus `runs/all_mags/` (the v1 distill);
- raw KOfam hits in `runs/chunk_NN/out/HMM_search_kofam_raw/`.

**Caveat: KO calls are not thresholded by DRAM2 beta17.** It keeps the best KOfam hit per gene at E ≤ 1e-5.
- Of 4,007,770 genes, 2,945,888 (74%) got a KO, and only 1,671,868 of those (57%) pass the KO's own KOfam score
  threshold.
- distill and product use the unfiltered calls, so module completeness here is optimistic.
- For KO-level work, filter the raw table on `kofam_passes_ko_threshold`. It is approximate for "domain"-type KOs;
  the raw hits allow exact re-thresholding.
