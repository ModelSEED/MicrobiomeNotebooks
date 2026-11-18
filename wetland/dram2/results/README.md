# DRAM2 results: 1,184 high-quality, multi-sample wetland MAG representatives (2026-11-17)

DRAM2 v2.0.0-beta17 was run natively on poplar (setup and scripts in `../`).
- **Inputs:** one representative per species cluster, CheckM2 >90% complete / <5% contaminated, cluster seen in ≥2 samples (`census/outputs/dram2_hq_multisample.tsv`).
- **Gene calling:** DRAM2's own Prodigal.
- **Databases:** KOfam, dbCAN, MEROPS, VOG, CAMPER, CANT-HYD, methyl, sulfur, FeGenie.
- **Runs:** 5 chunks, all verified, then one distill/product run over all 1,184 MAGs.

Files here:
- `DISTILL/distillate.xlsx`: DRAM2 distill (topics default, camper; ecosystems eng_sys, ag).
- `PRODUCT/product.tsv`: DRAM2 product (module/function completeness per MAG).
- **Coding and non-coding genes** (added 2026-11-18). DRAM2 writes them separately: coding genes in the raw
  annotations, non-coding genes in its tRNA/rRNA scans.
  - `noncoding/`: DRAM2's non-coding calls, namely `combined_trna_scan.tsv.gz` (tRNAscan-SE, per tRNA),
    `combined_rrna_scan.tsv` (barrnap 5S/16S/23S), and the per-MAG count matrices `collected_{t,r}rnas.tsv`.
  - `dram2_per_mag_genes.tsv`: per MAG counts of CDS, CDS with any DB hit, with a KO, with a threshold-passing KO,
    tRNA, distinct tRNA amino acids, and 5S/16S/23S rRNA.
  - **Totals (1,184 MAGs):**
    - 3,935,612 CDS, of which 2,946,753 have a DB hit;
    - 53,420 tRNA, present in every MAG (median 19 amino acids per MAG);
    - 925 5S, 748 16S and 525 23S rRNA.
  - DRAM2 calls no other non-coding RNA (no tmRNA, Rfam ncRNA or CRISPR). The KBDL Bakta builds carry those.
  - `dram2_genes_coding_noncoding.tsv.gz` (103 MB, local only): one row per gene with gene_type CDS/tRNA/rRNA,
    location and product. Built by `../merge_coding_noncoding.py`.
- **Residual outputs** (local only, `residual/`; listed in `RESIDUAL_MANIFEST.tsv` with sha256): everything else
  DRAM2 wrote, i.e. Prodigal gene calls (faa/fna/gff), per-database HMM and MMseqs2 hits, raw KOfam hits (for exact
  re-thresholding), per-chunk raw annotations, QUAST, processed inputs, run logs and configs. Bundled on poplar by
  `../bin/bundle_dram2_residual.sh` and pulled by `kbdl_build/auto_collect.sh`.

Full outputs on poplar, `/scratch1/afreiburger/wetland_dram2/`:
- `merged/raw-annotations.all.tsv` (1 GB; `mag` column = original MAG name);
- `runs/all_mags/out/{DISTILL,PRODUCT}`;
- raw KOfam hits in `runs/chunk_NN/out/HMM_search_kofam_raw/`.

**Caveat: KO calls are not thresholded by DRAM2 beta17.** It keeps the best KOfam hit per gene at E ≤ 1e-5.
- Of 3,935,612 genes, 2,896,210 (74%) got a KO, and only 1,643,919 of those (57%) pass the KO's own KOfam score threshold.
- distill and product use the unfiltered calls, so module completeness here is optimistic.
- For KO-level work, filter the merged table on `kofam_passes_ko_threshold`. It is approximate for "domain"-type KOs; the raw hits allow exact re-thresholding.
