# DRAM2 results: 1,184 high-quality, multi-sample wetland MAG representatives (2026-11-17)

DRAM2 v2.0.0-beta17 was run natively on poplar (setup and scripts in `../`).
- **Inputs:** one representative per species cluster, CheckM2 >90% complete / <5% contaminated, cluster seen in ≥2 samples (`census/outputs/dram2_hq_multisample.tsv`).
- **Gene calling:** DRAM2's own Prodigal.
- **Databases:** KOfam, dbCAN, MEROPS, VOG, CAMPER, CANT-HYD, methyl, sulfur, FeGenie.
- **Runs:** 5 chunks, all verified, then one distill/product run over all 1,184 MAGs.

Files here:
- `DISTILL/distillate.xlsx`: DRAM2 distill (topics default, camper; ecosystems eng_sys, ag).
- `PRODUCT/product.tsv`: DRAM2 product (module/function completeness per MAG).

Full outputs on poplar, `/scratch1/afreiburger/wetland_dram2/`:
- `merged/raw-annotations.all.tsv` (1 GB; `mag` column = original MAG name);
- `runs/all_mags/out/{DISTILL,PRODUCT}`;
- raw KOfam hits in `runs/chunk_NN/out/HMM_search_kofam_raw/`.

**Caveat: KO calls are not thresholded by DRAM2 beta17.** It keeps the best KOfam hit per gene at E ≤ 1e-5.
- Of 3,935,612 genes, 2,896,210 (74%) got a KO, and only 1,643,919 of those (57%) pass the KO's own KOfam score threshold.
- distill and product use the unfiltered calls, so module completeness here is optimistic.
- For KO-level work, filter the merged table on `kofam_passes_ko_threshold`. It is approximate for "domain"-type KOs; the raw hits allow exact re-thresholding.
