# DRAM 1.5.0 on the 558 species representatives: results, validation against KOfam, and the manuscript summary

*Track `dram`, results written 2026-11-13. Run set-up and history: `README.md` in this folder. Release: **KOfam 2026-10-26** (the same files as the KOfam track), Pfam 39.0, dbCAN V11, MEROPS 12.x; Prodigal 2.6.3 `meta`; contigs ≥1,000 bp. Nothing was committed and nothing on poplar was changed or deleted.*

## Review (2026-11-13)

*A second agent reviewed this file.*
- *The headline counts below were re-derived with a separate script from `merged/annotations.tsv.gz`, `distill/product.tsv`, `distill/genome_stats.tsv`, `metabolism_summary.xlsx` and the logs.*
- *The per-species KO attributions behind the product-row changes (items 3, 4 and 6) and the phylum breakdowns use the `ids:` columns of `species_558_dram_summary.tsv` and `summary/*.tsv`. `summarize_species.py` derives these from the same `ko_id` column.*

*This section lists what was checked and what was changed. No count in any table changed, and no TSV, script or file under `results/` was modified. Nothing was committed, and poplar was not touched.*

**Re-derived independently, all equal to the text:**
- **`merged/annotations.tsv.gz`:**
  - 1,645,407 rows, 558 `fasta` values and unique gene IDs.
  - 673,349 genes with `ko_id` (40.92 %; per-MAG median 40.81 %, range 11.82–64.75 %). Rank C/D/E is 673,349 / 613,410 / 358,648. Pfam is 75.57 % as a mean per MAG.
  - 511 bacteria and 47 archaea by `bin_taxonomy`.
- **`zcat merged/annotations.tsv.gz | md5sum`** = `e6f025947f9e82c07b67ab1f194aa423`, the poplar md5. The other 595 entries of `MD5SUMS_poplar.txt` pass `md5sum -c`.
- **Per-MAG logs:** 558 have `Exit status: 0`, none has a non-zero exit, and none mentions "database is locked", "disk I/O error" or a Traceback. The `genes=` values in `run_all.log` sum to 1,645,407. Wall time per MAG is 6.4–77.9 min (median 26.9) and 270.1 h summed; peak RSS is 25,033,448 kB.
- **Sheet KOs per species, from `ko_id`:**
  - *phnJ* 35, *grdH* 4, *grdI* 15, *mcrA* 1.
  - K14083: 136 species and 417 genes. 135 have no *mcrA*, 106 have more than one copy (up to 16), 23 are Desulfobacterota, and all 38 Rhodobacteraceae have it (113 genes).
  - Betaine reducers (DRAM): 15. TMA producers 0, one-sided candidates 27.
  - *mtgB* 16 (8 also reducers); *cutC* 4, *cntA* 23, *cutD* and *cntB* 0.
  - The "Salterns methane" sheet itself gives the same *phnJ*, *grdH*, *grdI* and K14083 counts.
- **`product.tsv`:**
  - Every non-zero count in §4: TMA→DMA 136 (135 without *mcrA*, 38 Rhodobacteraceae), *hao* 77 (41 Desulfobacterota; only `concoct.bin.031_sub` also has *amo*), nitrite→nitrate 48, *amo* 6, pmo 6 (the same 6 genomes), ammonia→nitrite 3, *fwdA* 16, acetate 298 / 285 / 55, SOX 94, *dsr* 51, *nif* 33.
  - Every CAZy row.
  - The modules: CO2→methane 201 non-zero and 1 complete, acetyl-CoA 39 and 2, Wood–Ljungdahl 443 and 0.
- **`genome_stats.tsv`:** 547 med and 11 high; 5S in 300, 16S in 235, 23S in 44.
- **Sec reconciliation (`validation/mag_marker_cells_differing_post_rescue.tsv`, recounted by hand):**
  - *grdH*: 12 MAGs, −12.
  - *grdA*: 19 MAGs, +28.
  - *grdB*: 13 MAGs, +11. The 12 positive *grdB* MAGs are exactly the 12 *grdH* MAGs.
  - 22 MAGs in all.
  - The 16 `betaine_reducer_final` species split by evidence as stated in §3: 10 / 4 / 1 / 1.

**Marker reconciliation logic** (`analysis/validate_vs_kofam.py`, `summarize_species.py`): sound.
- Genes are matched on MAG, contig, start, end and strand.
- The pre-rescue counts are cross-checked against two KOfam tables.
- `summarize_species.py` asserts that every DRAM guild equals the KOfam guild on the same MAG. For betaine it compares against `betaine_reducer_without_sec_rescue`.
- One gap was not stated: the post-rescue comparison reads `guilds_1324.tsv` for only 29 of the 68 markers. The other 39 equal the pre-rescue counts by construction. That is safe (all 100 replaced pieces in the 558 had a Sec best KO), and §2 now says so.

**The K14083 and Sec caveats** are stated in the short answer and in §3, §4 and §6.

**Changes made (wording and labels only):**
1. §1 distill row: "558 genomes × 98 functions" → "98 columns (32 module and ETC completeness fractions, 66 present/absent functions)". This is the dtype split of `product.tsv`.
2. §2: added the 29/39 note to the post-rescue row, citing `kofam/results/sec_readthrough/replaced_prodigal_piece_calls_new.tsv`.
3. §4 ammonia oxidation row: "(*amoABC*)" → the row's three KOs, with the note that only *amoB*/*amoC* are called (`ids:` columns of `species_558_dram_summary.tsv`: K10945 ×5, K10946 ×5, no K10944). BJGV01 rests on *amoC* plus *hao*, and *Nitrosopumilus* `maxbin2.Bin.029` on *amoB* alone. The count of 7 ammonia oxidizers stands on taxonomy plus these partial *amo* calls.
4. §4 "nitrite => nitric oxide (*nirK*/*nirS*)" → "(*nirK*/*nirS*, or *nirB* K00362)". DRAM's form row includes K00362 (`distill_forms/function_heatmap_form.20261111.tsv`), and 5 species are lit by it alone: Flavobacteriaceae `concoct.bin.170`, `metabat2.bin.149` and `metabat2.bin.280_sub`, and Oleiphilaceae `maxbin2.Bin.136` and `maxbin2.Bin.281_sub`.
5. §4 "nitric oxide => nitrous oxide (*norB*/*norC*)" → "(*norB* K04561)". The form row is K04561 + K15877, and K15877 is fungal P450nor (`kofam/results/hmm_info.tsv`), not *norC* (K02305). All 88 species are lit by K04561.
6. §4 nitrogen fixation: added that only 23 of the 33 have both K02586 and K02588 (6 *nifH* only, 4 *nifD* only).
7. §5/§6, the two 136s: "Report these as `mttB_family_without_mcrA`" now points at `k14083_caveat` (the 135). The `mttB_family_without_mcrA` column is 136 because it adds `maxbin2.Bin.023_sub` (M55B127, *mtgB* only) and excludes *Methanolobus*. A §5 bullet explains this.
8. §6 any-of bullet: added the *nirB*, *norB* and single-subunit *nif* cases.

**Overclaims to watch, unchanged in the text:**
- The 4 non-zero-E best-KO ties are explained as "most likely" DRAM's 2-digit E-value precision. That is plausible but not shown, and none is a marker.
- The bacterial *mtbB* calls clear their threshold by 34–39 bits (531/536 vs 497). That is moderate, not "comfortable"; the *mtmB* calls clear by 38–225 bits.
- The README's "`state/DONE` at 02:15" is the driver's clock. The marker files' mtimes (preserved by rsync) are about 4 min earlier (DONE 02:11, MERGE_OK 00:41 vs the logged 00:45), which looks like a constant host vs NFS clock offset. No step is out of order.

## Short answer

- **The run finished cleanly.** `state/DONE` exists; there is no `DONE_INCOMPLETE` or `FAILED`, and `state/missing.txt` is empty. All **558/558** MAGs were annotated in pass 1 with exit status 0, and **1,645,407 genes** were merged.
- **Everything needed is copied back and md5-checked:** 596 files are identical to poplar. The merged `annotations.tsv` is stored gzipped (662 MB → 88 MB), and its uncompressed md5 matches poplar.
- **DRAM's KO layer is identical to the KOfam track on the same 558 MAGs.**
  - The gene calls are identical: 1,645,407 genes in both, every gene matched by contig, start, end and strand. Both tracks used the same Prodigal build, `prodigal-2.6.3-h577a1d6_11`.
  - **All 3,785 marker genes have the same best KO.**
  - The 558 × 68 MAG × marker counts are **identical (0 of 37,944 cells differ)** to the KOfam counts before the Sec rescue.
  - Across all KOs, 50 of 673,349 genes with a KO get a different best KO. Both KOs pass in every one of them; these are E-value ties (46 at E = 0), and none involves a marker.
- **The only DRAM ≠ KOfam differences are the Sec read-through rescue.** 44 cells in 22 MAGs differ, all *grdH*, *grdA* or *grdB*:
  - DRAM has 4 *grdH* genes in 4 species, against 16 genes in 15 species with the rescue;
  - DRAM counts 28 extra *grdA* genes: both pieces of split loci are called;
  - DRAM counts 11 extra *grdB* genes: N-terminal GrdH pieces are called *grdB*.
- **Guild calls on the best MAG of each species (558).** Rules as in the KOfam track; the betaine column adds the Sec-rescued *grdH* from `guilds_1324.tsv`:

  | guild | DRAM | final (DRAM + Sec rescue) |
  |---|---|---|
  | MePn consumers (*phnJ*) | 35 | 35 |
  | betaine reducers (*grdH* or *grdI*) | **15** | **16** |
  | methanogens (*mcrA*) | 1 | 1 |
  | methylamine users | 1 | 1 |
  | TMA producers | **0** | 0 |
  | one-sided TMA candidates | 27 | 27 |

  - The 16th betaine reducer is *QPJE01* (Halanaerobiaceae) `metabat2.bin.356.fasta__R1_sediment`. Its only evidence is the Sec-joined *grdH*: no *grdI*, and DRAM calls the N-terminal piece *grdB*.
- **The K14083 caveat applies in full.** DRAM's product heatmap lights "trimethylamine => dimethylamine" in **136 species**. **135 of them have no *mcrA*.**
  - They include all 38 Rhodobacteraceae species and 23 Desulfobacterota species.
  - Only *Methanolobus* sp020055655 is a methylamine methanogen.
  - Read these 135 as the non-pyrrolysine MttB family, not as TMA methanogenesis.

## 1. Run completion and copy-back

| check | result |
|---|---|
| stage markers (`results/state/`) | `SETUP_A_OK`, `SETUP_OK`, `TEST_OK`, `ANNOTATE_OK`, `MERGE_OK`, `DISTILL_OK`, **`DONE`**; no `DONE_INCOMPLETE` or `FAILED`; `missing.txt` empty |
| driver log | `pass 1: missing 0` (2026-11-13 00:34), `merge OK: genes 1645407 MAGs in merged 558 MAG dirs 558` (00:45), `DONE: all MAGs annotated, merged, distilled` (02:15). The driver (PID 4007531) has exited; no DRAM processes are left |
| per-MAG logs (`logs/mags/*.log`) | 558 × `Exit status: 0`; no "database is locked", "disk I/O error" or Traceback |
| time | annotation from 2026-11-11 14:44 to 2026-11-13 00:34 (34 h wall at 8 parallel). Per MAG: median 27 min (6–78 min), 270 h summed, peak RSS 25 GB (Pfam MMseqs2). Merge took 11 min, distill 1 h, and the gene-names distill another 30 min |
| genes | 1,645,407 = sum of `genes=` in `run_all.log` = rows of `merged/annotations.tsv`, for 558 distinct `fasta` values |
| annotation depth | KO on 40.9 % of genes (per MAG median 40.8 %, range 11.8–64.8 %); Pfam on 75.6 % (mean per MAG); rank C 673,349, D 613,410, E 358,648 (no A/B, because KEGG and UniRef were not used) |
| distill | 9 sheets in `metabolism_summary.xlsx`: MISC, carbon utilization, Transporters, Energy, Organic Nitrogen, carbon utilization (Woodcroft), **Salterns methane** (65 rows × 558 genomes), rRNA, tRNA. `product.html`/`product.tsv` cover 558 genomes × 98 columns (32 module and ETC completeness fractions, 66 present/absent functions) |
| genome_stats | DRAM rates 547 MAGs "med" and 11 "high". 5S rRNA in 300 MAGs, 16S in 235, 23S in 44 |

**Copied** (`rsync -az`, then `md5sum -c` against md5s computed on poplar: **all 596 OK**; the list is in `results/MD5SUMS_poplar.txt`):

| local path (`results/`) | poplar source (`/scratch1/afreiburger/salterns_2026/dram/`) | size |
|---|---|---|
| `merged/annotations.tsv.gz` | `merged/annotations.tsv`, gzipped here. Check it with `zcat merged/annotations.tsv.gz \| md5sum` = `e6f025947f9e82c07b67ab1f194aa423` (the md5 of the poplar file) | 88 MB (662 MB raw) |
| `merged/trnas.tsv`, `merged/rrnas.tsv`, `merged/annotate.log` | `merged/` | 2.3 MB |
| `distill/` (`metabolism_summary.xlsx` with the "Salterns methane" sheet, `product.html`, `product.tsv`, `genome_stats.tsv`, `distill.log`) | `distill/` | 28 MB |
| `distill_gene_names/` (`metabolism_summary.xlsx` listing gene IDs instead of counts, `genome_stats.tsv`, `product.tsv`, `distill.log`) | `distill_gene_names/` | 11 MB |
| `env/`, `logs/` (incl. 558 per-MAG logs), `state/` | same | 5 MB |
| `distill_forms/` (DRAM's function heatmap, module step, genome summary and ETC forms, as installed) | `../dram_db/*.20261111.tsv` (md5 checked) | |

**Not copied (left on poplar):**
- `merged/genes.{faa,fna,gff}`, `merged/scaffolds.fna` and `merged/genbank/`;
- the per-MAG `dram_out/`;
- `distill_gene_names/product.html`. Its `product.tsv` is byte-identical to `distill/product.tsv`. The HTML differs only in the vega dataset hash names and in the order of IDs within tooltip `function_ids` strings (Python set order); this was checked on poplar by comparing the normalized token sets.

## 2. Validation against the KOfam track (`analysis/validate_vs_kofam.py` → `results/validation/`)

| check | expected | found |
|---|---|---|
| Prodigal build | same | both tracks: bioconda `prodigal-2.6.3-h577a1d6_11` (`env/dram15.explicit.txt`; KOfam `run_prodigal.sh` EmilyKin env) |
| genes per MAG, DRAM vs KOfam `genes_all.tsv.gz` | equal | **558/558 equal** (Δ = 0 everywhere); 1,645,407 = 1,645,407 |
| coordinate join (MAG, contig, start, end, strand) | all match | **1,645,407 matched**, 0 DRAM-only, 0 KOfam-only. Coordinates are the same 1-based inclusive values |
| genes with any KO | equal | 673,349 in both, and they are the same genes |
| best KO, all KOs | equal except E-value ties (README "Comparability") | **50 genes differ (0.007 %)**. In all 50 both KOs pass in both runs. 46 have E = 0 for the KOfam best KO, which DRAM breaks by profile order and KOfam by score. 4 are hits whose scores are within 0.6–2.3 bits of each other (KOfam E-values 1e-9 to 1e-176). DRAM's per-MAG E-values weren't kept, so the tie can't be shown directly; it is most likely at DRAM's 2-digit E-value precision. The pairs are paralog sets: K00164/K01616 ×12, K00123/K05299 ×11, K03296/K18138 ×6, K01537/K12953 ×4, … **None is a marker.** |
| best KO, 68 marker KOs | equal | **3,785/3,785 marker genes agree** (0 disagreements; `marker_gene_disagreements.tsv` is empty) |
| MAG × marker counts vs KOfam before the rescue (`all_mag_ko_counts.tsv` `n_best` = `all_mag_markers_guilds.tsv` `n_<marker>`) | equal | **0 of 37,944 cells differ** |
| MAG × marker counts vs KOfam after the rescue (`guilds_1324.tsv` `n_<marker>`) | differ only in the Sec KOs | **44 cells in 22 MAGs, all *grdH*, *grdA* or *grdB***; no non-Sec marker changes. This is tested directly for the 29 markers that have an `n_` column in `guilds_1324.tsv`. For the other 39 the post-rescue count is the pre-rescue count by construction; that is safe, because all 100 pieces the rescue replaced in these 558 MAGs had best KO *grdA* (61), *grdB* (37) or *grdH* (2) (`kofam/results/sec_readthrough/replaced_prodigal_piece_calls_new.tsv`) |
| KOs with threshold `-` (*mpnS* K18049, *phnY\** K23996, *phnZ1* K05915) | never called | 0 genes in DRAM and 0 as KOfam best KO. They are not in the custom sheet, which therefore has 65 of the 68 markers |
| "Salterns methane" sheet vs `ko_id` counts in `annotations.tsv` | equal | **0 of 36,270 cells differ** |
| product functions vs a re-implementation of DRAM 1.5 `make_functional_df` on `ko_id` + EC (non-CAZy rows) | equal | **0 mismatches** over 558 genomes; the contributing KOs per function are in `summary/product_functions_of_interest.tsv` |
| taxonomy/quality join (`genome_stats.tsv`) | the input tables | all 558 taxonomies equal `gtdb_taxonomy_558.tsv`, and completeness/contamination equal `checkm_quality_558.tsv`. None has DRAM's fallback values (0 / 100) |

**Every DRAM vs KOfam marker difference** (`validation/mag_marker_cells_differing_post_rescue.tsv`, DRAM − KOfam after the rescue). The reasons come from `kofam/RESULTS.md` §2.

| KO | MAGs | net genes | reason |
|---|---|---|---|
| *grdH* K21579 | 12 | −12 | Sec read-through not rescued by DRAM. Prodigal stops at the UGA codon, and neither piece reaches the 624.03 threshold. 11 species gain their first *grdH*; BM714 `maxbin2.Bin.114 R1_A_D2` gains a second copy |
| *grdA* K10670 | 19 | +28 | DRAM counts both pieces of a Sec-split *grdA* (the ~40-aa N-terminal piece and the C-terminal piece). One MAG, `concoct.bin.023 SF2_sediment` (UBA7803, Azotimanducaceae), is −1, because the rescue creates its only *grdA*/*grdB* call |
| *grdB* K10672 | 13 | +11 | The N-terminal piece of a Sec-truncated GrdH scores best as *grdB* (in 12 MAGs). The rescue replaces it with the joined *grdH*. The same `concoct.bin.023` is −1 (see *grdA*) |

There are no threshold or release differences: the two tracks share the KOfam files, the threshold rule (strict `>`, full or domain score) and the best-KO rule.

## 3. The "Salterns methane" sheet across the 558 species

Counts are genes with that best KO in each species' best MAG (`summary/salterns_methane_by_ko.tsv`; per species in `species_558_methane_guilds.tsv` and `species_558_salterns_methane_long.tsv`). Of the 558 species, 511 are bacteria and 47 archaea.

### MePn C-P lyase and other phosphonate catabolism

| gene | species | of which archaea | main families |
|---|---|---|---|
| *phnC* / *phnD* / *phnE* (transport) | 107 / 167 / 133 | 15 | Rhodobacteraceae 22–27, Haloferacaceae 8–9, Wenzhouxiangellaceae, Trueperaceae |
| *phnF* | 23 | 0 | Rhodobacteraceae 18 |
| *phnG* / *H* / *I* / ***J*** / *K* / *L* / *M* | 35 / 29 / 36 / **35** / 32 / 32 / 50 | 2–5 | Rhodobacteraceae 18–31, Haloferacaceae, Desulfobacteraceae |
| *phnN* / *phnO* / *phnP* | 37 / 1 / 219 | | *phnP* is not specific (219 species) |
| *phnX* / *phnW* / *phnY* / *phnZ* (2-AEP route) | 12 / 28 / 18 / 0 | 0 | SZUA-36, J153, Pirellulaceae, Rhodobacteraceae |
| *mpnS*, *phnY\**, *phnZ1* | not callable | | KOfam threshold `-`. The hits are in `kofam/results/all_no_threshold_hits.tsv.gz` for a manual cut-off |

- **35 MePn consumers (*phnJ*):**
  - 32 bacteria: Rhodobacteraceae 19 (*Roseovarius* 9, *Marivita* 2, *Roseibacterium*, *Lentibacter*, *Maritimibacter*, *Thermohalobaculum*, CAZRZN01, JAVJPJ01, JBKDZR01, DAINOF01), *Desulfotignum* 2, JAWLGC01, *Saliniramus*, JACKII01, CAJQQK01, JAAAOG01, CAMCWB01, WSZY01, DATHSN01, RI-231, UBA8366 and UBA8670;
  - 3 haloarchaea: *Haloquadratum*, *Halorubrum* and *Natronomonas*.
- **Locus completeness:**
  - 22 of the 35 carry all 7 *phnGHIJKLM*;
  - in 20 the locus is split over more than one contig;
  - 1 (JAWLGC01) has *phnJ* with fewer than 3 of the 7.
- **JAWLGC01** (`maxbin2.Bin.441.fasta__R2_sediment`, 99.4 % complete) gets its only *phnJ*/*phnK* from a 1,493 bp MaxBin2 contig. Its two cluster siblings have no *phn* gene. Treat this species' MePn call as a likely mis-binned contig (`kofam/RESULTS.md` §3).
- The per-cluster counts match KOfam exactly. In all 35 clusters the best MAG carries *phnJ*, so the 558-MAG count equals the any-MAG count.

### Betaine reductase and the Sec rescue

| gene | species (DRAM) | species (with the Sec rescue) | main families |
|---|---|---|---|
| *grdI* (guild rule) | 15 | 15 | M55B127 6, Dethiosulfatibacteraceae 2, JADQBU01, BM714, Acidaminobacteraceae, PWPR01, Halanaerobiaceae, Halarsenatibacteraceae, TCS61 |
| ***grdH*** (guild rule, Sec) | **4** | **15** | DRAM: BM714, DATFJG01, PWPR01, UBA8670 |
| *grdA* / *grdD* / *grdC* | 23 / 22 / 19 | 24 / 22 / 19 | the same Bacillota families, plus Desulfosalsimonadaceae (*grdA*) |
| *grdB* / *grdE* (glycine reductase paralogs) | 25 / 18 | 22 / 18 | |
| *grdF* / *grdG* (sarcosine reductase) | 0 / 0 | 0 | |
| *mtgB* (betaine:corrinoid methyltransferase) / *mtgA* | 16 / 0 | | M55B127 6, Desulfosalsimonadaceae 2, UBA11579, JADQBU01, UBA11574, JAAYZQ01; 8 of the 16 are also betaine reducers |
| betaine transport (*proVWX*, *opuA/BD/C*, *betT*, *opuD*, BCCT) | 192 / 193 / 218; 106 / 121 / 165; 117; 26; 70 | | widespread, not guild-defining |

- **Betaine reducers: 15 species by DRAM, 16 with the rescue.** All are Bacillota except DAOWES01 (Cloacimonadota TCS61):
  - M55B127 ×6;
  - UBA8670 ×2 (Dethiosulfatibacteraceae);
  - DATFJG01, BM714, PWPR01, CAYRVO01 (JADQBU01);
  - *Halanaerobium*, *QPJE01*, *Halonatronomonas*;
  - DAOWES01.
- **Evidence in the best MAG:**
  - *grdH*_sec_readthrough + *grdI*: 10;
  - *grdH* + *grdI*: 4;
  - *grdI* only: 1 (M55B127 `metabat2.bin.362_sub`);
  - Sec-rescued *grdH* only: 1 (*QPJE01* `metabat2.bin.356 R1_sediment`). Without the rescue this species is not a reducer in its best MAG.
- **Contig-edge *grdH*-like genes:** 2 of the 3 genes that KOfam still calls *grdB* are in these best MAGs: `metabat2.bin.362_sub R2_sediment` and `metabat2.bin.047 R1_A_D2`. Both are reducers through *grdI*, so no call depends on them.

### TMA production and consumption

| gene | species | notes |
|---|---|---|
| *cutC* / ***cutD*** | 4 / **0** | BM714 ×2, Desulfocapsaceae, Desulfobacteraceae (*Desulfotignum*). The *cutD*-like genes are PflA-type activases by best KO (`kofam/RESULTS.md` §3) |
| *cntA/yeaW* / ***cntB/yeaX*** | 23 / **0** | Rhodobacteraceae 6, JAKDMP01 5, J153, … |
| → TMA producers / one-sided candidates | **0 / 27** | 31 clusters if every MAG of a cluster counts (KOfam) |
| *torA* / *torZ* / *dmsA* (TMAO reduction) | 0 / 1 / 13 | DRAM's "TMAO reductase" product row (*torA*) is empty |
| *tmm* / *dmd-tmd* / *tdm* (TMA oxidation and demethylation) | 12 / 19 / 10 | 26 species with at least one; Rhodobacteraceae carry *tmm* in 8 species and *tdm* in 7 |

### Methanogenesis, methylamines, methanol, pyrrolysine

| gene | species | which |
|---|---|---|
| *mcrA* / *mcrB* / *mcrG* | **1** | ***Methanolobus* sp020055655** (`concoct.bin.101.fasta__R1_sediment`, 98.4 % complete; *mcrA* 983.3, *mcrB* 746.8, *mcrG* 423.5) |
| ***mttB*** K14083 | **136** | *Methanolobus* plus **135 species without *mcrA***: Rhodobacteraceae 38 (all 38 Rhodobacteraceae species), Desulfobacteraceae 6, M55B127 5, Woeseiaceae 4, Desulfosalsimonadaceae 4, UBA12077 4, …; 1 asgard archaeon (Promethearchaeaceae). 106 species carry more than one copy (up to 16). `summary/k14083_mttB_by_family.tsv` |
| *mttC* / *mtbC* / *mtbA* | 1 / 1 / 1 | *Methanolobus* only |
| *mtbB* K16178 | 2 | *Methanolobus*, and the bacterium **JAAYHV01** (Anaerovoracaceae, `metabat2.bin.069 R1_sediment`; 2 genes at 531 and 536 vs 497, with *pylBCD*) |
| *mtmB* K16176 | 4 | *Methanolobus*, and 3 bacteria: JAFDCJ01 (Desulfobacterota, 390, with *pylBCD*), UBA11579 and JAAYZQ01 (Chloroflexota, 298–484) |
| *mtaB* / *mtaC* / *mtaA* (methanol) | 3 / 2 / 2 | ***Methanohalophilus portucalensis*** (`metabat2.bin.187_sub R1_sediment`, 55 % complete, no *mcr* hit: *mtaB* ×2 at 730 and 791), plus 2 bacteria just above the 578.4 threshold (*Luteolibacter* 583.8, WJKA01 591.1) |
| *pylB/C/D/S* | 5 / 4 / 5 / 2 | *pylS* only in the 2 Methanosarcinaceae. *pylBCD* are also in JAAYHV01 and JAFDCJ01, and in *Desulfotignum* (`metabat2.bin.461`; *pylD* alone in `metabat2.bin.079`) |

- **Methanogens by the gene rule: 1 species** (*Methanolobus*); **methylamine users: 1** (the same).
- ***Methanohalophilus*** is a methanogen by taxonomy but not by the rule: an incomplete bin with the methanol branch and Pyl genes, but no *mcr*.
- **The bacterial *mtbB*/*mtmB* calls pass their KOfam thresholds comfortably**, and two of them sit with *pylBCD*. This looks like bacterial demethylation of methylamines (Pyl-containing methyltransferases), not methanogenesis. In DRAM's product they light "dimethylamine => monomethylamine" and "monomethylamine => ammonia" in genomes without *mcrA*.

## 4. DRAM product (`distill/product.html`, `product.tsv`): modules of interest

**How DRAM scores a function:** a function row is lit when *any* of its KOs or ECs is present, and a function is present when *all* its rows are. Modules (the left part of the heatmap) are a fraction of steps. Counts are species out of 558 (`summary/product_functions_of_interest.tsv`). The genomes behind every methylotrophy and monooxygenase call, with KOfam scores, are in `summary/notable_product_calls.tsv`.

### Methanogenesis and C1

| product row | species | driven by | reading |
|---|---|---|---|
| Methanogenesis: Key functional gene (*mcrABG*) | 1 | K00399/K00401/K00402 | *Methanolobus* |
| **trimethylamine => dimethylamine** | **136** | K14083 alone | **K14083 caveat:** 135 lack *mcrA*. The 38 Rhodobacteraceae species carry 113 copies. This is the non-Pyl MttB family (for example, glycine betaine or other quaternary-amine demethylation), not TMA methanogenesis. `product.html` shows these cells lit; read them with the *mcrA*/Pyl rule (`annotation_choice/README.md` §6.3). The column `k14083_caveat` in the species tables marks them |
| dimethylamine => monomethylamine | 2 | *mtbB* | *Methanolobus*; JAAYHV01 (bacterial, Pyl) |
| monomethylamine => ammonia | 4 | *mtmB* | *Methanolobus*; JAFDCJ01, UBA11579, JAAYZQ01 (bacteria) |
| methanol => methane | 3 | *mtaB* | *Methanohalophilus* (no *mcr*); *Luteolibacter* and WJKA01 at threshold +5 and +13 bits (weak) |
| putative but not defining CO2 => methane | 16 | *fwdA/fmdA* K00200 | Pirellulaceae 8, *Methanolobus*, *QPJE01*, *Methyloceanibacter*, … *fmd/fwd* also serves formaldehyde/formate oxidation, so this is not methanogenesis |
| acetate => methane, pt 1 / 2 / 3 | 298 / 285 / 55 | *acs* / *ack* / *pta* | generic acetate metabolism, **not methanogen-specific** |
| methane => methanol (pmo) / (mmo) | 6 / 0 | *amoB/C* K10945/K10946 | the same 6 genomes as "Bacterial/Archaeal ammonia oxidation" (ammonia oxidizers, below), not methanotrophs |
| module "Methanogenesis, CO2 => methane" | 201 non-zero, **1 complete** | module steps | complete only in *Methanolobus*. The 200 partial values are shared upstream steps |
| module "Acetyl-CoA pathway, CO2 => acetyl-CoA" | 39 non-zero, 2 complete | | complete in *Methanolobus* and *Methanohalophilus*; the partials are mostly Desulfobacterota (27) |
| module "Reductive acetyl-CoA pathway (Wood–Ljungdahl)" | 443 non-zero, 0 complete | | partial scores only (shared folate steps); don't quote them |
| "TMAO reductase" (*torA*) | 0 | | |

### Sulfur

| product row | species | main phyla / families |
|---|---|---|
| dissimilatory sulfate reduction (and oxidation), *dsrA*/*dsrB* | 51 | Desulfobacterota 34 (Desulfobacteraceae 8, Desulfocapsaceae 4, Desulfosalsimonadaceae 4, UBA11574 4), Pseudomonadota 9, Gemmatimonadota, Zixibacteria, Acidobacteriota. DRAM does not separate reductive from reverse (oxidative) Dsr; the 9 Pseudomonadota are candidates for rDsr |
| thiosulfate oxidation by SOX | 94 | Pseudomonadota 77 (Rhodobacteraceae 35 of 38) |
| thiosulfate => sulfite (*phsA*/*psrA*) | 64 | Desulfobacterota 27, Bacteroidota 16 |
| tetrathionate => thiosulfate (*ttrA*) | 2 | Nitrococcaceae, Geopsychrobacteraceae |

### Nitrogen

| product row | species | notes |
|---|---|---|
| Bacterial/Archaeal ammonia oxidation (row K10944/K10945/K10946; called here only as *amoB* K10945 and *amoC* K10946, with no K10944 *amoA* in any of the 558) | 6 | Nitrosopumilaceae 4, Nitrososphaeraceae 1 and the bacterium BJGV01 (Nitrosomonadaceae). BJGV01 (`concoct.bin.031_sub`) rests on *amoC* alone plus *hao*, and *Nitrosopumilus* `maxbin2.Bin.029` on *amoB* alone. The thaumarchaeal *amoA* is KEGG's newer K28504, which DRAM's form doesn't list. Through its EC (1.14.99.39) it lights "ammonia => nitrite" in 3 Nitrosopumilaceae, including a 7th ammonia oxidizer, *Nitrosomaritimum* `metabat2.bin.012_sub`. So there are **7 ammonia-oxidizer species** in total |
| **Bacterial (aerobic-specific) ammonia oxidation** (*hao* K10535) | **77** | **Caveat: *hao* alone.** Desulfobacterota 41, Pseudomonadota 10, Gemmatimonadota 6, …; KOfam scores 295–506 vs 279. Only BJGV01 also has *amo*. The rest are octaheme HAO-family proteins (for example, hydroxylamine or nitrite reductases of sulfate reducers), not nitrifiers |
| nitrite => nitrate (*nxrAB* = *narGH* KOs) | 48 | **Caveat:** these are the same KOs as nitrate reductase (K00370/K00371). They are not evidence of nitrite oxidation |
| nitrate => nitrite (*narGH*/*napAB*) | 69 | Pseudomonadota 34, Desulfobacterota 14 |
| nitrite => nitric oxide (*nirK*/*nirS*, or *nirB* K00362) | 117 | Pseudomonadota 51, Bacteroidota 23, Thermoproteota 6 (all Nitrososphaeria: thaumarchaeal *nirK*). **Caveat:** 5 species (3 Flavobacteriaceae, 2 Oleiphilaceae) are lit by *nirB* alone, the NADH nitrite reductase that makes ammonia. They are not denitrifiers |
| nitric oxide => nitrous oxide (*norB* K04561) | 88 | Pseudomonadota 52. All 88 are *norB*; the row's other KO, K15877, is fungal P450nor |
| nitrous oxide => nitrogen (*nosZ*) | 51 | Bacteroidota 28, Pseudomonadota 8, Gemmatimonadota 7 |
| DNRA (*nrfA*/*nrfH*) | 72 | Bacteroidota 27, Desulfobacterota 12 |
| nitrogen fixation (*nifH*/*nifD*) | 33 | Only 23 have both subunits; 6 have *nifH* only and 4 *nifD* only. Desulfobacterota 20, Planctomycetota 4; 1 archaeon (*Methanolobus*). Alternative (*vnfD*) 0; anammox (*hzs*) 0 |

### CAZy (product rows; present/absent)

Lit in at least one species: mixed-linkage glucans 151, amorphous cellulose 101, arabinan 89, xyloglucan 71, polyphenolics 13 and sulfated polysaccharides 11. These come mainly from Flavobacteriaceae, Balneolaceae, Prolixibacteraceae, Saprospiraceae and UBA12077.

13 rows are lit in no species: alpha-galactans, alpha-mannan, arabinose cleavage, beta-galactan, beta-mannan, chitin, crystalline cellulose, fucose, mucin, pectin, rhamnose, starch and xylans.

The full dbCAN calls per gene are in `annotations.tsv.gz` (`cazy_*` columns) and in the "carbon utilization" sheets of `metabolism_summary.xlsx`.

## 5. Guild table for the manuscript (best MAG of each of the 558 species; `summary/guilds_dram_vs_kofam_558.tsv`)

| guild (KOfam-track rules, `annotation_choice/README.md` §6.3) | DRAM (no rescue) | KOfam, same MAG (with rescue) | KOfam, any MAG of the cluster |
|---|---|---|---|
| methanogen (*mcrA*) | 1 | 1 | 1 |
| MePn consumer (*phnJ*) | 35 | 35 | 35 |
| betaine reducer (*grdH* or *grdI*) | **15** | **16** | 16 |
| methylamine user (*mttB/mtbB/mtmB* in an *mcrA* genome) | 1 | 1 | 1 |
| TMA producer (*cutC+cutD* or *cntA+cntB*) | 0 | 0 | 0 |
| one-sided TMA candidate | 27 | 27 | 31 |
| MttB family without *mcrA* (not a guild) | 136 | 136 | 145 |

- Every DRAM guild call equals the KOfam call on the same MAG, except the one betaine cluster that depends on the Sec rescue. The script asserts this.
- **The two 136s are different sets.** The "trimethylamine => dimethylamine" product row (136) is the 135 K14083 species without *mcrA* plus *Methanolobus*. The `mttB_family_without_mcrA` column (136) is any of *mttB*/*mtbB*/*mtmB*/*mtgB* without *mcrA*. It is the same 135 plus `maxbin2.Bin.023_sub.fasta__R2_sediment` (M55B127), which has *mtgB* but no K14083, and it excludes *Methanolobus*. The 135 product-row cases are the ones `k14083_caveat` flags.
- **For the paper:** use the KO layer from DRAM or KOfam (they are identical), and **add the Sec-rescued *grdH* from `kofam/results/guilds_1324.tsv`**. The column `betaine_reducer_final` in the species tables already does this.

## 6. Caveats

- **K14083 (*mttB*) over-call.** It lights "trimethylamine => dimethylamine" in 135 non-methanogen species. It is a family-level KO; the 2026-10 threshold (218.93) is passed by many non-Pyl MttB paralogs. Report these 135 as the non-Pyl MttB family, not as methylamine use; `k14083_caveat` flags exactly them. (`mttB_family_without_mcrA` is a related but different set of 136; see section 5.)
- **Sec-truncated *grdH*:**
  - DRAM can't rescue it. Prodigal ends the gene at the UGA, and the pieces stay below 624.03. So DRAM under-calls *grdH* (4 vs 15 species) and over-calls *grdA* and *grdB*.
  - The KOfam track's rescue (`kofam/scripts/sec_rescue.py`) is the correction. It depends on the HMM spanning the U and on beating the pieces; a SECIS is recorded but not required.
- **No thresholds for *mpnS*, *phnY\** or *phnZ1*.** Neither DRAM nor the KOfam track calls them. MePn *synthesis* and the oxidative MePn route therefore can't be read from this run.
- **Product rows are any-of-KOs.**
  - pmo = amo, *hao* ≠ nitrification, *nxr* = *nar*.
  - *nirB* alone lights "nitrite => nitric oxide" (5 species), the *nor* row is *norB* only, and 10 of the 33 nitrogen-fixation calls rest on one *nif* subunit. The species-table columns `product_amo` and `product_nifHD` are these any-of rows.
  - The acetate rows and *fwdA* are not methanogenesis.
  - Most module scores are partial (Wood–Ljungdahl 0 complete).
  - Quote module completeness only where it reaches 1, or give the KOs.
- **One MAG per species.** Cluster-level counts can be higher (TMA candidates 27 → 31, MttB family 136 → 145). JAWLGC01's MePn call rests on one short MaxBin2 contig.
- **Thresholds are release-specific.** Every call here uses KOfam 2026-10-26; record it in Methods.

## 7. Files

| path | what |
|---|---|
| `results/merged/`, `results/distill/`, `results/distill_gene_names/`, `results/env/`, `results/logs/`, `results/state/`, `results/distill_forms/` | copied from poplar (section 1) |
| `results/MD5SUMS_poplar.txt`, `results/MD5SUMS_local_gz.txt` | md5s computed on poplar (the 596 copied files; `merged/annotations.tsv` refers to the uncompressed file), and the md5 of the local `.gz` |
| `results/validation/` | `validation_report.txt`, `gene_counts_dram_vs_kofam.tsv`, `genes_unmatched_by_coordinates.tsv` (empty), `marker_genes_dram_vs_kofam.tsv.gz` (3,785 genes), `marker_gene_disagreements.tsv` (empty), `marker_ko_totals_dram_vs_kofam.tsv`, `mag_marker_cells_differing_post_rescue.tsv` (44 cells), `all_ko_best_ko_disagreements.tsv.gz` (50 genes) |
| **`results/species_558_methane_guilds.tsv`** | **the tidy per-species table for the manuscript.** One row per species cluster: best MAG, GTDB R232 lineage, CheckM, the MePn/betaine/TMA/methanogenesis gene counts, guild calls (`betaine_reducer_final` includes the Sec rescue), `k14083_caveat`, cluster-level counts, and the key product rows |
| `results/species_558_dram_summary.tsv` | the full per-species table (216 columns): all 65 sheet KO counts (`n_<gene>`), annotation depth, tRNA/rRNA, DRAM guild calls, all KOfam guild columns of `guilds_1324.tsv` (`kofam_*`) and `guilds_species_clusters_new.tsv` (`cluster_*`), every methanogenesis, C1, sulfur, nitrogen and CAZy product value (`product: …`), and the KOs behind each non-CAZy product call (`ids: …`) |
| `results/species_558_salterns_methane_long.tsv` | long form: species × 65 sheet KOs, count and DRAM gene IDs |
| `results/summary/` | `salterns_methane_by_ko.tsv`, `product_functions_of_interest.tsv`, `guilds_dram_vs_kofam_558.tsv`, `notable_product_calls.tsv`, `k14083_mttB_by_family.tsv` |
| `analysis/validate_vs_kofam.py`, `analysis/summarize_species.py` | the two scripts that write everything above (local, py_venv) |

**Re-run** from `annotation_runs/dram/`:
```bash
~/Documents/py_venv/bin/python -I analysis/validate_vs_kofam.py    # ~1 min
~/Documents/py_venv/bin/python -I analysis/summarize_species.py    # ~30 s
```

**Left on poplar** (`/scratch1/afreiburger/salterns_2026/dram/`), with nothing deleted:
- `dram_out/` (per-MAG outputs);
- `merged/genes.*`, `merged/scaffolds.fna` and `merged/genbank/`;
- `review/` (1.17 GB, deletable per `README.md`);
- `dram_db/` with the raw downloads and `pfam.mmsmsa.*`. The README notes these few hundred GB aren't needed after setup. Delete them only when you decide to.
