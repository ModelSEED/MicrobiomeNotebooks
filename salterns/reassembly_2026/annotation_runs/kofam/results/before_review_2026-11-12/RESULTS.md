# KOfam on the 1,324 re-assembly MAGs: results, Sec rescue, guilds, and comparison with the 310 originals

*Track `kofam`, results written 2026-11-12. Method and run details: `README.md`. Release: **KOfam 2026-10-26** (27,981 profiles); KO rule = DRAM 1.5 `sig_scores`; Prodigal 2.6.3 `-p meta`; contigs ≥1,000 bp. Nothing was committed.*

## Short answer

- **The run is complete and checks out.**
  - 40/40 chunks finished, covering **1,324/1,324 MAGs** and 3,911,146 genes. **1,578,156 genes (40.4 %) have a KO.**
  - The local guild code reproduces the poplar postprocessor exactly (0 differing MAGs in 16 columns).
- **Sec read-through rescue.** I implemented it for the four Grd selenoproteins: grdH K21579, grdB K10672, grdF K21584 and grdA K10670.
  - It rescues **all 35 bakta_test loci**, with joined scores identical to that track (715.7–823.1; max difference 0.0). K21579 is the best KO over all 27,981 profiles, and a SECIS_4 sits 3 nt after every TGA.
  - Across all 1,324 MAGs it adds **45 new grdH calls** (33 of the 35 bakta_test loci plus 12 more) and confirms 3 calls (the 2 BM714 loci and a PWPR01 locus). All 48 accepted grdH junctions have a SECIS.
  - It fixes the grdH→grdB mislabel and the grdA double count (grdA genes 288 → 175).
  - It creates no grdH call outside Bacillota and Cloacimonadota, and **no call of any kind in Rhodobacteraceae or haloarchaea**.
- **Guilds (all 1,324 MAGs / 558 species clusters):**

  | guild | MAGs | clusters |
  |---|---|---|
  | MePn consumers | 87 | 35 |
  | betaine reducers | **53** (49 without the rescue) | 16 |
  | methanogens | 3 | 1 |
  | methylamine users | 3 | 1 |
  | TMA producers | **0** | 0 |
  | one-sided TMA candidates | 62 | 31 |

  - 328 MAGs carry the non-pyrrolysine MttB family without *mcrA*, including every one of the 122 Rhodobacteraceae.
- **Against the originals (species-cluster level).** I ran KOfam 2026-10 on all 310 original MAGs with the same rules, so each flip can be attributed.
  - **MePn:** of the 22 RAST consumers, 13 are recovered with the same call and 6 sit in 5 clusters the re-assembly did not recover (5 *Roseovarius* MAGs and *Natrinema*).
    - 2 are lost through genome content: Roseovarius.2 (the new bins are 52–54 % complete with a partial *phn* locus and no *phnJ*) and Phaeodactylibacter.1 (a lone *phnJ*).
    - 1 is lost through the rule: **Albimonas.1**, whose 371-aa *phnJ* fragment scores 135 vs 194.6, with the same fragment in the new bins.
    - **19 new MePn species clusters** appear, 9 of them Rhodobacteraceae (4 *Roseovarius*), plus 3 gained in clusters the originals already had.
  - **Betaine:** of the 9 RAST reducers, 5 are recovered, 3 *Halanaerobium* sit in clusters not recovered, and 1 is lost through the rule: BM714.4, whose full-length GrdI scores 431 vs 506.
    - **11 new betaine-reducer clusters** appear, all Bacillota (M55B127 ×4, *Halanaerobium*, QPJE01, …) plus Cloacimonadota.
  - **Methanogens:** only *Methanolobus*, recovered in 3 of its 4 bins.
  - **TMA:** RAST's 4 "TMA producers" were matches to function strings. 3 of the 4 have no CutC or CntA under KOfam. Neither assembly has a complete CutC+CutD or CntA+CntB pair at the KOfam thresholds.

## 1. Results copied back and checked

- **Copied** with `rsync` from `poplar:/scratch1/afreiburger/salterns_2026/kofam/` to `results/` here (330 MB):
  - `all_genes_ko.tsv.gz`, `all_mag_ko_counts.tsv`, `all_mag_markers_guilds.tsv`, `all_marker_hits.tsv.gz` and `all_no_threshold_hits.tsv.gz`;
  - `per_mag/` (1,324 DRAM-like gene tables);
  - per-chunk `qc.txt`, `mag_markers_guilds.tsv` and `.postprocessed` markers, plus `combine_status.txt`;
  - the pooled gene → contig → binner table `genes_all.tsv.gz`;
  - `hmm_info.tsv` (thresholds and definitions) and `run_all.log`.
- **Left on poplar:**
  - the raw hmmsearch tables `hmm/` (8.0 GB);
  - `genes/` (3.0 GB of per-MAG faa/gff);
  - the chunk FASTAs;
  - the per-chunk copies of the combined tables.
- **md5** of the five combined tables and `genes_all.tsv.gz` is identical on poplar and here.

| check | expected (README) | found |
|---|---|---|
| chunks / MAGs | 40 / 1,324 | `combine_status.txt`: 40 chunks, `mags_covered 1324`; 1,324 `per_mag` files; MAG set = `mag_manifest.tsv` |
| run | forecast 17–23 h wall per chunk | 2026-11-11 17:36:39Z → last chunk 2026-11-12 10:10:11Z (**16.6 h**), 0 failed, combine done 10:11:23Z |
| genes | 3,911,146 (= `-Z`) | 3,911,146 = 2,826,115 searched representatives + 1,085,031 copied to identical proteins |
| genes with ≥1 KO | (42.3 % in chunk 01) | **1,578,156 (40.35 %)**; per chunk 34.8–46.3 % |
| tblout vs domtblout | passing pairs without a domain row = 0 in the partial test | 61,088,620 tblout pairs; 3,477,468 have no domain ≥14.7 (full score ≥14.7 from several weak domains). **12 of them pass a threshold** (18 genes after copying to identical proteins): K07112 ×13, K01154 ×3, K02004, K07497. **None is a marker**, so the only effect is that HMM coverage is reported as 0 for those calls |
| marker hits / no-threshold hits | — | 522,423 rows at score ≥14.7 / 1,537,058 hits at E ≤1e-5 |
| control bins (chunk-01 test) | see README | reproduced: both MaxBin2 *Roseovarius* `mepn_consumer` with 7/7 *phnGHIJKLM*; *Methanolobus* `concoct.bin.101` methanogen + methylamine user; all 5 *Albimonas* without *phnJ*; `maxbin2.Bin.031_sub R1_B_D1` MePn call lost without short MaxBin2 contigs |

## 2. Selenocysteine read-through rescue (`scripts/sec_rescue.py`, outputs in `results/sec_readthrough/`)

**What it does.** Sec KOs: all four Sec subunits in the marker table. The kofam track's poplar Prodigal genes are used unchanged (gene IDs are identical to the bakta_test pre-scan).
1. **Trigger:** any gene with a hit (score ≥14.7) to a Sec KO in `all_marker_hits`.
2. **Find the junction on the DNA**, not by joining proteins:
   - (a) the gene's own stop codon is TGA and the frame continues; or
   - (b) the first in-frame stop upstream of the gene's start is TGA. This is the *grdA* case, whose ~40-aa N-terminal piece Prodigal often does not call.

   The joined ORF runs from the upstream piece's Prodigal start, through the TGA (translated as U), to the next in-frame stop. If the upstream piece was not called, the start is the first ATG/GTG/TTG after the previous stop. The called partner gene is recorded as a pair when it is on the same strand and in the same frame, with no stop between. There is one record per junction, and each side must add ≥10 codons.
3. **Re-score** each joined protein:
   - against **all 27,981 profiles** on poplar (`scripts/sec_full_db.sh`, the kofam track's `-T 14.7 --domT 14.7 -Z 3911146`; about 2 min);
   - against the 68 markers locally, as a cross-check (the best KO agreed for 100 % of accepted junctions);
   - `cmscan` with Bakta's `ncRNA-regions` CMs (`--cut_tc -g --nohmmonly --rfam -Z 5`) on 603 nt around each TGA, for SECIS.
4. **Accept** a junction only if:
   - the best KO over all profiles (DRAM rule) is a Sec KO and passes its threshold;
   - an aligned domain for that KO spans the U, with ≥10 residues on each side;
   - the joined score beats every piece's score for that KO.

   A SECIS is recorded as evidence but not required.
5. **Status:**
   - `sec_readthrough`: accepted, and no piece already had that KO as its best KO. This is a new call, and it is the flag in the guild tables.
   - `confirmed`: accepted, and a piece already had the KO.
   - `rejected`: failed the acceptance gate.

   For accepted junctions, the pieces' Prodigal calls are replaced by the one joined call. Over all 356 accepted junctions, the pieces' previous best KOs (all-profile best for the re-assembly, `pieces_best_ko` in the scored table) were only K10670 ×291, K10672 ×184, K21579 ×5 and none ×207, so no other KO is dropped. The marker-gene calls removed from the guild counts are listed in `replaced_prodigal_piece_calls_new.tsv`.

**Funnel**

| | re-assembly (1,324) | originals (310) |
|---|---|---|
| trigger genes (any Sec-KO hit) | 991 in 245 MAGs | 190 in 57 MAGs |
| gene ends at TGA / upstream in-frame TGA | 533 / 538 | 109 / 106 |
| junctions (deduplicated; with a called partner pair) | 545 (413) | 112 (86) |
| best KO over all profiles is a Sec KO → of these, U spanned → of these, beats pieces | 324 → 304 → 304 | 55 → 52 → 52 |
| **accepted** (`sec_readthrough` + `confirmed`) | **304** (55 + 249) | **52** (7 + 45) |
| new calls by KO | **grdH 45**, grdB 9, grdA 1 | grdH 5, grdB 2 |
| SECIS (Rfam SECIS_4, starts 3 nt after the TGA's first base) | grdH 48/48, grdB 87/124, grdA 0/132; 2 of 241 rejected | grdH 7/7, grdB 13/21, grdA 0/24 |

**Checks.**
- **The 35 bakta_test loci.**
  - All 35 are found as junctions on the same Prodigal gene IDs, with the joined protein 394–440 aa and U at residue 324–354.
  - All 35 are accepted as grdH: 33 `sec_readthrough` and 2 `confirmed` (the two BM714 bins, which already passed at 635.3).
  - The joined grdH and grdB scores are identical to `bakta_test/results/interim/sec_readthrough_grdH_scores.tsv`, with a maximum difference of 0.0.
  - K21579 is the best KO over all 27,981 profiles in all 35, as the bakta_test `FOLLOWUP_DONE` run found.
  - A SECIS_4 is present in all 35.
- **13 further grdH junctions** are accepted (12 new calls, 1 confirmed), all with a SECIS (54.6–66.9). The bakta_test pre-scan looked only at the single best grdH gene of 184 bins with grdH ≥500, which is why it missed them:
  - 4 are a second Sec-grdH copy in bins already among the 35 (DATFJG01 `concoct.bin.001_sub`, `maxbin2.Bin.124_sub`; BM714 `maxbin2.Bin.114`, `.208`);
  - 1 is a Sec copy next to a full-length TAA copy (DATFJG01 `metabat2.bin.327`);
  - 3 are pre-scanned bins below the ≥500 filter (UBA8670 `maxbin2.Bin.008` 480.6 and `.013` 465.8; QPJE01 `metabat2.bin.356` 485.5);
  - 5 are in bins outside the pre-scanned taxa: CAYRVO01, PWPR01 (confirmed, piece 635.7), DAOWES01 ×2 (Cloacimonadota TCS61) and *Halonatronomonas*.
- Every one of the 48 accepted grdH junctions sits in a MAG that already had a passing *grdI*, *grdA*, *grdC* or *grdD*.
- **No calls in genomes without grd genes.**
  - Of 991 trigger genes, many are weak hits in Rhodobacteraceae, haloarchaea and other groups. Accepted junctions fall only in Bacillota (grdH 46 incl. confirmed), Synergistota (grdA/grdB), Cloacimonadota (grdH 2), SAR324 and Desulfobacterota (grdA/grdB only), and one Pseudomonadota MAG.
  - There are **none in Rhodobacteraceae or haloarchaea**.
  - The only accepted rescue in a MAG with no other passing non-Sec grd gene (grdI/A/C/D/E/G) is UBA7803 `concoct.bin.023 SF2_sediment` (Azotimanducaceae). It gets grdA 146 and grdB 292 (SECIS 62.4), not grdH, so it changes no guild.
- **Negative controls inside the run** (both sets together):
  - 23 junctions whose best KO is grdB were rejected because the U lies outside the HMM alignment:
    - 21 are a spurious upstream TGA 13–40 residues into an already complete GrdB, where the joined score changes by <1 bit;
    - 2 are mid-protein junctions whose joined score is lower than the piece's.
  - 278 junctions (221 re-assembly + 57 original) have a non-Sec best KO. 107 of these are D-proline reductase PrdB (K10794), another selenoprotein, so the read-through generalises, but it is not a guild KO and is not counted.
- **Desulfobacterota follow-up from bakta_test** (`metabat2.bin.289 R2_sediment`, the 345-aa grdB-like gene `k127_4725482_10`):
  - The HMM gate joins it (433 aa, U346, grdB 390.8 → 493.4, U spanned). It has no SECIS.
  - It already passed grdB, so it is `confirmed`.
  - The joined protein scores grdH 389.5 vs 624, so it stays glycine-reductase B and makes no betaine call.

**Effect on the calls (re-assembly):**

| quantity | without rescue | with rescue |
|---|---|---|
| grdH genes | 9 | 54 |
| MAGs with grdH | 9 | 49 |
| grdB genes (N-terminal GrdH pieces called *grdB* lose that call) | 190 | 154 |
| grdA genes (N and C pieces were counted twice in 113 of the 132 grdA loci) | 288 | 175 |
| **betaine reducers** | 49 | **53** |

The 4 MAGs gained are UBA8670 `metabat2.bin.079 R1_A_D2` and `maxbin2.Bin.014 R2_B_D1`, M55B127 `maxbin2.Bin.004 R1_C_D1` (the three without *grdI* named by bakta_test) and QPJE01 `metabat2.bin.356 R1_sediment`. All 4 are in species clusters that already had a betaine reducer, so the cluster count stays 16. On the originals the rescue turns 3 grdI-only reducers into grdH+grdI; the reducer count stays 7.

## 3. Guild calls (`scripts/build_guilds.py`)

**Files.**
- `results/guilds_1324.tsv`: one row per MAG, with sample, binner, species cluster, CheckM, GTDB R232 lineage, the guild calls, `*_anypass`, `*_short_maxbin2_only` flags, `betaine_reducer_evidence` (`grdH` / `grdI` / `grdH_sec_readthrough`), `betaine_reducer_without_sec_rescue`, `n_<marker>` (after the rescue), `n_<grd>_prodigal_only`, `n_<grd>_sec_readthrough`, `phn_accessory_counts`, `phnGHIJKLM_n_distinct/_n_contigs`, `pylBCDS_n_distinct`, and `marker_genes_on_short_maxbin2_contigs`.
- `results/guilds_species_clusters_new.tsv`: one row per species cluster (558), with the best MAG and, per guild, `_best_mag`, `_any_mag` and `_n_mags`.
  - The best MAG is the project's `best_of_cluster` flag (`annotation_choice/short_contig_fraction_per_mag.tsv`, highest completeness − 5 × contamination).
  - A max-quality/N50 tie-break would differ in 3 clusters, all exact quality ties with identical guild calls.

**Rules** (`annotation_choice/README.md` §6.3, on each gene's best KO; `_anypass` uses any passing KO and gives identical guild totals):
- methanogen = *mcrA*;
- MePn consumer = *phnJ*;
- betaine reducer = *grdH* (incl. Sec rescue) or *grdI*;
- methylamine user = *mttB/mtbB/mtmB* in an *mcrA* genome;
- TMA producer = *cutC*+*cutD* or *cntA*+*cntB*. *cutC*-only or *cntA*-only is `tma_candidate_one_sided`.

| guild | MAGs (of 1,324) | species clusters, best MAG / any MAG (of 558) | MAGs whose call needs a MaxBin2 contig <2,500 bp | clusters where every calling MAG needs one |
|---|---|---|---|---|
| methanogen | 3 | 1 / 1 | 0 | 0 |
| MePn consumer | 87 | 35 / 35 | 6 | **1** (JAWLGC01, see below) |
| betaine reducer | 53 (49 without the rescue) | 16 / 16 | 8 | 0 |
| methylamine user | 3 | 1 / 1 | 0 | 0 |
| TMA producer | 0 | 0 | — | — |
| TMA candidate (one-sided) | 62 (49 *cntA*-only, 13 *cutC*-only) | 27 / 31 | 1 | — |
| MttB family without *mcrA* (not a guild) | 328 | 136 / 145 | — | — |

**Notes.**
- **MePn.** 67 of the 87 consumers have their *phnGHIJKLM* genes split over more than one contig; 43 have all 7. Calling binners: MetaBAT2 47, CONCOCT 21, MaxBin2 19.
  - **7 clusters are heterogeneous** (*phnJ* in some MAGs only): 4 *Roseovarius*, *Lentibacter*, *Halorubrum* and JAWLGC01. The 14 non-calling MAGs are not less complete (median 78 % vs 74 % for calling MAGs). 10 of them carry 1–6 other *phn* genes, so most of this is a fragmented, partly binned locus rather than an absent operon. In all 35 MePn clusters, the 7 included, the best MAG carries *phnJ*. So the per-cluster count does not depend on using all 1,324 MAGs, but the per-MAG abundances do.
  - **JAWLGC01** (`maxbin2.Bin.441 R2_sediment`) gets its only *phnJ* (407) from a 1,493 bp MaxBin2 contig carrying only *phnJK*. Its two 98–99 % complete siblings and the original JAWLGC01.1 have no *phn* gene at all. **Treat this cluster's MePn call as a likely mis-binned contig.**
- **Methanogens.** The *Methanolobus* cluster has *mcrA* in 3 of 4 MAGs (983.3). `maxbin2.Bin.045 R1_C_D1` (61 % complete) has only a 28-aa *mcrA* fragment (33.3) on a 1,313 bp contig. It still has *mttB/mtbB/mtmB* (1/1/3), *pylS/pylD* and *mtbA*. By the strict rule it lands in `mttB_family_without_mcrA`, although it is a *Methanolobus*.
  - The new *Methanohalophilus* (`metabat2.bin.187_sub`, 55 %) has **no *mcr* hit at all**. It has *mtaB* ×2 (methanol), *pylS* and *pylB*, but no *mtt/mtb/mtm*. So it is a methanogen by taxonomy but not by the gene rule. Report it as an incomplete methanogen bin.
- **TMA.**
  - All 13 *cutC* genomes fail *cutD*: the best *cutD* is 330–346 vs 391.03 in BM714, *Desulfotignum* and JABDQA01, where *cutC* scores 1,036–1,086 vs 989.
  - All 49 *cntA* genomes fail *cntB*: the best are 386.5 and 377.7 (two *Roseovarius* bins) vs 424.6.
  - So "0 TMA producers" is a threshold statement about the activase/reductase subunit. The *cutD* and *cntB* thresholds look strict for these lineages; this is a candidate list, not an absence.

## 4. Original 310 vs re-assembly at species-cluster level (`scripts/compare_vs_original.py`, `results/comparison_vs_original/`)

**Three call sets.**
- **RAST:** `correlational/r232/guilds_r232.csv`, regular expressions on RAST strings. Totals after the join are 22 MePn, 9 betaine, 4 "TMA" and 1 methanogen. The 2020 manuscript sets are 27 MePn and 15 betaine.
- **KOfam-orig:** new here. `scripts/originals_markers.py` runs Prodigal plus the 68 marker profiles of KOfam 2026-10 on all 310 original MAGs (993,654 genes), with the same rules and Sec rescue: `results/originals_310/guilds_310_kofam.tsv`. Totals are 21 MePn, 7 betaine, 1 methanogen, 25 TMA candidates and 70 MttB family without *mcrA* (25 of them Rhodobacteraceae).
- **KOfam-new:** section 3.

The name join uses `key()` from `rerun_correlational_r232.py` and was asserted 310 ↔ 310.

**Validating KOfam-orig.**
- Against the earlier control test on the same 38 originals (2026-05 release, same local Prodigal; `scripts/check_originals_vs_control.py` → `results/originals_310/control_test_comparison.tsv`), **20 of 2,014 MAG × marker cells differ, and all come from the release change**:
  - *mttB* +1/+2 in 5 MAGs and *mtgB* −1/−2 in the same 5. These genes pass both KOs, and the 2026-10 *mttB* profile now has the lower E-value.
  - *grdB* and *grdE* +2/+3 (thresholds 412.6 → 278.0, …), *phnJ* +1 (246.6 → 194.6), *grdH*, *phnM*, *phnP* and *proX* ±1 (thresholds moved).
- Marker-only best KO (originals) vs best over all profiles (re-assembly):
  - on the 1,324 MAGs the two differ for **2 of 2,900 guild-KO genes**, both *cntA*;
  - this changes 2 MAGs' TMA-candidate flag and no guild call.
- The originals use the local Prodigal build. On the *Methanolobus* control it differs from poplar's by 2 short ORFs (README).

**Categories per (species cluster, guild).** Counted over the 604 clusters of both sets:

| category | meaning | MePn | betaine | methanogen | TMA (RAST vs KOfam candidate) |
|---|---|---|---|---|---|
| recovered_same | RAST member; ≥1 re-assembly MAG of the cluster calls it | 13 | 5 | 1 | 1 (J153.1) |
| not_recovered | RAST member; cluster absent from the re-assembly | 5 | 2 | 0 | 0 |
| lost_rule | RAST member; KOfam doesn't call it even in the original genome | 1 | 1 | 0 | 3 |
| lost_genome | RAST member; KOfam calls the original but no re-assembly MAG | 2 | 0 | 0 | 0 |
| gained_new_cluster | re-assembly call; no original MAG in the cluster | 19 | 11 | 0 | 10 |
| gained_rule | re-assembly call; KOfam also calls the original, RAST did not | 0 | 0 | 0 | 17 |
| gained_genome | re-assembly call; neither RAST nor KOfam in the original | 3 | 0 | 0 | 3 |
| kofam_orig_only | KOfam calls the original; RAST and the re-assembly don't | 0 | 0 | 0 | 6 |

### MePn consumers (RAST 22 MAGs in 21 clusters)

- **Recovered, same call (13 clusters):**
  - Roseovarius.1, .5, .8, .9 and .10;
  - Marivita.2 and .3;
  - Roseibacterium.1, JAVJPJ01.1 and Saliniramus.1;
  - Haloquadratum.4 (5/5 bins);
  - Desulfotignum.1 and .2.

  Scores in the original and new bins agree within about 1 point, except Marivita.3 (534 → 492). In 4 *Roseovarius* clusters only some bins call it (2/6, 3/4, 4/6, 3/5).
- **Not recovered (5 clusters, 6 MAGs):** Roseovarius.4, .6, .7 + .11, .12 and Natrinema.1. These are the original-only clusters in `comparison/README.md`.
- **Lost, genome content (2):**
  - **Roseovarius.2** (original 76 % complete, *phnJ* 531.5). The 2 new bins are 52–54 % complete, carry *phnG/H/I/K/L/M* hits but have **no *phnJ* hit at all**.
  - **Phaeodactylibacter.1** has a lone *phnJ* (322.5) with no other *phn* gene. Its 3 MaxBin2 bins have no *phnJ* hit, only weak *phnK/L/M*.
- **Lost, rule (1): Albimonas.1.**
  - RAST called the broken locus "C-P lyase". KOfam scores the 371-aa *phnJ* fragment 135.5 vs 194.6.
  - In the 5 new *Albimonas* bins, *phnJ* has no hit in 2; the others score 133.1 (the same 371-aa fragment, 36 % HMM coverage), 65.7 and 26.4.
  - **The rule and the genome agree: no assembly has a callable *phnJ*.**
- **Gained, genome (3):**
  - **CAZRZN01.1**: the original has 4 *phn* genes but no *phnJ* hit; both new bins have *phnJ* (534.8). It is in the 2020 "mepe" set of 27.
  - **Natronomonas.2**: no *phnJ* in the original; 2/2 new bins at 460.1.
  - **JAWLGC01.1**: one new MaxBin2 bin with *phnJ* on a short contig, flagged above.
- **Gained, new species clusters (19):**
  - 9 Rhodobacteraceae: 4 *Roseovarius*, *Lentibacter*, *Maritimibacter*, *Thermohalobaculum*, DAINOF01 and JBKDZR01;
  - 5 other Pseudomonadota: JACKII01, CAJQQK01, JAAAOG01, CAMCWB01 and WSZY01;
  - 1 UBA8366;
  - *Halorubrum* (haloarchaeon), RI-231 (Acidobacteriota), DATHSN01 (Actinomycetota) and UBA8670 (Bacillota).
- **Not RAST members:** the 5 other MAGs of the 2020 set of 27 have no *phnJ* under KOfam either: DAINOF01.1, Gaopeijia.1, Halorubrum.3, Thiolinea.1 and Albimonas.1 (1–4 *phn* genes).

### Betaine reducers (RAST 9 MAGs in 8 clusters)

- **Recovered, same call (5 clusters):** M55B127.1 (12/13 bins), M55B127.2 (7/7), DATFJG01.1 (8/8), BM714.2 (2/2) and UBA8670.2 (11/11).
  - In M55B127.1, M55B127.2 and UBA8670.2 the original's *grdH* is the Sec-truncated piece (596–606), so KOfam's original call rests on *grdI* plus the Sec rescue.
- **Not recovered (2 clusters):** Halanaerobium.1, and Halanaerobium.2 + .3.
  - Halanaerobium.3 (76 % complete) is not called by KOfam even in the original: *grdI* is a 225-aa fragment at 393 and *grdH* a 40-aa fragment, with no Sec junction for *grdH*.
  - The re-assembly's Halanaerobiaceae (*Halanaerobium* `metabat2.bin.206`, QPJE01) form **2 different species clusters**, both new betaine reducers.
- **Lost, rule (1): BM714.4** (25 % contamination; its cluster also holds BM714.1, a 2020-set member).
  - The original has a **full-length 428-aa GrdI scoring 431 vs 505.9** and *grdH* fragments of at most 361. The Sec-joined *grdH* reaches 595.8, below 624.
  - The 6 new bins of the cluster give the same numbers (*grdI* 427.7, Sec-joined *grdH* 595.8). RAST's "betaine reductase component B" string matched these genes.
  - A divergent GrdI just under the threshold in every assembly is a rule question (lower the threshold for this lineage or not), not missing data.
- **Gained, new species clusters (11):** M55B127 ×4, CAYRVO01 (JADQBU01), PWPR01, *Halanaerobium*, QPJE01, UBA8670, *Halonatronomonas* (Halarsenatibacteraceae) and DAOWES01 (Cloacimonadota TCS61).
  - All 11 have a passing *grdI* (594–848). 9 also have a *grdH*: 8 Sec-joined (716–817) and 1 full-length (781, UBA8670).
- **Not RAST members:** the 6 other members of the 2020 set of 15 (BM714.1, BM714.3, CAZFJG01.1/.2, DATFIM01.1 and UBA8670.1) have no *grdH* or *grdI* under KOfam either. They were model-curated in 2020, not gene-called.

### Methanogens and TMA

- **Methanogens:** Methanolobus.1 is recovered (3/4 new bins; *mcrA* 983.3 in both assemblies). There is no other *mcrA* genome in either set (the *Methanohalophilus* bin is discussed in section 3).
- **TMA, RAST's 4 candidates:**
  - **J153.1** is recovered as a KOfam candidate (*cntA* 397, *cntB* 348 vs 424.6).
  - **UBA11579.1, T1Sed10-7.1 and CAZFJG01.2** are `lost_rule`. RAST matched genes "clustered with choline TMA-lyase", and KOfam finds no *cutC* or *cntA* above 30 in either assembly.
- **TMA, KOfam candidates RAST missed:** 17 clusters where KOfam calls the original too, including:
  - Roseovarius.2, .3, .5 and .8 (CntA, which RAST calls "Rieske 2Fe-2S domain protein");
  - BM714.1/.3/.4 (*cutC* 1,077–1,086);
  - Desulfotignum.1 (*cutC* 1,037).

  There are also 10 new clusters, and 3 genome gains: Roseovarius.1, Roseovarius.9 and Marivita.2, where the original's best *cntA* hit scores 196–206 and the new bins score 455–459.

### Guild counts per taxon of interest (`taxa_of_interest.tsv`)

In each guild cell: RAST original MAGs / KOfam original MAGs / re-assembly MAGs (re-assembly species clusters).

| taxon | orig MAGs / clusters | new MAGs / clusters | clusters in both | MePn | betaine | methanogen | TMA candidate | MttB family w/o *mcrA* (new MAGs) |
|:--|:--|:--|--:|:--|:--|:--|:--|--:|
| Rhodobacteraceae | 25 / 24 | 122 / 38 | 19 | 16 / 15 / 58 (19) | 0 / 0 / 0 (0) | 0 / 0 / 0 (0) | 0 / 9 / 23 (10) | 122 |
| *Roseovarius* | 12 / 11 | 41 / 12 | 7 | 11 / 11 / 24 (9) | 0 / 0 / 0 (0) | 0 / 0 / 0 (0) | 0 / 7 / 17 (6) | 41 |
| *Marivita* | 3 / 3 | 17 / 3 | 3 | 2 / 2 / 16 (2) | 0 / 0 / 0 (0) | 0 / 0 / 0 (0) | 0 / 0 / 1 (1) | 17 |
| *Albimonas* | 1 / 1 | 5 / 1 | 1 | **1 / 0 / 0 (0)** | 0 / 0 / 0 (0) | 0 / 0 / 0 (0) | 0 / 0 / 0 (0) | 5 |
| *Halomonas* | 1 / 1 | 3 / 1 | 1 | 0 / 0 / 0 (0) | 0 / 0 / 0 (0) | 0 / 0 / 0 (0) | 0 / 0 / 0 (0) | 0 |
| *Methanolobus* | 1 / 1 | 4 / 1 | 1 | 0 / 0 / 0 (0) | 0 / 0 / 0 (0) | 1 / 1 / 3 (1) | 0 / 0 / 0 (0) | 1 (Bin.045) |
| *Methanohalophilus* | 0 / 0 | 1 / 1 | 0 | 0 | 0 | 0 / 0 / 0 (0), no *mcr* hit | 0 | 0 |
| Halanaerobiaceae | 3 / 2 | 3 / 2 | **0** | 0 / 0 / 0 (0) | 3 / 2 / 3 (2) | 0 | 0 | 0 |
| Desulfobacterota | 28 / 27 | 118 / 50 | 27 | 2 / 2 / 6 (3) (*Desulfotignum*, JAWLGC01) | 0 / 0 / 0 (0) | 0 | 0 / 1 / 4 (2) | 59 |
| haloarchaea (c. Halobacteria) | 23 / 23 | 75 / 28 | 17 | 2 / 2 / 8 (3) (*Haloquadratum* 5, *Natronomonas* 2, *Halorubrum* 1) | 0 | 0 | 0 | 0 |
| Bacillota (all) | 15 / 13 | 88 / 22 | 11 | 0 / 0 / 1 (1) | 9 / 7 / 51 (15) | 0 | 0 / 2 / 9 (2) | 45 |
| M55B127 (family) | 2 / 2 | 26 / 6 | 2 | 0 | 2 / 2 / 23 (6) | 0 | 0 | 26 |

- *Halomonas* has no guild call in either set. The second 2020 *Halomonas* is *Vreelandella* in R232 and has no re-assembly MAG.
- No betaine reducer is a Desulfobacterota. As in the bakta_test pre-scan, no Desulfobacterota bin passes *grdH* or *grdI*.

## Caveats

- **The originals were searched with the 68 marker profiles only**, so best KO is taken among the markers. On the new MAGs this differs from the all-profile best KO for 2 of 2,900 guild-KO genes (*cntA*). The originals also use the local Prodigal build.
- **Species-cluster comparisons inherit the re-assembly's fragmentation.** The median N50 is 6.4 kb vs 19.1 kb, and loci split over contigs (67 of 87 MePn consumers) can drop out of a bin. "Lost, genome content" is therefore a statement about the bins, not about the populations.
- The Sec rescue only checks the frame through a TGA. SECIS is recorded but not required, and *grdA* junctions never carry a SECIS_4 hit (0/156; their Sec is near the N-terminus). Their acceptance rests on the HMM gate and on the pieces already being *grdA*. *grdA* is support only and enters no guild rule.
- Thresholds are release-specific. Every number here is for KOfam 2026-10-26 and should be recorded in Methods.

## Files

| path | what |
|---|---|
| `results/` | copied from poplar (section 1) + the outputs below |
| `results/sec_readthrough/` | `{new,orig}.junctions.tsv` (all junctions), `sec_readthrough_junctions_scored.tsv` (scores, gate, SECIS, status: **the rescue table**), `{new,orig}.build_counts.txt`, `all.joined.faa`, `all.windows.fna`, the hmmsearch / cmscan outputs, `replaced_prodigal_piece_calls_new.tsv`, `fasta_map_{new,orig}.tsv` |
| `results/guilds_1324.tsv`, `results/guilds_species_clusters_new.tsv` | per-MAG and per-cluster guild calls (section 3) |
| `results/marker_only_vs_full_db_best_ko.tsv` | the 2 genes whose marker-only best KO differs from the all-profile one |
| `results/originals_310/` | `genes.tsv.gz`, `proteins.faa.gz`, `marker_hits.tsv.gz`, `markers.{tblout,domtblout}.gz`, `guilds_310_kofam.tsv`, `control_test_comparison.tsv`; `tmp/` (≈0.7 GB) holds the per-MAG Prodigal faa/gff, can be deleted (everything needed is in `genes.tsv.gz` / `proteins.faa.gz`), but `originals_markers.py` re-runs Prodigal without it |
| `results/comparison_vs_original/` | `clusters_guilds_old_vs_new.tsv` (604 clusters × guild: RAST / KOfam-orig / new n, category, best defining-gene scores), `original_guild_members.tsv`, `category_counts.tsv`, `taxa_of_interest.tsv` |
| `scripts/sec_rescue.py`, `sec_full_db.sh` | Sec rescue (build/score) and the poplar re-score (poplar copy: `kofam/sec_readthrough/`) |
| `scripts/originals_markers.py`, `check_originals_vs_control.py` | KOfam markers on the 310 originals; check against the 2026-05 control test |
| `scripts/build_guilds.py`, `compare_vs_original.py` | guild tables; comparison with the originals |
| `db/markers.hmm` | the 68 marker profiles of KOfam 2026-10 (copy of poplar `kofam/db/markers.hmm`) |

**Re-run** from `annotation_runs/kofam/`, with `~/Documents/py_venv/bin/python -I`:
1. `scripts/originals_markers.py`
2. `scripts/sec_rescue.py build new results/genes_all.tsv.gz results/all_marker_hits.tsv.gz results/sec_readthrough/fasta_map_new.tsv results/all_genes_ko.tsv.gz results/sec_readthrough`, and the same with `orig … results/originals_310/genes.tsv.gz results/originals_310/marker_hits.tsv.gz …fasta_map_orig.tsv - …`
3. `bash scripts/sec_full_db.sh` on poplar, plus a local `hmmsearch` against `db/markers.hmm`
4. `scripts/sec_rescue.py score results/sec_readthrough results/hmm_info.tsv <markers.domtblout> <full.tblout.gz> <full.domtblout.gz> <cmscan.tbl>`
5. `scripts/build_guilds.py`
6. `scripts/compare_vs_original.py`
