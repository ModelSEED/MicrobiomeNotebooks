# KOfam (DRAM's KO layer) + Sec read-through rescue + guild calls on the 38 kept co-assembly MAGs without KOfam calls — poplar run

*Track `kofam_coasm38`, set up and run 2026-11-19/20 (UTC). It is a re-pointed copy of the reviewed `../kofam_rollout/` track (935 Arm C rollout MAGs; `../kofam_rollout/README.md`), which is itself a copy of `../kofam/` (1,324 first-pass MAGs). The method is unchanged: Prodigal 2.6.3 `-p meta -g 11 -q` on contigs ≥ 1 kb; KOfam 2026-10-26 (all 27,981 profiles); whole-MAG protein chunks with one representative per identical protein; HMMER 3.4 `hmmsearch --cpu 1 --noali -T 14.7 --domT 14.7 -Z <total proteins>`; `postprocess.py` with the DRAM 1.5 best-KO rule; the Sec read-through rescue (from scratch on these 38); the guild rules (`build_guilds.py`). The `kofam_rollout` and `kofam` tracks were only read, here and on poplar. Nothing was committed and nothing was written to KBase.*

**Why these 38.** The final salterns MAG set is the 935 rollout MAGs plus the 520 kept MAGs of the four first-pass MEGAHIT site co-assemblies (R1/R2/R2A/SF2_sediment), re-scored with the rollout QC (`coasm_final/`). 482 of the 520 have KOfam calls from the first-pass track `../kofam/`; the other **38 had none** (`coasm_final/kofam_missing.tsv`: all MQ; 18 MaxBin2, 14 MetaBAT2, 6 CONCOCT; 34 bacteria, 4 archaea). This track gives them calls with the same method, so the caller can merge per-MAG tables into the final set.

## Status

**DONE 2026-11-20 ≈ 02:00Z** (`guilds_38.tsv` written 01:59:29Z). Search `ALL_DONE` at 01:55:32Z (16/16 chunks, 38 MAGs, no failed chunk or postprocess), validation PASSED, results copied back with md5, Sec rescue run from scratch, **`results/guilds_38.tsv`** written (38 rows, the 83 columns of `guilds_935.tsv`). See "Results".

- **`-Z 111,172`**, 16 production chunks + 1 separate validation chunk = 17 single-CPU hmmsearch jobs at nice 10 (≈ 110–113 % CPU each: worker + reader thread; max RSS 12–290 MB per chunk).
- Load before the launch: 1-min load 182 on 336 cores (`poplar_copies/logs/load_at_launch.txt`; 244 at setup, 00:48Z, with 2.36 TB RAM available); 126 at the end.
- Driver `run_all.sh` PID 494010 (own session, `nohup setsid`), chunk jobs and hmmsearch PIDs in `poplar_copies/logs/chunk_pids.tsv` / `hmmsearch_pids.tsv`; validation `run_validation.sh search` PID 473581 (own session). All exited normally.
- **First launch attempt failed instantly and was redone** (01:01:38Z → relaunch 01:02:11Z): `hmm/` did not exist. `run_all.sh` / `run_chunk.sh` create `results/chunks` and `logs` but not `hmm/` (in `kofam_rollout` it was created by hand at setup), so every hmmsearch exited with "Failed to open output file hmm/chunk_XX.out" before reading anything. Nothing was written; the logs of that attempt are in `poplar_copies/logs/launch1_failed_no_hmm_dir/`, the 0-chunk `combine` outputs were deleted, `hmm/` was created by hand, and `run_all.sh` was relaunched unchanged. **For a future copy: `mkdir -p hmm` before `run_all.sh`.**

## What changed versus `kofam_rollout`

| item | `kofam_rollout` (935 rollout MAGs) | `kofam_coasm38` (38 co-assembly MAGs) |
|---|---|---|
| MAG set | 935 FASTAs read in place from `rollout/mags_kept/*.fa` | **38 FASTAs read in place** from `coasm_qc/mags_kept/<mag_id>.fa`, listed in `meta/mags_38.txt` (= `mag_id` of `coasm_final/kofam_missing.tsv`). That dir holds all 520 kept co-assembly MAGs, so the scripts read the list, never a glob |
| `mag_manifest.tsv` | from `rollout/results/mags_kept.tsv` + GTDB-Tk | **same 22 columns**, from `coasm_final/coasm_bins_all.tsv` (`keep == yes`: CheckM2, `quality` → `mimag`, GUNC, contigs, size, N50), `coasm_final/gtdbtk/*.summary.tsv` (R232) and `mags.md5`; `quality_score` = completeness − 5 × contamination as in the rollout. Checked equal to `kofam_missing.tsv` for GTDB, quality, completeness, N50 |
| `run_prodigal.sh` | `ls rollout/mags_kept/*.fa`, `ALL_DONE` at 935, `xargs -P 40` | the 38 listed paths, **`ALL_DONE` at 38**, **`xargs -P 16`**, work dir `kofam_coasm38` |
| `-Z` (total proteins) | 2,732,576 | **111,172** (`chunks/chunk_stats.txt`; the 2 validation controls are not in it) |
| chunk grouping | 40 chunks of whole species clusters (LPT) | **16 chunks of whole MAGs** (LPT on residues; `make_chunks.py` without a `species_cluster` column = one group per MAG), as the brief asked. Identical proteins are rare across these 38 (1.4 %), so cluster grouping would save nothing |
| `run_all.sh`, `run_chunk.sh`, `progress.sh`, `pool_gene_tables.sh` | `W=…/kofam_rollout` | `W=…/kofam_coasm38`; comments only otherwise (chunk count and `-Z` come from `chunk_stats.txt`, as before) |
| `sec_full_db.sh` | `R=…/kofam_rollout` | `R=…/kofam_coasm38` (`-Z` 111,172 from this track's `chunk_stats.txt`) |
| Sec rescue marker cross-check (step 3d) | local `~/Documents/gtdbtk_venv/bin/hmmsearch` (HMMER 3.4) | **on poplar**, the EmilyKin env's `hmmsearch` (also HMMER 3.4), same flags and `-Z`; the local workstation is reserved for another session |
| `build_guilds.py` | `guilds_935.tsv`, set name `rollout` for both the Sec table and the cluster-table rows, per-cluster table `guilds_935_species_clusters.tsv` | **`guilds_38.tsv`**, `--n-mags 38`, Sec set name **`coasm38`** (junction ids `coasm38_secNNNN`) and a separate **`--cluster-set coassembly`**. `species_cluster` / `rollout_species_rep` are copied per MAG from the `coassembly|<mag>` rows of `comparison_rollout/mags_three_way.tsv`, so the columns equal `guilds_935.tsv`'s (83). **No species-cluster aggregation** (the caller merges per-MAG tables) |
| `run_validation.sh`, `compare_validation.py` | 8 first-pass controls × the 68 marker profiles vs the `kofam` marker test | **2 rollout MAGs × all 27,981 profiles** at this track's `-Z`, compared with `kofam_rollout`'s production output (see "Validation") |
| dropped | | `build_db.py`, `check_db.py` (the DB is not rebuilt; its md5 was re-verified instead) |
| byte-identical to `kofam_rollout` | | **`call_genes.py`, `make_chunks.py`, `postprocess.py`, `sec_rescue.py`** (md5 of the originals at copy time in `poplar_copies/kofam_rollout_scripts.md5`). Their docstrings still describe the rollout run; the code handles both naming schemes and needs no change for `<unit>__<binner>.bin.N[_sub]` co-assembly names |

**`build_guilds.py` regression test.** Run on scratch copies of `kofam_rollout`'s `all_genes_ko`, `all_marker_hits`, `all_mag_markers_guilds` and Sec table with `--manifest ../kofam_rollout/mag_manifest.tsv --n-mags 935 --set-name rollout --cluster-set rollout`, the modified script reproduces `kofam_rollout/results/guilds_935.tsv` exactly (935 × 83 cells, 0 differing, same column order), and `marker_only_vs_full_db_best_ko.tsv` and `replaced_prodigal_piece_calls_rollout.tsv` byte for byte.

## Inputs (md5)

| input | md5 / check |
|---|---|
| 38 MAG FASTAs `coasm_qc/mags_kept/<mag_id>.fa` | per file in `mags.md5` (md5 of that list `69e0ac32edcdddd4ebc7b6c8a7558bce`), 38 rows taken from `coasm_qc/mags_kept.md5` (`de1b3552c1eb680d256088ec3ba58ad3`; the local `coasm_final/mags_kept.md5` is identical). **`md5sum -c` 38/38 OK** on poplar (`poplar_copies/mags_md5_check.txt`) |
| contig lengths | **21,648 contigs, 105,677,914 bp, minimum 1,000 bp, 0 below 1,000** (`poplar_copies/contig_len_check.txt`), so the ≥ 1 kb gene-calling filter drops nothing (MaxBin2 bins go down to exactly 1,000 bp; MetaBAT2/CONCOCT ≥ 2,500) |
| `coasm_final/kofam_missing.tsv` | `991259c208293bff375b6efbb476de52` (38 rows) |
| `coasm_final/coasm_bins_all.tsv` | `1de20943f57eaa9ea61fad1bbe1c9bb8` (718 bins, 520 `keep == yes`); `coasm_final/DONE` = `2026-11-19T20:46:57Z units=4 bins=718 kept=520` |
| `coasm_final/gtdbtk/gtdbtk.bac120.summary.tsv` / `gtdbtk.ar53.summary.tsv` | `e50faf6d8e0f9ac26b19097d4c7cfee7` / `1e2c363fef56011a61e91dfe98fa720c` |
| `mag_manifest.tsv` (here and on poplar) | `6bb55b229171e56a77599b44f30f55b2` |
| `comparison_rollout/mags_three_way.tsv` (per-MAG `species_cluster` / `rollout_species_rep` only) | `685547e80727e8bcdd513366adaed107` (four-set version, 14:50 CST; 520 `coassembly` rows contain all 38; `rollout_species_rep` is False for all 520). The 935 rollout rows of this version give exactly the `species_cluster` / `rollout_species_rep` columns of `guilds_935.tsv` |
| `db/kofam_all.hmm` (→ `kofam/db/kofam_all.hmm`, symlink as in `kofam_rollout`) | **`6fcb36e48f6c0c3a8f34867a95731c6e`**, re-computed on poplar 2026-11-20 (`poplar_copies/kofam_all_hmm.md5`) = the `kofam_rollout` reference value; 7,685,073,092 B |
| `db/markers.hmm` / `db/hmm_info.tsv` / `meta/marker_ko_table.tsv` | `ebb733f9dc054ed8f53c6473a44d8dbe` / `88c4b270f57f7bd183fb07e9069317e0` / `05457c839419511ee9cbdbd19b1fdee7` (copies; = `kofam_rollout`) |
| validation controls `rollout/mags_kept/R2_A_H2O__metabat2.bin.034.fa` / `R1_C_D2__maxbin2.Bin.002_sub.fa` | `0a8ee3a0999f8e24df239e82313d714e` / `85f299154e3817be8fdb46a6d000b930` = `kofam_rollout/mags.md5`, `md5sum -c` OK |
| tools | Prodigal 2.6.3, HMMER 3.4 (Aug 2023), Infernal 1.1.5: EmilyKin env `544db00dd5c254ecfa7c6335967a046f_` on poplar, as in `kofam_rollout` |

## Gene calls (2026-11-20 00:56:24Z → 00:57:05Z)

- `run_prodigal.sh` (16 × `xargs -P`, nice 10), Prodigal 2.6.3 `-p meta -g 11 -q`. **38/38 MAGs**; contigs in = used = 21,648.
- **111,172 proteins** (31,663,451 residues in the `.faa`, counting the terminal `*` of complete genes). By binner: MaxBin2 58,948, MetaBAT2 38,537, CONCOCT 13,687 genes (18 / 14 / 6 MAGs). Binner column clean.
- Protein IDs `<unit>__<bin>|<contig>_<n>`; no contig name contains `|` (every header has exactly one).

## Chunking

`make_chunks.py` (unchanged) with a name-only manifest (`logs/chunk_manifest_names.tsv`) → one group per MAG, LPT on de-duplicated residues, 16 chunks, empty controls file (`poplar_copies/chunk_stats.txt`, `chunk_assignment.tsv`, `make_chunks.log`; `chunk_stats.txt` md5 `592b2aef67d64dacc73bbdff61d39faf`):

- `proteins 111172` (= `-Z`); 38 MAGs, each in exactly one chunk; 2–3 MAGs per chunk.
- 1.4 % of proteins are exact copies of a protein in another of the 38 MAGs; within-chunk de-duplication searches 111,078 (saves 0.08 %).
- Searched residues per chunk **1.86–2.16 M** (rollout: 13.9–14.3 M per chunk).

## Validation: PASS

**Design.** 2 rollout MAGs that `kofam_rollout` annotated, run as a separate control chunk (`validation/`; not in `genes/`, `chunks/`, `-Z`, `results/` or `guilds_38.tsv`):

- **`R2_A_H2O__metabat2.bin.034`**: *Natronomonas*, HQ, `mepn_consumer`, phnGHIJKLM 7/7 on one contig, MetaBAT2.
- **`R1_C_D2__maxbin2.Bin.002_sub`**: Bacillota *M55B127*, HQ, `betaine_reducer` (`grdH_sec_readthrough+grdI`, grdA/C/D), MaxBin2, so the MaxBin2 name parsing and the short-contig flag are exercised too.
- In `kofam_rollout` both sat in species-cluster chunks (39 and 12), where **1,949 of 2,970** and **1,897 of 2,258** of their proteins were not searched themselves but got the hits of an identical protein of another MAG (`kofam_rollout/chunks/chunk_39|12.members.tsv`). Here every protein is searched directly, so the comparison also checks the de-duplication copy path.

**Pipeline** (`scripts/run_validation.sh genes|search`): FASTA md5 = `kofam_rollout/mags.md5` → `call_genes.py` → `make_chunks.py` (1 chunk, 5,228 proteins, 5,220 searched) → `hmmsearch` × **all 27,981 profiles** at this track's **`-Z 111,172`**, same flags as `run_chunk.sh` → `postprocess.py chunk` + `combine`. 41.7 min (01:01:38Z → 01:43:30Z), 109 % CPU. Compared with `kofam_rollout`'s production output for the same 2 MAGs (`-Z 2,732,576`) by `scripts/compare_validation.py` → `validation/validation_comparison.tsv`, `validation/differing_cells.tsv` (empty).

| compared | cells | differing |
|---|---|---|
| gene calls: `.faa`, `.gff`, `.genes.tsv`, `.done` of both MAGs vs poplar `kofam_rollout/genes/` (md5) | 8 files | **0** (byte-identical) |
| per-MAG gene tables (5,228 genes), every column except `kofam_full_evalue` | 73,192 | **0** |
| … of which the best KO (`kofam_id`; 2,477 genes with a KO) | 5,228 | **0** |
| `genes_ko` (2,477 genes): best KO, score, threshold, score type, HMM coverage, n passing, all passing KOs with scores (in DRAM order) | 27,247 | **0** |
| `marker_hits` (1,027 gene × KO rows), every column except `full_evalue` | 19,513 | **0** |
| `mag_ko_counts` (1,988 MAG × KO rows: n best, n any-pass) | 3,976 | **0** |
| `mag_markers_guilds` (2 MAGs), every column except `chunk` | 182 | **0** |
| **total call cells** | **124,118** | **0** |
| E-values (`kofam_full_evalue`, `best_full_evalue`, marker `full_evalue`) | 2,477 / 2,477 / 1,027 | all non-zero ones differ, as expected: new/old ÷ (111,172 / 2,732,576 = 0.04068) has median 0.9999 (0.9994 for marker hits), range 0.945–1.057 (2-digit tblout rounding); 57 are 0 in both |

- The tie-break concern (E printed to 2 significant digits, scaled by 0.04 here) did not materialise: 0 best-KO or ordering differences in 2,477 genes.
- **Informational, Z-dependent by design:** `no_threshold_hits` keeps KOs without a threshold whose full E ≤ 1e-5, so the smaller `-Z` admits more rows: 2,626 here vs 2,409 in `kofam_rollout` for the 2 MAGs (all 2,409 present, 217 extra, 0 differing cells in common rows). No call depends on this table, but row counts of `all_no_threshold_hits.tsv.gz` are not comparable across tracks.

**Production checks.**
- In all 16 chunks, 0 passing (gene, KO) pairs lack a domtblout row (`results/chunks/chunk_XX.qc.txt`).
- `build_guilds.py`: the pre-rescue guild calls equal `postprocess.py`'s for all 38 MAGs in all 16 checked columns.
- Marker-only vs full-DB best KO: 52 guild-KO genes, 0 differing.

## Timing

| step | start → end (UTC) | wall |
|---|---|---|
| Prodigal (16 parallel) | 00:56:24 → 00:57:05 | 41 s |
| search, 16 chunks (1.86–2.16 M aa each) | 01:02:11 → 01:55:21 (last hmmsearch); `ALL_DONE` 01:55:32 | **41.1–53.2 min per chunk; 53.4 min to `ALL_DONE`** |
| validation chunk (2 MAGs, ≈ 1.5 M aa) | 01:01:38 → 01:43:30 | 41.7 min |
| pool gene tables, copy back, Sec rescue build, `sec_full_db.sh` (14 s), score, `build_guilds.py` | 01:56 → 01:59:29 | ≈ 4 min |

Per-chunk wall times are in `poplar_copies/logs/chunk_XX.log` (`/usr/bin/time -v`). The rate (21–26 min of wall per M searched aa at load 125–190) is faster than in `kofam_rollout` (≈ 33 min per M aa at load 250–330).

## After `ALL_DONE`: commands as run

All local commands run from `~/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/annotation_runs/kofam_coasm38` with `PY=~/Documents/py_venv/bin/python`. These are the `kofam_rollout` README steps 0–4, re-pointed. Step 5 (species aggregation) is skipped by design.

**0. Check completion (poplar).** `ALL_DONE`; `results/combine_status.txt` = `chunks_combined 16`, `mags_covered 38`; no `*.failed` / `*.postprocess_failed`; passing pairs without a domtblout row 0 in all 16 chunks; `results/per_mag` 38 files.

**1. Pool the gene tables (poplar).** `bash scripts/pool_gene_tables.sh` → `genes_all.tsv.gz`, 111,173 lines (seconds; run in the foreground).

**2. Copy back and check.**

```bash
P=poplar:/scratch1/afreiburger/salterns_2026/kofam_coasm38
rsync -a -e "ssh -o BatchMode=yes" $P/results/ results/
rsync -a -e "ssh -o BatchMode=yes" $P/genes_all.tsv.gz $P/db/hmm_info.tsv $P/logs/run_all.log results/
ssh -o BatchMode=yes poplar 'cd /scratch1/afreiburger/salterns_2026/kofam_coasm38 && md5sum results/all_*.tsv* results/per_mag/* genes_all.tsv.gz db/hmm_info.tsv' \
  | sed 's#results/##; s#db/##' | (cd results && md5sum -c --quiet -)      # OK: 5 all_* + 38 per_mag + genes_all + hmm_info
```

The raw hmmsearch tables `hmm/` (311 MB gzipped), `genes/` (79 MB) and `chunks/` (36 MB) stay on poplar.

**3. Sec read-through rescue, from scratch on these 38 (set name `coasm38`).**

```bash
mkdir -p results/sec_readthrough mags_sec_trigger
# 3a. MAGs with a Sec-KO trigger hit (grdH/grdB/grdF/grdA, score >= 14.7): 9 MAGs, 22 genes
$PY -I -c "import pandas as pd; h=pd.read_csv('results/all_marker_hits.tsv.gz',sep='\t',usecols=['mag','ko']); \
m=sorted(h[h.ko.isin(['K21579','K10672','K21584','K10670'])].mag.unique()); \
open('results/sec_readthrough/trigger_mags.txt','w').write(''.join(x+'.fa\n' for x in m)); print(len(m))"
rsync -a -e "ssh -o BatchMode=yes" --files-from=results/sec_readthrough/trigger_mags.txt \
  poplar:/scratch1/afreiburger/salterns_2026/coasm_qc/mags_kept/ mags_sec_trigger/
(cd mags_sec_trigger && awk 'NR==FNR{w[$1];next} ($2 in w)' ../results/sec_readthrough/trigger_mags.txt ../mags.md5 | md5sum -c --quiet -)   # OK (9)
sed 's/\.fa$//' results/sec_readthrough/trigger_mags.txt | while read m; do printf "%s\t%s\n" "$m" "$PWD/mags_sec_trigger/$m.fa"; done \
  > results/sec_readthrough/fasta_map_coasm38.tsv
# 3b. junctions
$PY -I scripts/sec_rescue.py build coasm38 results/genes_all.tsv.gz results/all_marker_hits.tsv.gz \
  results/sec_readthrough/fasta_map_coasm38.tsv results/all_genes_ko.tsv.gz results/sec_readthrough
cp results/sec_readthrough/coasm38.joined.faa results/sec_readthrough/all.joined.faa
cp results/sec_readthrough/coasm38.windows.fna results/sec_readthrough/all.windows.fna
# 3c. all 27,981 profiles + cmscan SECIS on poplar (-Z 111172), and 3d. the 68-marker cross-check, also on poplar
ssh -o BatchMode=yes poplar 'mkdir -p /scratch1/afreiburger/salterns_2026/kofam_coasm38/sec_readthrough'
rsync -a -e "ssh -o BatchMode=yes" results/sec_readthrough/all.joined.faa results/sec_readthrough/all.windows.fna \
  poplar:/scratch1/afreiburger/salterns_2026/kofam_coasm38/sec_readthrough/
ssh -o BatchMode=yes poplar 'cd /scratch1/afreiburger/salterns_2026/kofam_coasm38
  nohup setsid bash scripts/sec_full_db.sh > sec_readthrough/sec_full_db.log 2>&1 < /dev/null &
  cd sec_readthrough; A=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin
  nice -n 10 $A/hmmsearch -o /dev/null --tblout all.joined.markers.tblout --domtblout all.joined.markers.domtblout \
    --noali -T 14.7 --domT 14.7 -Z 111172 --cpu 1 ../db/markers.hmm all.joined.faa'
# SEC_DONE: "2026-11-19T19:59:00 done 39 chunks Z=111172" (poplar local time); copy back with md5 (OK):
S=poplar:/scratch1/afreiburger/salterns_2026/kofam_coasm38/sec_readthrough
rsync -a -e "ssh -o BatchMode=yes" $S/all.joined.full.tblout.gz $S/all.joined.full.domtblout.gz $S/all.windows.cmscan.tbl \
  $S/all.joined.markers.tblout $S/all.joined.markers.domtblout $S/SEC_DONE $S/SEC_START $S/sec_full_db.log results/sec_readthrough/
# 3e. score + acceptance gate -> results/sec_readthrough/sec_readthrough_junctions_scored.tsv
$PY -I scripts/sec_rescue.py score results/sec_readthrough results/hmm_info.tsv \
  results/sec_readthrough/all.joined.markers.domtblout results/sec_readthrough/all.joined.full.tblout.gz \
  results/sec_readthrough/all.joined.full.domtblout.gz results/sec_readthrough/all.windows.cmscan.tbl
```

If no MAG had had a trigger hit (or `build` had found 0 junctions; `sec_rescue.py build` would fail on an empty junction table), the rescue would have been skipped with `build_guilds.py --sec none`. That path prints "WARNING: no Sec rescue table (testing only)", which would then be the intended skip. Not needed here.

**4. Guild calls.** `$PY -I scripts/build_guilds.py` (defaults: `--set-name coasm38 --cluster-set coassembly --out guilds_38.tsv --n-mags 38`; log `results/build_guilds.log`) → `results/guilds_38.tsv`, `results/marker_only_vs_full_db_best_ko.tsv` (header only: 0 differing genes), `results/sec_readthrough/replaced_prodigal_piece_calls_coasm38.tsv`.

**5. Species-level aggregation: skipped** (the caller merges the per-MAG tables into the final set).

## Results

**Search.** 111,172 genes; 44,397 (39.9 %) have a passing KO (`all_genes_ko.tsv.gz`); 14,809 marker-hit rows. 411–5,634 genes per MAG.

**Sec read-through rescue** (`results/sec_readthrough/`, from scratch):

| | n |
|---|---|
| trigger MAGs (a Sec-KO hit ≥ 14.7) | **9** (22 trigger genes, 13 contigs) |
| junctions built | **14**; 9 with a called partner pair (`coasm38.build_counts.txt`: 14 trigger genes ending at TGA, 14 with an in-frame TGA upstream; a junction reached from both pieces counts once) |
| accepted | **5**, all `confirmed` (a piece already had that KO): grdA 3, grdB 2 |
| new read-through calls (`sec_readthrough`) | **0**, so **0 new grdH** |
| rejected | **9**: 5 have K10794 (D-proline reductase PrdB, a Sec protein of the same family, not a Grd KO) as best KO, like all 64 K10794 junctions in `kofam_rollout`; 2 have no significant KO; 2 keep a non-Sec best KO (K03366, K00053) whose alignment does not span the U |
| SECIS hits | 1 (SECIS_4 at `coasm38_sec0003`, grdB) |

- Smaller than the rollout (174 trigger MAGs, 397 junctions, 233 accepted, 30 new grdH), as expected for 38 MQ MAGs with 1 Bacillota.
- Four junctions have joined proteins identical to `kofam_rollout` junctions, with the same status (3 confirmed, 1 rejected).
- No guild call changes through the rescue (`betaine_reducer_without_sec_rescue` = `betaine_reducer`).

**Guilds** (`results/guilds_38.tsv`; best-KO calls after the rescue; `_anypass` counts are the same):

| guild | MAGs of 38 | which |
|---|---|---|
| methanogen | **0** | |
| `mepn_consumer` | **1** | `R1_sediment__metabat2.bin.275` *Roseovarius*: phnGHIJKLM 7/7 on 4 contigs |
| `betaine_reducer` | **1** | `R1_sediment__metabat2.bin.091_sub` Bacillota *M55B127*: evidence **grdI only** (grdA/C/D 2 of 3; grdB by Prodigal) |
| `methylamine_user` | **0** | |
| `tma_producer` | **0** | |
| `tma_candidate_one_sided` | **2** | `R1_sediment__metabat2.bin.275` (*Roseovarius*), `SF2_sediment__maxbin2.Bin.283_sub` (*Fodinibius*): cntA without cntB |
| `mttB_family_without_mcrA` | **7** | |

- No call depends only on genes on MaxBin2 contigs < 2,500 bp (all `_short_maxbin2_only` columns 0).
- **Consistency with the rollout.** 19 of the 38 share a species cluster (`mags_three_way.tsv`) with ≥ 1 rollout MAG.
  - Methanogen, MePn consumer, betaine reducer, methylamine user and TMA producer agree with the rollout cluster's any-MAG call in all 19.
  - The one-sided TMA candidate differs in 2: *Roseovarius* bin.275 is a candidate here, none of the 6 rollout MAGs of its cluster is; `SF2_sediment__maxbin2.Bin.307` is not, its 1 rollout cluster mate is.
- **Final-set arithmetic.** Adding these 38 to the 1,417 MAGs with calls in `kofam_rollout/README.md` "Final MAG set" gives 1,455/1,455 MAGs with calls.
  - By MAG counts: MePn consumer 108 + 1, betaine reducer 52 + 1, methanogen 4 + 0, methylamine user 4 + 0, TMA producer 0.
  - Species counts should not change for these guilds. The *Roseovarius* MAG's cluster (`rollout|R1_B_D1__metabat2.bin.017_sub`: 6 rollout, 2 co-assembly, 7 first-pass, 1 original MAG) and the *M55B127* MAG's cluster (`first_pass|maxbin2.Bin.023_sub.fasta__R2_sediment`: 1 rollout, 2 co-assembly, 2 first-pass MAGs) each already hold a rollout MAG with the same call. The caller's merge decides.

## Caveats and oddities

- **E-values use `-Z 111,172`**: ≈ 0.041 × `kofam_rollout`'s and ≈ 0.028 × `kofam`'s for the same hit. Calls do not depend on Z (validated: 0 of 124,118 call cells differ). But E-value columns and the row count of `all_no_threshold_hits.tsv.gz` (filter E ≤ 1e-5) are **not comparable across tracks**; this track's no-threshold table holds more weak hits.
- The *M55B127* betaine reducer has no grdH: the 5 junctions of this MAG are grdB/grdA (confirmed, 2 identical to rollout junctions) and 2 PrdB-type (rejected). The grdH-carrying contig of this species (accepted with joined scores 744–823 in rollout *M55B127* MAGs) is not in this bin. The call rests on grdI, as the rule allows.
- Inherited from `kofam_rollout` / `kofam`:
  - the `-T 14.7` full-score cut;
  - KOs with threshold `-` are never called;
  - the short-MaxBin2 flag covers MaxBin2 contigs < 2,500 bp. 18 of the 38 are MaxBin2 bins with contigs down to 1,000 bp; 2 of them have marker genes on such contigs, with no guild consequence.
- 17 concurrent jobs, not 16: 16 production chunks plus the separate control chunk.
- `compare_validation.py` had an indexing bug in its first run (it crashed before writing anything) and was fixed; only that local script changed.

## Files

| here | poplar `/scratch1/afreiburger/salterns_2026/kofam_coasm38/` |
|---|---|
| `scripts/` (identical copies on poplar; `poplar_copies/scripts_md5_local.txt`, `poplar_copies/logs/scripts_md5_launch.txt`) | `scripts/` |
| `mag_manifest.tsv`, `mags.md5`, `meta/mags_38.txt`, `meta/controls_validation.txt`, `meta/marker_ko_table.tsv`, `db/markers.hmm` | same files; `db/kofam_all.hmm` → `kofam/db/kofam_all.hmm`; `db/hmm_info.tsv`, `markers.hmm`, `db_check.txt`; `meta/empty_controls.txt` |
| `results/`: `guilds_38.tsv`, `all_*` (5 tables), `per_mag/` (38), `chunks/` (16 × per-chunk tables + `qc.txt`), `genes_all.tsv.gz`, `hmm_info.tsv`, `combine_status.txt`, `run_all.log`, `build_guilds.log`, `marker_only_vs_full_db_best_ko.tsv`, `sec_readthrough/` | `results/` (without the Sec / guild files), `genes_all.tsv.gz`, `sec_readthrough/` (full-DB re-score, cmscan, marker cross-check), `genes/` (38 × faa/gff/genes.tsv/done + `ALL_DONE`), `chunks/`, `hmm/` (raw tables, gzipped) |
| `mags_sec_trigger/` (9 trigger-MAG FASTAs, md5-checked) | — |
| `validation/`: `new/` (control-chunk results), `logs/` (md5 checks, gene-file comparison, `Z.txt`, `search.log`), `validation_comparison.tsv`, `differing_cells.tsv` (empty), `controls.md5`, `manifest.tsv` | `validation/` (+ `genes/`, `chunks/`, `hmm/`) |
| `poplar_copies/`: `logs/` (all poplar logs incl. `launch1_failed_no_hmm_dir/`), `chunk_stats.txt`, `chunk_assignment.tsv`, `make_chunks.log`, `mags_md5_check.txt`, `contig_len_check.txt`, `kofam_all_hmm.md5`, `genes_summary.txt`, `kofam_rollout_scripts.md5` | `logs/` |
| `README.md` (this file) | `README.md` (mirror) |
