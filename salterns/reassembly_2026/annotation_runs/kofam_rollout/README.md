# KOfam (DRAM's KO layer) + guild calls on the 935 kept MAGs of the Arm C rollout — poplar run

*Track `kofam_rollout`, set up and launched 2026-11-19. It is a re-pointed copy of the reviewed `../kofam/` track, which ran the 1,324 first-pass MAGs (`../kofam/README.md`, `RESULTS.md`; `reassembly_2026/WORKFLOW.md` §7.3). The method is unchanged: Prodigal 2.6.3 `-p meta -g 11 -q` on contigs ≥ 1 kb; KOfam 2026-10-26 (all 27,981 profiles); whole-MAG protein chunks with one representative per identical protein; 40 × HMMER 3.4 `hmmsearch --cpu 1 --noali -T 14.7 --domT 14.7 -Z <total proteins>`; `postprocess.py` with the DRAM 1.5 best-KO rule; the Sec read-through rescue; the guild rules. The original `kofam/` track was only read, here and on poplar. Nothing was committed and nothing was written to KBase.*

## Status

**Search RUNNING on poplar since 2026-11-19 14:25:55Z (08:25:55 CST).** 40 chunks, 935 MAGs, `-Z 2732576`. Validation PASSED before the launch (see "Validation").

- **Expected finish: 2026-11-19 ≈ 23:30Z to 2026-11-20 ≈ 02:30Z (17:30–20:30 CST)**, i.e. ≈ 9–12 h per chunk. The slowest chunk sets `ALL_DONE`; see "Timing and the ETA". `progress.sh` prints a live estimate.
- **Driver:** `run_all.sh`, **PID 3038311**, started as in the `kofam` track with `nohup setsid bash scripts/run_all.sh > logs/run_all.nohup 2>&1 < /dev/null &`.
  - The driver itself idles; every hmmsearch and postprocess runs under `nice -n 10`.
  - It leads its own session (SID 3038311, PPID 1), so an ssh drop does not affect it.
  - When the last chunk ends it runs `postprocess.py combine` and touches `ALL_DONE`.
- **Chunk jobs:** `run_chunk.sh` PIDs 3038346… in `logs/chunk_pids.tsv`. The hmmsearch PIDs 3038393…3038859 are in `logs/hmmsearch_pids.tsv` (columns chunk, pid, ppid, nice, %cpu, rss_kb). Copies of both are in `poplar_copies/`.
- **Load when launched:** 1-min load 247 on 336 cores (limit set for this run: ~600), 2.35 TB RAM available. Each hmmsearch runs at nice 10, ≈ 109 % CPU (one worker thread plus the reader thread) and ≈ 44 MB RSS. There are 40 single-CPU jobs, the maximum allowed for this run.
- Each chunk postprocesses itself when its search ends (`results/chunks/chunk_XX.postprocessed`) and then gzips its raw tables.

### Status commands

```bash
# live per-chunk progress: profiles done, fraction of total profile length (sum LENG), elapsed, projected total h per chunk
ssh -o BatchMode=yes poplar 'bash /scratch1/afreiburger/salterns_2026/kofam_rollout/scripts/progress.sh'

# short check: chunks done / failed, driver alive, ALL_DONE, load
ssh -o BatchMode=yes poplar 'cd /scratch1/afreiburger/salterns_2026/kofam_rollout
  ls hmm/*.done 2>/dev/null | wc -l; ls hmm/*.failed results/chunks/*.postprocess_failed 2>/dev/null
  kill -0 3038311 2>/dev/null && echo "driver alive" || echo "driver gone"
  tail -3 logs/run_all.log; ls ALL_DONE 2>/dev/null; uptime'
```

- The ETA column of `progress.sh` is the **projected total runtime of the chunk** (elapsed / fraction of ΣLENG done), counted from the chunk's start at 14:25:55Z. It is not the time remaining.
- Profiles are searched in K-number order. The low K numbers are large oxidoreductase families with many hits, so the early estimate is pessimistic; in the `kofam` track the 16-minute forecast was ≈ 15 % too long.
- `progress.sh` greps the growing `hmm/chunk_XX.out` files over NFS, so it slows down near the end. Run one copy at a time.

**Do not delete or rebuild `kofam/db/kofam_all.hmm` until `ALL_DONE` exists.** This track's `db/kofam_all.hmm` is a symlink to it, the 40 hmmsearch processes read it for the whole run, and deleting it from another NFS client would make every chunk fail (ESTALE). The same applies to `bakta_test/kofam/chunks/`, which `sec_full_db.sh` needs in step 3c. **Do not edit `scripts/run_all.sh` or `scripts/run_chunk.sh` on poplar until `ALL_DONE` exists.** bash reads a running script lazily by byte offset, so an in-place edit changes what the waiting processes execute (`../kofam/README.md`, "Independent review"). `postprocess.py` is safe to edit, since each chunk starts it as a new process; keep the poplar and local copies identical.

**Restart after a failure or a reboot:** `cd /scratch1/afreiburger/salterns_2026/kofam_rollout; nohup setsid bash scripts/run_all.sh > logs/run_all.nohup 2>&1 < /dev/null &`. Chunks with `hmm/chunk_XX.done` are skipped. Delete a killed chunk's partial `hmm/chunk_XX.*` first, because it restarts from the beginning. Re-run one chunk's postprocessing with `…/bin/python3 scripts/postprocess.py chunk $PWD XX`, then `postprocess.py combine $PWD`.

## What changed versus the `kofam` track

| item | `kofam` (1,324 first-pass MAGs) | `kofam_rollout` (935 rollout MAGs) |
|---|---|---|
| MAG set | 1,324 FASTAs copied to `salterns_2026/mags/<name>.fna` | **935 FASTAs read in place** from `salterns_2026/rollout/mags_kept/<sample>__<bin>.fa`. 935 files = 935 rows of `rollout/results/mags_kept.tsv` (all `keep == yes`) = the GTDB-Tk R232 summaries (831 bac120 + 104 ar53). md5 of every file in `mags.md5` |
| MAG names → binner / sample | `name.split('.')[0]`, `name.split('__')[1]` (scheme `<binner>.bin.N.fasta__<unit>`) | **`parse_mag()` in `call_genes.py` and `postprocess.py`**, which handles both schemes. Rollout names are `<sample>__<binner>.bin.N[_sub]`; the old parse would give binner `R1_A_D1__concoct` and unit `concoct.bin.007`, and would silently switch off the short-MaxBin2 flag. Checked: all 1,324 first-pass names parse exactly as before; all 935 rollout names parse to the `sample`/`bin` columns of `mags_kept.tsv` |
| `run_prodigal.sh` | `ls mags/*.fna`, `ALL_DONE` at 1,324 | `ls rollout/mags_kept/*.fa`, `basename .fa`, **`ALL_DONE` at 935**, work dir `kofam_rollout` |
| `-Z` (total proteins) | 3,911,146 | **2,732,576** (`chunks/chunk_stats.txt` → `run_chunk.sh`; the validation controls are not in it) |
| chunk grouping | 40 chunks of whole species clusters (558 first-pass skani clusters), LPT on de-duplicated residues; 8 control bins forced into chunk 01 | 40 chunks of **whole species clusters from `comparison_rollout/mags_three_way.tsv`** (373 rollout clusters), LPT, no forced controls. **This deviates from the brief**; see "Chunking" |
| `make_chunks.py` | needs a `species_cluster` column | `species_cluster` optional (missing → one group per MAG); empty controls file allowed |
| `run_all.sh`, `run_chunk.sh`, `pool_gene_tables.sh` | `W=…/kofam` | `W=…/kofam_rollout`; `run_all.sh` takes the chunk count from `chunk_stats.txt` (40) |
| `progress.sh` | `W=…/kofam`, fixed `/tmp/afr_kofam_cumleng.tsv`, 40 | `W=…/kofam_rollout`, private `mktemp` file, chunk count from `chunk_stats.txt` |
| `sec_full_db.sh` | `W=kofam/sec_readthrough`, `-Z 3911146` hard-coded | `W=kofam_rollout/sec_readthrough`, `-Z` read from `kofam_rollout/chunks/chunk_stats.txt` (2,732,576) |
| `build_guilds.py` | `guilds_1324.tsv`, `mags_old_vs_new.tsv`, asserts 558 clusters, originals-310 block | **`guilds_935.tsv`**, metadata from `mag_manifest.tsv`, no 1,324/558 asserts, no originals block; per-cluster table `guilds_935_species_clusters.tsv` only if the cluster table exists (best MAG = its `rollout_species_rep`); `--res/--manifest/--clusters/--sec/--out` options for tests |
| HMM DB | built by `build_db.py` into `kofam/db/` | `db/kofam_all.hmm` = **symlink** to `kofam/db/kofam_all.hmm` (7.7 GB, verified, not rebuilt); `hmm_info.tsv`, `markers.hmm`, `db_check.txt` copied |
| new scripts | | `check_db.py` (DB vs release tarball), `run_validation.sh`, `compare_validation.py` |
| dropped | | `originals_markers.py`, `compare_vs_original.py`, `check_originals_vs_control.py`, `compare_controls.py`, `run_tests.sh`, `make_test_hmms.py` (first-pass/originals only) |
| unchanged | | `sec_rescue.py`, `postprocess.py` call rule, `build_db.py` (kept for the record, not run) |

## Inputs (md5)

| input | md5 / check |
|---|---|
| 935 MAG FASTAs `rollout/mags_kept/*.fa` | per file in `mags.md5` (md5 of that list `57de4ae3e0289f6a8e4554bf95e91ce7`) |
| `rollout/results/mags_kept.tsv` | `53a255c5f647c4bf4935e854f924c3d6` |
| `rollout/gtdbtk/gtdbtk.bac120.summary.tsv` / `gtdbtk.ar53.summary.tsv` | `4e4319e132024e8a7b62e960a07588fd` / `3215493dce44c772b2ad4f75776df65f` (equal to `mags_kept.tsv`'s `gtdb_classification` for all 935) |
| `rollout/DONE` | `2026-11-19T12:27:06Z samples=39/39 kept_mags=935 gtdbtk=done optional_failed=0` |
| `mag_manifest.tsv` (here, and on poplar) | `7b90bbb31f15b439ea5db48c39919e94`: name, sample, bin, binner, CheckM2, quality score, MIMAG label, GUNC, contigs, size, N50, GTDB R232 + ranks, poplar path, md5 |
| KOfam release `kofam_db/ko_list.gz` / `profiles.tar.gz` | `697af5b3b22a919306a1e60f99cec1a6` / `bc68ad2a52ed48812b5496d15c215a7d`: **`md5sum -c kofam_db/md5sums.txt` OK** (the documented values, WORKFLOW.md §4) |
| `db/kofam_all.hmm` (→ `kofam/db/kofam_all.hmm`) | **`6fcb36e48f6c0c3a8f34867a95731c6e`** (7,685,073,092 B, 27,981 profiles; new reference value, no md5 was documented before) |
| `db/markers.hmm` | `ebb733f9dc054ed8f53c6473a44d8dbe` (68 profiles) |
| `db/hmm_info.tsv` | `88c4b270f57f7bd183fb07e9069317e0` |
| `meta/marker_ko_table.tsv` | `05457c839419511ee9cbdbd19b1fdee7` (copy of `kofam/meta/`) |
| `chunk_manifest.tsv` (chunk grouping only) | `77e9625ad675e25af85da94e9cdd4d33`, from `comparison_rollout/mags_three_way.tsv` md5 `804b38f7cb948bf65ccb579dbbfac340` (mtime 2026-11-19 08:03:08 CST; `chunk_manifest.source.txt`) |

**DB check** (`scripts/check_db.py`, output `poplar_copies/check_db.txt`). The release md5s match. Streaming `profiles.tar.gz`:
- every one of the 27,981 `K*.hmm` members is byte-identical to its range in `kofam_all.hmm`;
- the member set equals the profile set, the NAMEs are in sorted order, and there are no stray bytes;
- `markers.hmm` equals the 68 marker ranges;
- `hmm_info.tsv` regenerated from `kofam_all.hmm` + `ko_list.gz` is byte-identical (27,981 rows, ΣLENG 16,466,034);
- the 39 `bakta_test/kofam/chunks/chunk_*.hmm` profile chunks used by `sec_full_db.sh` concatenate to `kofam_all.hmm` byte-for-byte.

So the DB was **not rebuilt**.

## Gene calls (done 2026-11-19 14:11:07Z → 14:19:04Z)

- `run_prodigal.sh` (40 × `xargs -P`, nice 10), Prodigal 2.6.3 (bioconda `h577a1d6_11`, EmilyKin env) `-p meta -g 11 -q`.
- **935/935 MAGs. Contigs: 326,815 in = 326,815 used; minimum length 1,000 bp, 0 below 1,000** (`poplar_copies/contig_len_check.txt`, `genes_summary.txt`), so the ≥ 1 kb filter drops nothing, as expected from `asm/contigs.fa`.
- **2,732,576 proteins, 814,431,445 aa**; by binner: MetaBAT2 1,051,951, MaxBin2 951,962, CONCOCT 728,663 genes (371 / 315 / 249 MAGs). The binner column of every gene table is a clean `metabat2`/`maxbin2`/`concoct`.
- Protein IDs `<sample>__<bin>|<contig>_<n>`; contig names (`NODE_…`) contain no `|`.

## Chunking: grouped by species cluster (deviation from the brief)

The brief said to LPT-balance whole MAGs, because species clusters for the rollout did not exist yet. By the time the chunks were built, `comparison_rollout/mags_three_way.tsv` existed (08:03 CST, 935 rollout rows, 373 clusters), and the numbers differed a lot:

| grouping | searched proteins | searched aa | saved | aa per chunk (min–max) | expected wall per chunk |
|---|---|---|---|---|---|
| whole MAGs (per the brief; `logs/chunks_per_mag_not_used/`, copy in `poplar_copies/chunks_per_mag_not_used/`) | 2,687,357 | 802,127,748 | 1.5 % | 20.29–20.45 M | ≈ 16–17 h |
| **whole species clusters (used)** | **1,848,636** | **561,075,163** | **31.1 %** | **13.93–14.26 M** | **≈ 10–12 h** |

- 33.1 % of all proteins are exact copies of a protein in another MAG: the same strain re-assembled from replicate samples. Only cluster grouping puts those copies into one chunk.
- Grouping is exact either way. hmmsearch scores do not depend on the target database, and E-values depend on it only through the fixed `-Z`; hits are copied to identical proteins within each chunk (`chunks/chunk_XX.members.tsv`).
- So the switch saves ≈ 240 CPU-h and ≈ 5 h of wall time on a shared host and changes no result.
- The cluster labels are only a search grouping. They appear in `chunks/chunk_assignment.tsv` and in no result file. If `mags_three_way.tsv` is regenerated, nothing in the search changes; `build_guilds.py` re-reads the table when it runs.
- **Checks after chunking:**
  - `proteins 2732576` (= `-Z`);
  - 935 MAGs, each in exactly one chunk, and no species cluster split across chunks;
  - 17–31 MAGs and 8–10 clusters per chunk;
  - for chunks 01, 17 and 40, searched representatives + copied members = assigned proteins (66,419; 55,435; 71,625).

  Files: `poplar_copies/chunk_stats.txt`, `chunk_assignment.tsv`, `make_chunks.log`.

## Validation (before the launch): PASS

**Design.** The 8 first-pass control bins of `kofam/meta/controls_chunk01.txt` (2 MaxBin2 *Roseovarius*, the *Methanolobus* `concoct.bin.101`, 5 *Albimonas*) were run as a separate chunk. They were not in the production genes dir, its chunks or its `-Z`.
- **Inputs.** md5 of the 8 `salterns_2026/mags/<name>.fna` = `kofam/mags.md5` (`validation/logs/controls_md5_check.txt`).
- **Pipeline.** Each step used this track's re-pointed scripts (`scripts/run_validation.sh genes|search`):
  1. `call_genes.py`;
  2. `make_chunks.py` (1 chunk, one group per MAG; 27,413 proteins, 20,684 searched);
  3. `hmmsearch` with the 68 marker profiles at the production `-Z 2732576`;
  4. `postprocess.py chunk`.
- **Reference.** The original track's `kofam/test/markers` run, the same 68 profiles on chunk 01 of the 1,324 run at `-Z 3911146`. Compared with `scripts/compare_validation.py` → `validation/validation_comparison.tsv`.

| compared | cells | differing |
|---|---|---|
| gene calls: `.faa`, gene table, `.gff`, `.done` of the 8 bins vs `kofam/genes/` | 32 files | **0** (byte-identical) |
| per-MAG marker/guild table (`mag_markers_guilds`), every column except `chunk` | 8 × 91 = 728 | **0** |
| … of which the 68 `n_<marker>` best-KO counts | 544 | **0** |
| `kofam/test/control_comparison.tsv`, column `poplar_2026_10` | 544 | **0** |
| marker hits (5,348 gene × KO rows), all columns except `full_evalue` | 101,612 | **0** |
| genes with a passing KO (150 genes): best KO, score, threshold, coverage, all passing KOs with scores | 1,650 | **0** |
| per-MAG gene tables (27,413 genes), all columns except the E-value | 411,195 | **0** |
| E-values | 5,348 + 150 | all differ, as expected: new/old ÷ (2,732,576 / 3,911,146) has median 1.0004 and range 0.95–1.07 (the tblout prints 2 significant digits) |

The expected calls are reproduced:
- both *Roseovarius* are `mepn_consumer`;
- *Methanolobus* is `methanogen` + `methylamine_user`;
- no *Albimonas* has *phnJ*;
- 1 bin is a one-sided TMA candidate.

**Further checks on rollout names (not in the brief):**
- **Production chunk 01 × 68 markers** (`run_validation.sh prodparse 01` → `validation/prod_chunk01_markers/`, 23 MAGs, 51 s):
  - `binner` and `unit` equal the manifest's binner and sample for all 23 MAGs;
  - the short-MaxBin2 flag equals (binner = maxbin2 and contig < 2,500 bp) on all 6,374 marker-hit rows; 516 rows and 4 MAGs are flagged. The old parse would flag none.
- **Production chunk 21 × 68 markers** (22 MAGs, Bacillota with *grd* genes). A dry run of the post-search steps on rollout names, in the scratchpad, with marker-only scores standing in for the full-DB search:
  - `sec_rescue.py build` found 81 trigger genes in 14 MAGs and 36 junctions;
  - `score` ran, and the `build_guilds.py` path with a Sec table ran: 9 betaine reducers, 5 of them with Sec-read-through *grdH*;
  - the pre-rescue guild calls equal `postprocess.py`'s in all 16 checked columns.

  These are code-path tests, not results.
- `build_guilds.py` on the 8 validation controls (`--clusters none --sec none`) and on chunk 01 with a 23-MAG subset of the cluster table: it ran, and the pre-rescue calls equal `postprocess.py`'s for every MAG.

## Timing and the ETA

- **Rate.** In the `kofam` track a chunk of 20.2 M searched aa took 16.6 h of wall time at 85–93 % of a core (load 520–640), i.e. 13.9 CPU-h of worker time. Scaled by residues, a 14.26 M-aa chunk needs ≈ 9.8 CPU-h.
- **Load.** It is lower now: load 250–330, and each job gets ≈ 109 % (worker + reader).
- **Live forecast.** `progress.sh` 4 min after the launch: 7.8–10.0 h per chunk. At 14:33Z (7 min, 1.35 % of ΣLENG done on average): **7.1–10.3 h**. The early rate is pessimistic (see above).
- **Expected finish ≈ 9–12 h after 14:25:55Z, i.e. 2026-11-19 ≈ 23:30Z to 2026-11-20 ≈ 02:30Z**, later if poplar's load rises again.
- **Disk.** ≈ 100 GB of uncompressed raw tables at the peak, before the per-chunk gzip; 8.5 TB free on `/scratch1`.

## After `ALL_DONE`: exact commands

All local commands run from `~/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/annotation_runs/kofam_rollout` with `PY=~/Documents/py_venv/bin/python`. All are light; nothing heavy runs locally.

**0. Check completion (poplar).**

```bash
ssh -o BatchMode=yes poplar 'cd /scratch1/afreiburger/salterns_2026/kofam_rollout
  ls ALL_DONE && cat results/combine_status.txt          # chunks_combined 40 ..., mags_covered 935
  ls hmm/*.failed results/chunks/*.postprocess_failed 2>/dev/null   # must print nothing
  tail -4 logs/run_all.log
  awk -F"\t" "/passing_pairs/{print \$4}" results/chunks/chunk_*.qc.txt | sort | uniq -c   # passing pairs without a domtblout row: expect mostly 0 (first pass: 12 of 61 M pairs, no marker)
  ls results/per_mag | wc -l                              # 935'
```

**1. Pool the per-MAG gene tables (poplar; needed by the Sec rescue).**

```bash
ssh -o BatchMode=yes poplar 'cd /scratch1/afreiburger/salterns_2026/kofam_rollout && \
  nohup setsid bash scripts/pool_gene_tables.sh > logs/pool_gene_tables.log 2>&1 < /dev/null &'
# done when genes_all.tsv.gz exists (written as .part, then renamed); zcat genes_all.tsv.gz | wc -l = 2,732,577
```

**2. Copy the results back and check them.**

```bash
P=poplar:/scratch1/afreiburger/salterns_2026/kofam_rollout
rsync -a -e "ssh -o BatchMode=yes" $P/results/ results/
rsync -a -e "ssh -o BatchMode=yes" $P/genes_all.tsv.gz $P/db/hmm_info.tsv $P/logs/run_all.log results/
ssh -o BatchMode=yes poplar 'cd /scratch1/afreiburger/salterns_2026/kofam_rollout && md5sum results/all_*.tsv* genes_all.tsv.gz' \
  | sed 's#results/##' | (cd results && md5sum -c -)
ls results/per_mag | wc -l    # 935
```

The raw hmmsearch tables `hmm/` (expect ≈ 5–6 GB gzipped), `genes/` (2.1 GB) and `chunks/` (0.8 GB) stay on poplar.

**3. Sec read-through rescue, from scratch on this set.** Do not reuse the 1,324-set `kofam/results/sec_readthrough/` files.

```bash
mkdir -p results/sec_readthrough mags_sec_trigger
# 3a. FASTAs of the MAGs with a Sec-KO trigger hit (grdH/grdB/grdF/grdA, score >= 14.7), from poplar, md5-checked
$PY -I -c "import pandas as pd; h=pd.read_csv('results/all_marker_hits.tsv.gz',sep='\t',usecols=['mag','ko']); \
m=sorted(h[h.ko.isin(['K21579','K10672','K21584','K10670'])].mag.unique()); \
open('results/sec_readthrough/trigger_mags.txt','w').write(''.join(x+'.fa\n' for x in m)); print(len(m), 'trigger MAGs')"
rsync -a -e "ssh -o BatchMode=yes" --files-from=results/sec_readthrough/trigger_mags.txt \
  poplar:/scratch1/afreiburger/salterns_2026/rollout/mags_kept/ mags_sec_trigger/
(cd mags_sec_trigger && awk 'NR==FNR{w[$1];next} ($2 in w)' ../results/sec_readthrough/trigger_mags.txt ../mags.md5 | md5sum -c --quiet -) && echo md5 OK
sed 's/\.fa$//' results/sec_readthrough/trigger_mags.txt | while read m; do printf "%s\t%s\n" "$m" "$PWD/mags_sec_trigger/$m.fa"; done \
  > results/sec_readthrough/fasta_map_rollout.tsv

# 3b. build the junctions (set name 'rollout')
$PY -I scripts/sec_rescue.py build rollout results/genes_all.tsv.gz results/all_marker_hits.tsv.gz \
  results/sec_readthrough/fasta_map_rollout.tsv results/all_genes_ko.tsv.gz results/sec_readthrough
cp results/sec_readthrough/rollout.joined.faa results/sec_readthrough/all.joined.faa
cp results/sec_readthrough/rollout.windows.fna results/sec_readthrough/all.windows.fna

# 3c. re-score the joined proteins against all 27,981 profiles + cmscan SECIS (poplar, a few minutes, -Z 2732576)
ssh -o BatchMode=yes poplar 'mkdir -p /scratch1/afreiburger/salterns_2026/kofam_rollout/sec_readthrough'
rsync -a -e "ssh -o BatchMode=yes" results/sec_readthrough/all.joined.faa results/sec_readthrough/all.windows.fna \
  poplar:/scratch1/afreiburger/salterns_2026/kofam_rollout/sec_readthrough/
ssh -o BatchMode=yes poplar 'cd /scratch1/afreiburger/salterns_2026/kofam_rollout && \
  nohup setsid bash scripts/sec_full_db.sh > sec_readthrough/sec_full_db.log 2>&1 < /dev/null &'
# wait for sec_readthrough/SEC_DONE ("... done 39 chunks Z=2732576"), then:
S=poplar:/scratch1/afreiburger/salterns_2026/kofam_rollout/sec_readthrough
rsync -a -e "ssh -o BatchMode=yes" $S/all.joined.full.tblout.gz $S/all.joined.full.domtblout.gz $S/all.windows.cmscan.tbl \
  $S/SEC_DONE $S/SEC_START results/sec_readthrough/

# 3d. marker cross-check search (local, seconds; same flags as the kofam track, -Z of this set)
(cd results/sec_readthrough && nice -n 10 ~/Documents/gtdbtk_venv/bin/hmmsearch -o /dev/null \
  --tblout all.joined.markers.tblout --domtblout all.joined.markers.domtblout --noali -T 14.7 --domT 14.7 \
  -Z 2732576 --cpu 1 ../../db/markers.hmm all.joined.faa)

# 3e. score + acceptance gate -> results/sec_readthrough/sec_readthrough_junctions_scored.tsv
$PY -I scripts/sec_rescue.py score results/sec_readthrough results/hmm_info.tsv \
  results/sec_readthrough/all.joined.markers.domtblout results/sec_readthrough/all.joined.full.tblout.gz \
  results/sec_readthrough/all.joined.full.domtblout.gz results/sec_readthrough/all.windows.cmscan.tbl
```

Check the junction counts against the first pass (`../kofam/RESULTS.md` §2: 545 junctions, 304 accepted, 45 new *grdH*). The rollout should have the same order of magnitude, and accepted junctions should fall only in *grd*-carrying Bacillota/Synergistota/Cloacimonadota-type genomes.

**4. Guild calls.**

```bash
$PY -I scripts/build_guilds.py
```

- It writes `results/guilds_935.tsv`, one row per MAG: sample, binner, CheckM2, GTDB R232, guild calls, `_anypass`, `_short_maxbin2_only`, Sec-rescue columns and marker counts.
- It writes `results/marker_only_vs_full_db_best_ko.tsv` and `results/sec_readthrough/replaced_prodigal_piece_calls_rollout.tsv`.
- It asserts that the pre-rescue calls equal `postprocess.py`'s (0 differing MAGs).

**5. Species-level aggregation.**
- `build_guilds.py` writes `results/guilds_935_species_clusters.tsv` whenever `reassembly_2026/comparison_rollout/mags_three_way.tsv` exists. That file exists already (08:03 CST version, md5 `804b38f7…`; 935 rollout rows, 373 clusters, columns `genome`, `set`, `species_cluster`, `rollout_species_rep`).
- It uses only the `set == rollout` rows (`genome = rollout|<sample>__<bin>`) and asserts that they are exactly the 935 manifest MAGs.
- The best MAG per cluster is the table's `rollout_species_rep`, which equals the highest CheckM2 completeness − 5 × contamination among the cluster's rollout MAGs (checked for all 373). The script asserts one per cluster.
- Columns per cluster: `cluster_sets`, `n_mags`, `<guild>_n_mags`, `best_mag` + its quality and taxonomy, `<guild>_best_mag`, `<guild>_any_mag`, `<guild>_best_mag_short_maxbin2_only`, `betaine_reducer_evidence_best_mag`, `phnGHIJKLM_n_distinct_best_mag`.
- **If the cluster table is regenerated, re-run step 4 only.** No re-search is needed. If it has not been finalised when step 4 runs, `--clusters none` skips the aggregation, and this step stays as the documented next step.
- A three-way comparison of guild calls against the first pass and the originals at cluster level (as `../kofam/scripts/compare_vs_original.py` did for 1,324 vs 310) is a separate next step and is not set up here.

## Caveats

- **E-values use `-Z 2,732,576`**, so they are 0.70 × the `kofam` track's for the same hit. The best KO is the lowest E-value within a gene, which does not depend on Z in exact arithmetic. But the tblout prints E-values to 2 significant digits, and ties are broken by the higher score, so a near-tie between two passing KOs could rank differently from the `kofam` track. This affected 0 of 150 genes in the validation. It matters only when comparing gene-level calls across the two tracks, never within this run.
- Inherited from the `kofam` track:
  - the `-T 14.7` full-score cut (3 domain-type KOs with thresholds < 15.7 could miss a hit inside a 1-bit window; none is a marker);
  - KOs with threshold `-` are never called;
  - the short-MaxBin2 flag is about MaxBin2 contigs < 2,500 bp.
- The Sec rescue and `build_guilds.py` were only smoke-tested on rollout names (marker-only scores, chunks 01 and 21). The real run in step 3 uses the full-DB scores.
- poplar is shared. Another user's load can stretch the ETA.

## Files

| here | poplar `/scratch1/afreiburger/salterns_2026/kofam_rollout/` |
|---|---|
| `scripts/` (identical copies on poplar; md5 in `poplar_copies/scripts_md5_poplar.txt` at setup) | `scripts/` |
| `mag_manifest.tsv`, `mags.md5`, `chunk_manifest.tsv`, `chunk_manifest.source.txt` | same files (`mags.md5` written there) |
| `meta/marker_ko_table.tsv`, `meta/controls_validation.txt` | `meta/` (+ `empty_controls.txt`) |
| `db/markers.hmm` | `db/kofam_all.hmm` → `kofam/db/kofam_all.hmm`; `db/hmm_info.tsv`, `markers.hmm`, `db_check.txt` |
| `poplar_copies/`: `check_db.txt`, `check_db_release_md5.txt`, `contig_len_check.txt`, `genes_summary.txt`, `chunk_stats.txt`, `chunk_assignment.tsv`, `make_chunks.log`, `chunk_pids.tsv`, `hmmsearch_pids.tsv`, `chunks_per_mag_not_used/`, `original_scripts.md5` (md5 of `kofam/scripts/` at copy time) | `logs/` (also `run_prodigal.log`, `chunk_XX.log`, `run_all.log`) |
| — | `genes/` (935 × faa/gff/genes.tsv/done + `ALL_DONE`), `chunks/` (40 × faa/members.tsv), `hmm/`, `results/` |
| `validation/`: `new/` (this run, incl. `all_*` from `postprocess.py combine`), `orig/` (kofam `test/markers`, 8 controls), `validation_comparison.tsv`, `Z.txt`, `controls.md5`, `logs/`, `prod_chunk01_markers/` | `validation/` (+ `prod_chunk21_markers/`) |
| `README.md` (this file) | `README.md` (mirror) |

## Results (2026-11-19 23:00 UTC; steps 0–4 above run as written)

- **Search:** `ALL_DONE` at 22:08 UTC. 40/40 chunks combined, covering 935 MAGs; no failed chunk or postprocess. "Passing pairs without a domtblout row" is 0 in 33 chunks, 1 in 5 and 2 in 2.
- **Results:** copied back to `results/`; md5 OK for `all_*` and `genes_all.tsv.gz` (2,732,577 lines).
- **Sec rescue** (`results/sec_readthrough/`), run from scratch:
  - 174 trigger MAGs gave 397 junctions.
  - 233 were accepted:
    - 199 confirm existing calls (grdA 102, grdB 94, grdH 3);
    - 34 are new read-through calls (grdH 30, grdB 4).
  - 164 were rejected.
  - Same order as the first pass (545 junctions / 304 accepted / 45 new grdH over 1,324 MAGs).
- **Guilds:** `build_guilds.py` gives `results/guilds_935.tsv` and `results/guilds_935_species_clusters.tsv`, using the 4-set `comparison_rollout/mags_three_way.tsv` regenerated at 21:00 UTC. Pre-rescue calls equal `postprocess.py` in every column.
  - Counts are rollout MAGs; species are given as best MAG / any MAG.

    | Guild | MAGs | Species |
    |---|---|---|
    | methanogen (*Methanolobus*) | 2 | 1 |
    | MePn consumer | 78 | 21 / 23 |
    | betaine reducer | 32 | 8 / 9 |
    | methylamine user | 2 | 1 |
    | TMA producer | 0 | 0 |
    | one-sided TMA candidates | 60 | 23 / 26 |
    | mttB family without mcrA | 228 | |

- **Sanity check against the first pass** (`results/guild_agreement_vs_first_pass.tsv`), over the 309 species clusters both sets contain:
  - methanogen, betaine reducer, methylamine user and TMA calls agree in every shared species;
  - MePn consumer differs in 3: JAWLGC01 is called by the first pass only; *Roseovarius* and *Haloplanus* by the rollout only.
  - No pipeline problem.
- **Final MAG set** (rollout 935 + kept co-assembly 520, 585 species): co-assembly MAGs use their first-pass KOfam calls (same DB and rule).
  - 1,417 of 1,455 MAGs have calls. The 38 without are exactly `coasm_final/kofam_missing.tsv`, none HQ.

    | Guild | MAGs | Species |
    |---|---|---|
    | MePn consumer | 108 | 39 |
    | betaine reducer | 52 | 18 |
    | methanogen | 4 | 1 |
    | methylamine user | 4 | 1 |
    | TMA producer | 0 | 0 |
