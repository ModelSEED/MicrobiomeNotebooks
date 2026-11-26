# Salterns Arm C rollout: all 39 metagenomes, re-assembled and re-binned on poplar

Started 2026-11-13 on poplar. Nothing here is committed, and nothing is written to KBase.

**Approval.** Andrew approved this on 2026-11-13, after the binning test came back NO_GO only on SF2_C_D1's criterion d (`../binning_test/README.md`, "Status"):

- re-assemble and bin all 39 salterns metagenomes on poplar with the Arm C recipe, using the compute poplar has free;
- drop every MAG that fails GUNC or exceeds the CheckM2 contamination limit before it replaces a 2026 MAG.

The EmilyKin and GROW tests stay deferred.

**Status:** see "Launch record" at the end.

**R2A migration (2026-11-16).** 7 of the 9 R2A samples now run on frankensteinsmonster with this same driver and are imported back here before GTDB-Tk; R2A_A_H2O and R2A_C_H2O stay on poplar. See "R2A migrated to frankensteinsmonster" near the end.

**R2 migration (2026-11-17).** 8 of the 12 R2 samples (those ≤ 9 Gbp whose assembly had not started here) moved the same way. See "R2 migrated to frankensteinsmonster" near the end.

## What runs

One driver (`scripts/driver.py`, stdlib Python) runs every step as a detached bash task with a done marker. It is the reviewed binning-test driver (`../binning_test/scripts/driver.py`) with more task kinds, a load- and memory-aware scheduler, per-attempt SPAdes parameters, cleanup of intermediates, an append-only ledger and an end marker.

Per sample (all 39, per-sample assemblies only):

1. **Reads.** Download from NCBI's AWS open-data mirror, then `vdb-validate`. Pair and base counts must equal `manifest.tsv`, and every separator line must be a bare `+` after the rewrite. Then BBDuk (adapters, then artifacts and phiX) gives `clean_*`, and a quality trim gives `map_*` for coverage mapping.
2. **Arm C assembly.** `bbmerge-auto` (`rem extend2=50 k=62 ecct`), Tadpole correction (two passes), quality trim, then SPAdes 3.15.3 in **default mode** (`--only-assembler`, merged + unmerged reads).
3. **Binning, exactly as the binning test.**
   - bowtie2 maps of the site's coverage group (R1 9, R2 12, R2A 9, SF2 9 libraries) to the scaffolds ≥ 1 kb.
   - jgi depth tables, then MetaBAT2, MaxBin2 and CONCOCT, then DAS Tool 1.1.7.
4. **QC.** CheckM2 (primary), CheckM1 `--reduced_tree` (informational), GUNC and N50/size.
   - **KEEP = CheckM2 completeness ≥ 50 AND contamination < 10 AND GUNC pass.**
   - Kept MAGs are copied to `mags_kept/<sample>__<bin>.fa`.
5. **2026 side of the comparison (background).** CheckM2 + GUNC on all 2026 DAS Tool bins of the same unit.

At the end: GTDB-Tk 2.8.0 + GTDB R232 `classify_wf` on all kept MAGs, then the report, then `DONE` (or `FAILED` with reasons).

**Not done here.**

- **The 2026 site co-assemblies (R1/R2/R2A/SF2 sediment, MEGAHIT) are NOT redone.** Their 2026 MAGs remain as they are. Only the 39 per-sample assemblies are replaced.
- No KBase writes, and no DRAM (Andrew runs DRAM2 on h100).

### Reused, verified and not redone

| What | From | Check made by the driver |
|---|---|---|
| R1_A_D1 and SF2_C_D1 Arm C assemblies | `tadpole_pilot/<s>/armC/spades/` | `.done/C1_spades_default` exists; scaffolds ≥ 1 kb re-filtered here are byte-identical (`cmp`) to the binning test's `asm/contigs.fa` |
| Their binning and QC (NEW bins, CheckM2, CheckM1, GUNC; OLD CheckM2, GUNC) | `binning_test/` | the binning-test `.done` markers of `dastool`, `checkm2/gunc/checkm1:new`, `checkm2/gunc:old`; one QC row per bin (`rolib.py verify-import`: R1_A_D1 71 bins, SF2_C_D1 11) |
| BBDuk-clean reads of the 18 R1 + SF2 libraries (`clean_*` for assembly, `map_*` for mapping) | `binning_test/reads/<lib>/` (symlinked; the two pilot libraries resolve to `tadpole_pilot/<s>/qc/`) | `binning_test/.done/reads:<lib>` exists; the files are non-empty |

So only 21 libraries (R2 12, R2A 9) are downloaded. Nothing under `binning_test/` or `tadpole_pilot/` is modified or deleted.

## Exact commands

Each task is a bash file the driver writes once per attempt (`scripts/tasks/<task>.<attempt>.sh`, never edited afterwards). It runs with `set -euo pipefail`, `OMP/MKL/OPENBLAS_NUM_THREADS` = the task's threads, a short local `TMPDIR=/tmp/ro_<hash>` and a scratch temp `SCRTMP=$ROOT/tmp/<task>.<attempt>`. Both are removed on exit, or by the driver if the attempt was killed.

Every task inherits `nice 10` from the driver, which sets its own niceness to 10 at start. `$ROOT` = `/scratch1/afreiburger/salterns_2026/rollout`. Real-run thread counts are shown.

```
# dl:<lib>   (21 libraries; sra-tools 3.4.1)
curl -sS --retry 5 --retry-delay 20 -C - -o <srr>.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/<srr>/<srr>
vdb-validate <srr>.sra
fasterq-dump --split-files -e 6 -t fqd_tmp -O . <srr>.sra
awk 'NR%4==3{print "+" > out; next} {print > out} NR%4==2{n++; b+=length($0)} ...' > raw_<m>.fastq.tmp
#   pairs of both mates == manifest spots, total bases == manifest bases, 4 x spots lines per mate,
#   0 non-bare separator lines left; only then mv raw_<m>.fastq.tmp raw_<m>.fastq
# reads:<lib>  (BBTools 39.81, t=16; the 18 R1/SF2 libraries: symlinks to the binning test's files instead)
bbduk.sh -Xms2g -Xmx16g eoom t=16 in1=raw_1 in2=raw_2 out1=trim_1 out2=trim_2 ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo ftm=5 maq=5 minlen=45 mlf=0.333 maxns=0 ordered
bbduk.sh -Xms2g -Xmx16g eoom t=16 in1=trim_1 in2=trim_2 out1=clean_1 out2=clean_2 ref=artifacts,phix k=31 hdist=1 ordered
bbduk.sh -Xms2g -Xmx16g eoom t=16 in1=clean_1 in2=clean_2 out1=map_1 out2=map_2 qtrim=r trimq=10 minlen=45 ordered
#   chained Input/Result counts asserted; raw FASTQ deleted
# merge:<s> / tadpole:<s> / qtrimasm:<s>   (Arm C read preparation; heap -Xmx150g if the sample has <= 9 Gbp, else 200g;
#   250g on a retry after a Java OutOfMemoryError; -Xms8g so the heap is not reserved up front; eoom = exit on OOM)
bbmerge-auto.sh -Xms8g -Xmx<heap> eoom t=16 in1=clean_1 in2=clean_2 out=merged outu1=unmerged_1 outu2=unmerged_2 rem k=62 extend2=50 ecct prefilter=2 ihist=ihist.txt ordered
tadpole.sh -Xms8g -Xmx<heap> eoom t=16 mode=correct k=62 prefilter=2 in=unmerged_1 in2=unmerged_2 extra=merged out=ecc_u_1 out2=ecc_u_2 ordered
tadpole.sh -Xms8g -Xmx<heap> eoom t=16 mode=correct k=62 prefilter=2 in=merged extra=unmerged_1,unmerged_2 out=ecc_merged ordered
bbduk.sh -Xms2g -Xmx16g eoom t=16 in1=ecc_u_1 in2=ecc_u_2 out1=reads_1 out2=reads_2 qtrim=r trimq=10 minlen=45 ordered
#   merged/unmerged deleted after Tadpole; ecc_u_* after the trim
# asm:<s>   (SPAdes 3.15.3, python 3.9)
spades.py --only-assembler -k 33,55,77,99,127 -t 24 -m <mem> --phred-offset 33 -1 reads_1 -2 reads_2 --merged ecc_merged -o asm/spades
#   attempt 1: default mode, -m 400.  attempt 2: default mode, -m 800 if attempt 1 ran out of memory, else -m 400.
#   attempts 3-4 (default mode failed twice): the Arm A recipe, the same command with --meta, into asm/spades_meta,
#   marked mode=meta_fallback in asm/asm_meta.json and in the report (-m 800 after an out-of-memory failure).
#   asm_meta.json records mode, attempt, k list (from params.txt), -m, threads, peak memory (spades.log), version.
#   On success: K*/, tmp/, misc/ deleted (scaffolds, contigs, graphs, logs kept); prepared reads deleted.
# contigs:<s>
rolib.py asm-stats asm/final/scaffolds.fasta 1000 asm/contigs.fa asm/asm_stats.json   # length with N's, as KBase
# idx:<s> / map:<s>:<lib>   (bowtie2 2.5.5, samtools 1.25)
bowtie2-build --threads 16 --seed 1 asm/contigs.fa bt2/contigs
bowtie2 -x bt2/contigs -1 map_1 -2 map_2 --threads 14 | samtools view -F 0x04 -u - | samtools sort -@ 2 -m 2G -o bam/<lib>.bam -
samtools index bam/<lib>.bam
# depth:<s>   (jgi_summarize_bam_contig_depths from MetaBAT2 2.12.1; CONCOCT 1.1.0 helpers)
jgi_summarize_bam_contig_depths --outputDepth depth/<lib>.d2500.txt --minContigLength 2500 --minContigDepth 1 bam/<lib>.bam   # x n libs
jgi_summarize_bam_contig_depths --outputDepth depth/<lib>.d1000.txt --minContigLength 1000 --minContigDepth 1 bam/<lib>.bam   # x n libs
btlib.py merge-depth asm/contigs.fa depth/metabat_depth.txt depth/*.d2500.txt      # KBase combineDepthTablesUtil logic
btlib.py abund depth/<lib>.d1000.txt depth/abund_<lib>.txt                          # MaxBin2 abundance
cut_up_fasta.py asm/contigs.fa -c 10000 -o 0 --merge_last -b temp.bed > split_contigs.fa
concoct_coverage_table.py temp.bed bam/<libs>.bam > coverage_table.tsv
#   then the BAMs and the bowtie2 index are deleted (bowtie2 logs kept)
# metabat2:<s> / maxbin2:<s> / concoct:<s> / dastool:<s>   (binning-test commands and parameters, unchanged)
metabat2 --inFile asm/contigs.fa --outFile metabat2/bins/bin --abdFile depth/metabat_depth.txt --numThreads 16 --minContig 2500 --maxEdges 500 --unbinned
run_MaxBin.pl -contig asm/contigs.fa -out maxbin2/out/mb -abund_list depth/abund_list.txt -thread 16 -prob_threshold 0.8 -markerset 40 -min_contig_length 1000
concoct --composition_file split_contigs.fa -l 2500 -b out/ --coverage_file coverage_table.tsv -t 16 -k 4 -c 400 -i 500 --total_percentage_pca 90
merge_cutup_clustering.py out/clustering_gt2500.csv > out/clustering_merged.csv; extract_fasta_bins.py asm/contigs.fa out/clustering_merged.csv --output_path bins
DAS_Tool -i metabat2/contigs2bin.tsv,maxbin2/contigs2bin.tsv,concoct/contigs2bin.tsv -l metabat2,maxbin2,concoct -c asm/contigs.fa -o dastool/out/<s> \
  --search_engine diamond --score_threshold 0.5 --duplicate_penalty 0.6 --megabin_penalty 0.5 --write_bin_evals --write_bins -t 16
#   empty-binner rules as the binning test / pipeline.py (MaxBin2 "cannot be binned", DAS Tool "no bins with bin-score")
# checkm2|gunc|checkm1:new:<s>   and   checkm2|gunc:old:<s>  (old = the 2026 DAS Tool bins of the same unit)
checkm2 predict --input <dir> -x fa --output-directory .../checkm2/out --threads 16 --database_path CheckM2_database/uniref100.KO.1.dmnd --tmpdir $TMPDIR
gunc run --input_dir <dir> --file_suffix .fa --db_file gunc_db_progenomes2.1.dmnd --threads 16 --out_dir .../gunc --temp_dir $SCRTMP
CHECKM_DATA_PATH=binning_test/refdata/checkm1 checkm lineage_wf --reduced_tree -x fa -t 16 --pplacer_threads 4 --tab_table -f .../checkm1/qa.tsv --tmpdir $TMPDIR <dir> .../checkm1/out
# keep:<s>   (CheckM2 must have one row per bin, else the task fails and retries)
rolib.py keep --sample <s> --bins dastool/out/<s>_DASTool_bins --checkm2 quality_report.tsv --gunc GUNC.progenomes_2.1.maxCSS_level.tsv \
  --checkm1 qa.tsv --outdir <s>/keep --kept-dir mags_kept
# gtdbtk   (after every keep:<s> has finished or failed; GTDB-Tk 2.8.0, GTDB R232; flags as ~/gtdbtk_kbase/run_gtdbtk_run3.sh)
gtdbtk classify_wf --batchfile gtdbtk/batchfile.tsv --out_dir gtdbtk/out --cpus 48 --pplacer_cpus 8 --min_perc_aa 10 --tmpdir $TMPDIR
# report
rolib.py report --config scripts/config.json          # results/*.tsv + results/REPORT.md; then the driver writes DONE or FAILED
```

**KEEP rule details** (`rolib.py`, binning-test conventions).

- Comparisons are exact decimals on the printed CheckM2 values.
- GUNC pass means `pass.GUNC == True` in `GUNC.progenomes_2.1.maxCSS_level.tsv`, joined on GUNC's name (file name cut at the first `.fa`). A bin with no GUNC row, or a non-True value, does not pass.
- Drop reasons are listed per bin: `completeness<50`, `contamination>=10`, `gunc_fail`, `gunc_na`, `gunc_missing`, `checkm2_missing`.
- HQ = kept with CheckM2 ≥ 90 / < 5.

## Tools and versions

| Tool | Version | Where (poplar) |
|---|---|---|
| sra-tools | 3.4.1 | `/scratch1/afreiburger/metag-hybrid/env/miniforge/envs/sra` (read only) |
| BBTools | 39.81 (openjdk 25) | `/scratch1/afreiburger/emilykin/processed/reanalysis_norm/env_bbtools` (read only) |
| SPAdes | 3.15.3 (python 3.9.19) | `envs/spades3153` |
| bowtie2 / samtools / seqkit / driver python | 2.5.5 / 1.25 / 2.14.0 / 3.12 | `envs/bt_map` |
| MetaBAT2 | 2.12.1 | `envs/bt_metabat2` |
| MaxBin2 | 2.2.7 | `envs/bt_maxbin2` |
| CONCOCT | 1.1.0 | `envs/bt_concoct` |
| DAS Tool / diamond / prodigal | 1.1.7 / 2.2.8 / 2.6.3 | `envs/bt_dastool` |
| GUNC | 1.1.1, progenomes 2.1 db | `envs/bt_gunc`; db `binning_test/refdata/gunc` |
| CheckM1 | 1.1.3, data 2015_01_16 | `envs/bt_checkm1`; data `binning_test/refdata/checkm1` |
| CheckM2 | 1.1.0 | EmilyKin `bin_qc` env `/scratch1/afreiburger/emilykin/processed/.snakemake_envs/8c703f8abdc13cdc1ef374a2bcf37f14_` (read only) |
| GTDB-Tk | **2.8.0 (new env)** | `envs/gtdbtk280` (`gtdbtk280.create.log`, `gtdbtk280.list.txt`) |
| GTDB R232 package | md5 `25a59e0352b1fd150c589f56559767d4` | `refdata/gtdbtk_r232/extract.tmp/release232` (path in `refdata/gtdbtk_r232/.ok`) |

`envs/` and `refdata/` are under `/scratch1/afreiburger/salterns_2026/`.

**GTDB R232** (`scripts/setup_gtdbtk.sh`, log `logs/setup_gtdbtk.log`). It finished 2026-11-13 09:19 CST, before any assembly.

- The 12 split parts of the R232 package were downloaded in parallel from data.gtdb.ecogenomic.org. Each was checked against its md5 in GTDB's `MD5SUM.txt` (kept as `refdata/gtdbtk_r232/MD5SUM.txt`).
- The concatenation was checked against the full-package md5 `25a59e03…` (`full_md5.verified`).
- The package was streamed into `tar` and checked with `gtdbtk check_install` (all OK, `check_install.log`). The parts were then deleted.
- The data directory was left where `check_install` verified it (`extract.tmp/release232`) rather than moved.
- The 26 MB R232 mockup DB (`refdata/gtdbtk_r232_mockup`, md5-checked) is used only by the smoke test.

## Resource policy and scheduler

poplar has 336 cores and 3 TB RAM and is shared: the 1-min load is 300–600, and the CPU was 334/336 busy when measured. The scheduler recomputes everything every 2 min, and right after any task ends.

- **Threads.**
  - `others = load1 − (our threads in tasks running ≥ 2 min)`. Recent launches are not yet in the load average.
  - `free_cores = 336 − others`.
  - New tasks start only while `ours + task ≤ cap`.
    - The task's rule, literally: `cap_literal = min(256, max(64, free_cores + ours − 16))`.
    - That formula counts our threads twice, since `free_cores` already includes them. Below a load of about 350 it would feed back on itself up to 256.
    - So the driver also applies `cap_guard = max(64, free_cores − 16)`. Our total never exceeds what the other users leave free minus 16, with a floor of 64.
    - `cap = min(cap_literal, cap_guard)`, which is at most 256. Both values are logged in `STATE.json` and `ledger.jsonl`.
  - **At today's load (about 360–400 from other users) the cap sits at the 64-thread floor.** That is about 2 SPAdes lanes plus one 16-thread task. It rises only when other users' load drops below about 256.
- **Per-task threads (cost).** SPAdes 24; BBTools 16; bowtie2 map 17 (bowtie2 14 + samtools view 1 + sort 2); MetaBAT2/MaxBin2/CONCOCT/DAS Tool 16; CheckM2/GUNC/CheckM1 16; depth = number of libraries (9–12); download 6; GTDB-Tk 48 (`--pplacer_cpus 8`).
- **Memory.**
  - A task starts only if `MemAvailable − Σ(reserve − current RSS of our running tasks) − its reserve ≥ 400 GB`. Current RSS is summed over each task's session from `/proc`.
  - Reserves: merge/Tadpole = their `-Xmx` (150 or 200 GB; 250 after a Java OOM); SPAdes 80 GB (400 after an OOM retry; the pilot's peak was 37–55 GB); CheckM2/CheckM1 40; GUNC 30; CONCOCT 30; GTDB-Tk 170; others 2–20.
- **Order.** Downstream first: keep, QC, DAS Tool, binners, depth, maps, index, then reads/downloads, assembly and preparation. Within a kind, samples follow `config.json`: R2_restored_D_shore and SF2_A_H2O first (the discriminating batch from `../tadpole/RESULTS.md`), then the R1/SF2 samples whose reads are ready, then R2 and R2A.
  - **Exception: protected assembly lanes.** Assemblies are the long pole, so `min(cap // 64, lanes)` of them go ahead of downstream work, together with the read preparation that feeds them. Room is held for them when they don't fit yet. `lanes = min(6, cap // 40)`.
- **Limits.**
  - At most 6 assemblies, 3 read preparations (merge/Tadpole/trim), 4 downloads, 4 BBDuk preps and 6 QC tasks at once, with at most 3 each of CheckM2 and CheckM1 (they use `/tmp`).
  - A new sample's prep starts only while fewer than `2 + free lanes` samples are prepared and waiting for SPAdes.
  - The background QC of the 2026 bins holds at most 16 threads.
- **Disk.**
  - New downloads and preparations start only while /scratch1 has ≥ 1 TB free and the rollout tree is < 2.6 TB. A walker measures the tree every 30 min.
  - Deleted as soon as they are no longer needed:
    - raw FASTQ after BBDuk;
    - merged/unmerged reads after Tadpole;
    - `ecc_u_*` after the trim;
    - the SPAdes inputs and SPAdes `K*/`, `tmp/` and `misc/` after a successful assembly, plus a failed default-mode directory before the `--meta` fallback;
    - BAMs and the bowtie2 index after the depth/coverage tables;
    - CheckM1 per-bin and GUNC diamond intermediates.
  - Kept: `clean_*` and `map_*` of each library (about 12 GB per library; the 18 reused ones are symlinks), scaffolds/graphs (about 2–3 GB per sample), tables and bins. The expected peak is well under 1 TB.
- **/tmp** (26 GB, shared). Only small files, Python multiprocessing sockets and the CheckM2/CheckM1/GTDB-Tk temp directories go there. CheckM2 and GTDB-Tk set `tempfile.tempdir` from `--tmpdir`, and an NFS socket path fails on cleanup.
  - CheckM2, CheckM1 and GTDB-Tk start only while `/tmp` has ≥ 10 GB free.
  - GUNC, samtools sort, fasterq-dump and SPAdes keep their temporary files on /scratch1.
- Every process runs at nice 10. The driver only reads `/proc`, never signals a process it did not start, and never touches other users' jobs.

## Retries and failures

| Kind | Attempts | On failure |
|---|---|---|
| download | 3 | retry after ≥ 5 min |
| SPAdes | 4 | 1–2 default mode (2nd with `-m 800` if the 1st ran out of memory); 3–4 `--meta` (Arm A recipe), marked `meta_fallback` |
| merge / Tadpole | 2 | the retry gets `-Xmx250g` after a Java OutOfMemoryError |
| everything else | 2 | retry after ≥ 5 min |

- Out of memory means any of these in the previous attempt's logs: `Cannot allocate memory`, `unable to allocate`, `std::bad_alloc`, `OutOfMemoryError`, `OS return value: -9`, `Killed`, `Out of memory`. SPAdes 3.15.3 with too small an `-m` prints `mimalloc: error: unable to allocate OS memory` (tested).
- A task that dies without an exit code counts as a failed attempt.
- Optional tasks (CheckM1, and the 2026-side CheckM2/GUNC) never fail a sample.
- A gate that can never open is turned into a failure, so the run still ends: a GTDB-Tk setup `.FAIL` marker, or 2026 inputs never uploaded once all samples are finished.

## Decisions beyond the task text

1. **The thread cap.** The task's formula counts our threads twice. The driver applies it literally, plus a guard with the evident intent (our total ≤ what other users leave free − 16, floor 64). See "Resource policy".
2. **Protected assembly lanes and their read prep** are an exception to strict downstream-first. Without them, at the 64-thread floor, maps and binners would keep the assembly lane idle. The assemblies set the run's length.
3. **CheckM2 + GUNC for the 2026 side on all 39 units.** CheckM2 had existed only for the two binning-test units. The 2026 DAS Tool bins of the other 37 units were uploaded from `../../work/dastool/` (md5-checked) and are scored at background priority, so the comparison uses CheckM2 on both sides everywhere.
4. **Smoke GTDB-Tk on the real R232 instead of only a mock.** The mockup cannot place archaea.
5. **Timing of the real run.**
   - The download-only pre-run started (15:47) before the smoke had finished. The download step itself had already passed in the smoke.
   - The real driver launched (16:16) while the smoke's real-R232 GTDB-Tk step was still running. GTDB-Tk is the run's last step, days away, and that check then passed (16:23).
   - The smoke and real drivers overlapped for about an hour. That was at most about 30 smoke threads on top of the real run's cap.
6. **R232 location.** The data stays at `refdata/gtdbtk_r232/extract.tmp/release232`, where `check_install` verified it; it was not moved.
7. **Java heaps use `-Xms8g` and `eoom`.** Heaps are no longer reserved in full up front, and a Java OutOfMemoryError exits instead of hanging. Both flags were tested on BBTools 39.81's bbduk, bbmerge-auto and tadpole wrappers.

## Where everything is on poplar

Root: `/scratch1/afreiburger/salterns_2026/rollout/`

| Path | What |
|---|---|
| `scripts/` | `driver.py`, `rolib.py`, `btlib.py` (unchanged from the binning test), `config.json`, `smoke_config.json`, `make_config.py`, `status.sh`, `setup_gtdbtk.sh`, `smoke_prep.sh`, `upload_inputs.sh`; `scripts/tasks/` holds every generated task script |
| `reads/<lib>/` | `clean_{1,2}.fq.gz` (assembly input), `map_{1,2}.fq.gz` (mapping input), `read_counts.json`, `download.json`, `vdb-validate.log`, BBDuk logs; reused libraries hold symlinks and `reused.json` |
| `<s>/prep/` | BBMerge/Tadpole/trim logs and `ihist.txt` (the reads themselves are deleted after use) |
| `<s>/asm/` | `final` → `spades` (or `spades_meta`), `contigs.fa` (≥ 1 kb), `asm_stats.json`, `asm_meta.json` (mode, attempt, k list, `-m`, peak memory), `spades.<attempt>.{out,log}` |
| `<s>/bam/`, `<s>/depth/`, `<s>/{metabat2,maxbin2,concoct}/`, `<s>/dastool/out/<s>_DASTool_bins/` | bowtie2 logs, depth tables, binner outputs, NEW bins |
| `qc/{new,old}/<s>/{checkm2,checkm1,gunc}/` | `quality_report.tsv`, `qa.tsv`, `GUNC.progenomes_2.1.maxCSS_level.tsv` |
| `<s>/keep/` | `bins.tsv` (every bin: CheckM2, CheckM1, GUNC, size, contigs, N50, keep, quality, drop reasons), `kept.txt`, `summary.json` |
| `mags_kept/` | **kept MAG FASTAs**, `<sample>__<bin>.fa` |
| `gtdbtk/` | `batchfile.tsv`, `out/`, `gtdbtk.{bac120,ar53}.summary.tsv` |
| `results/` | `samples_summary.tsv`, `bins_all.tsv`, `mags_kept.tsv`, `bins_dropped.tsv`, `old2026_bins_checkm2_gunc.tsv`, `REPORT.md` |
| `inputs/` | `mags_old_vs_new.tsv`; `old/<unit>/` = the 2026 DAS Tool bins (md5-checked, `.md5`, `.source`, `.ok`); the two pilot units are copied from the binning test |
| `STATE.json` | rewritten every 20 s: scheduler snapshot, per-sample state, running tasks (RSS, params), every task's state |
| `ledger.jsonl` | append-only: driver starts, launches (params, reserve, scheduler numbers), ends (exit, elapsed, sampled peak RSS), scheduler snapshots, end result |
| `.done/`, `.failed/`, `.run/`, `.rc/`, `.state/` | markers: finished tasks, permanent failures, running attempts, exit codes, attempt counters |
| `logs/driver.log`, `logs/tasks/<task>.<attempt>.log` | driver log; one log per task attempt |
| **`DONE`** / **`FAILED`** (+ `FAILED.json`) | **end markers** |
| `smoke/` | the smoke test (same layout) |

## How to check status

```
ssh poplar 'bash /scratch1/afreiburger/salterns_2026/rollout/scripts/status.sh'
ssh poplar 'cat /scratch1/afreiburger/salterns_2026/rollout/DONE /scratch1/afreiburger/salterns_2026/rollout/FAILED 2>/dev/null'
ssh poplar 'cd /scratch1/afreiburger/salterns_2026/rollout && /scratch1/afreiburger/salterns_2026/envs/bt_map/bin/python -I scripts/driver.py report --config scripts/config.json'   # interim report
rsync -a -e "ssh -o BatchMode=yes" poplar:/scratch1/afreiburger/salterns_2026/rollout/results/ results/
```

`status.sh` prints:

- load, MemAvailable, the driver PID and whether it is alive, the end marker, STOP, the GTDB-Tk setup state and free disk;
- the scheduler snapshot: load, our threads, cap, memory, rollout usage, and why ready tasks are held back;
- task counts, each sample's state (last finished step, running steps), failed tasks and running tasks (elapsed, RSS, SPAdes mode/`-m`, heap);
- the last driver log lines.

`status.sh /scratch1/afreiburger/salterns_2026/rollout/smoke` shows the smoke test.

**Done / failed.**

- `$ROOT/DONE` = `<utc> samples=39/39 kept_mags=<n> gtdbtk=done optional_failed=<n>`.
- `$ROOT/FAILED` = `<utc> FAILED: incomplete samples ...; failed tasks ...` with details in `FAILED.json`.
- Either one is written only after GTDB-Tk and the report have run on what finished.

**Stop, resume, retry.**

- `touch $ROOT/STOP` makes the driver start nothing new; running tasks finish. `rm $ROOT/STOP` resumes.
- To restart the driver, e.g. after editing `driver.py`: `kill $(cat $ROOT/driver.pid)`. The tasks are detached and keep running. Then start it again with the launch command in "Launch record". It adopts running tasks (`.run/<task>.json` + `/proc`) and skips finished ones (`.done/`).
- New code applies to new launches only, because each attempt's script is written at launch.
- **Never run two drivers on the same root.** Before starting one, this must print 0: `pgrep -f '^/scratch1/afreiburger/salterns_2026/envs/bt_map/bin/python -I scripts/driver.py run --config scripts/config.json' | wc -l`. The anchor `^` keeps the pattern from matching the ssh command line itself.
- Retry a permanently failed task: `rm $ROOT/.failed/<task>` (a fresh round of attempts starts). If the run already ended, also `rm $ROOT/.done/gtdbtk $ROOT/.done/report $ROOT/FAILED`, then restart the driver.

## Smoke test (2026-11-13, passed: `DONE`)

Root `rollout/smoke/`. It ran the same driver code with `scripts/smoke_config.json`: 8-thread tasks, 8 GB heaps, and 20 s scheduling and retry waits. `scripts/smoke_prep.sh` prepared it. It ran from 15:29 to 16:23 UTC and ended with **`DONE` = `samples=3/3 kept_mags=8 gtdbtk=done optional_failed=2`**.

**Inputs.**

- **Library SMK_DL:** R2_B_H2O (SRR5209948), downloaded and verified in full, then cut to its first 300,000 pairs.
- **Library SMK_RU:** the first 300,000 pairs of the binning test's R1_B_D1 reads (the reuse path).
- **Samples:**
  - SMK_A: the normal path, with SPAdes `-m 1` on attempt 1 to force a real out-of-memory failure.
  - SMK_B: default-mode SPAdes forced to fail on attempts 1–2 with an invalid flag; no 2026 inputs uploaded.
  - SMK_I: the import path, on the binning test's SF2_C_D1 outputs.

**What ran and passed.**

| Path | Result |
|---|---|
| download + checks | vdb-validate "consistent"; 20,319,661 pairs and 6,136,537,622 bases = `manifest.tsv`; 0 non-bare separator lines; about 8 min |
| BBDuk prep / reuse links | chained BBDuk counts asserted; subsampled reuse copies |
| Arm C prep | BBMerge 46.2 % joined; Tadpole two passes; quality trim; merged/unmerged/`ecc_u_*` deleted after use |
| **SPAdes out-of-memory retry** | attempt 1 (`-m 1`) died with `mimalloc: error: unable to allocate OS memory … Cannot allocate memory`. The driver detected it and gave attempt 2 `-m 32` (the smoke's bumped value) and the bumped reserve. Default mode then finished: `asm_meta.json` has mode `default`, attempt 2, k `33,55,77,99,127`, peak 4.6 GB |
| **`--meta` fallback** | SMK_B attempts 1–2 failed. Attempt 3 ran `--meta` into `spades_meta/` and finished. `mode = meta_fallback` is in `asm_meta.json` and the report |
| binning chain | contigs ≥ 1 kb, bowtie2 index, 2 maps per sample, depth + CONCOCT coverage table, then BAMs and index deleted; MetaBAT2, MaxBin2 and CONCOCT; DAS Tool (1 bin per sample) |
| QC + KEEP | CheckM2, CheckM1, GUNC; KEEP table; copies in `mags_kept/` |
| import (SMK_I) | contigs byte-identical (`cmp`) to the binning test's; 11 bins with 11 CheckM2 / GUNC / CheckM1 rows. Kept 6 of 11 (2 HQ), matching the binning test (7 MQ+, 6/7 GUNC pass, 2 HQ) |
| 2026-side QC | CheckM2 + GUNC on SMK_A's 6 bins. SMK_B's never-uploaded inputs became optional failures once every sample had finished |
| report | `results/*.tsv` and `REPORT.md`, final and interim (`driver.py report`) |
| GTDB-Tk | see below; 8 / 8 kept MAGs classified |
| end markers | the `FAILED` path (round 1: GTDB-Tk on the mockup failed) and the `DONE` path (round 2) |
| resume | three driver restarts each adopted the tasks still running (dl + 2 QC; 4 QC/DAS Tool; later rounds) |
| killed task | `kill -TERM` of `gunc:new:SMK_B`'s process group gave `exit None`, which counted as a failed attempt |
| STOP | `STOP` held the retry; removing it launched attempt 2 |
| manual retry | removing `.failed/gtdbtk` + `.done/report` + `FAILED` and restarting re-ran GTDB-Tk, the report and the end marker |

**GTDB-Tk in the smoke.**

1. **First try: the R232 mockup DB** (26 MB, meant for CI). It ran the ANI screen, marker identification and alignment, then pplacer failed.
   - The mockup has no backbone/class-level split, so it needs `--full_tree`. With it, the 6 bacterial MAGs classified (rc 0).
   - Its archaeal package is empty, so it cannot place archaea. Two of the smoke MAGs are haloarchaea.
   - So the mockup cannot stand in for this data set.
2. **Second try: the production command on the real, verified R232** (8 CPUs, `--pplacer_cpus 2`, `--min_perc_aa 10`, no `--full_tree`).
   - 8 genomes in 10.5 min, with a peak RSS of 97.8 GB, inside the 170 GB reserve.
   - `/tmp/ro_*` held at most 32 MB.
   - Classes: 2 *Halorubrum* (Haloferacaceae), 6 Desulfobacterota/Gammaproteobacteria to species level.

**Problems found and fixed before or at launch.**

1. **`test -s a && test -s b` does not stop a `set -e` script when `a` fails** (only the last member of an `&&` list counts). The same holds for `cmp … && echo`. Found by review before the smoke: a `need` helper and separate lines replace them.
2. **Ledger write.** It collided on the `utc` key (`cycle error`; the launches still happened). Fixed.
3. **Download order.** Downloads started in alphabetical order. Library tasks now take the order of the sample whose reads they are.
4. **Leftover SPAdes directory.** A failed default-mode `spades/` directory would have stayed beside `spades_meta/`. It is now deleted before the fallback.
5. **Temp dirs of killed attempts.** An attempt killed by a signal skips its EXIT trap and left its temp dirs. The driver now removes them.
6. **Log overwrite on manual retry.** The first manual retry (still on the old code) reused attempt number 1 and overwrote `gtdbtk.1.{sh,log}` from round 1; round 1's attempt 2 log, with the same pplacer error, survives. Launch numbers are now absolute (files are never overwritten), and the attempt within the current round sets the SPAdes mode. Checked locally: 4 attempts default/default/meta/meta, then a fresh round from launch 5.
7. **Hold blocking tiny tasks** (found at the real launch). Room held for a protected task also blocked 1-thread `keep` tasks. Tasks of ≤ 2 threads are now exempt.

**Not exercised by the smoke:** a Java OutOfMemoryError retry with `-Xmx250g` (same detection code as SPAdes); a sample with zero DAS Tool bins end to end (each smoke sample had ≥ 1); the disk gates.

## Launch record

All times UTC on 2026-11-13; poplar's clock shows CST (UTC−6).

| When | What | PID |
|---|---|---|
| 14:30 | `setup_gtdbtk.sh`: env built 14:32; R232 parts 14:42; whole md5 14:45; extracted 15:07; `check_install` OK **15:19** | 1837739 (exited) |
| 15:29–16:23 | smoke test, ending `DONE` | 2672809 → 2688930 → 3036640 → 3362957 (exited) |
| 15:47 | **download-only pre-run** of the real driver (`--only-kinds dl`, log `logs/driver_prerun_dl.log`). The download path had already passed in the smoke. Restarted at 15:48 with the library-order fix | 2951006 → 2966497 (stopped 16:16) |
| 16:16:08 | **full driver launched**; it adopted the 4 running downloads | 3416667 (stopped 16:18) |
| **16:18:39** | **restarted with the tiny-task hold exemption**; it adopted 6 running tasks; 0 cycle errors | 3455692 (stopped 2026-11-16 17:08:55) |
| **2026-11-16 17:09:04** | **restarted for the R2A migration** (`scripts/restart_driver.sh`, see "R2A migrated to frankensteinsmonster"); it adopted 3/3 running tasks; 0 cycle errors | **973270** |

```
cd /scratch1/afreiburger/salterns_2026/rollout
nohup setsid nice -n 10 /scratch1/afreiburger/salterns_2026/envs/bt_map/bin/python -I scripts/driver.py run \
    --config scripts/config.json >> logs/driver.log 2>&1 < /dev/null &
```

**State at 16:24 UTC.**

- **R1_A_D1 complete** (imported): 71 NEW bins, **54 kept (26 HQ)**. Dropped: contamination ≥ 10 ×7, GUNC fail ×15 (8 by GUNC alone), completeness < 50 ×2. Kept median N50 16.5 kb, 155 Mbp.
- **SF2_C_D1 complete** (imported): 11 bins, **6 kept (2 HQ)**.
- The 16 other R1/SF2 libraries are linked to the binning test's reads.
- `merge:SF2_A_H2O` (heap 200g, 12.1 Gbp) and `merge:R2_restored_D_shore` (150g, 6.8 Gbp) are running, as the first discriminating batch.
- R2/R2A downloads: 7 of 21 done at 16:25, 4 running, the rest held behind the protected merges. CheckM2 on the 2026 bins of R2_restored_D_shore is done.
- **Scheduler:** load1 about 410–470, others about 395, **cap 64 (the floor)**. MemAvailable about 2.4 TB. /scratch1 9.9 TB free, /tmp 26 GB free. Rollout tree 269 GB, mostly raw R2/R2A FASTQ waiting for BBDuk.

**Expected timeline.** This is an estimate, at nice 10 on a host that is CPU-saturated by other users.

- At today's load the driver stays at the **64-thread floor**: about 2 SPAdes lanes when downstream work is idle, and 1 protected lane plus its read prep otherwise.
- Pilot default-mode SPAdes took 9.2 h (R1-type) and 13.3 h (SF2-type) at 12 threads. Assume about 6–9 h per sample at 24 threads.
- 37 assemblies plus about 1–3 h of binning and QC per sample: **about 7–12 days**, so `DONE` around 2026-11-20 to 11-25.
- If other users' load stays below about 250, the cap rises toward 256 (up to 6 assemblies) and it could finish in 3–5 days.
- The first per-sample results (R2_restored_D_shore, SF2_A_H2O) need all 12 R2 coverage libraries, or the 9 SF2 ones, plus their assemblies, so expect them in about 1–2 days.

**When `DONE` (or `FAILED`) exists:**

```
rsync -a -e "ssh -o BatchMode=yes" poplar:/scratch1/afreiburger/salterns_2026/rollout/results/ results/
rsync -a -e "ssh -o BatchMode=yes" poplar:/scratch1/afreiburger/salterns_2026/rollout/{DONE,FAILED,FAILED.json,STATE.json} results/ 2>/dev/null
```

Then:

1. Review `results/REPORT.md` and `samples_summary.tsv`: per sample, NEW vs 2026 on both CheckM2 and KBase CheckM, the drop reasons and the `meta_fallback` samples.
2. Bring Andrew the replacement decision. The KEEP rule already drops GUNC failures and CheckM2 contamination ≥ 10.

The kept MAG FASTAs are in `mags_kept/` on poplar (about 1–2 GB expected).

## Independent review (2026-11-13, 16:25–16:41 UTC)

A second agent checked the running rollout on poplar. **Nothing needed fixing**, so no job was killed or relaunched and no script was changed. PID 3455692 is still the driver.

What was checked:

- **Scripts.** The md5 of every file in `scripts/` is the same locally and on poplar. `btlib.py` is byte-identical (`cmp`) to `binning_test/scripts/btlib.py`, locally and on poplar.
- **Read prep and SPAdes vs the pilot (item 1).**
  - The `java`/`Executing` lines of the pilot's `R1_A_D1/logs/A1_bbmerge.log`, `A2a`/`A2b` Tadpole and `A3_qtrim` were compared with the running `SF2_A_H2O/prep/bbmerge.log` and the generated task scripts. They differ only in the intended ways:
    - `t=16` instead of 24;
    - `-Xms8g` plus `eoom` (`-XX:+ExitOnOutOfMemoryError` is on the running JVM's command line);
    - the heap by sample size (SF2_A_H2O, 12.1 Gbp: `-Xmx200g`; R2_restored_D_shore, 6.8 Gbp: `-Xmx150g`).
  - BBMerge (`rem k=62 extend2=50 ecct prefilter=2 ordered`), both Tadpole passes (`mode=correct k=62 prefilter=2`, `extra=`), the trim (`qtrim=r trimq=10 minlen=45`) and the BBDuk passes are unchanged.
  - The SPAdes command equals the pilot's Arm C `Command line` from `params.txt` (`--only-assembler -k 33,55,77,99,127 -m 400 --phred-offset 33`, `-1/-2` + `--merged`), except `-t 24`.
  - The smoke shows both retry paths in the files: `SMK_A` has `spades.1/2.{out,log}`, and the ledger records `oom_previous: mimalloc … Cannot allocate memory`. `SMK_B`'s `asm_meta.json` is `meta_fallback`, attempt 3, with `--meta` in its command line.
- **Rollout fixes.**
  - The separator fix runs `awk` under `set -e` and checks both mates: pairs, bases, line count and non-bare separator lines.
  - Heaps are sized as specified.
  - Intermediates are deleted: `reads/R2_restored_D_shore` has no `raw_*` or `trim_*` left.
  - Task logs are per launch. Tool logs are not (see the limitation below).
- **Reads (item 4).**
  - The pair and base counts of all 19 libraries prepared so far equal `manifest.tsv`, and their BBDuk counts chain.
  - The 18 reused libraries resolve to `binning_test/reads/<lib>/`, or to `tadpole_pilot/<s>/qc/` for the two pilot libraries. They were made with the same BBDuk parameters.
- **Binning and QC vs the binning test (item 2).** The commands match `binning_test/scripts/driver.py` line by line: bowtie2 `--threads 14` with `samtools view -F 0x04` and `sort -@ 2 -m 2G`; jgi 2500/1000; MetaBAT2, MaxBin2, CONCOCT and DAS Tool with the same flags; CheckM2, CheckM1 and GUNC. The only differences are the documented ones:
  - CONCOCT's cut-up and coverage table moved into the depth task;
  - GUNC's temp directory is on /scratch1;
  - map tasks are costed at 17 threads.
- **KEEP rule (item 3).** `rolib.py classify`: completeness ≥ 50, contamination < 10 (Decimal), `pass.GUNC == "True"`. A missing or `nan` GUNC row drops the bin.
- **Coverage groups and reuse (item 4).**
  - `make_config.py` groups by the manifest's `coverage_group`, the same column `pipeline.py` uses.
  - R2 holds 12 libraries, including R2_restored_C_black, R2_restored_D_shore and R2_restored_H2O. R1, R2A and SF2 hold 9 each.
  - The binning test's R1_A_D1 and SF2_C_D1 lists equal the rollout's R1 and SF2 groups.
  - The two imports are verified:
    - `cmp` of `contigs.fa` passed;
    - 71/71/71/71 rows for R1_A_D1 and 11/11/11/11 for SF2_C_D1;
    - 54 + 6 = 60 MAGs in `mags_kept/`.
- **Scheduler (item 5).**
  - All rollout processes run at nice 10.
  - A 20 s `/proc` tick measurement gave **20.3 cores actually used** against 56 booked threads (2 merges, 4 downloads) and a cap of 64. Load1 was about 445–490.
  - MemAvailable was 2.47 TB.
  - The 6 `/tmp/ro_*` directories (one per running task) are empty. /tmp has 26 GB free.
  - The driver sends no signal other than `kill -0`.
  - Two jobs of the same account (`biot` `cc_full.py` and `pc_full.py`, nice 0) are not part of this rollout and were not touched.
  - Booked threads exceed real use, so `others` is underestimated by that excess. In terms of real CPU use, the cap still keeps the rollout at or below (336 − other users − 16).
- **GTDB-Tk (item 6).**
  - All 12 parts are `ok` against `MD5SUM.txt`, and the whole package is `25a59e0352b1fd150c589f56559767d4`.
  - `check_install` reports OK for all 9 reference directories. Its "NOT FOUND" third-party lines appear only because it ran without the env on PATH. The task puts `envs/gtdbtk280/bin` on PATH, and pplacer, guppy, FastTree(MP), hmmer, prodigal and skani are all there.
  - The smoke's production command classified 8/8 genomes.
- **Resume (item 7).** `ledger.jsonl` and `driver.log` show three restarts adopting the running tasks (4, 4 and 6). No cycle errors.

**Known limitation (not fixed, documented).** Tool logs have fixed names and are not kept per launch:

- `prep/bbmerge.log`, `tadpole_{pairs,merged}.log` and `qtrim.log`;
- `reads/<lib>/bbduk_*.log`;
- the binner and QC logs inside directories that the next attempt removes with `rm -rf`, such as `metabat2/`, `qc/new/<s>/checkm2/` and `dastool/`.

A retry therefore overwrites the previous attempt's error text, and the per-launch task log keeps only `tail`s from the success paths.

- The out-of-memory detection is not affected: `plan_params` reads the previous attempt's logs when the retry is launched, before the retry runs.
- After a permanent failure, the last attempt's tool log survives.
- What can be lost is why attempt 1 failed when attempt 2 then succeeded.
- This is forensic debt, not an error, so the driver was not restarted for it. If a retry happens, read the ledger's `oom_previous`.

**To watch.**

1. **GTDB-Tk at about 1,000 MAGs.** `--tmpdir` is on /tmp (26 GB, shared), and the 10 GB gate is checked only at start. The smoke used 32 MB for 8 genomes.
2. **CheckM2 `--tmpdir` on /tmp,** with up to 3 runs at once. The binning test ran a 71-bin set this way without trouble.
3. **The `--meta` fallback follows any two failed default-mode attempts,** including a SPAdes run killed by a reboot or an administrator. Every such sample is marked `meta_fallback` in `asm_meta.json` and in `samples_summary.tsv`. Check why its default-mode attempts failed and, if it was not the recipe, rerun it in default mode (`rm .failed/…`, see "Stop, resume, retry").
4. **The launch summary's "CheckM2 + GUNC on both sides for all 39 units" is the plan, not yet the state.** At review, the 2026 side was scored for R1_A_D1 and SF2_C_D1 (imported) and for R2_restored_D_shore. The other 36 units are queued at background priority (16 threads).
5. **Clock skew.** /scratch1 file times run about 3.5 min behind poplar's clock: `.done/dl:R2_restored_D_shore` shows 10:11, while its log ends at 10:14:45 CST. The driver never uses file times.

## R2A migrated to frankensteinsmonster

**Approval.** On 2026-11-16 Andrew approved moving the R2A samples from poplar to frankensteinsmonster, to finish sooner. That machine is a Threadripper 3970X (64 threads), 251 GB RAM, shared with other sessions.

**What runs where.**

| Samples | Where |
|---|---|
| R2A_A_D1, R2A_A_D2, R2A_B_D1, R2A_B_D2, R2A_C_D1, R2A_C_D2, R2A_B_H2O (6.65–8.80 Gbp) | **frankensteinsmonster**, root `~/Documents/salterns_rollout_local/rollout/` (its README: `~/Documents/salterns_rollout_local/README.md`) |
| **R2A_A_H2O (9.65 Gbp), R2A_C_H2O (10.20 Gbp)** | **still poplar**, as ordinary samples |
| the other 30 samples; GTDB-Tk on all kept MAGs; the report; `DONE` | poplar, unchanged |

**Why two R2A samples stay on poplar.**

- The recipe gives samples above 9 Gbp a 200g Java heap.
- The migration rule was: start a merge or Tadpole only when MemAvailable ≥ heap + 40 GB. For a 200g heap that is 240 GB.
- frankensteinsmonster has 251 GB in total, and its driver scope is capped at 220G, so that condition can never hold. A smaller heap would change a recipe parameter, so it was not used.
- R2A_A_H2O and R2A_C_H2O need no reads copied back: all 9 R2A libraries stayed on poplar.

**Status when this section was written (17:10 UTC).** No R2A task had started on poplar beyond the downloads, BBDuk and the 2026-side QC of six R2A units. All 9 `dl:`/`reads:` markers were done. No `merge:R2A_*` had run, because `prep_buffer` held them and the samples were near the end of the order.

### Exactly what changed on poplar

All under `/scratch1/afreiburger/salterns_2026/rollout/`.

1. **Backups**, as `scripts/<file>.pre_r2a_20261116`:

   | File | md5 before |
   |---|---|
   | `driver.py` | `f99883b9…` |
   | `rolib.py` | `db74f74c…` |
   | `config.json` | `a6e6e228…` |
   | `make_config.py` | `c68b13b9…` |
   | `status.sh` | `b0cc2a74…` |

2. **`scripts/driver.py`** (now `507e1670…`; the same file runs on frankensteinsmonster). Everything new is switched on by config keys only.
   - **Staged import**, a second `import_from` variant: `{"staged": "<root>/imports/<s>", "marker": "READY", "source": "..."}`.
     - `gate()` holds `import:<s>` as **`import_wait`** until `imports/<s>/READY` exists.
     - The import body runs `md5sum -c MANIFEST.md5` and checks that `asm/contigs.fa`, `asm_meta.json`, `asm_stats.json`, `final/scaffolds.fasta` and the DAS Tool bins directory exist.
     - It copies `imports/<s>/<s>/` to `<s>/` and `imports/<s>/qc/new/<s>/` to `qc/new/<s>/`.
     - `rolib.py verify-import` checks one CheckM2/GUNC row per bin (and CheckM1 rows when present).
     - It then re-applies the KEEP rule in its scratch directory. The result must equal the remote `keep/kept.txt` (`cmp`), and the kept FASTAs must be byte-identical to the remote `mags_kept/` (md5).
     - It writes `<s>/IMPORTED.json`.
     - After that, `keep:<s>` runs here unchanged and fills `mags_kept/`.
     - Unlike the binning-test imports, staged samples **keep their 2026-side `checkm2:old`/`gunc:old` tasks** (the R2A_B_H2O pair is still pending). Libraries used only by staged samples are not "used" here. That changes nothing in this run, because R2A_A_H2O and R2A_C_H2O still use all 9 R2A libraries.
   - **`fail_stuck`:** if the other host writes `imports/<s>/FAILED`, `import:<s>` is marked failed and the run ends `FAILED` for that sample instead of waiting forever.
   - **Used only by frankensteinsmonster's config:**
     - a `push:<s>` task kind (below);
     - `asm_after_all_prep`, `protect_lanes`, `scope_mem_gb` (the sum of max(reserve, RSS) of running tasks must stay below the scope cap);
     - `skip_gtdbtk` and `skip_old_qc`. `write_end` and `fail_stuck` treat GTDB-Tk as required only where it is a task, so poplar's behaviour is unchanged.
   - **Normal samples are unchanged.** Every generated task script of R2A_B_D1 under the old config (28 scripts) is identical between the old and the new driver.
3. **`scripts/rolib.py`** (now `89425e47…`). Report only:
   - a new `processed_on` column in `samples_summary.tsv`;
   - `imported_from_binning_test` and the 2026-side source say "binning test" only for the two binning-test imports;
   - `REPORT.md` lists the samples imported from frankensteinsmonster.
4. **`scripts/make_config.py`** (now `7e9a2575…`) gains `MIGRATED` (the 7 samples, each given the staged `import_from`) and `host: "poplar"`.
   - The regenerated **`scripts/config.json`** (now `e924414a…`) differs from the old one only by `"host": "poplar"` and the 7 `import_from` blocks.
   - The R2A samples keep their `lib`/`group`, so the report still counts 9 coverage libraries.
   - `smoke_config.json` regenerated byte-identical.
5. **`scripts/status.sh`** prints a "staged imports" block (READY / imported per sample).
6. **New files:**
   - `scripts/restart_driver.sh`: STOP, wait for the acknowledgement in `STATE.json`, SIGTERM only while the main thread is in `hrtimer_nanosleep`, check that no driver runs, relaunch with the README command, count adoptions, remove STOP;
   - `scripts/make_config_local.py` and `scripts/local_config.json` (copies; they run on frankensteinsmonster);
   - `imports/README.txt`;
   - `README_rollout.md`, a copy of this README in the rollout root, so the record is also on poplar.
7. **Driver restart:** `bash scripts/restart_driver.sh $ROOT scripts/config.json --rm-stop`.
   - 17:08:43 STOP was created and acknowledged at 17:08:55.
   - PID 3455692 was stopped while sleeping.
   - **New PID 973270** started at 17:09:04 and adopted 3/3 running tasks: `asm:R1_B_H2O`, `tadpole:SF2_B_H2O`, `map:R1_A_H2O:R1_A_D2`.
   - STOP was removed at 17:09:09. 0 cycle errors.
   - There are now 910 tasks, down from 1064: the 7 samples' 161 chain tasks are gone and 7 imports were added.
   - The scheduler shows `import_wait: 7`, and the other samples continue.
8. **Checked before the restart.**
   - A dry run of the old and new driver on a copy of poplar's markers (twice, the last time at 17:08): the only task-set changes are those above.
   - No `.done` marker was orphaned, and no sample state changed.
   - The ready set differs only by the 7 R2A merges (now gone) and the 7 waiting imports.

**`DONE` waits for R2A.** The chain is the existing one:

- `keep:<s>` depends on `import:<s>`;
- `gtdbtk` waits (soft) for every `keep`;
- `report` waits (soft) for every task.

Nothing turns a waiting import into a failure, except `imports/<s>/FAILED`. So poplar writes `DONE`/`FAILED` only after all 7 imports, their keeps and GTDB-Tk on all kept MAGs, R2A included.

### Import-back (automatic)

1. **On frankensteinsmonster**, `push:<s>` runs after `keep:<s>`.
   - It writes `MANIFEST.md5` (md5 of every file of `<s>/`, `qc/new/<s>/` and `mags_kept/<s>__*.fa`, with root-relative paths).
   - It runs `ssh poplar rm -f imports/<s>/READY` and creates `imports/<s>.part/`.
   - It rsyncs the files there with their relative paths, plus `MANIFEST.md5` and `PUSHED.json`.
   - It then runs `md5sum -c` **on poplar** and only then renames `.part` to `imports/<s>/` and touches `READY`.
   - ssh uses `BatchMode=yes` through the existing ControlMaster. There are 4 attempts, at least 5 min apart.
2. **On poplar**, the waiting `import:<s>` starts within about 2 min of `READY` and runs as described above. Then `keep:<s>` runs, and GTDB-Tk / report / `DONE` wait for all of it.

**Manual fallback.** Use this if a push fails permanently, for example because the ssh master is gone. On frankensteinsmonster:

```
cd ~/Documents/salterns_rollout_local/rollout
ssh -o BatchMode=yes poplar true && rm .failed/push:<s>   # the driver retries the push (fresh round)
# or run the same script by hand:  bash scripts/tasks/push__<s>.<n>.sh   (the last generated push script)
```

**If a sample cannot be finished on frankensteinsmonster** (for example a Java OutOfMemoryError at 150g, whose retry there also gets 150g):

- Either write one line to `imports/<s>/FAILED` on poplar; the run then ends `FAILED` naming it.
- Or run it on poplar instead:
  1. Remove it from `MIGRATED` in `make_config.py` and regenerate `config.json`.
  2. Restart with `restart_driver.sh`. It becomes an ordinary sample again: merge → … → keep, and its reads are already here.
  3. Stop it on frankensteinsmonster (`touch STOP` there) if it is still running.

### Test (2026-11-16, passed)

**Local smoke** (`~/Documents/salterns_rollout_local/smoke/`). It used the same driver in its own scope, `salterns-r2a-smoke-*`, with MemoryMax 64G and the memory guard.

- **Round 1** used 300k pairs of R2A_B_D1 sediment reads. It assembled to a single contig ≥ 1 kb (0 ≥ 2.5 kb), so MetaBAT2 and CONCOCT had nothing to bin. Kept as `smoke_round1_300k_R2A_B_D1/`.
- **Round 2** used 1M pairs of R2A_B_H2O plus R2A_C_H2O (water) and went through the whole chain:
  - reuse/subsample, merge, Tadpole, trim, SPAdes default mode (0.36 h), contigs, index, 2 maps, depth;
  - MetaBAT2, MaxBin2, CONCOCT, DAS Tool (2 bins);
  - CheckM2, GUNC and CheckM1, then keep (1 kept, HQ), push, report;
  - `DONE samples=1/1 kept_mags=1 optional_failed=0`.
- The first round-2 CheckM1 failed. The local env's `checkm/DATA_CONFIG` pointed at `~/.checkm`. poplar's env points at its own data copy (`envs/bt_checkm1/checkm_data`), md5-identical to `binning_test/refdata/checkm1` (32/32 files). The local `DATA_CONFIG` now points at the verified local copy. The rerun passed, and the sample was re-pushed.

**Poplar side** (`rollout/smoke_import/`). It used the same sample as a staged import, with no GTDB-Tk and the 2026-side QC on a copy of the six SF2_C_D1 bins.

- The driver waited in `import_wait` while the 2026-side QC ran.
- **Bug found and fixed:** the import checked the empty `READY` with `need` (non-empty test). It now uses `test -e`. The import was then retried with the fixed driver:
  - md5 OK for 186 files;
  - verify-import gave 2/2/2/2 rows;
  - the KEEP rule recomputed on poplar equalled the remote one;
  - the kept FASTA md5 `e33d1403…` is identical on both hosts;
  - keep → `mags_kept/`, report (`processed_on = frankensteinsmonster:…`), `DONE samples=1/1`.
- The same root rehearsed `restart_driver.sh` (STOP, kill while sleeping, relaunch, remove STOP).
- Two more import checks, run locally with `driver.py run-one` on a copy of the staged test data:
  - a sample pushed **without CheckM1 rows** (CheckM1 is optional) imports and keeps (2/2 CheckM2/GUNC rows, same kept set);
  - a staged file changed by one byte makes the import fail at `md5sum -c`.

**Task commands, poplar vs frankensteinsmonster** (`~/Documents/salterns_rollout_local/task_script_diff_R2A_B_D1.txt`, every task of R2A_B_D1).

- After normalising paths, the only differences are:
  - SPAdes `-t 24` → `-t 16` and `-m 400` → `-m 120` (the `--meta` fallback likewise);
  - bowtie2 `--threads 14` → `13` (map task 17 → 16 threads);
  - `OMP_NUM_THREADS` and the reserve shown in each task's start line.
- Java heaps (`-Xms8g -Xmx150g`), every BBTools, SPAdes, binner, DAS Tool, CheckM2/CheckM1/GUNC and KEEP parameter, and the k list are identical.
- The tool envs there were built from poplar's exact package URLs and md5s, except DAS Tool, which reuses `~/Documents/dastool_env`: the same das_tool, diamond, prodigal and blast builds, but a different r-base build.

### Status, stop, resume

```
ssh poplar 'bash /scratch1/afreiburger/salterns_2026/rollout/scripts/status.sh'      # poplar (shows import_wait + staged imports)
bash ~/Documents/salterns_rollout_local/status_local.sh rollout                       # frankensteinsmonster
```

- **poplar:** as before (`touch STOP` / `rm STOP`). Restart with `bash scripts/restart_driver.sh /scratch1/afreiburger/salterns_2026/rollout scripts/config.json --rm-stop`.
- **frankensteinsmonster:**
  - `touch ~/Documents/salterns_rollout_local/rollout/STOP` / `rm` it.
  - To restart: `kill $(cat rollout/driver.pid)` while it sleeps, then `bash ~/Documents/salterns_rollout_local/launch_local.sh rollout`. Each start gets a new scope plus a guard; the old scope lives on with its adopted tasks until they end.
  - Never `systemctl --user stop` a scope.

### ETA (estimate, 2026-11-16 17:10 UTC)

**poplar**

- Since the first assembly, 13 assemblies have finished in 71.4 h, about 5.5 h each (about 1.6 lanes at the 64-thread floor).
- 17 assemblies remain: R1_B_H2O (running), R1_C_H2O, SF2_B_H2O, SF2_C_H2O, 11 × R2, R2A_A_H2O and R2A_C_H2O. That is about 90 h, so the last assembly finishes around **2026-11-20 12:00 UTC**.
- Then come binning and QC of the last sample (about 3–4 h) and GTDB-Tk on about 1,000 MAGs. GTDB-Tk takes several hours to a day on this shared host, and its `/tmp` gate must be open. So **`DONE` around 2026-11-21 (to 11-22)**.
- With all 24 assemblies (before the migration) it would have been about 2026-11-22/23. The 7 migrated samples save about 7 × 5.5 h ≈ **1.6 days**.

**frankensteinsmonster** (`~/Documents/salterns_rollout_local/README.md`)

- **Read preparation**, one sample at a time and all seven before any assembly: about 6–8 h, provided MemAvailable reaches 190 GB.
- **Assemblies:** 7 SPAdes runs at 16 threads, two at a time, about 5–8 h each. That is about 20–30 h, with binning and QC overlapping in the third 16-thread lane.
- **Push** of each sample as soon as it is kept.
- **All 7 imported** about 1.5 days after the first merge starts, around **2026-11-18**, well before poplar's own assemblies end.

### Risks

1. **The local memory gate depends on other sessions.**
   - Each merge or Tadpole run needs MemAvailable ≥ 190 GB, counting our own outstanding reserves.
   - Over 6 h on 2026-11-15/16 (idle, before launch) that held 95% of the time (median 209 GB).
   - At launch (17:10 UTC) it was 178 GB because of another session's 31 GB job. The first merge started at 17:28, when MemAvailable reached 209 GB.
   - A long spell below 190 GB delays only the frankensteinsmonster samples, and they have about 2 days of slack against poplar.
2. **A Java OutOfMemoryError at 150g cannot be retried with a bigger heap there.** That sample goes back to poplar (above). poplar's 150g preps peaked at 72–110 GB RSS.
3. **The push needs the ssh ControlMaster to poplar.** If it is gone, the pushes fail after 4 attempts and need the manual fallback. The local results are unaffected.
4. **poplar `/tmp`** was 5.4 GB free on 2026-11-16 04:40 UTC and 22 GB at 17:00, against the 10 GB gate for CheckM2/CheckM1/GTDB-Tk ("To watch" 1 and 2). It is not ours.
5. **The local machine is shared.** If it runs out of memory, the kernel prefers our tasks (oom_score_adj +300), and the guard freezes our scope below 12 GB. A frozen task keeps its RAM, and a killed one is retried (and counted as an OOM).
6. **Memory booking.** The local SPAdes reserve is 60 GB, and the scope-memory check books max(reserve, RSS) at launch.
   - Two SPAdes (60 + 60) plus CheckM1 (80) pass the 210 GB check.
   - If both assemblies then grow towards `-m 120` while CheckM1's pplacer is near 75 GB, the 220G cgroup cap kills the largest process, which is then retried.
   - poplar's peaks for samples of this size were 29–56 GB. If the first local `asm_meta.json` shows a peak above 60 GB, raise `reserve_gb.asm` in `make_config_local.py` and restart the local driver.

## R2 migrated to frankensteinsmonster

**Approval.** Andrew, 2026-11-17 ~20:00 UTC: "migrate all 7 and maybe the 8th if it has not started on poplar yet".

**Why.** poplar's load had climbed back to ~550–590 on 336 threads (others ~470), which put the scheduler at its 64-thread floor: one new assembly lane, with assemblies taking 5–17 h there. frankensteinsmonster had just assembled the R2A sediment samples in 4.5–6 h each at 16 threads (two lanes), and its R2A work was ending.

**What moved and what stayed.**

| Samples | Where |
|---|---|
| R2_B_D1 (7.55 Gbp), R2_C_D1 (8.34), R2_C_D2 (6.41), R2_restored_C_black (8.58), R2_A_H2O (8.33), R2_B_H2O (6.14), R2_C_H2O (7.70), R2_restored_H2O (8.44) | **frankensteinsmonster**, same root as the R2A samples, local order 7–14 (after the R2A samples) |
| R2_A_D1, R2_A_D2 (assembling here since ~15:00 UTC), SF2_C_H2O (assembling), SF2_B_H2O (assembled) | poplar |
| R2_B_D2 (12.44 Gbp), R2A_A_H2O (9.65), R2A_C_H2O (10.20) | poplar: above 9 Gbp, so a 200g heap that frankensteinsmonster cannot give |

- **R2_B_D1** was the "8th": its merge and Tadpole had run here (17:21 UTC), but not its trim or assembly. It is redone there from its BBDuk reads with the same commands. Its `.done/merge:R2_B_D1` and `.done/tadpole:R2_B_D1` stay here, unused (no task in the new task set reads them), and so do its Tadpole outputs under `R2_B_D1/prep/` (`ecc_*`, 4.2 GB; not deleted).
- **All 12 R2 libraries stay here** as coverage for R2_A_D1, R2_A_D2 and R2_B_D2. Their 2026-side `checkm2:old`/`gunc:old` tasks for the 8 samples stay here too, as for R2A.

### Exactly what changed on poplar

All under `/scratch1/afreiburger/salterns_2026/rollout/`. `driver.py` (`507e1670…`), `rolib.py` (`89425e47…`) and `btlib.py` are **unchanged**: the staged-import machinery of the R2A migration covers this.

1. **Backups:** `scripts/config.json.pre_r2_20261117` (`e924414a…`), `scripts/make_config.py.pre_r2_20261117` (`7e9a2575…`).
2. **`scripts/make_config.py`** (now `1cb39882…`): `MIGRATED_R2` (the 8 samples), appended to `MIGRATED`.
3. **`scripts/config.json`** (now `d9d754cb…`): differs from the old one only by the 8 `import_from` (staged) blocks. `smoke_config.json` regenerated byte-identical.
4. **New helper scripts** (also in this repo's `scripts/`): `dryrun_migration.sh` (below) and `gen_task_scripts.py` (writes a sample's attempt-1 task scripts from a config without running anything).
5. **Driver restart**, `bash scripts/restart_driver.sh $ROOT scripts/config.json --rm-stop`:
   - 20:26:32 UTC STOP; acknowledged 20:26:53; driver 973270 stopped while sleeping.
   - **New PID 2919716**, started 20:27:01 with **710 tasks** (was 910); adopted 3/3 running tasks (`asm:R2_A_D1`, `asm:R2_A_D2`, `asm:SF2_C_H2O`); STOP removed 20:27:04.
   - The 8 samples show as `import_wait` until their `imports/<s>/READY` appears.
6. **`tmp/r2_md5/`**: poplar's own md5sum of each R2 library's `clean_*`, `map_*`, `read_counts.json` and `download.json` (nice 19, idle I/O), the reference for the copy.

### Checked before the restarts

- **Dry run, poplar** (`scripts/dryrun_migration.sh`: `driver.py plan` with the old and the new config on two copies of the markers):
  - 208 tasks removed, all `waiting` except the two R2_B_D1 markers above: 26 per sample (merge … checkm1); 8 `import:<s>` added.
  - 0 common tasks changed state, 0 sample-state changes; the 3 running tasks are in the new set.
  - Markers used by the old task set but not the new: only `merge:R2_B_D1`, `tadpole:R2_B_D1`.
- **Dry run, frankensteinsmonster** (same script, its markers):
  - 248 tasks added = 8 samples × 28 (merge … keep, push; 12 maps each) + `dl:`/`reads:` of the 12 R2 libraries; 0 removed.
  - 0 common tasks changed state; the R2A sample entries and the R2A libraries are byte-identical in the new config.
- **Task commands, poplar vs frankensteinsmonster, R2_C_D2** (`gen_task_scripts.py` on both hosts' dry-run copies; `~/Documents/salterns_rollout_local/task_script_diff_R2_C_D2{,.normalized}.txt`):
  - after normalising paths, the only differences are SPAdes `-t 24 -m 400` → `-t 16 -m 120`, bowtie2 `--threads 14` → `13`, `OMP_NUM_THREADS`, the reserves in the start lines, poplar's two 2026-side QC tasks and the local `push`;
  - Java heaps (`-Xms8g -Xmx150g`), k list, every BBTools, binner, DAS Tool, CheckM2/CheckM1/GUNC and KEEP parameter identical — the same result as R2A_B_D1's diff.

### What changed on frankensteinsmonster

Root `~/Documents/salterns_rollout_local/rollout/` (its README: "R2 migration").

1. **`scripts/make_config_local.py`** (this repo; backup `*.pre_r2_20261117`): `MIGRATED_R2`, the R2 libraries, R2 samples after the R2A ones (local order 7–14). **`local_config.json`** (`572f5738…`, installed as `rollout/scripts/config.json`; old one kept there as `config.json.pre_r2_20261117`, `5b015767…`): the R2A sample and library entries are byte-identical; it adds the 12 R2 libraries and the 8 R2 samples. The smoke configs regenerated byte-identical.
2. **Reads:** the 12 R2 libraries (`clean_*`, `map_*`, logs; 121 GB) copied 20:22–20:49 UTC with `logs/copy_inputs_r2.sh` (+ a second stream, `copy_inputs_r2_stream2.sh`; ~25–75 MB/s, poplar-I/O bound), then `logs/verify_inputs.sh <12 libs>` (`logs/verify_inputs_r2.log`): md5 = poplar's `tmp/r2_md5/` lists for all 12, pairs/bases = `manifest.tsv`, BBDuk chain counts link up; moved to `rollout/reads/<lib>/` with `.done/dl:<lib>` / `.done/reads:<lib>` provenance markers.
3. **Restart** with the new `~/Documents/salterns_rollout_local/restart_local.sh --rm-stop` (the local twin of `restart_driver.sh`: STOP, acknowledge, SIGTERM while sleeping, `launch_local.sh rollout`, adoption check, remove STOP):
   - 20:53:19 UTC STOP; driver 211753 stopped while sleeping at 20:53:30.
   - **New PID 863284**, scope **`salterns-r2a-local-20261117t205330z.scope`** (MemoryMax 220G), mem_guard **863316**; 442 tasks; adopted 1/1 (`asm:R2A_B_H2O`, which stays in the old scope `…20261116t171024z` with its old guard 211961 until it ends); STOP removed 20:53:34.
   - First scheduler cycle: all 14 ready tasks (6 R2A `checkm1:new`, 8 `merge:R2_*`) held by the memory gate (MemAvailable 138.5 GB, outstanding reserve 37 GB).

### ETA and risk (2026-11-17 21:00 UTC)

- **The local memory gate decides it.** Other sessions' jobs held ~80 GB, so with our last R2A assembly gone MemAvailable is about 171 GB: enough for the 6 waiting R2A CheckM1 runs (80 + 40 GB), not for a 150g-heap merge/Tadpole (150 + 40 = 190 GB). The gate is not lowered: poplar's 150g–200g Tadpole runs reached 144–157 GB RSS. The R2 preps start when those jobs release ~20 GB (MemAvailable ranged 118–207 GB over the previous 24 h).
- **Once the preps start:** 8 × ~0.5 h preps, then 8 assemblies two at a time (4.5–6 h each here) and ~2 h of binning/QC tail: ~28–33 h, pushes as each sample is kept.
- **poplar:** 3 running assemblies (SF2_C_H2O, R2_A_D1, R2_A_D2), 3 waiting (R2_B_D2, R2A_A_H2O, R2A_C_H2O) and SF2_B_H2O's binning, at one new lane while its load stays ~550: last assembly ~2026-11-20, `DONE` ~2026-11-20/21 if the local samples are back by then (was ~11-22/23 with all 11 here).
- **Fallback:** if no `merge:R2_*` has started locally by 2026-11-18 18:00 UTC, send some R2 samples back (remove them from `MIGRATED_R2` in both generators, regenerate, `dryrun_migration.sh`, restart both drivers) — with Andrew's OK.

## Local files

`scripts/` here holds copies of every script on poplar (`driver.py`, `rolib.py`, `btlib.py`, `make_config.py`, `config.json`, `smoke_config.json`, `status.sh`, `setup_gtdbtk.sh`, `smoke_prep.sh`, `upload_inputs.sh`, and since the R2A migration `restart_driver.sh`, `imports_README.txt` = poplar's `imports/README.txt`). `make_config.py` and `upload_inputs.sh` run locally.

For the R2A migration, `make_config_local.py` (run locally) writes `local_config.json` (the 7 samples on frankensteinsmonster, installed there as `rollout/scripts/config.json`), `local_smoke_config.json` and `smoke_import_config.json` (the poplar half of the import test, root `rollout/smoke_import/`). `driver.py`, `rolib.py` and `btlib.py` are md5-identical here, on poplar and on frankensteinsmonster.

For the R2 migration (2026-11-17): `dryrun_migration.sh` (old vs new config `plan` on copies of a root's markers) and `gen_task_scripts.py` (a sample's attempt-1 task scripts from a config); both are also in poplar's and frankensteinsmonster's `rollout/scripts/`. Backups of the edited files: `*.pre_r2_20261117`.

### Update 2026-11-18 16:30 UTC: load distribution decided

- Other sessions' jobs released memory at ~11:45 UTC; FM ran all 8 R2 read preparations 11:46–15:26 UTC (Tadpole peaks 84–112 GB) and started `asm:R2_B_D1` and `asm:R2_C_D1` at 15:26 (two lanes; 6 queued).
- poplar: 28/39 complete (all 7 R2A imported); its last assembly, R2A_C_H2O, was at K77 after 7.4 h; recent poplar assemblies took 7.8–11.4 h in its single lane (cap 64; others' load ~260).
- **Andrew asked for the distribution that finishes soonest. All 8 R2 samples stay on FM:** FM assembles in ~5–6 h with two lanes; poplar's one lane frees only when R2A_C_H2O ends and takes ~8–11 h per sample, so moving two of FM's last-round samples there would end ~5–8 h later than leaving them, and moving one saves nothing (FM still needs four rounds). FM cannot take a third lane within its gates (a third SPAdes needs 100 GB of room; ~90 GB at best).
- ETA: FM's last assembly ~2026-11-19 12:00–15:30 UTC, last push ~15:00–19:00 UTC; then import, keep and GTDB-Tk R232 on poplar → `DONE` ~2026-11-19 late to 2026-11-20.
- The one-off 18:05 UTC deadline check was cancelled (preps started).
