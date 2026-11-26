# Binning test: do Arm C assemblies give substantially better MAGs than the 2026 pipeline?

Started 2026-11-13 on poplar. Nothing here is committed, and nothing is written to KBase.

**Question.** The tadpole pilot (`../tadpole/RESULTS.md`) found that Arm C (BBTools merge + Tadpole correction + default-mode SPAdes 3.15.3) regenerates the Tringe-lab "(iii)" assembly of R1_A_D1 at 0.93× / 1.02× of its Mbp in contigs ≥ 2.5 kb / ≥ 10 kb. This test asks whether that assembly gives substantially better MAGs when binned exactly like the 2026 KBase pipeline, on R1_A_D1 and on SF2_C_D1 (the worst-recovered sample; its Arm C assembly was still running at the start).

If the decision is GO, the full salterns rollout (all 39 metagenomes through Arm C and this binning) starts on poplar. The EmilyKin and GROW tests wait until this is resolved.

**Status:** see "Status" at the end.

## Pre-registered decision rule

Fixed before any result was seen and implemented literally in `scripts/decide.py`, which writes `DECISION.json` (every metric, its threshold, pass/fail) and `REPORT.md`.

Per sample S in {R1_A_D1, SF2_C_D1}, compare **NEW** (the Arm C assembly binned like 2026) with **OLD** (the 2026 pipeline's DAS Tool output for the same unit). OLD is *all* DAS Tool bins of the 2026 unit, not only the KBase-CheckM-selected ones, and NEW is all DAS Tool bins of the new binning, so both sides are filtered by the same CheckM2 rule. Both are scored here with the same tools: CheckM2 (primary), CheckM1 `lineage_wf --reduced_tree` (secondary, to connect to the KBase CheckM values) and GUNC.

| | Criterion | Threshold |
|---|---|---|
| a | MQ+ MAGs (CheckM2 completeness ≥ 50, contamination < 10) | NEW ≥ 1.2 × OLD |
| b | HQ MAGs (≥ 90, < 5) | NEW ≥ 1.2 × OLD **or** NEW ≥ OLD + 5 |
| c | median N50 of MQ+ MAGs | NEW ≥ 1.5 × OLD |
| d | median CheckM2 contamination of MQ+ MAGs; GUNC pass fraction of MQ+ MAGs | NEW ≤ OLD + 1.0 point **and** NEW ≥ OLD − 0.05 |
| e | original species recovered: skani ANI ≥ 95 and max AF ≥ 50 between an MQ+ MAG and an original derep310 MAG of sample S | NEW ≥ OLD |

- S passes if a, b, c, d and e all pass. CheckM1-based a and b are reported too (informational).
- **GO** for the full salterns rollout only if **both** samples pass.
- Only R1_A_D1 passes: **NO_GO**; the SF2-type decision needs Andrew.
- SF2_C_D1's Arm C assembly fails, or has no done marker by **2026-11-14T06:18:00Z** (24 h after the task was issued): decide on R1_A_D1 alone, as **R1_ONLY_PASS** (or NO_GO if R1_A_D1 fails), never GO.

Conventions that the rule text does not fix (stated in `DECISION.json` as well):

- A median or fraction over an empty MQ+ set is undefined. An undefined OLD with a defined NEW passes; an undefined NEW fails.
- GUNC pass means `pass.GUNC == True` in `GUNC.progenomes_2.1.maxCSS_level.tsv`. A genome missing from GUNC's output counts as not passing; one missing from the CheckM2 report is not MQ+.
- A sample whose processing failed (a required step failed after its retries) is INCOMPLETE, and so is the decision (needs Andrew). The exception is the SF2 assembly failure or timeout above.
- Originals for e): R1_A_D1 has 11 derep310 MAGs; SF2_C_D1 has one (concoct_out.13, which is below MQ in the 2026 re-QC: 95.6 / 15.1).
- Added at review, before any decision (see "Review" below):
  - Every threshold is evaluated in exact rational arithmetic on the values as printed (CheckM2 values as decimal strings, GUNC pass fractions as ratios, N50 medians). A value exactly at a threshold therefore gets the inclusive inequality the rule states.
  - GUNC rows are joined to genomes on GUNC's own name, the file name cut at its first `.fa`. This changes only originals named `*.fa_assembly`, which are informational.

`scripts/decide.py` was unit-tested on synthetic tables before launch for every branch: GO, NO_GO (R1 passes and SF2 fails), R1_ONLY_PASS (SF2 timeout), NO_GO (SF2 assembly failed and R1 fails), INCOMPLETE, empty OLD or NEW sets, and a crash (it then writes `DECISION.json` with `decision: ERROR` and exits 1; the driver retries once and still writes the done marker).

## What runs, and how it matches the 2026 pipeline

One driver (`scripts/driver.py`, stdlib Python) runs every step on poplar as a detached bash task with a done marker. The two samples are two lanes that run side by side: R1_A_D1 starts at once, and SF2_C_D1 starts when its Arm C done marker appears.

| Step | Here | 2026 KBase pipeline (`../../pipeline.py`, `../../local/`) | Match? |
|---|---|---|---|
| Assembly | Arm C `scaffolds.fasta` (tadpole pilot), sequences ≥ 1,000 bp (length with N's) | metaSPAdes 3.15.3 `scaffolds.fasta` saved with `min_contig_length 1000` | the variable under test |
| Reads for coverage | the 9 libraries of the sample's site (R1: R1_{A,B,C}_{D1,D2,H2O}; SF2 likewise), from SRA (NCBI ODP mirror) | the same 9 libraries per site (`coverage_group`) | same libraries |
| Read prep for coverage | Arm C's BBDuk steps without merging or Tadpole: adapter pass, artifact/phiX pass, then `qtrim=r trimq=10 minlen=45` | Trimmomatic 0.39 PE | changed with the assembly recipe, as the task specifies |
| Mapping | bowtie2 2.5.5 defaults (end-to-end, `--sensitive`), `-1/-2`, unmapped dropped (`samtools view -F 0x04`), sorted and indexed BAMs; one BAM per library, shared by all binners | CONCOCT and MaxBin2 apps: bowtie2 defaults (2.3.0 / 2.2.3); MetaBAT2 app: BBMap 37.90 `fast` | **deviation for MetaBAT2** (bowtie2 instead of BBMap fast), by design: one mapper for all |
| MetaBAT2 depth | `jgi_summarize_bam_contig_depths --minContigLength 2500 --minContigDepth 1` per BAM, merged with KBase's own logic (`combineDepthTablesUtil.py`: zero-fill, totalAvgDepth = sum, assembly order; `btlib.py merge-depth`) | the same | same |
| MetaBAT2 | **2.12.1** (the KBase app's binary version): `--minContig 2500 --maxEdges 500 --unbinned`; unbinned/lowDepth/tooShort not used | `metabat_wrapper.sh`: the same flags | same version and flags |
| MaxBin2 | **2.2.7**: `run_MaxBin.pl -abund_list <9 files> -thread 16 -prob_threshold 0.8 -markerset 40 -min_contig_length 1000` | 2.2.4 with `-reads_list` (MaxBin maps each library itself with bowtie2 2.2.3) and `-thread 32` | **deviation**: abundance = per-library mean depth from `jgi_summarize --minContigLength 1000`, not MaxBin's internal read mapping; 2.2.4 → 2.2.7 is crash fixes only |
| CONCOCT | **1.1.0**: `cut_up_fasta.py -c 10000 -o 0 --merge_last -b temp.bed`; `concoct_coverage_table.py temp.bed <9 BAMs>`; `concoct -l 2500 -k 4 -c 400 -i 500 --total_percentage_pca 90 -t 16`; `merge_cutup_clustering.py clustering_gt2500.csv`; `extract_fasta_bins.py` | `ConcoctUtil.py` (kb_concoct 1.3.4, CONCOCT 1.1.0 + 3 later commits): the same commands and values | same |
| DAS Tool | **1.1.7** (diamond 2.2.8, prodigal 2.6.3, as `~/Documents/dastool_env`): `--search_engine diamond --score_threshold 0.5 --duplicate_penalty 0.6 --megabin_penalty 0.5 --write_bin_evals --write_bins -t 16`, labels `metabat2,maxbin2,concoct` | `local/dastool_worker.py`: the same, `-t 8` (thread count does not change the result: `../../binning_experiments/dastool_min_contig/`) | same |
| Empty binners | MaxBin2 "dataset cannot be binned" → empty, DAS Tool runs on the others; DAS Tool "no bins with bin-score" → NEW has 0 bins | the same in `pipeline.py` / `dastool_worker.py` | same |

QC, identical for NEW, OLD and the sample's originals:

- **CheckM2 1.1.0** `predict -x fa --database_path .../CheckM2_database/uniref100.KO.1.dmnd`.
- **CheckM1 1.1.3** (new env `bt_checkm1`; the EmilyKin env's 1.2.3 cannot import, see the smoke test) `lineage_wf --reduced_tree -x fa -t 16 --pplacer_threads 4 --tab_table`, data `checkm_data_2015_01_16` (md5 631012fa598c43fdeb88c619ad282c4d) through `CHECKM_DATA_PATH`. KBase used CheckM 1.1.2 (kb_Msuite) on the same reduced tree; `DECISION.json` reports how the OLD bins' KBase values compare with CheckM1 here.
- **GUNC 1.1.1** `run --file_suffix .fa` with the progenomes_2.1 database (the `gunc download_db` default).
- **skani 0.3.2** `triangle -E` over NEW + OLD + the sample's originals (derep310 and nonderep140).
- N50, contig count and size from the FASTA files (`decide.py`).

CheckM2 and skani come read-only from the EmilyKin `bin_qc` env, which is not modified (`/scratch1/afreiburger/emilykin/processed/.snakemake_envs/8c703f8abdc13cdc1ef374a2bcf37f14_`).

**OLD** = all DAS Tool bins of the 2026 unit, copied from `../../work/dastool/<unit>/attempt1-*/out/<unit>_DASTool_bins/` (R1_A_D1 60, SF2_C_D1 6; md5 lists in `inputs/old/<unit>.md5` on poplar). **Originals** = the sample's MAGs from `salterns/MAG_contigs/` (R1_A_D1: 11 derep310 + 12 nonderep140; SF2_C_D1: 1 + 4), md5-checked against `index.tsv`.

## Exact commands

Each task is a bash file the driver writes once per attempt (`scripts/tasks/<task>.<attempt>.sh`, never edited afterwards) with `set -euo pipefail`, `OMP/MKL/OPENBLAS_NUM_THREADS` = the task's thread count and its own `TMPDIR`. Everything inherits `nice -n 10` from the driver. `$ROOT` = `/scratch1/afreiburger/salterns_2026/binning_test`.

```
# dl:<lib>  (16 libraries; the two pilot samples reuse the pilot's verified reads)
curl -sS --retry 5 --retry-delay 20 -C - -o <srr>.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/<srr>/<srr>
vdb-validate <srr>.sra                                    # sra-tools 3.4.1
fasterq-dump --split-files -e 6 -t fqd_tmp -O . <srr>.sra
awk 'NR%4==3{print "+"; next} {print}' > raw_<m>.fastq.tmp (counts pairs and bases in the same pass) && mv
#   pairs of both mates must equal manifest.tsv 'spots', total bases 'bases', line counts 4 x spots
# reads:<lib>  (BBTools 39.81, t=16; pilot samples: only the 3rd pass, on tadpole_pilot/<sample>/qc/clean_*)
bbduk.sh -Xmx16g in1=raw_1 in2=raw_2 out1=trim_1 out2=trim_2 ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo ftm=5 maq=5 minlen=45 mlf=0.333 maxns=0 ordered
bbduk.sh -Xmx16g in1=trim_1 in2=trim_2 out1=clean_1 out2=clean_2 ref=artifacts,phix k=31 hdist=1 ordered
bbduk.sh -Xmx16g in1=clean_1 in2=clean_2 out1=map_1 out2=map_2 qtrim=r trimq=10 minlen=45 ordered
#   read_counts.json: Input/Result reads and bases of each pass; chained counts are asserted
# asm:<S>
btlib.py filter-fasta tadpole_pilot/<S>/armC/spades/scaffolds.fasta 1000 <S>/asm/contigs.fa
# idx:<S> / map:<S>:<lib>  (bowtie2 2.5.5, samtools 1.25)
bowtie2-build --threads 16 --seed 1 asm/contigs.fa bt2/contigs
bowtie2 -x bt2/contigs -1 map_1.fq.gz -2 map_2.fq.gz --threads 14 | samtools view -F 0x04 -u - | samtools sort -@ 2 -m 2G -o bam/<lib>.bam -
samtools index bam/<lib>.bam
# depth:<S>  (jgi_summarize_bam_contig_depths from MetaBAT2 2.12.1)
jgi_summarize_bam_contig_depths --outputDepth depth/<lib>.d2500.txt --minContigLength 2500 --minContigDepth 1 bam/<lib>.bam   # x 9
jgi_summarize_bam_contig_depths --outputDepth depth/<lib>.d1000.txt --minContigLength 1000 --minContigDepth 1 bam/<lib>.bam   # x 9
btlib.py merge-depth asm/contigs.fa depth/metabat_depth.txt depth/*.d2500.txt          # KBase combineDepthTablesUtil logic
btlib.py abund depth/<lib>.d1000.txt depth/abund_<lib>.txt                              # MaxBin2 abundance, x 9
# metabat2:<S>
metabat2 --inFile asm/contigs.fa --outFile metabat2/bins/bin --abdFile depth/metabat_depth.txt --numThreads 16 --minContig 2500 --maxEdges 500 --unbinned
# maxbin2:<S>
run_MaxBin.pl -contig asm/contigs.fa -out maxbin2/out/mb -abund_list depth/abund_list.txt -thread 16 -prob_threshold 0.8 -markerset 40 -min_contig_length 1000
# concoct:<S>
cut_up_fasta.py asm/contigs.fa -c 10000 -o 0 --merge_last -b temp.bed > split_contigs.fa
concoct_coverage_table.py temp.bed bam/<9 libs>.bam > coverage_table.tsv
concoct --composition_file split_contigs.fa -l 2500 -b out/ --coverage_file coverage_table.tsv -t 16 -k 4 -c 400 -i 500 --total_percentage_pca 90
merge_cutup_clustering.py out/clustering_gt2500.csv > out/clustering_merged.csv
extract_fasta_bins.py asm/contigs.fa out/clustering_merged.csv --output_path bins
# dastool:<S>  (contig2bin labels metabat2.bin.NNN / maxbin2.Bin.NNN / concoct.bin.NNN, as in 2026)
DAS_Tool -i metabat2/contigs2bin.tsv,maxbin2/contigs2bin.tsv,concoct/contigs2bin.tsv -l metabat2,maxbin2,concoct -c asm/contigs.fa -o dastool/out/<S> \
  --search_engine diamond --score_threshold 0.5 --duplicate_penalty 0.6 --megabin_penalty 0.5 --write_bin_evals --write_bins -t 16
# checkm2|checkm1|gunc:<set>:<S>   set = new (dastool/out/<S>_DASTool_bins), old (inputs/old/<S>), orig (inputs/orig/<S>)
checkm2 predict --input <dir> -x fa --output-directory qc/<set>/<S>/checkm2/out --threads 16 --database_path <CheckM2_database/uniref100.KO.1.dmnd>
CHECKM_DATA_PATH=refdata/checkm1 checkm lineage_wf --reduced_tree -x fa -t 16 --pplacer_threads 4 --tab_table -f qc/<set>/<S>/checkm1/qa.tsv <dir> qc/<set>/<S>/checkm1/out
gunc run --input_dir <dir> --file_suffix .fa --db_file refdata/gunc/gunc_db_progenomes2.1.dmnd --threads 16 --out_dir qc/<set>/<S>/gunc
# skani:<S>
skani triangle -l genomes.txt -t 8 -E -o qc/skani/<S>/skani_triangle.tsv          # NEW + OLD + originals of S
# decide  (when both samples are complete, incomplete, or SF2's assembly failed / timed out)
python -I scripts/decide.py --config scripts/config.json --state STATE.json --out $ROOT   # -> DECISION.json, REPORT.md, mags.tsv
#   then the driver writes $ROOT/DONE
```

**Scheduling.** At most 96 threads in all (`thread_cap`), counting each task's thread setting; at most 4 downloads and 6 read preparations at once; OLD/original QC (priority ≥ 60) never holds more than 32 threads, so the NEW chain keeps most of the budget. Failed tasks are retried once after ≥ 5 min (downloads twice). A restarted driver adopts tasks that are still running (`.run/<task>.json` + `/proc` check) and skips finished ones (`.done/<task>`).

**The SF2 wait.** `asm:SF2_C_D1` starts only when `tadpole_pilot/SF2_C_D1/.done/C1_spades_default` and `armC/spades/scaffolds.fasta` exist. The lane is abandoned (`.state/SF2_C_D1.abandoned`) as `asm_timeout` if the marker is missing at 2026-11-14T06:18:00Z, or as `asm_failed` if SPAdes (PID 3634657) and its two `run_sample.sh` parents (3634245, 3634387) have all been gone for 15 min without the marker. The driver only watches these processes; it never signals them.

## Where everything is on poplar

Root: `/scratch1/afreiburger/salterns_2026/binning_test/`

| Path | What |
|---|---|
| `scripts/` | copies of `scripts/` here: `driver.py`, `decide.py`, `test_decide.py` (review tests of decide.py), `btlib.py`, `config.json`, `make_config.py`, `status.sh`, `setup_envs.sh`, `setup_checkm1.sh`, `fetch_gunc_db.sh`, `smoke_prep.sh`, `smoke_config.json`; `scripts/tasks/` holds every generated task script |
| `inputs/old/<S>/` | OLD: the 2026 DAS Tool bins (`.fa`), `inputs/old/<S>.md5` and `.source.txt` |
| `inputs/orig/<S>/` | the sample's original MAGs and `originals.tsv` (mag, population, file, md5) |
| `inputs/kbase_checkm/` | the 2026 KBase CheckM tables (`checkm_<S>.tsv` = all DAS Tool bins; `selected_<S>.tsv`) |
| `reads/<lib>/` | `map_{1,2}.fq.gz` (mapping input), `clean_{1,2}.fq.gz` (BBDuk output before the quality trim; kept for a rollout), `read_counts.json`, `download.json`, BBDuk logs and stats, `vdb-validate.log` |
| `<S>/asm/` | `contigs.fa` (Arm C scaffolds ≥ 1 kb), `asm_stats.json` |
| `<S>/bt2/`, `<S>/bam/` | bowtie2 index; one sorted, indexed BAM per library plus its bowtie2 log (alignment rate) |
| `<S>/depth/` | per-BAM jgi tables, `metabat_depth.txt`, `abund_<lib>.txt`, `abund_list.txt` |
| `<S>/{metabat2,maxbin2,concoct}/` | binner outputs and `contigs2bin.tsv` |
| `<S>/dastool/out/` | DAS Tool 1.1.7 output; NEW bins in `<S>_DASTool_bins/` |
| `qc/{new,old,orig}/<S>/{checkm2,checkm1,gunc}/` | QC outputs (`quality_report.tsv`, `qa.tsv`, `GUNC.progenomes_2.1.maxCSS_level.tsv`) |
| `qc/skani/<S>/` | `genomes.txt`, `skani_triangle.tsv` |
| `refdata/checkm1/`, `refdata/gunc/` | CheckM1 data, GUNC progenomes_2.1 db (each with `.ok`) |
| `STATE.json` | driver state, rewritten every 20 s: per-sample state, per-task state, running tasks, threads |
| `.done/`, `.failed/`, `.run/`, `.rc/`, `.state/` | markers: finished tasks, permanent failures, running attempts, exit codes, attempts / abandon records |
| `logs/driver.log`, `logs/tasks/<task>.<attempt>.log` | driver log; one log per task attempt |
| `DECISION.json`, `REPORT.md`, `mags.tsv` | the decision, its report, per-MAG metrics (NEW, OLD, originals) |
| **`DONE`** | **done marker**, written after `DECISION.json`: `<utc> decision=<GO / NO_GO / R1_ONLY_PASS / INCOMPLETE / ERROR>` |
| `smoke/` | the smoke test (same layout) |

Envs (new, under `/scratch1/afreiburger/salterns_2026/envs/`, each with `bt_<name>.create.log` and `bt_<name>.list.txt`; mamba 2.6.0, conda-forge + bioconda, package cache `envs/pkgs`):

| Env | Key packages |
|---|---|
| `bt_map` | bowtie2 2.5.5, samtools 1.25, seqkit 2.14.0, pigz 2.8, python 3.12.16 (runs the driver) |
| `bt_metabat2` | metabat2 2.12.1 (with its jgi_summarize_bam_contig_depths) |
| `bt_maxbin2` | maxbin2 2.2.7, hmmer 3.4, FragGeneScan 1.32 |
| `bt_concoct` | concoct 1.1.0, numpy 1.26.4, pandas 3.0.6, scikit-learn 1.9.1, samtools 1.25 |
| `bt_dastool` | das_tool 1.1.7, diamond 2.2.8, prodigal 2.6.3, R 4.4.3 |
| `bt_gunc` | gunc 1.1.1 |
| `bt_checkm1` | checkm-genome 1.1.3, python 3.8.20, hmmer 3.4, prodigal 2.6.3, pplacer 1.1.alpha19, setuptools 69.5.1 (`setup_checkm1.sh`) |

Also used, read-only: BBTools 39.81 (`/scratch1/afreiburger/emilykin/processed/reanalysis_norm/env_bbtools`), sra-tools 3.4.1 (`/scratch1/afreiburger/metag-hybrid/env/miniforge/envs/sra`), and the EmilyKin `bin_qc` env (CheckM2 1.1.0 and skani 0.3.2; its checkm-genome 1.2.3 is not used) with the CheckM2 database in `/scratch1/afreiburger/emilykin/processed/refdata/checkm2/`.

## How to check status

```
ssh poplar 'bash /scratch1/afreiburger/salterns_2026/binning_test/scripts/status.sh'
ssh poplar 'cat /scratch1/afreiburger/salterns_2026/binning_test/DONE'          # exists only when finished
rsync -a -e "ssh -o BatchMode=yes" poplar:/scratch1/afreiburger/salterns_2026/binning_test/{DECISION.json,REPORT.md,mags.tsv,STATE.json} results/
```

`status.sh` prints the host load, the driver PID and whether it is alive, `DONE`, the SF2_C_D1 Arm C assembly state (marker, SPAdes PID, last log line), per-sample states, task counts, failed/blocked tasks, running tasks and the last driver log lines. `status.sh <root>` works for the smoke test too.

**Stopping and restarting.** `touch $ROOT/STOP` makes the driver start nothing new. To restart, kill the driver PID (`driver.pid`) while running tasks keep going (they are detached), then start it again with the launch command below: it adopts the running tasks and skips finished ones. Never run two drivers on the same root. Delete `.failed/<task>` to give a permanently failed task another try.

## Smoke test (2026-11-13, passed)

Run with the same driver code on `smoke/` (`scripts/smoke_config.json`, `scripts/smoke_prep.sh`):

- **Inputs:** sample SMOKE_A = the 215 longest R1_A_D1 Arm C scaffolds (30.05 Mbp) with 2 libraries (4 M R1_A_D1 + 2 M SF2_C_D1 cleaned read pairs from the pilot), 6 OLD bins and 4 originals of R1_A_D1. Sample SMOKE_B never gets an assembly and has empty OLD/original sets, with a 20-minute deadline.
- **What ran and passed:** pilot-type read prep (quality-trim pass, counts), asm filter, bowtie2-build, 2 maps, per-BAM jgi depth (both thresholds) and the KBase-style merge, MetaBAT2, MaxBin2 (`-abund_list`), CONCOCT (pandas 3 / scikit-learn 1.9 work), DAS Tool (7 bins from 3 binners), CheckM2, CheckM1, GUNC, skani, the EMPTY paths (SMOKE_B's empty sets), the SF2 timeout path (SMOKE_B → `asm_timeout`), `decide.py` on real outputs, the `DONE` marker, and driver restarts that adopted running tasks.
- **Result:** `DONE` = `decision=NO_GO` (SMOKE_B timed out → decided on SMOKE_A alone, which failed b, d and e). The numbers mean nothing with 2 libraries on 30 Mbp. One useful check: **CheckM1 1.1.3 here reproduces the KBase CheckM 1.1.2 values** of the 6 OLD bins to a median 0.003 points of completeness and 0.001 of contamination.
- **Separately tested at full scale before launch:** the download step (4 R1 runs: pair and base counts equal `manifest.tsv`, e.g. R1_A_D2 29,528,975 pairs / 8,917,750,450 bases) and the non-pilot read prep (R1_A_H2O: 54.35 M → 52.87 M reads, chained counts asserted).

**Problems the smoke test found, all fixed before the real launch:**

1. **CheckM1 in the EmilyKin env cannot import** (`No module named 'pkg_resources'`: setuptools 82). A new env `bt_checkm1` has checkm-genome **1.1.3** (python 3.8, setuptools 69.5.1; closest to KBase's 1.1.2), so the EmilyKin env stays untouched.
2. **CheckM2's `mp.Manager()` died with EOFError** because the task TMPDIR under `/scratch1/.../binning_test/tmp/...` made the AF_UNIX socket path 114 characters (limit 107). Every task now uses a short `TMPDIR=/tmp/bt_<hash>` (26 GB free), removed on exit; big temporary files go to explicit `-T` / `-t` directories.
3. **`zcat | head` under `pipefail` dies of SIGPIPE** (smoke prep); `seqkit head` is used instead, and the task bodies were checked for the same pattern (`cat DECISION.json | head` → `head -c`).
4. **The GUNC database download from zenodo ran at 1–3 MB/s.** It is fetched from gunc's second (EMBL) URL in 8 byte ranges instead; the mirror answers 403 to some concurrent requests, so each range is retried until it has its exact size. Both gunc 1.1.1 md5s match (gz `bc93a855…`, dmnd `447c9330…`).
5. GUNC and CheckM1 tasks now wait for their reference data and env markers (`refdata/gunc/.ok`, `refdata/checkm1/.ok`, `envs/.bt_setup/checkm1.ok`) instead of failing early.

## Review (2026-11-13, 07:45–14:10 UTC, before the decision)

This is an independent review of the launched test, against `../../pipeline.py`, the KBase app sources and the ledger.

**Timing.** The review began at 07:45 UTC. By the time the decide.py fix was ready (13:59 UTC), every task except two SF2 QC tasks and `decide` had finished. There were no failures, no retries and no abandoned lane. The fix was deployed before `decide` ran. To check the outputs first, `STOP` held new launches from 13:59:31 to 14:00:23 UTC. Nothing was waiting to start then except `decide`.

**1. Binner and mapping parameters vs the 2026 pipeline.**
- **Verified equal.** These were checked against the 2026 ledger parameters and the app sources (kbase-metabat at both the release caf4e3c0 and the beta bd317234 that 2026 actually ran; kb_concoct 936c9bfe; kb_maxbin 3574f0b2):
  - MetaBAT2 2.12.1 with `--minContig 2500 --maxEdges 500 --unbinned`. `metabat_wrapper.sh` is identical in release and beta.
  - Per-BAM `jgi_summarize --minContigLength 2500 --minContigDepth 1`, merged by `combineDepthTablesUtil.py` (identical in release and beta). `btlib.py merge-depth` reproduces its logic, and the merged table rows equal the contigs ≥ 2.5 kb (R1 31,890; SF2 7,965).
  - CONCOCT: `bowtie2 -1/-2` defaults, `samtools view -F 0x04`, cut-up 10000/0/`--merge_last`, `-l 2500 -k 4 -c 400 -i 500 --total_percentage_pca 90`, no `--no_cov_normalization`/`--no_total_coverage`.
  - MaxBin2: `-prob_threshold 0.8 -markerset 40 -min_contig_length 1000`, no plotmarker.
  - DAS Tool 1.1.7, diamond 2.2.8 and prodigal 2.6.3: the same packages as `~/Documents/dastool_env`, with the same flags.
- **Unavoidable differences.** All of them apply to NEW only, so they are confounded with the assembly. Expected direction:

  | Difference | Expected direction |
  |---|---|
  | MetaBAT2 coverage from bowtie2 2.5.5 end-to-end BAMs, not BBMap 37.90 `fast`. The same jgi 97 %-identity and edge filters follow. | small, not predictable: neutral |
  | MaxBin2 abundance is jgi mean depth (97 % identity, edges excluded), not MaxBin 2.2.4's own bowtie2 2.2.3 mapping plus read counting. MaxBin 2.2.7 vs 2.2.4 is crash fixes only. | small; if anything slightly toward NEW (cleaner per-library abundance) |
  | Coverage reads: BBDuk (Q10 right trim, minlen 45) instead of Trimmomatic (SLIDINGWINDOW 4:15, MINLEN 50). Slightly more bases are kept, evenly across contigs. | neutral |
  | CONCOCT's bowtie2 is 2.5.5, not Debian 2.3.0; the command is the same. | neutral |
  | Thread counts: MaxBin 16 vs 32, DAS Tool 16 vs 8. | none |
  | MetaBAT2 has no fixed seed on either side. | run-to-run noise on both |

  None of these can plausibly produce the R1 gaps seen (MQ+ 62 vs 47, HQ 27 vs 11). Criterion c (MAG N50 15.4 vs 4.7 kb) reflects contig length, which the binning cannot create.

**2. Scoring is identical for OLD and NEW.**
- One task template (`driver.py body()`) is used for NEW, OLD and the originals: CheckM2, CheckM1, GUNC and skani. Only the input directory differs.
- Row counts equal bin counts in every QC table: R1 NEW 71, OLD 60, originals 23; SF2 NEW 11, OLD 6, originals 5.
- CheckM1 1.1.3 here reproduces KBase CheckM 1.1.2 on all 60 R1 OLD bins: median |Δ| 0.003 completeness and 0.002 contamination, and MQ+ 45 = 45.
- OLD and NEW bin names contain no inner `.fa`, so GUNC's name truncation does not touch them.

**3. decide.py vs the pre-registered rule: one class of error, fixed.**
- **The error.** decide.py compared floats. At an exact tie, the inclusive inequalities in d came out false:
  - `1.36 <= 0.36 + 1.0` and `1.57 <= 0.57 + 1.0` are False in float;
  - `3/20 >= 1/5 - 0.05` (0.15 vs 0.15000000000000002) is False too.
  - CheckM2 prints two decimals and SF2's OLD set is tiny (fractions like k/5), so such ties are reachable.
- **The fix.** All thresholds now use exact `Fraction` arithmetic:
  - `5·NEW >= 6·OLD` for a and b (and CheckM1 a/b);
  - `Fraction(3,2)` for c;
  - decimal-exact contamination medians for d1;
  - GUNC ratios for d2.
- **Second fix.** GUNC rows are joined on GUNC's own name (`basename.split('.fa')[0]`). The originals named `*.fa_assembly` had no GUNC value before; this is informational only.
- Both conventions were added to the decide.py docstring, so they land in `DECISION.json`. The rule text itself is unchanged.
- **Tests.** `scripts/test_decide.py` builds 38 synthetic roots and compares decide.py with an independent exact oracle: every criterion at its tie and just past it, empty and undefined sets, missing CheckM2/GUNC rows, the skani ANI/AF/MQ+/nonderep140 rules, and all decision branches.
  - The pre-review decide.py failed 3 of them (the d1/d2 ties, each flipping GO to NO_GO). The fixed version passes 38/38, locally and on poplar.
  - A NaN still gives `decision: ERROR`.
  - On the real smoke outputs, the fixed version gives the same decision and criteria as before; only the GUNC value of the `.fa_assembly` original changes.
- **Deployment.** The fix was deployed on poplar by `scp` to `decide.py.tmp` and then `mv`, at 13:59 UTC before `decide` ran. Checksums: `d1c47ac2…` before review (kept as `logs/decide.py.prereview`), `b2f78405…` after. Local and poplar `scripts/` are byte-identical.

**4. Reads, samples and originals: correct.**
- Coverage libraries are R1 ×9 and SF2 ×9, the same samples and SRRs as `manifest.tsv` and the 2026 ledger `reads_list`s. All 18 were Trimmomatic outputs in 2026.
- Every download matched the manifest pair and base counts.
- OLD is all DAS Tool bins of the only 2026 DAS Tool attempt of each unit: R1 60 (`275645/513/1`) and SF2 6 (`275645/1168/1`). Both are md5-identical to `../../work/dastool/<unit>/attempt1-*/out/`.
- Originals are R1 11 derep310 + 12 nonderep140 and SF2 1 + 4, with md5s equal to `MAG_contigs/index.tsv`. Only derep310 enter e.
- The 2026 assemblies of both units were metaSPAdes, not the MEGAHIT fallback.

**5. Budget and isolation: within limits.**
- At 07:45 UTC, 88 of 96 threads were in use, and every task ran at nice 10.
- A map task really uses about 17 threads (bowtie2 14, samtools view 1, sort 2) against its cost of 16. The true peak can therefore exceed 96 by about 6 threads while 6 maps run. This is minor and was left as is.
- The DRAM driver (PID 4007531) and the SF2 SPAdes (PID 3634657) were never signalled. The driver only uses `kill(pid, 0)` and reads `/proc/<pid>/cmdline`.
- `/scratch1` has 9.8 TB free; `/tmp` has 26 GB.

## Status

**Finished 2026-11-13 14:04:25 UTC: `DONE` = `decision=NO_GO`.** The driver (PID 3607621) exited normally. All 87 tasks finished on their first attempt, and the decision was made by the reviewed decide.py. `DECISION.json`, `REPORT.md`, `mags.tsv`, `STATE.json` and `DONE` are copied to `results/`.

| | MQ+ NEW / OLD | HQ | median MAG N50 (bp) | median contamination | GUNC pass (MQ+) | originals recovered | result |
|---|---|---|---|---|---|---|---|
| R1_A_D1 | 62 / 47 | 27 / 11 | 15,416 / 4,656 | 3.34 / 3.48 | 0.87 / 0.79 | 10 / 8 of 11 | **PASS** (a–e) |
| SF2_C_D1 | 7 / 4 | 2 / 0 | 9,133 / 4,852 | **4.67 / 1.70** | **0.86 (6/7) / 1.00 (4/4)** | 0 / 0 of 1 | **FAIL** on d only (d1 threshold 2.70, d2 threshold 0.95) |

- **Rule outcome.** Only R1_A_D1 passes, so the rule gives **NO_GO**. The full rollout is not started, and the SF2-type decision is Andrew's.
- **Not a tie or an artifact.**
  - The pre-review decide.py gives the same NO_GO on the same outputs.
  - Every QC table has one row per bin.
  - No MQ+ MAG lacks a GUNC row.
- **Context for that decision (not part of the rule).**
  - SF2's d failure rests on 7 vs 4 MQ+ MAGs.
  - OLD's MQ+ MAGs are mostly about 64 % complete: median 64.0 vs NEW 88.2. Less complete bins tend to show lower CheckM2 contamination.
  - NEW has twice the MQ+ Mbp (23.8 vs 11.3).
  - The one GUNC failure is `metabat2.bin.003_sub` (63.6 / 5.85).
  - Neither set recovers SF2's single derep310 original, which re-scores here at 92.3 / 10.7.

**Launched 2026-11-13 07:35:56 UTC (01:35 CST) on poplar. Driver PID 3607621** (since 07:40:41 UTC; `$ROOT/driver.pid`; log `logs/driver.log`).

```
cd /scratch1/afreiburger/salterns_2026/binning_test
nohup setsid nice -n 10 /scratch1/afreiburger/salterns_2026/envs/bt_map/bin/python -I scripts/driver.py run \
    --config scripts/config.json >> logs/driver.log 2>&1 < /dev/null &
```

How it got there:

- 06:45 UTC: a download-only pre-run of the same driver (`--only-kinds dl --cap 30`, log `logs/driver_prerun_dl.log`), started while the smoke test was being built.
- 07:09 UTC: widened to `dl,reads,asm,idx --cap 80`.
- 07:35 UTC: replaced by the full driver (PID 3544261), which adopted the 7 running download and read-prep tasks.
- 07:40 UTC: restarted once more (PID 3607621) with two hardening changes, adopting the 8 running tasks (4 maps, 4 downloads):
  - each scheduling cycle is wrapped in an exception guard (an NFS error is logged as `cycle error` instead of killing the driver);
  - an exit-code file that is still empty is read again on the next cycle rather than counted as a failure.

At launch:

- R1_A_D1's assembly filter and bowtie2 index were done (Arm C ≥ 1 kb: 129,813 sequences, 393.9 Mbp), its reads and SF2_C_D1's were prepared, 4 R1 downloads were verified, and the first R1 mapping had started.
- SF2_C_D1's Arm C SPAdes (PID 3634657) had no done marker yet. *Correction at review:* it was in **K77**, not K127. `spades-core` was running `K77/configs/config.info`; K33 took 2 h 18 min and K55 2 h 49 min. It finished anyway at 12:52 UTC (`.done/C1_spades_default`), well before the deadline.

**Expected timeline (at launch; actual: SF2 Arm C marker at 12:52 UTC, `DONE` at 14:04 UTC).** These are estimates at nice 10 on a host at load 430–550.

- Each read prep took about 20 min (16 threads).
- R1_A_D1's 9 mappings, binning, DAS Tool and QC should finish around 13:00–16:00 UTC on 2026-11-13.
- SF2_C_D1 adds about 5–7 h after its assembly finishes, so the decision (`DONE`) is expected about **2026-11-13 18:00–23:00 UTC**.
- If the SF2 assembly never finishes, the decision comes at the deadline, **2026-11-14 06:18 UTC**, as R1_ONLY_PASS or NO_GO.

**Watch points not reached by the smoke test** (both resolved: MetaBAT2 binned 108 / 17 bins with no abundance complaint, its depth rows equal the contigs ≥ 2.5 kb, and no QC task failed)**:**

- **MetaBAT2 on a mixed-length depth table.** The smoke assembly had no contig under 2.5 kb, so `metabat2 --minContig 2500` with the merged jgi table of an assembly that also has 1–2.5 kb contigs was not exercised. KBase ran the same binary and flags, so it should pass. When `depth:R1_A_D1` finishes, check `R1_A_D1/metabat2/metabat2.log` for an abundance-file complaint, and that the rows of `depth/metabat_depth.txt` match the number of contigs ≥ 2.5 kb.
- **Space in `/tmp`.** All QC temporary files now live in `/tmp` (26 GB free, shared). If a CheckM2, CheckM1 or GUNC task fails after the binning chain succeeded, check `df -h /tmp` first.

**When `DONE` exists:** copy `DECISION.json`, `REPORT.md` and `mags.tsv` here (see the rsync above) and act on the decision:

- **GO:** start the full 39-sample salterns rollout on poplar (Arm C assembly, then this binning).
- **NO_GO / R1_ONLY_PASS / INCOMPLETE / ERROR:** report to Andrew. The SF2-type decision is his.
