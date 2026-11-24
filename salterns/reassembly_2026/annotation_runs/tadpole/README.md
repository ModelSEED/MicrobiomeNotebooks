# Tadpole pilot: can the (iii) MAG-source assemblies be regenerated?

Started 2026-11-11 on poplar. Nothing here is committed, and nothing was written to KBase.

**Results (2026-11-12): see `RESULTS.md`.** Arm C (merge + Tadpole + default-mode SPAdes) passes the go rule on R1_A_D1 (0.93 / 1.02 of (iii)). Arm C on SF2_C_D1 was added and is still running.

**Question.** All 450 original MAGs were binned from the Tringe-lab "(iii)" assemblies (Tadpole + SPAdes, about 2018, unpublished). The 2026 KBase metaSPAdes assemblies hold only about 57 % of (iii)'s sequence in contigs ≥ 2.5 kb and about 25–30 % in contigs ≥ 10 kb (`../../binning_experiments/assembly_vs_jgi/README.md`). That README ranks the likely causes as:

1. unmerged pairs;
2. BayesHammer instead of Tadpole correction;
3. the SPAdes mode;
4. quality trimming.

This pilot reruns two samples with the BBTools-style read preparation described in Zhou et al. 2022 to see how much of the gap it closes.

**Decision rule** (from the assembly_vs_jgi README): if an arm reaches about **≥ 85 %** of the (iii) values (Mbp in contigs ≥ 2.5 kb and ≥ 10 kb), roll it out to all 39 samples. If none does, the remaining gap is in something unrecorded, and asking the Tringe lab for (iii) is the only route.

## Targets

From `binning_experiments/assembly_vs_jgi/assembly_vs_jgi.tsv`. The (iii) values are reconstructed from MAG contig names and lean slightly low. The 2026 values are KBase metaSPAdes 3.15.3 `scaffolds.fasta` ≥ 1 kb.

| Sample | SRA run | Assembly | Mbp ≥ 1 kb | Mbp ≥ 2.5 kb | Mbp ≥ 5 kb | Mbp ≥ 10 kb | Mbp ≥ 20 kb | N50 of ≥ 2.5 kb |
|---|---|---|---|---|---|---|---|---|
| R1_A_D1 | SRR5209943 | (iii) target | – | **268.9** | 189.0 | **131.4** | 87.2 | 9,579 |
| | | 2026 KBase | 354.1 | 162.8 | 78.6 | 32.6 | 12.4 | 4,843 |
| | | 85 % of (iii) | – | 228.6 | – | 111.7 | – | – |
| SF2_C_D1 | SRR5210468 | (iii) target | – | **59.3** | 31.1 | **18.2** | 9.0 | 5,344 |
| | | 2026 KBase | 76.9 | 18.3 | 6.7 | 2.0 | 0.5 | 4,029 |
| | | 85 % of (iii) | – | 50.4 | – | 15.5 | – | – |

SF2_C_D1 is the worst-recovered sample: 2026 / (iii) = 0.31 at ≥ 2.5 kb. R1_A_D1 sits near the median: 0.61, against a median of 0.57 over the 39 samples.

## Design

| Arm | Read preparation | Assembler |
|---|---|---|
| A: merge + Tadpole | BBDuk (adapters, then artifacts and phiX) → `bbmerge-auto.sh rem k=62 extend2=50 ecct` → `tadpole.sh mode=correct k=62` on merged and unmerged reads together → BBDuk `qtrim=r trimq=10 minlen=45` on the unmerged pairs | metaSPAdes 3.15.3 `--meta --only-assembler -k 33,55,77,99,127`, `-1/-2` unmerged, `--merged` merged |
| B: Tadpole only | the same BBDuk output → `tadpole.sh mode=correct k=62` on the pairs → BBDuk `qtrim=r trimq=10 minlen=45` | the same command without `--merged` |
| C: SPAdes mode (R1_A_D1 only) | Arm A's reads | `spades.py` default multicell mode, the same flags without `--meta`; starts when Arm A's assembly ends |

- **A vs B isolates the merging.** The two arms share the BBDuk output, the Tadpole settings, the post-correction quality trim and every SPAdes flag except `--merged`.
  - Strictly, it isolates "BBMerge as run here". `ecct` gives the pairs that fail to merge an extra Tadpole correction inside BBMerge.
  - Arm A's Tadpole table also counts each overlap region once (merged) rather than twice (pairs).
  - Both are part of what Zhou et al. describe, so they are not separated.
- **C vs A isolates the SPAdes mode** (cause 3).
- **What separates A and B from the 2026 KBase run.** They differ in read preparation (BBDuk instead of Trimmomatic, plus Tadpole) **and** in corrector: `--only-assembler` skips BayesHammer, because Tadpole replaces it. This is deliberate. Running both would double-correct the reads and add about 3 h of BayesHammer per run, and it matches JGI's practice of using metaSPAdes `--only-assembler` after BBTools correction. The cost is that a gain in B over KBase cannot be split between causes 2 and 4.

**Decisions and their reasons**

- **SPAdes 3.15.3, not poplar's 4.2.0.** The existing assembly env (`metag-hybrid/env/miniforge/envs/assembly`) has SPAdes 4.2.0. Using it would add a version change to the comparison with KBase, so a new env was built: `/scratch1/afreiburger/salterns_2026/envs/spades3153` (bioconda `spades=3.15.3`, python 3.9, pigz).
- **BBDuk settings are the rqcfilter2 defaults** from the installed BBTools 39.81 help, without the human and microbe removal, clumpify and sketch steps:
  - adapter trimming `ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo ftm=5`;
  - `maq=5 minlen=45 mlf=0.333 maxns=0`;
  - then `ref=artifacts,phix k=31 hdist=1`.
  - The raw reads are 151 bp, so `ftm=5` trims the 151st base, as JGI does. rqcfilter2 does no quality trimming (`qtrim=f`).
- **The quality trim comes after merging, not before.** The task suggested `qtrim=r trimq=10` in the first BBDuk pass. It was moved after merging and correction, as in BBTools' own `assemblyPipeline.sh`, because quality trimming before BBMerge removes overlap and lowers the merge rate. Merged reads are not quality-trimmed. Arm B's pairs get the identical trim, so the A vs B comparison stays clean.
- **Tadpole uses k = 62,** the value in the Tadpole usage example for `mode=correct`, and the same k as the BBMerge extension. Both Tadpole and BBMerge use `prefilter=2`, which keeps singleton k-mers out of the main table to save memory.
- **Every Java step has an explicit `-Xmx`.** bbmerge-auto.sh otherwise takes 84 % of RAM. The limits are 200 GB for BBMerge and Tadpole and 16 GB for BBDuk.
- **Threads: at most 24 per sample and 48 for the track.** Preprocessing runs one step at a time at 24 threads. Each SPAdes run gets 12 threads (`-t 12 -m 400`). Arm A's assembly starts as soon as its reads are ready, and Arm B's preparation (12 threads) runs beside it. Everything runs at `nice -n 10`.

## Where everything is on poplar

Work root: `/scratch1/afreiburger/salterns_2026/tadpole_pilot/`

| Path | What |
|---|---|
| `scripts/` | copies of `scripts/` here: `download.sh`, `run_sample.sh`, `eval_assemblies.py`, `targets.tsv`, `status.sh`, `smoke.sh` |
| `reads/` | `SRR5209943_{1,2}.fastq` (R1_A_D1) and `SRR5210468_{1,2}.fastq` (SF2_C_D1); `ena_filereport.tsv` (ENA portal API: fastq_md5, read_count, base_count); `download.log`; `*.vdb-validate.log`; markers `<run>.ok`, `<run>.plusfixed` |
| `<sample>/qc/` | BBDuk output `clean_{1,2}.fq.gz` (shared by both arms), `adapt_stats.txt`, `filter_stats.txt` |
| `<sample>/armA/` | `merged.fq.gz`, `unmerged_{1,2}.fq.gz`, `ihist.txt`; Tadpole-corrected `ecc_*`; SPAdes input `reads_{1,2}.fq.gz` + `ecc_merged.fq.gz`; assembly in `spades/` |
| `<sample>/armB/` | `ecc_{1,2}.fq.gz`, `reads_{1,2}.fq.gz`, `spades/` |
| `R1_A_D1/armC/spades/` | default-mode SPAdes on the Arm A reads |
| `<sample>/logs/<step>.log` | one log per step; `.done/<step>` markers; `STATUS`; `driver.log`; `driver.pid`; `armA_chain.pid` |
| `smoke/` | smoke test (200,000 R1_A_D1 pairs, k = 33,55) |
| `../envs/spades3153/` | SPAdes 3.15.3 env (`spades3153.create.log`) |

Tools and versions:

- BBTools 39.81 (openjdk 25), from `/scratch1/afreiburger/emilykin/processed/reanalysis_norm/env_bbtools` (read only)
- SPAdes 3.15.3 with Python 3.9.19, from `/scratch1/afreiburger/salterns_2026/envs/spades3153`
- sra-tools 3.4.1, from `/scratch1/afreiburger/metag-hybrid/env/miniforge/envs/sra`

## Exact commands

**Download** (`scripts/download.sh`). ENA https ran at 0.01–1.3 MB/s from poplar, so the reads come from NCBI's AWS open-data mirror (about 50 MB/s):

```
curl -C - -o <run>.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/<run>/<run>
vdb-validate <run>.sra                         # internal MD5/CRC checks
fasterq-dump --split-files -e 6 -t tmp_fqd -O . <run>.sra
```

The pair count and base count must equal ENA's `read_count` / `base_count`. ENA's `fastq_md5` cannot be compared, because fasterq-dump writes different header lines.

**Per sample** (`scripts/run_sample.sh <sample> <run> [armC]`; every command is prefixed with `nice -n 10`):

```
# 0b  separator lines -> bare "+" (BBTools 39.81 needs it; see Gotchas)
awk 'NR%4==3{print "+"; next} {print}'
# 1   BBDuk, rqcfilter2 defaults (no qtrim)
bbduk.sh -Xmx16g t=24 in1= in2= out1=qc/trim_1 out2=qc/trim_2 ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo ftm=5 maq=5 minlen=45 mlf=0.333 maxns=0 ordered
bbduk.sh -Xmx16g t=24 in1=qc/trim_1 in2=qc/trim_2 out1=qc/clean_1 out2=qc/clean_2 ref=artifacts,phix k=31 hdist=1 ordered
# Arm A
bbmerge-auto.sh -Xmx200g t=24 in1=qc/clean_1 in2=qc/clean_2 out=merged outu1=unmerged_1 outu2=unmerged_2 rem k=62 extend2=50 ecct prefilter=2 ihist=ihist.txt ordered
tadpole.sh -Xmx200g t=24 mode=correct k=62 prefilter=2 in=unmerged_1 in2=unmerged_2 extra=merged out=ecc_u_1 out2=ecc_u_2 ordered
tadpole.sh -Xmx200g t=24 mode=correct k=62 prefilter=2 in=merged extra=unmerged_1,unmerged_2 out=ecc_merged ordered
bbduk.sh -Xmx16g t=24 in1=ecc_u_1 in2=ecc_u_2 out1=reads_1 out2=reads_2 qtrim=r trimq=10 minlen=45 ordered
spades.py --meta --only-assembler -k 33,55,77,99,127 -t 12 -m 400 --phred-offset 33 -1 reads_1 -2 reads_2 --merged ecc_merged -o armA/spades      # background
# Arm B (12 threads, beside Arm A's assembly)
tadpole.sh -Xmx200g t=12 mode=correct k=62 prefilter=2 in=qc/clean_1 in2=qc/clean_2 out=ecc_1 out2=ecc_2 ordered
bbduk.sh -Xmx16g t=12 in1=ecc_1 in2=ecc_2 out1=reads_1 out2=reads_2 qtrim=r trimq=10 minlen=45 ordered
spades.py --meta --only-assembler -k 33,55,77,99,127 -t 12 -m 400 --phred-offset 33 -1 reads_1 -2 reads_2 -o armB/spades
# Arm C (R1_A_D1 only; starts when Arm A's assembly ends, in its 12-thread slot)
spades.py --only-assembler -k 33,55,77,99,127 -t 12 -m 400 --phred-offset 33 -1 armA/reads_1 -2 armA/reads_2 --merged armA/ecc_merged -o armC/spades
```

**Evaluation** (`scripts/eval_assemblies.py`, stdlib only). Run it when the assemblies finish, or at any time for partial results:

```
cd /scratch1/afreiburger/salterns_2026/tadpole_pilot
/scratch1/afreiburger/salterns_2026/envs/spades3153/bin/python -I scripts/eval_assemblies.py \
    --root $PWD --targets scripts/targets.tsv --out-prefix $PWD/eval/pilot
```

- **Output files:**
  - `eval/pilot_assemblies.tsv`: for each arm, `scaffolds.fasta` and `contigs.fasta`. It gives n and Mbp in sequences ≥ 1, 2.5, 5, 10 and 20 kb; N50 of sequences ≥ 1 kb and ≥ 2.5 kb; max length; the length-weighted median SPAdes k-mer coverage of sequences ≥ 2.5 kb; the SPAdes insert size; and ratios to (iii) and to 2026 KBase. `go_85pct_of_iii` = yes when Mbp ≥ 2.5 kb and Mbp ≥ 10 kb both reach 85 % of (iii).
  - `eval/pilot_prep.tsv`: BBDuk retention, BBMerge join %, ambiguous % and average insert, Tadpole errors detected and corrected, and step times.
- **Which file to compare.** The primary comparison is `scaffolds.fasta`, which is what KBase saved. `contigs.fasta` is reported as well, because whether (iii) used contigs or scaffolds is unknown.
- **`targets.tsv`** holds all 39 samples, taken from `assembly_vs_jgi.tsv`, so the same script works for a full rollout.

## How to check status

From any shell that can `ssh poplar`, including a fresh login in tmux. The jobs are detached with `setsid nohup`, so they do not depend on any ssh connection staying up.

```
ssh poplar 'bash /scratch1/afreiburger/salterns_2026/tadpole_pilot/scripts/status.sh'
```

`status.sh` prints the following:

- the host load;
- the download markers;
- for each sample: whether the driver PID is alive, the current step, the finished steps and the last driver lines;
- for each arm: whether SPAdes has finished, and its last log line;
- the live SPAdes and Java processes, with CPU and RSS.

A sample is complete when `<sample>/.done/ALL` exists (`STATUS` = `ALL_DONE`). A failed step shows `STATUS` = `FAILED <step>`; read `<sample>/logs/<step>.log`. The drivers can be resumed: rerunning the same `nohup setsid bash scripts/run_sample.sh ...` command skips finished steps.

## Gotchas found while setting up

1. **The ENA download was too slow.** ENA's https mirror stalled at about 8 KB/s and measured 1.3 MB/s at best. NCBI's open-data mirror gives about 50 MB/s, so the reads come from NCBI and are checked with vdb-validate plus exact ENA read and base counts.
2. **BBTools 39.81 rejects fasterq-dump's separator line.** Its multithreaded FASTQ reader (`stream/FASTQ.java`, `quadToReadVec`) asserts that the separator line is exactly `+`. fasterq-dump writes `+SRR5209943.1 1 length=151`.
   - The worker threads die, and BBDuk then hangs at 0 % CPU instead of exiting.
   - The 200,000-pair smoke file passed, apparently because a small input takes a different reader path.
   - Fix: `run_sample.sh` rewrites the separator lines once per run (step `fix_plus_<run>`). First attempt: `R1_A_D1/logs/bbduk_adapt.try1.missing_plus.failed`.
3. **Tadpole's documented one-call form fails.** `tadpole.sh in=u1,merged in2=u2,null` dies with `Can't read file 'null'` in 39.81. It is replaced by two passes with `extra=`; each pass builds its k-mer table from all the reads.
4. **SPAdes 3.15.3 logs "Files with merged reads were ignored.".** This only concerns the read-length estimate (`support.get_primary_max_reads_length(..., ["merged reads"])`). The merged reads are assembled: the smoke log shows "Converting merged reads", an average read length of 176.6 (Arm B: 148.7), and 597,740 processed reads = 2 × (2 × 105,818 unmerged + 87,234 merged).
5. **bbmerge-auto.sh sizes its heap itself.** It takes 84 % of RAM, so every Java step passes `-Xmx`.

## Status (2026-11-11, 11:13 CST on poplar)

| Item | State |
|---|---|
| SPAdes 3.15.3 env | built: `envs/spades3153` |
| R1_A_D1 reads (SRR5209943) | verified. vdb-validate "consistent"; 29,912,586 pairs and 9,033,600,972 bases = ENA |
| SF2_C_D1 reads (SRR5210468) | verified. vdb-validate "consistent"; 30,085,210 pairs and 9,085,733,420 bases = ENA |
| Smoke test (200k pairs, k = 33,55) | **passed end to end**: BBDuk, BBMerge (45.2 % joined, mean insert 245 bp, as in `merge_rate_R1_A_D1.txt`), both Tadpole passes, quality trim, metaSPAdes A and B, Arm C default mode, eval |
| R1_A_D1 driver | **PID 3919366** (`R1_A_D1/driver.pid`), started 10:58:37 with Arm C. Separator fix done (0 non-bare lines left); BBDuk adapter pass done 11:11; artifact/phiX pass running; BBMerge next |
| SF2_C_D1 driver | **PID 3919367** (`SF2_C_D1/driver.pid`), started 10:58:37 (no Arm C). Separator fix done 11:11 (0 non-bare lines left); BBDuk running cleanly since 11:11 |
| Download job | PID 3582415, finished 11:05:26 with fail=0 (`reads/download.done`); `.sra` files deleted after the checks |

**Expected runtime.** These are estimates. Load on poplar is about 450–560, and every step runs at nice 10.

| Step | Estimate | Note |
|---|---|---|
| BBDuk, two passes | about 10–20 min | |
| bbmerge-auto with Tadpole extension and ecct | about 1–2 h | slowest preprocessing step |
| Tadpole correction, two passes (Arm A) | about 0.5–1 h each | |
| Tadpole correction (Arm B) | about 0.5–1 h | |
| metaSPAdes `--only-assembler`, 12 threads | **about 6–12 h per arm per sample** | KBase took about 2.5 h of assembly at 64 threads |

- Arm A and Arm B for both samples: done roughly 8–15 h after the start, i.e. overnight 2026-11-11/12.
- Arm C on R1_A_D1 starts after its Arm A and needs about another 6–12 h. Default multicell mode on a metagenome may be slower or hit the `-m 400` cap; it is optional, so if it fails, the A/B result still stands.

**Results so far.** Full data, R1_A_D1 BBDuk adapter pass:

- input 59,825,172 reads, 9,033,600,972 bases;
- 3.36 % of reads adapter-trimmed;
- 3.62 % of reads removed (maq / minlen / maxns);
- **95.57 % of bases kept**, against 93.61 % for the 2026 Trimmomatic step.

The read-preparation figures below come from the 200k-pair smoke test only. The full-data numbers will be in `eval/pilot_prep.tsv`.

| R1_A_D1, 200k pairs | Value |
|---|---|
| BBDuk adapters | 95.7 % of bases kept (3.1 % of reads adapter-trimmed) |
| BBDuk artifacts and phiX | 4 reads removed |
| BBMerge `rem k=62 extend2=50 ecct` | 45.19 % joined, 54.80 % ambiguous, 0.02 % no solution, mean insert 245.2 bp. This is a lower bound: extension depends on depth |
| Tadpole k = 62, Arm B | errors found in 0.08 % of reads; 300 corrected. Few corrections, because depth on a subsample is low |
| Quality trim after correction | 99.0 % (A) and 99.3 % (B) of bases kept |

No assembly results yet; the smoke assemblies are too small to mean anything.

**When the assemblies finish:**

1. Run the evaluation command above.
2. Copy `eval/` back here: `rsync -a -e "ssh -o BatchMode=yes" poplar:/scratch1/afreiburger/salterns_2026/tadpole_pilot/eval/ eval/`.
3. Apply the decision rule. Arm A ≥ 85 % of (iii) on both samples means roll out to 39 samples with the same driver. A gain in A over B is the merging effect. Arm C vs A is the SPAdes-mode effect.

## Independent review (2026-11-11, 11:15–11:21 CST)

A second agent checked the setup on poplar before the long steps. **Nothing needed fixing**, so no job was killed or relaunched. PIDs 3919366 (R1_A_D1) and 3919367 (SF2_C_D1) are still the drivers.

What was checked:

- **Scripts.** The md5 of every file in `scripts/` is the same locally and on poplar. `run_sample.sh` matches the commands above, and the running BBDuk and BBMerge command lines (from `ps` and the step logs) match it too.
- **Reads.** `download.log` gives exact ENA pair and base counts for both runs, and vdb-validate reports "consistent" with QUALITY md5 ok. The quality strings vary (B–G, `.`, `9`, `=`), so these are original qualities, not SRA Lite.
  - The separator rewrite lost no records. BBDuk read exactly ENA's totals after the fix: R1_A_D1 59,825,172 reads / 9,033,600,972 bases; SF2_C_D1 60,170,420 / 9,085,733,420.
  - No `*.plusfix.tmp` files are left.
- **BBDuk parameters vs `rqcfilter2.sh` help (39.81).** trimk 23, mink 11, trimhdist 1, maq 5, minlength 45, mlf 0.333, maxns 0, qtrim f and filterk 31 all match. rqcfilter2's `trimpolyg=6` is left out, which does not matter for these 2015 HiSeq reads.
- **Full-data QC so far.**
  - R1_A_D1: the adapter pass kept 95.57 % of bases; the artifact/phiX pass removed 914 reads (0.00 %).
  - SF2_C_D1: the adapter pass kept 85.86 % of bases, with 12.93 % of reads adapter-trimmed. That is consistent with the 2026 Trimmomatic step keeping 84.23 % for this sample.
- **Tadpole `extra=`.** In the 39.81 help it means "files for use as kmer data, but not for error-correction". So both Arm A passes build their k-mer table from all the reads, as described.
- **Smoke-test SPAdes logs.**
  - SPAdes 3.15.3 runs on Python 3.9.19 from the env (`#!/usr/bin/env python` resolves through the driver's PATH).
  - Arm A and Arm C convert and count the merged reads: 597,740 processed = 2 × (2 × 105,818 + 87,234).
  - Max read length without merged reads is 150, so k = 127 is valid.
  - None of the BBMerge unmerged reads is longer than 151 bp, so the `extend2` extensions are not carried into the unmerged output.
- **Targets.** `targets.tsv` for R1_A_D1 and SF2_C_D1 equals `binning_experiments/assembly_vs_jgi/assembly_vs_jgi.tsv`. The 85 % go rule in `eval_assemblies.py` uses `orig_Mbp_ge2500_est` and `orig_Mbp_ge10000_est`, as the decision rule says.
- **CPU, measured from `/proc` tick counts, including bgzip children.**

  | Window | R1_A_D1 | SF2_C_D1 | Total |
  |---|---|---|---|
  | 10 s, both samples in BBDuk | 13.8 cores | 23.0 cores | 36.8 |
  | 30 s, R1_A_D1 in BBMerge | 8.1 cores | 11.7 cores | 19.8 |

  Each JVM has 213 G1 GC threads by default on this 336-core host, but GC used about 0 CPU seconds. Memory: 1.55 TB free; BBMerge RSS was 40 GB at 2 min.
- **Processes.** No stray tadpole-pilot processes are left from try 1, the smoke test or the download.

**Points to watch. None justified a relaunch.**

1. **Thread count while both arms run.** The 24-per-sample budget does not count BBTools' bgzip writer threads (`bgzip -@ 24`, about 0.8 core each as measured). While Arm A's SPAdes (12) runs beside Arm B's Tadpole (12 plus about 2 bgzip) on both samples, the track can briefly reach about 50 cores. If `status.sh` shows a sustained total above 48, lower `TASM` or `TPRE` for the rollout.
2. **Memory.** The Java heaps are fixed at 200 GB (`-Xms` = `-Xmx` in the BBTools wrappers), so at its peak each sample can hold about 200 GB for Arm B's Tadpole plus about 100–150 GB for SPAdes.
3. **A weakness in `fix_plus` to fix before the 39-sample rollout.** It replaces the original FASTQ whenever the rewritten file is smaller, without checking awk's exit status. It also checks for leftover separator lines only in R1. Here the exact BBDuk input counts show that nothing was lost. For the rollout, use `awk ... > tmp && mv` and check both mates.
   - Do **not** edit `scripts/run_sample.sh` on poplar while a driver runs. bash reads the script as it executes, so editing the file would corrupt the running driver.
