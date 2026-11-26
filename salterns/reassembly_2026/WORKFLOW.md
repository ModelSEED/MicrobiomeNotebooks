# Salterns metagenomes: processing workflow, start to finish (reproduction guide)

Written 2026-11-18 for a future session, human or Claude, that has to rebuild the environment and re-run or extend this work without the conversation that produced it.

Every version, checksum, path and parameter below was read from the run's own scripts, logs and environment exports. Where a fact is not recorded, the text says so rather than guessing.

**The final workflow is the "Arm C rollout"** (section 4): per-sample re-assembly and re-binning of all 39 salterns metagenomes.

- **Status on 2026-11-18:** 28/39 samples complete on poplar. The last 8, all R2, are assembling on frankensteinsmonster. `DONE` is expected around 2026-11-19/20.
- **What came before** (section 8): the original 2020 MAGs and a first-pass 2026 KBase pipeline. Both are kept only as comparators.
- **Downstream analyses** (section 6) were developed on the first-pass MAGs. They are to be re-applied to the rollout MAGs with the same methods.

## Contents

1. [Data](#1-data)
2. [Hosts and conventions](#2-hosts-and-conventions)
3. [Software environments (exact)](#3-software-environments-exact)
4. [Reference databases (exact)](#4-reference-databases-exact)
5. [The per-sample workflow (Arm C rollout)](#5-the-per-sample-workflow-arm-c-rollout)
6. [Orchestration: driver, scheduling, two hosts](#6-orchestration-driver-scheduling-two-hosts)
7. [Downstream analyses](#7-downstream-analyses)
8. [History: why this recipe (originals, KBase first pass, pilots)](#8-history-why-this-recipe)
9. [Gotchas that cost time](#9-gotchas-that-cost-time)
10. [Where everything is](#10-where-everything-is)
11. [Minimal reproduction recipe](#11-minimal-reproduction-recipe)

---

## 1. Data

**39 shotgun metagenomes** from solar-saltern ponds: Illumina paired-end reads, JGI/IMG.

- **Manifest:** `reassembly_2026/manifest.tsv`, one row per sample (39 rows + header). Columns: `sample, type, site, coverage_group, layer, srr, img, gold, bioproject, biosample, spots, bases, sra_release, odp_url, odp_bytes`.
- **Source:** NCBI SRA run accessions (`srr`), downloaded from the NCBI SRA open-data mirror on AWS (`odp_url`, e.g. `https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5209943/SRR5209943`).
- **Sizes:** 6.1–12.4 Gbp per sample. Samples above 9 Gbp get a larger Java heap in read preparation (section 5.3).
- **Coverage groups** (`coverage_group`): each sample is binned with coverage from every library of its site.

| Group | Libraries | Samples |
|---|---|---|
| R1 | 9 | R1_{A,B,C}_{D1,D2,H2O} |
| R2 | 12 | R2_{A,B,C}_{D1,D2,H2O} + R2_restored_C_black, R2_restored_D_shore, R2_restored_H2O |
| R2A | 9 | R2A_{A,B,C}_{D1,D2,H2O} |
| SF2 | 9 | SF2_{A,B,C}_{D1,D2,H2O} |

D1/D2 are sediment layers, H2O is water, and "restored" marks the restored-pond samples. The manifest's `type`/`layer` columns are authoritative.

**Comparators.**

- **Original MAGs (2020 study):** 450 in total, 310 dereplicated plus 140 non-dereplicated.
  - FASTA: `salterns/MAG_contigs/` (`index.tsv`); provenance in memory note "Salterns MAG inventory" and `reassembly_2026/provenance/`.
  - They came from an unpublished Tringe-lab per-sample Tadpole + SPAdes re-assembly, reconstructed from the SPAdes `NODE_` names in the original contigs.
- **First-pass 2026 MAGs:** 1,324 MAGs from the KBase pipeline (section 8.2).

## 2. Hosts and conventions

| Host | CPU / RAM / OS | Role | Access |
|---|---|---|---|
| **poplar** (shared cluster node) | AMD EPYC 9634, 336 threads, ~3 TB RAM, Ubuntu 24.04.4 | the rollout driver; GTDB-Tk; KOfam and DRAM runs | `ssh -o BatchMode=yes poplar` through an ssh ControlMaster (login-gce, Duo 2FA). A person must log in once in a persistent tmux session (`ssh -MN -o ControlPersist=yes -o ClearAllForwardings=yes poplar`); non-interactive commands then reuse the master. Never log in interactively from a script. |
| **frankensteinsmonster** ("FM", the workstation) | AMD Threadripper 3970X, 64 threads, 251 GB RAM, Ubuntu 24.04.5 | 15 samples of the rollout (7 R2A, 8 R2) by the same driver; local analyses; first-pass GTDB-Tk | local |

**Conventions** (from the project's CLAUDE.md files and the user's rules):

- **Python on FM:** `~/Documents/py_venv/bin/python` (Python 3.12.3, uv-managed). Install with `uv pip install`. Never use system `python`/`pip`.
- **Conda envs:** built with **micromamba 2.10.0** on FM (`~/Documents/dastool_env/bin/micromamba`) and **mamba 2.6.0** on poplar (`/scratch1/afreiburger/metag-hybrid/env/miniforge/bin/mamba`), with channels `conda-forge` and `bioconda` only.
- **Running jobs:**
  - Every heavy job runs at `nice 10`, detached (`nohup setsid …`).
  - On FM a job runs inside a transient **systemd user scope** with `MemoryMax`, `MemorySwapMax=0`, `OOMPolicy=continue` and `ManagedOOMPreference=avoid`, plus `oom_score_adj +300`. A memory guard (`reassembly_2026/local/mem_guard.sh`) freezes only that scope below 12 GB MemAvailable and thaws it above 30 GB.
  - Never `systemctl --user stop` such a scope: that kills its tasks.
- **Process hygiene:**
  - Never signal processes you did not start.
  - Never `pkill -f`/`pgrep -f` with a pattern that also appears in your own command line (it matches itself).
  - Anchor patterns with `^`.
- **Secrets:** the KBase token is in `~/.kbase/token`; never print it.
- **Writes:** KBase writes happen only on explicit request. Commit and push only when asked.

## 3. Software environments (exact)

**Every environment can be recreated package-for-package** from the explicit lists in `reassembly_2026/workflow_envs/`.

- Each list has every package's URL and md5 (`@EXPLICIT` format); `SPECS.md5` holds checksums of the spec files themselves.
- Most were exported from poplar's working environments. `bbtools.txt`, `checkm2.txt` and `sra.txt` describe pre-existing read-only envs from other projects (EmilyKin, metag-hybrid) that the rollout uses.

```bash
MM=~/Documents/dastool_env/bin/micromamba            # or any micromamba/mamba >= 2
export MAMBA_ROOT_PREFIX=<somewhere>/mamba_root
$MM create -y -p <prefix>/<env> -f reassembly_2026/workflow_envs/<env>.txt
$MM list -p <prefix>/<env> --explicit --md5 > <env>.list.txt   # record what you got
```

(That is what `~/Documents/salterns_rollout_local/envs/create_envs.sh` did on FM for every env except DAS Tool and GTDB-Tk.)

**What has been tested.**

- **Built from these specs and in use on FM:** `bbtools`, `spades3153`, `bt_map`, `bt_metabat2`, `bt_maxbin2`, `bt_concoct`, `bt_gunc`, `bt_checkm1` and `checkm2` (`create_envs.sh`; `~/Documents/salterns_rollout_local/envs/*.create.log`).
- **Exported from working poplar envs, not yet used to build:** `bt_dastool`, `sra`, `gtdbtk280` (exported 2026-11-18 with `mamba list --explicit --md5`) and `dram15` (DRAM's own explicit export).
- **All 13** pass `micromamba create --dry-run -f <spec>` (micromamba 2.10.0, 2026-11-18). Every transaction resolves exactly the spec's packages (53/53, 211/211, 136/136, 67/67 for those four).

| Spec file | Core tool(s), exact build | Other pins that matter | Used for |
|---|---|---|---|
| `sra.txt` | sra-tools 3.4.1 (`2_linux_64`) | | download, `vdb-validate`, `fasterq-dump` |
| `bbtools.txt` | **BBTools/bbmap 39.81** (`h9b5c0a0_1`) | openjdk 25.0.2 | BBDuk, BBMerge, Tadpole |
| `spades3153.txt` | **SPAdes 3.15.3** (`h95f258a_0`) | python 3.9.19 | assembly |
| `bt_map.txt` | bowtie2 2.5.5 (`h63e9258_2`), samtools 1.25, seqkit 2.14.0 | htslib 1.25, python 3.12.16 (**also the driver's Python on poplar**), pigz | mapping, driver |
| `bt_metabat2.txt` | **MetaBAT2 2.12.1** (`1`) | | MetaBAT2, `jgi_summarize_bam_contig_depths` |
| `bt_maxbin2.txt` | **MaxBin2 2.2.7** (`h503566f_8`) | FragGeneScan 1.32, HMMER 3.4, IDBA 1.1.3, bowtie2 2.5.5, perl 5.32.1 | MaxBin2 |
| `bt_concoct.txt` | **CONCOCT 1.1.0** (`py312hb1d17a5_9`) | python 3.12.16, scikit-learn 1.9.1, numpy 1.26.4, pandas 3.0.6, biopython 1.88 | CONCOCT and its helper scripts |
| `bt_dastool.txt` | **DAS Tool 1.1.7** (`r44hdfd78af_1`) | diamond 2.2.8, prodigal 2.6.3, blast 2.17.0, r-base 4.4.3 (`hd767689_12`), pullseq 1.0.2, ruby 4.0.9 | DAS Tool |
| `bt_gunc.txt` | **GUNC 1.1.1** (`pyhdfd78af_0`) | diamond 2.1.24, python 3.14.10, prodigal 2.6.3 | chimerism QC |
| `bt_checkm1.txt` | **CheckM 1.1.3** (`checkm-genome py_1`) | python 3.8.20, pplacer 1.1.alpha19, HMMER 3.4, prodigal 2.6.3, numpy 1.24.4 (built as `checkm-genome=1.1.3 python=3.8 setuptools<70`) | informational QC |
| `checkm2.txt` | **CheckM2 1.1.0** (`pyh7e72e81_1`) | python 3.12.13, diamond 2.1.11, prodigal 2.6.3, tensorflow 2.17.0 (CUDA build; CheckM2 sets `CUDA_VISIBLE_DEVICES=-1`, so it runs on CPU), **skani 0.3.2** | primary QC; skani recovery checks on poplar |
| `gtdbtk280.txt` | **GTDB-Tk 2.8.0** (`pyhdfd78af_0`) | python 3.13.18, pplacer 1.1.alpha19, skani 0.4.0, FastANI 1.34, FastTree 2.2.0, HMMER 3.4, prodigal 2.6.3, numpy 2.5.4, pydantic 1.10.22 (built as `gtdbtk=2.8.0`, conda-forge + bioconda) | taxonomy |
| `dram15.txt` | **DRAM 1.5.0** (`pyhdfd78af_0`) | Python 3.10.22, pandas 1.5.3, SQLAlchemy 1.4.54, numpy 1.26.4, scikit-bio 0.5.8, setuptools 80.10.2, MMseqs2 18.8cc5c, HMMER 3.4, tRNAscan-SE 2.0.13, barrnap 1.10.5 | DRAM annotation of species representatives |

**Notes per environment.**

- **DAS Tool on FM** reuses `~/Documents/dastool_env/env`. It has the same das_tool, diamond, prodigal, blast and pullseq builds as `bt_dastool.txt`; only r-base (4.4.3 `h502d0c9_11`), fonts, libtiff and pango differ.
- **CheckM1** needs `checkm/DATA_CONFIG` inside the env to point at the data directory (or `CHECKM_DATA_PATH` set). With the default `~/.checkm` it fails. On FM the original file is `DATA_CONFIG.orig_default`. The EmilyKin env's checkm-genome 1.2.3 cannot import, because setuptools ≥ 82 has no `pkg_resources`.
- **DRAM 1.5.0** can be built two ways:
  - **from `dram15.txt`:** this is the env after its fix (setuptools 80.10.2), so nothing else is needed;
  - **fresh**, with `mamba create … "dram=1.5.0" "python=3.10" "pandas>=1.5,<2" "sqlalchemy<2"`: then run `mamba install "setuptools<81" "pandas>=1.5,<2" "sqlalchemy<2" "dram=1.5.0"`, because DRAM imports `pkg_resources`.
  - Either way, always run DRAM with `PYTHONNOUSERSITE=1` (`annotation_runs/dram/scripts/dram_env.sh`).
- **Driver scripts** (`driver.py`, `rolib.py`, `btlib.py`) are stdlib-only Python and run under any Python ≥ 3.10.
  - poplar: `bt_map`'s python 3.12.16. FM: `~/Documents/py_venv`.
  - md5 on all hosts: `driver.py` `507e1670ff8002c55fc427b1d64d3420`, `rolib.py` `89425e472d3bd8355c88146d27c84ff1`, `btlib.py` `15d12b482943d7fea594c2ce1653b158`.
- **Other tools used downstream** (section 7):
  - skani 0.4.0 in `~/Documents/gtdbtk_venv` (first-pass species clustering);
  - HMMER 3.4 and Prodigal 2.6.3 (`h577a1d6_11`) and Python 3.11 from the EmilyKin annotation env (KOfam track);
  - Infernal 1.1.5 `cmscan` and Bakta 1.12.0 (EmilyKin Bakta env, Bakta DB v6.0 full 2025-02-24; Sec-rescue SECIS check and the Bakta test).

## 4. Reference databases (exact)

| Database | Version | Source | Integrity | Location (poplar unless noted) |
|---|---|---|---|---|
| **GTDB R232** for GTDB-Tk | release 232 (`232.0`) | `https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/split_package/gtdbtk_r232_data.tar.gz.part_*` (12 parts) | every part against GTDB's `MD5SUM.txt`; whole package md5 **`25a59e0352b1fd150c589f56559767d4`**; `gtdbtk check_install` OK | `/scratch1/afreiburger/salterns_2026/refdata/gtdbtk_r232/extract.tmp/release232` (`GTDBTK_DATA_PATH`; left where `check_install` verified it). Built by `annotation_runs/rollout/scripts/setup_gtdbtk.sh`. FM has its own R232 copy for the local GTDB-Tk runs (`~/gtdbtk_kbase/env.sh`). |
| GTDB R232 mockup | 26 MB test DB | same release, `…/mockup_db/r232_mockup.tar.gz` | md5 from `MD5SUM.txt` | `refdata/gtdbtk_r232_mockup` (smoke test only; it cannot place archaea) |
| **CheckM2 DIAMOND DB** | `uniref100.KO.1.dmnd` | not recorded in this project (the copy predates it, in the EmilyKin `refdata/`). `checkm2 database --download` is CheckM2's own route; check the md5 | md5 **`4db58455c51ff15804b6b8e19ff99f58`** | `/scratch1/afreiburger/emilykin/processed/refdata/checkm2/CheckM2_database/`; FM `~/Documents/salterns_rollout_local/refdata/checkm2/` |
| **GUNC proGenomes 2.1** | `gunc_db_progenomes2.1.dmnd` (the `gunc download_db` default) | `https://black.embl.de/~fullam/gunc/gunc_db_progenomes2.1.dmnd.gz` (GUNC 1.1.1's second URL; zenodo was 1–3 MB/s), 8 byte ranges (`annotation_runs/binning_test/scripts/fetch_gunc_db.sh`) | gz **`bc93a855e0760aad5c4e5f2d0e26da46`**, dmnd **`447c9330056b02f29f30fe81fe4af4eb`** (GUNC 1.1.1's own md5s) | `/scratch1/afreiburger/salterns_2026/binning_test/refdata/gunc/`; FM `…/salterns_rollout_local/refdata/gunc/` |
| **CheckM1 data** | `checkm_data_2015_01_16` | download URL not recorded in this project; check the md5 | md5 **`631012fa598c43fdeb88c619ad282c4d`**; 32 files md5-identical across copies | `/scratch1/afreiburger/salterns_2026/binning_test/refdata/checkm1/`; FM `…/salterns_rollout_local/refdata/checkm1/` (33 files incl. `.ok`) |
| **BBTools resources** | `adapters`, `artifacts`, `phix` refs bundled with bbmap 39.81 | inside the bbtools env | env spec | (named `ref=adapters`, `ref=artifacts,phix`) |
| **KOfam** | **2026-10-26 release**: `ko_list.gz` (Last-Modified 2026-10-26 05:26 GMT, 917,028 B, md5 **`697af5b3b22a919306a1e60f99cec1a6`**); `profiles.tar.gz` (Last-Modified 2026-10-27 01:41 GMT, 1,567,491,980 B, md5 **`bc68ad2a52ed48812b5496d15c215a7d`**; 27,981 profiles) | `https://www.genome.jp/ftp/db/kofam/` | md5 + `gzip -t` | `/scratch1/afreiburger/salterns_2026/kofam_db/`; local copy of the ko_list `annotation_choice/sources/kegg/kofam_ko_list_2026-10-26.gz`. **Use this pair for every KO call.** |
| **DRAM 1.5 DBs** | KOfam as above; **Pfam 39.0** (`Pfam-A.full.gz` md5 `1ca35677412407da16a340bdd293e773`, `Pfam-A.hmm.dat.gz` md5 `13c095351722642b3c51b56d7ecf9a0a`); **dbCAN V11** (`dbCAN-HMMdb-V11.txt`, `CAZyDB.08062022.fam-activities.txt`, `CAZyDB.08062022.fam.subfam.ec.txt`); **MEROPS** `pepunit.lib` (2023-02-22); **VOGDB 237** (built, not searched); RefSeq viral `viral.1.protein.faa.gz` (2026-09-04, unused); no UniRef90, no KEGG Genes | EBI (Pfam, MEROPS), `https://pro.unl.edu/dbCAN2/download/Databases/V11/` (DRAM's `bcb.unl.edu` URLs now return HTML), `fileshare.lisc.univie.ac.at/vog/latest/`, NCBI | `dram_db/raw/MANIFEST.tsv` (URL, Last-Modified, bytes, md5) | `/scratch1/afreiburger/salterns_2026/dram_db/`; setup details in `annotation_runs/dram/README.md` |
| Bakta DB (test only) | v6.0 full, 2025-02-24 | | `version.json` | `/scratch1/afreiburger/emilykin/processed/refdata/bakta/db` |

## 5. The per-sample workflow (Arm C rollout)

**Approval and scope.** Andrew approved the rollout on 2026-11-13 ("All 39 + QC filter") after the binning test's NO_GO (section 8.4). All 39 samples go through Arm C, and only MAGs passing the KEEP rule are kept.

**Overview** (one sample `<s>` with site coverage group `G`):

```
SRA run ──download+validate──▶ raw FASTQ ──BBDuk (adapters; artifacts+phiX)──▶ clean_* (assembly) + map_* (mapping, also q-trimmed)
clean_* ──bbmerge-auto rem extend2 ecct──▶ merged + unmerged ──Tadpole correct ×2──▶ ecc_merged + ecc_u_* ──q-trim──▶ reads_1/2
reads_1/2 + ecc_merged ──SPAdes 3.15.3 DEFAULT mode (not --meta), --only-assembler, k 33..127──▶ scaffolds ──≥1 kb──▶ contigs.fa
map_* of every library in G ──bowtie2──▶ BAMs ──▶ depth tables ──▶ MetaBAT2 · MaxBin2 · CONCOCT ──DAS Tool 1.1.7──▶ bins
bins ──CheckM2 1.1.0 · GUNC 1.1.1 · CheckM1 1.1.3──▶ KEEP = CheckM2 ≥50 % complete AND <10 % contaminated AND GUNC pass
all kept MAGs ──GTDB-Tk 2.8.0 classify_wf, GTDB R232──▶ taxonomy ──▶ report, DONE
```

### 5.1–5.9 Exact commands

These are the task bodies the driver generates (`annotation_runs/rollout/scripts/driver.py`), as recorded in `annotation_runs/rollout/README.md` "Exact commands".

- **Thread counts shown are poplar's.** On FM: SPAdes `-t 16 -m 120` (`-m 160` after an out-of-memory failure) and bowtie2 `--threads 13`. Every other parameter is identical, as the task-script diffs `~/Documents/salterns_rollout_local/task_script_diff_{R2A_B_D1,R2_C_D2}*.txt` show.
- **Every task script** runs with `set -euo pipefail`, `OMP/MKL/OPENBLAS_NUM_THREADS` = the task's threads, a short local `TMPDIR=/tmp/ro_<hash>` and a scratch temp `SCRTMP=$ROOT/tmp/<task>.<attempt>`.

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

**KEEP rule details** (`rolib.py keep`, binning-test conventions).

- Comparisons are exact decimals on the printed CheckM2 values.
- GUNC pass means `pass.GUNC == True` in `GUNC.progenomes_2.1.maxCSS_level.tsv`, joined on GUNC's genome name: **the file name cut at the first `.fa`** (see Gotchas).
- A bin with no GUNC row, or a non-`True` value, does not pass.
- Drop reasons recorded per bin: `completeness<50`, `contamination>=10`, `gunc_fail`, `gunc_na`, `gunc_missing`, `checkm2_missing`.
- **HQ** = kept with CheckM2 ≥ 90 % complete and < 5 % contaminated. **MQ+** = ≥ 50 / < 10.
- CheckM1 is informational only and never decides KEEP.

**Retries.**

| Kind | Attempts | On failure |
|---|---|---|
| download | 3 | retry after ≥ 5 min |
| SPAdes | 4 | attempts 1–2 in default mode (the 2nd with more memory after an OOM); 3–4 use `--meta` (Arm A recipe) into `asm/spades_meta`, marked `meta_fallback` in `asm_meta.json` |
| merge / Tadpole | 2 | the retry gets `-Xmx250g` after a Java OutOfMemoryError (FM: 150g, so the sample goes back to poplar) |
| everything else | 2 | retry after ≥ 5 min |

- OOM detection looks for these strings in the logs: `Cannot allocate memory`, `unable to allocate`, `std::bad_alloc`, `OutOfMemoryError`, `OS return value: -9`, `Killed`, `Out of memory`.
- SPAdes 3.15.3 with too small an `-m` prints `mimalloc: error: unable to allocate OS memory`.
- As of 2026-11-18 no sample has needed `meta_fallback` and none has failed.

**Comparator QC.** CheckM2 + GUNC with the same commands on:

- the first-pass 2026 DAS Tool bins of the same unit (task `checkm2:old`/`gunc:old`, background priority; inputs uploaded to `rollout/inputs/old/<unit>/`, md5-checked);
- the 450 original MAGs (`/scratch1/afreiburger/salterns_2026/vs2020/run.sh`, `rerun.sh`, `rest.sh`; inputs `vs2020/originals_all/`).

Everything is then compared under one rule.

**Reused, not re-run, in this rollout.** For a fresh reproduction, run these steps for all 39 samples with the commands above.

- **BBDuk reads of the 18 R1 + SF2 libraries** come from the binning test (`/scratch1/afreiburger/salterns_2026/binning_test/reads/<lib>`, symlinked; `reuse` in `config.json`).
  - They were made with the same three BBDuk commands (without `-Xms8g eoom`, which do not change the output).
  - The two pilot libraries resolve to `tadpole_pilot/<s>/qc/`, made with the same rqcfilter2-default settings.
- **R1_A_D1 and SF2_C_D1** reuse their Arm C assemblies from the Tadpole pilot and their binning and QC from the binning test.
  - The driver checks that the ≥ 1 kb scaffolds are byte-identical (`cmp`) to the binning test's `asm/contigs.fa`, and that there is one QC row per bin.
- So only the 21 R2 + R2A libraries were downloaded (the `dl:` comment above).

## 6. Orchestration: driver, scheduling, two hosts

All code is in `reassembly_2026/annotation_runs/rollout/scripts/`, with live copies in each run root.

| File | Role |
|---|---|
| `driver.py` | the scheduler: builds the task DAG from `config.json`, writes each task as a bash script under `scripts/tasks/<task>.<attempt>.sh` (never edited after), launches detached tasks, keeps markers and a ledger. Commands: `run`, `plan` (task states from markers, launches nothing), `run-one --task`, `report`, and `run --only-kinds dl` |
| `rolib.py` | `asm-stats`, `spades-meta`, `keep`, `verify-import`, `report` (results/*.tsv + REPORT.md) |
| `btlib.py` | `merge-depth` (KBase combineDepthTablesUtil logic), `abund` (MaxBin2 abundance files), `filter-fasta`, `bbduk-counts`, `das-inputs`, `fasta-stats` |
| `make_config.py` | writes poplar's `config.json` from `manifest.tsv`. Groups by `coverage_group`; sample order puts R2_restored_D_shore and SF2_A_H2O first, then R1/SF2, then R2/R2A. Holds `MIGRATED` / `MIGRATED_R2`, the samples run on FM as staged imports |
| `make_config_local.py` | FM's config (`local_config.json`) for the migrated samples, with the FM resource policy |
| `setup_gtdbtk.sh` | GTDB-Tk 2.8.0 env + R232 download/verification/extraction |
| `status.sh` | status of a run root (scheduler line, holds, samples, running tasks, staged imports) |
| `restart_driver.sh` | poplar restart without disturbing tasks: STOP → acknowledged → SIGTERM only while the driver sleeps → relaunch → check adoptions → remove STOP. FM twin: `~/Documents/salterns_rollout_local/restart_local.sh` |
| `dryrun_migration.sh`, `gen_task_scripts.py` | checks for any config change: `plan` old vs new config on copies of a root's markers; write a sample's task scripts from a config without running them, to diff hosts |
| `upload_inputs.sh`, `smoke_prep.sh`, `smoke_config.json`, `local_smoke_config.json`, `smoke_import_config.json` | inputs of the 2026-side comparison; smoke tests |

**Launch on poplar** (root `/scratch1/afreiburger/salterns_2026/rollout`):

```bash
cd /scratch1/afreiburger/salterns_2026/rollout
# exactly one driver per root; this must print 0 first:
pgrep -f '^/scratch1/afreiburger/salterns_2026/envs/bt_map/bin/python -I scripts/driver.py run --config scripts/config.json' | wc -l
nohup setsid nice -n 10 /scratch1/afreiburger/salterns_2026/envs/bt_map/bin/python -I scripts/driver.py run \
    --config scripts/config.json >> logs/driver.log 2>&1 < /dev/null &
ssh -o BatchMode=yes poplar 'bash /scratch1/afreiburger/salterns_2026/rollout/scripts/status.sh'    # status
```

**Launch on FM** (root `~/Documents/salterns_rollout_local/rollout`): `bash ~/Documents/salterns_rollout_local/launch_local.sh rollout`.

- It opens a new scope `salterns-r2a-local-<utc>.scope` (MemoryMax 220G) plus a memory guard.
- Status: `bash ~/Documents/salterns_rollout_local/status_local.sh`.
- Restart: `bash ~/Documents/salterns_rollout_local/restart_local.sh --rm-stop`.

**Run-root layout** (both hosts):

| Path | Content |
|---|---|
| `reads/<lib>/` | `clean_*`, `map_*`, `read_counts.json`, `download.json`, logs |
| `<s>/prep/`, `<s>/asm/` | `final` → `spades`, `contigs.fa`, `asm_stats.json`, `asm_meta.json` |
| `<s>/bam`, `<s>/depth`, `<s>/{metabat2,maxbin2,concoct}`, `<s>/dastool/out/<s>_DASTool_bins/` | mapping, depth and binning outputs |
| `qc/{new,old}/<s>/{checkm2,checkm1,gunc}/` | QC outputs |
| `<s>/keep/{bins.tsv,kept.txt,summary.json}`, `mags_kept/<s>__<bin>.fa` | KEEP results and the kept MAGs |
| `gtdbtk/`, `results/` | `samples_summary.tsv`, `bins_all.tsv`, `mags_kept.tsv`, `bins_dropped.tsv`, `REPORT.md` |
| `STATE.json`, `ledger.jsonl` | live state (rewritten every 20 s); append-only ledger |
| `.done/ .failed/ .run/ .rc/ .state/` | markers |
| `logs/driver.log`, `logs/tasks/` | driver log; one log per task attempt |
| end markers | **`DONE`** or **`FAILED` (+`FAILED.json`)** |

**Scheduling policy (poplar, shared).**

- **Thread cap:** `cap = min(min(256, max(64, free_cores + ours − 16)), max(64, free_cores − 16))`, where `free_cores = 336 − (load1 − our threads)`. At other users' load ≥ 256 the cap sits at the 64-thread floor.
- **Assembly lanes:** `min(6, cap // 40)`, protected ahead of downstream work.
- **Memory:** a task starts only if `MemAvailable − Σ(reserve − RSS of running tasks) − its reserve ≥ 400 GB`. Reserves: merge/Tadpole = their heap; SPAdes 80; CheckM2/CheckM1 40; GUNC 30; CONCOCT 30; GTDB-Tk 170.
- **Limits:** ≤ 3 read preparations at a time, ≤ 6 QC tasks, `prep_buffer 2`.
- **Disk gates:** ≥ 1 TB free on /scratch1; ≥ 10 GB `/tmp` for CheckM2/CheckM1/GTDB-Tk.
- **Intermediates** are deleted as soon as they are no longer needed (raw FASTQ, merged reads, SPAdes `K*/`, BAMs…).
- **Priority:** downstream first (keep, QC, binning) before new samples.

**FM policy** (`make_config_local.py`):

- **Threads:** cap fixed at 48. SPAdes runs 16 threads, 2 lanes, plus one 16-thread lane for binning/QC.
- **Read prep:** one at a time, and only when `MemAvailable − outstanding − heap ≥ 40 GB` (so ≥ 190 GB free for a 150g heap). `asm_after_all_prep`: all preps finish before the first assembly.
- **Memory booking:** SPAdes `-m 120`, reserve 60; CheckM1 reserve 80. The sum of max(reserve, RSS) stays ≤ 210 GB.
- **Not on FM:** no GTDB-Tk and no 2026-side QC.
- **Do not loosen the prep gate:** 150g/200g-heap Tadpole runs reached 144–157 GB RSS on poplar.

**Two hosts.** A sample listed in `MIGRATED` runs end to end on FM.

- FM's `push:<s>` task sends `<s>/`, `qc/new/<s>/` and `mags_kept/<s>__*.fa` to poplar's `imports/<s>.part/` with an md5 manifest. It verifies the manifest on poplar, then renames the folder to `imports/<s>/` and touches `READY`.
- poplar's `import:<s>` waits (`import_wait`) for `READY`, re-verifies md5, copies into place, re-applies KEEP and requires byte-identical kept FASTAs. `keep:<s>`, GTDB-Tk, the report and `DONE` then follow on poplar.
- If FM writes `imports/<s>/FAILED` on poplar, the run ends `FAILED` for that sample.

Migration history:

- **2026-11-16:** 7 R2A samples ≤ 9 Gbp (all imported 2026-11-17/18).
- **2026-11-17:** 8 R2 samples (R2_B_D1, R2_C_D1, R2_C_D2, R2_restored_C_black, R2_A_H2O, R2_B_H2O, R2_C_H2O, R2_restored_H2O).
- Samples above 9 Gbp (200g heap) stay on poplar.

Procedure for any move (both times):

1. Edit `MIGRATED` in both generators and regenerate.
2. `dryrun_migration.sh` on both hosts: removed/added tasks only, 0 state changes.
3. `gen_task_scripts.py` diff of one sample: only thread/memory parameters differ.
4. Copy the site's read libraries with poplar-side md5 lists (`~/Documents/salterns_rollout_local/logs/{copy_inputs*.sh,verify_inputs.sh}`).
5. Restart both drivers with their restart scripts.

Full records: `annotation_runs/rollout/README.md`, sections "R2A migrated…" and "R2 migrated…".

## 7. Downstream analyses

**Status of the downstream steps.**

- **Applied so far** to the first-pass 1,324 MAGs (section 8.2), to the originals, and to the rollout's kept MAGs per completed sample (QC and recovery numbers).
- **After `DONE`**, the same steps are to be run on the complete rollout MAG set:
  1. collect kept MAGs + GTDB-Tk R232;
  2. dereplicate;
  3. compare with the 2026 and original MAGs, including the 30 originals the first pass lost;
  4. KOfam/guild calls;
  5. refresh the manuscript.

### 7.1 Recovery of original MAGs and the manuscript numbers

- **skani recovery** (poplar, `vs2020/skani_completed.sh`, re-runnable): skani **0.3.2** `triangle -l genomes_completed.txt -t 16 -E` over the 450 originals plus the kept MAGs of every completed sample.
  - An original counts as recovered when a kept MAG **of the same sample** has ANI ≥ 95 % and the larger aligned fraction ≥ 50 %.
- **`reassembly_2026/manuscript_numbers.py`** (run on FM with `~/Documents/py_venv/bin/python`; it reads poplar over ssh):
  - Per complete sample: NEW (kept, HQ, median N50), 2026 first pass (same rule), originals (same rule; CheckM2/GUNC from `vs2020/`), originals recovered, and recovery of the 30 "lost" originals.
  - Writes `manuscript_numbers_<date>.json` and `table_s7_partial_<date>.tsv`.
  - **Run `vs2020/skani_completed.sh` on poplar first.**

### 7.2 Species clustering / dereplication

`reassembly_2026/compare_vs_original.py` → `comparison/` (first pass vs the 310 originals):

1. skani **0.4.0** (`~/Documents/gtdbtk_venv/bin/skani`) `triangle -t 48 -E --min-af 30` over all genomes.
2. Edges with ANI ≥ 95 % and max(aligned fraction) ≥ 50 % (GTDB's species thresholds).
3. Greedy centroid clustering, best genome first, ranked by `completeness − 5 × contamination`.
4. Best MAG per cluster = the dereplicated set.

The first pass gave 558 species clusters. The same rule applies to the rollout set; record which skani version is used, since 0.3.2 and 0.4.0 both exist in these envs.

### 7.3 Functional annotation for the guilds: KOfam (DRAM 1.5 KO rule)

Track `annotation_runs/kofam/` (README, RESULTS, `scripts/`); work dir `poplar:/scratch1/afreiburger/salterns_2026/kofam`.

1. **Genes:** Prodigal 2.6.3 (`h577a1d6_11`) `prodigal -p meta -g 11 -q` on contigs ≥ 1,000 bp. Protein IDs `<MAG>|<contig>_<n>` (`scripts/call_genes.py`, `run_prodigal.sh`).
2. **HMM DB:** all 27,981 KOfam 2026-10-26 profiles concatenated in K order (`scripts/build_db.py`, which also writes `hmm_info.tsv` with threshold and `score_type`).
3. **Search:** whole-MAG protein chunks (`make_chunks.py`, chunks of whole species clusters, one representative per identical protein), then `hmmsearch --cpu 1 --noali -T 14.7 --domT 14.7 -Z <total proteins> --tblout --domtblout` with **HMMER 3.4**.
   - 14.7 is the lowest numeric KOfam threshold. `-Z` is fixed so E-values do not depend on chunking.
4. **Call rule = DRAM 1.5** (`postprocess.py`):
   - significant when score > threshold, using the full or domain score per `score_type`;
   - KOs with threshold `-` are never called;
   - best KO per gene = the significant KO with the lowest full E-value, ties broken by higher score, then KO id;
   - all passing KOs are also kept.
5. **Selenocysteine read-through rescue** for GrdH K21579, GrdB K10672, GrdF K21584 and GrdA K10670 (`scripts/sec_rescue.py build|score`, `sec_full_db.sh`).
   - Prodigal stops at the Sec TGA, so the pieces are joined through the TGA (translated as U) and re-scored against all profiles.
   - Accept when: the best KO is the Sec KO and passes its threshold; an aligned domain spans the U ±10 aa; and the joined score beats every piece.
   - A SECIS element is recorded as evidence, not required: Infernal 1.1.5 `cmscan --cut_tc -g --nohmmonly --rfam -Z 5` with Bakta's ncRNA-region CMs on 603 nt around the TGA.
6. **Guild rules per MAG** (`build_guilds.py`):

| Guild | Rule |
|---|---|
| methanogen | *mcrA* K00399 |
| MePn consumer | *phnJ* K06163 |
| betaine reducer | *grdH* K21579 **or** *grdI* K21578 |
| methylamine user | *mttB* K14083 / *mtbB* K16178 / *mtmB* K16176 **only in an *mcrA* genome**; K14083 is over-called in non-methanogens, reported as `mttB_family_without_mcrA` |
| TMA producer | *cutC*+*cutD* (K20038+K20037) **or** *cntA*+*cntB* (K22443+K22444) |

### 7.4 DRAM 1.5 (species representatives)

`annotation_runs/dram/` (README with exact setup); work dir `poplar:/scratch1/afreiburger/salterns_2026/dram`. Run on the 558 best-per-species first-pass MAGs:

```bash
PYTHONNOUSERSITE=1 nice -n 10 DRAM.py annotate -i mags/<MAG>.fna -o dram_out/<MAG> --min_contig_size 1000 --threads 6 \
  --gtdb_taxonomy inputs/gtdb_taxonomy_558.tsv --checkm_quality inputs/checkm_quality_558.tsv
DRAM.py merge_annotations -i "$W/dram_out/*" -o $W/merged
nice -n 10 DRAM.py distill -i $W/merged/annotations.tsv -o $W/distill --trna_path $W/merged/trnas.tsv \
  --rrna_path $W/merged/rrnas.tsv --custom_distillate inputs/custom_distillate_salterns.dram1.tsv
# then the same distill with --distillate_gene_names into $W/distill_gene_names (non-fatal); see scripts/driver.sh
```

- Settings: DRAM defaults otherwise, with no `--use_vogdb` and no UniRef/KEGG Genes.
- Database setup bypasses `DRAM-setup.py prepare_databases`: it calls `prepare_databases()` directly. That avoids the dead dbCAN URLs and VOG 237's `hmm/` subfolder (see the README).
- DRAM2 (planned on h100) was not reachable and was not used.

### 7.5 Bakta test (not part of the workflow)

Bakta 1.12.0 (DB v6.0 full) was tested on 9 MAGs to decide whether Bakta CDS should replace Prodigal for the KO layer. The decision was to keep Prodigal + the Sec rescue, and to use Bakta only for deposition (`annotation_runs/bakta_test/README.md`, `annotation_choice/README.md`).

### 7.6 Manuscript refresh

The manuscript review document is regenerated with:

```
manuscript_numbers.py → clauded_insights/make_r232_review_docx.py → docx validate
```

- The validation step is `python docx/scripts/office/validate.py` from the docx skill.
- The master text is Susannah's Google Doc (view-only); the local mirror and tracked-changes rules are in memory note "Salterns Google Doc mirror".

## 8. History: why this recipe

### 8.1 Original 2020 MAGs

- **Set:** 310 dereplicated + 140 non-dereplicated, from per-sample DAS Tool binning with single-sample coverage.
- **Assembly:** NODE names show a Tringe-lab per-sample Tadpole + SPAdes re-assembly, not the JGI SPAdes assemblies.
- **Re-scored here:** with CheckM2 1.1.0 + GUNC 1.1.1. GUNC pass: 264/310 dereplicated, 115/140 non-dereplicated.

### 8.2 First pass, 2026 (KBase narrative 275645; `reassembly_2026/pipeline.py`, README)

| Stage | App / tool |
|---|---|
| trim | Trimmomatic 0.39 (`kb_trimmomatic`): ILLUMINACLIP TruSeq3-PE-2 2:30:10, LEADING:3, TRAILING:3, SLIDINGWINDOW:4:15, MINLEN:50 |
| assemble | metaSPAdes 3.15.3 (`kb_SPAdes run_metaSPAdes`, k 33–127; MEGAHIT 1.2.9 fallback) |
| co-assemble | MEGAHIT 1.2.9 site co-assemblies of sediment cores |
| binning | MetaBAT2 / MaxBin2 / CONCOCT KBase apps with site multi-sample coverage |
| refinement | DAS Tool **1.1.7 locally** (KBase's 1.1.2 has a scoring bug) |
| QC | KBase CheckM 1.1.2 `lineage_wf` reduced tree; kept ≥ 50 / < 10 |
| taxonomy | GTDB-Tk 2.8.0 / R232 locally |

- **Result:** 1,324 MAGs in 558 species clusters. They were fragmented: median MAG N50 about 5 kb versus about 16 kb for the originals.
- **Taxonomy writeback** to KBase: `reassembly_2026/local/writeback_gtdbtk_r232.py`.

### 8.3 Tadpole pilot (`annotation_runs/tadpole/`, RESULTS.md)

**Question:** could read preparation explain the fragmentation? Tested on R1_A_D1 and SF2_C_D1.

| Arm | Read prep | SPAdes mode |
|---|---|---|
| A | bbmerge `rem extend2 ecct` + Tadpole | `--meta` |
| B | Tadpole only | `--meta` |
| C | as A | **default (non-meta)** |

**Arm C regenerated the original assembly's size and contiguity** on R1_A_D1: 0.93× / 1.02× of the original Mbp in contigs ≥ 2.5 kb / ≥ 10 kb.

Effects, each measured on top of the one before:

- Tadpole correction alone: about 0.
- Merging with Tadpole extension: 2.5×.
- Default mode: a further 1.6×.

### 8.4 Binning test (`annotation_runs/binning_test/`; pre-registered rule `scripts/decide.py`)

**Design:** Arm C assemblies binned exactly as the first pass, scored with CheckM2 + GUNC. Criteria a–e: MQ+ ≥ 1.2×; HQ ≥ 1.2× or +5; median N50 ≥ 1.5×; contamination/GUNC no worse; original species recovered ≥.

| Sample | MQ+ | HQ | Result |
|---|---|---|---|
| R1_A_D1 | 62 vs 47 | 27 vs 11 | **passed** |
| SF2_C_D1 | 7 vs 4 | | failed only criterion d: contamination 4.67 vs 1.70 |

- **Outcome:** NO_GO by the pre-registered rule.
- **Decision:** Andrew chose **"All 39 + QC filter"**, which became the KEEP rule above.

## 9. Gotchas that cost time

1. **GUNC genome names** are the file name cut at the **first** `.fa`: `x.contigs.fa_assembly.fa` → `x.contigs`. Join GUNC to CheckM2 with the same rule (`rolib.py gunc_key`). Copies made as `.fna` and linked as `.fa` need the same treatment.
2. **AF_UNIX socket paths are limited to 107 characters:** CheckM2's multiprocessing failed with `EOFError` under a long NFS TMPDIR. Use a short local `TMPDIR` (`/tmp/ro_<hash>`).
3. **CheckM1 `DATA_CONFIG`** must point at the data; checkm-genome 1.2.3 needs `setuptools<81` (no `pkg_resources` in newer versions).
4. **CheckM2/GUNC "diamond not found"** when called by absolute path: put the env's `bin/` first on `PATH`.
5. **Java:** use `-Xms8g` (don't reserve the heap up front) and `eoom` (exit on OOM instead of hanging). Heap size can change BBTools' prefilter sizing, so keep heaps identical across hosts (150g ≤ 9 Gbp, 200g above).
6. **The SPAdes OOM message** is `mimalloc: error: unable to allocate OS memory`. `--meta` is a fallback only, so record the mode.
7. **DAS Tool:** KBase's 1.1.2 has a silent scoring bug; use 1.1.7.
8. **Two skani versions** (0.3.2 in the CheckM2 env, 0.4.0 in GTDB-Tk/gtdbtk_venv): state which one produced a number.
9. **KOfam `-` thresholds** (*mpnS*, *phnY\**, *phnZ1*…) are never callable under DRAM 1.5 rules. *phnL*'s threshold changed between releases, so pin the release.
10. **Prodigal builds differ** slightly (2 extra ORFs on one control MAG between bioconda and venv builds). Pin the build.
11. **DRAM 1.5 setup:** dead `bcb.unl.edu` dbCAN URLs return HTML with HTTP 200; VOG 237 moved its HMMs into `hmm/`; `~/.local` site-packages on poplar break imports without `PYTHONNOUSERSITE=1`.
12. **poplar `/tmp`** is small (≈ 22–26 GB) and shared: gate on ≥ 10 GB free.
13. **Restarts:**
    - Stop a driver only while it sleeps (`/proc/<pid>/wchan` = `hrtimer_nanosleep`), after `STOP` is acknowledged.
    - Tasks are detached and are adopted by the next driver.
    - Never run two drivers on one root.
    - Reboots kill tasks; on relaunch the driver counts them as failed attempts and retries.
14. **The poplar ssh master dies with its tmux session.** Keep `ssh -MN` running in a persistent tmux session. On `Permission denied`, ask the person to log in again; never automate Duo.

## 10. Where everything is

| What | Where |
|---|---|
| This guide, env specs | `salterns/reassembly_2026/WORKFLOW.md`, `salterns/reassembly_2026/workflow_envs/` |
| Rollout code + full run record | `salterns/reassembly_2026/annotation_runs/rollout/{README.md,scripts/}` (poplar copy `rollout/README_rollout.md`) |
| Rollout run roots | poplar `/scratch1/afreiburger/salterns_2026/rollout/`; FM `~/Documents/salterns_rollout_local/rollout/` (+ its `README.md`) |
| poplar tool envs / refdata | `/scratch1/afreiburger/salterns_2026/envs/`, `/scratch1/afreiburger/salterns_2026/refdata/`, `…/binning_test/refdata/` |
| FM tool envs / refdata | `~/Documents/salterns_rollout_local/envs/` (specs, create logs, `.list.txt`), `…/refdata/` |
| Pilots / tests | `annotation_runs/{tadpole,binning_test,bakta_test}/` |
| Annotation | `annotation_runs/{kofam,dram}/`; decisions in `annotation_choice/README.md` |
| Comparisons | `reassembly_2026/comparison/`, `compare_vs_original.py`, `manuscript_numbers.py`; poplar `/scratch1/afreiburger/salterns_2026/vs2020/` |
| First-pass pipeline | `reassembly_2026/{pipeline.py,README.md,local/,provenance/}`; KBase narrative 275645 |
| Inputs | `reassembly_2026/manifest.tsv`; originals `salterns/MAG_contigs/` |

**Not committed as of 2026-11-18:** `annotation_runs/`, `workflow_envs/`, `manuscript_numbers*.py/json`, `comparison/` additions and `local/writeback_gtdbtk_r232.py`. Commit them, ideally together with this file, so a fresh clone has the code.

## 11. Minimal reproduction recipe

On a Linux host with ≥ 64 threads, ~2–3 TB scratch, and RAM as follows:

- **≤ 9 Gbp samples:** ≥ 256 GB. Read prep uses a 150g Java heap with peaks up to ~155 GB RSS, and SPAdes peaks at 30–67 GB.
- **The 17 samples > 9 Gbp:** ≥ ~300 GB, because of the 200g heap. Or run them on a larger machine, as was done here.

1. **Envs:** create each env from `workflow_envs/<env>.txt` (section 3).
   - Fix CheckM1's `DATA_CONFIG`.
   - Run DRAM with `PYTHONNOUSERSITE=1`.
2. **Databases:** fetch and md5-check GTDB R232, the CheckM2 DB, the GUNC proGenomes 2.1 DB, the CheckM1 2015_01_16 data and the KOfam 2026-10-26 pair (section 4).
   - `setup_gtdbtk.sh` and `fetch_gunc_db.sh` show the exact procedure.
   - Run `gtdbtk check_install`.
3. **Config:** copy `annotation_runs/rollout/scripts/` to a new root and edit the host paths in `make_config.py`. Those are `SCR`, `ROOT`, `ENVS` and the `TOOLS` dict; the binning-test reuse paths are only needed if you reuse reads.
   - Set `MIGRATED = []` for a single host, or use `make_config_local.py` for a second host.
   - Adjust the resource policy (`cores`, `cap_min/max`, `mem_floor_gb`, `kind_limits`, `reserve_gb`, `spades.mem_gb`) to the host, then generate `config.json`.
4. **Check:**
   - `driver.py plan --config scripts/config.json` shows the task DAG;
   - run the smoke configs on a small subsample (`smoke_prep.sh`);
   - `gen_task_scripts.py` lets you diff commands against `annotation_runs/rollout/README.md` "Exact commands".
5. **Run:** launch the driver detached at nice 10 (section 6) and watch `status.sh`. A sample takes about 6–20 h on a contended host, dominated by SPAdes.
6. **Afterwards:** section 7 in order, i.e. recovery/QC tables, species clustering, KOfam + Sec rescue + guilds, DRAM on species representatives, then the manuscript numbers.
