# DRAM annotation of all 12,611 wetland MAGs (KBase narrative 272312) on a server

**For:** a Claude instance with shell access to a Linux server, working for Andrew Freiburger.
**Written:** 2026-10-31 by the Claude session that classified these MAGs (GTDB-Tk 2.8.0 / GTDB R232) and runs the narrative's MAG workflow. The facts below were checked on that date. Commands for the other server are marked *not run here* where this session could not test them (no conda on Andrew's analysis server).
**Goal:** DRAM functional annotation and DRAM distillation for every MAG, handed back to Andrew. KBase stays read-only.

Read the whole file before acting. Each section ends with a **Check**; if a check fails, stop and report the exact error and what you tried.

---

## 1. Background

### The MAGs

KBase narrative **272312** ("wetland bins", workspace `freiburgermsu:narrative_1790959438780`) holds **12,611 MAGs**.
- They are medium-quality-or-better DAS Tool bins (DAS Tool SCG estimates: ≥ 50% complete, < 10% contaminated) from wetland metagenomes. The MEGAHIT assemblies are in narrative 261700.
- Each is a KBase Assembly object named `<bin id>__<SRR accession>`, for example `bin.001.fasta__SRR10127268`.
- Size: 0.27–10.5 Mbp, about 39 Gbp in total, so expect roughly **35–40 million genes**.

Already done:
- **Taxonomy:** GTDB R232 lineages for every MAG, from a local GTDB-Tk 2.8.0 run. These are in `taxonomy/`, and the same lineages are in KBase.
- **Quality:** KBase CheckM (lineage_wf, reduced tree) on every MAG, as 127 jobs of up to 100 MAGs. Results are in `quality/`; see 5c if it says PENDING.
- **DRAM on 5 MAGs:** KBase's DRAM app annotated the 5 "phylogenetically distinct" MAGs the narrative selected (`validation/selected_MAGs.txt`) and their 4 source metagenomes. Their output is your parity reference.

### Why this runs on a server

KBase's DRAM app was too slow for 12,611 MAGs:
- It annotates a job's MAGs **one after another** on a 32-CPU "extreme" node, with 30 threads.
- The **KOfam** search alone took **24 min for a 0.66 Mb MAG and 45 min for a 5.96 Mb MAG**, roughly 21 min plus 4 min per Mb.
- 12,611 MAGs would therefore need more than 6,000 node-hours on a scarce KBase queue.

A server can run many DRAM processes at once.

### What KBase's DRAM run did (the parity target)

- **Software:** the kb_DRAM app 0.1.2 (commit `935498b6`) with **DRAM 1.4.0**. The app's container installs prodigal, mmseqs2 (≠ 10.6d92c), HMMER (≠ 3.3.1), tRNAscan-SE ≥ 2, barrnap and scipy 1.8.1 via conda.
- **Annotation call:** `annotate_bins(fastas, out, min_contig_size=2500, trans_table='11', bit_score_threshold=60, rbh_bit_score_threshold=350, low_mem_mode=True, threads=30)`.
  - `low_mem_mode` means KOfam (not KEGG Genes) and no UniRef90.
  - prodigal ran in DRAM's default `meta` mode.
  - DRAM 1.4.0 renames genes `<MAG>_<contig>_<n>` regardless of `rename_bins`.
- **Distillation:** `summarize_genomes(annotations, trnas, rrnas, groupby_column='fasta')`.
- **Databases:** built by DRAM-setup on **2022-09-27**: KOfam, Pfam (mmseqs profile DB), dbCAN **v11**, MEROPS; UniRef none. VOGDB and RefSeq viral were configured but not searched for MAGs.

## 2. Ground rules

- **KBase is read-only for you.** Never save KBase objects or submit KBase jobs. A KBase writeback of DRAM results (KBase's app saves annotated Genome objects) is a separate decision for Andrew.
- **The KBase token is a password.** You only need one if you download the MAGs from KBase (5b). Never print, echo or log it, never put it on a command line, never copy it elsewhere. The tools read `~/.kbase_token` themselves.
  - If it's missing, ask Andrew to create it on this server: log in at narrative.kbase.us, copy the `kbase_session` cookie value, write it as the only line of `~/.kbase_token`, then `chmod 600 ~/.kbase_token`. He must never paste it into a chat.
- **One environment.** Create one dedicated conda/mamba environment. Don't change other environments or system packages, and don't use `sudo` without Andrew's OK.
  - If this server has its own instructions (a `CLAUDE.md`, admin rules) that restrict new environments, ask Andrew first.
- **Shared servers:**
  - Check `uptime`, `free -g` and `ps` before heavy work.
  - Run DRAM under `nice -n 10`.
  - Size the parallelism so other users keep headroom, and tell Andrew if you can't.
- Use one work folder, `~/dram_wetland`, plus the DRAM database folder you choose in section 7.

## 3. Choose a profile

| | **P: parity with KBase** | **R: recommended (default)** |
|---|---|---|
| DRAM | 1.4.0, KBase's version | **1.5.0**, the last DRAM 1 release (bug fixes since 1.4.0); DRAM 1 is no longer developed |
| KO assignment | KOfam (`--low_mem_mode`) | KEGG Genes **if Andrew has a KEGG FTP licence**, else KOfam |
| UniRef90 | no | no, by default (see optional) |
| Databases | `prepare_databases --skip_uniref`, today's releases (dbCAN v11 default) | same |
| Taxonomy / quality columns | none, as KBase | `--gtdb_taxonomy` (R232) and `--checkm_quality` (KBase CheckM) |
| Translation table | 11 for all | 11, plus table 4 for the 4 MAGs GTDB-Tk called table 4 |
| Throughput | — | many DRAM processes in parallel, 25 MAGs each |

**Optional extras for R,** only with Andrew's agreement:
- **KEGG Genes:** `prepare_databases --kegg_loc kegg.pep`, then annotate without `--low_mem_mode`. Uses about 100 GB RAM per annotate process instead of under 50 GB.
- **UniRef90** (`--use_uniref`):
  - setup needs ~512 GB RAM and ~500 GB disk;
  - annotation needs ~220 GB RAM per process and is "drastically" slower.
  - For 12,611 MAGs it's rarely worth it; consider it for a subset only.
- **VOGDB** (`--use_vogdb`): viral orthologous groups. That's another ~25,000-model HMM search per MAG, roughly doubling the time.
- **CAMPER** (`--use_camper`, DRAM 1.5.0): polyphenol metabolism, relevant for wetlands. DRAM's own code labels it beta. Setup doesn't download it: get `CAMPER v1.0.0-beta.1` from `https://github.com/WrightonLabCSU/CAMPER/archive/refs/tags/v1.0.0-beta.1.tar.gz`, then run `DRAM-setup.py set_database_locations --camper_tar_gz_loc <file>`.

Keep these exactly as KBase and DRAM's defaults have them: prodigal mode (`meta`), `--min_contig_size 2500`, `--bit_score_threshold 60`, `--rbh_bit_score_threshold 350`.

**Why not DRAM2?**
- DRAM2, the Nextflow rewrite, parallelizes well, but its output (`raw-annotations.tsv`) differs from DRAM 1's `annotations.tsv`.
- The narrative's cell C and Andrew's DRAM analyses read DRAM 1 columns (`ko_id`/`kegg_id`, `cazy_*`, `peptidase_family`).
- If this server already has DRAM2 (KBUtilLib's `DRAM2Utils` is pinned to an install on a host called "h100"), mention it to Andrew, but use DRAM 1 unless he says otherwise.

**Check:** you know which profile you're running, and Andrew has agreed to any optional extra.

## 4. Preflight

```bash
mkdir -p ~/dram_wetland && cd ~/dram_wetland
nproc; free -g; df -h ~ .; uptime
command -v mamba || command -v conda || echo "no conda"
for u in https://pypi.org https://conda.anaconda.org https://raw.githubusercontent.com https://bcb.unl.edu \
         https://ftp.ebi.ac.uk https://www.genome.jp; do printf '%-36s ' $u; curl -s -o /dev/null -m 20 -w '%{http_code}\n' $u; done
curl -s -m 30 --list-only ftp://ftp.genome.jp/pub/db/kofam/ | head -3        # KOfam is fetched over FTP
```

| Resource | Needed |
|---|---|
| CPU | Any; more is faster. Throughput comes from running many DRAM processes in parallel. |
| RAM | Database setup ~128 GB with `--skip_uniref` (~512 GB with UniRef90). Annotation < 50 GB per process with KOfam (~100 GB with KEGG Genes, ~220 GB with UniRef90). Measure in the pilot. |
| Disk | Databases ~30 GB (with UniRef90 ~500 GB); MAG FASTAs 37 GB; DRAM output ~15–20 MB per MAG, ~200–250 GB in total; merged copies extra. Put databases and output on local disk. |
| Network | HTTPS to PyPI, conda channels, GitHub and bcb.unl.edu (dbCAN); **FTP** to ftp.genome.jp (KOfam) and ftp.ebi.ac.uk (Pfam, MEROPS). |

**If there's no conda,** install Miniforge in your home folder; no root is needed:
```bash
curl -L -o /tmp/miniforge.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash /tmp/miniforge.sh -b -p ~/miniforge3 && source ~/miniforge3/bin/activate
```

**Check:** resources meet the table (or Andrew accepts a reduced setup), and the hosts answer, FTP included. If FTP is blocked, see the Troubleshooting table.

## 5. Inputs

### 5a. The handoff bundle (small)

On Andrew's analysis server the bundle is `~/gtdbtk_kbase/dram_handoff.tar.gz`, with `.sha256` beside it. Copy both here; use `scp`/`rsync` if this server can reach that one, otherwise ask Andrew to transfer them. Then:
```bash
cd ~/dram_wetland && sha256sum -c dram_handoff.tar.gz.sha256
tar -xzf dram_handoff.tar.gz --strip-components=1 && sha256sum -c SHA256SUMS
```
It contains:
- `RUNBOOK.md` (this file);
- `mags/all_mags.tsv` (for a KBase download);
- `mags/mag_manifest.tsv`: name, KBase set, size, contigs, FASTA md5, KBase ref;
- `taxonomy/`: GTDB R232 summaries and `translation_table_4.txt`;
- `quality/`;
- `validation/`;
- `tools/`.

### 5b. The MAG FASTAs (37 GB)

- **Option A (preferred if reachable):** rsync from Andrew's analysis server, `~/gtdbtk_kbase/genomes/`. That folder holds exactly these 12,611 files as `<name>.fna`.
- **Option B:** download from KBase with the bundled tools. This needs `~/.kbase_token`; see section 2.
  - `tools/kbase_mags_kbaseus.py` wraps `kbase_mags.py` and fetches from `kbase.us`. Use it if `services.kbase.us` is unreachable from here, as it is from Andrew's analysis server.
  ```bash
  python3 tools/kbase_mags.py fetch --list mags/all_mags.tsv --genome-dir genomes --manifest genomes_manifest.tsv
  ```
  `fetch` checks every file's md5 against KBase and resumes if interrupted. To go faster, split `mags/all_mags.tsv` into 4 lists (keep the header line) and run 4 fetches in parallel.

Either way, verify against the manifest:
```bash
python3 - <<'EOF'
import csv, hashlib, os
bad = []
for r in csv.DictReader(open("mags/mag_manifest.tsv"), delimiter="\t"):
    p = f"genomes/{r['name']}.fna"
    if not os.path.isfile(p) or hashlib.md5(open(p, "rb").read()).hexdigest() != r["md5"]:
        bad.append(r["name"])
print(f"{len(bad)} missing/mismatched FASTAs", bad[:5])
EOF
```
**Check:** 0 missing or mismatched, and `ls genomes/*.fna | wc -l` = 12611.

### 5c. CheckM quality table

`quality/checkm_quality_for_DRAM.tsv` has the columns `Bin Id`, `Completeness` and `Contamination` that DRAM's `--checkm_quality` reads.

If you find `quality/PENDING.txt` instead, KBase's CheckM jobs were still running when the bundle was built. Either:
- ask Andrew for a new bundle;
- or, with a KBase token, run `python3 tools/collect_checkm.py --wait` in `~/dram_wetland/quality`. It is read-only and waits until all 127 jobs are final.

DRAM can start without it: leave `CHECKM_Q` empty, and add the two columns to the merged annotations later with a join on `fasta`.

**Check:** the quality table has 12,611 rows, or you have a plan for adding it later.

## 6. Build the DRAM environment

**Profile R** (*not run here*): the Bioconda recipe for DRAM 1.5.0 pins `pandas >=1.5,<2`, `scikit-bio >=0.5.8,<0.6`, `mmseqs2 >10.6d92c` and `trnascan-se >=2`.
```bash
mamba create -y -n dram-1.5.0 -c conda-forge -c bioconda "dram=1.5.0" "sqlalchemy<2"
eval "$(conda shell.bash hook)" && conda activate dram-1.5.0
```

**Profile P:** Bioconda's DRAM 1.4.0 recipe pins only `pandas >=1` and `scipy <=1.8.1`. Add the pins yourself, or pandas 2 gets installed and DRAM breaks:
```bash
mamba create -y -n dram-1.4.0 -c conda-forge -c bioconda "dram=1.4.0" python=3.10 "pandas>=1.5,<2" \
    "scikit-bio>=0.5.7,<0.6" "sqlalchemy<2"
```
- If either solve fails, create the environment from the dependencies and install DRAM without letting pip change them:

  ```bash
  mamba create -y -n dram-1.4.0 -c conda-forge -c bioconda python=3.10 "pandas>=1.5,<2" "numpy<2" "scikit-bio>=0.5.8,<0.6" "sqlalchemy<2" "altair>=4,<5" networkx openpyxl prodigal "mmseqs2>10.6d92c" "hmmer>=3.3.2" "trnascan-se>=2" barrnap parallel ruby wget curl
  conda activate dram-1.4.0 && pip install --no-deps DRAM-bio==1.4.0
  ```

These pins matter:
- DRAM 1.4.0 and 1.5.0 call `DataFrame.iteritems()`, which pandas 2 removed.
- Their database code uses SQLAlchemy 1.x's `declarative_base`.

Record the environment:
```bash
mkdir -p env && conda env export > env/dram_env.yml && conda list --explicit > env/dram_env.explicit.txt
DRAM-setup.py version; python -c "import pandas, sqlalchemy; print(pandas.__version__, sqlalchemy.__version__)"
for t in prodigal mmseqs hmmsearch tRNAscan-SE barrnap; do printf '%-12s ' $t; command -v $t || echo MISSING; done
ls -l /usr/bin/time || echo "no GNU time: run_dram_batch.sh still works, without the per-batch time/RSS summary"
```
**Check:** `DRAM-setup.py version` prints 1.5.0 (or 1.4.0 for P). pandas is below 2, SQLAlchemy is below 2, and every tool is found.

## 7. Databases

```bash
export DRAMDB=/path/on/local/disk/dram_db && mkdir -p "$DRAMDB"
nice -n 10 DRAM-setup.py prepare_databases --output_dir "$DRAMDB" --skip_uniref --threads 16 --verbose 2>&1 | tee env/dram_setup.log
#   with a KEGG licence add:  --kegg_loc /path/to/kegg.pep        (UniRef90: leave out --skip_uniref; see section 3)
DRAM-setup.py print_config | tee env/dram_config.txt
```
Setup can take several hours. DRAM's README says up to 5 h with UniRef90, and about 2 h of processing on 10 cores when nothing needs downloading. If it fails midway, re-run the same command.

**Check:** `print_config` shows paths for these entries, with download dates:
- KOfam db and KOfam KO list (or KEGG db);
- Pfam db and Pfam hmm dat;
- dbCAN db with family activities and subfamily EC;
- MEROPS peptidase db;
- the Description db;
- the five distillation sheets (genome summary, module step, ETC module, function heatmap, AMG).

UniRef should show `None` unless you chose it. Save `env/dram_config.txt` for the report.

## 8. Pilots (before the full run)

Write the settings file. Every later step sources it.
```bash
cd ~/dram_wetland && cat > dram_env.sh <<'EOF'
eval "$(conda shell.bash hook)" && conda activate dram-1.5.0      # dram-1.4.0 for profile P
export DRAM_WORK=$HOME/dram_wetland
export BATCH_DIR=$DRAM_WORK/batches DRAM_OUT=$DRAM_WORK/dram_out DRAM_LOGS=$DRAM_WORK/logs
export PROFILE=recommended                     # "parity" adds --low_mem_mode (KOfam, no UniRef) as KBase did
export THREADS=2                               # threads per DRAM process: set from the timing pilot
export TAX_BAC=$DRAM_WORK/taxonomy/gtdbtk.bac120.summary.tsv TAX_ARC=$DRAM_WORK/taxonomy/gtdbtk.ar53.summary.tsv
export CHECKM_Q=$DRAM_WORK/quality/checkm_quality_for_DRAM.tsv
export DRAM_EXTRA=""                           # optional extras from section 3, e.g. "--use_camper"
export PY=python3
EOF
source dram_env.sh
python3 tools/make_batches.py --manifest mags/mag_manifest.tsv --genomes genomes --out batches --batch-size 25
```
- **Profile P:** set `PROFILE=parity`, `TAX_BAC=` and `CHECKM_Q=` (empty), so the run matches KBase's call.
- **Profile R** uses KOfam unless you set up KEGG Genes. If you did, keep `PROFILE=recommended` and don't add `--low_mem_mode`.

`make_batches.py` creates, all as symlinks into `genomes/`:
- `batches/partNNN_kk/`: about 510 batches of up to 25 MAGs, grouped like KBase's `MAG_all.partNNN` sets;
- `batches/tt4/`: the 4 table-4 MAGs;
- `batches/pilot_selected5/` and `batches/pilot_timing/` (24 MAGs across the size range).

**8a. Validation against KBase.** `tools/run_dram_batch.sh` runs `DRAM.py annotate` with the profile's settings on one batch folder:
```bash
nice -n 10 bash tools/run_dram_batch.sh pilot_selected5
```
If `validation/kbase_dram_selected5_annotations.tsv` exists, compare. Gene IDs are `<MAG>_<contig>_<n>` in both runs, so with the same prodigal version the genes should be identical:
```bash
python3 - <<'EOF'
import csv, re, sys
from collections import defaultdict
csv.field_size_limit(sys.maxsize)
KO, CZ = re.compile(r"K\d{5}"), re.compile(r"\b(?:GH|GT|PL|CE|AA|CBM)\d+\b")
def load(p):
    g = {}
    for r in csv.DictReader(open(p), delimiter="\t"):
        gid = r.get("") or next(iter(r.values()))
        g[gid] = (r.get("fasta", ""), set(KO.findall((r.get("ko_id") or "") + " " + (r.get("kegg_id") or ""))),
                  set(CZ.findall(" ".join(r.get(c) or "" for c in ("cazy_ids", "cazy_hits", "cazy_id")))),
                  (r.get("peptidase_family") or "").strip())
    return g
kb, me = load("validation/kbase_dram_selected5_annotations.tsv"), load("dram_out/pilot_selected5/annotations.tsv")
per = defaultdict(lambda: [0, 0, 0, 0, 0])
for gid in set(kb) | set(me):
    m = (kb.get(gid) or me.get(gid))[0]
    s = per[m]; s[0] += gid in kb; s[1] += gid in me
    if gid in kb and gid in me:
        s[2] += kb[gid][1] == me[gid][1]; s[3] += kb[gid][2] == me[gid][2]; s[4] += kb[gid][3] == me[gid][3]
for m, (a, b, ko, cz, pep) in sorted(per.items()):
    both = min(a, b) or 1
    print(f"{m}: genes KBase {a} / here {b}; same KO {ko/both:.1%}, same CAZy {cz/both:.1%}, same MEROPS {pep/both:.1%}")
EOF
```
- **Expect:** identical gene counts and very high agreement, typically above 90%. Differences come from database releases (2022 versus now) and, for profile R, from any KEGG Genes or UniRef90 use.
- **If gene counts differ,** compare prodigal versions (KBase's is in its container; yours is `prodigal -v`) before anything else.
- If the KBase file isn't in the bundle yet, KBase's DRAM job on the 5 MAGs was still running. Note that and continue; validate when Andrew sends it.

**8b. Timing.**
```bash
nice -n 10 bash tools/run_dram_batch.sh pilot_timing
grep -E "Elapsed|Maximum resident|Percent of CPU" logs/pilot_timing.log
```
Take minutes per MAG = elapsed ÷ 24, and peak RSS. If cores allow, repeat on a copy of `pilot_timing` with `THREADS=4` to see whether more threads per process helps; usually it doesn't help much. Then compute:
- **Parallel processes:** `P = min(floor(usable cores ÷ THREADS), floor(0.8 × free RAM GB ÷ peak RSS GB))`
- **Estimated days:** `12,611 × minutes per MAG ÷ P ÷ 1440`

For reference, KBase's 30-thread run took 26 and 48 min for the first two MAGs; KOfam dominates.

**Check:**
- The 5-MAG comparison is acceptable, or explained.
- You know P, THREADS and the estimated duration.
- **Report these numbers to Andrew before the full run** if the estimate is over about 2 weeks, or if the run would take most of a shared machine.

## 9. Full run

```bash
cd ~/dram_wetland && source dram_env.sh
ls batches | grep -v '^pilot_' > batch_list.txt && wc -l batch_list.txt          # ~511 batches
P=<from 8b>
nohup bash -c "xargs -a batch_list.txt -P $P -I{} nice -n 10 bash tools/run_dram_batch.sh {}" >> run_all.log 2>&1 &
echo $! > run_all.pid
```
- `run_dram_batch.sh` writes one line per batch to `run_all.log`: done, already done, or FAILED with the log path.
- Each batch's DRAM log, with a `/usr/bin/time -v` summary, is in `logs/<batch>.log`.
- **Monitor** every 30 min or so (not with `tail -f`):

  ```bash
  echo "done $(grep -c ': done' run_all.log)  failed $(grep -c FAILED run_all.log)  of $(wc -l < batch_list.txt)"; free -g | sed -n 2p
  ```
- **Failures:** read the batch log, fix the cause, and re-run the same `xargs` command. Finished batches are skipped, and a failed batch's partial output is removed before its retry.
- **Memory:** if the machine runs short, stop new launches (`kill $(cat run_all.pid)`; running batches finish), then restart with a smaller P.

**Check:**
- Every name in `batch_list.txt` has `dram_out/<batch>/annotations.tsv`.
- The "MAGs with genes" counts in `run_all.log` add up to 12,611. DRAM drops contigs under 2,500 bp, so a MAG with only short contigs would have no genes; list any.

## 10. Merge and summarize

DRAM's own `merge_annotations` and `distill` load every annotation into pandas at once. For ~37 million genes that needs a very large amount of RAM. So merge by streaming, and distill per batch.
```bash
cd ~/dram_wetland && source dram_env.sh && mkdir -p merged distill
dirs=$(sed 's#^#dram_out/#' batch_list.txt)           # the real batches only: the pilot outputs would duplicate MAGs
for f in annotations trnas rrnas; do python3 tools/merge_tsv.py merged/$f.tsv $(for d in $dirs; do echo "$d/$f.tsv"; done); done
for f in genes.faa genes.fna scaffolds.fna; do for d in $dirs; do cat "$d/$f"; done > merged/$f; done   # IDs carry the MAG name
for d in $dirs; do cat "$d/genes.gff"; done | awk 'NR > 1 && /^##gff-version/ {next} 1' > merged/genes.gff
# distill refuses an existing output folder: remove distill/<batch> before re-running a failed one
xargs -a batch_list.txt -P "$P" -I{} sh -c 'DRAM.py distill -i dram_out/{}/annotations.tsv -o distill/{} --groupby_column fasta \
    $( [ -s dram_out/{}/trnas.tsv ] && echo --trna_path dram_out/{}/trnas.tsv ) \
    $( [ -s dram_out/{}/rrnas.tsv ] && echo --rrna_path dram_out/{}/rrnas.tsv ) > logs/distill_{}.log 2>&1 || echo "distill FAILED {}"'
for f in genome_stats product; do python3 tools/merge_tsv.py merged/$f.tsv $(sed "s#^#distill/#; s#\$#/$f.tsv#" batch_list.txt); done
```
- `metabolism_summary.xlsx` and `product.html` stay per batch, under `distill/<batch>/`.
- If this server has at least about 5 times the size of `merged/annotations.tsv` in free RAM, you may also run one `DRAM.py distill` on the merged files.

**Check:**
- `merged/annotations.tsv` has one distinct `fasta` value per MAG with genes.
- `merged/genome_stats.tsv` has the same number of rows.
- No "distill FAILED" lines were printed.

## 11. Deliverables and handback

Assemble `~/dram_wetland/dram_results/`:
- **`merged/`:** `annotations.tsv.gz`, `trnas.tsv`, `rrnas.tsv`, `genes.faa.gz`, `genes.fna.gz`, `genes.gff.gz`, `genome_stats.tsv`, `product.tsv`.
- **`distill/`:** per-batch `metabolism_summary.xlsx` and `product.html`.
- **`logs/`:** `run_all.log`, the per-batch logs and the setup log.
- **`env/`:** the environment export, `dram_config.txt` and the DRAM version.
- **`validation_report.md`:** the 8a comparison and the 8b timing, plus the actual run time.
- **`SHA256SUMS`.**

Per-batch `genbank/` folders and `scaffolds.fna` are large; keep them on this server unless Andrew wants them. Transfer `dram_results/` to Andrew's analysis server (`~/gtdbtk_kbase/dram_results/`) if it's reachable, otherwise where he says. Don't delete inputs or outputs until he confirms receipt.

## 12. What to report to Andrew

- **Setup:** the server (cores, RAM), the profile and any extras, the DRAM version, and the database releases with dates from `print_config`.
- **Validation:** gene counts and percent same KO, CAZy and MEROPS per MAG, against KBase.
- **Timing:** minutes per MAG, P and THREADS, the estimate, and the actual wall time.
- **Results:**
  - MAGs annotated and genes total;
  - the share of genes with a KO;
  - the share with any of KO, CAZy, MEROPS or Pfam;
  - MAGs without genes;
  - failed batches and how they were resolved.
- **Location:** where the results are, and their total size.

## Troubleshooting

| Symptom | Fix |
|---|---|
| `'DataFrame' object has no attribute 'iteritems'` | pandas 2 got installed. Pin `pandas<2` (DRAM 1.4/1.5 need it); install DRAM with `pip --no-deps` if pip upgraded pandas. |
| SQLAlchemy errors about `declarative_base` or sessions | Pin `sqlalchemy<2`. |
| `mmseqs: Illegal instruction` | The mmseqs2 build doesn't match the CPU. Install the build for its instruction set (SSE4.1 or AVX2). |
| `prepare_databases` download errors, or FTP blocked | Re-run with the same `--output_dir`. If FTP is blocked, download KOfam (`profiles.tar.gz`, `ko_list.gz`), Pfam and MEROPS manually over HTTPS mirrors, and pass `--kofam_hmm_loc`, `--kofam_ko_list_loc` and so on. Ask Andrew before using third-party mirrors. |
| DRAM: output directory already exists | `run_dram_batch.sh` removes partial output itself. If you run DRAM by hand, remove the directory first. |
| Processes killed for memory | Lower P. Check peak RSS in `logs/*.log`. Without KEGG Genes or UniRef, a process should stay well under 50 GB. |
| tRNAscan-SE / infernal errors | Use tRNAscan-SE ≥ 2 from Bioconda in the same environment. |
| A batch keeps failing on one MAG | Move that MAG into its own batch folder, run it alone, and report it with its log. |
