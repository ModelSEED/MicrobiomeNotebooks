# GTDB-Tk 2.8.0 / GTDB R232 environment, and KBUtilLib's GTDB skani support

**For:** a Claude instance working for Andrew Freiburger on Linux. You can be on his server, where everything below already exists, or on a new machine.
**Written:** 2026-10-30, by the Claude session that built and used this environment. Every command and expected output below was checked on Andrew's server on that date, except where a step is marked *not re-run*.

Read this whole file before acting. Each section ends with a **Check**; if a check fails, stop and report the exact output instead of working around it.

---

## 1. What exists

| Item | Location (Andrew's server) | Notes |
|---|---|---|
| GTDB-Tk environment | `~/Documents/gtdbtk_venv` | uv venv, CPython 3.12.3, `gtdbtk==2.8.0` + 5 pinned deps; 4 Bioconda tools copied into `bin/` |
| GTDB R232 reference package | `~/Documents/gtdbtk_data/release232` | 94 GiB unpacked; `skani/` alone is 80 GiB (199,923 species-representative sketches) |
| Package tarball | `~/Documents/gtdbtk_data/gtdbtk_r232_data.tar.gz` | 60,806,405,195 bytes, md5 `25a59e0352b1fd150c589f56559767d4` (matches GTDB's `MD5SUM.txt`). Can be deleted; GTDB-Tk uses only the unpacked folder. |
| Original build script | `salterns/run_gtdbtk_local.sh` in the MicrobiomeNotebooks repo (tracked) | Built the venv and data on 2026-10-20; resumable |
| KBUtilLib GTDB skani code | PR [cshenry/KBUtilLib#105](https://github.com/cshenry/KBUtilLib/pull/105); branch `skani-gtdb-database` on the fork `freiburgermsu/KBUtilLib`, commit `c3a279f` | Checked out at `~/Documents/KBUtilLib_skani_gtdb` (git worktree, remote `fork`) |
| skani registry entry | `~/.kbutillib/skani_databases.json`, key `gtdb_r232` | Already registered on Andrew's server |

What it has been used for:
- **Salterns, 2026-10-20:** 310 dereplicated MAGs, in `salterns/MAG_QC/gtdbtk_local_r232/`.
- **Wetland, 2026-10-29:** all 2,449 MAGs of KBase narrative 272312, in `~/gtdbtk_kbase/gtdbtk_all/`. The lineages were then written back to KBase.
- **Why locally:** KBase's own GTDB-Tk app failed on these sets. Its release version crashes while relaying GTDB-Tk's progress-bar output, and its beta version has a fastANI glibc error.

## 2. Ground rules

- **Venvs need Andrew's approval.** His global `~/.claude/CLAUDE.md` says the only Python environment is `~/Documents/py_venv`, and creating any other venv needs his approval first.
  - `gtdbtk_venv` is the approved exception: approved 2026-10-20 for salterns and reused for wetland on 2026-10-29.
  - If you need to rebuild it, or build it at another path or on another machine, ask Andrew first.
- **Never install `gtdbtk` into `py_venv`.** GTDB-Tk 2.8.0 pins numpy 1.x and pydantic 1.x, while `py_venv` has numpy 2.x and pydantic 2.x, and KBUtilLib requires `pydantic>=2`.
  - The two stay separate on purpose. KBUtilLib runs in `py_venv` and calls the `skani` binary from `gtdbtk_venv` as an external program.
- **Never put `gtdbtk_venv/bin` on PATH globally,** whether in `.bashrc` or anywhere else. It contains `python`, `python3` and pip, which would shadow `py_venv`.
  - Prepend it only in the shell or command that runs `gtdbtk` (see section 4d).
- **Check memory before big jobs.** Andrew's server is shared: 64 cores, 251 GB RAM, and many Claude sessions run jobs there.
  - Check `free -g` and the largest processes (`ps -eo pid,rss,args --sort=-rss | head`) before starting anything big.
  - Run GTDB-Tk inside a memory-capped systemd scope (section 5).
  - Tell busy peer sessions (`ListAgents`, then `SendMessage`) before using more than about 100 GB.
- **The KBase token is a password.** Some commands below don't need it at all. Never print it, and never put it on a command line or in a file.

## 3. On Andrew's server: verify instead of rebuilding

```bash
export PATH="$HOME/Documents/gtdbtk_venv/bin:$PATH"
export GTDBTK_DATA_PATH="$HOME/Documents/gtdbtk_data/release232"
gtdbtk --version | head -1
uv pip freeze --python ~/Documents/gtdbtk_venv/bin/python
for b in prodigal hmmsearch hmmalign pplacer guppy skani; do printf '%-9s ' $b; command -v $b; ldd "$(command -v $b)" | grep -c 'not found'; done
gtdbtk check_install
```

**Check:**
- The first line reads `gtdbtk: version 2.8.0 …`.
- The freeze shows exactly these packages:

  ```
  dendropy==4.6.4  gtdbtk==2.8.0  numpy==1.26.4  pydantic==1.10.26  tqdm==4.70.1  typing-extensions==4.16.0
  ```
- Every binary resolves to `~/Documents/gtdbtk_venv/bin/…` with `0` unresolved libraries.
  - `pplacer` and `guppy` also print `not a dynamic executable`. That's expected, because they are static builds.
- `check_install` (about 70 s) reports every reference folder `OK` and ends with `Done.`.
  - `FastTree` and `FastTreeMP` show `NOT FOUND`. That's expected: only `de_novo_wf` needs them.

If all of these hold, skip to section 5 or 6.

## 4. Rebuild from scratch (needs Andrew's approval; see section 2)

These are the steps `salterns/run_gtdbtk_local.sh` ran on 2026-10-20. The commands were *not re-run* today. Their inputs were re-checked today: all four Bioconda URLs return HTTP 200, the R232 URL returns the expected size, and the package versions resolve.

**Prerequisites:**
- Linux x86_64. The build was tested on Ubuntu 24.04, glibc 2.39.
- `uv`, `curl`, `tar`, `bzip2`, `unzip` and `zstd`. No conda is needed.
- Disk: about 152 GiB at peak (57 GiB tarball plus 94 GiB unpacked), or about 94 GiB if you stream-extract (section 4c, option B).
- RAM: at least 64 GB to run anything. About 150 GB free for large batches at `--pplacer_cpus 8` (section 5).

### 4a. Python environment

```bash
VENV="$HOME/Documents/gtdbtk_venv"
uv venv "$VENV" --python 3.12
# exact pins = the working environment; unpinned, pydantic would resolve to 1.10.27 today
uv pip install --python "$VENV/bin/python" gtdbtk==2.8.0 dendropy==4.6.4 numpy==1.26.4 \
    pydantic==1.10.26 tqdm==4.70.1 typing-extensions==4.16.0
```

### 4b. External tools: Bioconda packages, unpacked without conda

```bash
VENV="$HOME/Documents/gtdbtk_venv"; BIOCONDA=https://conda.anaconda.org/bioconda/linux-64
for f in prodigal-2.6.3-h516909a_2.tar.bz2 hmmer-3.4-hdbdd923_2.tar.bz2 \
         pplacer-1.1.alpha19-h9ee0642_2.tar.bz2 skani-0.4.0-hec9b1f2_0.conda; do
  name="${f%%-*}"; dir="$VENV/pkgs/$name"
  [ -f "$dir/.done" ] && continue
  rm -rf "$dir" && mkdir -p "$dir"
  curl -sSfL -o "$dir/$f" "$BIOCONDA/$f"
  case "$f" in
    *.tar.bz2) tar -xjf "$dir/$f" -C "$dir" ;;
    *.conda)   mkdir -p "$dir/z" && unzip -q -o "$dir/$f" -d "$dir/z"
               for t in "$dir"/z/pkg-*.tar.zst; do zstd -dq --stdout "$t" | tar -x -C "$dir"; done ;;
  esac
  cp -f "$dir"/bin/* "$VENV/bin/"
  touch "$dir/.done"
done
```

Keep pplacer at **1.1.alpha19**, the static build. 1.1.alpha23 needs `libgsl.so.25`, which Ubuntu 24.04 doesn't have.

### 4c. GTDB R232 reference package (60.8 GB download)

```bash
DATA="$HOME/Documents/gtdbtk_data"; mkdir -p "$DATA" && cd "$DATA"
URL=https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz
# Option A (resumable; ~152 GiB peak): download, verify, extract
curl -sSfL -C - -o gtdbtk_r232_data.tar.gz "$URL"
stat -c %s gtdbtk_r232_data.tar.gz                                        # must be 60806405195
echo "25a59e0352b1fd150c589f56559767d4  gtdbtk_r232_data.tar.gz" | md5sum -c -   # must say OK (~3 min)
tar -xzf gtdbtk_r232_data.tar.gz                                          # creates release232/
# Option B (~94 GiB; not resumable): stream-extract and hash on the way through
# curl -sSfL "$URL" | tee >(md5sum > r232.md5) | tar -xz ; cat r232.md5   # must be 25a59e0352b1fd150c589f56559767d4
```

The md5 comes from GTDB's release checksum file, `https://data.gtdb.ecogenomic.org/releases/release232/232.0/MD5SUM.txt`. There is no checksum file next to the tarball. Delete the tarball only after section 4e passes, and only if the disk space is needed.

### 4d. How to activate (per shell, never globally)

```bash
export PATH="$HOME/Documents/gtdbtk_venv/bin:$PATH"
export GTDBTK_DATA_PATH="$HOME/Documents/gtdbtk_data/release232"
```

### 4e. Verify

Run all of section 3, then the bundled test, which takes about 6 min. Its archaeal pplacer step uses about 44 GB RAM at 8 threads.

```bash
gtdbtk test --out_dir /tmp/gtdbtk_selftest --cpus 8 ; echo "exit=$?"
```

**Check:** the output contains `Test has successfully finished.` and ends with `exit=0`.

## 5. Running `classify_wf`

Batchfile format: one genome per line, `<absolute FASTA path><TAB><genome name>`, with no header.

```bash
RUN=myrun; OUT="$PWD/gtdbtk_$RUN"; TMP=/tmp/gtdbtk_$RUN; mkdir -p "$TMP"
setsid nohup systemd-run --user --scope --unit=gtdbtk-$RUN -p MemoryMax=170G -p MemorySwapMax=0 \
  bash -c "export PATH=$HOME/Documents/gtdbtk_venv/bin:\$PATH GTDBTK_DATA_PATH=$HOME/Documents/gtdbtk_data/release232 TMPDIR=$TMP;
           gtdbtk classify_wf --batchfile batchfile.tsv --out_dir $OUT --cpus 32 --pplacer_cpus 8 --min_perc_aa 10 --tmpdir $TMP" \
  > gtdbtk_$RUN.log 2>&1 < /dev/null &
```

Results are in `$OUT/gtdbtk.bac120.summary.tsv` and `$OUT/gtdbtk.ar53.summary.tsv`, which are links into `classify/`. `$OUT/gtdbtk.json` records the version, command and steps. Check the log with `tail`, not `tail -f`.

Lessons learned, all of them hit in real runs:
- **`--tmpdir` must already exist,** or GTDB-Tk refuses it.
- **Keep the tmpdir path short,** for example `/tmp/gtdbtk_<run>`. `multiprocessing` creates Unix sockets inside it, and socket paths are capped at 107 bytes. A tmpdir nested in a long output path made Prodigal die with `OSError: AF_UNIX path too long`.
- **2.8.0 always runs a skani ANI screen.** There is no `--mash_db` and no `--skip_ani_screen` (those belong to 2.3.x), and fastANI and Mash are not used. Genomes whose closest skani hit passes AF ≥ 0.5 and the species' ANI radius are classified at that step; the rest go to tree placement.
- **Memory is set by pplacer.** Measured on the 2,449-MAG run at `--pplacer_cpus 8`, sampled every 30 s by `mem_guard.sh`:

  | Stage | pplacer RSS | Whole job (anonymous) |
  |---|---|---|
  | Archaeal tree | 42 GB | 51 GB |
  | Bacterial backbone | 29 GB | 38 GB |
  | Class-level trees | **133 GB** | **145 GB** |

  - The 160 GB cap held, but with only about 15 GB to spare. Use `MemoryMax=170G` or more if the machine has room.
  - Otherwise lower `--pplacer_cpus`, since memory grows with pplacer threads, and add `--scratch_dir <dir>`. This is slower but uses far less RAM.
  - If pplacer is killed, rerun with `--pplacer_cpus 1–2` and `--scratch_dir`.
- **Runtime is mostly a fixed cost.** The 8 bacterial class-level trees take about 12–22 min each, almost regardless of how many genomes go into each.
  - 2,449 MAGs took 3.4 h with 32 CPUs (1,103 classified at the ANI screen).
  - 310 MAGs took about 3 h with 48 CPUs and 16 pplacer threads.

  One combined run is much cheaper than several small ones.
- **R232 names differ from KBase's R214.** KBase's GTDB-Tk app uses GTDB-Tk 2.3.2 with R214, so the same organism can carry a different name. In the wetland MAGs:
  - 974 of 2,249 lineages were identical, and 1,531 kept the same phylum and genus;
  - Polyangia moved from `p__Myxococcota` to `p__Myxococcota_A`;
  - `p__Desulfobacterota_E` was merged into `p__Desulfobacterota`;
  - many placeholder genera were retired.

  Don't mix R214 and R232 labels in one analysis without saying so.

## 6. KBUtilLib: GTDB species assignment from the package's skani database

### What the code does

The code lives in `kbutillib.domains.genome.skani_utils.SKANIUtils`. Upstream moved it there; the old `kbutillib.skani_utils` still works but is a deprecation shim.

- **`register_gtdb_database(package_dir, database_name=None, overwrite=False)`**
  - Checks that the package has `skani/database/{sketches.db,index.db,markers.bin}`, `taxonomy/gtdb_taxonomy.tsv` and `radii/gtdb_radii.tsv`.
  - Reads the release from `metadata/metadata.txt`.
  - Records a `"type": "gtdb"` entry in `~/.kbutillib/skani_databases.json`. The default name is `gtdb_<release>`, e.g. `gtdb_r232`. Nothing is copied.
  - Returns `registered: False` if the entry already exists; pass `overwrite=True` if the package has moved.
- **`classify_genomes_gtdb(fastas, database_name="gtdb_r232", threads=1, af_threshold=0.5, max_hits=5, timeout=None)`**
  - Runs `skani search` against the package database, then applies GTDB-Tk's ANI-screen rule:
    1. AF = max(AF_ref, AF_query).
    2. Among hits with AF ≥ 0.5, take the highest ANI.
    3. Assign that representative's species if ANI ≥ its radius in `gtdb_radii.tsv`.
  - Returns `{genome_id: {...}}`, where `genome_id` is the file name minus its FASTA extension and `.gz`. Each record has:
    - `status`: `species assigned`, `below radius`, `low AF` or `no hit`;
    - `classification`: the GTDB lineage, set only when a species is assigned;
    - `reference`, `ani` (percent), `af` (fraction), `radius` and `reference_taxonomy`;
    - `hits`: the closest `max_hits` hits.
- **`canonical_gtdb_id(gid)`:** GTDB-Tk's accession canonicalization (`GB_GCA_…`, `RS_GCF_…` and `…_genomic.fna.gz` all map to `G…`).
- **`SKANIUtils(skani_executable=...)`:** points at a skani binary that isn't on PATH, such as the one in `gtdbtk_venv`.

What it is **not**: it is not a replacement for `classify_wf`. It reproduces only GTDB-Tk's first step. Genomes with any status other than `species assigned` need tree placement (section 5). In the wetland set, 1,103 of 2,449 MAGs were assigned at the ANI screen.

Validated on 2026-10-30 against the wetland `classify_wf` run:
- 30 of 30 randomly sampled ANI-screen MAGs got GTDB-Tk's exact species lineage.
- 30 of 30 tree-placed MAGs were correctly left unassigned (4 below radius, 11 low AF, 15 no hit).

Cost grows with the number of queries: 60 genomes took 3.5 min and 42 GB RSS at 16 threads; 2 genomes took 26 s and 15 GB with the database already in page cache. Batch genomes into one call.

### Getting the code

- **On Andrew's server:** use the worktree `~/Documents/KBUtilLib_skani_gtdb` (branch `skani-gtdb-database`).
  - Do **not** use `~/Documents/KBUtilLib`. It is on an old pre-reorg `main` that doesn't have this code, and other notebooks and sessions import it. Update it only with Andrew's OK.
  - Its `git stash` holds an earlier, pre-reorg version of this same change. Don't pop it onto current upstream.
- **On another machine:** run `git clone --branch skani-gtdb-database https://github.com/freiburgermsu/KBUtilLib.git`.
  - Once PR #105 is merged, use `https://github.com/cshenry/KBUtilLib` `main` instead.
- **Python environment:** KBUtilLib needs its own dependencies (cobra, pandas, PyYAML, requests, `pydantic>=2`, …) in the analysis venv, never in `gtdbtk_venv`.
  - On Andrew's server, `py_venv` already has them, and KBUtilLib is used by putting its `src/` on `sys.path`, not by installing it.
  - On a new machine, `uv pip install -e <checkout>` into the analysis venv is the normal route (*not tested here*).

### Usage (run as written on 2026-10-30, from outside the repo, with `~/Documents/py_venv/bin/python`)

```python
import os, sys
sys.path.insert(0, os.path.expanduser("~/Documents/KBUtilLib_skani_gtdb/src"))   # PR #105 branch checkout
from kbutillib.domains.genome.skani_utils import SKANIUtils

u = SKANIUtils(skani_executable=os.path.expanduser("~/Documents/gtdbtk_venv/bin/skani"))
u.register_gtdb_database("~/Documents/gtdbtk_data/release232")   # once per package; no-op if already registered

res = u.classify_genomes_gtdb(["/path/to/bin1.fna", "/path/to/bin2.fna"], "gtdb_r232", threads=16)
for genome, r in res.items():
    print(genome, r["status"], r["reference"], r["ani"], r["af"], r["radius"], r["classification"])
```

Expected, using two wetland MAGs:
- `bin.001.fasta__SRR10127268`: `species assigned GCA_036392515.1 99.72 0.9835 95.0 d__Bacteria;…;s__JABFXC01 sp036392515`, identical to GTDB-Tk.
- `bin.001.fasta__SRR10156750`: `no hit`. GTDB-Tk placed it in the tree as `g__DATHCB01`.

`SKANIUtils` inherits KBUtilLib's config and token handling: it reads `~/.kbutillib/config.yaml`, or else a `config.yaml` in the working directory or its parents. To make it hermetic, as the tests do, pass `config_file=False, token_file=None, kbase_token_file=None`; no KBase token is needed. Don't create `~/.kbutillib/config.yaml` just to set the skani path. KBUtilLib loads a single config file, so that file would take precedence over, not merge with, any other `config.yaml`. Pass `skani_executable` instead.

### Tests and follow-up changes

```bash
cd ~/Documents/KBUtilLib_skani_gtdb
~/Documents/py_venv/bin/python -m pytest tests/domains/test_genome_skani_gtdb.py tests/domains/test_genome_skani_hardening.py tests/domains/test_genome_skani_build_sketch_database.py -q
~/Documents/py_venv/bin/ruff check src/kbutillib/domains/genome/skani_utils.py tests/domains/test_genome_skani_gtdb.py
```

**Check:** 48 passed (13 of them new and offline: skani is mocked, and a fake package is built in `tmp_path`), and ruff prints `All checks passed!`.

- **Pushing:** Andrew's GitHub account (`freiburgermsu`) can **read but not push** to `cshenry/KBUtilLib`; a push to `origin` returns 403.
  - Commit follow-ups on `skani-gtdb-database` and run `git push fork skani-gtdb-database`. PR #105 updates by itself.
  - Commit and push only when Andrew asks.
- **Repo conventions:** commit messages are `feat(scope): …` style, and the repo's own rules are in `.claude/CLAUDE.md` there. Stage files by name, and never use `--no-verify`.
- **Known gap:** the PR doesn't include a fallback for when `skani` isn't on PATH. One possible follow-up: `register_gtdb_database()` records the skani binary path in the registry entry, and `SKANIUtils` falls back to it. The fallback has not been implemented; it was only offered to Andrew.

## 7. Related tooling from the wetland run (`~/gtdbtk_kbase`, not in any repo)

| File | Purpose |
|---|---|
| `run_gtdbtk_all.sh` | The 2,449-MAG `classify_wf` launcher (batchfile from `all_manifest.tsv`) |
| `mem_guard.sh` | Watches a GTDB-Tk systemd scope: logs memory every 30 s, alerts, and freezes only its own cgroup if `MemAvailable` drops below 12 GB |
| `kbase_mags.py` (+ `kbase_mags_kbaseus.py`) | Lists, fetches (md5-checked) and collects for KBase narrative 272312. The wrapper rewrites blob downloads to `kbase.us`, because `services.kbase.us` is unreachable from Andrew's server. |
| `writeback_r232.py` | Writes R232 lineages back to KBase MAG objects, re-points the `MAG_pool` sets, and saves a KBase report. **Gotcha:** for Assembly objects, KBase extracts the `GTDB_*` metadata from `data.std_lineages.gtdb`; user metadata with the same keys is silently overwritten. |
| `post_retry_watch.sh` | Writes the 200 held MAGs (MAG_pool part07/part19) once KBase's 5th GTDB-Tk attempt on those chunks is final. Log: `post_retry_watch.log`, deadline 2026-11-03. |
