> **Status (2026-11-12): historical — this R214 path was not used.** The wetland MAGs were classified with
> GTDB-Tk 2.8.0 / GTDB R232 instead (all 12,611 MAGs; the 248 from the DAS Tool runs repeated on 2026-11-03 follow in run 3, started 2026-11-12), and written to
> `data.std_lineages.gtdb`, not to user metadata. The current hand-off is `gtdbtk_r232_and_kbutillib_skani.md`.
> Do not run `kbase_mags.py writeback`: it labels results as GTDB-Tk 2.3.2 / R214 (`GTDB_source_ver` 214.1) and
> writes user metadata, which KBase overwrites from `std_lineages`.

# Local GTDB-Tk run for the KBase MAGs that KBase couldn't classify

**For:** Claude running on Andrew's Linux server, with shell access.
**Context:** written by Claude (Cowork), which runs Andrew's KBase workflow in narrative 272312.
**Goal:** give GTDB taxonomy to the MAGs whose KBase GTDB-Tk jobs keep failing. Use a dedicated, tested GTDB-Tk environment that matches KBase's run, then hand the results back.

Work through the sections in order. Each section ends with a check; if a check fails, stop and report.

---

## 1. Background

- KBase narrative **272312** ("wetland bins", workspace `freiburgermsu:narrative_1790959438780`, owner `freiburgermsu`) holds **2,449 MAGs**. They are medium-quality or better (≥ 50% complete, < 10% contaminated) and stored as KBase Assembly objects.
  - MAG objects are named `<bin id>__<SRR accession>`, for example `bin.001.fasta__SRR10127268`.
  - The MAGs are split into 25 AssemblySets of 100 MAGs each: `MAG_pool.part01.AssemblySet` to `MAG_pool.part25.AssemblySet`.
- KBase ran **GTDB-Tk 2.3.2** with the **GTDB R214** reference package on each chunk. It used the defaults: `--min_perc_aa 10`, the split (non-full) tree, and the ANI screen.
- When KBase's run succeeds, it saves a new version of each MAG object and adds three metadata fields:
  - `GTDB_lineage`: the full `d__…;s__` string
  - `GTDB_taxon_id`: the lowest named rank
  - `GTDB_source_ver`: `214.1`

  As of Oct 29, 2026, 2,049 MAGs carry these fields.
- Four chunks keep failing: `part04`, `part07`, `part16` and `part19`, 400 MAGs in total.
  - GTDB-Tk itself finishes. KBase's job runner then crashes while reading GTDB-Tk's progress-bar output (`'utf-8' codec can't decode byte 0xe2 in position 16414`), and KBase reports this as "No output generated".
  - `part07` and `part19` crashed at the same byte in 4 of 4 attempts, so more KBase retries are unlikely to help. KBase allows 5 attempts per chunk.
- Andrew's workflow selects MAGs by reading `GTDB_lineage` from each MAG's metadata. The optional last step writes your results into that same field, so his workflow continues unchanged.

## 2. Ground rules

- Keep your files in two places: a work folder, `~/gtdbtk_kbase`, and the reference-data folder chosen in section 5. Don't change other conda environments, system packages or files outside those two folders. Don't use `sudo` without Andrew's OK.
- Treat KBase as **read-only**. The only exception is the final `writeback --apply` in section 9, and that needs Andrew's explicit approval in your conversation.
- The KBase token is a password.
  - Never print it, echo it, log it or put it on a command line.
  - Never copy it into another file, or into a chat.
  - Don't use `set -x` in any shell that touches it.
  - The helper script reads it from `~/.kbase_token` by itself.
- If a check fails, stop and report the exact error and what you tried. Don't work around a failed test.

## 3. Preflight

```bash
mkdir -p ~/gtdbtk_kbase && cd ~/gtdbtk_kbase
nproc; free -g; df -h ~ /data 2>/dev/null
command -v mamba || command -v conda || echo "no conda"
curl -sI https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz | head -1
curl -s -o /dev/null -w '%{http_code}\n' https://kbase.us/services/ws      # any HTTP code means reachable
stat -c '%a %s' ~/.kbase_token                                              # expect: 600 and a non-zero size
```

| Resource | Needed |
|---|---|
| RAM | **≥ 64 GB free.** pplacer needs about 65 GB for the bacterial tree at default settings. With less, stop and ask Andrew; `--scratch_dir` cuts memory but is much slower. |
| Disk | **About 170 GB** on the data volume: a 77.6 GB archive plus about 85 GB unpacked. About 90 GB is enough if you stream-extract (section 5). Add about 10 GB for the Mash sketch and the work folder. |
| CPUs | Any number; more is faster. GTDB's guide is about 1 hour per 1,000 genomes on 64 CPUs. |
| Network | HTTPS to `data.gtdb.ecogenomic.org` and `kbase.us` |

**If there is no conda**, install Miniforge in your home folder; no root is needed:
```bash
curl -L -o /tmp/miniforge.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash /tmp/miniforge.sh -b -p ~/miniforge3 && source ~/miniforge3/bin/activate
```

**If `~/.kbase_token` is missing,** ask Andrew to create it on the server himself:
1. Log in at narrative.kbase.us in a browser.
2. Copy the value of the `kbase_session` cookie (DevTools, Application tab, Cookies, `https://narrative.kbase.us`).
3. On the server, write it as the only line of `~/.kbase_token`, then run `chmod 600 ~/.kbase_token`.

He should not paste the token into a chat. It is his login token, so it stops working if he logs out of KBase. If KBase later answers with a token or authorization error, ask him to refresh the file.

**Check:** RAM, disk and network meet the table; the token file exists with mode 600.

## 4. Build the GTDB-Tk 2.3.2 environment

Pin the exact version KBase used, in its own environment:

```bash
mamba create -y -n gtdbtk-2.3.2 -c conda-forge -c bioconda gtdbtk=2.3.2     # or: conda create ... (same arguments)
eval "$(conda shell.bash hook)" && conda activate gtdbtk-2.3.2              # the hook is needed in non-interactive shells
gtdbtk --version            # must report 2.3.2
for t in prodigal hmmsearch pplacer guppy fastANI mash; do printf '%-9s ' $t; command -v $t || echo MISSING; done
mkdir -p ~/gtdbtk_kbase/env
conda env export -n gtdbtk-2.3.2 > ~/gtdbtk_kbase/env/gtdbtk-2.3.2.env.yml
conda list -n gtdbtk-2.3.2 --explicit > ~/gtdbtk_kbase/env/gtdbtk-2.3.2.explicit.txt   # exact rebuild list
```

Always give the version: an unpinned `gtdbtk` can resolve to an ancient release (1.0.2). Keep `conda-forge` before `bioconda`.

**Check:** `gtdbtk --version` says 2.3.2, and every tool in the loop is found.

## 5. Reference data: GTDB R214 package

GTDB-Tk 2.3.2 supports R214, the release KBase used; it does not support R220 or later. Pick a data folder on a volume with enough space, here called `$DB`, for example `/data/gtdbtk`.

```bash
export DB=/data/gtdbtk && mkdir -p "$DB" && cd "$DB"
URL=https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz
# Option A, if disk allows (~170 GB): download, verify, extract
wget -c "$URL"
echo "630745840850c532546996b22da14c27  gtdbtk_r214_data.tar.gz" | md5sum -c -     # must say OK
tar -xzf gtdbtk_r214_data.tar.gz
# Option B, if disk is tight (~90 GB): stream-extract and check the md5 on the way through
# wget -qO- "$URL" | tee >(md5sum > r214.md5) | tar -xz ; cat r214.md5   # must be 630745840850c532546996b22da14c27
ls "$DB"        # note the unpacked folder name (expected: release214); use it below if different
mkdir -p "$DB/mash"
conda env config vars set -n gtdbtk-2.3.2 GTDBTK_DATA_PATH="$DB/release214"
# One file that sets everything up; every later section starts with: source ~/gtdbtk_kbase/env.sh
cat > ~/gtdbtk_kbase/env.sh <<EOF
eval "\$(conda shell.bash hook)"
conda activate gtdbtk-2.3.2
export DB="$DB"
export GTDBTK_DATA_PATH="$DB/release214"
EOF
source ~/gtdbtk_kbase/env.sh && echo "$GTDBTK_DATA_PATH" && ls "$GTDBTK_DATA_PATH"
```

Delete the archive only after the md5 check and the tests in section 6 pass, and only if disk is needed.

**Check:** the md5 matches `630745840850c532546996b22da14c27`, and `GTDBTK_DATA_PATH` points at the unpacked folder (the one containing `markers/`, `metadata/`, `fastani/` and so on).

## 6. Test the environment

**6a. Install check and the bundled test.**
```bash
source ~/gtdbtk_kbase/env.sh
gtdbtk check_install
gtdbtk test --out_dir ~/gtdbtk_kbase/selftest --cpus 8 ; echo "exit=$?"
```
The test runs `classify_wf` on bundled genomes and exits non-zero if anything is missing.
**Check:** both commands succeed (`exit=0`).

**6b. Known-answer test on 20 of Andrew's own MAGs that KBase already classified.** This proves your setup reproduces KBase's results before you touch the real MAGs. First save the helper script from the appendix as `~/gtdbtk_kbase/kbase_mags.py` and verify it:

```bash
source ~/gtdbtk_kbase/env.sh && cd ~/gtdbtk_kbase
sha256sum kbase_mags.py       # must be a210e9950ac3b8ea642aad445cf21fe88d272d753fc556261c53183d3a3327d9
python3 kbase_mags.py list --controls 20      # read-only; also writes mags_to_classify.tsv (used in section 7)
python3 kbase_mags.py fetch --list controls_to_check.tsv --genome-dir controls --manifest controls_manifest.tsv
gtdbtk classify_wf --genome_dir controls -x fna --out_dir gtdbtk_controls \
    --cpus "$(nproc)" --pplacer_cpus 2 --min_perc_aa 10 --mash_db "$DB/mash/gtdb_r214.msh"
python3 kbase_mags.py collect --manifest controls_manifest.tsv --out-dir gtdbtk_controls --out controls_lineages.tsv
```

The first `classify_wf` run builds the Mash sketch of the reference genomes at `$DB/mash/gtdb_r214.msh`; this is a one-time step that can take up to about an hour. Later runs reuse it. In GTDB-Tk 2.3.2 `classify_wf` requires either `--mash_db` or `--skip_ani_screen`; use `--mash_db` to match KBase's run.

**Check:** `collect` prints `CONTROLS: x/20 identical …; 20/20 agree on phylum and genus.` Phylum and genus must agree for all 20. If any species-level differences appear, list them in your report. If phylum or genus differs for any control, stop and report.

## 7. Get the MAGs to classify

```bash
source ~/gtdbtk_kbase/env.sh && cd ~/gtdbtk_kbase
python3 kbase_mags.py list
```

`list` writes `mags_to_classify.tsv` with every MAG that has no `GTDB_lineage` yet. It leaves out chunks whose KBase GTDB-Tk job is still queued or running, and prints how many times each listed chunk has failed in KBase.

- Proceed if Andrew asked you to run now, or if every listed chunk shows **5** KBase failures.
- If a listed chunk shows fewer than 5 failures, KBase will retry it on Andrew's next daily run, at 8:59 a.m. Chicago time. Ask Andrew whether to wait or go ahead.
- `list --include-active` also lists chunks that are still running in KBase. Use it only if Andrew says so.

```bash
python3 kbase_mags.py fetch          # downloads into genomes/, checks each file's md5, writes genomes_manifest.tsv
ls genomes | wc -l                   # should equal the number of rows in mags_to_classify.tsv
```

`fetch` resumes where it left off if interrupted.
**Check:** the file count matches and `fetch` reported no errors.

## 8. Run GTDB-Tk

```bash
source ~/gtdbtk_kbase/env.sh && cd ~/gtdbtk_kbase
nohup gtdbtk classify_wf --genome_dir genomes -x fna --out_dir gtdbtk_out \
    --cpus "$(nproc)" --pplacer_cpus 2 --min_perc_aa 10 --mash_db "$DB/mash/gtdb_r214.msh" \
    > gtdbtk_run.log 2>&1 &
echo $! > gtdbtk.pid
# Poll every few minutes (don't use tail -f): 400 MAGs should take well under an hour once the Mash sketch exists
kill -0 "$(cat gtdbtk.pid)" 2>/dev/null && echo running || echo finished; tail -n 5 gtdbtk_run.log
```

- **Memory:** pplacer uses more memory with more threads. With about 64–96 GB of RAM keep `--pplacer_cpus 1–2`; with 128 GB or more you can raise it. If pplacer is killed for lack of memory, rerun with `--pplacer_cpus 1` and add `--scratch_dir ~/gtdbtk_kbase/pplacer_scratch`.
- **Settings:** don't add `--full_tree` or change `--min_perc_aa`; both would make the results differ from KBase's.

**Check:** the log ends with `Done.`, and `gtdbtk_out/` contains `gtdbtk.bac120.summary.tsv` (and `gtdbtk.ar53.summary.tsv` if any archaea are present), either at the top level or under `classify/`.

## 9. Collect results, report, and write back to KBase (only with approval)

```bash
source ~/gtdbtk_kbase/env.sh && cd ~/gtdbtk_kbase
python3 kbase_mags.py collect        # writes mag_gtdb_lineages.tsv and prints domain, phylum and novelty counts
python3 kbase_mags.py writeback      # DRY RUN: shows how many MAGs would get a lineage, with 3 examples
```

Report to Andrew, then wait. If he approves writing to KBase:

```bash
python3 kbase_mags.py writeback --apply
```

`--apply` saves a new version of each classified MAG object. Each version holds the same data plus `GTDB_lineage`, `GTDB_taxon_id`, `GTDB_source_ver` = `214.1` and `GTDB_classified_by` = `local GTDB-Tk 2.3.2, GTDB R214 package`. Each version also gets a provenance note naming this local run. This mirrors what KBase's own GTDB-Tk app writes.

The script never writes:
- the control MAGs
- MAGs that came out "Unclassified"
- MAGs that KBase classified in the meantime

It finishes by re-reading the objects to confirm the lineage is there. After that, Andrew's next run of the "MAG annotation B" cell, or his daily scheduled run, uses these lineages to pick the 5 MAGs.

There is one exception. If the narrative already contains an object named `selected_MAGs.AssemblySet`, the 5 MAGs were picked before your results arrived, and the workflow won't redo the pick on its own. Tell Andrew in that case.

**Check:** the dry run's count equals `collect`'s "with a lineage" count, minus any MAGs KBase classified meanwhile. After `--apply`, the final line reports all checked objects carrying `GTDB_lineage`.

## 10. What to report to Andrew

- Environment: GTDB-Tk version, `GTDBTK_DATA_PATH`, the md5 check result, and the `check_install` and `test` outcomes.
- The control result line (x/20 identical; 20/20 phylum and genus) and any differences.
- The MAGs classified: how many in total; domains; top phyla; novelty (no GTDB genus / genus but no species / species).
- Any MAGs filtered or unclassified, by name, with GTDB-Tk's reason. Look in `gtdbtk_out/align/*filtered.tsv`, `gtdbtk_out/identify/*failed*` and the `note` and `warnings` columns.
- Where the files are: `~/gtdbtk_kbase/mag_gtdb_lineages.tsv`, `gtdbtk_out/`, `env/`.
- Whether you wrote back to KBase, and the verification line.

## Troubleshooting

| Symptom | Fix |
|---|---|
| `You must specify a path to the mash database with --mash_db` | Add `--mash_db "$DB/mash/gtdb_r214.msh"` (or `--skip_ani_screen`, but then mention it in the report). |
| pplacer killed, or `MemoryError` | Use `--pplacer_cpus 1` and `--scratch_dir`; check that nothing else is using RAM. |
| `check_install` reports missing or malformed files | Re-verify the md5 and re-extract; make sure `GTDBTK_DATA_PATH` points at the unpacked folder itself. |
| The helper stops with a KBase error mentioning a token or authorization | Ask Andrew to refresh `~/.kbase_token`. Never ask for the token in chat. |
| `fetch` reports an md5 mismatch or a network error | Rerun `fetch`; it retries and skips files already verified. |
| MAGs missing from the summary | GTDB-Tk filtered them (too few marker genes, below `--min_perc_aa 10`) or couldn't find markers. List them; don't change settings to force a call. |
| Conda solve is very slow or fails | Use `mamba`; keep `conda-forge` listed before `bioconda`; try adding `--strict-channel-priority`. Don't relax the `gtdbtk=2.3.2` pin. |
| `conda activate` fails with "Run 'conda init' before 'conda activate'" | Run `eval "$(conda shell.bash hook)"` first (`env.sh` does this). |

---

## Appendix: `kbase_mags.py`

Save this exactly as `~/gtdbtk_kbase/kbase_mags.py`. Its sha256 should be `a210e9950ac3b8ea642aad445cf21fe88d272d753fc556261c53183d3a3327d9`. It uses only the Python standard library, so it runs with any Python 3.8 or later, including the one in the GTDB-Tk environment.

```python
#!/usr/bin/env python3
"""Move KBase MAGs that KBase's GTDB-Tk could not classify through a local GTDB-Tk run.

Subcommands (run in order):
  list       find MAG Assembly objects in narrative 272312 that have no GTDB lineage yet
             (--controls N also picks N MAGs KBase already classified, as a known-answer test set)
  fetch      download FASTA files (md5-checked) for a list
  collect    read local GTDB-Tk output for a manifest; compares with KBase's lineages for control sets
  writeback  add the lineage to each MAG's KBase object metadata (dry run unless --apply)

Standard library only. The KBase token is read from $KB_AUTH_TOKEN or ~/.kbase_token and is never printed.
"""
import argparse, csv, glob, hashlib, json, os, random, re, socket, sys, time, urllib.error, urllib.request
from collections import Counter, defaultdict

WSID = 272312                                   # narrative 272312 ("wetland bins")
WS = "https://kbase.us/services/ws"
HS = "https://kbase.us/services/handle_service"
EE2 = "https://kbase.us/services/ee2"
GTDBTK_APP = "kb_gtdbtk/run_kb_gtdbtk_classify_wf"
ACTIVE = {"created", "estimating", "queued", "running"}
MAG_NAME = re.compile(r"__SRR\d+$")             # MAG objects are named <bin id>__<run accession>
CHUNK = re.compile(r"^MAG_pool\.part(\d\d)\.AssemblySet$")
LIST_TSV, MANIFEST, LINEAGES = "mags_to_classify.tsv", "genomes_manifest.tsv", "mag_gtdb_lineages.tsv"
CONTROLS_TSV = "controls_to_check.tsv"


def token():
    t = os.environ.get("KB_AUTH_TOKEN", "").strip()
    if not t:
        p = os.path.expanduser("~/.kbase_token")
        if os.path.exists(p):
            t = open(p).read().strip()
    if not t:
        sys.exit("No KBase token: set KB_AUTH_TOKEN or put it in ~/.kbase_token (chmod 600).")
    return t


def rpc(url, method, params, tries=4):
    """KBase JSON-RPC 1.1 call. Service errors stop the script; network/gateway errors are retried."""
    body = json.dumps({"version": "1.1", "id": "1", "method": method, "params": params}).encode()
    for k in range(tries):
        req = urllib.request.Request(url, data=body, headers={"Authorization": token(),
                                                              "Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                res = json.load(r)
        except urllib.error.HTTPError as e:      # KBase reports service errors as JSON with HTTP 500
            try:
                res = json.load(e)
            except Exception:
                res = None                         # gateway/proxy error: retry
        except (urllib.error.URLError, TimeoutError, socket.timeout, ConnectionError):
            res = None                             # network problem: retry
        if res is not None:
            if res.get("error"):
                err = res["error"]
                sys.exit(f"{method} failed: {str(err.get('message') or err)[:400]}")
            return res["result"]
        if k == tries - 1:
            sys.exit(f"{method}: no usable response from {url} after {tries} tries")
        time.sleep(10 * (k + 1))


def upa(i):
    return f"{i[6]}/{i[0]}/{i[4]}"


def list_objects(type_str, meta=False):
    out, minid = [], 1
    while True:
        batch = rpc(WS, "Workspace.list_objects", [{"ids": [WSID], "type": type_str, "minObjectID": minid,
                                                     "includeMetadata": int(meta), "limit": 5000}])[0]
        out += batch
        if len(batch) < 5000:
            return out
        minid = max(i[0] for i in batch) + 1


def write_tsv(path, rows, cols):
    if any("kbase_lineage" in r for r in rows) and "kbase_lineage" not in cols:
        cols = cols + ["kbase_lineage"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ------------------------------------------------------------------ list
def cmd_list(a):
    mags = {i[1]: i for i in list_objects("KBaseGenomeAnnotations.Assembly", meta=True) if MAG_NAME.search(i[1])}
    chunk_of = {}
    for s in list_objects("KBaseSets.AssemblySet"):
        if CHUNK.match(s[1]):
            items = rpc(WS, "Workspace.get_objects2", [{"objects": [{"ref": upa(s)}]}])[0]["data"][0]["data"]["items"]
            for it in items:
                chunk_of[it.get("label") or ""] = s[1]
    jobs = rpc(EE2, "execution_engine2.check_workspace_jobs",
               [{"workspace_id": WSID, "return_list": 1,
                 "exclude_fields": ["job_output", "authstrat", "condor_job_ads", "scheduler_type", "scheduler_id"]}])[0]
    jobs = jobs["job_states"] if isinstance(jobs, dict) else jobs
    by_chunk = defaultdict(Counter)
    for j in jobs:
        ji = j.get("job_input") or {}
        if ji.get("app_id") != GTDBTK_APP:
            continue
        base = (ji.get("params") or [{}])[0].get("output_tree_basename", "")
        st = j.get("status")
        by_chunk[base.replace(".GTDB_tree", ".AssemblySet")]["active" if st in ACTIVE else st] += 1
    rows, skipped = [], Counter()
    for name, i in sorted(mags.items()):
        if (i[10] or {}).get("GTDB_lineage"):
            continue                                    # KBase already classified it
        chunk = chunk_of.get(name, "")
        st = by_chunk.get(chunk, Counter())
        if st["active"] and not a.include_active:
            skipped[chunk] += 1
            continue
        rows.append({"name": name, "ref": upa(i), "chunk": chunk, "kbase_gtdbtk_failures": st["error"],
                     "size_bp": (i[10] or {}).get("Size", ""), "contigs": (i[10] or {}).get("N Contigs", "")})
    write_tsv(LIST_TSV, rows, ["name", "ref", "chunk", "kbase_gtdbtk_failures", "size_bp", "contigs"])
    print(f"{len(mags)} MAG objects in narrative {WSID}; "
          f"{sum(1 for i in mags.values() if (i[10] or {}).get('GTDB_lineage'))} already have a GTDB lineage.")
    print(f"{len(rows)} MAGs to classify locally -> {LIST_TSV}")
    for c, n in sorted(Counter(r["chunk"] for r in rows).items()):
        print(f"   {c or '(not in a MAG_pool chunk)'}: {n} MAGs, KBase GTDB-Tk failures: {by_chunk[c]['error']}")
    for c, n in sorted(skipped.items()):
        print(f"   skipped {c}: {n} MAGs, a KBase GTDB-Tk job for this chunk is still queued/running")
    if a.controls:
        done = sorted(n for n, i in mags.items() if (i[10] or {}).get("GTDB_lineage"))
        pick = sorted(random.Random(214).sample(done, min(a.controls, len(done))))
        write_tsv(CONTROLS_TSV, [{"name": n, "ref": upa(mags[n]), "chunk": chunk_of.get(n, ""),
                                  "kbase_lineage": mags[n][10]["GTDB_lineage"]} for n in pick],
                  ["name", "ref", "chunk", "kbase_lineage"])
        print(f"{len(pick)} already-classified MAGs picked as known-answer controls -> {CONTROLS_TSV}")


# ------------------------------------------------------------------ fetch
def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def cmd_fetch(a):
    rows = read_tsv(a.list)
    os.makedirs(a.genome_dir, exist_ok=True)
    manifest, done = [], 0
    for k in range(0, len(rows), 50):
        part = rows[k:k + 50]
        objs = rpc(WS, "Workspace.get_objects2", [{"objects": [{"ref": r["ref"], "included": ["/fasta_handle_ref"]}
                                                               for r in part]}])[0]["data"]
        hids = [o["data"]["fasta_handle_ref"] for o in objs]
        handles = rpc(HS, "AbstractHandle.hids_to_handles", [hids])[0]
        hmap = {h["hid"]: h for h in handles}
        for r, hid in zip(part, hids):
            h = hmap[hid]
            path = os.path.join(a.genome_dir, r["name"] + ".fna")
            if not (os.path.exists(path) and md5sum(path) == h.get("remote_md5")):
                url = f'{h["url"]}/node/{h["id"]}?download_raw'
                tmp = path + ".part"
                for attempt in range(4):
                    try:
                        req = urllib.request.Request(url, headers={"Authorization": "OAuth " + token()})
                        with urllib.request.urlopen(req, timeout=600) as resp, open(tmp, "wb") as out:
                            while True:
                                block = resp.read(1 << 20)
                                if not block:
                                    break
                                out.write(block)
                        if h.get("remote_md5") and md5sum(tmp) != h["remote_md5"]:
                            raise IOError("md5 mismatch")
                        os.replace(tmp, path)
                        break
                    except Exception as e:
                        if attempt == 3:
                            sys.exit(f"download failed for {r['name']}: {e}")
                        time.sleep(10 * (attempt + 1))
            with open(path) as fh:
                if not fh.read(1) == ">":
                    sys.exit(f"{path} is not FASTA")
            manifest.append(dict(r, file=path, md5=md5sum(path)))
            done += 1
        print(f"  {done}/{len(rows)} MAGs downloaded and md5-checked")
    write_tsv(a.manifest, manifest, ["name", "ref", "chunk", "file", "md5", "size_bp", "contigs"])
    print(f"Wrote {len(manifest)} FASTA files to {a.genome_dir}/ and {a.manifest}")


# ------------------------------------------------------------------ collect
def first_real(r, *keys):
    return next((r[k] for k in keys if r.get(k) not in (None, "", "N/A")), "")


def cmd_collect(a):
    want = {r["name"]: r for r in read_tsv(a.manifest)}
    got, files = {}, sorted(set(os.path.realpath(p) for p in
                                glob.glob(os.path.join(a.out_dir, "**", "gtdbtk.*.summary.tsv"), recursive=True)))
    files = [f for f in files if re.search(r"gtdbtk\.(bac120|ar53|ar122)\.summary\.tsv$", f)]
    if not files:
        sys.exit(f"No gtdbtk.*.summary.tsv under {a.out_dir}; did classify_wf finish?")
    for f in files:
        for r in read_tsv(f):
            g = r.get("user_genome")
            if g in want and r.get("classification"):
                got[g] = r
    rows, missing = [], []
    for name, m in sorted(want.items()):
        r = got.get(name)
        if not r:
            missing.append(name)
            continue
        row = {"name": name, "ref": m["ref"], "chunk": m["chunk"], "classification": r["classification"],
               "classification_method": r.get("classification_method", ""),
               "closest_reference": first_real(r, "fastani_reference", "closest_placement_reference"),
               "ani": first_real(r, "fastani_ani", "closest_placement_ani"),
               "red_value": r.get("red_value", ""), "warnings": r.get("warnings", ""), "note": r.get("note", "")}
        if m.get("kbase_lineage"):
            row["kbase_lineage"] = m["kbase_lineage"]
        rows.append(row)
    write_tsv(a.out, rows, ["name", "ref", "chunk", "classification", "classification_method", "closest_reference",
                            "ani", "red_value", "warnings", "note"])
    ok = [r for r in rows if not r["classification"].startswith("Unclassified")]
    ranks = lambda c: dict(x.split("__", 1) for x in c.split(";") if "__" in x)
    print(f"Read {len(files)} summary file(s). {len(rows)}/{len(want)} MAGs in the output, {len(ok)} with a lineage.")
    print("Domains:", dict(Counter(ranks(r["classification"]).get("d", "") for r in ok)))
    print("Phyla:", dict(Counter(ranks(r["classification"]).get("p", "") for r in ok).most_common()))
    nov = Counter("no genus" if not ranks(r["classification"]).get("g") else
                  "no species" if not ranks(r["classification"]).get("s") else "species" for r in ok)
    print("Novelty:", dict(nov))
    if missing:
        print(f"{len(missing)} MAGs are not in any summary (filtered for too few markers, or failed). "
              f"See {a.out_dir}/align/*filtered.tsv and {a.out_dir}/identify/*failed*: {', '.join(missing[:10])}"
              + (" ..." if len(missing) > 10 else ""))
    ctrl = [r for r in rows if r.get("kbase_lineage")]
    if ctrl:                                  # known-answer controls: compare with KBase's GTDB-Tk lineages
        exact = [r for r in ctrl if r["classification"] == r["kbase_lineage"]]
        genus = [r for r in ctrl if ranks(r["classification"]).get("g") == ranks(r["kbase_lineage"]).get("g")
                 and ranks(r["classification"]).get("p") == ranks(r["kbase_lineage"]).get("p")]
        print(f"CONTROLS: {len(exact)}/{len(ctrl)} identical to KBase's lineage; {len(genus)}/{len(ctrl)} agree "
              "on phylum and genus.")
        for r in ctrl:
            if r["classification"] != r["kbase_lineage"]:
                print(f"   differs: {r['name']}\n      local: {r['classification']}\n      KBase: {r['kbase_lineage']}")
    print(f"Wrote {a.out}")


# ------------------------------------------------------------------ writeback
def taxon_id(classification):
    parts = [x for x in classification.split(";") if "__" in x and x.split("__", 1)[1]]
    return parts[-1] if parts else ""


def cmd_writeback(a):
    rows = [r for r in read_tsv(a.lineages)
            if r["classification"] and not r["classification"].startswith("Unclassified") and not r.get("kbase_lineage")]
    host, now = socket.gethostname(), time.strftime("%Y-%m-%dT%H:%M:%S+0000", time.gmtime())
    plan, skip = [], []
    for k in range(0, len(rows), 25):
        part = rows[k:k + 25]
        names = [{"wsid": WSID, "name": r["name"]} for r in part]    # latest version of each MAG
        infos = rpc(WS, "Workspace.get_object_info3", [{"objects": names, "includeMetadata": 1}])[0]["infos"]
        for r, info in zip(part, infos):
            if (info[10] or {}).get("GTDB_lineage"):
                skip.append(r["name"])                           # KBase classified it in the meantime
            else:
                plan.append((r, info))
    print(f"{len(plan)} MAGs would get a GTDB lineage; {len(skip)} skipped (already have one in KBase).")
    for r, info in plan[:3]:
        print(f"   {r['name']} ({upa(info)}): {r['classification']}")
    if not a.apply:
        print("Dry run only. Re-run with --apply after Andrew approves writing to KBase.")
        return
    saved = 0
    for k in range(0, len(plan), 25):
        part = plan[k:k + 25]
        data = rpc(WS, "Workspace.get_objects2", [{"objects": [{"ref": upa(i)} for _, i in part]}])[0]["data"]
        objs = []
        for (r, info), d in zip(part, data):
            meta = dict(info[10] or {})
            meta.update({"GTDB_lineage": r["classification"], "GTDB_taxon_id": taxon_id(r["classification"]),
                         "GTDB_source_ver": "214.1", "GTDB_classified_by": "local GTDB-Tk 2.3.2, GTDB R214 package"})
            objs.append({"type": info[2], "objid": info[0], "data": d["data"], "meta": meta,
                         "provenance": [{"time": now, "script": "kbase_mags.py writeback", "script_ver": "1",
                                         "input_ws_objects": [upa(info)],
                                         "description": f"GTDB-Tk 2.3.2 classify_wf (GTDB R214) run locally on {host}; "
                                                        "lineage added to object metadata"}]})
        rpc(WS, "Workspace.save_objects", [{"id": WSID, "objects": objs}])
        saved += len(objs)
        print(f"  saved {saved}/{len(plan)}")
    check = rpc(WS, "Workspace.get_object_info3",
                [{"objects": [{"wsid": WSID, "name": r["name"]} for r, _ in plan[:50]], "includeMetadata": 1}])[0]["infos"]
    print(f"Verified: {sum(1 for i in check if (i[10] or {}).get('GTDB_lineage'))}/{len(check)} checked objects "
          "now carry GTDB_lineage.")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("list")
    s.add_argument("--include-active", action="store_true",
                   help="also list MAGs whose chunk still has a KBase GTDB-Tk job queued or running")
    s.add_argument("--controls", type=int, default=0, metavar="N",
                   help=f"also write N already-classified MAGs to {CONTROLS_TSV} as a known-answer test")
    s = sub.add_parser("fetch")
    s.add_argument("--list", default=LIST_TSV)
    s.add_argument("--genome-dir", default="genomes")
    s.add_argument("--manifest", default=MANIFEST)
    s = sub.add_parser("collect")
    s.add_argument("--manifest", default=MANIFEST)
    s.add_argument("--out-dir", default="gtdbtk_out", help="GTDB-Tk classify_wf --out_dir")
    s.add_argument("--out", default=LINEAGES)
    s = sub.add_parser("writeback")
    s.add_argument("--lineages", default=LINEAGES)
    s.add_argument("--apply", action="store_true", help="actually save new object versions in KBase")
    a = p.parse_args()
    {"list": cmd_list, "fetch": cmd_fetch, "collect": cmd_collect, "writeback": cmd_writeback}[a.cmd](a)


if __name__ == "__main__":
    main()
```
