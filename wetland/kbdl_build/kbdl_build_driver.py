#!/usr/bin/env python3
"""Build every wetland MAG in KBDL with Bakta (handoff Phase 2; per Andrew 2026-11-15: gene_caller=bakta, all MAGs).

Runs ON POPLAR against the KBDL job service at http://localhost:8791 (Chris Henry's service). Inputs:
  archive_plan.tsv  517 archives (canary 2 MAGs, pilot 23, giant 1, bulk 25 each) with each MAG's predicted KBDL content
                    hash (= fasta_object_id); MAG FASTAs in MAGS_DIR/<MAG>.fna.
Per job: KBDLBuildGenome {archive, gene_caller='bakta', call_non_coding=true, skani_db=--skani-db, bakta_db=6.0-full}.
--skani-db must be the GTDB R232 SKANIDB object id (per Andrew 2026-11-15: R232, not r220); the canary and pilot ran on
r220 (ref:cee67d7f…) before that decision. Each job's result is checked against the skani db in its own payload.
These are permanent: KBDL dedups builds by content hash across owners, so the first build fixes caller and taxonomy.

v2 (2026-11-15) applies the 22 confirmed findings of an adversarial review (wf_8d5f3d4a-1f3). Safety:
  * one driver per ledger (flock); every payload is asserted right before POST /jobs;
  * POST /jobs is never blindly retried: an intent row is written first; any unresolved intent (lost response, crash,
    5xx) parks its archive as 'lost' and is reconciled against GET /jobs by submit time (adopt / never_created after
    a grace period); nothing new is submitted while an intent is unresolved;
  * failed jobs are always marked harvested; worker_lost/internal_error/store_contention retry the archive (max 3
    attempts), anything else marks the archive 'failed' and halts;
  * skani-failure members are queued as a subset retry and submitted by the main loop under the in-flight cap;
  * uploads: deterministic .tar.xz, GET /objects/{sha} first; 202 is not "stored": upload job result == local sha256,
    GET /objects/{id} (GenomeArchive, private), /files names+sizes == members; max 2 re-uploads, then halt;
  * per member: in genomes xor errors; fasta_object_id == predicted hash; gene_caller == bakta; non_coding_source ==
    called; Genome JSON: genome_id == hash[:16], contigs re-rendered hash == predicted hash, >= 1 CDS, contigs/bp,
    non-coding count, provenance (09_expected_genome_provenance.json incl. cds_overlap_rule/fraction), cds_discarded;
    a temporarily unreadable Genome JSON is re-checked later, not failed;
  * Bakta edge-interval failures (KBDL's strict linear-interval check rejects a Bakta feature ending past a contig end,
    a deterministic per-MAG outcome; first seen 2026-11-17 on bin.001.fasta__SRR11599182) are recorded as member errors
    and listed in kv bakta_edge_mags; halt only when > BAKTA_EDGE_PER_JOB in one job or more than BAKTA_EDGE_RATE of the
    harvested MAGs (floor BAKTA_EDGE_FLOOR); observed ~1% (2 of ~250 by 2026-11-17 21:11);
  * circuit breaker (stop submitting/uploading, keep polling + harvesting) on any other 'Bakta:' error, caller/hash/content/
    provenance mismatch, failed archive, 4xx on submit, or a per-job quality breaker (> half of the >= 1 Mb members
    without tRNA, or > 2 members below 500 CDS/Mb);
  * watchdog ALERTs only (never cancel: KBDL keeps building after a cancel): running > 3x predicted + 30 min, queued
    > 12 h, no harvest for 6 h while jobs are in flight;
  * phases canary -> pilot -> bulk; a phase with failed/partial archives blocks the next one explicitly (STOPPED with
    a GATE message) unless the operator sets kv accept_phase:<phase>=1.
Operator: after a poplar reboot or crash just restart the same command (the ledger resumes; running builds that the
reboot killed come back as worker_lost and are retried). To clear a halt: fix the cause, DELETE FROM kv WHERE k='halted'.
Ledger: SQLite on local disk (WORK, never NFS). Token: read from stdin, then the process detaches; never written to disk.
Usage (poplar): python3 -I kbdl_build_driver.py --through canary|pilot|bulk [--max-in-flight 3] [--dry-run] < token
"""
import argparse, calendar, fcntl, hashlib, io, json, lzma, os, sqlite3, sys, tarfile, time
import requests

BASE = "http://localhost:8791"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = "/scratch/afreiburger/wetland_kbdl_build"
MAGS_DIR = "/scratch1/afreiburger/wetland_coverm_pilot/all_mags"
GTDB_R220 = "ref:cee67d7f33f57c07d569f9fb0a0537a0263094f380fa0b10faa56aafdfcc56ae"
BAKTA_DB = "ref:2fdfb17daa745a4fe41065a38a9a1d367469fe85fcb7a1b1023f8142d586f577"
PHASES = ["canary", "pilot", "bulk"]
ALLOWED_KEYS = {"archive", "gene_caller", "call_non_coding", "skani_db", "bakta_db"}
TERMINAL = {"completed", "failed"}
RETRYABLE = {"worker_lost", "internal_error", "store_contention"}
MAX_ATTEMPTS = 3
RECONCILE_GRACE = 900          # s before an intent with no visible job is declared never_created
BAKTA_EDGE_PER_JOB = 2         # Bakta edge-interval member failures tolerated per job ...
BAKTA_EDGE_RATE = 0.03         # ... and in all as a share of harvested MAGs (observed ~1%); above it: systematic, halt
BAKTA_EDGE_FLOOR = 10          # never halt on the total below this many
RECHECK_EVERY = 1800           # s between re-checks of members whose Genome JSON was temporarily unreadable

p = argparse.ArgumentParser()
p.add_argument("--through", choices=PHASES, default="canary")
p.add_argument("--max-in-flight", type=int, default=3)
p.add_argument("--prefetch", type=int, default=6, help="archives uploaded+verified ahead of submission")
p.add_argument("--dry-run", action="store_true", help="build archives + hashes, print payloads, no HTTP writes")
p.add_argument("--foreground", action="store_true")
p.add_argument("--skani-db", default=None, help="SKANIDB object id for new submissions (the GTDB R232 ref)")
a = p.parse_args()

TOKEN = sys.stdin.read().strip()
assert len(TOKEN) >= 20, "token expected on stdin"
os.makedirs(os.path.join(WORK, "archives"), exist_ok=True)
os.makedirs(os.path.join(WORK, "results"), exist_ok=True)
if not a.foreground and not a.dry_run:
    if os.fork():
        sys.exit(0)
    os.setsid()
    log_fh = open(os.path.join(WORK, "driver.log"), "a", buffering=1)
    sys.stdout = sys.stderr = log_fh
if not a.dry_run:                                     # single instance per ledger
    LOCK = open(os.path.join(WORK, "driver.lock"), "a+")
    try:
        fcntl.flock(LOCK, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        LOCK.seek(0)
        print(f"another driver holds {WORK}/driver.lock (pid {LOCK.read().strip()}); exiting", flush=True)
        sys.exit(1)
    LOCK.seek(0)
    LOCK.truncate()
    LOCK.write(str(os.getpid()))
    LOCK.flush()
S = requests.Session()
S.headers["Authorization"] = "Bearer " + TOKEN
EXPECT_PROV = {k: v for k, v in json.load(open(os.path.join(HERE, "09_expected_genome_provenance.json"))).items()
               if not k.startswith("_")}
assert EXPECT_PROV.get("cds_overlap_rule") == "applied" and EXPECT_PROV.get("cds_overlap_drop_fraction") == 0.5


def ev(msg):
    line = f"[{time.strftime('%F %T')}] {msg}"
    print(line, flush=True)
    with open(os.path.join(WORK, "events.log"), "a") as fh:
        fh.write(line + "\n")


# ---------------- ledger ----------------
db = sqlite3.connect(os.path.join(WORK, "ledger.sqlite"), timeout=60)
db.execute("PRAGMA journal_mode=WAL")
db.executescript("""
CREATE TABLE IF NOT EXISTS archives(name TEXT PRIMARY KEY, phase TEXT, n INT, pred_s REAL, path TEXT, sha256 TEXT,
  upload_job TEXT, object_id TEXT, state TEXT DEFAULT 'planned', updated TEXT);
CREATE TABLE IF NOT EXISTS members(mag TEXT PRIMARY KEY, archive TEXT, predicted_hash TEXT, bp INT, n_contigs INT,
  job_id TEXT, genome_id TEXT, genome_object_id TEXT, fasta_object_id TEXT, prefix TEXT, gene_caller TEXT,
  non_coding_source TEXT, counts_json TEXT, taxonomy_json TEXT, cds INT, n_noncoding INT, error_code TEXT,
  error_msg TEXT, flags TEXT, verified INT DEFAULT 0, updated TEXT);
CREATE TABLE IF NOT EXISTS intents(id INTEGER PRIMARY KEY AUTOINCREMENT, archive TEXT, payload TEXT, created REAL,
  posted REAL, job_id TEXT, outcome TEXT);
CREATE TABLE IF NOT EXISTS jobs(job_id TEXT PRIMARY KEY, archive TEXT, genome_ids TEXT, submitted REAL, state TEXT,
  started REAL, finished REAL, failure TEXT, harvested INT DEFAULT 0, alerted INT DEFAULT 0, attempt INT DEFAULT 1);
CREATE TABLE IF NOT EXISTS kv(k TEXT PRIMARY KEY, v TEXT);
""")


def now():
    return time.strftime("%F %T")


def kv_get(k, d=None):
    r = db.execute("SELECT v FROM kv WHERE k=?", (k,)).fetchone()
    return r[0] if r else d


def kv_set(k, v):
    db.execute("INSERT OR REPLACE INTO kv VALUES(?,?)", (k, str(v)))
    db.commit()


def kv_del(k):
    db.execute("DELETE FROM kv WHERE k=?", (k,))
    db.commit()


def halted():
    return kv_get("halted") == "1"


def halt(msg):
    if not halted():
        ev(f"HALT {msg}")
    kv_set("halted", "1")
    kv_set("halt_reason", msg)
    return None


plan = [l.rstrip("\n").split("\t") for l in open(os.path.join(HERE, "archive_plan.tsv"))][1:]
if not db.execute("SELECT count(*) FROM archives").fetchone()[0]:
    order = {}
    for name, phase, mag, member, h, bp, nc, pred in plan:
        ph = phase if phase in ("canary", "pilot") else "bulk"
        order.setdefault(name, [ph, 0, 0.0])
        order[name][1] += 1
        order[name][2] += float(pred)
        db.execute("INSERT INTO members(mag, archive, predicted_hash, bp, n_contigs) VALUES(?,?,?,?,?)",
                   (mag, name, h, int(bp), int(nc)))
    for name, (ph, n, ps) in order.items():
        db.execute("INSERT INTO archives(name, phase, n, pred_s, state, updated) VALUES(?,?,?,?,'planned',?)",
                   (name, ph, n, ps, now()))
    db.commit()
n_m = db.execute("SELECT count(*), count(DISTINCT mag) FROM members").fetchone()
assert n_m == (12859, 12859), n_m
SKANI_DB = a.skani_db or kv_get("skani_db")
if SKANI_DB and kv_get("skani_db") and kv_get("skani_db") != SKANI_DB:
    sys.exit(f"--skani-db {SKANI_DB} differs from the one this run was started with ({kv_get('skani_db')}); refusing")
if SKANI_DB == GTDB_R220:
    sys.exit("refusing --skani-db r220: new builds must use the GTDB R232 SKANIDB (per Andrew 2026-11-15)")
if SKANI_DB and not a.dry_run:
    kv_set("skani_db", SKANI_DB)


def members_of(name):
    return [r[0] for r in db.execute("SELECT mag FROM members WHERE archive=? ORDER BY mag", (name,))]


# ---------------- HTTP ----------------
def http(method, path, retry=None, **kw):
    """retry: GETs (and read-only POST /jobs/status) retry on connection errors and 5xx; other POSTs never do."""
    retry = (method == "GET") if retry is None else retry
    timeout = kw.pop("timeout", 120)
    for i in range(8):
        try:
            r = S.request(method, BASE + path, timeout=timeout, **kw)
        except requests.RequestException as e:
            if not retry:
                raise
            ev(f"WARN {method} {path}: {e!s:.150}; retry")
            time.sleep(min(300, 15 * 2 ** i))
            continue
        if r.status_code == 401:
            if not retry:
                return r
            ev("ALERT 401 from KBDL (token?); pausing 30 min")
            time.sleep(1800)
            continue
        if r.status_code in (500, 502, 503, 504) and retry:
            time.sleep(min(300, 15 * 2 ** i))
            continue
        return r
    raise RuntimeError(f"{method} {path}: giving up")


# ---------------- archives / uploads ----------------
def build_archive(name):
    path = os.path.join(WORK, "archives", name + ".tar.xz")
    mems = members_of(name)
    raw = io.BytesIO()
    with tarfile.open(fileobj=raw, mode="w", format=tarfile.GNU_FORMAT) as tf:
        for m in mems:
            data = open(os.path.join(MAGS_DIR, m + ".fna"), "rb").read()
            ti = tarfile.TarInfo(m + ".fna")
            ti.size, ti.mtime, ti.mode, ti.uid, ti.gid, ti.uname, ti.gname = len(data), 0, 0o644, 0, 0, "", ""
            tf.addfile(ti, io.BytesIO(data))
    xz = lzma.compress(raw.getvalue(), format=lzma.FORMAT_XZ, preset=6)
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:
        fh.write(xz)
    os.replace(tmp, path)
    sha = hashlib.sha256(xz).hexdigest()
    db.execute("UPDATE archives SET path=?, sha256=?, state='built', updated=? WHERE name=?", (path, sha, now(), name))
    db.commit()
    return path, sha


def verify_object(oid, mems):
    r = http("GET", f"/objects/{oid}")
    if r.status_code != 200:
        return f"GET /objects -> {r.status_code}"
    o = r.json()
    if o.get("object_type") != "GenomeArchive" or o.get("visibility") != "private":
        return f"object type/visibility {o.get('object_type')}/{o.get('visibility')}"
    r = http("GET", f"/objects/{oid}/files")
    if r.status_code != 200:
        return f"GET files -> {r.status_code}"
    files = {f["name"]: f.get("size_bytes") for f in r.json()}
    want = {m + ".fna": os.path.getsize(os.path.join(MAGS_DIR, m + ".fna")) for m in mems}
    if files != want:
        return f"file listing mismatch: {len(files)} listed vs {len(want)} expected"
    return None


def ensure_uploaded(name):
    """Upload (idempotently) + verify one archive. Returns object_id when verified, else None."""
    state, path, sha, ujob = db.execute("SELECT state, path, sha256, upload_job FROM archives WHERE name=?",
                                        (name,)).fetchone()
    mems = members_of(name)
    if state == "verified":
        return sha
    if state not in ("planned", "built", "uploading"):
        return None
    if state == "planned" or not path or not os.path.exists(path):
        path, sha = build_archive(name)
    if a.dry_run:
        ev(f"[dry run] archive {name}: {len(mems)} members, sha256 {sha[:16]}…, {os.path.getsize(path)/1e6:.1f} MB")
        return sha
    if ujob is None:
        r = http("GET", f"/objects/{sha}")
        if r.status_code != 200:
            with open(path, "rb") as fh:
                try:
                    r = http("POST", "/objects", files={"file": (name + ".tar.xz", fh)}, timeout=1800,
                             data={"object_type": "GenomeArchive", "name": name, "visibility": "private",
                                   "transform": "none"})
                except requests.RequestException as e:
                    ev(f"WARN upload POST {name}: {e!s:.150}; will retry (GET /objects/{{sha}} makes it idempotent)")
                    return None
            if r.status_code == 202:
                j = r.json()
                if j.get("expected_object_id") != sha:
                    return halt(f"upload {name}: expected_object_id {j.get('expected_object_id')} != local {sha}")
                ujob = j["job_id"]
                db.execute("UPDATE archives SET upload_job=?, state='uploading', updated=? WHERE name=?", (ujob, now(), name))
                db.commit()
            elif r.status_code != 201:
                return halt(f"upload {name}: HTTP {r.status_code} {r.text[:200]}")
    if ujob:
        st = http("GET", f"/jobs/{ujob}").json()
        if st["state"] not in TERMINAL:
            return None
        if st["state"] == "failed":
            fail = st.get("failure") or {}
            n = int(kv_get(f"upload_fail:{name}", 0)) + 1
            kv_set(f"upload_fail:{name}", n)
            db.execute("UPDATE archives SET upload_job=NULL, state='built', updated=? WHERE name=?", (now(), name))
            db.commit()
            if fail.get("code") in ("invalid_input", "validation_failed", "cancelled") or n > 2:
                return halt(f"upload job {ujob} for {name} failed (attempt {n}): {fail}; inspect, then clear kv "
                            f"'halted' (and 'upload_fail:{name}' to allow more re-uploads) and restart")
            ev(f"WARN upload job {ujob} for {name} failed (attempt {n}): {fail}; will re-upload")
            return None
        res = http("GET", f"/jobs/{ujob}/result").json()
        if res.get("object_id") != sha:
            return halt(f"upload {name}: result object_id {res.get('object_id')} != {sha}")
    err = verify_object(sha, mems)
    if err:
        return halt(f"upload {name}: {err}")
    db.execute("UPDATE archives SET object_id=?, state='verified', updated=? WHERE name=?", (sha, now(), name))
    db.commit()
    ev(f"ARCHIVE {name} verified ({len(mems)} members, {os.path.getsize(path)/1e6:.1f} MB)")
    return sha


# ---------------- submission ----------------
def payload_for(oid, genome_ids=None):
    if not SKANI_DB or not SKANI_DB.startswith("ref:") and len(SKANI_DB) != 64:
        raise SystemExit("no valid --skani-db (GTDB R232 SKANIDB object id) configured; not submitting")
    arch = {"object_id": oid}
    if genome_ids:
        arch["genome_ids"] = list(genome_ids)
    pl = {"schema_version": "1", "job_type": "KBDLBuildGenome",
          "params": {"archive": arch, "gene_caller": "bakta", "call_non_coding": True,
                     "skani_db": {"object_id": SKANI_DB}, "bakta_db": {"object_id": BAKTA_DB}}}
    pr = pl["params"]
    assert set(pr) <= ALLOWED_KEYS and pr["gene_caller"] == "bakta" and pr["call_non_coding"] is True
    assert pr["skani_db"] == {"object_id": SKANI_DB} and SKANI_DB != GTDB_R220 and pr["bakta_db"] == {"object_id": BAKTA_DB}
    assert set(arch) <= {"object_id", "genome_ids"} and len(oid) == 64
    return pl


def submit(name, oid, genome_ids=None):
    """POST one build job. Returns job_id, or None (dry run, lost -> archive parked as 'lost', or halted)."""
    pl = payload_for(oid, genome_ids)
    if a.dry_run:
        ev(f"[dry run] would POST /jobs {json.dumps(pl)}")
        return None
    attempt = 1 + db.execute("SELECT count(*) FROM jobs WHERE archive=? AND genome_ids IS NULL", (name,)).fetchone()[0] \
        if not genome_ids else 2
    cur = db.execute("INSERT INTO intents(archive, payload, created) VALUES(?,?,?)", (name, json.dumps(pl), time.time()))
    iid = cur.lastrowid
    db.execute("UPDATE archives SET state='lost', updated=? WHERE name=?", (now(), name))   # until the job is recorded
    db.commit()
    try:
        r = http("POST", "/jobs", json=pl, retry=False)
    except requests.RequestException as e:
        db.execute("UPDATE intents SET outcome='unknown' WHERE id=?", (iid,))
        db.commit()
        ev(f"WARN POST /jobs for {name} lost ({e!s:.120}); will reconcile")
        return None
    if r.status_code >= 500:
        db.execute("UPDATE intents SET outcome='unknown' WHERE id=?", (iid,))
        db.commit()
        ev(f"WARN POST /jobs for {name} -> {r.status_code}; will reconcile")
        return None
    if r.status_code != 202:
        db.execute("UPDATE intents SET outcome='rejected' WHERE id=?", (iid,))
        db.execute("UPDATE archives SET state='failed' WHERE name=?", (name,))
        db.commit()
        return halt(f"POST /jobs {name} -> {r.status_code} {r.text[:300]}")
    try:
        jid = r.json()["job_id"]
        db.execute("UPDATE intents SET posted=?, job_id=?, outcome='ok' WHERE id=?", (time.time(), jid, iid))
        db.execute("INSERT INTO jobs(job_id, archive, genome_ids, submitted, state, attempt) VALUES(?,?,?,?,?,?)",
                   (jid, name, json.dumps(genome_ids) if genome_ids else None, time.time(), "queued", attempt))
        db.execute("UPDATE archives SET state='submitted', updated=? WHERE name=?", (now(), name))
        if genome_ids:
            db.execute("DELETE FROM kv WHERE k=?", (f"subset_pending:{name}",))
        db.commit()
    except Exception:
        db.rollback()
        db.execute("UPDATE intents SET outcome='unknown' WHERE id=?", (iid,))
        db.commit()
        raise
    ev(f"SUBMIT {name} -> job {jid}" + (f" (subset of {len(genome_ids)})" if genome_ids else f" (attempt {attempt})"))
    return jid


def reconcile_intents():
    """Resolve intents whose POST outcome is unknown (lost response, crash between POST and commit, 5xx)."""
    pend = db.execute("SELECT id, archive, created, payload FROM intents WHERE outcome='unknown' "
                      "OR (outcome IS NULL AND job_id IS NULL) ORDER BY created").fetchall()
    if not pend:
        return 0
    lst = http("GET", "/jobs").json()
    known = {r[0] for r in db.execute("SELECT job_id FROM jobs")}
    t0 = pend[0][2] - 5

    def ts(j):
        try:
            return calendar.timegm(time.strptime(j["submitted_at"][:19], "%Y-%m-%dT%H:%M:%S"))
        except Exception:
            return 0
    orphans = sorted((j for j in lst if j.get("job_type") == "KBDLBuildGenome" and j["job_id"] not in known
                      and ts(j) >= t0), key=ts)
    if orphans and len(orphans) != len(pend):
        halt(f"RECONCILE ambiguous: orphan jobs {[j['job_id'] for j in orphans]} vs intents {[x[0] for x in pend]}. "
             "Resolve by hand: INSERT each job into jobs(job_id, archive, genome_ids, attempt), set intents.outcome="
             "'adopted', set archives.state='submitted', then clear kv halted")
        return len(pend)
    unresolved = 0
    for i, (iid, name, created, payload) in enumerate(pend):
        gids = json.loads(payload)["params"]["archive"].get("genome_ids")
        if orphans:
            jid = orphans[i]["job_id"]
            att = 2 if gids else 1 + db.execute("SELECT count(*) FROM jobs WHERE archive=? AND genome_ids IS NULL",
                                                (name,)).fetchone()[0]
            db.execute("UPDATE intents SET job_id=?, outcome='adopted' WHERE id=?", (jid, iid))
            db.execute("INSERT INTO jobs(job_id, archive, genome_ids, submitted, state, attempt) VALUES(?,?,?,?,?,?)",
                       (jid, name, json.dumps(gids) if gids else None, created, "queued", att))
            db.execute("UPDATE archives SET state='submitted' WHERE name=?", (name,))
            if gids:
                db.execute("DELETE FROM kv WHERE k=?", (f"subset_pending:{name}",))
            ev(f"RECONCILE adopted orphan job {jid} for {name}" + (f" (subset of {len(gids)})" if gids else ""))
        elif time.time() - created >= RECONCILE_GRACE:
            db.execute("UPDATE intents SET outcome='never_created' WHERE id=?", (iid,))
            db.execute("UPDATE archives SET state=? WHERE name=?", ("partial" if gids else "verified", name))
            ev(f"RECONCILE intent {iid} for {name}: no job was created; " +
               ("subset retry stays pending" if gids else "archive will be resubmitted"))
        else:
            unresolved += 1
    db.commit()
    return unresolved


# ---------------- verification ----------------
def rehash_contigs(contigs):
    lines = []
    for c in contigs:
        lines.append(">" + str(c.get("source_id")))
        s = str(c.get("sequence") or "").upper()
        lines += [s[i:i + 60] for i in range(0, len(s), 60)]
    return hashlib.sha256(("\n".join(lines) + "\n").encode()).hexdigest() if lines else ""


def check_genome(mag, e):
    """Verify the stored Genome JSON. Returns (cds, n_nc, flags, problems); problems == None means 'unreadable now'."""
    try:
        r = http("GET", f"/objects/{e['genome_object_id']}/content", timeout=600)
    except RuntimeError:
        return None, None, [], None
    if r.status_code != 200:
        return None, None, [], None if r.status_code >= 500 else [f"genome content HTTP {r.status_code}"]
    g = r.json()
    bp, nc_, pred = db.execute("SELECT bp, n_contigs, predicted_hash FROM members WHERE mag=?", (mag,)).fetchone()
    feats, ncf, prov = g.get("features") or [], g.get("non_coding_features") or [], g.get("non_coding_provenance") or {}
    contigs = g.get("contigs") or []
    problems, flags = [], []
    if g.get("genome_id") != pred[:16]:
        problems.append(f"Genome JSON genome_id {g.get('genome_id')} != hash[:16]")
    if rehash_contigs(contigs) != pred:
        problems.append("Genome JSON contigs do not re-hash to the MAG's content hash")
    if not feats:
        problems.append("0 CDS (gene-free build)")
    elif len(feats) / max(bp / 1e6, 1e-9) < 500:
        flags.append(f"low CDS density {len(feats)/(bp/1e6):.0f}/Mb")
    if len(contigs) != nc_ or sum(int(c.get("length") or 0) for c in contigs) != bp:
        problems.append(f"contigs {len(contigs)}/{nc_} or bp mismatch")
    if len(ncf) != sum(int(v) for v in (e.get("non_coding_counts") or {}).values()):
        problems.append(f"non_coding_features {len(ncf)} != counts sum")
    for k, v in EXPECT_PROV.items():
        if prov.get(k) != v:
            problems.append(f"provenance {k}={prov.get(k)!r} expected {v!r}")
    cd = prov.get("cds_discarded")
    if not (isinstance(cd, dict) and set(cd) <= {"overlap", "spurious", "length"}
            and all(isinstance(x, int) and x >= 0 for x in cd.values())):
        problems.append(f"provenance cds_discarded {cd!r}")
    ko = prov.get("cds_kept_overlapping_rna")
    if not (ko is None or (isinstance(ko, int) and ko >= 0)):
        problems.append(f"provenance cds_kept_overlapping_rna {ko!r}")
    return len(feats), len(ncf), flags, problems


def record_member(mag, jid, e, cds, nnc, flags, problems):
    tax = e.get("taxonomy") or {}
    unreadable = problems is None
    db.execute("""UPDATE members SET job_id=?, genome_id=?, genome_object_id=?, fasta_object_id=?, prefix=?,
                  gene_caller=?, non_coding_source=?, counts_json=?, taxonomy_json=?, cds=?, n_noncoding=?,
                  flags=?, verified=?, error_code=?, error_msg=?, updated=? WHERE mag=?""",
               (jid, e.get("genome_id"), e.get("genome_object_id"), e.get("fasta_object_id"), e.get("prefix"),
                e.get("gene_caller"), e.get("non_coding_source"), json.dumps(e.get("non_coding_counts")),
                json.dumps(tax), cds, nnc, ";".join(flags) or None, 0 if (unreadable or problems) else 1,
                "content_unavailable" if unreadable else None,
                None if unreadable else ("; ".join(problems)[:500] or None), now(), mag))


def harvest(jid, name, state, failure, genome_ids):
    if state == "failed":
        code = (failure or {}).get("code")
        att = db.execute("SELECT attempt FROM jobs WHERE job_id=?", (jid,)).fetchone()[0]
        db.execute("UPDATE jobs SET harvested=1, failure=? WHERE job_id=?", (json.dumps(failure), jid))
        ev(f"JOBFAIL {name} job {jid} (attempt {att}): {failure}")
        if code in RETRYABLE and att < MAX_ATTEMPTS:
            if genome_ids:
                kv_set(f"subset_pending:{name}", genome_ids)
                db.execute("UPDATE archives SET state='partial' WHERE name=?", (name,))
            else:
                db.execute("UPDATE archives SET state='verified' WHERE name=?", (name,))   # resubmit whole archive
        else:
            db.execute("UPDATE archives SET state='failed' WHERE name=?", (name,))
            halt(f"job {jid} ({name}) failed: {failure}")
        db.commit()
        return
    r = http("GET", f"/jobs/{jid}/result", timeout=600)
    if r.status_code != 200:
        ev(f"WARN GET /jobs/{jid}/result -> {r.status_code}; will re-poll")
        return
    res = r.json()
    with open(os.path.join(WORK, "results", f"{jid}.json"), "w") as fh:
        json.dump(res, fh)
    want = json.loads(genome_ids) if genome_ids else members_of(name)
    gen, errs = res.get("genomes") or {}, res.get("errors") or {}
    if set(gen) | set(errs) != set(want) or set(gen) & set(errs):
        db.execute("UPDATE jobs SET harvested=1 WHERE job_id=?", (jid,))
        db.execute("UPDATE archives SET state='failed' WHERE name=?", (name,))
        db.commit()
        return halt(f"job {jid}: result keys do not partition the submitted members")
    pl = db.execute("SELECT payload FROM intents WHERE job_id=?", (jid,)).fetchone()
    job_db = json.loads(pl[0])["params"]["skani_db"]["object_id"] if pl else GTDB_R220
    retry_subset, flag_count, no_trna_big, low_cds, edge_job = [], 0, 0, 0, 0
    n_big = 0
    for mag in want:
        if mag in errs:
            msg = str(errs[mag].get("message"))
            db.execute("UPDATE members SET job_id=?, error_code=?, error_msg=?, verified=0, updated=? WHERE mag=?",
                       (jid, errs[mag].get("code"), msg[:500], now(), mag))
            if msg.startswith("Bakta:") and "returned an out-of-range interval" in msg:
                edge = json.loads(kv_get("bakta_edge_mags", "[]"))
                if mag not in edge:
                    edge.append(mag)
                    kv_set("bakta_edge_mags", json.dumps(edge))
                edge_job += 1
                ev(f"BAKTA-EDGE {mag} ({len(edge)} in all): {msg[:160]}")
                nharv = db.execute("SELECT count(*) FROM members WHERE job_id IS NOT NULL").fetchone()[0]
                if edge_job > BAKTA_EDGE_PER_JOB or len(edge) > max(BAKTA_EDGE_FLOOR, BAKTA_EDGE_RATE * nharv):
                    halt(f"Bakta edge-interval failures: {edge_job} in job {jid}, {len(edge)} of {nharv} harvested")
            elif msg.startswith("Bakta:"):
                halt(f"Bakta error for {mag}: {msg[:200]}")
            elif "skani failed" in msg:
                retry_subset.append(mag)
            ev(f"MEMBER-ERROR {mag}: {errs[mag].get('code')} {msg[:160]}")
            continue
        e = gen[mag]
        pred = db.execute("SELECT predicted_hash FROM members WHERE mag=?", (mag,)).fetchone()[0]
        pre = []
        if e.get("fasta_object_id") != pred:
            pre.append(f"fasta_object_id {e.get('fasta_object_id')} != predicted {pred}")
        if e.get("gene_caller") != "bakta":
            pre.append(f"gene_caller {e.get('gene_caller')}")
        if e.get("non_coding_source") != "called":
            pre.append(f"non_coding_source {e.get('non_coding_source')}")
        tax = e.get("taxonomy") or {}
        if tax.get("ani") is not None and tax.get("source_database") != job_db:
            pre.append(f"taxonomy source_database {tax.get('source_database')}")
        cds, nnc, flags, problems = check_genome(mag, e)
        if problems is not None:
            problems = pre + problems
        elif pre:
            problems = pre
        bp = db.execute("SELECT bp FROM members WHERE mag=?", (mag,)).fetchone()[0]
        if int((e.get("non_coding_counts") or {}).get("tRNA") or 0) == 0:
            flags.append("no tRNA")
            if bp >= 1e6:
                no_trna_big += 1
        if bp >= 1e6:
            n_big += 1
        if any(f.startswith("low CDS density") for f in flags):
            low_cds += 1
        flag_count += bool(flags)
        record_member(mag, jid, e, cds, nnc, flags, problems)
        if problems:
            halt(f"verification failed for {mag}: {'; '.join(problems)[:300]}")
    db.execute("UPDATE jobs SET harvested=1 WHERE job_id=?", (jid,))
    nver = db.execute("SELECT count(*) FROM members WHERE archive=? AND verified=1", (name,)).fetchone()[0]
    ntot = db.execute("SELECT count(*) FROM members WHERE archive=?", (name,)).fetchone()[0]
    if retry_subset:
        if genome_ids:          # a subset retry failed skani again: give up on those members
            ev(f"WARN {name}: skani failed again for {len(retry_subset)} member(s) after the subset retry: {retry_subset}")
        else:
            kv_set(f"subset_pending:{name}", json.dumps(retry_subset))
    db.execute("UPDATE archives SET state=? WHERE name=?", ("done" if nver == ntot else "partial", name))
    db.commit()
    ev(f"HARVEST {name} job {jid}: {len(gen)} built, {len(errs)} errors, {nver}/{ntot} verified"
       + (f"; flags on {flag_count} (no tRNA in {no_trna_big} of {n_big} MAGs >= 1 Mb, low CDS density {low_cds})"
          if flag_count else ""))
    if (n_big and no_trna_big > n_big / 2) or low_cds > 2:
        halt(f"quality breaker on {name}: no tRNA in {no_trna_big}/{n_big} MAGs >= 1 Mb, {low_cds} below 500 CDS/Mb")


def recheck_unreadable():
    if time.time() - float(kv_get("last_recheck", 0)) < RECHECK_EVERY:
        return
    kv_set("last_recheck", time.time())
    for mag, jid, goid, name in db.execute("SELECT m.mag, m.job_id, m.genome_object_id, m.archive FROM members m "
                                           "WHERE m.error_code='content_unavailable'").fetchall():
        res = json.load(open(os.path.join(WORK, "results", f"{jid}.json")))
        e = res["genomes"][mag]
        cds, nnc, flags, problems = check_genome(mag, e)
        if problems is None:
            continue
        record_member(mag, jid, e, cds, nnc, flags, problems)
        if problems:
            halt(f"verification failed for {mag} on re-check: {'; '.join(problems)[:300]}")
        nver = db.execute("SELECT count(*) FROM members WHERE archive=? AND verified=1", (name,)).fetchone()[0]
        ntot = db.execute("SELECT count(*) FROM members WHERE archive=?", (name,)).fetchone()[0]
        if nver == ntot:
            db.execute("UPDATE archives SET state='done' WHERE name=?", (name,))
        db.commit()
        ev(f"RECHECK {mag}: {'verified' if not problems else 'FAILED'}")


def poll():
    act = db.execute("SELECT job_id, archive, genome_ids, submitted, alerted FROM jobs WHERE harvested=0").fetchall()
    if act:
        r = http("POST", "/jobs/status", json={"job_ids": [x[0] for x in act]}, retry=True)
        st = r.json().get("jobs", {})
        for jid, name, gids, sub, alerted in act:
            s = st.get(jid)
            if not s:
                continue
            state = s.get("state")
            db.execute("UPDATE jobs SET state=? WHERE job_id=?", (state, jid))
            db.commit()
            if state in TERMINAL:
                harvest(jid, name, state, s.get("failure"), gids)
                continue
            if alerted:
                continue
            pred = db.execute("SELECT pred_s FROM archives WHERE name=?", (name,)).fetchone()[0] or 0
            started = s.get("started_at")
            if state == "running" and started:
                try:
                    t0 = calendar.timegm(time.strptime(started[:19], "%Y-%m-%dT%H:%M:%S"))
                except Exception:
                    t0 = sub
                if time.time() - t0 > 3 * pred * 1.11 + 1800 + 600:
                    ev(f"ALERT job {jid} ({name}) running {(time.time()-t0)/3600:.1f} h vs predicted {pred/3600:.2f} h: "
                       f"ask Chris to inspect the Bakta container (do NOT cancel)")
                    db.execute("UPDATE jobs SET alerted=1 WHERE job_id=?", (jid,))
            elif state == "queued" and time.time() - sub > 12 * 3600:
                ev(f"ALERT job {jid} ({name}) queued for {(time.time()-sub)/3600:.1f} h")
                db.execute("UPDATE jobs SET alerted=1 WHERE job_id=?", (jid,))
        db.commit()
    inflight = db.execute("SELECT count(*) FROM jobs WHERE harvested=0 AND state NOT IN ('completed','failed')").fetchone()[0]
    if inflight:
        last = float(kv_get("last_harvest", time.time()))
        if time.time() - last > 6 * 3600 and time.time() - float(kv_get("last_progress_alert", 0)) > 6 * 3600:
            ev(f"ALERT no harvest for {(time.time()-last)/3600:.1f} h with {inflight} job(s) in flight")
            kv_set("last_progress_alert", time.time())
    return inflight


def phase_status(ph):
    rows = dict(db.execute("SELECT state, count(*) FROM archives WHERE phase=? GROUP BY state", (ph,)).fetchall())
    return rows


def gate_open(ph):
    """True if phase ph may start: every archive of the previous phase is done (or the operator accepted it)."""
    i = PHASES.index(ph)
    if i == 0:
        return True
    prev = PHASES[i - 1]
    st = phase_status(prev)
    if set(st) <= {"done"} or kv_get(f"accept_phase:{prev}") == "1":
        return True
    return False


def work_items(limit):
    """Next archives to upload/submit within --through, respecting phase gates."""
    out = []
    for ph in PHASES[:PHASES.index(a.through) + 1]:
        if not gate_open(ph):
            break
        out += db.execute("SELECT name FROM archives WHERE phase=? AND state IN ('planned','built','uploading',"
                          "'verified') ORDER BY name", (ph,)).fetchall()
        if len(out) >= limit:
            break
    return [r[0] for r in out[:limit]]


# ---------------- main loop ----------------
ev(f"START driver v2 pid {os.getpid()} through={a.through} max_in_flight={a.max_in_flight} dry_run={a.dry_run}")
if a.dry_run:
    for name in work_items(3):
        oid = ensure_uploaded(name)
        if oid:
            submit(name, oid)
    sys.exit(0)
if not db.execute("SELECT 1 FROM kv WHERE k='last_harvest'").fetchone():
    kv_set("last_harvest", time.time())
while True:
    try:
        unresolved = reconcile_intents()
        n_before = db.execute("SELECT count(*) FROM jobs WHERE harvested=1").fetchone()[0]
        inflight = poll()
        if db.execute("SELECT count(*) FROM jobs WHERE harvested=1").fetchone()[0] > n_before:
            kv_set("last_harvest", time.time())
        recheck_unreadable()
        inflight += unresolved
        # 1. pending subset retries first (under the in-flight cap)
        for k, gids in db.execute("SELECT k, v FROM kv WHERE k LIKE 'subset_pending:%'").fetchall():
            if halted() or unresolved or inflight >= a.max_in_flight:
                break
            name = k.split(":", 1)[1]
            oid = db.execute("SELECT object_id FROM archives WHERE name=?", (name,)).fetchone()[0]
            if submit(name, oid, json.loads(gids)):
                inflight += 1
        # 2. keep a few archives uploaded ahead
        for name in work_items(a.prefetch + a.max_in_flight):
            if halted():
                break
            ensure_uploaded(name)
        # 3. submit verified archives
        for name in work_items(a.max_in_flight * 2):
            if halted() or unresolved or inflight >= a.max_in_flight:
                break
            st, oid = db.execute("SELECT state, object_id FROM archives WHERE name=?", (name,)).fetchone()
            if st == "verified" and oid:
                if submit(name, oid):
                    inflight += 1
        # 4. exit conditions
        active = db.execute("SELECT count(*) FROM jobs WHERE harvested=0").fetchone()[0]
        pending = db.execute("SELECT count(*) FROM kv WHERE k LIKE 'subset_pending:%'").fetchone()[0]
        open_intents = db.execute("SELECT count(*) FROM intents WHERE outcome='unknown' "
                                  "OR (outcome IS NULL AND job_id IS NULL)").fetchone()[0]
        if not active and not open_intents and not pending:
            if halted():
                ev(f"STOPPED (halted: {kv_get('halt_reason')}); fix the cause, DELETE FROM kv WHERE k='halted', restart")
                break
            if not work_items(1):
                blocked = [ph for ph in PHASES[1:PHASES.index(a.through) + 1] if not gate_open(ph)]
                if blocked:
                    prev = PHASES[PHASES.index(blocked[0]) - 1]
                    ev(f"STOPPED GATE: phase '{blocked[0]}' blocked because '{prev}' is not all done "
                       f"({phase_status(prev)}); inspect the errors, then set kv accept_phase:{prev}=1 to proceed")
                else:
                    ev(f"DONE through={a.through}: " + json.dumps(dict(db.execute(
                        "SELECT state, count(*) FROM archives GROUP BY state").fetchall())))
                break
    except Exception as e:
        ev(f"WARN loop error: {e!r:.300}")
    time.sleep(60)
