#!/usr/bin/env python3
"""Build every wetland MAG in KBDL with Bakta (handoff Phase 2; per Andrew 2026-11-15: gene_caller=bakta, all MAGs).

Runs ON POPLAR against the KBDL job service at http://localhost:8791 (Chris Henry's service). Inputs:
  archive_plan.tsv  517 archives (canary 2 MAGs, pilot 23, giant 1, bulk 25 each) with each MAG's predicted KBDL content
                    hash (= fasta_object_id); MAG FASTAs in MAGS_DIR/<MAG>.fna.
Per job: KBDLBuildGenome {archive, gene_caller='bakta', call_non_coding=true, skani_db=GTDB r220 ref, bakta_db=6.0-full}.
These are permanent: KBDL dedups builds by content hash across owners, so the first build fixes caller and taxonomy.

Safety (from the 2026-11-15 contract review, ~/gtdbtk_kbase/next_kbdl-build-contract.md):
  * every payload is asserted (key set, gene_caller, call_non_coding, refs) right before POST /jobs;
  * POST /jobs is never blindly retried: an intent row is written first; a lost response is reconciled from GET /jobs;
  * uploads: deterministic .tar.xz, GET /objects/{sha} first (idempotent), 202 is not "stored": upload job result ==
    local sha256, GET /objects/{id} (GenomeArchive, private), /files names+sizes == members;
  * per member: in genomes xor errors; fasta_object_id == predicted hash; gene_caller == bakta; non_coding_source ==
    called; tRNA > 0 (flag); Genome JSON content: >= 1 CDS (flag < 500/Mb), contigs/bp match, non-coding count matches,
    provenance matches 09_expected_genome_provenance.json;
  * circuit breaker (stop submitting, keep polling) on any 'Bakta:' error, caller/hash/provenance mismatch, a whole-job
    invalid_input / dependency_unavailable / tool_failure;
  * watchdog: a job running > 3x its predicted time + 30 min raises an ALERT only (never cancel: KBDL keeps building).
Stages: canary -> pilot -> bulk; the driver stops after --through (default: canary) so results can be inspected.
Ledger: SQLite on local disk (WORK, never NFS). Token: read from stdin, then the process detaches; never written to disk.
Usage (poplar): python3 -I kbdl_build_driver.py --through canary|pilot|bulk [--max-in-flight 3] [--dry-run] < token
"""
import argparse, calendar, hashlib, io, json, lzma, os, sqlite3, sys, tarfile, time
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

p = argparse.ArgumentParser()
p.add_argument("--through", choices=PHASES, default="canary")
p.add_argument("--max-in-flight", type=int, default=3)
p.add_argument("--prefetch", type=int, default=6, help="archives uploaded+verified ahead of submission")
p.add_argument("--dry-run", action="store_true", help="build archives + hashes, print payloads, no HTTP writes")
p.add_argument("--foreground", action="store_true")
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
S = requests.Session()
S.headers["Authorization"] = "Bearer " + TOKEN
EXPECT_PROV = {k: v for k, v in json.load(open(os.path.join(HERE, "09_expected_genome_provenance.json"))).items()
               if not k.startswith("_")}


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


# ---------------- HTTP helpers ----------------
def http(method, path, retry=None, **kw):
    retry = (method == "GET") if retry is None else retry
    for i in range(8):
        try:
            r = S.request(method, BASE + path, timeout=kw.pop("timeout", 120), **kw)
        except requests.RequestException as e:
            if not retry:
                raise
            ev(f"WARN {method} {path}: {e!s:.150}; retry")
            time.sleep(min(300, 15 * 2 ** i))
            continue
        if r.status_code == 401:
            ev("ALERT 401 from KBDL (token?); pausing 30 min")
            time.sleep(1800)
            continue
        if r.status_code == 503 and retry:
            time.sleep(min(300, 15 * 2 ** i))
            continue
        return r
    raise RuntimeError(f"{method} {path}: giving up")


# ---------------- archives ----------------
def build_archive(name):
    path = os.path.join(WORK, "archives", name + ".tar.xz")
    mems = [r[0] for r in db.execute("SELECT mag FROM members WHERE archive=? ORDER BY mag", (name,))]
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
    return path, sha, mems


def verify_object(name, oid, mems):
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
    """Upload (idempotently) + verify one archive. Returns object_id or None (and halts on mismatch)."""
    row = db.execute("SELECT state, path, sha256, upload_job FROM archives WHERE name=?", (name,)).fetchone()
    state, path, sha, ujob = row
    mems = [r[0] for r in db.execute("SELECT mag FROM members WHERE archive=? ORDER BY mag", (name,))]
    if state == "verified":
        return sha
    if state == "planned" or not path or not os.path.exists(path):
        path, sha, mems = build_archive(name)
    if a.dry_run:
        ev(f"[dry run] archive {name}: {len(mems)} members, sha256 {sha[:16]}…, {os.path.getsize(path)/1e6:.1f} MB")
        return sha
    if ujob is None:
        r = http("GET", f"/objects/{sha}")
        if r.status_code != 200:
            with open(path, "rb") as fh:
                r = http("POST", "/objects", files={"file": (name + ".tar.xz", fh)},
                         data={"object_type": "GenomeArchive", "name": name, "visibility": "private", "transform": "none"},
                         timeout=1800)
            if r.status_code == 201:
                pass
            elif r.status_code == 202:
                j = r.json()
                if j.get("expected_object_id") != sha:
                    return halt(f"upload {name}: expected_object_id {j.get('expected_object_id')} != local {sha}")
                ujob = j["job_id"]
                db.execute("UPDATE archives SET upload_job=?, state='uploading', updated=? WHERE name=?", (ujob, now(), name))
                db.commit()
            else:
                return halt(f"upload {name}: HTTP {r.status_code} {r.text[:200]}")
    if ujob:
        st = http("GET", f"/jobs/{ujob}").json()
        if st["state"] not in TERMINAL:
            return None
        if st["state"] == "failed":
            db.execute("UPDATE archives SET upload_job=NULL, state='built' WHERE name=?", (name,))
            db.commit()
            ev(f"WARN upload job {ujob} for {name} failed: {st.get('failure')}; will re-upload")
            return None
        res = http("GET", f"/jobs/{ujob}/result").json()
        if res.get("object_id") != sha:
            return halt(f"upload {name}: result object_id {res.get('object_id')} != {sha}")
    err = verify_object(name, sha, mems)
    if err:
        return halt(f"upload {name}: {err}")
    db.execute("UPDATE archives SET object_id=?, state='verified', updated=? WHERE name=?", (sha, now(), name))
    db.commit()
    ev(f"ARCHIVE {name} verified ({len(mems)} members, {os.path.getsize(path)/1e6:.1f} MB)")
    return sha


# ---------------- jobs ----------------
def halt(msg):
    if kv_get("halted") != "1":
        ev(f"HALT {msg}")
    kv_set("halted", "1")
    kv_set("halt_reason", msg)
    return None


def payload_for(oid, genome_ids=None):
    arch = {"object_id": oid}
    if genome_ids:
        arch["genome_ids"] = list(genome_ids)
    pl = {"schema_version": "1", "job_type": "KBDLBuildGenome",
          "params": {"archive": arch, "gene_caller": "bakta", "call_non_coding": True,
                     "skani_db": {"object_id": GTDB_R220}, "bakta_db": {"object_id": BAKTA_DB}}}
    pr = pl["params"]
    assert set(pr) <= ALLOWED_KEYS and pr["gene_caller"] == "bakta" and pr["call_non_coding"] is True
    assert pr["skani_db"] == {"object_id": GTDB_R220} and pr["bakta_db"] == {"object_id": BAKTA_DB}
    assert set(arch) <= {"object_id", "genome_ids"}
    return pl


def submit(name, oid, genome_ids=None, attempt=1):
    pl = payload_for(oid, genome_ids)
    if a.dry_run:
        ev(f"[dry run] would POST /jobs {json.dumps(pl)}")
        return None
    cur = db.execute("INSERT INTO intents(archive, payload, created) VALUES(?,?,?)", (name, json.dumps(pl), time.time()))
    iid = cur.lastrowid
    db.commit()
    try:
        r = S.post(BASE + "/jobs", json=pl, timeout=120)
    except requests.RequestException as e:
        db.execute("UPDATE intents SET outcome='unknown' WHERE id=?", (iid,))
        db.commit()
        ev(f"WARN POST /jobs for {name} lost ({e!s:.120}); will reconcile")
        return None
    if r.status_code != 202:
        db.execute("UPDATE intents SET outcome='rejected' WHERE id=?", (iid,))
        db.commit()
        return halt(f"POST /jobs {name} -> {r.status_code} {r.text[:300]}")
    jid = r.json()["job_id"]
    db.execute("UPDATE intents SET posted=?, job_id=?, outcome='ok' WHERE id=?", (time.time(), jid, iid))
    db.execute("INSERT INTO jobs(job_id, archive, genome_ids, submitted, state, attempt) VALUES(?,?,?,?,?,?)",
               (jid, name, json.dumps(genome_ids) if genome_ids else None, time.time(), "queued", attempt))
    db.execute("UPDATE archives SET state='submitted', updated=? WHERE name=?", (now(), name))
    db.commit()
    ev(f"SUBMIT {name} -> job {jid}" + (f" (subset {len(genome_ids)}, attempt {attempt})" if genome_ids else ""))
    return jid


def reconcile_intents():
    pend = db.execute("SELECT id, archive, created FROM intents WHERE outcome='unknown'").fetchall()
    if not pend:
        return
    lst = http("GET", "/jobs").json()
    known = {r[0] for r in db.execute("SELECT job_id FROM jobs")}
    orphans = [j for j in lst if j.get("job_type") == "KBDLBuildGenome" and j["job_id"] not in known]
    for iid, name, created in pend:
        if len(orphans) == 1 and len(pend) == 1:
            jid = orphans[0]["job_id"]
            db.execute("UPDATE intents SET job_id=?, outcome='adopted' WHERE id=?", (jid, iid))
            db.execute("INSERT INTO jobs(job_id, archive, submitted, state) VALUES(?,?,?,?)", (jid, name, created, "queued"))
            db.execute("UPDATE archives SET state='submitted' WHERE name=?", (name,))
            ev(f"RECONCILE adopted orphan job {jid} for {name}")
        elif not orphans:
            db.execute("UPDATE intents SET outcome='never_created' WHERE id=?", (iid,))
            ev(f"RECONCILE intent {iid} for {name}: no job was created; archive will be resubmitted")
        else:
            halt(f"RECONCILE ambiguous: {len(orphans)} orphan build jobs for {len(pend)} lost intents")
    db.commit()


def check_genome(mag, e):
    """Fetch the Genome JSON and verify CDS / contigs / non-coding / provenance. Returns (cds, n_nc, flags, problems)."""
    r = http("GET", f"/objects/{e['genome_object_id']}/content", timeout=600)
    if r.status_code != 200:
        return None, None, [], [f"genome content HTTP {r.status_code}"]
    g = r.json()
    bp, nc_ = db.execute("SELECT bp, n_contigs FROM members WHERE mag=?", (mag,)).fetchone()
    feats, ncf, prov = g.get("features") or [], g.get("non_coding_features") or [], g.get("non_coding_provenance") or {}
    contigs = g.get("contigs") or []
    problems, flags = [], []
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
    return len(feats), len(ncf), flags, problems


def harvest(jid, name, state, failure, genome_ids):
    if state == "failed":
        code = (failure or {}).get("code")
        ev(f"JOBFAIL {name} job {jid}: {failure}")
        att = db.execute("SELECT attempt FROM jobs WHERE job_id=?", (jid,)).fetchone()[0]
        if code in ("worker_lost", "internal_error", "store_contention") and att < 3:
            db.execute("UPDATE archives SET state='retry' WHERE name=?", (name,))
        else:
            halt(f"job {jid} ({name}) failed: {failure}")
        return
    res = http("GET", f"/jobs/{jid}/result", timeout=600).json()
    with open(os.path.join(WORK, "results", f"{jid}.json"), "w") as fh:
        json.dump(res, fh)
    want = json.loads(genome_ids) if genome_ids else \
        [r[0] for r in db.execute("SELECT mag FROM members WHERE archive=? ORDER BY mag", (name,))]
    gen, errs = res.get("genomes") or {}, res.get("errors") or {}
    if set(gen) | set(errs) != set(want) or set(gen) & set(errs):
        return halt(f"job {jid}: result keys do not partition the submitted members")
    retry_subset = []
    for mag in want:
        if mag in errs:
            msg = str(errs[mag].get("message"))
            db.execute("UPDATE members SET job_id=?, error_code=?, error_msg=?, updated=? WHERE mag=?",
                       (jid, errs[mag].get("code"), msg[:500], now(), mag))
            if msg.startswith("Bakta:"):
                halt(f"Bakta error for {mag}: {msg[:200]}")
            elif "skani failed" in msg:
                retry_subset.append(mag)
            ev(f"MEMBER-ERROR {mag}: {errs[mag].get('code')} {msg[:160]}")
            continue
        e = gen[mag]
        pred = db.execute("SELECT predicted_hash FROM members WHERE mag=?", (mag,)).fetchone()[0]
        problems = []
        if e.get("fasta_object_id") != pred:
            problems.append(f"fasta_object_id {e.get('fasta_object_id')} != predicted {pred}")
        if e.get("gene_caller") != "bakta":
            problems.append(f"gene_caller {e.get('gene_caller')}")
        if e.get("non_coding_source") != "called":
            problems.append(f"non_coding_source {e.get('non_coding_source')}")
        tax = e.get("taxonomy") or {}
        if tax.get("ani") is not None and tax.get("source_database") != GTDB_R220:
            problems.append(f"taxonomy source_database {tax.get('source_database')}")
        cds, nnc, flags, p2 = check_genome(mag, e)
        problems += p2
        if int((e.get("non_coding_counts") or {}).get("tRNA") or 0) == 0:
            flags.append("no tRNA")
        db.execute("""UPDATE members SET job_id=?, genome_id=?, genome_object_id=?, fasta_object_id=?, prefix=?,
                      gene_caller=?, non_coding_source=?, counts_json=?, taxonomy_json=?, cds=?, n_noncoding=?,
                      flags=?, verified=?, error_code=NULL, error_msg=?, updated=? WHERE mag=?""",
                   (jid, e.get("genome_id"), e.get("genome_object_id"), e.get("fasta_object_id"), e.get("prefix"),
                    e.get("gene_caller"), e.get("non_coding_source"), json.dumps(e.get("non_coding_counts")),
                    json.dumps(tax), cds, nnc, ";".join(flags) or None, 0 if problems else 1,
                    "; ".join(problems)[:500] or None, now(), mag))
        if problems:
            halt(f"verification failed for {mag}: {'; '.join(problems)[:300]}")
    db.execute("UPDATE jobs SET harvested=1 WHERE job_id=?", (jid,))
    nver = db.execute("SELECT count(*) FROM members WHERE archive=? AND verified=1", (name,)).fetchone()[0]
    ntot = db.execute("SELECT count(*) FROM members WHERE archive=?", (name,)).fetchone()[0]
    db.execute("UPDATE archives SET state=? WHERE name=?", ("done" if nver == ntot else "partial", name))
    db.commit()
    ev(f"HARVEST {name} job {jid}: {len(gen)} built, {len(errs)} errors, {nver}/{ntot} verified")
    if retry_subset and kv_get("halted") != "1":
        att = db.execute("SELECT attempt FROM jobs WHERE job_id=?", (jid,)).fetchone()[0]
        if att < 2:
            submit(name, db.execute("SELECT object_id FROM archives WHERE name=?", (name,)).fetchone()[0],
                   retry_subset, attempt=att + 1)


def poll():
    act = db.execute("SELECT job_id, archive, genome_ids, submitted, alerted FROM jobs WHERE harvested=0").fetchall()
    if not act:
        return 0
    r = http("POST", "/jobs/status", json={"job_ids": [x[0] for x in act]}, retry=True)
    st = r.json().get("jobs", {})
    inflight = 0
    for jid, name, gids, sub, alerted in act:
        s = st.get(jid)
        if not s:
            continue
        state = s.get("state")
        db.execute("UPDATE jobs SET state=? WHERE job_id=?", (state, jid))
        if state in TERMINAL:
            db.commit()
            harvest(jid, name, state, s.get("failure"), gids)
            continue
        inflight += 1
        if state == "running" and not alerted:
            pred = db.execute("SELECT pred_s FROM archives WHERE name=?", (name,)).fetchone()[0] or 0
            started = s.get("started_at")
            try:
                t0 = calendar.timegm(time.strptime(started[:19], "%Y-%m-%dT%H:%M:%S")) if started else sub
            except Exception:
                t0 = sub
            if time.time() - t0 > 3 * pred * 1.11 + 1800 + 600:
                ev(f"ALERT job {jid} ({name}) running {(time.time()-t0)/3600:.1f} h vs predicted {pred/3600:.2f} h: "
                   f"ask Chris to inspect the Bakta container (do NOT cancel)")
                db.execute("UPDATE jobs SET alerted=1 WHERE job_id=?", (jid,))
    db.commit()
    return inflight


def phase_done(ph):
    return not db.execute("SELECT count(*) FROM archives WHERE phase=? AND state NOT IN ('done')", (ph,)).fetchone()[0]


def next_archives(limit):
    allowed = PHASES[:PHASES.index(a.through) + 1]
    out = []
    for ph in allowed:
        if PHASES.index(ph) > 0 and not phase_done(PHASES[PHASES.index(ph) - 1]):
            break
        rows = db.execute("SELECT name, state FROM archives WHERE phase=? AND state IN ('planned','built','uploading',"
                          "'verified','retry') ORDER BY name", (ph,)).fetchall()
        out += rows
        if len(out) >= limit:
            break
    return out[:limit]


# ---------------- main loop ----------------
ev(f"START driver pid {os.getpid()} through={a.through} max_in_flight={a.max_in_flight} dry_run={a.dry_run}")
if a.dry_run:
    for name, _ in next_archives(3):
        oid = ensure_uploaded(name)
        submit(name, oid)
    sys.exit(0)
while True:
    try:
        reconcile_intents()
        inflight = poll()
        halted = kv_get("halted") == "1"
        cands = next_archives(a.prefetch + a.max_in_flight)
        for name, state in cands:                     # keep a few archives uploaded ahead
            if halted:
                break
            ensure_uploaded(name)
        for name, state in next_archives(a.max_in_flight * 2):
            if halted or inflight >= a.max_in_flight:
                break
            st = db.execute("SELECT state, object_id FROM archives WHERE name=?", (name,)).fetchone()
            if st[0] in ("verified", "retry") and st[1]:
                att = 1 + db.execute("SELECT count(*) FROM jobs WHERE archive=?", (name,)).fetchone()[0]
                if submit(name, st[1], attempt=att):
                    inflight += 1
                halted = kv_get("halted") == "1"
        if not db.execute("SELECT count(*) FROM jobs WHERE harvested=0").fetchone()[0]:
            if halted:
                ev(f"STOPPED (halted: {kv_get('halt_reason')}); fix, clear kv 'halted', restart")
                break
            if not next_archives(1):
                ev(f"DONE through={a.through}: " + json.dumps(dict(db.execute(
                    "SELECT state, count(*) FROM archives GROUP BY state").fetchall())))
                break
    except Exception as e:
        ev(f"WARN loop error: {e!r:.300}")
    time.sleep(60)
