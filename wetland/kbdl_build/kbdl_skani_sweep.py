#!/usr/bin/env python3
"""skani every wetland MAG against GTDB R232 AND the novel-genome representatives DB in KBDL (per Andrew, 2026-11-17:
"ensure that the novel representatives DB is also run, in addition to R232").

BATCH SIZE MATTERS: at most 100 MAGs per job. skani search builds an approximate marker-hash index when it is given
> 100 query files ("--no-marker-index ... default: load index if > 100 query files"), and that index LOSES hits, more
so the more queries a call has (KBUtilLib passes every query of a job to one skani call). Measured 2026-11-17 on
novel_all_representatives: same 50 MAGs, 131 hits brute-force (50 queries) vs 124 indexed (101 queries; all 7 lost
hits from the 31.6 kb canary MAG); 2,000 queries indexed: 2,476 rows (-n 5) vs 2,372 for the same MAGs inside the
census's 12,859-query run, which also missed bin.001.fasta__SRR10968274's only hits. So archives are grouped up to
MAX_QUERIES = 100 MAGs (4 bulk archives), which keeps every job brute-force. Batches 2-4 of the first launch had 125
MAGs (indexed); they were superseded and re-run.

Why separate KBDLSKANI jobs: KBDLBuildGenome takes one skani DB for the taxonomy it attaches (R232 here) and does not
expose its raw hits. As in Chris Henry's own R232 load (KBDLSKANI over batches of 5 GenomeArchives), the hits for the
clearinghouse landing come from KBDLSKANI jobs with skani_dbs = [gtdb_r232_consolidated, novel_all_representatives].
This driver follows the build driver's ledger READ-ONLY and sweeps every archive once it is uploaded (object_id set),
including the canary and pilot archives built earlier. It never touches the build driver, its ledger or its archives
(delete_archive_on_completion stays false).

KBDLSKANI facts this relies on (KBDLJobRunningPrototype src/kbdl_service/job_types/skani.py, read 2026-11-17):
  results: {query stem: [{genome_id: ref-file stem, distance: 1 - ani (ani a 0-1 fraction), database: db object_id}]},
  one flat list per query over all DBs, best first; a zero-hit query from a searched DB contributes nothing (or []).
  errors: {"<query>::<db object_id>": {...}} per failed entry, "archive::<oid>" for archive problems.
  databases_queried: the requested DBs (a DB that failed wholesale fails the job).  Query stem = archive member
  "<mag>.fna" minus ".fna" = the MAG name.
Per job this driver asserts databases_queried == both DBs (else halt). MAGs with a per-entry error (or in an archive
with an archive:: error) are left out of that batch's hits and retried once as a genome_ids subset batch; a MAG that
fails again is recorded in table failed_mags. A failed job is retried once; a batch that fails twice is recorded the
same way. Halt only above MAX_FAILED_MAGS failed MAGs or MAX_FAILED_BATCHES failed batches (something systematic).
The tail (the last group of < MAX_QUERIES MAGs) is flushed, and the sweep ends, once EVERY archive in the build
ledger (all phases, giant included) is uploaded (object_id set): the sweep needs uploads, not builds, and since
2026-11-18 a separate --upload-only process uploads all archives ahead of the build.

CAVEAT for any downstream use: KBDLSKANI hits carry no alignment fraction. ANI >= 0.95 here is NOT a species match by
the census rule (ANI >= 95 and AF >= 0.5): e.g. bin.009.fasta__SRR16547334 has a 0.9912 R232 hit at AF 0.484, which
GTDB-Tk rejected. Species calls stay with the census; these hits are for the clearinghouse distance table.

Outputs on poplar (WORK): skani_ledger.sqlite; skani_results/<job>.json (raw result); skani_hits/batch_<id>.tsv with
columns mag, archive, db, ref_genome, ani: every hit, plus one row with empty ref_genome/ani for each MAG x DB pair
without a hit (so every MAG has >= 1 row per DB); skani_events.log.
Runs ON POPLAR against http://localhost:8791. Token from stdin; the process then detaches (unless --foreground); the
token is never written to disk.
Usage: python3 -I kbdl_skani_sweep.py [--max-in-flight 4] [--only A,B,...] [--foreground] < token
  --only restricts the sweep to the named archives and lets their last batch be partial (for a test run).
"""
import argparse, calendar, fcntl, json, os, sqlite3, sys, time
import requests

BASE = "http://localhost:8791"
WORK = "/scratch/afreiburger/wetland_kbdl_build"
BUILD_LEDGER = os.path.join(WORK, "ledger.sqlite")
R232 = "ref:eac24851ac1263d5bede542ffafe92aa18f98930542277eb5ebeb8de4a22da5f"
NOVEL = "ref:1be8415b2a773d2632949cc46f1b68d3956c4ce67eab4088a9446fadf5157404"
DBS = [R232, NOVEL]
DBNAME = {R232: "gtdb_r232_consolidated", NOVEL: "novel_all_representatives"}
POLL_S = 120
MAX_QUERIES = 100          # MAGs per KBDLSKANI job: > 100 query files switches skani to its lossy marker index
MAX_FAILED_MAGS = 50
MAX_FAILED_BATCHES = 3
BUILD_PHASES = ("canary", "pilot", "bulk")

p = argparse.ArgumentParser()
p.add_argument("--max-in-flight", type=int, default=4)
p.add_argument("--only", default="")
p.add_argument("--foreground", action="store_true")
a = p.parse_args()
ONLY = {x for x in a.only.split(",") if x}
TOKEN = sys.stdin.read().strip()
assert len(TOKEN) >= 20, "token on stdin"
for d in ("skani_results", "skani_hits"):
    os.makedirs(os.path.join(WORK, d), exist_ok=True)
if not a.foreground:
    if os.fork():
        sys.exit(0)
    os.setsid()
    sys.stdout = sys.stderr = open(os.path.join(WORK, "skani_driver.log"), "a", buffering=1)
LOCK = open(os.path.join(WORK, "skani_sweep.lock"), "a+")
try:
    fcntl.flock(LOCK, fcntl.LOCK_EX | fcntl.LOCK_NB)
except BlockingIOError:
    sys.exit("another skani sweep holds the lock")
LOCK.seek(0); LOCK.truncate(); LOCK.write(str(os.getpid())); LOCK.flush()
S = requests.Session()
S.headers["Authorization"] = "Bearer " + TOKEN


def ev(msg):
    line = f"[{time.strftime('%F %T')}] {msg}"
    print(line, flush=True)
    with open(os.path.join(WORK, "skani_events.log"), "a") as fh:
        fh.write(line + "\n")


db = sqlite3.connect(os.path.join(WORK, "skani_ledger.sqlite"), timeout=60)
db.executescript("""
CREATE TABLE IF NOT EXISTS batches(id INTEGER PRIMARY KEY AUTOINCREMENT, archives TEXT, job_id TEXT, state TEXT,
  created REAL, attempt INT DEFAULT 0, note TEXT);
CREATE TABLE IF NOT EXISTS archives(name TEXT PRIMARY KEY, object_id TEXT, batch_id INT, state TEXT);
CREATE TABLE IF NOT EXISTS kv(k TEXT PRIMARY KEY, v TEXT);
CREATE TABLE IF NOT EXISTS failed_mags(mag TEXT PRIMARY KEY, archive TEXT, batch_id INT, reason TEXT);
""")
if "subset" not in [r[1] for r in db.execute("PRAGMA table_info(batches)")]:
    db.execute("ALTER TABLE batches ADD COLUMN subset TEXT")   # {archive: [mag]} for a retry of errored MAGs
    db.commit()


def halted():
    r = db.execute("SELECT v FROM kv WHERE k='halted'").fetchone()
    return bool(r and r[0] == "1")


def halt(msg):
    if not halted():
        ev(f"HALT {msg}")
    db.execute("INSERT OR REPLACE INTO kv VALUES('halted','1')")
    db.execute("INSERT OR REPLACE INTO kv VALUES('halt_reason',?)", (msg[:500],))
    db.commit()


def build_view():
    """({archive: object_id} for uploaded archives, {archive: [mag]}, is every archive uploaded?)"""
    b = sqlite3.connect(f"file:{BUILD_LEDGER}?mode=ro", uri=True, timeout=60)
    try:
        arch = dict(b.execute("SELECT name, object_id FROM archives WHERE object_id IS NOT NULL").fetchall())
        mem = {}
        for mag, name in b.execute("SELECT mag, archive FROM members"):
            mem.setdefault(name, []).append(mag)
        left = b.execute("SELECT count(*) FROM archives WHERE object_id IS NULL").fetchone()[0]
    finally:
        b.close()
    if ONLY:
        missing = ONLY - set(arch)
        assert not missing, f"--only archives not uploaded: {sorted(missing)}"
        arch = {n: o for n, o in arch.items() if n in ONLY}
        return arch, mem, True
    return arch, mem, left == 0


def dbid(x):
    return x.get("object_id") if isinstance(x, dict) else x


def submit(bid, names, oids, attempt, subset=None):
    pl = {"schema_version": "1", "job_type": "KBDLSKANI",
          "params": {"skani_dbs": [{"object_id": d} for d in DBS],
                     "archives": [{"object_id": o, **({"genome_ids": subset[n]} if subset else {})}
                                  for n, o in zip(names, oids)]}}
    assert len(oids) == len(names) and all(oids)
    assert sum(len(subset[n]) for n in names) <= MAX_QUERIES if subset else True
    db.execute("UPDATE batches SET state='unknown', created=?, attempt=? WHERE id=?", (time.time(), attempt, bid))
    db.commit()
    try:
        r = S.post(BASE + "/jobs", json=pl, timeout=120)   # never retried blindly: 'unknown' is reconciled
    except requests.RequestException as e:
        ev(f"WARN POST lost for batch {bid} ({e!s:.100}); will reconcile")
        return
    if r.status_code != 202:
        if r.status_code >= 500:
            ev(f"WARN POST batch {bid} -> {r.status_code}; will reconcile")
            return
        db.execute("UPDATE batches SET state='rejected', note=? WHERE id=?", (r.text[:300], bid))
        db.commit()
        return halt(f"POST batch {bid} -> {r.status_code} {r.text[:200]}")
    jid = r.json()["job_id"]
    db.execute("UPDATE batches SET job_id=?, state='queued' WHERE id=?", (jid, bid))
    db.commit()
    ev(f"SUBMIT batch {bid} ({len(names)} archives: {names[0]} .. {names[-1]}"
       + (f"; subset of {sum(map(len, subset.values()))} MAGs" if subset else "") + f") -> {jid} (attempt {attempt})")


def record_failed(mags_by_archive, bid, reason):
    for n, mags in mags_by_archive.items():
        for m in mags:
            db.execute("INSERT OR REPLACE INTO failed_mags VALUES(?,?,?,?)", (m, n, bid, reason[:500]))
    db.commit()
    nf = db.execute("SELECT count(*) FROM failed_mags").fetchone()[0]
    ev(f"FAILED-MAGS batch {bid}: {sum(map(len, mags_by_archive.values()))} recorded ({nf} in all): {reason[:200]}")
    if nf > MAX_FAILED_MAGS:
        halt(f"{nf} MAGs failed (> {MAX_FAILED_MAGS})")


def reconcile():
    pend = db.execute("SELECT id, created FROM batches WHERE state='unknown' ORDER BY created").fetchall()
    if not pend:
        return 0
    known = {r[0] for r in db.execute("SELECT job_id FROM batches WHERE job_id IS NOT NULL")}
    lst = S.get(BASE + "/jobs", timeout=120).json()

    def ts(j):
        try:
            return calendar.timegm(time.strptime(j["submitted_at"][:19], "%Y-%m-%dT%H:%M:%S"))
        except Exception:
            return 0
    orph = sorted((j for j in lst if j.get("job_type") == "KBDLSKANI" and j["job_id"] not in known
                   and ts(j) >= pend[0][1] - 5), key=ts)
    if orph and len(orph) != len(pend):
        halt(f"RECONCILE ambiguous: orphan KBDLSKANI jobs {[j['job_id'] for j in orph]} vs batches {[x[0] for x in pend]}")
        return len(pend)
    unresolved = 0
    for i, (bid, created) in enumerate(pend):
        if orph:
            db.execute("UPDATE batches SET job_id=?, state='queued' WHERE id=?", (orph[i]["job_id"], bid))
            ev(f"RECONCILE adopted {orph[i]['job_id']} for batch {bid}")
        elif time.time() - created > 900:
            db.execute("UPDATE batches SET state='retry_same' WHERE id=?", (bid,))
            ev(f"RECONCILE batch {bid}: no job was created; will resubmit")
        else:
            unresolved += 1
    db.commit()
    return unresolved


def harvest(bid, jid, names, mem, subset):
    r = S.get(f"{BASE}/jobs/{jid}/result", timeout=600)
    if r.status_code != 200:
        ev(f"WARN result {jid} -> {r.status_code}; re-poll")
        return
    res = r.json()
    with open(os.path.join(WORK, "skani_results", f"{jid}.json"), "w") as fh:
        json.dump(res, fh)
    queried = sorted(dbid(x) for x in (res.get("databases_queried") or []))
    if queried != sorted(DBS):
        return halt(f"batch {bid} job {jid}: databases_queried {queried} != both DBs")
    owner = {m: n for n in names for m in (subset[n] if subset else mem[n])}
    if len(owner) > MAX_QUERIES:
        return halt(f"batch {bid} has {len(owner)} MAGs > {MAX_QUERIES} (indexed skani); not harvesting it")
    results = res.get("results") or {}
    errs = res.get("errors") or {}
    stray = set(results) - set(owner)
    bad_db = {h.get("database") for hs in results.values() for h in hs} - set(DBS)
    if stray or bad_db:
        return halt(f"batch {bid} job {jid}: unexpected query keys {sorted(stray)[:3]} / databases {bad_db}")
    oid2name = {db.execute("SELECT object_id FROM archives WHERE name=?", (n,)).fetchone()[0]: n for n in names}
    bad = set()
    for k in errs:
        if k.startswith("archive::"):
            bad |= {m for m, n in owner.items() if n == oid2name.get(k.split("::", 1)[1])}
        else:
            bad.add(k.rsplit("::", 1)[0])
    unknown = bad - set(owner)
    if unknown:
        return halt(f"batch {bid} job {jid}: error keys name non-members {sorted(unknown)[:3]}")
    if errs:
        ev(f"BATCH-ERRORS batch {bid} job {jid}: {len(errs)} error entries on {len(bad)} MAGs, "
           f"e.g. {json.dumps(dict(list(errs.items())[:2]))[:400]}")
        by_arch = {}
        for m in sorted(bad):
            by_arch.setdefault(owner[m], []).append(m)
        if subset:      # already a retry: final
            record_failed(by_arch, bid, "per-entry errors on the subset retry: " + json.dumps(dict(list(errs.items())[:1])))
        else:
            cur = db.execute("INSERT INTO batches(archives, state, created, subset, note) VALUES(?,?,?,?,?)",
                             (json.dumps(sorted(by_arch)), "pending", time.time(), json.dumps(by_arch),
                              f"retry of errored MAGs from batch {bid}"))
            ev(f"RETRY-SUBSET batch {cur.lastrowid} for {len(bad)} MAGs of batch {bid}")
    owner = {m: n for m, n in owner.items() if m not in bad}
    rows, nohit = [], {d: 0 for d in DBS}
    for mag in sorted(owner):
        hits = results.get(mag) or []
        for d in DBS:
            mine = [h for h in hits if h.get("database") == d]
            for h in mine:
                rows.append(f"{mag}\t{owner[mag]}\t{DBNAME[d]}\t{h['genome_id']}\t{1.0 - float(h['distance']):.6f}\n")
            if not mine:
                rows.append(f"{mag}\t{owner[mag]}\t{DBNAME[d]}\t\t\n")
                nohit[d] += 1
    tmp = os.path.join(WORK, "skani_hits", f".batch_{bid}.tsv.tmp")
    with open(tmp, "w") as fh:
        fh.write("mag\tarchive\tdb\tref_genome\tani\n")
        fh.writelines(rows)
    os.replace(tmp, os.path.join(WORK, "skani_hits", f"batch_{bid}.tsv"))
    db.execute("UPDATE batches SET state='done' WHERE id=?", (bid,))
    if not subset:
        for n in names:
            db.execute("UPDATE archives SET state='done' WHERE name=?", (n,))
    db.commit()
    nh = sum(len(results.get(m) or []) for m in owner)
    ev(f"HARVEST batch {bid} job {jid}: {len(owner)} MAGs, {nh} hits; no hit: "
       + ", ".join(f"{DBNAME[d]} {nohit[d]}" for d in DBS)
       + f"; queries absent from results: {len(set(owner) - set(results))}" + (f"; held out {len(bad)}" if bad else ""))


ev(f"START skani sweep pid {os.getpid()} dbs={[DBNAME[d] for d in DBS]} max_queries={MAX_QUERIES} "
   f"max_in_flight={a.max_in_flight} only={sorted(ONLY) or 'all'}")
while True:
    try:
        unresolved = reconcile()
        arch, mem, build_finished = build_view()
        for n, o in arch.items():
            db.execute("INSERT OR IGNORE INTO archives(name, object_id, state) VALUES(?,?,'new')", (n, o))
        db.commit()
        act = db.execute("SELECT id, job_id, archives, subset FROM batches WHERE state IN ('queued','running')").fetchall()
        if act:
            st = S.post(BASE + "/jobs/status", json={"job_ids": [x[1] for x in act]}, timeout=120).json().get("jobs", {})
            for bid, jid, names, sub in act:
                s = (st.get(jid) or {}).get("state")
                names, sub = json.loads(names), (json.loads(sub) if sub else None)
                if s == "completed":
                    harvest(bid, jid, names, mem, sub)
                elif s in ("failed", "cancelled"):
                    att = db.execute("SELECT attempt FROM batches WHERE id=?", (bid,)).fetchone()[0]
                    fail = json.dumps((st.get(jid) or {}).get("failure"))[:300]
                    ev(f"JOBFAIL batch {bid} job {jid} ({s}, attempt {att}): {fail}")
                    db.execute("UPDATE batches SET state=? WHERE id=?", ("retry" if att < 2 else "failed", bid))
                    db.commit()
                    if att >= 2:
                        record_failed(sub or {n: mem[n] for n in names}, bid, f"job failed twice: {fail}")
                        nb = db.execute("SELECT count(*) FROM batches WHERE state='failed'").fetchone()[0]
                        if nb >= MAX_FAILED_BATCHES:
                            halt(f"{nb} batches failed twice")
                elif s in ("queued", "running"):
                    db.execute("UPDATE batches SET state=? WHERE id=?", (s, bid))
            db.commit()
        # new batches from fresh archives, up to MAX_QUERIES MAGs each: a group is closed when it is full or the next
        # archive would overflow it; the last open group goes out once the build (or the --only set) is complete
        fresh = [r[0] for r in db.execute("SELECT name FROM archives WHERE state='new' ORDER BY name")]
        groups, cur_g, nq = [], [], 0
        for n in fresh:
            k = len(mem[n])
            assert 0 < k <= MAX_QUERIES, (n, k)
            if nq + k > MAX_QUERIES:
                groups.append(cur_g); cur_g, nq = [], 0
            cur_g.append(n); nq += k
            if nq == MAX_QUERIES:
                groups.append(cur_g); cur_g, nq = [], 0
        if cur_g and build_finished:
            groups.append(cur_g)
        for take in groups:
            cur = db.execute("INSERT INTO batches(archives, state, created) VALUES(?,?,?)",
                             (json.dumps(take), "pending", time.time()))
            for n in take:
                db.execute("UPDATE archives SET state='batched', batch_id=? WHERE name=?", (cur.lastrowid, n))
            db.commit()
        inflight = unresolved + db.execute("SELECT count(*) FROM batches WHERE state IN ('queued','running')").fetchone()[0]
        for bid, names, att, state, sub in db.execute("SELECT id, archives, attempt, state, subset FROM batches WHERE "
                                                      "state IN ('pending','retry','retry_same') ORDER BY id").fetchall():
            if halted() or unresolved or inflight >= a.max_in_flight:
                break
            names = json.loads(names)
            if ONLY and not set(names) <= ONLY:
                continue
            oids = [db.execute("SELECT object_id FROM archives WHERE name=?", (n,)).fetchone()[0] for n in names]
            submit(bid, names, oids, att if state == "retry_same" else att + 1, json.loads(sub) if sub else None)
            inflight += 1
        open_jobs = db.execute("SELECT count(*) FROM batches WHERE state IN ('queued','running','unknown')").fetchone()[0]
        todo = [json.loads(r[0]) for r in
                db.execute("SELECT archives FROM batches WHERE state IN ('pending','retry','retry_same')")]
        if ONLY:
            todo = [t for t in todo if set(t) <= ONLY]
        new = db.execute("SELECT count(*) FROM archives WHERE state='new'").fetchone()[0]
        if not open_jobs and (halted() or (build_finished and not todo and not new)):
            ev("STOPPED (halted: see kv halt_reason)" if halted() else
               f"DONE {'--only set' if ONLY else 'all archives'} swept against both databases")
            break
    except Exception as e:
        ev(f"WARN loop error: {e!r:.300}")
    time.sleep(POLL_S)
