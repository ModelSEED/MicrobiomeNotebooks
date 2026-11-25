#!/usr/bin/env python3
"""TRANSYT (KBDLGenomeAnnotation, tool TRANSYT_LOCAL) on the 1,206 DRAM2/Bakta representatives (per Andrew,
2026-11-25). KBDL runs TRANSYT natively, so no KBDL change is needed.

One job per genome. Inputs:
  - genome_ref = the representative's verified KBDL Genome, from the build ledger, read only;
  - tax_id = rep_taxids.tsv (make_taxids.py: the NCBI taxid of the GTDB R232 species representative of the MAG's
    deepest named rank; the KBDL genomes carry GTDB taxonomy without an NCBI id).
Representatives with no verified KBDL genome (Bakta edge-interval failures) are listed and skipped.
Result per MAG -> WORK/transyt_results/<mag>.json ({gene_id: {"TRANSYT": {TC, MSRXN, MSCPD, TRANSYT_RXN}}}).
Ledger WORK/transyt_ledger.sqlite; events WORK/transyt_events.log. A failed job is retried once. A POST whose outcome
is unknown is never re-sent blindly: it is matched to a new job by time first.
Runs ON POPLAR against http://localhost:8791. Token on stdin; detaches unless --foreground.
Usage: python3 -I transyt_sweep.py [--max-in-flight 8] < token
"""
import argparse, calendar, csv, fcntl, json, os, sqlite3, sys, time
import requests

BASE = "http://localhost:8791"
WORK = "/scratch/afreiburger/wetland_kbdl_build"
HERE = os.path.dirname(os.path.abspath(__file__))
p = argparse.ArgumentParser(); p.add_argument("--max-in-flight", type=int, default=8); p.add_argument("--foreground", action="store_true")
a = p.parse_args()
TOKEN = sys.stdin.read().strip(); assert len(TOKEN) >= 20
os.makedirs(os.path.join(WORK, "transyt_results"), exist_ok=True)
if not a.foreground:
    if os.fork(): sys.exit(0)
    os.setsid(); sys.stdout = sys.stderr = open(os.path.join(WORK, "transyt_driver.log"), "a", buffering=1)
LOCK = open(os.path.join(WORK, "transyt.lock"), "a+")
try: fcntl.flock(LOCK, fcntl.LOCK_EX | fcntl.LOCK_NB)
except BlockingIOError: sys.exit("another TRANSYT sweep holds the lock")
LOCK.seek(0); LOCK.truncate(); LOCK.write(str(os.getpid())); LOCK.flush()
S = requests.Session(); S.headers["Authorization"] = "Bearer " + TOKEN
def ev(m):
    line = f"[{time.strftime('%F %T')}] {m}"; print(line, flush=True)
    with open(os.path.join(WORK, "transyt_events.log"), "a") as fh: fh.write(line + "\n")
db = sqlite3.connect(os.path.join(WORK, "transyt_ledger.sqlite"), timeout=60)
db.execute("CREATE TABLE IF NOT EXISTS t(mag TEXT PRIMARY KEY, genome TEXT, tax_id TEXT, job_id TEXT, state TEXT, attempt INT DEFAULT 0, created REAL, note TEXT)")
taxa = {r["mag"]: r["tax_id"] for r in csv.DictReader(open(os.path.join(HERE, "rep_taxids.tsv")), delimiter="\t")}
b = sqlite3.connect(f"file:{WORK}/ledger.sqlite?mode=ro", uri=True)
gen = {m: g for m, v, g in b.execute("SELECT mag, verified, genome_object_id FROM members") if v == 1 and g}
b.close()
skipped = sorted(m for m in taxa if m not in gen)
for m, t in taxa.items():
    if m in gen:
        db.execute("INSERT OR IGNORE INTO t(mag, genome, tax_id, state) VALUES(?,?,?, 'pending')", (m, gen[m], t))
db.commit()
ev(f"START TRANSYT sweep pid {os.getpid()}: {len(taxa)} reps, {len(gen.keys() & taxa.keys())} with a KBDL genome; "
   f"skipped (no verified KBDL genome): {skipped}")

def submit(mag, genome, tax, att):
    pl = {"schema_version": "1", "job_type": "KBDLGenomeAnnotation",
          "params": {"genome_ref": {"object_id": genome}, "tools": ["TRANSYT_LOCAL"], "tax_id": tax}}
    db.execute("UPDATE t SET state='unknown', created=?, attempt=? WHERE mag=?", (time.time(), att, mag)); db.commit()
    try: r = S.post(BASE + "/jobs", json=pl, timeout=120)
    except requests.RequestException as e: return ev(f"WARN POST lost for {mag}: {e!s:.80}")
    if r.status_code != 202:
        db.execute("UPDATE t SET state='rejected', note=? WHERE mag=?", (r.text[:300], mag)); db.commit()
        return ev(f"REJECTED {mag}: {r.status_code} {r.text[:200]}")
    db.execute("UPDATE t SET job_id=?, state='queued' WHERE mag=?", (r.json()["job_id"], mag)); db.commit()

def reconcile():
    pend = db.execute("SELECT mag, created FROM t WHERE state='unknown'").fetchall()
    if not pend: return 0
    known = {r[0] for r in db.execute("SELECT job_id FROM t WHERE job_id IS NOT NULL")}
    ts = lambda j: calendar.timegm(time.strptime(j["submitted_at"][:19], "%Y-%m-%dT%H:%M:%S"))
    orph = [j for j in S.get(BASE + "/jobs", timeout=120).json() if j.get("job_type") == "KBDLGenomeAnnotation"
            and j["job_id"] not in known and ts(j) >= min(c for _, c in pend) - 5]
    if len(pend) == 1 and len(orph) == 1:
        db.execute("UPDATE t SET job_id=?, state='queued' WHERE mag=?", (orph[0]["job_id"], pend[0][0])); ev(f"RECONCILE adopted for {pend[0][0]}")
    elif not orph:
        for m, c in pend:
            if time.time() - c > 900: db.execute("UPDATE t SET state='pending' WHERE mag=?", (m,))
    else:
        ev(f"WARN {len(orph)} unmatched jobs for {len(pend)} unknown POSTs; waiting"); return len(pend)
    db.commit(); return 0

while True:
    try:
        unresolved = reconcile()
        act = db.execute("SELECT mag, job_id, attempt FROM t WHERE state IN ('queued','running')").fetchall()
        if act:
            st = S.post(BASE + "/jobs/status", json={"job_ids": [x[1] for x in act]}, timeout=120).json().get("jobs", {})
            for mag, jid, att in act:
                s = (st.get(jid) or {}).get("state")
                if s == "completed":
                    r = S.get(f"{BASE}/jobs/{jid}/result", timeout=300)
                    if r.status_code == 200:
                        res = r.json(); n = sum(1 for v in res.values() if (v.get("TRANSYT") or {}))
                        json.dump(res, open(os.path.join(WORK, "transyt_results", f"{mag}.json"), "w"))
                        db.execute("UPDATE t SET state='done', note=? WHERE mag=?", (f"{n} of {len(res)} genes", mag))
                elif s in ("failed", "cancelled"):
                    ev(f"JOBFAIL {mag} ({s}, attempt {att}): {json.dumps(st[jid].get('failure'))[:200]}")
                    db.execute("UPDATE t SET state=? WHERE mag=?", ("pending" if att < 2 else "failed", mag))
                elif s in ("queued", "running"):
                    db.execute("UPDATE t SET state=? WHERE mag=?", (s, mag))
            db.commit()
        inflight = unresolved + db.execute("SELECT count(*) FROM t WHERE state IN ('queued','running','unknown')").fetchone()[0]
        for mag, g, tax, att in db.execute("SELECT mag, genome, tax_id, attempt FROM t WHERE state='pending' ORDER BY mag").fetchall():
            if inflight >= a.max_in_flight or unresolved: break
            submit(mag, g, tax, att + 1); inflight += 1
        c = dict(db.execute("SELECT state, count(*) FROM t GROUP BY state").fetchall())
        if not any(c.get(k) for k in ("pending", "queued", "running", "unknown")):
            ev(f"DONE {c}"); break
    except Exception as e:
        ev(f"WARN loop error: {e!r:.200}")
    time.sleep(60)
