"""Import the GROWdb raw metagenomes (runs.tsv: 171 SRA runs = the 133 samples behind the 3,825 MAGs of KBase
ws 145226) into KBase narrative 275929 as interleaved PairedEndLibraries.

App: kb_uploadmethods/import_sra_as_reads_from_web (beta), the same app/params the salterns reassembly used
(salterns/reassembly_2026/pipeline.py). KBase fetches each run itself from the AWS ODP mirror of SRA
(runs.tsv:odp_url), so nothing is uploaded from this machine.

One library per SRA run, named <sample>.<SRR>.reads. A library is accepted only if read_count == 2 x SRA spots
and total_bases == SRA bases. Multi-run samples are NOT merged (distinct sequencing libraries).

Usage (KBUtilLib from ~/Documents/KBUtilLib/src; token from ~/.kbase/token, never printed):
    python kbase_import.py status
    python kbase_import.py plan                       # what advance would submit now; submits nothing
    python kbase_import.py advance                    # poll, validate, submit up to the cap
    python kbase_import.py watch --interval 900       # advance in a loop until every run is done/failed
The largest run is a pilot: nothing else is submitted until it has completed and validated.
"""
import argparse, csv, fcntl, json, os, sys, time
from collections import Counter
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.expanduser("~/Documents/KBUtilLib/src"))
from kbutillib.kb_app_runner.nms import NMSSpecCache  # noqa: E402
from kbutillib.kb_app_runner.runner import AppRunner  # noqa: E402
from kbutillib.kb_genome_utils import KBGenomeUtils  # noqa: E402
from kbutillib.kb_job_utils import KBJobUtils  # noqa: E402

WS_ID = 275929
WS_NAME = "freiburgermsu:narrative_1794152940209"
APP, TAG = "kb_uploadmethods/import_sra_as_reads_from_web", "beta"
RUNS = os.path.join(HERE, "runs.tsv")
LEDGER = os.path.join(HERE, "kbase_import_ledger.json")
CAP = 10                      # live njs jobs from this driver (salterns watcher has its own njs cap of 15)
MAX_TRIES = 3                 # counted (non-cancelled) attempts per run
LIVE = {"created", "estimating", "queued", "running"}
ORG = "grow"                  # KBase Organization "Genome Resolved Open Watersheds database" (user is a Member)
GROUPS = "https://kbase.us/services/groups"
WS_DESCRIPTION = ("Raw metagenomes (NCBI SRA) of the 133 GROWdb samples the 3,825 GROWdb MAGs (KBase ws 145226) were "
                  "binned from: 65 NASQAN + 2 Yukon (BioProject PRJNA946291) and 66 WHONDRS-2019 surface-water samples "
                  "(JGI per-sample BioProjects). One interleaved PairedEndLibrary per SRA run, <sample>.<SRR>.reads, "
                  "imported from the AWS ODP mirror; read and base counts verified against SRA.")


def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def log(msg):
    print(f"{now()} {msg}", flush=True)


def retry_rpc(fn, *a, tries=4, **kw):
    for k in range(tries):
        try:
            return fn(*a, **kw)
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(10 * (k + 1))


def load_runs():
    rows = list(csv.DictReader(open(RUNS), delimiter="\t"))
    for r in rows:
        r["name"] = f"{r['sample']}.{r['Run']}.reads"
    return rows


def pilot_run(rows):
    return max(rows, key=lambda r: int(r["bases"]))["Run"]


def load_ledger():
    if os.path.exists(LEDGER):
        return json.load(open(LEDGER))
    return {"workspace_id": WS_ID, "workspace_name": WS_NAME, "app": APP, "tag": TAG, "created_utc": now(), "runs": {}}


def save_ledger(L):
    tmp = LEDGER + ".tmp"
    json.dump(L, open(tmp, "w"), indent=1)
    os.replace(tmp, LEDGER)


def params_for(r):
    return [{"workspace_name": WS_NAME, "download_type": "Direct Download",
             "sra_urls_to_add": [{"file_url": r["odp_url"], "name": r["name"], "sequencing_tech": "Illumina",
                                  "single_genome": 0}]}]


class KB:
    def __init__(self, submit=True):
        self.u = KBGenomeUtils(kb_version="prod")
        self.jobs = KBJobUtils(self.u, kb_version="prod")
        self.ws = self.u.ws_client()
        self.ee2 = self.jobs.ee2_client
        self.nms = NMSSpecCache()
        self.runner = AppRunner(ws=self.u, kb_version="prod", job_store=self.jobs, nms_cache=self.nms) if submit else None

    def states(self, ids):
        out = {}
        for i in range(0, len(ids), 100):
            res = retry_rpc(self.ee2.check_jobs, {"job_ids": ids[i:i + 100], "exclude_fields": ["job_input"]})
            for s in (res.get("job_states", res) if isinstance(res, dict) else res):
                if isinstance(s, dict) and s.get("job_id"):
                    out[s["job_id"]] = s
        return out

    def info(self, name):
        infos = retry_rpc(self.ws.get_object_info3, {"objects": [{"wsid": WS_ID, "name": name}], "ignoreErrors": 1,
                                                    "includeMetadata": 1})["infos"]
        return infos[0] if infos else None

    def find_job(self, name, since_utc, exclude=()):
        """EE2 job of this app in WS_ID whose params name `name`, created after since_utc - 2 min (a lost submit)."""
        since_ms = int(datetime.fromisoformat(since_utc).timestamp() * 1000) - 120_000
        res = retry_rpc(self.ee2.check_workspace_jobs, {"workspace_id": WS_ID, "return_list": 1,
                                                        "exclude_fields": ["job_output"]})
        for j in sorted(res.get("job_states", []), key=lambda j: j.get("created") or 0, reverse=True):
            ji = j.get("job_input") or {}
            if (ji.get("app_id") == APP and (j.get("created") or 0) >= since_ms and j.get("job_id") not in exclude
                    and j.get("status") in LIVE | {"completed"} and json.dumps(name) in json.dumps(ji.get("params"))):
                return j
        return None


def validate(info, r):
    """(ok, reason) for an imported library against SRA's spot/base counts."""
    if not info:
        return False, "no output object"
    if not info[2].startswith("KBaseFile.PairedEndLibrary"):
        return False, f"output is {info[2]}, expected a PairedEndLibrary"
    m = info[10] or {}
    got = (int(m.get("read_count") or -1), int(m.get("total_bases") or -1))
    exp = (2 * int(r["spots"]), int(r["bases"]))
    return (True, None) if got == exp else (False, f"read_count/total_bases {got} != SRA {exp}")


def err_text(s):
    e = s.get("error")
    return str((e.get("message") if isinstance(e, dict) else None) or s.get("errormsg") or "")[:2000]


def counted(a):
    return a.get("status") not in ("terminated",) and not a.get("transient")


def poll(kb, L, runs):
    live = {a["job_id"]: (srr, a) for srr, t in L["runs"].items() for a in t["attempts"]
            if a.get("status") in LIVE and a.get("job_id")}
    st = kb.states(list(live)) if live else {}
    for jid, (srr, a) in live.items():
        s = st.get(jid)
        if not s:
            continue
        t, status = L["runs"][srr], s["status"]
        if status in LIVE:
            a["status"] = status
            continue
        a["status"], a["finished_utc"] = status, now()
        if status == "completed":
            ok, why = validate(kb.info(runs[srr]["name"]), runs[srr])
            if ok:
                i = kb.info(runs[srr]["name"])
                t.update(state="done", ref=f"{i[6]}/{i[0]}/{i[4]}")
                log(f"{srr}: done -> {t['ref']} ({runs[srr]['name']})")
            else:
                a.update(status="error", error=why, rejected=True)
                t["state"] = "retry"
                log(f"{srr}: rejected: {why}")
        else:
            a["error"] = err_text(s)
            transient = any(k in a["error"] for k in ("timed out", "Connection", "503", "502", "504", "Temporary"))
            a["transient"] = transient and status == "error"
            t["state"] = "retry"
            log(f"{srr}: {status}: {a['error'][:300]}")


def advance(kb, L, rows, dry_run=False):
    runs = {r["Run"]: r for r in rows}
    for r in rows:
        L["runs"].setdefault(r["Run"], {"sample": r["sample"], "name": r["name"], "state": "todo", "attempts": []})
    poll(kb, L, runs)
    save_ledger(L)
    pilot = pilot_run(rows)
    fan_out = L["runs"][pilot]["state"] == "done"
    live = sum(1 for t in L["runs"].values() for a in t["attempts"] if a.get("status") in LIVE)
    order = [pilot] + [r["Run"] for r in sorted(rows, key=lambda r: (r["source"], r["sample"], r["Run"])) if r["Run"] != pilot]
    planned = []
    for srr in order:
        t, r = L["runs"][srr], runs[srr]
        if t["state"] in ("done", "failed") or any(a.get("status") in LIVE for a in t["attempts"]):
            if not fan_out:
                break                                   # pilot in flight: nothing else yet
            continue
        if sum(map(counted, t["attempts"])) >= MAX_TRIES:
            t["state"] = "failed"
            log(f"{srr}: failed after {MAX_TRIES} attempts")
            continue
        if live >= CAP:
            break
        info = kb.info(r["name"])                       # an object already saved (lost ledger write): adopt it
        if info:
            ok, why = validate(info, r)
            if ok:
                t.update(state="done", ref=f"{info[6]}/{info[0]}/{info[4]}", adopted=True)
                log(f"{srr}: adopted existing {t['ref']}")
                continue
        if t.get("submitting"):                         # a submit that raised or was interrupted: did EE2 take it?
            j = kb.find_job(r["name"], t["submitting"], exclude={a.get("job_id") for a in t["attempts"]})
            t.pop("submitting")
            if j:
                t["attempts"].append({"job_id": j["job_id"], "submitted_utc": now(), "status": j["status"],
                                      "reconciled": True})
                t["state"] = "active"
                live += 1
                save_ledger(L)
                log(f"{srr}: adopted orphan job {j['job_id']}")
                continue
        planned.append(srr)
        live += 1
        if dry_run:
            if not fan_out:
                break
            continue
        t["submitting"] = now()
        save_ledger(L)                                  # write-ahead
        try:
            jid = kb.runner.run_app(APP, params_for(r), WS_ID, tag=TAG, meta={"stage": "import", "unit": srr}).job_id
        except Exception as e:
            log(f"{srr}: submit raised {type(e).__name__}: {str(e)[:300]} (reconciled next cycle)")
            save_ledger(L)
            continue
        t.pop("submitting", None)
        t["attempts"].append({"job_id": jid, "submitted_utc": now(), "status": "queued", "params": params_for(r)})
        t["state"] = "active"
        save_ledger(L)
        log(f"{srr}: submitted {jid} ({r['name']}, {int(r['bases'])/1e9:.1f} Gbp){' [pilot]' if srr == pilot else ''}")
        if not fan_out:
            break
    save_ledger(L)
    return planned


def org_add(kb, L, rows):
    """Once every run is imported and validated: describe the narrative and ask to attach it to the grow Org.
    As an Org Member (not admin) this files a request the Org admins approve. Returns True when settled."""
    o = L.setdefault("org_add", {})
    if o.get("status") in ("requested", "added"):
        return True
    if not (len(L["runs"]) == len(rows) and all(t["state"] == "done" for t in L["runs"].values())):
        return False
    import requests
    tok = open(os.path.expanduser("~/.kbase/token")).read().strip()
    try:
        retry_rpc(kb.ws.set_workspace_description, {"id": WS_ID, "description": WS_DESCRIPTION})
        r = requests.post(f"{GROUPS}/group/{ORG}/resource/workspace/{WS_ID}", headers={"Authorization": tok}, timeout=120)
        body = r.json() if r.headers.get("content-type", "").startswith("application/json") else {"text": r.text[:500]}
        if not r.ok:
            raise RuntimeError(f"HTTP {r.status_code}: {json.dumps(body)[:500]}")
    except Exception as e:
        o.update(status="error", error=str(e)[:600], utc=now())
        log(f"org add to '{ORG}' failed ({str(e)[:300]}); retried next cycle")
        return False
    o.update(status="added" if body.get("complete") else "requested", request_id=body.get("id"), response=body, utc=now())
    log(f"org '{ORG}': {'narrative attached' if o['status'] == 'added' else 'add request filed for admin approval'} "
        f"(request {o.get('request_id')})")
    return True


def status(L, rows):
    c = Counter(t["state"] for t in L["runs"].values())
    c["todo"] += len(rows) - len(L["runs"])
    live = Counter(a["status"] for t in L["runs"].values() for a in t["attempts"] if a.get("status") in LIVE)
    done = [r for r in rows if L["runs"].get(r["Run"], {}).get("state") == "done"]
    print(f"narrative {WS_ID} ({WS_NAME}); {len(rows)} runs / {len({r['sample'] for r in rows})} samples")
    print(f"states: {dict(c)}   live jobs: {dict(live)}")
    print(f"imported: {len(done)} runs, {sum(int(r['bases']) for r in done)/1e12:.3f} of "
          f"{sum(int(r['bases']) for r in rows)/1e12:.3f} Tbp; pilot {pilot_run(rows)} -> "
          f"{L['runs'].get(pilot_run(rows), {}).get('state', 'todo')}")
    o = L.get("org_add") or {}
    print(f"org '{ORG}' add: {o.get('status', 'pending (runs once all runs are done)')}"
          f"{' request ' + str(o.get('request_id')) if o.get('request_id') else ''}")
    for srr, t in L["runs"].items():
        if t["state"] in ("retry", "failed"):
            print(f"  {t['state']:6} {srr} {t['name']}: {(t['attempts'][-1].get('error') or '')[:200] if t['attempts'] else ''}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["status", "plan", "advance", "watch"])
    ap.add_argument("--interval", type=int, default=900)
    a = ap.parse_args()
    rows = load_runs()
    if a.cmd == "status":
        status(load_ledger(), rows)
        return
    lock = open(LEDGER + ".lock", "w")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)    # one writer at a time
    kb = KB(submit=a.cmd != "plan")
    while True:
        L = load_ledger()
        planned = advance(kb, L, rows, dry_run=a.cmd == "plan")
        if a.cmd == "plan":
            print("would submit:", planned)
            return
        if a.cmd == "advance":
            status(L, rows)
            return
        if all(t["state"] in ("done", "failed") for t in L["runs"].values()) and len(L["runs"]) == len(rows):
            if all(t["state"] == "done" for t in L["runs"].values()):
                settled = org_add(kb, L, rows)
                save_ledger(L)
                if not settled:                         # transient org/ws error: try again next cycle
                    time.sleep(a.interval)
                    continue
            else:
                log(f"{sum(t['state'] == 'failed' for t in L['runs'].values())} runs failed: NOT adding the narrative "
                    f"to org '{ORG}' (needs a decision)")
            status(L, rows)
            log("all runs settled; watcher exiting")
            return
        time.sleep(a.interval)


if __name__ == "__main__":
    main()
