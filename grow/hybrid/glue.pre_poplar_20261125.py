"""Hybrid local/KBase glue shared by the two GROW drivers (user decision 2026-11-23, both BRIEFs).

Each driver passes itself as P (its module: task, get_task, output_of, log, now, save_ledger, validate_output,
on_completed, WORK, WS_ID, WS_NAME, LIVE, ...). For the stages in P.HYBRID the driver no longer submits KBase jobs:
  * a pending/retry task is run locally by hw.py (one process per step, detached, result.json at the end);
  * an active task whose KBase job is still QUEUED is taken over when a local slot is free at that moment: the job
    is re-checked, cancelled (attempt status 'terminated', takeover=True: not counted against the unit), and the
    step starts locally in the same cycle. Jobs already RUNNING in KBase are left to finish and adopted as before;
  * the first local step of every stage is a pilot: one at a time until one has finished 'done' (hybrid_pilot_ok).
Local attempts are recorded in the ledger as {"local": True, "host": ..., "workdir", "pid", "versions", ...}.
"""
import hashlib
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gate as G  # noqa: E402

HW = os.path.join(HERE, "hw.py")
PY = os.path.expanduser("~/Documents/py_venv/bin/python")
READS_CACHE = os.path.join(HERE, "cache", "reads")
HOST = os.uname().nodename
# threads, memory GB, extra disk GB (reads cache + intermediates) per local step
RES = {"assemble": (24, 100, 250), "coassemble": (32, 150, 400), "metabat": (16, 40, 150), "maxbin": (16, 24, 200),
       "concoct": (16, 40, 200), "extract": (8, 32, 20), "checkm": (8, 32, 20)}
MAX_LIVE = {"assemble": 1, "coassemble": 1}
TF_STAGES = {"checkm": (8, 16), "extract": (8, 16)}    # TF threads / GB for the CheckM part
TF_FM_FOOTPRINT = (2, 8)                                # what such a step still uses on FM (fetch, extract, saves)            # MEGAHIT runs one at a time (memory)


def script_ver():
    return hashlib.md5(open(HW, "rb").read()).hexdigest()[:12]


def queued_kbase(t):
    a = t["attempts"][-1] if t["attempts"] else None
    return bool(a and not a.get("local") and a.get("job_id") and a.get("status") in ("created", "estimating", "queued"))


def local_live(ledger, stage=None):
    return [(k, a) for k, t in ledger["tasks"].items() for a in t["attempts"]
            if a.get("local") and a.get("hybrid") and a.get("status") == "running" and (stage is None or t["stage"] == stage)]


def pilot_ok(ledger, stage):
    return bool(ledger.get("hybrid_pilot_ok", {}).get(stage))


class Cycle:
    """Per-advance state: the gate is consulted per launch; 'waiting' counts ready steps the gate refused."""

    def __init__(self, P, ledger, driver):
        self.P, self.ledger, self.driver = P, ledger, driver
        self.waiting = 0
        self.waiting_tf = 0
        self.reasons = {}
        with G.Gate() as g:                                # this driver's live DAS Tool workers count too (8 threads)
            for t in ledger["tasks"].values():
                if t["stage"] != "dastool":
                    continue
                for a in t["attempts"]:
                    if a.get("local") and a.get("status") in P.LIVE and a.get("workdir") and a.get("pid") \
                            and a["workdir"] not in g.d["workers"]:
                        g.add(a["workdir"], driver, "dastool", t["sample"], 8, 16, a["pid"])

    def try_run(self, kb, stage, u, t, variant, params, out_name, unit_key, dry_run=False):
        """Called by advance() for a hybrid stage once params are known. Returns True if the step is handled
        (launched or deliberately waiting); the caller then `continue`s."""
        P, ledger = self.P, self.ledger
        res_stage = "coassemble" if stage == "coassemble" else stage
        thr, mem, disk = RES[res_stage]
        if stage not in ("assemble", "coassemble"):        # fit the driver's reserved share (2026-11-24); threads
            thr = min(thr, G.SHARE[self.driver][0])        # only change speed (MaxBin uses min(32, threads))
        live_stage = local_live(ledger, stage)
        if stage in MAX_LIVE and len(live_stage) >= MAX_LIVE[stage]:
            return True
        if not pilot_ok(ledger, stage) and live_stage:
            return True                                    # stage pilot: one local step until one is done
        if not pilot_ok(ledger, stage):                    # a failed pilot holds the stage (no unit after unit)
            failed = [k for k, tt in ledger["tasks"].items() if tt["stage"] == stage for x in tt["attempts"]
                      if x.get("hybrid") and x.get("status") == "error" and not x.get("transient")
                      and not x.get("pilot_cleared")]
            if failed:
                held = ledger.setdefault("hybrid_held", {})
                if held.get(stage) != failed:
                    held[stage] = failed
                    P.log(f"hybrid {stage}: local pilot failed ({failed}) - stage held for inspection; set "
                          f"'pilot_cleared': true on the failed attempt(s) after fixing")
                return True
        with G.Gate() as g:
            ok, why = g.can(self.driver, thr, mem, disk, P.WORK)
            remote = None
            if not ok and stage in TF_STAGES:              # CheckM (and extract's CheckM) can run on TF instead
                fm_ok, _ = g.can(self.driver, TF_FM_FOOTPRINT[0], TF_FM_FOOTPRINT[1], disk, P.WORK)
                tf_ok, tf_why = g.can_tf(self.driver, TF_STAGES[stage][0], TF_STAGES[stage][1]) if fm_ok else (False, "")
                if fm_ok and tf_ok:
                    ok, remote = True, "TF"
                    thr, mem = TF_FM_FOOTPRINT
                elif tf_why:
                    why += f"; TF: {tf_why}"
            if not ok:
                self.waiting += 1
                if stage in TF_STAGES:
                    self.waiting_tf += 1
                self.reasons[stage] = why
                return True
            if dry_run:
                return True
            takeover = queued_kbase(t)
            if takeover:
                a0 = t["attempts"][-1]
                try:
                    s = kb.states([a0["job_id"]]).get(a0["job_id"]) or {}
                except Exception as e:  # noqa: BLE001
                    P.log(f"{stage}:{unit_key}: takeover check failed ({type(e).__name__}); next cycle")
                    return True
                if s.get("status") not in ("created", "estimating", "queued"):
                    a0["status"] = s.get("status", a0["status"])
                    return True                            # started meanwhile (or finished): left to KBase
                try:
                    kb.ee2.cancel_job({"job_id": a0["job_id"]})
                except Exception as e:  # noqa: BLE001
                    P.log(f"{stage}:{unit_key}: cancel of {a0['job_id']} failed ({type(e).__name__}: {str(e)[:150]}); "
                          f"not taken over")
                    return True
                a0.update(status="terminated", takeover=True, finished_utc=P.now(),
                          error="cancelled while queued: taken over by the local hybrid run (user decision 2026-11-23)")
                P.log(f"{stage}:{unit_key}: KBase job {a0['job_id']} cancelled while queued - taken over locally")
                variant = P.choose_variant(stage, t) or variant
                got = P.params_for(variant, u, ledger, kb, *self.extra) if hasattr(self, "extra") else None
                if got and got not in ("skip",):
                    params, out_name = got
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            wd = os.path.join(P.WORK, "hybrid", stage, unit_key, f"attempt{len(t['attempts']) + 1}-{stamp}")
            os.makedirs(wd, exist_ok=True)
            p0 = params[0] if isinstance(params, list) else params
            p0 = dict(p0)
            if stage == "extract":
                p0["binned_contig_name"] = self.extract_name(u)
            reads = ([p0["read_library_ref"]] if "read_library_ref" in p0 else list(p0.get("reads_list") or []))
            spec = {"stage": "assemble" if stage in ("assemble", "coassemble") else stage, "unit": unit_key,
                    "ws_id": P.WS_ID, "ws_name": P.WS_NAME, "params": p0, "out_name": out_name, "threads": thr,
                    "mem_gb": int(mem * 0.9), "reads_cache": READS_CACHE, "reads_refs": reads,
                    "script_ver": script_ver(), "variant": variant}
            if remote:
                spec.update(remote=remote, remote_threads=TF_STAGES[stage][0])
            json.dump(spec, open(os.path.join(wd, "spec.json"), "w"), indent=1)
            a = {"local": True, "hybrid": True, "host": HOST, "pid": None, "workdir": wd, "variant": variant,
                 "params": params, "output_name": out_name, "status": "running", "submitted_utc": P.now(),
                 "threads": thr, "mem_gb": mem, "script_ver": spec["script_ver"]}
            if remote:
                a["remote"] = remote
            t["attempts"].append(a)
            t["state"] = "active"
            g.add(wd, self.driver, stage, unit_key, thr, mem)
            if remote:
                g.add(wd + "#TF", self.driver, stage, unit_key, TF_STAGES[stage][0], TF_STAGES[stage][1], host="TF")
            P.save_ledger(ledger)                          # write-ahead
            pr = subprocess.Popen(["nice", "-n", "5", PY, "-u", HW, wd], stdout=open(os.path.join(wd, "worker.log"), "a"),
                                  stderr=subprocess.STDOUT, start_new_session=True, cwd=wd)
            a["pid"] = pr.pid
            g.set_pid(wd, pr.pid)
            if remote:
                g.set_pid(wd + "#TF", pr.pid)
            P.save_ledger(ledger)
            P.log(f"{stage}:{unit_key}: started locally ({variant}, {thr} threads / {mem} GB, pid {pr.pid})"
                  + (f" - CheckM on TF ({TF_STAGES[stage][0]} threads)" if remote else ""))
        return True

    def finish(self):
        with G.Gate() as g:
            g.waiting(self.driver, self.waiting, self.waiting_tf)
        if self.waiting:
            self.P.log(f"hybrid: {self.waiting} local steps ready but waiting for resources: {self.reasons}")


def poll(P, kb, ledger):
    """Finished local hybrid steps -> ledger (the driver's own rules: floors, ladders, extract checks)."""
    for key, t in ledger["tasks"].items():
        for a in t["attempts"]:
            if not (a.get("local") and a.get("hybrid") and a.get("status") == "running"):
                continue
            rp = os.path.join(a["workdir"], "result.json")
            if not os.path.exists(rp):
                if P.worker_alive(a.get("pid"), a["workdir"]):
                    age_h = (time.time() - datetime.fromisoformat(a["submitted_utc"]).timestamp()) / 3600
                    if age_h > 72 and not a.get("stale_logged"):
                        a["stale_logged"] = P.now()
                        P.log(f"{key}: local step {a['pid']} running for {age_h:.0f} h - see {a['workdir']}/worker.log")
                    continue
                a.update(status="error", finished_utc=P.now(), error="local worker exited without a result")
                t["state"] = "retry"
                P.log(f"{key}: local worker {a.get('pid')} exited without a result - retried")
                continue
            r = json.load(open(rp))
            a.update(finished_utc=r.get("finished_utc") or P.now(), versions=r.get("versions"), host=r.get("host"))
            stage = t["stage"]
            if r["status"] == "empty":
                a["status"] = "completed"
                t.update(state="empty", reason=r.get("error"))
                P.log(f"{key}: empty (local) - {r.get('error')}")
                continue
            if r["status"] != "done":
                a.update(status="error", error=(r.get("error") or "")[:2000])
                if r.get("intrinsic"):
                    a["intrinsic"] = True
                t["state"] = "retry"
                P.log(f"{key}: local step failed{' (intrinsic)' if r.get('intrinsic') else ''}: {a['error'][:200]}")
                continue
            ref = r["ref"]
            a["status"] = "completed"
            if stage == "extract":
                fake = {"job_output": {"result": [{"assembly_ref_list": r["bin_refs"],
                                                   "assembly_set_ref": ref if len(r["bin_refs"]) > 1 else None}]}}
                ok, out, reason = P.on_completed(kb, ledger, key, t, a, fake)
            elif stage == "checkm":
                ok, out, reason = True, ref, None
            else:
                ok, reason = P.validate_output(kb, stage, ref, a["params"])
                out = ref
            if not ok:
                a.update(status="error", error=reason, rejected=True, rejected_ref=ref)
                if stage in ("assemble", "coassemble"):
                    a["intrinsic"] = True
                t["state"] = "retry"
                P.log(f"{key}: local result rejected: {reason}")
                continue
            t.update(state="done", output_ref=out, local=True)
            ledger.setdefault("hybrid_pilot_ok", {}).setdefault(stage, key)
            P.log(f"{key}: done locally -> {out}")
            if stage == "checkm":
                t["collect_pending"] = True
            if stage == "extract" and r.get("checkm_report_ref"):
                unit = t["sample"]
                ck = P.task(ledger, "checkm", unit)
                if ck["state"] not in ("done",):
                    ck["attempts"].append({"local": True, "hybrid": True, "host": r.get("host"), "variant": "checkm",
                                           "via": "extract (local CheckM 1.1.2 in the same step)",
                                           "status": "completed", "submitted_utc": a["submitted_utc"],
                                           "finished_utc": a["finished_utc"], "versions": r.get("versions"),
                                           "workdir": a["workdir"], "n_pass_50_10": r.get("n_pass_50_10")})
                    ck.update(state="done", output_ref=r["checkm_report_ref"], collect_pending=True, local=True)
                    ledger.setdefault("hybrid_pilot_ok", {}).setdefault("checkm", f"checkm:{unit}")
                    P.log(f"checkm:{unit}: done locally (with extract) -> {r['checkm_report_ref']}")


def write_needs(driver, refs):
    p = os.path.join(HERE, f"needs_{driver}.json")
    with open(p + ".tmp", "w") as f:
        json.dump(sorted(set(refs)), f)
    os.replace(p + ".tmp", p)
