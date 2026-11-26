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
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "assembly_binning_2026", "local")))
from kbrpc import rpc  # noqa: E402  (read-only)

HW = os.path.join(HERE, "hw.py")
PY = os.path.expanduser("~/Documents/py_venv/bin/python")
READS_CACHE = os.path.join(HERE, "cache", "reads")
HOST = os.uname().nodename
# threads, memory GB, extra disk GB (reads cache + intermediates) per local step
RES = {"assemble": (24, 100, 250), "coassemble": (32, 150, 400), "metabat": (16, 40, 150), "maxbin": (16, 24, 200),
       "concoct": (16, 40, 200), "extract": (8, 32, 20), "checkm": (8, 32, 20)}
MAX_LIVE = {"assemble": 1, "coassemble": 1}            # per route (FM / poplar)
# user decision 2026-11-25: MEGAHIT on FM only for small libraries. providencecreek_0007 (23.5 GB gz, 46.6 Gbp) passed
# 161 GB RSS in its k29 iterate (-m caps only the SdBG build), i.e. >= ~7 GB RAM per GB gz and still rising. With a
# 100 GB per-job budget, 10 GB gz allows ~10x; a memory watchdog (hw.py) stops anything that grows past 100 GB RSS
# or drives MemAvailable below 20 GB and requeues it to poplar, which resumes from the checkpoints.
FM_MAX_LIB_GB = 10.0
FM_RSS_FLOOR_GB = 20
POPLAR_RES = (32, 300)                                  # poplar threads / GB (-m 270 GB for the SdBG builds)
POPLAR_FM_FOOTPRINT = (1, 4)                            # the FM-side worker (copy, poll, save)
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
        self.waiting_poplar = 0
        self.fm_assemble_waiting = False
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
        route = self.route_of(stage, t, params) if stage in ("assemble", "coassemble") else "FM"
        if route is None:
            return True
        pkey = stage + ("@poplar" if route == "poplar" else "")
        if stage in ("metabat", "maxbin", "concoct") and self.fm_assemble_waiting:
            self.waiting += 1                              # an FM assembly waits for room: it goes first
            self.reasons[stage] = "an FM assembly is next in line"
            return True
        if stage not in ("assemble", "coassemble"):        # fit the driver's reserved share (2026-11-24); threads
            thr = min(thr, G.SHARE[self.driver][0])        # only change speed (MaxBin uses min(32, threads))
        live_stage = local_live(ledger, stage)
        if stage in ("assemble", "coassemble"):
            live_stage = [x for x in live_stage if (x[1].get("remote") == "poplar") == (route == "poplar")]
        if stage in MAX_LIVE and len(live_stage) >= MAX_LIVE[stage]:
            return True
        if not pilot_ok(ledger, pkey) and live_stage:
            return True                                    # stage pilot: one local step until one is done
        if not pilot_ok(ledger, pkey):                     # a failed pilot holds the stage (no unit after unit)
            failed = [k for k, tt in ledger["tasks"].items() if tt["stage"] == stage for x in tt["attempts"]
                      if x.get("hybrid") and x.get("status") == "error" and not x.get("transient")
                      and not x.get("pilot_cleared")
                      and (stage not in ("assemble", "coassemble") or (x.get("remote") == "poplar") == (route == "poplar"))]
            if failed:
                held = ledger.setdefault("hybrid_held", {})
                if held.get(pkey) != failed:
                    held[pkey] = failed
                    P.log(f"hybrid {stage}: local pilot failed ({failed}) - stage held for inspection; set "
                          f"'pilot_cleared': true on the failed attempt(s) after fixing")
                return True
        with G.Gate() as g:
            remote = None
            if route == "poplar":
                need = self.poplar_need_gb(t, variant)
                ok, why = g.can(self.driver, *POPLAR_FM_FOOTPRINT, 0, P.WORK, share=False)
                if ok:
                    ok, why = g.can_poplar(self.driver, POPLAR_RES[0], POPLAR_RES[1], need)
                if not ok:
                    self.waiting_poplar += 1
                    self.reasons[pkey] = why
                    return True
                remote = "poplar"
            else:
                ok, why = g.can(self.driver, thr, mem, disk, P.WORK)
            if not ok and stage in ("assemble", "coassemble"):
                self.fm_assemble_waiting = True
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
            if remote == "TF":
                spec.update(remote=remote, remote_threads=TF_STAGES[stage][0])
            elif remote == "poplar":
                spec.update(remote="poplar", threads=POPLAR_RES[0], mem_gb=int(POPLAR_RES[1] * 0.9))
                prev = [x for x in t["attempts"] if x.get("overflow") and x.get("workdir") and x.get("variant") == variant]
                ck = prev and os.path.join(prev[-1]["workdir"], "tmp", "megahit_out", "checkpoints.txt")
                if ck and os.path.exists(ck) and not t.get("no_resume"):
                    spec["resume_from"] = os.path.join(prev[-1]["workdir"], "tmp")
            elif stage in ("assemble", "coassemble"):
                spec.update(rss_limit_gb=mem, mem_floor_gb=FM_RSS_FLOOR_GB)
            json.dump(spec, open(os.path.join(wd, "spec.json"), "w"), indent=1)
            a = {"local": True, "hybrid": True, "host": HOST, "pid": None, "workdir": wd, "variant": variant,
                 "params": params, "output_name": out_name, "status": "running", "submitted_utc": P.now(),
                 "threads": thr, "mem_gb": mem, "script_ver": spec["script_ver"]}
            if remote:
                a["remote"] = remote
            if spec.get("resume_from"):
                a["resume_from"] = spec["resume_from"]
            t["attempts"].append(a)
            t["state"] = "active"
            if remote == "poplar":
                a.update(threads=POPLAR_RES[0], mem_gb=POPLAR_RES[1], compute_host="poplar")
                thr, mem = POPLAR_FM_FOOTPRINT
            g.add(wd, self.driver, stage, unit_key, thr, mem)
            if remote == "TF":
                g.add(wd + "#TF", self.driver, stage, unit_key, TF_STAGES[stage][0], TF_STAGES[stage][1], host="TF")
            elif remote == "poplar":
                g.add(wd + "#poplar", self.driver, stage, unit_key, POPLAR_RES[0], POPLAR_RES[1], host="poplar")
            P.save_ledger(ledger)                          # write-ahead
            pr = subprocess.Popen(["nice", "-n", "5", PY, "-u", HW, wd], stdout=open(os.path.join(wd, "worker.log"), "a"),
                                  stderr=subprocess.STDOUT, start_new_session=True, cwd=wd)
            a["pid"] = pr.pid
            g.set_pid(wd, pr.pid)
            if remote:
                g.set_pid(wd + "#" + remote, pr.pid)
            P.save_ledger(ledger)
            P.log(f"{stage}:{unit_key}: started locally ({variant}, {thr} threads / {mem} GB, pid {pr.pid})"
                  + (f" - CheckM on TF ({TF_STAGES[stage][0]} threads)" if remote == "TF" else "")
                  + (f" - MEGAHIT on poplar ({POPLAR_RES[0]} threads"
                     + (", resuming from FM checkpoints)" if spec.get("resume_from") else ")") if remote == "poplar" else ""))
        return True

    def lib_gb(self, t, params):
        if t.get("lib_gz_gb") is None:
            p0 = params[0] if isinstance(params, list) else params
            ref = p0.get("read_library_ref")
            d = rpc("ws", "Workspace.get_objects2", [{"objects": [{"ref": ref, "included": ["/lib1/size"]}]}])[0]
            t["lib_gz_gb"] = round(d["data"][0]["data"]["lib1"]["size"] / 1e9, 2)
        return t["lib_gz_gb"]

    def route_of(self, stage, t, params):
        """FM for small per-sample libraries; poplar for large ones, co-assemblies and anything the FM memory
        watchdog stopped (user decision 2026-11-25). None if the size cannot be read now."""
        if stage == "coassemble" or t.get("route") == "poplar":
            return "poplar"
        try:
            return "poplar" if self.lib_gb(t, params) > FM_MAX_LIB_GB else "FM"
        except Exception as e:  # noqa: BLE001
            self.P.log(f"{stage}: library size lookup failed ({type(e).__name__}); next cycle")
            return None

    def poplar_need_gb(self, t, variant):
        prev = [x for x in t["attempts"] if x.get("overflow") and x.get("workdir") and x.get("variant") == variant]
        if prev and not t.get("no_resume"):
            return 250                                     # copied checkpoints (~90-120 GB) plus growth
        return 8 * (t.get("lib_gz_gb") or 30)

    def finish(self):
        with G.Gate() as g:
            g.waiting(self.driver, self.waiting, self.waiting_tf, self.waiting_poplar)
        if self.waiting or self.waiting_poplar:
            self.P.log(f"hybrid: {self.waiting} local + {self.waiting_poplar} poplar steps ready but waiting: {self.reasons}")


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
            if r.get("overflow"):                           # FM memory watchdog: not a tool failure; poplar next
                a.update(status="terminated", overflow=True, error=(r.get("error") or "")[:2000])
                t.update(state="retry", route="poplar")
                P.log(f"{key}: stopped by the FM memory watchdog ({a['error'][:160]}) - requeued to poplar "
                      f"(resumes from the checkpoints)")
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
            ledger.setdefault("hybrid_pilot_ok", {}).setdefault(
                stage + ("@poplar" if a.get("remote") == "poplar" else ""), key)
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
