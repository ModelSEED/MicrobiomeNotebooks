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
MAX_LIVE = {"assemble": 1, "coassemble": 1}            # on FM; poplar is limited by gate.can_poplar (2026-11-25, 2nd)
# user decision 2026-11-25: MEGAHIT on FM only for small libraries. providencecreek_0007 (23.5 GB gz, 46.6 Gbp) passed
# 161 GB RSS in its k29 iterate (-m caps only the SdBG build), i.e. >= ~7 GB RAM per GB gz and still rising. With a
# 100 GB per-job budget, 10 GB gz allows ~10x; a memory watchdog (hw.py) stops anything that grows past 100 GB RSS
# or drives MemAvailable below 20 GB and requeues it to poplar, which resumes from the checkpoints.
# Second decision 2026-11-25: FM takes per-sample libraries up to 16 GB gz, one assembly at a time; the gate is asked
# for the expected peak (7 GB per GB gz + 10, at most the per-sample share of 120 GB) and the watchdog stops MEGAHIT at
# the share (120 GB RSS) or MemAvailable < 20 GB. Every remaining assembly goes to the first route that can take it:
# FM (small libraries), poplar (concurrent jobs), else a KBase job (fallback queue; taken over while still queued).
FM_MAX_LIB_GB = 16.0
FM_GB_PER_GZ_GB = 7
FM_RSS_FLOOR_GB = 20
POPLAR_RES = (32, 300)                                  # poplar threads / GB (-m 270 GB for the SdBG builds)
POPLAR_FM_FOOTPRINT = (1, 4)                            # the FM-side worker (copy, poll, save)
TF_STAGES = {"checkm": (8, 16), "extract": (8, 16),    # TF threads / GB for the CheckM part
             "maxbin": (10, 8), "concoct": (10, 20)}     # run_MaxBin.pl / concoct clustering (2026-11-25)
TF_FM = {"maxbin": (1, 2), "concoct": (1, 2)}           # their FM part only fetches reads, ships them and polls: the
                                                        # mapping runs on TF too (FM is often loaded past its cores)
BINNERS = ("metabat", "maxbin", "concoct")
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
        (launched or deliberately waiting); the caller then `continue`s. False (assemblies only): no local route can
        take it now - the caller submits the KBase job (fallback queue, user decision 2026-11-25)."""
        P, ledger = self.P, self.ledger
        res_stage = "coassemble" if stage == "coassemble" else stage
        thr, mem, disk = RES[res_stage]
        if t.get("dedup_hold"):                            # a superseded local result waits for a person (no-duplicate rule)
            return True
        is_asm = stage in ("assemble", "coassemble")
        routes = self.routes_of(stage, t, params) if is_asm else ["FM"]
        if routes is None:
            return True
        if is_asm and queued_kbase(t) and not any(self.route_open(stage, t, r) for r in routes):
            if "poplar" in routes:                         # still waiting for poplar: keeps per-sample's priority there
                self.waiting_poplar += 1
            return True                                    # its KBase job waits in the queue (fallback route)
        fm_held = stage in BINNERS and self.fm_assemble_waiting   # an FM assembly waits for room: it gets FM first;
                                                                  # the binner may still go to TF or KBase
        if not is_asm:                                     # fit the driver's reserved share (2026-11-24); threads
            thr = min(thr, G.SHARE[self.driver][0])        # only change speed (MaxBin uses min(32, threads))
            live_stage = local_live(ledger, stage)
            if stage in MAX_LIVE and len(live_stage) >= MAX_LIVE[stage]:
                return True
            if not self.pilot_open(stage, stage, live_stage):
                return True
        why_all = []
        with G.Gate() as g:
            remote = None
            ok = False
            work_fs = None
            for route in routes:
                if is_asm:
                    if not self.route_open(stage, t, route):
                        why_all.append(f"{route}: busy / pilot")
                        continue
                if route == "poplar":
                    need = self.poplar_need_gb(t, variant)
                    ok, why = g.can(self.driver, *POPLAR_FM_FOOTPRINT, 0, P.WORK, share=False)
                    if ok:
                        ok, why, work_fs = g.can_poplar(self.driver, POPLAR_RES[0], POPLAR_RES[1], need)
                    if ok:
                        remote = "poplar"
                        break
                    why_all.append(f"poplar: {why}")
                    continue
                if is_asm:                                 # FM MEGAHIT: ask for the expected peak, watchdog at the share
                    mem = min(G.SHARE[self.driver][1], int(FM_GB_PER_GZ_GB * self.lib_gb(t, params)) + 10)
                if fm_held:
                    ok, why = False, "an FM assembly is next in line"
                    break
                ok, why = g.can(self.driver, thr, mem, disk, P.WORK)
                if ok:
                    break
                if is_asm:
                    why_all.append(f"FM: {why}")
                    # the binners wait (this assembly gets FM's next room) only if yielding our own steps would
                    # make that room; memory that other users hold is not freed by holding our binners
                    ot, om = g.used("FM", self.driver)
                    tt, tm = g.used("FM")
                    if (G.mem_available_gb() + om >= mem + G.MEM_MARGIN_GB and tm - om + mem <= G.TOTAL_MEM_GB
                            and tt - ot + thr <= G.TOTAL_THREADS):
                        self.fm_assemble_waiting = True
            if is_asm and not ok:
                if "poplar" in routes:
                    self.waiting_poplar += 1
                if "FM" in routes:
                    self.waiting += 1
                self.reasons[f"{stage}:{unit_key}"] = "; ".join(why_all)
                if self.kbase_ok(t) and not dry_run:
                    t["hybrid_route_kbase"] = P.now()      # the driver submits the KBase job now (fallback route)
                    return False
                return True
            if not ok and stage in TF_STAGES:              # CheckM, extract's CheckM, MaxBin2 and CONCOCT can
                fthr, fmem = TF_FM.get(stage, TF_FM_FOOTPRINT)    # compute on TF; their FM part stays small
                fm_ok, fm_why = g.can(self.driver, fthr, fmem, disk, P.WORK, share=stage not in BINNERS)
                if fm_ok and stage in BINNERS and not self.tf_binner_open(stage):
                    fm_ok, fm_why = False, f"{stage}@TF pilot"
                tf_ok, tf_why = g.can_tf(self.driver, TF_STAGES[stage][0], TF_STAGES[stage][1]) if fm_ok else (False, "")
                if fm_ok and tf_ok:
                    ok, remote = True, "TF"
                    thr, mem = fthr, fmem
                else:
                    why += f"; TF: {tf_why or 'FM part: ' + str(fm_why)}"
            if not ok:
                self.waiting += 1
                if stage in TF_STAGES:
                    self.waiting_tf += 1
                self.reasons[stage] = why
                if stage in BINNERS and self.kbase_ok(t) and not dry_run:
                    t["hybrid_route_kbase"] = P.now()      # first route that can take it: KBase (fallback queue,
                    return False                           # user decision 2026-11-25; taken over while queued)
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
                try:                                       # no-duplicate rule (2026-11-25): EE2 must confirm the end
                    s2 = (kb.states([a0["job_id"]]).get(a0["job_id"]) or {}).get("status")
                except Exception as e:  # noqa: BLE001
                    s2 = None
                if s2 in ("running", "completed"):         # it started or finished meanwhile: KBase keeps it
                    a0["status"] = "running"               # (poll() follows it; a completed one is adopted there)
                    P.log(f"{stage}:{unit_key}: KBase job {a0['job_id']} is {s2} after the cancel - left to KBase")
                    return True
                if s2 not in ("terminated", "error"):
                    P.log(f"{stage}:{unit_key}: cancel of {a0['job_id']} not confirmed yet (EE2: {s2}) - no local start")
                    return True
                a0.update(status="terminated", takeover=True, finished_utc=P.now(), ee2_after_cancel=s2,
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
                spec.update(remote="poplar", threads=POPLAR_RES[0], mem_gb=int(POPLAR_RES[1] * 0.9),
                            disk_gb=self.poplar_need_gb(t, variant), poplar_work_fs=work_fs)
                prev = [x for x in t["attempts"] if x.get("overflow") and x.get("workdir") and x.get("variant") == variant]
                ck = prev and os.path.join(prev[-1]["workdir"], "tmp", "megahit_out", "checkpoints.txt")
                if ck and os.path.exists(ck) and not t.get("no_resume"):
                    spec["resume_from"] = os.path.join(prev[-1]["workdir"], "tmp")
            elif stage in ("assemble", "coassemble"):
                spec.update(rss_limit_gb=G.SHARE[self.driver][1], mem_floor_gb=FM_RSS_FLOOR_GB)
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
                g.add(wd + "#poplar", self.driver, stage, unit_key, POPLAR_RES[0], POPLAR_RES[1], host="poplar",
                      disk_gb=spec["disk_gb"], work_fs=work_fs)
            P.save_ledger(ledger)                          # write-ahead
            pr = subprocess.Popen(["nice", "-n", "5", PY, "-u", HW, wd], stdout=open(os.path.join(wd, "worker.log"), "a"),
                                  stderr=subprocess.STDOUT, start_new_session=True, cwd=wd)
            a["pid"] = pr.pid
            g.set_pid(wd, pr.pid)
            if remote:
                g.set_pid(wd + "#" + remote, pr.pid)
            P.save_ledger(ledger)
            P.log(f"{stage}:{unit_key}: started locally ({variant}, {thr} threads / {mem} GB, pid {pr.pid})"
                  + (f" - {'run' if stage in BINNERS else 'CheckM'} on TF ({TF_STAGES[stage][0]} threads)" if remote == "TF" else "")
                  + (f" - MEGAHIT on poplar ({POPLAR_RES[0]} threads"
                     + (", resuming from FM checkpoints" if spec.get("resume_from") else "")
                     + f", work dir on {work_fs})" if remote == "poplar" else ""))
        return True

    def lib_gb(self, t, params):
        if t.get("lib_gz_gb") is None:
            p0 = params[0] if isinstance(params, list) else params
            ref = p0.get("read_library_ref")
            d = rpc("ws", "Workspace.get_objects2", [{"objects": [{"ref": ref, "included": ["/lib1/size"]}]}])[0]
            t["lib_gz_gb"] = round(d["data"][0]["data"]["lib1"]["size"] / 1e9, 2)
        return t["lib_gz_gb"]

    def routes_of(self, stage, t, params):
        """Local routes in order of preference (user decisions 2026-11-25): FM for per-sample libraries <= 16 GB gz,
        then poplar; poplar only for larger ones, co-assemblies and anything the FM watchdog stopped. t['route']
        pins one route (the FM pilot; an overflow -> poplar). None if the library size cannot be read now."""
        try:
            gz = self.lib_gb(t, params)
        except Exception as e:  # noqa: BLE001
            self.P.log(f"{stage}: library size lookup failed ({type(e).__name__}); next cycle")
            return None
        if t.get("route") in ("FM", "poplar"):
            return [t["route"]]
        if stage == "coassemble" or gz > FM_MAX_LIB_GB:
            return ["poplar"]
        return ["FM", "poplar"]

    def pilot_open(self, stage, pkey, live):
        """Stage pilot: one local step until one is done; a failed pilot holds the stage (no unit after unit)."""
        P, ledger = self.P, self.ledger
        if pilot_ok(ledger, pkey):
            return True
        if live:
            return False
        poplar = pkey.endswith("@poplar")
        failed = [k for k, tt in ledger["tasks"].items() if tt["stage"] == stage for x in tt["attempts"]
                  if x.get("hybrid") and x.get("status") == "error" and not x.get("transient")
                  and not x.get("pilot_cleared")
                  and (stage not in ("assemble", "coassemble") or (x.get("remote") == "poplar") == poplar)]
        if failed:
            held = ledger.setdefault("hybrid_held", {})
            if held.get(pkey) != failed:
                held[pkey] = failed
                P.log(f"hybrid {pkey}: local pilot failed ({failed}) - held for inspection; set "
                      f"'pilot_cleared': true on the failed attempt(s) after fixing")
            return False
        return True

    def route_open(self, stage, t, route):
        """FM: one MEGAHIT at a time (and the pilot rule). poplar: concurrency is the gate's (can_poplar); until a
        fresh poplar run has been verified, one fresh run at a time (a resume from checkpoints is not counted)."""
        live = [x for x in local_live(self.ledger, stage) if (x[1].get("remote") == "poplar") == (route == "poplar")]
        if route == "FM":
            return len(live) < MAX_LIVE[stage] and self.pilot_open(stage, stage, live)
        fresh = [x for x in live if not x[1].get("resume_from")]
        return self.pilot_open(stage, stage + "@poplar", fresh)

    def tf_binner_open(self, stage):
        """TF binner pilot (MaxBin2 / CONCOCT on TF, 2026-11-25): one at a time until one is done; a failed one holds
        the TF route of that stage (FM and KBase still take the stage) until 'pilot_cleared' is set on it."""
        if pilot_ok(self.ledger, stage + "@TF"):
            return True
        tf_att = [x for tt in self.ledger["tasks"].values() if tt["stage"] == stage for x in tt["attempts"]
                  if x.get("hybrid") and x.get("remote") == "TF"]
        if any(x.get("status") == "running" for x in tf_att):
            return False
        bad = [x for x in tf_att if x.get("status") == "error" and not x.get("transient") and not x.get("pilot_cleared")]
        return not bad

    def kbase_ok(self, t):
        """KBase is the fallback queue for an assembly no local route can take now - not for a pinned route (the FM
        pilot; an overflow that resumes from its checkpoints) and not while its last attempt is still in flight."""
        a = t["attempts"][-1] if t["attempts"] else None
        return not t.get("route") and not (a and a.get("status") in ("created", "estimating", "queued", "running"))

    def poplar_need_gb(self, t, variant):
        """Projected size of MEGAHIT's work dir (output + tmp; the reads are staged on /scratch1 since 2026-11-25, 3rd
        decision): ~7 GB per GB gz (providencecreek, 23.5 GB gz: 170 GB with its reads at k29); 250 for a resume."""
        prev = [x for x in t["attempts"] if x.get("overflow") and x.get("workdir") and x.get("variant") == variant]
        if prev and not t.get("no_resume"):
            return 250
        return 7 * (t.get("lib_gz_gb") or 30)

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
            if r["status"] == "superseded":                 # another version arrived first: keep that one
                a.update(status="terminated", dedup=True, superseded_by=r.get("superseded_by"),
                         error=(r.get("error") or "")[:1000])
                if stage in ("extract", "checkm"):          # no name-based adoption for these: a person decides
                    t.update(state="retry", dedup_hold=r.get("superseded_by"))
                    P.log(f"{key}: local result NOT pushed - {a['error'][:200]}; task held (dedup_hold) for inspection")
                else:                                       # the driver adopts the existing object after validating it
                    t["state"] = "retry"
                    P.log(f"{key}: local result NOT pushed - {a['error'][:200]}; the existing version is checked next")
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
                stage + ("@poplar" if a.get("remote") == "poplar" else
                         "@TF" if a.get("remote") == "TF" and stage in BINNERS else ""), key)
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


# ------------------------------------------------------------------------------------------------------ no duplicates
# User decision 2026-11-25 (relayed by the coordinating session): no work is ever duplicated between KBase and
# FM / TF / poplar, in either direction. Enforced here for both drivers:
#  1. one live attempt per task across all routes; a local start after a KBase job only once EE2 confirms the
#     cancel (try_run's takeover); 2. no KBase (re)submit while any attempt is live (the drivers' submit guard; a local
#     result is pushed by its worker before result.json, so "pending push" is a live local attempt);
#  3. every cycle: dedup() - a task with several live attempts keeps the one furthest along (running > queued; ties:
#     local, then the older) and the others are cancelled / killed, recorded status 'terminated' + dedup (not counted);
#     EE2 jobs live in the driver's workspace that the ledger does not know are adopted if their params match a task
#     that has nothing live, otherwise cancelled (only the drivers' own apps); jobs the ledger records as finished but
#     EE2 shows live are cancelled if the task is done or has another live attempt, else re-adopted;
#  4. hw.py checks before each push that no other version of the output arrived meanwhile (result 'superseded').
LOCAL_LIVE = ("running", "queued", "started")


def live_attempts(t, LIVE):
    return [a for a in t["attempts"] if (a.get("local") and a.get("status") in LOCAL_LIVE)
            or (not a.get("local") and a.get("job_id") and a.get("status") in LIVE)]


def _desc(a):
    return (f"local:{a.get('compute_host') or a.get('remote') or a.get('host', 'FM')}" if a.get("local")
            else f"kbase:{a.get('job_id')} ({a.get('status')})")


def _cancel_job(P, kb, jid):
    """Cancel an EE2 job and confirm: returns its status after the cancel (terminated/error/completed/...)."""
    try:
        kb.ee2.cancel_job({"job_id": jid})
    except Exception as e:  # noqa: BLE001
        P.log(f"dedup: cancel of {jid} raised {type(e).__name__}: {str(e)[:150]}")
    try:
        return (kb.states([jid]).get(jid) or {}).get("status")
    except Exception:  # noqa: BLE001
        return None


def _kill_local(P, key, a):
    """Stop a local attempt on FM and on its remote host (TF / poplar): processes whose command line holds the
    attempt's workdir / run id. Returns True if nothing of it is left running on FM."""
    import signal
    wd = a.get("workdir") or ""
    if not wd:
        return True
    for pid in os.listdir("/proc"):
        if not pid.isdigit() or int(pid) == os.getpid():
            continue
        try:
            if wd.encode() in open(f"/proc/{pid}/cmdline", "rb").read():
                os.kill(int(pid), signal.SIGTERM)
        except OSError:
            pass
    host = a.get("remote") or a.get("compute_host")
    if host in ("TF", "poplar"):
        rid = os.path.basename(wd.rstrip("/")) + "_" + key.split(":", 1)[1]
        pat = f"[{rid[0]}]{rid[1:]}"                       # the bracket keeps the ssh shell itself out of the match
        subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", host,
                        f"for p in $(pgrep -u $(id -u) -f '{pat}'); do kill $p; done; true"],
                       capture_output=True, text=True, timeout=90)
    time.sleep(5)
    return not P.worker_alive(a.get("pid"), wd)


def dedup(P, kb, ledger, dry_run=False):
    """Rule 3: one live attempt per task; EE2 orphans adopted or cancelled. Run every cycle right after poll()."""
    if dry_run:
        return
    # (a) several live attempts on one task
    for key, t in ledger["tasks"].items():
        la = live_attempts(t, P.LIVE)
        if len(la) < 2:
            continue
        idx = {id(a): i for i, a in enumerate(t["attempts"])}
        keep = max(la, key=lambda a: (2 if a.get("local") or a.get("status") == "running" else 1,
                                      1 if a.get("local") else 0, -idx[id(a)]))
        for a in la:
            if a is keep:
                continue
            if a.get("local"):
                ok, what = _kill_local(P, key, a), "killed"
            else:
                st = _cancel_job(P, kb, a["job_id"])
                ok, what = st in ("terminated", "error", "completed"), f"cancelled (EE2: {st})"
                if st == "completed":                      # finished first after all: poll() handles the output;
                    a["status"] = "running"                 # the local one would then be the duplicate next cycle
                    P.log(f"{key}: dedup - KBase job {a['job_id']} completed meanwhile; its output is checked next")
                    continue
            if not ok:
                P.log(f"{key}: dedup - could not stop {_desc(a)} yet; retried next cycle")
                continue
            gone = _desc(a)
            a.update(status="terminated", dedup=True, finished_utc=P.now(),
                     error=f"duplicate of {_desc(keep)}: {what} by the no-duplicate check (user decision 2026-11-25)")
            P.log(f"{key}: dedup - kept {_desc(keep)}, {what} {gone} (not counted)")
        t["state"] = "active"
    # (b) EE2 jobs live in this workspace that the ledger does not know, or records as finished
    try:
        res = kb.ee2.check_workspace_jobs({"workspace_id": P.WS_ID, "return_list": 1,
                                           "exclude_fields": ["job_output"]})
    except Exception as e:  # noqa: BLE001
        P.log(f"dedup: EE2 workspace job list failed ({type(e).__name__}); next cycle")
        return
    known = {a["job_id"]: (k, t, a) for k, t in ledger["tasks"].items() for a in t["attempts"] if a.get("job_id")}
    our_apps = {v[0] for v in P.APPS.values()}
    for j in res.get("job_states", []):
        jid, st = j.get("job_id"), j.get("status")
        if st not in P.LIVE:
            continue
        if jid in known:
            k, t, a = known[jid]
            if a.get("status") in P.LIVE:
                continue
            others = [x for x in live_attempts(t, P.LIVE) if x is not a]
            if t["state"] in ("done", "empty", "failed") or others:
                s2 = _cancel_job(P, kb, jid)
                a.update(dedup=True, dedup_note=f"EE2 still showed it {st} although the ledger had it "
                                                 f"{a.get('status')}: cancelled (EE2: {s2})")
                P.log(f"{k}: dedup - KBase job {jid} was {st} in EE2 but finished in the ledger: cancelled ({s2})")
            else:
                a.update(status=st, readopted=P.now())
                t["state"] = "active"
                P.log(f"{k}: dedup - KBase job {jid} is {st} in EE2 (ledger: finished): re-adopted")
            continue
        ji = j.get("job_input") or {}
        app, ptxt = ji.get("app_id"), json.dumps(ji.get("params"))
        match = None
        for k, t in ledger["tasks"].items():
            if t.get("submitting") or t.get("submit_errors"):
                continue                                   # the driver's own lost-submit reconciliation
            for a in reversed(t["attempts"]):
                v = a.get("variant")
                if v in P.APPS and P.APPS[v][0] == app:
                    needle = P.recon_key(t["stage"], a.get("params"), a.get("output_name"))
                    if needle and json.dumps(needle) in ptxt:
                        match = (k, t, a)
                    break
            if match:
                break
        if match and not live_attempts(match[1], P.LIVE) and match[1]["state"] not in ("done", "empty", "failed"):
            k, t, a = match
            t["attempts"].append({"job_id": jid, "variant": a.get("variant"), "app_id": app, "params": ji.get("params"),
                                  "output_name": a.get("output_name"), "submitted_utc": P.now(), "status": st,
                                  "reconciled": True, "dedup_adopted": True})
            t["state"] = "active"
            P.log(f"{k}: dedup - adopted EE2 orphan {jid} ({app}, {st})")
        elif match or app in our_apps:
            s2 = _cancel_job(P, kb, jid)
            P.log(f"dedup - EE2 orphan {jid} ({app}, {st}) "
                  + (f"duplicates {match[0]}" if match else "matches no task") + f": cancelled (EE2: {s2})")
            ledger.setdefault("dedup_cancelled", []).append({"job_id": jid, "app_id": app, "utc": P.now(),
                                                              "task": match[0] if match else None, "ee2": s2})
        else:
            if jid not in ledger.setdefault("dedup_foreign_logged", []):
                ledger["dedup_foreign_logged"].append(jid)
                P.log(f"dedup - EE2 job {jid} ({app}) in ws {P.WS_ID} is not one of this driver's apps: left alone")
