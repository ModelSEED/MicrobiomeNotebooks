#!/usr/bin/env python3
"""Watch the 11 MEGAHIT retry jobs (jobs.json) every 30 min. Logs state changes to watch.events; flags STALL when a job
has run > 72 h or its log has been silent > 24 h (the 2026-07 attempts hung silently for ~2 weeks). Read-only on KBase.
Exits when every job is terminal (writes ALLDONE) or after 21 days. Usage: KB_AUTH_TOKEN=... setsid nohup python -I -B watch_megahit.py &"""
import json, os, sys, time
sys.path.insert(0, "/home/freiburger/gtdbtk_kbase")
import kbase_mags as k
EE2 = "https://kbase.us/services/ee2"
HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(HERE, "watch.events")
def ev(msg):
    open(EV, "a").write(f"[{time.strftime('%F %T')}] {msg}\n")
jobs = json.load(open(os.path.join(HERE, "jobs.json")))
last, flagged, t_end = {}, set(), time.time() + 21 * 86400
ev(f"START watching {len(jobs)} MEGAHIT jobs (pid {os.getpid()})")
while time.time() < t_end:
    try:
        st = k.rpc(EE2, "execution_engine2.check_jobs", [{"job_ids": [v["job_id"] for v in jobs.values()],
                                                         "exclude_fields": ["job_input", "job_output"]}])[0]["job_states"]
    except Exception as e:
        ev(f"WARN check_jobs failed: {e!s:.200}"); time.sleep(1800); continue
    by = {s["job_id"]: s for s in st}
    for name, v in jobs.items():
        s = by.get(v["job_id"], {})
        status = s.get("status")
        if last.get(name) != status:
            extra = f" error={str((s.get('error') or {}).get('message') or s.get('errormsg') or '')[:160]}" if status == "error" else ""
            ev(f"{name} {last.get(name)} -> {status}{extra}")
            last[name] = status
        if status == "running" and name not in flagged:
            run_h = (time.time() * 1000 - (s.get("running") or time.time() * 1000)) / 3.6e6
            try:
                lg = k.rpc(EE2, "execution_engine2.get_job_logs", [{"job_id": v["job_id"], "skip_lines": 0}])[0]
                lines = lg.get("lines") or []
                ts = max((l.get("ts") or 0) for l in lines) / 1000 if lines else None
                silent_h = (time.time() - ts) / 3600 if ts else None
            except Exception:
                silent_h = None
            if run_h > 72 or (silent_h is not None and silent_h > 24):
                ev(f"STALL {name}: running {run_h:.0f} h, log silent {silent_h if silent_h is None else round(silent_h)} h (job {v['job_id']})")
                flagged.add(name)
    if all(last.get(n) in ("completed", "error", "terminated") for n in jobs):
        ev("ALLDONE " + json.dumps({n: last[n] for n in jobs}))
        break
    time.sleep(1800)
