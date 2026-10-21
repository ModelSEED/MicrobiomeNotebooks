#!/usr/bin/env python
"""List the genomes / read libraries that failed in each project.

Salterns (ws 155805): the 200 MAGs queued through
KBDatalakeApps/build_genome_datalake_tables. A genome counts as done only if its
"<name>_.dlt" output object exists -- job state is not trusted, because EE2
reports 'running' for jobs that have been dead for days.

Wetland (ws 261700): MEGAHIT assemblies. A read library counts as failed if its
job errored/terminated, or is a zombie (still 'running' but no log output for
--silent-hours), AND no assembly object for it exists in the workspace.

Writes:
    salterns/BERDL_MAGs/failed_genomes.csv
    wetland/failed_assemblies.csv
"""

import argparse
import csv
import json
import os
import re
import sys
from collections import Counter
from datetime import datetime, timezone

_KBUTILLIB_SRC = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "KBUtilLib", "src")
)
if os.path.isdir(_KBUTILLIB_SRC) and _KBUTILLIB_SRC not in sys.path:
    sys.path.insert(0, _KBUTILLIB_SRC)

from kbutillib.kb_job_utils import KBJobUtils  # noqa: E402
from kbutillib.kb_ws_utils import KBWSUtils  # noqa: E402

SALTERNS_WS = 155805
WETLAND_WS = 261700
TS = re.compile(r"(20\d\d-\d\d-\d\d[ T]\d\d:\d\d:\d\d)")


def flatten(states):
    while states and isinstance(states[0], list):
        states = [x for s in states for x in s]
    return states


def last_log_age_h(ee2, job_id, now):
    """Hours since the job's last timestamped log line (None if unknown)."""
    try:
        lines = ee2.get_job_logs({"job_id": job_id}).get("lines", [])
    except Exception:
        return None, None
    for l in reversed(lines):
        m = TS.search(l.get("line", ""))
        if m:
            when = datetime.strptime(m.group(1).replace("T", " "),
                                     "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)
            return (now - when).total_seconds() / 3600, m.group(1)
    return None, None


def salterns_report(ws, ee2, now, ledger_path, out_path, silent_hours):
    ledger = json.load(open(ledger_path))
    objs = ws.list_ws_objects(SALTERNS_WS, include_metadata=False)
    jobs = {j["job_id"]: j for j in ledger["jobs"]}
    states = {j["job_id"]: j for j in flatten(
        ee2.check_jobs({"job_ids": list(jobs), "exclude_fields": ["job_input"]})
        .get("job_states", []))}

    # Classify each job once, then attribute to its genomes.
    verdicts = {}
    for jid, entry in jobs.items():
        st = states.get(jid, {})
        status = st.get("status", "unknown")
        silent, last = (None, None)
        if status == "running":
            silent, last = last_log_age_h(ee2, jid, now)
            if silent is not None and silent >= silent_hours:
                status, reason = "zombie", f"no log output for {silent:.0f}h"
            else:
                reason = "still running"
        elif status == "error":
            reason = st.get("errormsg") or "job error"
        elif status == "queued":
            reason = "not started"
        elif status == "completed":
            reason = "job completed"
        else:
            reason = status
        verdicts[jid] = (status, str(reason)[:200], last)

    rows = []
    for jid, entry in jobs.items():
        status, reason, last = verdicts[jid]
        for name in entry["genomes"]:
            produced = f"{name}_.dlt" in objs
            rows.append({
                "genome": name,
                "genome_ref": entry["genome_refs"][entry["genomes"].index(name)]
                if name in entry["genomes"] else "",
                "batch_output": entry["output"],
                "job_id": jid,
                "job_status": status,
                "output_produced": "yes" if produced else "no",
                "outcome": "succeeded" if produced else (
                    "pending" if status == "queued" else "failed"),
                "reason": "" if produced else reason,
                "last_log": last or "",
            })

    rows.sort(key=lambda r: (r["outcome"] != "failed", r["batch_output"], r["genome"]))
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


def wetland_report(ws, ee2, now, out_path, silent_hours):
    objs = ws.list_ws_objects(WETLAND_WS, include_metadata=False)
    states = flatten(ee2.check_workspace_jobs({"workspace_id": WETLAND_WS})["job_states"])
    mega = [j for j in states
            if "megahit" in str((j.get("job_input") or {}).get("app_id", "")).lower()]

    # Which libraries ever produced an assembly? (a later retry may have worked)
    succeeded_out = set()
    for j in mega:
        p = ((j.get("job_input") or {}).get("params") or [{}])[0]
        if j.get("status") == "completed" and p.get("output_contigset_name"):
            succeeded_out.add(p["output_contigset_name"])

    rows = []
    for j in mega:
        status = j.get("status")
        p = ((j.get("job_input") or {}).get("params") or [{}])[0]
        out_name = p.get("output_contigset_name") or ""
        silent, last = (None, None)
        if status == "running":
            silent, last = last_log_age_h(ee2, j["job_id"], now)
            if silent is not None and silent >= silent_hours:
                status = "zombie"
            else:
                continue                      # genuinely still running
        elif status not in ("error", "terminated"):
            continue                          # completed / queued

        elapsed = ((j.get("finished") or now.timestamp() * 1000) - j["running"]) / 3.6e6 \
            if j.get("running") else 0
        rows.append({
            "read_library_ref": p.get("read_library_ref", ""),
            "output_contigset_name": out_name,
            "job_id": j["job_id"],
            "job_status": status,
            "elapsed_hours": f"{elapsed:.1f}",
            "silent_hours": f"{silent:.0f}" if silent is not None else "",
            "last_log": last or "",
            "assembly_exists": "yes" if out_name in objs else "no",
            "recovered_by_retry": "yes" if (out_name in succeeded_out or out_name in objs)
                                  else "no",
            "errormsg": str(j.get("errormsg") or "")[:160],
        })

    rows.sort(key=lambda r: (r["recovered_by_retry"] == "yes", r["job_status"],
                             r["output_contigset_name"]))
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--silent-hours", type=float, default=24.0,
                    help="a 'running' job with no log output for this long is a zombie")
    ap.add_argument("--ledger", default=os.path.join(here, "BERDL_MAGs", "datalake_jobs.json"))
    ap.add_argument("--salterns-out", default=os.path.join(here, "BERDL_MAGs",
                                                           "failed_genomes.csv"))
    ap.add_argument("--wetland-out", default=os.path.normpath(
        os.path.join(here, "..", "wetland", "failed_assemblies.csv")))
    args = ap.parse_args()

    now = datetime.now(timezone.utc)
    ws = KBWSUtils(kb_version="prod")
    ee2 = KBJobUtils(ws, kb_version="prod").ee2_client

    s = salterns_report(ws, ee2, now, args.ledger, args.salterns_out, args.silent_hours)
    print(f"[salterns] {len(s)} genomes: {dict(Counter(r['outcome'] for r in s))}")
    for out, n in sorted(Counter(r["batch_output"] for r in s
                                 if r["outcome"] == "failed").items()):
        print(f"    {n:3d} failed in {out}")
    print(f"[salterns] -> {args.salterns_out}")

    os.makedirs(os.path.dirname(args.wetland_out), exist_ok=True)
    w = wetland_report(ws, ee2, now, args.wetland_out, args.silent_hours)
    print(f"\n[wetland] {len(w)} failed/stuck assemblies: "
          f"{dict(Counter(r['job_status'] for r in w))}")
    print(f"    recovered by a later retry: "
          f"{sum(1 for r in w if r['recovered_by_retry'] == 'yes')}")
    print(f"    still missing an assembly : "
          f"{sum(1 for r in w if r['recovered_by_retry'] == 'no')}")
    print(f"[wetland] -> {args.wetland_out}")


if __name__ == "__main__":
    main()
