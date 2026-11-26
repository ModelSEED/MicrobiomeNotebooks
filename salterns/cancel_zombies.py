#!/usr/bin/env python
"""Cancel zombie EE2 jobs: still reported 'running' but silent for --silent-hours.

Spans both projects -- salterns datalake (ws 155805) and wetland MEGAHIT
(ws 261700). Deliberately targeted: it only touches jobs that are BOTH 'running'
AND log-silent past the threshold. Queued jobs, completed jobs and jobs that are
still emitting log lines are never cancelled.

It also archives salterns ledger entries whose job has reached a dead state
(cancelled/error/terminated) into ledger["cancelled"], so those genomes stop
counting as "covered by a live job" and become re-queueable by
queue_datalake_jobs.py.

Usage:
    python cancel_zombies.py --dry-run     # list what would be cancelled
    python cancel_zombies.py               # cancel them
"""

import argparse
import json
import os
import re
import sys
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
DEAD_STATES = ("error", "terminated", "cancelled")
TS = re.compile(r"(20\d\d-\d\d-\d\d[ T]\d\d:\d\d:\d\d)")


def flatten(states):
    while states and isinstance(states[0], list):
        states = [x for s in states for x in s]
    return states


def silence_hours(ee2, job_id, now):
    """Hours since the job's last timestamped log line; None if none found."""
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


def find_zombies(ee2, states, now, threshold, label_of):
    out = []
    for j in states:
        if j.get("status") != "running":
            continue
        silent, last = silence_hours(ee2, j["job_id"], now)
        elapsed = (now.timestamp() * 1000 - j["running"]) / 3.6e6 if j.get("running") else 0
        if silent is not None and silent >= threshold:
            out.append({"job_id": j["job_id"], "label": label_of(j),
                        "elapsed_h": elapsed, "silent_h": silent, "last_log": last})
        else:
            print(f"  keeping {j['job_id']} {label_of(j)[:44]} "
                  f"(elapsed {elapsed:.0f}h, silent "
                  f"{'unknown' if silent is None else f'{silent:.0f}h'}) -- below threshold")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--silent-hours", type=float, default=24.0,
                    help="cancel 'running' jobs with no log output for this long")
    ap.add_argument("--ledger", default=os.path.join(here, "BERDL_MAGs", "datalake_jobs.json"))
    ap.add_argument("--dry-run", action="store_true", help="list only; cancel nothing")
    ap.add_argument("--skip-wetland", action="store_true")
    ap.add_argument("--skip-salterns", action="store_true")
    args = ap.parse_args()

    now = datetime.now(timezone.utc)
    ws = KBWSUtils(kb_version="prod")
    jobs_util = KBJobUtils(ws, kb_version="prod")
    ee2 = jobs_util.ee2_client

    ledger = json.load(open(args.ledger))
    byid = {j["job_id"]: j for j in ledger["jobs"]}
    targets = []

    if not args.skip_salterns:
        print(f"[salterns ws {SALTERNS_WS}] checking {len(byid)} ledger jobs")
        st = flatten(ee2.check_jobs({"job_ids": list(byid),
                                     "exclude_fields": ["job_input"]}).get("job_states", []))
        z = find_zombies(ee2, st, now, args.silent_hours,
                         lambda j: byid.get(j["job_id"], {}).get("output", "?"))
        for e in z:
            e["project"] = "salterns"
        targets += z

    if not args.skip_wetland:
        print(f"[wetland ws {WETLAND_WS}] checking MEGAHIT jobs")
        st = flatten(ee2.check_workspace_jobs({"workspace_id": WETLAND_WS})["job_states"])
        mega = [j for j in st
                if "megahit" in str((j.get("job_input") or {}).get("app_id", "")).lower()]
        z = find_zombies(ee2, mega, now, args.silent_hours,
                         lambda j: (((j.get("job_input") or {}).get("params") or [{}])[0]
                                    .get("output_contigset_name") or j["job_id"]))
        for e in z:
            e["project"] = "wetland"
        targets += z

    print(f"\n=== {len(targets)} zombie job(s) to cancel "
          f"(running, silent >= {args.silent_hours:.0f}h) ===")
    for e in sorted(targets, key=lambda e: (e["project"], e["label"])):
        print(f"  {e['project']:9} {e['job_id']} {e['label'][:44]:44} "
              f"elapsed={e['elapsed_h']:6.0f}h silent={e['silent_h']:6.0f}h "
              f"last={e['last_log']}")

    if args.dry_run:
        print("\n[dry-run] nothing cancelled.")
        return
    if not targets:
        print("nothing to cancel.")
        return

    ok = fail = 0
    for e in targets:
        try:
            jobs_util.cancel_job(e["job_id"])
            ok += 1
            print(f"  cancelled {e['job_id']}  {e['label'][:50]}")
        except Exception as exc:
            fail += 1
            print(f"  FAILED {e['job_id']}: {type(exc).__name__}: {exc}", file=sys.stderr)
    print(f"\n[cancel] cancelled={ok} failed={fail}")

    # Ledger hygiene: any salterns entry whose job is now dead should stop
    # blocking its genomes from being re-queued.
    st = flatten(ee2.check_jobs({"job_ids": list(byid),
                                 "exclude_fields": ["job_input"]}).get("job_states", []))
    state_of = {j["job_id"]: j.get("status") for j in st}
    live, archived = [], []
    for entry in ledger["jobs"]:
        s = state_of.get(entry["job_id"])
        if s in DEAD_STATES:
            entry["final_state"] = s
            entry["archived_utc"] = now.isoformat()
            archived.append(entry)
        else:
            live.append(entry)
    if archived:
        ledger.setdefault("cancelled", []).extend(archived)
        ledger["jobs"] = live
        tmp = args.ledger + ".part"
        with open(tmp, "w") as fh:
            json.dump(ledger, fh, indent=2)
        os.replace(tmp, args.ledger)
        freed = sum(len(e.get("genomes") or []) for e in archived)
        print(f"[ledger] archived {len(archived)} dead entr(ies) covering {freed} genomes; "
              f"{len(live)} still live -> those genomes are re-queueable")


if __name__ == "__main__":
    main()
