#!/usr/bin/env python
"""Queue KBDatalakeApps/build_genome_datalake_tables (beta) for the genomes in
narrative 155805 that have NOT been through it yet -- one EE2 job per genome.

Companion to download_berdl_mags.py. That script pulled the 110 genomes updated
in the last two weeks; those are the *outputs* of earlier runs of this app (they
carry the "_.dlt" suffix). This script finds the original genomes that still have
no ".dlt" counterpart and submits one job each.

Selection (both conditions must hold):
  * the genome's latest version was saved more than --days ago, and
  * no object named "<genome_name>_<suffix>" exists in the workspace
    -- that name is what the app produces, so its presence means "already done".

Parameters replicate the earlier runs exactly, read from the provenance of the
objects the app already produced:
    suffix=".dlt"  export_databases=1  everything else 0
    output="<genome_name>.dlt"    (KBaseFBA.GenomeDataLakeTables)
Jobs are BATCHED (--batch-size, default 22 genomes per job). Every job of this
app reserves an entire clientgroup='extreme' node (32 CPU, 200 GB RAM) no matter
how many genomes it carries, so batching is simply good practice. The completed
2026-07-17 runs averaged 19 genomes in 0.92 h (~6 min/genome), which is what a
batch of ~20 reproduces. Do NOT collapse everything into one job: the
26-GenomeSet single-job attempts (May-Jul 2026) errored or hung for 197 h and
1443 h before being terminated.

Batch size does NOT control scheduling latency. This user gets ~10 concurrent EE2
jobs cluster-wide, and if something else (e.g. a MEGAHIT bigmem backlog in
another workspace) already holds those slots, these jobs wait regardless of
whether there are 10 of them or 200. Check with
ee2.check_workspace_jobs({'workspace_id': N}) across the user's workspaces before
concluding a job is misconfigured -- 'queued' for days usually means the cap is
full, not that anything is wrong.

Submitted job ids are appended to BERDL_MAGs/datalake_jobs.json, so re-running
never double-submits a genome. KBUtilLib also records every job in its own
SQLite ledger at ~/.kbjobs/kbjobs.db.

Usage:
    python queue_datalake_jobs.py --dry-run        # show what would be queued
    python queue_datalake_jobs.py --limit 1        # submit one, validate it
    python queue_datalake_jobs.py                  # submit the rest
    python queue_datalake_jobs.py --status         # poll submitted jobs
"""

import argparse
import json
import os
import sys
import time
from collections import Counter
from datetime import datetime, timedelta, timezone

_KBUTILLIB_SRC = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "KBUtilLib", "src")
)
if os.path.isdir(_KBUTILLIB_SRC) and _KBUTILLIB_SRC not in sys.path:
    sys.path.insert(0, _KBUTILLIB_SRC)

from kbutillib.kb_app_runner.runner import AppRunner  # noqa: E402
from kbutillib.kb_job_utils import KBJobUtils  # noqa: E402
from kbutillib.kb_ws_utils import KBWSUtils  # noqa: E402

APP_ID = "KBDatalakeApps/build_genome_datalake_tables"
APP_TAG = "beta"
GENOME_TYPE_PREFIX = "KBaseGenomes.Genome"
SUFFIX = ".dlt"


def parse_ws_date(ts: str) -> datetime:
    for fmt in ("%Y-%m-%dT%H:%M:%S%z", "%Y-%m-%dT%H:%M:%S.%f%z"):
        try:
            return datetime.strptime(ts, fmt)
        except ValueError:
            continue
    raise ValueError(f"unrecognized workspace timestamp: {ts!r}")


def load_ledger(path: str) -> dict:
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return {"app_id": APP_ID, "tag": APP_TAG, "suffix": SUFFIX, "jobs": []}


def save_ledger(path: str, ledger: dict) -> None:
    tmp = path + ".part"
    with open(tmp, "w") as fh:
        json.dump(ledger, fh, indent=2)
    os.replace(tmp, path)


def build_params(ws: int, batch, suffix: str, output: str) -> dict:
    """UI-shape params; AppRunner remaps them through the NMS input mapping.

    *batch* is a list of workspace object-info tuples processed by one job.
    Every job reserves a whole clientgroup='extreme' node (32 CPU / 200 GB),
    so batch size is what controls demand on the scarcest resource class --
    not job count.
    """
    return {
        "input_refs": [f"{ws}/{o[0]}/{o[4]}" for o in batch],
        "suffix": suffix,
        "output": output,
        "save_models": 0,
        "export_genome_data": 0,
        "export_pangenome_data": 0,
        "export_databases": 1,
        "export_folder_models": 0,
        "export_folder_phenotypes": 0,
    }


def ledger_genomes(ledger: dict) -> set:
    """Genome names already covered by a live (non-cancelled) job."""
    out = set()
    for j in ledger["jobs"]:
        out.update(j.get("genomes") or ([j["genome"]] if j.get("genome") else []))
    return out


def cmd_cancel(jobs_util, ledger, path) -> None:
    """Cancel every live job in the ledger and archive its entry."""
    live = ledger["jobs"]
    if not live:
        print("[cancel] ledger has no live jobs.")
        return
    print(f"[cancel] cancelling {len(live)} job(s)...")
    ok = fail = 0
    for j in live:
        try:
            jobs_util.cancel_job(j["job_id"])
            ok += 1
        except Exception as e:
            fail += 1
            print(f"  FAILED {j['job_id']}: {type(e).__name__}: {e}", file=sys.stderr)
        j["cancelled_utc"] = datetime.now(timezone.utc).isoformat()
    ledger.setdefault("cancelled", []).extend(live)
    ledger["jobs"] = []
    save_ledger(path, ledger)
    print(f"[cancel] cancelled={ok} failed={fail}; ledger cleared "
          f"({len(ledger['cancelled'])} archived)")


def cmd_status(jobs_util, ledger) -> None:
    ids = [j["job_id"] for j in ledger["jobs"]]
    if not ids:
        print("[status] no jobs submitted yet.")
        return
    print(f"[status] polling {len(ids)} jobs...")
    records = jobs_util.check_jobs(ids)
    by_state = Counter(str(r.state.value if hasattr(r.state, "value") else r.state)
                       for r in records.values())
    for state, n in by_state.most_common():
        print(f"  {n:5d}  {state}")
    byid = {j["job_id"]: j for j in ledger["jobs"]}
    bad = [(jid, r) for jid, r in records.items()
           if str(getattr(r.state, "value", r.state)) in ("error", "terminated")]
    if bad:
        print(f"\n[status] {len(bad)} job(s) not successful:")
        for jid, r in bad[:20]:
            e = byid.get(jid, {})
            print(f"  {jid}  {e.get('output', e.get('genome', '?'))[:55]} "
                  f"({e.get('n_genomes', 1)} genomes)")
            print(f"      {r.error_message or (r.ee2_raw or {}).get('errormsg') or '(no message)'}"[:300])
        if len(bad) > 20:
            print(f"  ... and {len(bad) - 20} more")
        print("\n  Re-queue the failures by removing them from the ledger and re-running:")
        print("    python queue_datalake_jobs.py            # skips genomes still in the ledger")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--ws", type=int, default=155805, help="workspace id (default 155805)")
    ap.add_argument("--days", type=int, default=14,
                    help="genomes saved within this many days count as recent (default 14)")
    ap.add_argument("--suffix", default=SUFFIX, help=f"app output suffix (default {SUFFIX!r})")
    ap.add_argument("--ledger", default=os.path.join(here, "BERDL_MAGs", "datalake_jobs.json"),
                    help="submitted-job ledger (default BERDL_MAGs/datalake_jobs.json)")
    ap.add_argument("--limit", type=int, default=None, help="submit at most N jobs")
    ap.add_argument("--batch-size", type=int, default=22,
                    help="genomes per job (default 22). Each job holds a whole "
                         "32-CPU/200GB 'extreme' node, so this -- not job count -- "
                         "sets cluster demand. ~20-70 matches the runs that completed.")
    ap.add_argument("--batch-prefix", default="BERDL_unprocessed",
                    help="output object name prefix for batch jobs")
    ap.add_argument("--sleep", type=float, default=0.25,
                    help="seconds between submissions (default 0.25)")
    ap.add_argument("--dry-run", action="store_true", help="select and print; submit nothing")
    ap.add_argument("--status", action="store_true", help="poll already-submitted jobs and exit")
    ap.add_argument("--cancel", action="store_true",
                    help="cancel every live job in the ledger, archive them, and exit")
    ap.add_argument("--include-done", action="store_true",
                    help="also (re)submit genomes that already have a suffixed output")
    args = ap.parse_args()

    ws_util = KBWSUtils(kb_version="prod")
    if not ws_util.get_token(namespace="kbase"):
        sys.exit("[auth] no KBase token -- put one in ~/.kbase/token or set $KB_AUTH_TOKEN.")
    jobs_util = KBJobUtils(ws_util, kb_version="prod")
    ledger = load_ledger(args.ledger)

    if args.status:
        cmd_status(jobs_util, ledger)
        return

    if args.cancel:
        cmd_cancel(jobs_util, ledger, args.ledger)
        return

    cutoff = datetime.now(timezone.utc) - timedelta(days=args.days)
    all_objs = ws_util.list_ws_objects(args.ws, include_metadata=False)
    names = set(all_objs)
    genomes = {n: o for n, o in all_objs.items() if o[2].startswith(GENOME_TYPE_PREFIX)}
    older = {n: o for n, o in genomes.items() if parse_ws_date(o[3]) < cutoff}
    # Never feed the app its own output back in. A genome named "<x>_<suffix>"
    # IS a previous run's product; it has no "<x>_<suffix>_<suffix>" of its own,
    # so without this guard it looks unprocessed and gets re-queued forever --
    # which only shows up once those outputs age past the --days window.
    outputs = {n for n in older if n.endswith(f"_{args.suffix}")}
    # The app writes "<input name>_<suffix>"; that object existing means done.
    done = {n for n in older if f"{n}_{args.suffix}" in names}
    todo = {n: o for n, o in older.items()
            if n not in outputs and (args.include_done or n not in done)}

    already_submitted = ledger_genomes(ledger)
    todo = {n: o for n, o in todo.items() if n not in already_submitted}

    ordered = sorted(todo.values(), key=lambda o: o[0])
    size = max(1, args.batch_size)
    batches = [ordered[i : i + size] for i in range(0, len(ordered), size)]
    if args.limit:
        batches = batches[: args.limit]

    print(f"[ws] {args.ws}: {len(all_objs)} objects, {len(genomes)} genomes")
    print(f"[select] {len(genomes) - len(older)} recent (< {args.days}d, skipped as up to date)")
    print(f"[select] {len(older)} older; {len(done)} already have '<name>_{args.suffix}' output; "
          f"{len(outputs)} are themselves '_{args.suffix}' outputs (excluded)")
    print(f"[select] {len(already_submitted)} already covered by a live job")
    print(f"[select] -> {len(ordered)} genome(s) in {len(batches)} job(s) "
          f"of up to {size} genomes each")

    # Batch numbering restarts at 1 every run, so on a resubmission it would
    # collide with outputs an earlier run already produced (overwriting a good
    # result with a retry's) or with names still-queued jobs will claim. Skip
    # any name already taken in the workspace or promised by a ledger entry.
    taken = set(names)
    taken |= {j.get("output") for j in ledger["jobs"]}
    taken |= {j.get("output") for j in ledger.get("cancelled", [])}

    def allocate_names(count):
        out, i = [], 1
        while len(out) < count:
            cand = f"{args.batch_prefix}_batch{i:02d}{args.suffix}"
            if cand not in taken:
                out.append(cand)
                taken.add(cand)
            i += 1
        return out

    batch_names = allocate_names(len(batches))

    def out_name(i):
        return batch_names[i - 1]

    if args.dry_run:
        for i, b in enumerate(batches, 1):
            print(f"  job {i:02d}: {len(b):>3} genomes -> output={out_name(i)}")
            print(f"           first={b[0][1][:58]}")
        if batches:
            print("\n[dry-run] example UI params (job 01):")
            p = build_params(args.ws, batches[0], args.suffix, out_name(1))
            p["input_refs"] = p["input_refs"][:3] + [f"... {len(batches[0])} total"]
            print(json.dumps(p, indent=2))
        print("[dry-run] nothing submitted.")
        return

    if not batches:
        print("[done] nothing to submit.")
        return

    runner = AppRunner(ws=ws_util, kb_version="prod", job_store=jobs_util)
    submitted = failed = 0
    for i, batch in enumerate(batches, 1):
        names = [o[1] for o in batch]
        output = out_name(i)
        params = build_params(args.ws, batch, args.suffix, output)
        try:
            handle = runner.run_app(
                APP_ID, params, args.ws, tag=APP_TAG,
                meta={"n_genomes": len(batch), "output": output},
            )
            ledger["jobs"].append({
                "job_id": handle.job_id, "genomes": names,
                "genome_refs": params["input_refs"], "output": output,
                "n_genomes": len(batch),
                "submitted_utc": datetime.now(timezone.utc).isoformat(),
            })
            save_ledger(args.ledger, ledger)  # persist per job: a crash loses nothing
            submitted += 1
            print(f"[{i}/{len(batches)}] queued {handle.job_id}  "
                  f"{len(batch):>3} genomes -> {output}")
        except Exception as e:
            failed += 1
            print(f"[{i}/{len(batches)}] FAILED {output}: {type(e).__name__}: {e}",
                  file=sys.stderr)
        if args.sleep:
            time.sleep(args.sleep)

    print(f"\n[done] submitted={submitted} failed={failed} "
          f"covering {sum(len(b) for b in batches)} genomes")
    print(f"[ledger] {args.ledger}  (total {len(ledger['jobs'])} jobs)")
    print(f"[monitor] python {os.path.basename(__file__)} --status")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
