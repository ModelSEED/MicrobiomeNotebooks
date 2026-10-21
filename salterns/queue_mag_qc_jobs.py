#!/usr/bin/env python
"""Re-run GTDB-Tk and CheckM on the 310 dereplicated salterns MAGs in KBase.

Target: narrative "Metabolism 310 Bins" (ws 77580, cliff_bueno), where the 310 Assemblies
live; input is its AssemblySet `SSFB_Assembly_Set` (77580/312/1), verified to hold exactly
the 310. Parameters mirror John-Paul Faria's 2023 runs on the 140 unretained MAGs (ws 163264),
which completed in ~4 h (GTDB-Tk) and ~2.5 h (CheckM) on the `njs` client group:

    kb_gtdbtk/run_kb_gtdbtk_classify_wf   GTDB-Tk 2.3.2, GTDB r214 (newest KBase offers)
    kb_Msuite/run_checkM_lineage_wf       CheckM 1.0.18 lineage_wf, reduced tree, plots saved

GTDB-Tk writes its taxonomy back onto the input Assemblies, so each gets a new version in
ws 77580; earlier versions are kept and every downstream assembly_ref is version-pinned.

Jobs go through EE2 (outputs and reports land in the narrative's data panel; no app cells
are added to the Narrative). Ledger: MAG_QC/jobs.json.

Usage:
    python queue_mag_qc_jobs.py --dry-run     # print the job params
    python queue_mag_qc_jobs.py               # submit any app not already in the ledger
    python queue_mag_qc_jobs.py --status      # state, client group, timings, log silence
    python queue_mag_qc_jobs.py --collect     # download finished reports' files to MAG_QC/<app>/
    python queue_mag_qc_jobs.py --gtdbtk-parts 3   # GTDB-Tk on 3 ~100-member AssemblySets
    python queue_mag_qc_jobs.py --gtdbtk-parts 3 --tag beta --replace
                                              # cancel unfinished part jobs; resubmit on beta

The whole-set GTDB-Tk job died when KBase's JobRunner crashed decoding GTDB-Tk's UTF-8
progress-bar output ("No output generated"); --gtdbtk-parts archives that ledger entry,
saves SSFB_Assembly_Set_part<i>of<n> sets in ws 77580 (reused on re-runs), and submits one
GTDB-Tk job per part so each job emits less of the log that triggers the crash.
"""

import argparse
import json
import os
import re
import sys
from datetime import datetime, timezone

import requests

_KBUTILLIB_SRC = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "KBUtilLib", "src")
)
if os.path.isdir(_KBUTILLIB_SRC) and _KBUTILLIB_SRC not in sys.path:
    sys.path.insert(0, _KBUTILLIB_SRC)

from kbutillib.kb_app_runner.nms import NMSSpecCache  # noqa: E402
from kbutillib.kb_app_runner.runner import AppRunner  # noqa: E402
from kbutillib.kb_genome_utils import KBGenomeUtils  # noqa: E402
from kbutillib.kb_job_utils import KBJobUtils  # noqa: E402

WS = 77580
WS_NAME = "cliff_bueno:narrative_1605678421558"
INPUT_REF = "77580/312/1"  # SSFB_Assembly_Set: the 310 dereplicated Assemblies
APPS = {
    "gtdbtk": ("kb_gtdbtk/run_kb_gtdbtk_classify_wf", [{
        "workspace_id": WS, "input_object_ref": INPUT_REF, "output_tree_basename": "GTDB_Tree",
        "copy_proximals": "0", "save_trees": "0", "min_perc_aa": 10, "db_ver": "214",
        "keep_intermediates": "0", "overwrite_tax": "0", "dendrogram_report": "1"}]),
    "checkm": ("kb_Msuite/run_checkM_lineage_wf", [{
        "workspace_name": WS_NAME, "input_ref": INPUT_REF, "reduced_tree": "1",
        "save_plots_dir": "1", "threads": "4"}]),
}
LOG_TS = re.compile(r"(20\d\d-\d\d-\d\d[ T]\d\d:\d\d:\d\d)")


def utc(ms):
    return datetime.fromtimestamp(ms / 1000, timezone.utc) if ms else None


def load(path):
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return {"workspace": WS, "input_ref": INPUT_REF, "jobs": {}}


def save(path, ledger):
    tmp = path + ".part"
    with open(tmp, "w") as fh:
        json.dump(ledger, fh, indent=2)
    os.replace(tmp, path)


def job_states(ee2, ids):
    res = ee2.check_jobs({"job_ids": ids})
    return {s["job_id"]: s for s in res.get("job_states", res)}


def last_log(ee2, job_id):
    """Timestamp of the job's last log line (None if no log yet).

    EE2 stamps each line with `ts` (epoch ms); fall back to a timestamp in the text."""
    try:
        lines = ee2.get_job_logs({"job_id": job_id}).get("lines", [])
    except Exception:
        return None
    if lines and lines[-1].get("ts"):
        return utc(lines[-1]["ts"])
    for line in reversed(lines):
        m = LOG_TS.search(line.get("line", ""))
        if m:
            return datetime.strptime(m[1].replace("T", " "), "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)
    return None


def cmd_status(ee2, ledger):
    if not ledger["jobs"]:
        print("[status] nothing submitted yet.")
        return
    now = datetime.now(timezone.utc)
    states = job_states(ee2, [j["job_id"] for j in ledger["jobs"].values()])
    for key, j in ledger["jobs"].items():
        s = states[j["job_id"]]
        q, r, f = utc(s.get("queued")), utc(s.get("running")), utc(s.get("finished"))
        cg = (s.get("job_input") or {}).get("requirements", {}).get("clientgroup", "?")
        line = f"{key:7s} {j['job_id']} {s['status']:9s} clientgroup={cg}"
        if q:
            line += f"  queued {((r or now) - q).total_seconds() / 3600:.1f} h"
        if r:
            line += f"  ran {((f or now) - r).total_seconds() / 3600:.1f} h"
        if s["status"] == "running":
            ll = last_log(ee2, j["job_id"])
            line += f"  last log {((now - ll).total_seconds() / 3600):.1f} h ago" if ll else "  no log yet"
        if s["status"] == "error":
            line += f"\n        error: {(s.get('errormsg') or '')[:300]}"
        print(line)


def cmd_collect(u, ee2, ledger, out_dir):
    """Download every file the finished reports link to (plus GTDB-Tk's output archive)."""
    token = u.get_token(namespace="kbase")
    states = job_states(ee2, [j["job_id"] for j in ledger["jobs"].values()])
    for key, j in ledger["jobs"].items():
        s = states[j["job_id"]]
        if s["status"] != "completed":
            print(f"[collect] {key}: {s['status']} -- skipped")
            continue
        result = ((s.get("job_output") or {}).get("result") or [{}])[0]
        dest = os.path.join(out_dir, key)
        os.makedirs(dest, exist_ok=True)
        report = u.ws_client().get_objects2({"objects": [{"ref": result["report_ref"]}]})["data"][0]["data"]
        with open(os.path.join(dest, "report.json"), "w") as fh:
            json.dump({"result": result, "report": report}, fh, indent=2)
        links = (report.get("file_links") or []) + (report.get("html_links") or [])
        if result.get("archive_shock_id"):
            links.append({"name": f"{key}_archive", "shock_id": result["archive_shock_id"]})
        for link in links:
            name = os.path.basename(link.get("name") or link.get("label") or "file")
            path = os.path.join(dest, name)
            if os.path.exists(path):
                continue
            if link.get("handle"):
                got = u.download_blob_file(link["handle"], path)
            else:
                node = link.get("shock_id") or link.get("URL", "").rstrip("/").split("/")[-1]
                url = f"{u.shock_url}/node/{node}?download_raw"
                with requests.get(url, headers={"Authorization": "OAuth " + token}, stream=True) as r:
                    got = r.ok
                    if r.ok:
                        with open(path, "wb") as fh:
                            for chunk in r.iter_content(1 << 20):
                                fh.write(chunk)
            print(f"[collect] {key}: {'ok  ' if got else 'FAIL'} {name}")


def make_partitions(ws, n, dry_run=False):
    """Split INPUT_REF's items into n contiguous AssemblySets in WS; reuse any already saved."""
    items = ws.get_objects2({"objects": [{"ref": INPUT_REF}]})["data"][0]["data"]["items"]
    existing = {o[1]: f"{o[6]}/{o[0]}/{o[4]}"
                for o in ws.list_objects({"ids": [WS], "type": "KBaseSets.AssemblySet"})}
    size = -(-len(items) // n)
    parts = []
    for i in range(n):
        name = f"SSFB_Assembly_Set_part{i + 1}of{n}"
        chunk = items[i * size:(i + 1) * size]
        if name not in existing and not dry_run:
            info = ws.save_objects({"id": WS, "objects": [{
                "type": "KBaseSets.AssemblySet-1.1", "name": name,
                "data": {"description": f"Part {i + 1} of {n} of SSFB_Assembly_Set ({INPUT_REF}): "
                                        f"{len(chunk)} assemblies, split for GTDB-Tk", "items": chunk},
                "provenance": [{"script": "queue_mag_qc_jobs.py", "input_ws_objects": [INPUT_REF],
                                "description": f"contiguous {n}-way split of SSFB_Assembly_Set"}],
            }]})[0]
            existing[name] = f"{info[6]}/{info[0]}/{info[4]}"
        parts.append((name, existing.get(name, "(not saved: dry run)"), len(chunk)))
    return parts


def cmd_gtdbtk_parts(u, jobs, ledger, ledger_path, n, dry_run, tag=None, replace=False):
    if "gtdbtk" in ledger["jobs"] and not dry_run:
        failed = ledger["jobs"].pop("gtdbtk")
        failed["archived_reason"] = "KBase JobRunner UTF-8 log-decode crash; re-run as partitions"
        ledger.setdefault("failed", []).append(failed)
        save(ledger_path, ledger)
    if replace:
        live = [k for k in ledger["jobs"] if k.startswith("gtdbtk_p")]
        states = job_states(jobs.ee2_client, [ledger["jobs"][k]["job_id"] for k in live]) if live else {}
        for k in live:
            j = ledger["jobs"][k]
            if states[j["job_id"]]["status"] == "completed":
                continue
            if dry_run:
                print(f"[dry-run] would cancel {k} ({j['job_id']}, {states[j['job_id']]['status']})")
                continue
            jobs.cancel_job(j["job_id"])
            j.update(cancelled_utc=datetime.now(timezone.utc).isoformat(),
                     archived_reason=f"replaced by a resubmission on tag={tag or 'release'}")
            ledger.setdefault("replaced", []).append(ledger["jobs"].pop(k))
            save(ledger_path, ledger)
            print(f"[cancel] {k}: {j['job_id']}")
    parts = make_partitions(u.ws_client(), n, dry_run)
    runner = AppRunner(ws=u, kb_version="prod", job_store=jobs)
    app_id, base = APPS["gtdbtk"]
    service_ver = NMSSpecCache().get(app_id, tag).service_ver
    for i, (name, ref, count) in enumerate(parts, 1):
        key = f"gtdbtk_p{i}"
        params = [{**base[0], "input_object_ref": ref}]
        if key in ledger["jobs"]:
            print(f"[skip] {key}: already in ledger as {ledger['jobs'][key]['job_id']}")
        elif dry_run:
            print(f"[dry-run] {key}: {name} ({count} assemblies) ref={ref} tag={tag or 'release'} ver={service_ver[:8]}")
        else:
            handle = runner.run_app(app_id, params, WS, tag=tag, meta={"input": ref, "set": name})
            ledger["jobs"][key] = {"app_id": app_id, "job_id": handle.job_id, "params": params,
                                   "set": name, "n_assemblies": count, "tag": tag or "release",
                                   "service_ver": service_ver,
                                   "submitted_utc": datetime.now(timezone.utc).isoformat()}
            save(ledger_path, ledger)
            print(f"[submit] {key}: {name} ({count} assemblies) tag={tag or 'release'} "
                  f"ver={service_ver[:8]} -> job {handle.job_id}")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(here, "MAG_QC"), help="ledger + results dir (default ./MAG_QC)")
    ap.add_argument("--dry-run", action="store_true", help="print params; submit nothing")
    ap.add_argument("--status", action="store_true", help="poll submitted jobs and exit")
    ap.add_argument("--collect", action="store_true", help="download finished reports' files and exit")
    ap.add_argument("--gtdbtk-parts", type=int, metavar="N",
                    help="run GTDB-Tk as N jobs on N partition AssemblySets, then exit")
    ap.add_argument("--tag", choices=("release", "beta", "dev"), default="release",
                    help="app release tag for GTDB-Tk partition jobs (default release)")
    ap.add_argument("--replace", action="store_true",
                    help="with --gtdbtk-parts: cancel unfinished partition jobs and resubmit them")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    ledger_path = os.path.join(args.out, "jobs.json")
    ledger = load(ledger_path)
    u = KBGenomeUtils(kb_version="prod")
    jobs = KBJobUtils(u, kb_version="prod")
    if args.status:
        return cmd_status(jobs.ee2_client, ledger)
    if args.collect:
        return cmd_collect(u, jobs.ee2_client, ledger, args.out)
    if args.gtdbtk_parts:
        tag = None if args.tag == "release" else args.tag
        return cmd_gtdbtk_parts(u, jobs, ledger, ledger_path, args.gtdbtk_parts, args.dry_run,
                                tag=tag, replace=args.replace)

    runner = AppRunner(ws=u, kb_version="prod", job_store=jobs)
    for key, (app_id, params) in APPS.items():
        if key in ledger["jobs"]:
            print(f"[skip] {key}: already in ledger as {ledger['jobs'][key]['job_id']}")
            continue
        if args.dry_run:
            print(f"[dry-run] {key}: {app_id}\n          {json.dumps(params)}")
            continue
        handle = runner.run_app(app_id, params, WS, meta={"input": INPUT_REF, "set": "SSFB_Assembly_Set"})
        ledger["jobs"][key] = {"app_id": app_id, "job_id": handle.job_id, "params": params,
                               "submitted_utc": datetime.now(timezone.utc).isoformat()}
        save(ledger_path, ledger)  # persist per job: a crash loses nothing
        print(f"[submit] {key}: {app_id} -> job {handle.job_id}")
    if not args.dry_run:
        print(f"[ledger] {ledger_path}\n[monitor] python {os.path.basename(__file__)} --status")


if __name__ == "__main__":
    main()
