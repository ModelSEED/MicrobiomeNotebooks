#!/usr/bin/env python
"""Run KBase's DRAM app on the 310 dereplicated salterns MAGs, as size-balanced parallel jobs.

Target: narrative "Metabolism 310 Bins" (ws 77580), where the 310 Assemblies live and where the 2026 CheckM and
GTDB-Tk re-runs were written (queue_mag_qc_jobs.py). Input: its AssemblySet `SSFB_Assembly_Set` (77580/312/1).

KBase's DRAM (kb_DRAM/run_kb_dram_annotate, DRAM 1.4.0 at service version 935498b6; KOfam, Pfam, dbCAN, MEROPS,
no UniRef) annotates the MAGs of one job one after another, about 21 min plus 4 min per Mb each, so the 310 MAGs
(1,018 Mb) are split into --parts AssemblySets of nearly equal total length (largest-first greedy) and run as one
job per part. App parameters are the app defaults also used for the wetland MAGs: min_contig_size 2500 (every
contig of these MAGs is >= 2,500 bp, so nothing is dropped), trans_table 11, bitscore 60, rbh_bitscore 350.

Jobs go through EE2; outputs (DRAM-annotated Genomes, one GenomeSet per part, reports) land in ws 77580.
Ledger: prfba_augmentation/dram/mag_jobs.json.

Usage:
    python queue_dram_jobs.py --dry-run      # show the partition and the job params; save/submit nothing
    python queue_dram_jobs.py                # save the part sets (reused if present) and submit missing jobs
    python queue_dram_jobs.py --status       # state, client group, queue/run time, last log line age
    python queue_dram_jobs.py --collect      # download finished reports' files to prfba_augmentation/dram/mags/
"""

import argparse
import json
import os
import sys
from datetime import datetime, timezone

_KBUTILLIB_SRC = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "KBUtilLib", "src")
)
if os.path.isdir(_KBUTILLIB_SRC) and _KBUTILLIB_SRC not in sys.path:
    sys.path.insert(0, _KBUTILLIB_SRC)

from kbutillib.kb_app_runner.nms import NMSSpecCache  # noqa: E402
from kbutillib.kb_app_runner.runner import AppRunner  # noqa: E402
from kbutillib.kb_genome_utils import KBGenomeUtils  # noqa: E402
from kbutillib.kb_job_utils import KBJobUtils  # noqa: E402

import queue_mag_qc_jobs as qc  # noqa: E402  (ledger I/O, job status and report download)

WS = 77580
INPUT_REF = "77580/312/1"  # SSFB_Assembly_Set: the 310 dereplicated Assemblies
APP_ID = "kb_DRAM/run_kb_dram_annotate"
DRAM_PARAMS = dict(min_contig_size=2500, trans_table=11, bitscore=60, rbh_bitscore=350, is_metagenome=0)
SET_NAME = "SSFB_Assembly_Set_DRAM_part{:02d}of{:02d}"
OUT_NAME = "SSFB_DRAM_part{:02d}of{:02d}.GenomeSet"


def balanced_parts(ws, n):
    """Split INPUT_REF's items into n sets of nearly equal total DNA size (largest first, into the lightest set)."""
    items = ws.get_objects2({"objects": [{"ref": INPUT_REF}]})["data"][0]["data"]["items"]
    refs = [it["ref"] for it in items]
    infos = ws.get_object_info3({"objects": [{"ref": r} for r in refs], "includeMetadata": 1})["infos"]
    size = {r: int((info[10] or {}).get("Size") or (info[10] or {}).get("dna_size") or info[9]) for r, info in zip(refs, infos)}
    bins = [[] for _ in range(n)]
    load = [0] * n
    for r in sorted(refs, key=lambda x: -size[x]):
        k = min(range(n), key=lambda i: (load[i], len(bins[i])))
        bins[k].append(r)
        load[k] += size[r]
    by_ref = {it["ref"]: it for it in items}
    return [([by_ref[r] for r in b], sum(size[r] for r in b)) for b in bins], len(refs)


def save_parts(ws, parts, dry_run):
    existing = {o[1]: f"{o[6]}/{o[0]}/{o[4]}" for o in ws.list_objects({"ids": [WS], "type": "KBaseSets.AssemblySet"})}
    n, out = len(parts), []
    for i, (chunk, bp) in enumerate(parts, 1):
        name = SET_NAME.format(i, n)
        if name not in existing and not dry_run:
            info = ws.save_objects({"id": WS, "objects": [{
                "type": "KBaseSets.AssemblySet-1.1", "name": name,
                "data": {"description": f"Part {i} of {n} of SSFB_Assembly_Set ({INPUT_REF}): {len(chunk)} assemblies, "
                                        f"{bp / 1e6:.1f} Mb, size-balanced for DRAM", "items": chunk},
                "provenance": [{"script": "queue_dram_jobs.py", "input_ws_objects": [INPUT_REF],
                                "description": f"size-balanced {n}-way split of SSFB_Assembly_Set for kb_DRAM"}],
            }]})[0]
            existing[name] = f"{info[6]}/{info[0]}/{info[4]}"
        out.append((name, existing.get(name, "(not saved: dry run)"), len(chunk), bp))
    return out


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(here, "prfba_augmentation", "dram"), help="ledger + results dir")
    ap.add_argument("--parts", type=int, default=31, help="number of jobs (default 31, ~10 MAGs each)")
    ap.add_argument("--dry-run", action="store_true", help="show the partition and params; save/submit nothing")
    ap.add_argument("--status", action="store_true", help="poll submitted jobs and exit")
    ap.add_argument("--collect", action="store_true", help="download finished reports' files and exit")
    ap.add_argument("--cancel", action="store_true", help="cancel every unfinished job in the ledger and exit")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    ledger_path = os.path.join(args.out, "mag_jobs.json")
    ledger = qc.load(ledger_path)
    u = KBGenomeUtils(kb_version="prod")
    jobs = KBJobUtils(u, kb_version="prod")
    if args.status:
        return qc.cmd_status(jobs.ee2_client, ledger)
    if args.collect:
        return qc.cmd_collect(u, jobs.ee2_client, ledger, os.path.join(args.out, "mags"))
    if args.cancel:
        states = qc.job_states(jobs.ee2_client, [j["job_id"] for j in ledger["jobs"].values()])
        for key in list(ledger["jobs"]):
            j = ledger["jobs"][key]
            st = states[j["job_id"]]["status"]
            if st in ("completed", "error", "terminated"):
                print(f"[keep] {key}: {st}")
                continue
            jobs.cancel_job(j["job_id"])
            j.update(cancelled_utc=datetime.now(timezone.utc).isoformat(), state_at_cancel=st,
                     archived_reason="cancelled at the user's request (2026-11-04): DRAM2 to be run outside KBase instead")
            ledger.setdefault("cancelled", []).append(ledger["jobs"].pop(key))
            qc.save(ledger_path, ledger)
            print(f"[cancel] {key}: {j['job_id']} (was {st})")
        return None

    ws = u.ws_client()
    parts, n_items = balanced_parts(ws, args.parts)
    assert n_items == 310, f"{INPUT_REF} holds {n_items} assemblies, expected 310"
    saved = save_parts(ws, parts, args.dry_run)
    service_ver = NMSSpecCache().get(APP_ID, None).service_ver
    runner = None if args.dry_run else AppRunner(ws=u, kb_version="prod", job_store=jobs)
    for i, (name, ref, count, bp) in enumerate(saved, 1):
        key = f"dram_p{i:02d}"
        params = dict(assembly_input_ref=ref, output_name=OUT_NAME.format(i, len(saved)),
                      desc=f"DRAM on {name} ({count} MAGs, {bp / 1e6:.1f} Mb)", **DRAM_PARAMS)
        if key in ledger["jobs"]:
            print(f"[skip] {key}: already in ledger as {ledger['jobs'][key]['job_id']}")
        elif args.dry_run:
            print(f"[dry-run] {key}: {name} {count} MAGs {bp / 1e6:.1f} Mb ver={service_ver[:8]} {json.dumps(params)}")
        else:
            handle = runner.run_app(APP_ID, params, WS, meta={"input": ref, "set": name})
            ledger["jobs"][key] = {"app_id": APP_ID, "job_id": handle.job_id, "params": params, "set": name,
                                   "n_assemblies": count, "bp": bp, "service_ver": service_ver,
                                   "submitted_utc": datetime.now(timezone.utc).isoformat()}
            qc.save(ledger_path, ledger)  # persist per job: a crash loses nothing
            print(f"[submit] {key}: {name} ({count} MAGs, {bp / 1e6:.1f} Mb) -> job {handle.job_id}")
    if not args.dry_run:
        print(f"[ledger] {ledger_path}\n[monitor] python {os.path.basename(__file__)} --status")


if __name__ == "__main__":
    main()
