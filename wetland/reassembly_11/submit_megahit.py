#!/usr/bin/env python3
"""Retry MEGAHIT on the 11 wetland read libraries that were never assembled (per Andrew, 2026-11-15).

Each had one MEGAHIT attempt in KBase narrative 261700: 10 went log-silent on 2026-07-30 and were terminated 08-14, and
SRR13021855 errored ("Unexpected Job Error"). The retry uses exactly the configuration of all 770 completed assemblies:
MEGAHIT.run_megahit, service_ver e2e6d557da8a92d8b9a6cb4893a2a7293a3068bf, min_contig_length 2000, other params default.
Skips a library whose output Assembly already exists or that has a queued/running MEGAHIT job.
Dry run unless --apply. Job ids -> jobs.json. Usage: KB_AUTH_TOKEN=... python -I -B submit_megahit.py [--apply]
"""
import csv, json, os, sys, time
sys.path.insert(0, "/home/freiburger/gtdbtk_kbase")
import kbase_mags as k

EE2 = "https://kbase.us/services/ee2"
WSID, WSNAME = 261700, "freiburgermsu:narrative_1780514279067"
SERVICE_VER = "e2e6d557da8a92d8b9a6cb4893a2a7293a3068bf"
HERE = os.path.dirname(os.path.abspath(__file__))
apply = "--apply" in sys.argv
rows = [r for r in csv.DictReader(open("/home/freiburger/Documents/MicrobiomeNotebooks/wetland/failed_assemblies.csv"))
        if r["recovered_by_retry"] == "no"]
assert len(rows) == 11, len(rows)
existing = {i[1] for i in k.rpc(k.WS, "Workspace.list_objects", [{"ids": [WSID], "type": "KBaseGenomeAnnotations.Assembly"}])[0]}
res = k.rpc(EE2, "execution_engine2.check_workspace_jobs", [{"workspace_id": WSID, "return_list": 1,
                                                             "exclude_fields": ["job_output", "condor_job_ads"]}])[0]
states = res["job_states"] if isinstance(res, dict) else res
active = {(j.get("job_input") or {}).get("params", [{}])[0].get("read_library_ref") for j in states
          if j.get("status") in ("created", "estimating", "queued", "running")
          and (j.get("job_input") or {}).get("app_id") == "MEGAHIT/run_megahit"}
ledger = json.load(open(os.path.join(HERE, "jobs.json"))) if os.path.exists(os.path.join(HERE, "jobs.json")) else {}
for r in rows:
    ref, out = r["read_library_ref"], r["output_contigset_name"]
    info = k.rpc(k.WS, "Workspace.get_object_info3", [{"objects": [{"ref": ref}]}])[0]["infos"][0]
    if out in existing or ref in active:
        print(f"skip {out}: {'assembly exists' if out in existing else 'MEGAHIT job active'}")
        continue
    job = {"method": "MEGAHIT.run_megahit", "app_id": "MEGAHIT/run_megahit", "service_ver": SERVICE_VER, "wsid": WSID,
           "source_ws_objects": [ref],
           "params": [{"workspace_name": WSNAME, "read_library_ref": ref, "output_contigset_name": out,
                       "min_contig_length": 2000, "megahit_parameter_preset": None, "min_count": None, "k_min": None,
                       "k_max": None, "k_step": None, "k_list": []}]}
    if not apply:
        print(f"[dry run] {out} <- {ref} ({info[1]}, {info[2]})")
        continue
    jid = k.rpc(EE2, "execution_engine2.run_job", [job], tries=1)[0]
    ledger[out] = {"job_id": jid, "read_library_ref": ref, "submitted": time.strftime("%F %T"), "previous_job": r["job_id"]}
    json.dump(ledger, open(os.path.join(HERE, "jobs.json"), "w"), indent=1)
    print(f"+ {out} <- {ref}: job {jid}")
