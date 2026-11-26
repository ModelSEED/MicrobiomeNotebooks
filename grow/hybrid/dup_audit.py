"""Read-only duplicate audit for the GROW per-sample (ws 276048) and sitecov (ws 277632) drivers.

Flags: (1) a task with more than one live attempt (KBase and/or local); (2) KBase jobs live in EE2 that the
ledger does not know (orphans); (3) KBase jobs live in EE2 that the ledger records as finished/cancelled.
Exit code 1 if anything is flagged. Never cancels anything itself.
usage: ~/Documents/py_venv/bin/python grow/hybrid/dup_audit.py
"""
import collections, json, os, sys

sys.path.insert(0, os.path.expanduser("~/Documents/KBUtilLib/src"))
from kbutillib.kb_genome_utils import KBGenomeUtils  # noqa: E402
from kbutillib.kb_job_utils import KBJobUtils  # noqa: E402

GROW = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIVE = {"created", "estimating", "queued", "running"}
RUNS = [("per_sample", "assembly_binning_2026/ledger.json", 276048),
        ("sitecov", "multisample_binning_2026/ledger.json", 277632)]

ee2 = KBJobUtils(KBGenomeUtils(kb_version="prod"), kb_version="prod").ee2_client
bad = 0
for name, rel, ws in RUNS:
    tasks = json.load(open(os.path.join(GROW, rel)))["tasks"]
    known, live = {}, collections.defaultdict(list)
    for key, t in tasks.items():
        for a in t.get("attempts", []):
            if a.get("job_id"):
                known[a["job_id"]] = (key, a.get("status"))
            if a.get("status") in LIVE or (a.get("local") and a.get("status") in ("running", "queued", "started")):
                live[key].append(f"local:{a.get('host', 'FM')}" if a.get("local") else f"kbase:{a.get('job_id')}")
    res = ee2.check_workspace_jobs({"workspace_id": ws, "return_list": 1, "exclude_fields": ["job_input", "job_output"]})
    ee2_live = {j["job_id"]: j["status"] for j in res.get("job_states", []) if j["status"] in LIVE}
    dup = {k: v for k, v in live.items() if len(v) > 1}
    orphans = [j for j in ee2_live if j not in known]
    stale = {j: known[j] for j in ee2_live if j in known and known[j][1] not in LIVE}
    bad += bool(dup) + bool(orphans) + bool(stale)
    print(f"{name} (ws {ws}): EE2 live {len(ee2_live)}, ledger live tasks {len(live)}")
    print(f"  duplicate live attempts: {dup or 'none'}")
    print(f"  EE2 orphans (not in ledger): {orphans or 'none'}")
    print(f"  EE2 live but ledger says finished: {stale or 'none'}")
sys.exit(1 if bad else 0)
