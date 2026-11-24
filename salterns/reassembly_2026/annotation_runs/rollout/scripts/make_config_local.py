#!/usr/bin/env python3
"""Configs for the part of the rollout that runs on frankensteinsmonster (Andrew: R2A 2026-11-16, R2 2026-11-17), with the
same driver.py:

  local_config.json        7 R2A + 8 R2 samples (<= 9 Gbp) -> ~/Documents/salterns_rollout_local/rollout/scripts/config.json
  local_smoke_config.json  1 tiny sample on 1M pairs of two R2A water libraries, pushed to poplar's smoke_import root
  smoke_import_config.json poplar side of that smoke: the same sample as a staged import (root rollout/smoke_import)

Task commands and parameters are the poplar ones (make_config.py); only tool paths, thread counts and the resource
policy of this 64-thread / 251 GB shared machine differ:
  - at most 48 of our threads (cap 48 fixed), SPAdes at 16 threads (2 lanes + one 16-thread lane for binning/QC),
    bowtie2 13 threads (map task = 16), everything nice 10 (the driver renices itself);
  - one Java read preparation (merge/Tadpole/trim) at a time, started only when MemAvailable - outstanding reserves
    - heap >= 40 GB (mem_floor_gb 40, reserve = heap); all samples are prepared before the first assembly starts
    (asm_after_all_prep: a 150g-heap prep and a SPAdes run do not fit together here);
  - SPAdes -m 120 (160 after an out-of-memory failure), at most 2 at once, reserve 60 GB (poplar peaks 29-67 GB);
  - sum of max(reserve, RSS) of our tasks <= 210 GB (scope_mem_gb; the driver's scope has MemoryMax=220G);
  - no GTDB-Tk and no 2026-side QC here (poplar does both); each finished sample is pushed to poplar (push:<s>).
Run:  ~/Documents/py_venv/bin/python -I make_config_local.py
"""
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.normpath(os.path.join(HERE, "..", "..", "..", "manifest.tsv"))
HOME = os.path.expanduser("~")
LR = f"{HOME}/Documents/salterns_rollout_local"
ROOT = f"{LR}/rollout"
ENVS = f"{LR}/envs"
POPLAR_ROOT = "/scratch1/afreiburger/salterns_2026/rollout"
MIGRATED_R2A = ["R2A_A_D1", "R2A_A_D2", "R2A_B_D1", "R2A_B_D2", "R2A_C_D1", "R2A_C_D2", "R2A_B_H2O"]   # = make_config.py
# R2 migration (Andrew, 2026-11-17): the 8 R2 samples <= 9 Gbp whose assembly had not started on poplar (= make_config.py
# MIGRATED_R2). They come after the R2A samples (local order 7-14), so the R2A entries stay as they were; all 12 R2
# libraries are copied here as coverage (copy_inputs_r2.sh, verify_inputs.sh).
MIGRATED_R2 = ["R2_B_D1", "R2_C_D1", "R2_C_D2", "R2_restored_C_black", "R2_A_H2O", "R2_B_H2O", "R2_C_H2O",
               "R2_restored_H2O"]
MIGRATED = MIGRATED_R2A + MIGRATED_R2
GROUP_SIZES = {"R2A": 9, "R2": 12}
TOOLS = {
    "py": f"{HOME}/Documents/py_venv/bin/python",
    "sra": "/nonexistent/sra-tools-not-installed-here (reads copied from poplar)",
    "bbtools": f"{ENVS}/bbtools/bin",
    "spades": f"{ENVS}/spades3153/bin",
    "map": f"{ENVS}/bt_map/bin",
    "metabat2": f"{ENVS}/bt_metabat2/bin",
    "maxbin2": f"{ENVS}/bt_maxbin2/bin",
    "concoct": f"{ENVS}/bt_concoct/bin",
    "dastool": f"{HOME}/Documents/dastool_env/env/bin",
    "gunc": f"{ENVS}/bt_gunc/bin",
    "qcenv": f"{ENVS}/checkm2/bin",
    "checkm2_db": f"{LR}/refdata/checkm2/CheckM2_database/uniref100.KO.1.dmnd",
    "checkm1": f"{ENVS}/bt_checkm1/bin",
    "checkm1_env_ok": f"{ENVS}/.setup/bt_checkm1.ok",
    "checkm1_data": f"{LR}/refdata/checkm1",
    "gunc_db_dir": f"{LR}/refdata/gunc",
}
THREADS = {"dl": 6, "reads": 16, "merge": 16, "tadpole": 16, "qtrimasm": 16, "asm": 16, "contigs": 1, "idx": 16,
           "map": 15, "metabat2": 16, "maxbin2": 16, "concoct": 16, "dastool": 16, "checkm2": 16, "gunc": 16,
           "checkm1": 16, "keep": 1, "import": 2, "push": 1, "report": 1}
# observed on poplar (ledger peaks): map 5, depth 6, binners <= 1, DAS Tool 1, CheckM2 16, GUNC 10, CheckM1 75 GB
RESERVE = {"dl": 4, "reads": 20, "merge": 150, "tadpole": 150, "qtrimasm": 20, "asm": 60, "contigs": 2, "idx": 10,
           "map": 10, "depth": 10, "metabat2": 5, "maxbin2": 5, "concoct": 10, "dastool": 5, "checkm2": 20,
           "gunc": 15, "checkm1": 80, "keep": 2, "import": 4, "push": 2, "report": 4}
BASE = {
    "host": "frankensteinsmonster",
    "cores": 64, "cap_max": 48, "cap_min": 48, "headroom": 16, "mem_floor_gb": 40, "settle_s": 120,
    "poll_s": 20, "sched_s": 120, "retry_wait_s": 300,
    "min_free_scratch_gb": 300, "max_rollout_gb": 1000, "tmp_min_free_gb": 10,
    "kind_limits": {"dl": 1, "reads": 4, "asm": 2, "checkm2": 2, "checkm1": 1, "import": 2, "push": 1},
    "group_limits": {"prep": {"kinds": ["merge", "tadpole", "qtrimasm"], "max": 1},
                     "qc": {"kinds": ["checkm2", "gunc", "checkm1"], "max": 3}},
    "prep_buffer": 9, "asm_lane_threads": 16, "protect_lanes": 2, "background_cap": 16,
    "asm_after_all_prep": True, "scope_mem_gb": 210, "skip_gtdbtk": True, "skip_old_qc": True,
    "attempts": {"default": 2, "dl": 3, "asm": 4, "push": 4},
    "spades": {"klist": "33,55,77,99,127", "mem_gb": [120, 160], "reserve_gb": [60, 160]},
    # 150g as on poplar for <= 9 Gbp (all 7 samples here); a Java OutOfMemoryError retry also gets 150g, because a
    # bigger heap cannot be given on this machine: such a sample goes back to poplar (README "R2A migrated")
    "heap": {"xms": "8g", "le9": "150g", "gt9": "200g", "oom": "150g", "small": "16g", "split_bases": 9000000000},
    "tools": TOOLS, "binning_test_root": "/nonexistent/binning_test (not used on this host)", "old_table": "",
}


def main():
    rows = list(csv.DictReader(open(MANIFEST), delimiter="\t"))
    libs, groups = {}, {}
    for g, n in GROUP_SIZES.items():
        gr = [r for r in rows if r["coverage_group"] == g]
        assert len(gr) == n, (g, len(gr))
        groups[g] = sorted(r["sample"] for r in gr)
        libs.update({r["sample"]: {"srr": r["srr"], "spots": int(r["spots"]), "bases": int(r["bases"]), "group": g}
                     for r in gr})
    samples = {}
    # poplar's order (make_config.py) within each migrated batch; R2A first (orders 0-6 as before), then R2
    poplar_order = json.load(open(os.path.join(HERE, "config.json")))["samples"]
    batches = [sorted(b, key=lambda x: poplar_order[x]["order"]) for b in (MIGRATED_R2A, MIGRATED_R2)]
    for i, s in enumerate(batches[0] + batches[1]):
        g = libs[s]["group"]
        assert libs[s]["bases"] <= BASE["heap"]["split_bases"], s
        assert poplar_order[s]["group"] == groups[g] and poplar_order[s]["import_from"]["staged"] == f"{POPLAR_ROOT}/imports/{s}"
        samples[s] = {"lib": s, "group": groups[g], "group_name": g, "order": i, "bases": libs[s]["bases"], "unit": s,
                      "push_to": {"host": "poplar", "staged": f"{POPLAR_ROOT}/imports/{s}", "marker": "READY"}}
    cfg = dict(BASE, root=ROOT, name="rollout_r2a_local", threads=THREADS, reserve_gb=RESERVE, libraries=libs,
               samples=samples)
    json.dump(cfg, open(os.path.join(HERE, "local_config.json"), "w"), indent=1)

    # ---- local smoke: 1M pairs of R2A_B_H2O (own reads) + R2A_C_H2O (coverage), via the reuse + subsample path on
    #      the verified copies in ROOT/reads; the whole chain incl. QC/keep, then a real push to poplar's smoke_import
    S = f"{LR}/smoke"
    sth = {k: min(v, 8) for k, v in THREADS.items()}
    smoke = dict(BASE, root=S, name="smoke_r2a_local", threads=sth,
                 reserve_gb={k: min(v, 16) for k, v in RESERVE.items()}, sched_s=20, retry_wait_s=20, settle_s=60,
                 scope_mem_gb=56, spades={"klist": "33,55,77,99,127", "mem_gb": [32, 48], "reserve_gb": [8, 16]},
                 heap={"xms": "1g", "le9": "8g", "gt9": "8g", "oom": "12g", "small": "4g", "split_bases": 9000000000},
                 binning_test_root=ROOT)
    # round 1 (300k pairs of the R2A_B_D1 sediment library) assembled to one contig >= 1 kb, too little to bin;
    # round 2 uses 1M pairs of the two water libraries (fewer, more abundant genomes, as poplar's smoke on R2_B_H2O)
    smoke["libraries"] = {
        "SMK_1": {"srr": "-", "spots": 0, "bases": 0, "group": "SMK", "reuse": f"{ROOT}/reads/R2A_B_H2O",
                  "reuse_lib": "R2A_B_H2O", "subsample_pairs": 1000000},
        "SMK_2": {"srr": "-", "spots": 0, "bases": 0, "group": "SMK", "reuse": f"{ROOT}/reads/R2A_C_H2O",
                  "reuse_lib": "R2A_C_H2O", "subsample_pairs": 1000000}}
    smoke["samples"] = {"SMK_L": {"lib": "SMK_1", "group": ["SMK_1", "SMK_2"], "group_name": "SMK", "order": 0,
                                  "bases": libs["R2A_B_H2O"]["bases"], "unit": "R2A_B_H2O",
                                  "push_to": {"host": "poplar", "staged": f"{POPLAR_ROOT}/smoke_import/imports/SMK_L",
                                              "marker": "READY"}}}
    json.dump(smoke, open(os.path.join(HERE, "local_smoke_config.json"), "w"), indent=1)

    # ---- poplar side of the smoke (root rollout/smoke_import): SMK_L as a staged import, 2026-side QC on a copy of
    #      the binning test's 6 SF2_C_D1 bins (as smoke_prep.sh did for SMK_A); no GTDB-Tk (unchanged code, and
    #      poplar's /tmp is near the gate); the report and the end marker wait for the import like gtdbtk does.
    pc = json.load(open(os.path.join(HERE, "smoke_config.json")))
    PI = f"{POPLAR_ROOT}/smoke_import"
    imp = dict({k: v for k, v in pc.items() if k not in ("libraries", "samples")}, root=PI, name="smoke_import",
               host="poplar", skip_gtdbtk=True, tmp_min_free_gb=2, old_table="")
    imp["libraries"] = {k: {"srr": "-", "spots": 0, "bases": 0, "group": "SMK"} for k in ("SMK_1", "SMK_2")}
    imp["samples"] = {"SMK_L": {"lib": "SMK_1", "group": ["SMK_1", "SMK_2"], "group_name": "SMK", "order": 0,
                                "bases": libs["R2A_B_H2O"]["bases"], "unit": "R2A_B_H2O",
                                "import_from": {"staged": f"{PI}/imports/SMK_L", "marker": "READY",
                                                "source": f"frankensteinsmonster:{S}"}}}
    json.dump(imp, open(os.path.join(HERE, "smoke_import_config.json"), "w"), indent=1)
    print(f"wrote local_config.json ({len(samples)} samples: {' '.join(samples)}), local_smoke_config.json, "
          f"smoke_import_config.json")


if __name__ == "__main__":
    main()
