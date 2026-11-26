#!/usr/bin/env python3
"""Write config.json (the real test) and smoke_config.json for driver.py from ../../../manifest.tsv (run locally).

Usage: ~/Documents/py_venv/bin/python -I make_config.py --deadline 2026-11-14T06:18:00Z --smoke-deadline <UTC>
"""
import argparse
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.normpath(os.path.join(HERE, "..", "..", "..", "manifest.tsv"))
SCR = "/scratch1/afreiburger/salterns_2026"
ROOT = f"{SCR}/binning_test"
PILOT = f"{SCR}/tadpole_pilot"
ENVS = f"{SCR}/envs"
TOOLS = {
    "py": f"{ENVS}/bt_map/bin/python",
    "sra": "/scratch1/afreiburger/metag-hybrid/env/miniforge/envs/sra/bin",
    "bbtools": "/scratch1/afreiburger/emilykin/processed/reanalysis_norm/env_bbtools/bin",
    "map": f"{ENVS}/bt_map/bin",
    "metabat2": f"{ENVS}/bt_metabat2/bin",
    "maxbin2": f"{ENVS}/bt_maxbin2/bin",
    "concoct": f"{ENVS}/bt_concoct/bin",
    "dastool": f"{ENVS}/bt_dastool/bin",
    "gunc": f"{ENVS}/bt_gunc/bin",
    "qcenv": "/scratch1/afreiburger/emilykin/processed/.snakemake_envs/8c703f8abdc13cdc1ef374a2bcf37f14_/bin",
    "checkm2_db": "/scratch1/afreiburger/emilykin/processed/refdata/checkm2/CheckM2_database/uniref100.KO.1.dmnd",
    "checkm1": f"{ENVS}/bt_checkm1/bin",
    "checkm1_env_ok": f"{ENVS}/.bt_setup/checkm1.ok",
    "checkm1_data": f"{ROOT}/refdata/checkm1",
    "gunc_db_dir": f"{ROOT}/refdata/gunc",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deadline", required=True)
    ap.add_argument("--smoke-deadline", required=True)
    a = ap.parse_args()
    rows = list(csv.DictReader(open(MANIFEST), delimiter="\t"))
    libs = {}
    for r in rows:
        if r["coverage_group"] not in ("R1", "SF2"):
            continue
        libs[r["sample"]] = {"srr": r["srr"], "spots": int(r["spots"]), "bases": int(r["bases"]), "site": r["coverage_group"]}
    libs["R1_A_D1"]["pilot_qc"] = f"{PILOT}/R1_A_D1"
    libs["SF2_C_D1"]["pilot_qc"] = f"{PILOT}/SF2_C_D1"
    site = lambda g: sorted(k for k, v in libs.items() if v["site"] == g)
    assert len(site("R1")) == 9 and len(site("SF2")) == 9, (site("R1"), site("SF2"))
    cfg = {
        "root": ROOT, "thread_cap": 96, "kind_limits": {"dl": 4, "reads": 6}, "min_free_gb": 400,
        "background_prio": 60, "background_cap": 32,
        "sf2_deadline_utc": a.deadline, "tools": TOOLS, "libraries": libs,
        "samples": {
            "R1_A_D1": {"role": "primary", "site": "R1", "libs": site("R1"),
                        "asm_fasta": f"{PILOT}/R1_A_D1/armC/spades/scaffolds.fasta",
                        "asm_marker": f"{PILOT}/R1_A_D1/.done/C1_spades_default", "asm_pids": []},
            "SF2_C_D1": {"role": "secondary", "site": "SF2", "libs": site("SF2"),
                         "asm_fasta": f"{PILOT}/SF2_C_D1/armC/spades/scaffolds.fasta",
                         "asm_marker": f"{PILOT}/SF2_C_D1/.done/C1_spades_default",
                         "asm_pids": [{"pid": 3634657, "needle": "spades"}, {"pid": 3634245, "needle": "run_sample.sh"},
                                      {"pid": 3634387, "needle": "run_sample.sh"}]},
        },
    }
    json.dump(cfg, open(os.path.join(HERE, "config.json"), "w"), indent=1)
    S = f"{ROOT}/smoke"
    smoke = {
        "root": S, "thread_cap": 48, "kind_limits": {"dl": 1, "reads": 2}, "min_free_gb": 400,
        "background_prio": 60, "background_cap": 32,
        "sf2_deadline_utc": a.smoke_deadline, "tools": TOOLS,
        "libraries": {"SMK1": {"srr": "-", "spots": 0, "bases": 0, "site": "SMK", "pilot_qc": f"{S}/input/pilotlike/SMK1"},
                      "SMK2": {"srr": "-", "spots": 0, "bases": 0, "site": "SMK", "pilot_qc": f"{S}/input/pilotlike/SMK2"}},
        "samples": {
            "SMOKE_A": {"role": "primary", "site": "SMK", "libs": ["SMK1", "SMK2"], "asm_fasta": f"{S}/input/smoke_asm.fa",
                        "asm_marker": f"{S}/input/asm.done", "asm_pids": []},
            "SMOKE_B": {"role": "secondary", "site": "SMK", "libs": ["SMK1", "SMK2"], "asm_fasta": f"{S}/input/never.fa",
                        "asm_marker": f"{S}/input/never.done", "asm_pids": []},
        },
    }
    json.dump(smoke, open(os.path.join(HERE, "smoke_config.json"), "w"), indent=1)
    print("wrote config.json and smoke_config.json")


if __name__ == "__main__":
    main()
