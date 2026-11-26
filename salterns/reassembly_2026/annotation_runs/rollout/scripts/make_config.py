#!/usr/bin/env python3
"""Write config.json (the 39-sample rollout) and smoke_config.json for driver.py from ../../../manifest.tsv.
Run locally:  ~/Documents/py_venv/bin/python -I make_config.py
"""
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.normpath(os.path.join(HERE, "..", "..", "..", "manifest.tsv"))
SCR = "/scratch1/afreiburger/salterns_2026"
ROOT = f"{SCR}/rollout"
BT = f"{SCR}/binning_test"
PILOT = f"{SCR}/tadpole_pilot"
ENVS = f"{SCR}/envs"
TOOLS = {
    "py": f"{ENVS}/bt_map/bin/python",
    "sra": "/scratch1/afreiburger/metag-hybrid/env/miniforge/envs/sra/bin",
    "bbtools": "/scratch1/afreiburger/emilykin/processed/reanalysis_norm/env_bbtools/bin",
    "spades": f"{ENVS}/spades3153/bin",
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
    "checkm1_data": f"{BT}/refdata/checkm1",
    "gunc_db_dir": f"{BT}/refdata/gunc",
}
GTDBTK = {"env_bin": f"{ENVS}/gtdbtk280/bin", "env_ok": f"{ENVS}/.ro_setup/gtdbtk280.ok",
          "data_ok": f"{SCR}/refdata/gtdbtk_r232/.ok",
          "fail_markers": [f"{ENVS}/.ro_setup/gtdbtk280.FAIL", f"{SCR}/refdata/gtdbtk_r232/.FAIL"],
          "cpus": 48, "pplacer_cpus": 8, "min_perc_aa": 10}
THREADS = {"dl": 6, "reads": 16, "merge": 16, "tadpole": 16, "qtrimasm": 16, "asm": 24, "contigs": 1, "idx": 16,
           "map": 16, "metabat2": 16, "maxbin2": 16, "concoct": 16, "dastool": 16, "checkm2": 16, "gunc": 16,
           "checkm1": 16, "keep": 1, "import": 2, "gtdbtk": 48, "report": 1}
RESERVE = {"dl": 4, "reads": 20, "merge": 150, "tadpole": 150, "qtrimasm": 20, "asm": 80, "contigs": 2, "idx": 10,
           "map": 10, "depth": 10, "metabat2": 10, "maxbin2": 15, "concoct": 30, "dastool": 20, "checkm2": 40,
           "gunc": 30, "checkm1": 40, "keep": 2, "import": 4, "gtdbtk": 170, "report": 4}
BASE = {
    "cores": 336, "cap_max": 256, "cap_min": 64, "headroom": 16, "mem_floor_gb": 400, "settle_s": 120,
    "poll_s": 20, "sched_s": 120, "retry_wait_s": 300,
    "min_free_scratch_gb": 1000, "max_rollout_gb": 2600, "tmp_min_free_gb": 10,
    "kind_limits": {"dl": 4, "reads": 4, "asm": 6, "checkm2": 3, "checkm1": 3, "gtdbtk": 1, "import": 2},
    "group_limits": {"prep": {"kinds": ["merge", "tadpole", "qtrimasm"], "max": 3},
                     "qc": {"kinds": ["checkm2", "gunc", "checkm1"], "max": 6}},
    "prep_buffer": 2, "asm_lane_threads": 40, "background_cap": 16,
    "attempts": {"default": 2, "dl": 3, "asm": 4},
    "spades": {"klist": "33,55,77,99,127", "mem_gb": [400, 800], "reserve_gb": [80, 400]},
    "heap": {"xms": "8g", "le9": "150g", "gt9": "200g", "oom": "250g", "small": "16g", "split_bases": 9000000000},
    "tools": TOOLS, "binning_test_root": BT,
}
FIRST = ["R2_restored_D_shore", "SF2_A_H2O"]   # tadpole RESULTS.md section 3: the discriminating first batch
# R2A migration (Andrew, 2026-11-16): these 7 R2A samples (<= 9 Gbp, Java heap 150g) run on frankensteinsmonster with
# this same driver (~/Documents/salterns_rollout_local/, make_config_local.py). Here each one is a staged import: its
# import task waits for imports/<s>/READY, which that host's push:<s> task writes after an md5-verified copy; keep:<s>,
# GTDB-Tk, the report and DONE then follow as for every other sample. R2A_A_H2O (9.65 Gbp) and R2A_C_H2O (10.20 Gbp)
# need a 200g heap, which that 251 GB machine cannot give, so they stay here as ordinary samples.
# R2 migration (Andrew, 2026-11-17; poplar's load ~550 had left one assembly lane at the 64-thread floor): the 8 R2
# samples <= 9 Gbp whose assembly had not started here move the same way. R2_B_D1's merge and Tadpole had run here; it is
# redone there from its BBDuk reads (its merge/tadpole markers here are left as they are, no task uses them). R2_A_D1,
# R2_A_D2 (assembling) and R2_B_D2 (12.44 Gbp, 200g heap) stay; all 12 R2 libraries stay here as coverage.
MIGRATED_R2 = ["R2_B_D1", "R2_C_D1", "R2_C_D2", "R2_restored_C_black", "R2_A_H2O", "R2_B_H2O", "R2_C_H2O",
               "R2_restored_H2O"]
MIGRATED = ["R2A_A_D1", "R2A_A_D2", "R2A_B_D1", "R2A_B_D2", "R2A_C_D1", "R2A_C_D2", "R2A_B_H2O"] + MIGRATED_R2
REMOTE = "frankensteinsmonster:~/Documents/salterns_rollout_local/rollout"


def main():
    rows = list(csv.DictReader(open(MANIFEST), delimiter="\t"))
    assert len(rows) == 39
    libs, groups = {}, {}
    for r in rows:
        g = r["coverage_group"]
        L = {"srr": r["srr"], "spots": int(r["spots"]), "bases": int(r["bases"]), "group": g}
        if g in ("R1", "SF2"):      # BBDuk-clean + map reads already made by the binning test
            L["reuse"] = f"{BT}/reads/{r['sample']}"
            L["reuse_lib"] = r["sample"]
        libs[r["sample"]] = L
        groups.setdefault(g, []).append(r["sample"])
    sizes = {g: len(v) for g, v in groups.items()}
    assert sizes == {"R1": 9, "R2": 12, "R2A": 9, "SF2": 9}, sizes
    order = FIRST + [r["sample"] for r in rows if r["coverage_group"] in ("R1", "SF2") and r["sample"] not in FIRST] + \
        [r["sample"] for r in rows if r["coverage_group"] in ("R2", "R2A") and r["sample"] not in FIRST]
    assert sorted(order) == sorted(libs), "order must list every sample once"
    samples = {}
    for i, s in enumerate(order):
        g = libs[s]["group"]
        S = {"lib": s, "group": sorted(groups[g]), "group_name": g, "order": i, "bases": libs[s]["bases"], "unit": s}
        if s in ("R1_A_D1", "SF2_C_D1"):
            S["import_from"] = {"bt_sample": s, "pilot_dir": f"{PILOT}/{s}"}
        if s in MIGRATED:
            assert libs[s]["bases"] <= BASE["heap"]["split_bases"], s
            S["import_from"] = {"staged": f"{ROOT}/imports/{s}", "marker": "READY", "source": REMOTE}
        samples[s] = S
    cfg = dict(BASE, root=ROOT, name="rollout", host="poplar", threads=THREADS, reserve_gb=RESERVE, gtdbtk=GTDBTK,
               old_table=f"{ROOT}/inputs/mags_old_vs_new.tsv", libraries=libs, samples=samples)
    json.dump(cfg, open(os.path.join(HERE, "config.json"), "w"), indent=1)

    # ---- smoke: 300k pairs of R2_B_H2O (downloaded and verified in full first), 300k pairs of R1_B_D1's reused reads,
    #      three samples: SMK_A (normal; SPAdes -m 1 on attempt 1 -> real OOM -> bumped retry), SMK_B (default mode
    #      forced to fail twice -> --meta fallback; no 2026 inputs -> optional old QC fails), SMK_I (import path on the
    #      binning test's SF2_C_D1 outputs); GTDB-Tk 2.8.0 on the real R232 (the mockup cannot place archaea).
    S = f"{ROOT}/smoke"
    r2bh = next(r for r in rows if r["sample"] == "R2_B_H2O")
    sth = {k: (8 if v >= 8 else v) for k, v in THREADS.items()}
    sth.update({"gtdbtk": 8, "reads": 8})
    smoke = dict(BASE, root=S, name="smoke", threads=sth, reserve_gb={k: min(v, 16) for k, v in RESERVE.items()},
                 sched_s=20, retry_wait_s=20, settle_s=60, old_table=f"{ROOT}/inputs/mags_old_vs_new.tsv",
                 spades={"klist": "33,55,77,99,127", "mem_gb": [32, 48], "reserve_gb": [8, 16]},
                 heap={"xms": "1g", "le9": "8g", "gt9": "8g", "oom": "12g", "small": "4g", "split_bases": 9000000000},
                 # GTDB-Tk: the real, verified R232 with the production flags (8 CPUs). The R232 mockup DB was tried
                 # first: it cannot place archaea (empty ar53 refpkg) and needs --full_tree for bacteria.
                 gtdbtk={"env_bin": GTDBTK["env_bin"], "env_ok": GTDBTK["env_ok"], "data_ok": GTDBTK["data_ok"],
                         "fail_markers": GTDBTK["fail_markers"], "cpus": 8, "pplacer_cpus": 2, "min_perc_aa": 10})
    smoke["reserve_gb"]["gtdbtk"] = RESERVE["gtdbtk"]
    smoke["libraries"] = {
        "SMK_DL": {"srr": r2bh["srr"], "spots": int(r2bh["spots"]), "bases": int(r2bh["bases"]), "group": "SMK",
                   "subsample_pairs": 300000},
        "SMK_RU": {"srr": "-", "spots": 0, "bases": 0, "group": "SMK", "reuse": f"{BT}/reads/R1_B_D1",
                   "reuse_lib": "R1_B_D1", "subsample_pairs": 300000}}
    grp = ["SMK_DL", "SMK_RU"]
    smoke["samples"] = {
        "SMK_A": {"lib": "SMK_DL", "group": grp, "group_name": "SMK", "order": 0, "bases": int(r2bh["bases"]),
                  "unit": "R2_B_H2O", "spades_mem_gb": [1, 32]},
        "SMK_B": {"lib": "SMK_DL", "group": grp, "group_name": "SMK", "order": 1, "bases": int(r2bh["bases"]),
                  "unit": "R2_B_H2O", "spades_force_fail_default": True},
        "SMK_I": {"lib": "SMK_RU", "group": [], "group_name": "SMK", "order": 2, "bases": 0, "unit": "SF2_C_D1",
                  "import_from": {"bt_sample": "SF2_C_D1", "pilot_dir": f"{PILOT}/SF2_C_D1"}}}
    json.dump(smoke, open(os.path.join(HERE, "smoke_config.json"), "w"), indent=1)
    print(f"wrote config.json ({len(samples)} samples, {sum(1 for L in libs.values() if not L.get('reuse'))} downloads, "
          f"{sum(1 for L in libs.values() if L.get('reuse'))} reused libraries) and smoke_config.json")
    print("order:", " ".join(order))


if __name__ == "__main__":
    main()
