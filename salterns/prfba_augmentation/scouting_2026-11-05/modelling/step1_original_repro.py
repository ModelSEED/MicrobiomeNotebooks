"""Step 1: re-simulate the ORIGINAL Roseovarius.9 MePn model (salterns/models/MePn_consumers/) with the
analysis(1).ipynb cell-27 protocol on GLPK and the as-run (22f4629) abundances; compare with
correlational/saltern_fluxes.json; then FVA of methane / MePn / Pi exchange at the lactate optimum.
Output: step1_original_repro.json (scratch only)."""
import json
import os
import sys
import time

os.chdir(os.path.dirname(os.path.abspath(__file__)))
import cobra  # noqa: E402

import common as C  # noqa: E402

cobra.Configuration().solver = "glpk"
ab = C.abundances()
stored = json.load(open(C.STORED))
mid = C.MAG + ".contigs"
model = cobra.io.load_json_model(C.MEPN_MODELS / f"{mid}.json")
res = {"model_file": str(C.MEPN_MODELS / f"{mid}.json"), "n_reactions": len(model.reactions),
       "solver": "glpk", "abundances_22f4629": ab.to_dict(), "by_multiple": {}}

# growth on the simulation medium with the stored objective (cell 24 'dontGrow' check)
m0 = model.copy()
m0.medium = {k: v for k, v in C.MEDIA.items() if k in m0.reactions}
m0.objective = "bio1"
res["max_bio1_on_cell5_medium"] = m0.slim_optimize()
m0.medium = {k: v for k, v in C.M_GENOME_SCALE.items() if k in m0.reactions}
res["max_bio1_on_fliu_gapfill_medium_noPi"] = m0.slim_optimize()

for mult in (0.07, 0.33):
    t0 = time.perf_counter()
    sim = C.simulate_mepn(model, ab, mult)
    t_sim = time.perf_counter() - t0
    st = {k: stored[str(mult)][k][mid] for k in ("methane", "lactate", "phosphate", "MePn")}
    cmp = {}
    for core in C.CORES:
        s = sim[core]
        cmp[core] = {"status": s["status"],
                     "methane_rerun": s["methane"], "methane_stored": st["methane"][core],
                     "lactate_rerun": s["lactate"], "lactate_stored": st["lactate"][core],
                     "MePn_rerun": s["MePn"], "Pi_rerun": s["phosphate"],
                     "methane_per_abundance_rerun": (s["methane"] / s["abundance"]) if s["methane"] is not None else None,
                     "lactate_per_abundance_rerun": (s["lactate"] / s["abundance"]) if s["lactate"] is not None else None}
    t0 = time.perf_counter()
    fva = C.fva_mepn(model, ab, mult)
    t_fva = time.perf_counter() - t0
    for core in C.CORES:
        a = ab[core]
        cmp[core]["FVA"] = fva[core]
        cmp[core]["FVA_methane_per_abundance"] = [v / a for v in fva[core]["EX_cpd01024_e0"]]
    res["by_multiple"][str(mult)] = {"cores": cmp, "time_sim_s": round(t_sim, 2), "time_fva_s": round(t_fva, 2),
                                     "max_abs_diff_methane": max(abs(cmp[c]["methane_rerun"] - cmp[c]["methane_stored"])
                                                                 for c in C.CORES if cmp[c]["methane_rerun"] is not None)}

json.dump(res, open("step1_original_repro.json", "w"), indent=1)
for mult, d in res["by_multiple"].items():
    print(f"multiple {mult}: sim {d['time_sim_s']} s, FVA {d['time_fva_s']} s, max|dCH4| {d['max_abs_diff_methane']:.3g}")
    for core, c in d["cores"].items():
        print(f"  {core} {c['status']} CH4/a rerun {c['methane_per_abundance_rerun']:.6f} | lact rerun {c['lactate_rerun']:.6g} "
              f"stored {c['lactate_stored']:.6g} | CH4/a FVA [{c['FVA_methane_per_abundance'][0]:.6f}, "
              f"{c['FVA_methane_per_abundance'][1]:.6f}] | MePn FVA {[round(x, 7) for x in c['FVA']['EX_cpd25960_e0']]} "
              f"| Pi FVA {[round(x, 7) for x in c['FVA']['EX_cpd00009_e0']]}")
print("max bio1 on cell-5 medium:", res["max_bio1_on_cell5_medium"], "| on fliu Pi-free medium:",
      res["max_bio1_on_fliu_gapfill_medium_noPi"])
