"""Step 2c: does MePn act as an energy source (an energy-generating cycle) in the original vs the rebuilt model?
Max ATPM_c0 with bio1 = 0 on: (A) 'nothing' = H2O, H+, O2 only; (B) MePn + H2O + H+ + O2 + the cell-5 ions, no
carbon source, no Pi; (C) as B with Pi.  Then the pFBA flux path of (B) in the original model.  GLPK.
Output: step2c_energy_check.json"""
import json
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))
import cobra  # noqa: E402

import common as C  # noqa: E402

cobra.Configuration().solver = "glpk"
IONS = {k: v for k, v in C.MEDIA.items() if k not in ("EX_cpd00159_e0", "EX_cpd25960_e0", "EX_cpd00009_e0", "EX_cpd00305_e0")}
TESTS = {"A_water_H_O2_only": {"EX_cpd00001_e0": 1000, "EX_cpd00067_e0": 1000, "EX_cpd00007_e0": 1000},
         "B_MePn_noC_noPi": {**IONS, "EX_cpd25960_e0": 10},
         "C_MePn_noC_withPi": {**IONS, "EX_cpd25960_e0": 10, "EX_cpd00009_e0": 100},
         "D_lactate1_Pi_noMePn": {**IONS, "EX_cpd00159_e0": 1, "EX_cpd00009_e0": 100}}


def atp_tests(model):
    out = {}
    for lab, md in TESTS.items():
        m = model.copy()
        m.medium = {k: v for k, v in md.items() if k in m.reactions}
        m.reactions.bio1.bounds = (0, 0)
        m.objective = "ATPM_c0"
        out[lab] = m.slim_optimize(error_value=float("nan"))
    return out


def path(model, lab="B_MePn_noC_noPi"):
    m = model.copy()
    m.medium = {k: v for k, v in TESTS[lab].items() if k in m.reactions}
    m.reactions.bio1.bounds = (0, 0)
    m.objective = "ATPM_c0"
    if not m.slim_optimize(error_value=0) > 1e-6:
        return None
    s = cobra.flux_analysis.pfba(m)
    rows = []
    for r in m.reactions:
        v = s.fluxes[r.id]
        if abs(v) > 1e-6:
            rows.append((r.id, round(v, 4), r.build_reaction_string(use_metabolite_names=False)[:140],
                         "GPR" if r.gene_reaction_rule else "no GPR"))
    return sorted(rows, key=lambda x: -abs(x[1]))


res = {}
for lab, f in (("original", C.MEPN_MODELS / f"{C.MAG}.contigs.json"), ("rebuilt_v7", "arm_a_v7_pifree_gapfilled.json")):
    m = cobra.io.load_json_model(f)
    res[lab] = {"ATPM_max": atp_tests(m), "pfba_path_B": path(m)}
json.dump(res, open("step2c_energy_check.json", "w"), indent=1)
for lab, d in res.items():
    print(lab, d["ATPM_max"])
    if d["pfba_path_B"]:
        print("  pFBA flux path, test B (", len(d["pfba_path_B"]), "reactions ):")
        for row in d["pfba_path_B"][:45]:
            print("   ", row)
