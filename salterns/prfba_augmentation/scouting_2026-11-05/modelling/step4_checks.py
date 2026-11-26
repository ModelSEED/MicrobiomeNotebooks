"""Step 4: (1) closed form on the 20 stored MAGs (saltern_fluxes.json, R1_A): MePn/Pi ratio, total P per biomass,
r/(1+r)*(Pi+MePn) vs methane; (2) does arm b1 grow without MePn; (3) which bio1 reactant the original R9 cannot
make with MePn closed (demand test).  GLPK.  Output: step4_checks.json"""
import json, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import cobra, pandas as pd
from cobra import Reaction
import common as C
cobra.Configuration().solver = "glpk"
st = json.load(open(C.STORED)); ab = pd.read_csv(C.AB_AS_RUN, index_col=0)
ids = json.load(open(C.CORR / "iterativeIDs_r232.json")); key = lambda n: n.replace(".contigs", "").replace("bins.", "bins_")
IDK = {key(k): v for k, v in ids.items()}
out = {"closed_form": {}}
for r in ("0.07", "0.33"):
    for mag in st[r]["methane"]:
        c = "R1_A"; a = ab.loc[mag.replace(".contigs", ""), c]
        ch4, mepn, pi = st[r]["methane"][mag][c], -st[r]["MePn"][mag][c], -st[r]["phosphate"][mag][c]
        P = (pi + mepn) / a; rr = float(r)
        out["closed_form"][f"{IDK.get(key(mag), mag)}|{r}"] = {
            "MePn_over_Pi": mepn / pi, "P_per_biomass": P, "CH4_per_biomass": ch4 / a,
            "closed_form_r/(1+r)*P": rr / (1 + rr) * P, "ratio_bound": abs(mepn / pi - rr) < 1e-6}
cf = out["closed_form"]
out["closed_form_summary"] = {
    "n": len(cf), "n_ratio_bound": sum(v["ratio_bound"] for v in cf.values()),
    "max_abs_err_closed_form_if_bound": max(abs(v["CH4_per_biomass"] - v["closed_form_r/(1+r)*P"]) for v in cf.values() if v["ratio_bound"]),
    "P_per_biomass_range": [min(v["P_per_biomass"] for v in cf.values()), max(v["P_per_biomass"] for v in cf.values())],
    "not_bound": {k: round(v["MePn_over_Pi"], 4) for k, v in cf.items() if not v["ratio_bound"]}}
# (2) b1 / b2 with MePn closed
for lab, f in (("b1", "arm_b1_mock_v7_gapfilled.json"), ("b2", "arm_b2_mock_v7_gapfilled.json"), ("arm_a", "arm_a_v7_pifree_gapfilled.json")):
    m = cobra.io.load_json_model(f); m.medium = {k: v for k, v in C.MEDIA.items() if k in m.reactions}; m.objective = "bio1"
    g = m.slim_optimize(); m.reactions.EX_cpd25960_e0.lower_bound = 0
    out[f"{lab}_growth_with_and_without_MePn"] = [g, m.slim_optimize(error_value=float("nan"))]
# (3) demand test in the original R9 with MePn closed
orig = cobra.io.load_json_model(C.MEPN_MODELS / f"{C.MAG}.contigs.json")
orig.medium = {k: v for k, v in C.MEDIA.items() if k in orig.reactions}; orig.reactions.EX_cpd25960_e0.lower_bound = 0
orig.reactions.bio1.bounds = (0, 0)
blocked = []
for met, coef in orig.reactions.bio1.metabolites.items():
    if coef >= 0: continue
    with orig:
        dm = Reaction("DMtest", lower_bound=0, upper_bound=1000); dm.add_metabolites({met: -1}); orig.add_reactions([dm])
        orig.objective = "DMtest"
        v = orig.slim_optimize(error_value=0.0)
        if not v > 1e-9: blocked.append((met.id, met.name))
out["original_R9_bio1_reactants_not_producible_without_MePn"] = blocked
json.dump(out, open("step4_checks.json", "w"), indent=1, default=str)
print(json.dumps({k: v for k, v in out.items() if k != "closed_form"}, indent=1, default=str))
