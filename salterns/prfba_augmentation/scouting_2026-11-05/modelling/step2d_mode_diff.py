"""Step 2d: what does the original model gain from MePn above the ratio minimum?  R1_A, ratio 0.07, original model:
pFBA at the lactate optimum with the ratio free (MePn = 0.47 a) vs with MePn pinned to 0.07 x Pi (equality).
Lists the largest flux differences.  GLPK.  Output: step2d_mode_diff.json"""
import json, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import cobra
import common as C
cobra.Configuration().solver = "glpk"
ab = C.abundances(); a = float(ab["R1_A"])
orig = cobra.io.load_json_model(C.MEPN_MODELS / f"{C.MAG}.contigs.json")
res = {}
for lab, pin in (("ratio_free", False), ("MePn_pinned_to_0.07Pi", True)):
    m = C.prepare_mepn(orig, 0.07)
    if pin:
        m.constraints["Phosphate_ratio"].ub = 0
    C.set_core(m, a)
    opt = m.slim_optimize()
    m.reactions.EX_cpd00159_e0.lower_bound = opt * (1 + 1e-9) if opt < 0 else opt
    m.objective = "bio1"
    s = cobra.flux_analysis.pfba(m)
    res[lab] = {"lactate_per_biomass": -opt / a, "fluxes": s.fluxes.to_dict()}
f1, f2 = res["ratio_free"]["fluxes"], res["MePn_pinned_to_0.07Pi"]["fluxes"]
diff = sorted(((r, (f1[r] - f2[r]) / a) for r in f1), key=lambda x: -abs(x[1]))
rows = [(r, round(d, 4), orig.reactions.get_by_id(r).build_reaction_string(use_metabolite_names=True)[:150],
         bool(orig.reactions.get_by_id(r).gene_reaction_rule)) for r, d in diff[:40] if abs(d) > 1e-6]
out = {"lactate_per_biomass": {k: v["lactate_per_biomass"] for k, v in res.items()}, "top_flux_differences_per_biomass": rows}
json.dump(out, open("step2d_mode_diff.json", "w"), indent=1)
print(out["lactate_per_biomass"])
for r in rows: print(r)
