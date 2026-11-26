"""Step 2e: is MePn the original model's only PRPP source?  (i) PRPP-producing reactions present; (ii) min MePn
uptake per biomass with bio1 fixed (R1_A), lactate free, ratio constraint 0.07; (iii) same with rxn00770 (ribose-
phosphate diphosphokinase, R5P + ATP -> PRPP + AMP) added from the rebuilt model.  GLPK.  Output: step2e_prpp.json"""
import json, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import cobra
import common as C
cobra.Configuration().solver = "glpk"
a = float(C.abundances()["R1_A"])
orig = cobra.io.load_json_model(C.MEPN_MODELS / f"{C.MAG}.contigs.json")
reb = cobra.io.load_json_model("arm_a_v7_pifree_gapfilled.json")
out = {}
for lab, m in (("original", orig), ("rebuilt_v7", reb)):
    prpp = [(r.id, r.name[:50], r.build_reaction_string(True)[:110], bool(r.gene_reaction_rule))
            for r in m.metabolites.cpd00103_c0.reactions if r.id.startswith("rxn")]
    out[lab] = {"PRPP_reactions": prpp, "has_rxn00770": "rxn00770_c0" in m.reactions}
    mm = C.prepare_mepn(m, 0.07)
    mm.reactions.bio1.bounds = (a, a)
    mm.reactions.EX_cpd00159_e0.bounds = (-1000, 0)
    mm.objective = mm.problem.Objective(mm.reactions.EX_cpd25960_e0.flux_expression, direction="max")  # max (least negative) = min uptake
    out[lab]["min_MePn_uptake_per_biomass_lactate_free"] = -mm.slim_optimize(error_value=float("nan")) / a
mm = C.prepare_mepn(orig, 0.07)
r770 = reb.reactions.rxn00770_c0.copy() if "rxn00770_c0" in reb.reactions else None
if r770 is not None and "rxn00770_c0" not in mm.reactions:
    mm.add_reactions([r770]); r770.gene_reaction_rule = ""
    mm.reactions.bio1.bounds = (a, a); mm.reactions.EX_cpd00159_e0.bounds = (-1000, 0)
    mm.objective = mm.problem.Objective(mm.reactions.EX_cpd25960_e0.flux_expression, direction="max")
    out["original_plus_rxn00770"] = {"min_MePn_uptake_per_biomass_lactate_free": -mm.slim_optimize(error_value=float("nan")) / a}
    C.set_core(mm, a); s = mm.optimize()
    out["original_plus_rxn00770"]["cell27_CH4_per_biomass"] = s.fluxes["EX_cpd01024_e0"] / a
    out["original_plus_rxn00770"]["cell27_lactate_per_biomass"] = -s.fluxes["EX_cpd00159_e0"] / a
    out["rxn00770_reaction"] = r770.build_reaction_string(True)
json.dump(out, open("step2e_prpp.json", "w"), indent=1)
print(json.dumps(out, indent=1))
