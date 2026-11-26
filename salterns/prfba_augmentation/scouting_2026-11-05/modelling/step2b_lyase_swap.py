"""Step 2b: is the original k = 0.469737 a property of fliu's hand-written C-P lyase (rxs43120: NADH -> NAD+,
CH4 straight to e0) rather than of the genome?  On the rebuilt arm-(a) v7 model: (i) as built (template
rxn26459, H+, CH4 in c0); (ii) rxn26459 closed and fliu's rxs43120 added with the same GPR; (iii) the
reverse on the ORIGINAL model (rxs43120 closed, template rxn26459 stoichiometry added).  Also the lactate
cost per unit biomass.  GLPK, cell-27 protocol, as-run abundances.  Output: step2b_lyase_swap.json"""
import json
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))
import cobra  # noqa: E402
from cobra import Reaction  # noqa: E402

import common as C  # noqa: E402

cobra.Configuration().solver = "glpk"
ab = C.abundances()


def k(model, mult):
    sim = C.simulate_mepn(model, ab, mult)
    a = sim["R1_A"]
    if a["status"] != "optimal":
        return {"status": {c: sim[c]["status"] for c in C.CORES}}
    return {"status": a["status"], "CH4_per_biomass": a["methane"] / a["abundance"],
            "MePn_per_biomass": -a["MePn"] / a["abundance"], "Pi_per_biomass": -a["phosphate"] / a["abundance"],
            "lactate_per_biomass": -a["lactate"] / a["abundance"],
            "all_cores_same_k": max(abs(sim[c]["methane"] / sim[c]["abundance"] - a["methane"] / a["abundance"]) for c in C.CORES) < 1e-9}


def add_fliu_lyase(m, gpr=""):
    mt = m.metabolites
    r = Reaction("rxs43120_c0", "rxs43120_c0", "MePn Degradation", 0, 1000)
    r.add_metabolites({mt.cpd00017_c0: -1, mt.cpd26519_c0: -1, mt.cpd01024_e0: 1, mt.cpd03091_c0: 1,
                       mt.cpd00060_c0: 1, mt.cpd26508_c0: 1, mt.cpd00003_c0: 1, mt.cpd00004_c0: -1})
    m.add_reactions([r])
    r.gene_reaction_rule = gpr
    return r


def add_template_lyase(m):
    mt = m.metabolites
    r = Reaction("rxn26459_c0", "RXN0-6734.c (template stoichiometry)", "", 0, 1000)
    if "cpd01024_c0" not in mt:
        from cobra import Metabolite
        m.add_metabolites([Metabolite("cpd01024_c0", "CH4", "Methane [c0]", 0, "c0")])
    r.add_metabolites({mt.cpd00017_c0: -1, mt.cpd00067_c0: -1, mt.cpd26519_c0: -1, mt.cpd00060_c0: 1,
                       mt.get_by_id("cpd01024_c0"): 1, mt.cpd03091_c0: 1, mt.cpd26508_c0: 1})
    tr = Reaction("rxn10471_c0", "Methane transport", "", -1000, 1000)
    tr.add_metabolites({mt.cpd01024_e0: -1, mt.get_by_id("cpd01024_c0"): 1})
    m.add_reactions([r] + ([tr] if "rxn10471_c0" not in m.reactions else []))
    return r


out = {}
rebuilt = cobra.io.load_json_model("arm_a_v7_pifree_gapfilled.json")
orig = cobra.io.load_json_model(C.MEPN_MODELS / f"{C.MAG}.contigs.json")
out["rxs43120_mass_balance_in_original"] = orig.reactions.rxs43120_c0.check_mass_balance()
for mult in (0.07, 0.33):
    d = out.setdefault(str(mult), {})
    d["rebuilt_as_built_template_lyase"] = k(rebuilt, mult)
    m = rebuilt.copy()
    gpr = m.reactions.rxn26459_c0.gene_reaction_rule
    m.reactions.rxn26459_c0.bounds = (0, 0)
    add_fliu_lyase(m, gpr)
    d["rebuilt_with_fliu_rxs43120"] = k(m, mult)
    d["original_as_stored"] = k(orig, mult)
    m = orig.copy()
    m.reactions.rxs43120_c0.bounds = (0, 0)
    add_template_lyase(m)
    d["original_with_template_rxn26459"] = k(m, mult)
json.dump(out, open("step2b_lyase_swap.json", "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
