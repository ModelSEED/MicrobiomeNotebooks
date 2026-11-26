"""Step 3b: per-reaction probabilities in simulation.  Imports modelseedpy from the ModelSEEDpy-prfba worktree
(branch fix-prfba-single-model, 09856f5; NOT the venv's editable ~/Documents/ModelSEEDpy).  R1_A core, ratio 0.07.
 T1 MSProbability.prFBA on the saltern objective (max EX_cpd00159_e0, bio1 fixed)      -> expected ValueError
 T2 same with min_growth=None                                                           -> 95 % floor of a negative optimum
 T3 prFBA(biomass_id='bio1') with bio1 fixed                                            -> floor redundant; min sum(1-p)|v| + 100|EX|
 T4 lexicographic: min-lactate LP, lactate held within tol of its optimum, then the prFBA coefficients
    (msprobability.py lines 268-277) minimised; ex_weight 0 or 100; methane FVA at the same lactate bound
 T5 MSProbability.apply_threshold(model, 0.5) and a function-level threshold (drop the p = 0.1 mock reaction)
Models: original (MAG genes p = 1, gene-less 0), arm (a) v7, arm (b1) mock.  GLPK.  Output: step3b_prfba.json"""
import json
import os
import sys

sys.path.insert(0, "/home/freiburger/Documents/ModelSEEDpy-prfba")
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import modelseedpy  # noqa: E402
assert modelseedpy.__file__.startswith("/home/freiburger/Documents/ModelSEEDpy-prfba/"), modelseedpy.__file__
import mscommunity  # noqa: E402,F401  (msprobability.py line 3 imports it at module level)
import cobra  # noqa: E402
from cobra.flux_analysis import flux_variability_analysis  # noqa: E402
from optlang.symbolics import Zero  # noqa: E402
from modelseedpy.core.msprobability import MSProbability  # noqa: E402
from modelseedpy.core.msmodelutl import MSModelUtil  # noqa: E402

import common as C  # noqa: E402

cobra.Configuration().solver = "glpk"
a = float(C.abundances()["R1_A"])
MULT = 0.07


def with_p(model):
    m = model.copy()
    for r in m.reactions:
        if "probability" not in r.notes:
            r.notes["probability"] = "1.0" if r.genes else "0.0"
    return m


def prepared(model):
    m = C.prepare_mepn(model, MULT)
    C.set_core(m, a)
    return m


def summary(m, sol):
    if sol is None or sol.status != "optimal":
        return {"status": None if sol is None else sol.status}
    f = sol.fluxes
    return {"status": sol.status, "CH4_per_biomass": f["EX_cpd01024_e0"] / a, "lactate_per_biomass": -f["EX_cpd00159_e0"] / a,
            "MePn_per_biomass": -f["EX_cpd25960_e0"] / a, "Pi_per_biomass": -f["EX_cpd00009_e0"] / a,
            "bio1": f["bio1"], "objective": sol.objective_value}


def lexicographic(model, tol, ex_weight, min_prob=0.01, prob_exp=1, missing=0.0):
    m = prepared(model)
    v = m.slim_optimize()
    m.reactions.EX_cpd00159_e0.lower_bound = v * (1 + tol) - 1e-9 * abs(v)
    # direct min/max of CH4 export under the lactate bound (cobra FVA's floor lb = fraction * optimum is wrong
    # for a negative maximised objective: fraction 0 forces EX_cpd00159_e0 >= 0)
    rng = []
    for sense in ("min", "max"):
        with m:
            m.objective = m.problem.Objective(m.reactions.EX_cpd01024_e0.flux_expression, direction=sense)
            rng.append(m.slim_optimize(error_value=float("nan")) / a)
    coef = {}
    for r in m.reactions:                                  # msprobability.py prFBA coefficient loop
        if r.id[:3] == "rxn":
            w = max(min_prob, 1 - float(r.notes.get("probability", missing)) ** prob_exp)
            coef[r.forward_variable] = w
            coef[r.reverse_variable] = w
        elif r.id[:3] == "EX_" and ex_weight:
            coef[r.forward_variable] = ex_weight
            coef[r.reverse_variable] = ex_weight
    MSModelUtil.get(m).add_objective(Zero, "min", coef)
    out = summary(m, m.optimize())
    out["stage1_lactate_per_biomass"] = -v / a
    out["CH4_range_per_biomass_at_this_lactate_bound"] = rng
    return out


models = {"original": with_p(cobra.io.load_json_model(C.MEPN_MODELS / f"{C.MAG}.contigs.json")),
          "arm_a_v7": with_p(cobra.io.load_json_model("arm_a_v7_pifree_gapfilled.json")),
          "arm_b1_mock": with_p(cobra.io.load_json_model("arm_b1_mock_v7_gapfilled.json"))}
res = {"modelseedpy_file": modelseedpy.__file__, "abundance_R1_A": a, "ratio": MULT}
for name, model in models.items():
    d = res[name] = {}
    d["plain_FBA_cell27"] = summary(*(lambda m: (m, m.optimize()))(prepared(model)))
    try:                                                             # T1
        MSProbability.prFBA(prepared(model))
        d["T1_prFBA_lactate_objective"] = "no error"
    except Exception as e:
        d["T1_prFBA_lactate_objective"] = f"{type(e).__name__}: {str(e)[:220]}"
    try:                                                             # T2
        m = prepared(model)
        d["T2_prFBA_min_growth_None"] = summary(m, MSProbability.prFBA(m, min_growth=None))
        d["T2_floor_lb"] = m.constraints["min_value"].lb
    except Exception as e:
        d["T2_prFBA_min_growth_None"] = f"{type(e).__name__}: {str(e)[:220]}"
    try:                                                             # T3
        m = prepared(model)
        d["T3_prFBA_biomass_id_bio1"] = summary(m, MSProbability.prFBA(m, biomass_id="bio1"))
    except Exception as e:
        d["T3_prFBA_biomass_id_bio1"] = f"{type(e).__name__}: {str(e)[:220]}"
    for tol in (0.0, 0.05):                                          # T4
        for exw in (0.0, 100.0):
            d[f"T4_lexi_tol{tol}_exw{exw:g}"] = lexicographic(model, tol, exw)
    m = model.copy()                                                 # T5
    try:
        MSProbability.apply_threshold(m, 0.5)
        mm = m.copy()
        mm.medium = {k: v for k, v in C.MEDIA.items() if k in mm.reactions}
        mm.objective = "bio1"
        d["T5_apply_threshold_0.5"] = {"n_reactions_closed": sum(1 for r in m.reactions if r.bounds == (0, 0)),
                                       "of": len(m.reactions), "max_bio1_cell5": mm.slim_optimize(error_value=float("nan"))}
    except Exception as e:
        d["T5_apply_threshold_0.5"] = f"{type(e).__name__}: {str(e)[:200]}"
    if name == "arm_b1_mock":
        for rid in ("rxn00844_c0", "rxn05602_c0", "rxn26447_c0"):
            m = model.copy()
            if rid in m.reactions:
                m.reactions.get_by_id(rid).bounds = (0, 0)
                d[f"T5_plainFBA_without_{rid}_p{model.reactions.get_by_id(rid).notes['probability']}"] = \
                    summary(*(lambda x: (x, x.optimize()))(prepared(m)))
json.dump(res, open("step3b_prfba.json", "w"), indent=1, default=str)
for name in models:
    print("=====", name)
    for k, v in res[name].items():
        print(" ", k, ":", {kk: (round(vv, 6) if isinstance(vv, float) else vv) for kk, vv in v.items()} if isinstance(v, dict) else v)
