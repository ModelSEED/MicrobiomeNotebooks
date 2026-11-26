"""Step 2b: under (b) [cell-5 medium, MePn uptake blocked], which single metabolites restore growth of the
MePn-obligate models when supplied (free source -> X, <= 1000)?  Every cytosolic metabolite is tested.
Also: the MePn-dependent subnetwork = reactions that can carry flux under (a) but not under (b) (own blocked-
reaction finder: max/min each reaction, every LP's support is credited, GLPK, 20 s per LP, timeouts reported).
Output: step2b_rescue_scan.json"""
import json
import sys
import time
from multiprocessing import Pool

sys.dont_write_bytecode = True
sys.path.insert(0, "/tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns/a24fbb2e-9d33-5e05-aeea-a07d777ccde7/scratchpad/mepn_verify/analysis")
from cobra import Reaction  # noqa: E402

import vcommon as V  # noqa: E402

OBLIGATE = ["Roseovarius.9", "Roseovarius.11", "Albimonas.1"]
TOL = 1e-7


def unblocked(model, candidates=None):
    """reactions that can carry |flux| > TOL; returns (set, timeouts)."""
    model.solver.configuration.timeout = 20
    cand = [r.id for r in model.reactions] if candidates is None else list(candidates)
    active, timeouts = set(), []
    for rid in cand:
        if rid in active:
            continue
        r = model.reactions.get_by_id(rid)
        for direction in ("max", "min"):
            if direction == "max" and r.upper_bound <= 0:
                continue
            if direction == "min" and r.lower_bound >= 0:
                continue
            model.objective = model.problem.Objective(r.flux_expression, direction=direction)
            model.solver.optimize()
            st = model.solver.status
            if st != "optimal":
                if st not in ("infeasible",):
                    timeouts.append((rid, direction, st))
                continue
            v = model.solver.objective.value
            if abs(v) > TOL:
                prim = model.solver.primal_values
                for rr in model.reactions:
                    if abs(prim[rr.forward_variable.name] - prim[rr.reverse_variable.name]) > TOL:
                        active.add(rr.id)
                active.add(rid)
                break
    return active, timeouts


def scan(iid):
    t0 = time.time()
    inv = {V.iterative_id(m): m for m in V.model_ids()}
    mid = inv[iid]
    model = V.load(mid)
    V.set_medium_strict(model)
    model.reactions.EX_cpd25960_e0.lower_bound = 0
    model.objective = "bio1"
    rescuers = {}
    src = Reaction("SRCtest", lower_bound=0, upper_bound=0)
    model.add_reactions([src])
    model.solver.configuration.timeout = 20
    for met in model.metabolites:
        if met.compartment != "c0":
            continue
        src.add_metabolites({met: 1})
        src.upper_bound = 1000
        v = model.slim_optimize(error_value=0.0)
        if v > 1e-6:
            rescuers[met.id] = {"name": met.name, "bio1": v}
        src.upper_bound = 0
        src.subtract_metabolites({met: 1})
    model.remove_reactions([src])
    t1 = time.time()
    m_a = V.load(mid)
    V.set_medium_strict(m_a)
    act_a, to_a = unblocked(m_a)
    m_b = V.load(mid)
    V.set_medium_strict(m_b)
    m_b.reactions.EX_cpd25960_e0.lower_bound = 0
    act_b, to_b = unblocked(m_b, candidates=sorted(act_a))
    dep = sorted(act_a - act_b)
    return {"model": mid, "iterative_id": iid, "single_metabolite_rescuers": rescuers,
            "n_unblocked_with_MePn": len(act_a), "n_unblocked_without_MePn": len(act_b),
            "lp_timeouts": {"a": to_a, "b": to_b},
            "MePn_dependent_reactions": {r: {"name": m_a.reactions.get_by_id(r).name,
                                             "reaction": m_a.reactions.get_by_id(r).build_reaction_string(True),
                                             "bounds": m_a.reactions.get_by_id(r).bounds,
                                             "GPR": m_a.reactions.get_by_id(r).gene_reaction_rule}
                                         for r in dep},
            "unblocked_only_without_MePn": sorted(act_b - act_a),
            "time_s": {"rescue": round(t1 - t0, 1), "blocked": round(time.time() - t1, 1)}}


if __name__ == "__main__":
    with Pool(3) as p:
        out = {r["iterative_id"]: r for r in p.map(scan, OBLIGATE)}
    json.dump(out, open(V.OUT / "step2b_rescue_scan.json", "w"), indent=1)
    for iid, r in out.items():
        print(f"== {iid}: {len(r['single_metabolite_rescuers'])} single-metabolite rescuers; unblocked rxns with MePn "
              f"{r['n_unblocked_with_MePn']}, without {r['n_unblocked_without_MePn']}; MePn-dependent reactions "
              f"{len(r['MePn_dependent_reactions'])}; LP timeouts {len(r['lp_timeouts']['a'])}/{len(r['lp_timeouts']['b'])}; {r['time_s']}")
        for k, v in sorted(r["single_metabolite_rescuers"].items(), key=lambda x: -x[1]["bio1"]):
            print(f"   rescuer {k:14s} {v['name'][:50]:50s} bio1 {v['bio1']:.3f}")
