"""Step 4a: exhaustive size-1 repair scan for the MePn-obligate models.  Universal set = local ModelSEED GramNeg V6
template (the family fliu gapfilled with: gapfill.ipynb cells 62/72, 'GramNegModelTemplateV6_test', KBase ws 12218,
not available locally), with fliu's cell-72 template edits applied (4 reactions removed, lower bound 0 on 20).
Condition (b): cell-5 medium, EX_cpd25960_e0 lower bound 0, maximise bio1.
  (i)  every template reaction absent from the model is added alone, with the bounds its GapfillDirection allows
       ('=' both, '>' forward, '<' reverse); new extracellular metabolites get an export-only exchange (as
       MSBuilder.add_exchanges_to_model would add, with uptake closed by the medium);
  (ii) every model reaction that is irreversible in the model is opened in the other direction alone, when the
       template's GapfillDirection allows that direction.
A candidate 'repairs' the model if max bio1 > 1e-3.  RAST support from annot.support().
Output: step4a_single_repair_scan.json"""
import json
import sys
import time
from multiprocessing import Pool

sys.dont_write_bytecode = True
sys.path.insert(0, "/tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns/a24fbb2e-9d33-5e05-aeea-a07d777ccde7/scratchpad/mepn_verify/analysis")
from cobra import Reaction  # noqa: E402

import annot as A  # noqa: E402
import vcommon as V  # noqa: E402

OBLIGATE = ["Roseovarius.9", "Roseovarius.11", "Albimonas.1"]
FLIU_REMOVE = ['rxn23850_c', 'rxn09193_c', 'rxn05315_c', 'rxn05206_c']                 # gapfill.ipynb cell 72
FLIU_LB0 = ['rxn09188_c', 'rxn00929_c', 'rxn00145_c', 'rxn00499_c', 'rxn00146_c', 'rxn08934_c', 'rxn08941_c',
            'rxn08935_c', 'rxn08942_c', 'rxn11942_c', 'rxn11941_c', 'rxn43329_c', 'rxn11940_c', 'rxn40559_c',
            'rxn04159_c', 'rxn04705_c', 'rxn14415_c', 'sul00003_c', 'rxn48579_c', 'rxn39860_c']
THRESH = 1e-3


def load_template():
    from modelseedpy.core.mstemplate import MSTemplateBuilder
    tj = A.template_json()
    tpl = MSTemplateBuilder.from_dict(tj).build()
    return tpl, {r["id"]: r for r in tj["reactions"]}


def gf_bounds(tr_json):
    d = tr_json.get("GapfillDirection", tr_json["direction"])
    lb = -1000 if d in ("=", "<") else 0
    ub = 1000 if d in ("=", ">") else 0
    if tr_json["id"] in FLIU_LB0:
        lb = 0
    return lb, ub


def scan(iid):
    t0 = time.time()
    inv = {V.iterative_id(m): m for m in V.model_ids()}
    mid = inv[iid]
    model = V.load(mid)
    tpl, tj = load_template()
    V.set_medium_strict(model)
    model.reactions.EX_cpd25960_e0.lower_bound = 0
    model.objective = "bio1"
    base = model.slim_optimize(error_value=0.0)
    cands, new_rxns, gfb = [], [], {}
    for trx in tpl.reactions:
        if trx.id in FLIU_REMOVE:
            continue
        rid = trx.id + "0"
        if rid in model.reactions:
            continue
        r = trx.to_reaction(model, "0")
        gfb[rid] = gf_bounds(tj[trx.id])
        if gfb[rid] == (0, 0):
            continue
        r.bounds = (0, 0)
        new_rxns.append(r)
        cands.append(rid)
    model.add_reactions(new_rxns)
    # export-only exchanges for extracellular metabolites that have none
    ex_added = []
    for met in model.metabolites:
        if met.compartment == "e0" and not any(rr.boundary for rr in met.reactions):
            ex = Reaction(f"EX_{met.id}", lower_bound=0, upper_bound=1000)
            ex.add_metabolites({met: -1})
            ex_added.append(ex)
    model.add_reactions(ex_added)
    model.solver.configuration.timeout = 30
    check = model.slim_optimize(error_value=0.0)
    hits = {}
    for rid in cands:
        r = model.reactions.get_by_id(rid)
        r.bounds = gfb[rid]
        v = model.slim_optimize(error_value=0.0)
        if v > THRESH:
            sol_flux = model.solver.primal_values[r.forward_variable.name] - model.solver.primal_values[r.reverse_variable.name]
            hits[rid] = {"bio1": v, "flux_direction": "forward" if sol_flux > 0 else ("reverse" if sol_flux < 0 else "zero"),
                         "equation": r.build_reaction_string(True), "name": r.name,
                         "template_gapfill_bounds": gfb[rid], "RAST": A.support(mid, rid[:-1])}
        r.bounds = (0, 0)
    # (ii) reversal of existing irreversible reactions
    rev_hits = {}
    for r in list(model.reactions):
        if r.id in cands or r.id.startswith(("EX_", "SK_", "DM_", "bio")) or r.id == "ATPM_c0":
            continue
        tid = r.id[:-1]
        if tid not in tj:
            continue
        lb, ub = gf_bounds(tj[tid])
        old = r.bounds
        if r.lower_bound >= 0 and lb < 0:
            r.lower_bound = -1000
        elif r.upper_bound <= 0 and ub > 0:
            r.upper_bound = 1000
        else:
            continue
        v = model.slim_optimize(error_value=0.0)
        if v > THRESH:
            rev_hits[r.id] = {"bio1": v, "model_bounds": old, "opened_to": r.bounds, "equation": r.build_reaction_string(True),
                              "GPR": r.gene_reaction_rule.replace(mid + "__.RAST.", ""), "RAST": A.support(mid, tid)}
        r.bounds = old
    return {"model": mid, "iterative_id": iid, "bio1_without_MePn_before": base, "bio1_check_after_adding_closed": check,
            "n_candidates": len(cands), "n_export_exchanges_added": len(ex_added), "single_additions": hits,
            "single_reversals": rev_hits, "time_s": round(time.time() - t0, 1)}


if __name__ == "__main__":
    with Pool(3) as p:
        out = {r["iterative_id"]: r for r in p.map(scan, OBLIGATE)}
    json.dump(out, open(V.OUT / "step4a_single_repair_scan.json", "w"), indent=1, default=str)
    for iid, r in out.items():
        print(f"== {iid}: before {r['bio1_without_MePn_before']:.3g} (check {r['bio1_check_after_adding_closed']:.3g}); "
              f"{r['n_candidates']} candidates; {len(r['single_additions'])} single additions and "
              f"{len(r['single_reversals'])} single reversals restore growth; {r['time_s']} s")
        for rid, h in sorted(r["single_additions"].items(), key=lambda x: -x[1]["bio1"]):
            print(f"   ADD {rid:12s} bio1 {h['bio1']:.3f} {h['flux_direction']:7s} RAST={h['RAST'].get('rast_supported')} "
                  f"{h['name'][:60]} | {h['equation'][:110]}")
        for rid, h in sorted(r["single_reversals"].items(), key=lambda x: -x[1]["bio1"]):
            print(f"   REV {rid:12s} bio1 {h['bio1']:.3f} {h['model_bounds']}->{h['opened_to']} RAST={h['RAST'].get('rast_supported')} "
                  f"GPR={'yes' if h['GPR'] else 'EMPTY'} | {h['equation'][:110]}")
