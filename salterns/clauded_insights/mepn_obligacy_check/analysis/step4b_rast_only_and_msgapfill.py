"""Step 4b: (i) can ANY set of RAST-supported template reactions restore growth without MePn?  All RAST-supported
template reactions absent from the model are added at once (GapfillDirection bounds, fliu's template edits) together
with all RAST-supported reversals of irreversible model reactions; LP feasibility under (b).  If bio1 stays 0, no
repair built only from annotated roles exists, and the minimal repair must contain an unannotated reaction.
(ii) ModelSEEDpy MSGapfill cross-check: local GramNeg V6 template (fliu's edits applied), cell-5 medium without
MePn, target bio1, minimum_obj 0.1, reaction_scores = 1 for RAST-supported template reactions (penalty x1 instead of
x2), binary_check=True, prefilter=False (no ATP test conditions).
Output: step4b_rast_only_and_msgapfill.json"""
import json
import sys
import time
from multiprocessing import Pool

sys.dont_write_bytecode = True
sys.path.insert(0, "/tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns/a24fbb2e-9d33-5e05-aeea-a07d777ccde7/scratchpad/mepn_verify/analysis")
from cobra import Reaction  # noqa: E402

import annot as A  # noqa: E402
import vcommon as V  # noqa: E402
from step4a_single_repair_scan import FLIU_LB0, FLIU_REMOVE, gf_bounds, load_template  # noqa: E402

OBLIGATE = ["Roseovarius.9", "Roseovarius.11", "Albimonas.1"]
MEDIA_NO_MEPN = {k: v for k, v in V.MEDIA.items() if k != "EX_cpd25960_e0"}


def rast_only(iid):
    inv = {V.iterative_id(m): m for m in V.model_ids()}
    mid = inv[iid]
    model = V.load(mid)
    tpl, tj = load_template()
    added, rev = [], []
    for trx in tpl.reactions:
        if trx.id in FLIU_REMOVE or trx.id + "0" in model.reactions:
            continue
        if not A.support(mid, trx.id).get("rast_supported"):
            continue
        r = trx.to_reaction(model, "0")
        r.bounds = gf_bounds(tj[trx.id])
        added.append(r)
    model.add_reactions(added)
    for met in model.metabolites:
        if met.compartment == "e0" and not any(rr.boundary for rr in met.reactions):
            ex = Reaction(f"EX_{met.id}", lower_bound=0, upper_bound=1000)
            ex.add_metabolites({met: -1})
            model.add_reactions([ex])
    for r in model.reactions:
        tid = r.id[:-1]
        if tid in tj and r not in added and A.support(mid, tid).get("rast_supported"):
            lb, ub = gf_bounds(tj[tid])
            if r.lower_bound >= 0 and lb < 0:
                r.lower_bound = -1000
                rev.append(r.id)
            elif r.upper_bound <= 0 and ub > 0:
                r.upper_bound = 1000
                rev.append(r.id)
    # reversals of fliu's deliberate loop-breaking bounds are excluded
    for rid in FLIU_LB0:
        if rid + "0" in model.reactions and rid + "0" in rev:
            model.reactions.get_by_id(rid + "0").lower_bound = 0
            rev.remove(rid + "0")
    V.set_medium_strict(model)
    model.reactions.EX_cpd25960_e0.lower_bound = 0
    model.objective = "bio1"
    v = model.slim_optimize(error_value=0.0)
    return {"n_rast_supported_additions": len(added), "n_rast_supported_reversals": len(rev), "bio1_without_MePn": v}


def msgapfill(iid):
    from modelseedpy import MSGapfill, MSMedia
    inv = {V.iterative_id(m): m for m in V.model_ids()}
    mid = inv[iid]
    model = V.load(mid)
    tpl, tj = load_template()
    for rid in FLIU_REMOVE:
        if rid in tpl.reactions:
            tpl.reactions.get_by_id(rid).lower_bound = 0
            tpl.reactions.get_by_id(rid).upper_bound = 0
    for rid in FLIU_LB0:
        if rid in tpl.reactions:
            tpl.reactions.get_by_id(rid).lower_bound = 0
    scores = {}
    for trx in tpl.reactions:
        sup = A.support(mid, trx.id)
        if sup.get("rast_supported"):
            genes = sorted({g for gl in sup["rast_supported_roles"].values() for g in gl})
            scores[trx.id[:-2] if trx.id.endswith("_c") else trx.id.split("_")[0]] = {g: 1.0 for g in genes}
    model.objective = "bio1"
    media = MSMedia.from_dict({k[3:-3]: v for k, v in MEDIA_NO_MEPN.items()})
    t0 = time.time()
    gf = MSGapfill(model, default_gapfill_templates=[tpl], test_conditions=[], default_target="bio1",
                   reaction_scores=scores)
    sol = gf.run_gapfilling(media=media, target="bio1", minimum_obj=0.1, binary_check=True, prefilter=False)
    res = {"time_s": round(time.time() - t0, 1), "n_scored_reactions": len(scores)}
    if not sol:
        res["solution"] = None
        return res
    res["solution"] = {"new": dict(sol.get("new", {})), "reversed": dict(sol.get("reversed", {}))}
    res["new_RAST_supported"] = {k: A.support(mid, k[:-1]).get("rast_supported") for k in sol.get("new", {})}
    res["reversed_RAST_supported"] = {k: A.support(mid, k[:-1]).get("rast_supported") for k in sol.get("reversed", {})}
    return res


def run(iid):
    out = {"iterative_id": iid, "rast_only": rast_only(iid)}
    try:
        out["msgapfill"] = msgapfill(iid)
    except Exception as e:  # report, do not hide
        import traceback
        out["msgapfill"] = {"error": f"{type(e).__name__}: {e}", "trace": traceback.format_exc()[-2000:]}
    return out


if __name__ == "__main__":
    with Pool(3) as p:
        res = {r["iterative_id"]: r for r in p.map(run, OBLIGATE)}
    json.dump(res, open(V.OUT / "step4b_rast_only_and_msgapfill.json", "w"), indent=1, default=str)
    for iid, r in res.items():
        print(f"== {iid}: RAST-only universe: {r['rast_only']}")
        print(f"   MSGapfill: {json.dumps(r['msgapfill'], default=str)[:1500]}")
