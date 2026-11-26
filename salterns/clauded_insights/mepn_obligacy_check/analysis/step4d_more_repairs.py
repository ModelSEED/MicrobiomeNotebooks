"""Step 4d: remaining repair variants.
  (1) RAST-supported template ADDITIONS only (no reversals) all at once: does growth without MePn return?
  (2) the RAST-supported single REVERSALS found by step4a (rxn00330 malate synthase for Roseovarius.9/Albimonas.1,
      rxn00139 adenine phosphoribosyltransferase for Roseovarius.11): growth, ATPM(no C), cell 27, FVA;
  (3) the MSGapfill solutions of step4b (Roseovarius.9: rxn00740 + rxn00691; Albimonas.1: rxn00336;
      Roseovarius.11: rxn00770) integrated as a set: same evaluation;
  (4) step4c 'alternatives' that are export transporters: re-evaluated WITH an export-only exchange for the new
      extracellular metabolite (as MSBuilder.add_exchanges_to_model would add), since step4c added none.
Output: step4d_more_repairs.json"""
import json
import sys
from multiprocessing import Pool

sys.dont_write_bytecode = True
sys.path.insert(0, "/tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns/a24fbb2e-9d33-5e05-aeea-a07d777ccde7/scratchpad/mepn_verify/analysis")
from cobra import Reaction  # noqa: E402

import annot as A  # noqa: E402
import vcommon as V  # noqa: E402
from step1b_fva import fva_model  # noqa: E402
from step4a_single_repair_scan import FLIU_REMOVE, gf_bounds, load_template  # noqa: E402
from step4c_repair_and_rerun import atpm_no_carbon, build_repaired  # noqa: E402

INV = {V.iterative_id(m): m for m in V.model_ids()}


def add_export_exchanges(model):
    added = []
    for met in model.metabolites:
        if met.compartment == "e0" and not any(rr.boundary for rr in met.reactions):
            ex = Reaction(f"EX_{met.id}", lower_bound=0, upper_bound=1000)
            ex.add_metabolites({met: -1})
            model.add_reactions([ex])
            added.append(ex.id)
    return added


def cell27_summary(mid, model):
    row = V.abundances_as_run().loc[V.mag_of(mid)]
    out = {}
    for r in V.RATIOS:
        sim = V.cell27(model, row, r)
        out[str(r)] = {"cores": sim, "statuses": {c: sim[c]["status"] for c in sim},
                       "k_methane_per_abundance": {c: sim[c]["methane"] / sim[c]["abundance"] for c in sim},
                       "MePn_over_Pi": {c: sim[c]["MePn"] / sim[c]["phosphate"] if sim[c]["phosphate"] else None for c in sim}}
    return out


def evaluate(model, mid, label, extra=None):
    res = {"label": label, "model": mid, "iterative_id": V.iterative_id(mid)}
    res.update(V.growth_triplet(model))
    res["ATPM_max_no_carbon"] = atpm_no_carbon(model)
    res["cell27"] = cell27_summary(mid, model)
    res["fva"] = fva_model(mid, model=model, tag=label)["fva"]
    if extra:
        res.update(extra)
    return res


def job(args):
    kind, iid, payload = args
    mid = INV[iid]
    if kind == "rast_additions_only":
        model = V.load(mid)
        tpl, tj = load_template()
        add = []
        for trx in tpl.reactions:
            if trx.id in FLIU_REMOVE or trx.id + "0" in model.reactions or not A.support(mid, trx.id).get("rast_supported"):
                continue
            r = trx.to_reaction(model, "0")
            r.bounds = gf_bounds(tj[trx.id])
            add.append(r)
        model.add_reactions(add)
        add_export_exchanges(model)
        g = V.growth_triplet(model)
        return {"kind": kind, "iterative_id": iid, "n_added": len(add), "added": [r.id for r in add], **g}
    if kind == "reversal":
        model = V.load(mid)
        r = model.reactions.get_by_id(payload)
        old = r.bounds
        r.bounds = (-1000, 1000)
        return {"kind": kind, **evaluate(model, mid, f"reverse-open {payload}", {"old_bounds": old})}
    if kind == "set":
        rxns, dirs = payload
        model, added = build_repaired(mid, rxns, dirs)
        ex = add_export_exchanges(model)
        return {"kind": kind, **evaluate(model, mid, "+".join(rxns), {"added": added, "export_exchanges_added": ex})}


if __name__ == "__main__":
    scan = json.load(open(V.OUT / "step4a_single_repair_scan.json"))
    gfs = json.load(open(V.OUT / "step4b_rast_only_and_msgapfill.json"))
    jobs = [("rast_additions_only", iid, None) for iid in scan]
    for iid, r in scan.items():
        for rid in r["single_reversals"]:
            jobs.append(("reversal", iid, rid))
        sol = gfs[iid]["msgapfill"].get("solution") or {}
        new = sol.get("new", {})
        if new:
            jobs.append(("set", iid, (sorted(new), {k: ("forward" if v == ">" else "reverse") for k, v in new.items()})))
        for rid, h in r["single_additions"].items():
            if "transport" in h["name"].lower() or "TRANS-" in h["name"] or "permease" in h["name"].lower() or "uptake" in h["name"].lower() \
                    or "transporter" in h["name"].lower():
                jobs.append(("set", iid, ([rid], {rid: h["flux_direction"]})))
    with Pool(8) as p:
        res = p.map(job, jobs)
    json.dump(res, open(V.OUT / "step4d_more_repairs.json", "w"), indent=1, default=str)
    for x in res:
        if x["kind"] == "rast_additions_only":
            print(f"{x['iterative_id']:15s} RAST-supported additions only ({x['n_added']}): growth a/b/c "
                  f"{x['a_as_is']:.3f}/{x['b_MePn_blocked']:.3f}/{x['c_Pi_blocked_MePn_only_P']:.3f}")
            continue
        k7 = list(x["cell27"]["0.07"]["k_methane_per_abundance"].values())
        k33 = list(x["cell27"]["0.33"]["k_methane_per_abundance"].values())
        print(f"{x['iterative_id']:15s} {x['kind']:8s} {x['label'][:40]:40s} a/b/c {x['a_as_is']:.3f}/{x['b_MePn_blocked']:.3f}/"
              f"{x['c_Pi_blocked_MePn_only_P']:.3f} ATPM(noC) {x['ATPM_max_no_carbon']:.3g} k07 {min(k7):.6f}-{max(k7):.6f} "
              f"k33 {min(k33):.6f}-{max(k33):.6f} st {set(x['cell27']['0.07']['statuses'].values())}")
