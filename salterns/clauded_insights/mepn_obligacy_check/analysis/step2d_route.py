"""Step 2d: the MePn-dependent route, by flux comparison.  For each obligate model, the cell-27 problem at bio1 = 1
(r = 0.07, cell-5 medium, lactate uptake minimised, then pFBA at that optimum) is solved for the ORIGINAL model and
for the PRIMARY repaired model (step4c); reactions whose flux differs by > 0.01 (per unit biomass) are listed with
GPR and provenance (draft = in fliu's pre-gapfill draft; gapfilled = cell-93 gapfill solution; module = hand-added
MePn module (gapfill.ipynb cells 47/50/76); drain = SK_ sink added by _add_dm (cell 92)).
Output: step2d_route.json"""
import json
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, "/tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns/a24fbb2e-9d33-5e05-aeea-a07d777ccde7/scratchpad/mepn_verify/analysis")
import cobra  # noqa: E402
from cobra.flux_analysis import pfba  # noqa: E402

import vcommon as V  # noqa: E402

MODULE = {"rxs00001_c0", "rxn26456_c0", "rxn26457_c0", "rxs43120_c0", "rxn26447_c0", "rxn04651_c0"}


def solve(model):
    m = model.copy()
    V.add_ratio_constraint(m, 0.07)
    V.set_medium_strict(m)
    m.reactions.bio1.bounds = (1, 1)
    ex = m.reactions.EX_cpd00159_e0
    ex.bounds = (-1000, 0)
    m.objective = ex
    opt = m.slim_optimize()
    ex.bounds = (opt - 1e-9, min(0, opt + 1e-9))
    return opt, pfba(m).fluxes


def main():
    prov = json.load(open(V.OUT / "step0_provenance.json"))
    rep = json.load(open(V.OUT / "step4c_repair_and_rerun.json"))
    prim = json.load(open(V.OUT / "step4_primary_repairs.json"))
    inv = {V.iterative_id(m): m for m in V.model_ids()}
    out = {}
    for iid, rxns in prim.items():
        mid = inv[iid]
        orig = V.load(mid)
        repm = cobra.io.load_json_model(rep[iid]["+".join(rxns)]["saved_model"])
        draft_ids = {r["id"] for r in json.load(open(V.DRAFTS / f"{mid}__.RAST.json"))["reactions"]}
        lo, fo = solve(orig)
        lr, fr = solve(repm)
        rows = []
        for r in orig.reactions:
            d = fo[r.id] - fr.get(r.id, 0.0)
            if abs(d) > 0.01:
                src = ("module" if r.id in MODULE else "gapfilled" if r.id in prov[mid]["gapfilled"] else
                       "drain" if r.id.startswith("SK_") else "exchange" if r.id.startswith("EX_") else
                       "draft" if r.id in draft_ids else "other_added")
                rows.append({"id": r.id, "name": r.name, "flux_original": fo[r.id], "flux_repaired": fr.get(r.id, 0.0),
                             "difference": d, "provenance": src, "GPR_empty": not r.gene_reaction_rule,
                             "GPR": r.gene_reaction_rule.replace(mid + "__.RAST.", "")[:200],
                             "equation": r.build_reaction_string(True)[:200]})
        for rid in rxns:
            rows.append({"id": rid, "name": repm.reactions.get_by_id(rid).name, "flux_original": 0.0,
                         "flux_repaired": fr[rid], "difference": -fr[rid], "provenance": "repair",
                         "GPR_empty": not repm.reactions.get_by_id(rid).gene_reaction_rule, "GPR": "",
                         "equation": repm.reactions.get_by_id(rid).build_reaction_string(True)[:200]})
        rows.sort(key=lambda x: -abs(x["difference"]))
        out[iid] = {"lactate_per_biomass": {"original": lo, "repaired": lr}, "rows": rows,
                    "n_by_provenance": {p: sum(1 for x in rows if x["provenance"] == p) for p in
                                        ("module", "gapfilled", "drain", "draft", "exchange", "repair", "other_added")}}
    json.dump(out, open(V.OUT / "step2d_route.json", "w"), indent=1, default=float)
    for iid, d in out.items():
        print(f"== {iid}: lactate/biomass original {d['lactate_per_biomass']['original']:.4f} repaired "
              f"{d['lactate_per_biomass']['repaired']:.4f}; {len(d['rows'])} reactions differ by > 0.01; {d['n_by_provenance']}")
        for x in d["rows"][:45]:
            print(f"   {x['id']:15s} {x['flux_original']:+8.4f} -> {x['flux_repaired']:+8.4f} {x['provenance']:9s} "
                  f"GPR={'EMPTY' if x['GPR_empty'] else 'yes  '} {x['equation'][:115]}")


if __name__ == "__main__":
    main()
