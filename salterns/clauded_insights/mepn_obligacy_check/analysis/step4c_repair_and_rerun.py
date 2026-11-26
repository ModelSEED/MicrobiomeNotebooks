"""Step 4c: build the repaired models (original model + the repair reaction(s) from the local GramNeg V6 template,
direction locked to the one used, as fliu's _integrate_solution did; GPR = the MAG genes whose RAST function matches a
template role of the reaction, else empty) and re-run the analysis(1).ipynb cell-27 protocol on them (GLPK, as-run
abundances, notebook core order), with sanity checks:
  growth (a)/(b)/(c); max ATPM_c0 on the cell-5 medium without lactate (no carbon source) vs the original model;
  MePn/Pi ratio binding; FVA of CH4 at the lactate optimum.
Every single-reaction repair found by step4a is evaluated; PRIMARY[iid] is the one carried into step 5.
Outputs: step4c_repair_and_rerun.json, repaired_models/<iterativeID>__<repair>.json"""
import json
import sys
from multiprocessing import Pool

sys.dont_write_bytecode = True
sys.path.insert(0, "/tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns/a24fbb2e-9d33-5e05-aeea-a07d777ccde7/scratchpad/mepn_verify/analysis")
import cobra  # noqa: E402

import annot as A  # noqa: E402
import vcommon as V  # noqa: E402
from step1b_fva import fva_model  # noqa: E402
from step4a_single_repair_scan import gf_bounds, load_template  # noqa: E402

PRIMARY = json.load(open(V.OUT / "step4_primary_repairs.json"))   # {iid: [rxn ids]} written after reviewing step4a


def build_repaired(mid, rxn_ids, directions):
    model = V.load(mid)
    tpl, tj = load_template()
    added = []
    for rid in rxn_ids:
        trx = tpl.reactions.get_by_id(rid[:-1])
        r = trx.to_reaction(model, "0")
        lb, ub = gf_bounds(tj[rid[:-1]])
        d = directions.get(rid, "both")
        if d == "forward":
            lb = 0
        elif d == "reverse":
            ub = 0
        r.bounds = (lb, ub)
        sup = A.support(mid, rid[:-1])
        genes = sorted({g for gl in sup.get("rast_supported_roles", {}).values() for g in gl})
        model.add_reactions([r])
        if genes:
            r.gene_reaction_rule = " or ".join(f"{mid}__.RAST.{g}" for g in genes)
        r.notes["repair"] = "added from local GramNeg V6 template to remove MePn obligacy (mepn_verify step 4)"
        added.append({"id": rid, "bounds": r.bounds, "equation": r.build_reaction_string(True),
                      "GPR": r.gene_reaction_rule.replace(mid + "__.RAST.", ""), "RAST_roles": sup.get("rast_supported_roles"),
                      "template_roles": sup.get("template_roles")})
    model.id = mid
    return model, added


def atpm_no_carbon(model):
    med = {k: v for k, v in V.MEDIA.items() if k not in ("EX_cpd00159_e0",)}
    with model:
        model.medium = med
        model.objective = "ATPM_c0"
        return model.slim_optimize(error_value=float("nan"))


def evaluate(args):
    iid, rxn_ids, directions, tag = args
    inv = {V.iterative_id(m): m for m in V.model_ids()}
    mid = inv[iid]
    orig = V.load(mid)
    model, added = build_repaired(mid, rxn_ids, directions)
    res = {"iterative_id": iid, "model": mid, "repair": rxn_ids, "tag": tag, "added": added}
    res.update({k: v for k, v in V.growth_triplet(model).items()})
    res["original_growth"] = V.growth_triplet(orig)
    res["ATPM_max_no_carbon"] = {"original": atpm_no_carbon(orig), "repaired": atpm_no_carbon(model)}
    row = V.abundances_as_run().loc[V.mag_of(mid)]
    res["cell27"] = {}
    for r in V.RATIOS:
        sim = V.cell27(model, row, r)
        res["cell27"][str(r)] = {
            "cores": sim, "statuses": {c: sim[c]["status"] for c in sim},
            "k_methane_per_abundance": {c: sim[c]["methane"] / sim[c]["abundance"] for c in sim},
            "MePn_over_Pi": {c: sim[c]["MePn"] / sim[c]["phosphate"] if sim[c]["phosphate"] else None for c in sim},
            "P_per_biomass": {c: -(sim[c]["MePn"] + sim[c]["phosphate"]) / sim[c]["abundance"] for c in sim},
            "lactate_per_abundance": {c: sim[c]["lactate"] / sim[c]["abundance"] for c in sim}}
    res["fva"] = fva_model(mid, model=model, tag=tag)["fva"]
    out_dir = V.OUT / "repaired_models"
    out_dir.mkdir(exist_ok=True)
    fn = out_dir / f"{iid}__{'+'.join(x[:-3] for x in rxn_ids)}.json"
    cobra.io.save_json_model(model, str(fn))
    res["saved_model"] = str(fn)
    return res


if __name__ == "__main__":
    scan = json.load(open(V.OUT / "step4a_single_repair_scan.json"))
    jobs = []
    for iid, r in scan.items():
        for rid, h in r["single_additions"].items():
            jobs.append((iid, [rid], {rid: h["flux_direction"]}, "primary" if PRIMARY.get(iid) == [rid] else "alternative"))
        for rid in PRIMARY.get(iid, []):
            if rid not in r["single_additions"]:
                raise SystemExit(f"primary repair {rid} for {iid} is not a single-addition hit")
    with Pool(8) as p:
        res = p.map(evaluate, jobs)
    out = {}
    for x in res:
        out.setdefault(x["iterative_id"], {})["+".join(x["repair"])] = x
    json.dump(out, open(V.OUT / "step4c_repair_and_rerun.json", "w"), indent=1, default=str)
    for iid, d in out.items():
        print(f"== {iid}")
        for rep, x in sorted(d.items(), key=lambda kv: kv[1]["tag"] != "primary"):
            k7 = list(x["cell27"]["0.07"]["k_methane_per_abundance"].values())
            k33 = list(x["cell27"]["0.33"]["k_methane_per_abundance"].values())
            rat7 = list(x["cell27"]["0.07"]["MePn_over_Pi"].values())
            fv = [c["CH4_per_abundance_range"] for c in x["fva"]["0.07"].values() if "CH4_per_abundance_range" in c]
            print(f"   {x['tag']:11s} {rep:14s} growth a/b/c {x['a_as_is']:.3f}/{x['b_MePn_blocked']:.3f}/{x['c_Pi_blocked_MePn_only_P']:.3f} "
                  f"ATPM(no C) orig {x['ATPM_max_no_carbon']['original']:.3g} rep {x['ATPM_max_no_carbon']['repaired']:.3g} | "
                  f"k07 {min(k7):.6f}-{max(k7):.6f} k33 {min(k33):.6f}-{max(k33):.6f} MePn/Pi07 {min(rat7):.4f}-{max(rat7):.4f} "
                  f"FVA07 CH4/a [{min(f[0] for f in fv):.6f},{max(f[1] for f in fv):.6f}] st {set(x['cell27']['0.07']['statuses'].values())}")
