"""Step 1 (all 24 models in salterns/models/MePn_consumers/):
  growth on the cell-5 medium (a) as is, (b) MePn uptake blocked, (c) Pi blocked (MePn the only P source);
  obligate := (b) < 1e-6 while (a) > 0;
  cell-27 re-run at r = 0.07 and 0.33 (GLPK, as-run abundances, notebook core order) vs saltern_fluxes.json;
  MePn/Pi ratio and whether the ratio constraint binds (FVA at the lactate optimum: step1b_fva.py).
Output: step1_growth_and_repro.json"""
import json
import sys
from multiprocessing import Pool

sys.dont_write_bytecode = True
sys.path.insert(0, "/tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns/a24fbb2e-9d33-5e05-aeea-a07d777ccde7/scratchpad/mepn_verify/analysis")
import vcommon as V  # noqa: E402


def fva_at_optimum(model_in, row, multiple):
    from cobra.flux_analysis import flux_variability_analysis
    model = model_in.copy()
    V.add_ratio_constraint(model, multiple)
    V.set_medium_strict(model)
    out = {}
    for sample, a in row.items():
        model.reactions.bio1.bounds = (a, a)
        model.reactions.EX_cpd00159_e0.bounds = (-1000, 0)
        model.objective = "EX_cpd00159_e0"
        try:
            r = flux_variability_analysis(model, [V.CH4_EX, V.MEPN_EX, V.PI_EX], fraction_of_optimum=1.0, processes=1)
            out[sample] = {k: [float(r.loc[k, "minimum"]), float(r.loc[k, "maximum"])] for k in r.index}
        except Exception as e:  # infeasible
            out[sample] = {"error": str(e)[:120]}
        model.reactions.bio1.bounds = (0, 0)
    return out


def run(mid):
    model = V.load(mid)
    ab = V.abundances_as_run()
    st = V.stored()
    mag = V.mag_of(mid)
    res = {"model": mid, "iterative_id": V.iterative_id(mid), "n_reactions": len(model.reactions),
           "medium_keys_missing": [k for k in V.MEDIA if k not in model.reactions]}
    res.update(V.growth_triplet(model))
    res["in_stored"] = mid in st["0.07"]["methane"]
    res["cell27"] = {}
    if not res["medium_keys_missing"] and mag in ab.index:
        row = ab.loc[mag]
        for r in V.RATIOS:
            sim = V.cell27(model, row, r)
            d = {"cores": sim}
            if res["in_stored"]:
                diffs = {c: abs(sim[c]["methane"] - st[str(r)]["methane"][mid][c]) for c in sim}
                d["max_abs_diff_methane_vs_stored"] = max(diffs.values())
                d["max_abs_diff_lactate_vs_stored"] = max(abs(sim[c]["lactate"] - st[str(r)]["lactate"][mid][c]) for c in sim)
                d["max_abs_diff_MePn_vs_stored"] = max(abs(sim[c]["MePn"] - st[str(r)]["MePn"][mid][c]) for c in sim)
                d["max_abs_diff_Pi_vs_stored"] = max(abs(sim[c]["phosphate"] - st[str(r)]["phosphate"][mid][c]) for c in sim)
                d["max_rel_diff_methane_vs_stored"] = max(diffs[c] / abs(st[str(r)]["methane"][mid][c]) for c in sim
                                                          if st[str(r)]["methane"][mid][c] != 0)
                d["stored_k_methane_per_abundance"] = {c: st[str(r)]["methane"][mid][c] / sim[c]["abundance"] for c in sim}
            ok = [c for c in sim if sim[c]["status"] == "optimal"]
            d["statuses"] = {c: sim[c]["status"] for c in sim}
            d["k_methane_per_abundance"] = {c: sim[c]["methane"] / sim[c]["abundance"] for c in sim}
            d["MePn_over_Pi"] = {c: (sim[c]["MePn"] / sim[c]["phosphate"]) if sim[c]["phosphate"] else None for c in sim}
            d["ratio_binds"] = {c: (abs(sim[c]["MePn"] / sim[c]["phosphate"] - r) < 1e-6) if sim[c]["phosphate"] else None
                                for c in sim}
            d["P_per_biomass"] = {c: -(sim[c]["MePn"] + sim[c]["phosphate"]) / sim[c]["abundance"] for c in sim}
            d["closed_form_r_over_1pr_x_P"] = {c: r / (1 + r) * d["P_per_biomass"][c] for c in sim}
            d["n_optimal"] = len(ok)
            res["cell27"][str(r)] = d
    return res


if __name__ == "__main__":
    mids = V.model_ids()
    with Pool(8) as p:
        results = p.map(run, mids)
    out = {r["model"]: r for r in results}
    json.dump(out, open(V.OUT / "step1_growth_and_repro.json", "w"), indent=1)
    print(f"{'iterativeID':22s} {'a':>8s} {'b_noMePn':>9s} {'c_MePnOnlyP':>11s} oblig stored  "
          f"maxdCH4(0.07) maxdCH4(0.33)  k(0.07)     k(0.33)   binds07 binds33  MePn/Pi07")
    for mid in mids:
        r = out[mid]
        line = f"{r['iterative_id']:22s} {r['a_as_is']:8.4f} {r['b_MePn_blocked']:9.4f} {r['c_Pi_blocked_MePn_only_P']:11.4f} " \
               f"{str(r['obligate']):5s} {str(r['in_stored']):6s}"
        if r["in_stored"]:
            c7, c33 = r["cell27"]["0.07"], r["cell27"]["0.33"]
            k7 = list(c7["k_methane_per_abundance"].values())
            k33 = list(c33["k_methane_per_abundance"].values())
            line += f" {c7['max_abs_diff_methane_vs_stored']:.2e}     {c33['max_abs_diff_methane_vs_stored']:.2e}    " \
                    f"{min(k7):.6f}-{max(k7):.6f} {min(k33):.6f}-{max(k33):.6f} " \
                    f"{sum(bool(x) for x in c7['ratio_binds'].values())}/6 {sum(bool(x) for x in c33['ratio_binds'].values())}/6 " \
                    f"{min(v for v in c7['MePn_over_Pi'].values() if v is not None):.3f}-{max(v for v in c7['MePn_over_Pi'].values() if v is not None):.3f}" \
                    f" st07={set(c7['statuses'].values())} st33={set(c33['statuses'].values())}"
        print(line)
    allmax = {r: max(out[m]["cell27"][r]["max_abs_diff_methane_vs_stored"] for m in mids if out[m]["in_stored"]) for r in ("0.07", "0.33")}
    print("max |dCH4| over all 20 stored models x 6 cores:", allmax)
