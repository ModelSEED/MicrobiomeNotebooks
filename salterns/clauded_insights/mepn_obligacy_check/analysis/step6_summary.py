"""Step 6: collate steps 0-5 into one per-model table (all 24 models) + the corrected totals.
Outputs: summary_per_model.csv, summary.json"""
import csv
import json
import statistics
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, "/tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns/a24fbb2e-9d33-5e05-aeea-a07d777ccde7/scratchpad/mepn_verify/analysis")
import vcommon as V  # noqa: E402

O = V.OUT
s1 = json.load(open(O / "step1_growth_and_repro.json"))
s1b = json.load(open(O / "step1b_fva.json"))
s4c = json.load(open(O / "step4c_repair_and_rerun.json"))
s5 = json.load(open(O / "step5_consistent_basis.json"))
prim = json.load(open(O / "step4_primary_repairs.json"))
prov = json.load(open(O / "step0_provenance.json"))
stored = V.stored()
rows = []
for mid in V.model_ids():
    d = s1[mid]
    iid = d["iterative_id"]
    row = {"model": mid, "iterative_id": iid, "obligate": d["obligate"], "grows_cell5_with_mepn": d["a_as_is"],
           "grows_cell5_without_mepn": d["b_MePn_blocked"], "grows_pi_free_mepn_only": d["c_Pi_blocked_MePn_only_P"],
           "in_saltern_fluxes_json": d["in_stored"], "n_gapfilled_reactions_cell93": len(prov[mid]["gapfilled"])}
    if d["in_stored"]:
        for r in ("0.07", "0.33"):
            ks = [stored[r]["methane"][mid][c] / d["cell27"][r]["cores"][c]["abundance"] for c in V.abundances_as_run().columns]
            row[f"ch4_per_abundance_stored_r{r}"] = statistics.mean(ks)
            row[f"ch4_per_abundance_stored_r{r}_range"] = [min(ks), max(ks)]
            row[f"max_abs_diff_rerun_vs_stored_r{r}"] = d["cell27"][r]["max_abs_diff_methane_vs_stored"]
            row[f"ratio_binds_all_cores_r{r}"] = all(bool(x) for x in d["cell27"][r]["ratio_binds"].values())
            row[f"MePn_over_Pi_r{r}"] = statistics.mean(v for v in d["cell27"][r]["MePn_over_Pi"].values())
        if iid in prim:
            x = s4c[iid]["+".join(prim[iid])]
            row["repair_reactions"] = "+".join(prim[iid])
            row["repaired_growth_a_b_c"] = [x["a_as_is"], x["b_MePn_blocked"], x["c_Pi_blocked_MePn_only_P"]]
            row["repaired_ATPM_no_carbon"] = x["ATPM_max_no_carbon"]["repaired"]
            for r in ("0.07", "0.33"):
                row[f"ch4_per_abundance_repaired_r{r}"] = statistics.mean(x["cell27"][r]["k_methane_per_abundance"].values())
                row[f"MePn_over_Pi_repaired_r{r}"] = statistics.mean(x["cell27"][r]["MePn_over_Pi"].values())
        else:
            for r in ("0.07", "0.33"):
                row[f"ch4_per_abundance_repaired_r{r}"] = row[f"ch4_per_abundance_stored_r{r}"]
    if iid in s1b:
        row["fva_ch4_per_abundance_r0.07"] = [min(c["CH4_per_abundance_range"][0] for c in s1b[iid]["fva"]["0.07"].values()),
                                              max(c["CH4_per_abundance_range"][1] for c in s1b[iid]["fva"]["0.07"].values())]
    rows.append(row)
keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("iterative_id", "model", "obligate"), k))
with open(O / "summary_per_model.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys)
    w.writeheader()
    for r in rows:
        w.writerow({k: (json.dumps(v) if isinstance(v, (list, dict)) else v) for k, v in r.items()})
tot = {r: {lab: {k: s5["areal"][r][lab][k] for k in ("model_over_measured_total", "model_over_measured_min",
                                                      "model_over_measured_max", "spearman_rho", "spearman_p",
                                                      "pearson_r2", "share_of_three", "sim_by_core")}
           for lab in ("before", "after")} for r in ("0.07", "0.33")}
json.dump({"per_model": rows, "consistent_basis_areal": tot, "regression": s5["regression"],
           "alternatives_range_areal": s5["alternatives_range_areal"],
           "relab_only": {r: {lab: {k: s5["relab_only"][r][lab][k] for k in ("model_over_measured_total", "spearman_rho",
                                                                             "share_of_three")}
                              for lab in ("before", "after")} for r in ("0.07", "0.33")}},
          open(O / "summary.json", "w"), indent=1, default=float)
for r in rows:
    print(f"{r['iterative_id']:22s} oblig={str(r['obligate'])[0]} a/b/c {r['grows_cell5_with_mepn']:.3f}/{r['grows_cell5_without_mepn']:.3f}/"
          f"{r['grows_pi_free_mepn_only']:.3f} stored={str(r['in_saltern_fluxes_json'])[0]} "
          + (f"k07 {r['ch4_per_abundance_stored_r0.07']:.6f} ({r['ch4_per_abundance_stored_r0.07_range'][0]:.6f}-{r['ch4_per_abundance_stored_r0.07_range'][1]:.6f}) "
             f"k33 {r['ch4_per_abundance_stored_r0.33']:.6f} dmax {max(r['max_abs_diff_rerun_vs_stored_r0.07'], r['max_abs_diff_rerun_vs_stored_r0.33']):.1e} "
             f"-> repaired k07 {r['ch4_per_abundance_repaired_r0.07']:.6f} k33 {r['ch4_per_abundance_repaired_r0.33']:.6f} {r.get('repair_reactions', '')}"
             if r["in_saltern_fluxes_json"] else ""))
print(json.dumps(tot, indent=1, default=float)[:3000])
