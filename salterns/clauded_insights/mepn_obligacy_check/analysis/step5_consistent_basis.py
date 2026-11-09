"""Step 5: simulated vs measured methane on the consistent areal basis, before vs after repairing the three
MePn-obligate models.  Reuses simulated(), fit(), biomass() from
salterns/clauded_insights/modelling_consistent_basis_r232.py (imported by path, no bytecode written; its main() is
not run).  'before' = stored saltern_fluxes.json; 'after' = the same 20 models with the primary repaired models'
cell-27 methane (step4c, PRIMARY) substituted for Roseovarius.9, Roseovarius.11 and Albimonas.1.
Regression check: 'before' must reproduce modelling_consistent_basis_r232.json exactly.
Output: step5_consistent_basis.json"""
import importlib.util
import json
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, "/tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns/a24fbb2e-9d33-5e05-aeea-a07d777ccde7/scratchpad/mepn_verify/analysis")
import vcommon as V  # noqa: E402

CB_PATH = V.S / "clauded_insights" / "modelling_consistent_basis_r232.py"
spec = importlib.util.spec_from_file_location("mcb", str(CB_PATH))
mcb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mcb)

THREE = ("Roseovarius.9", "Roseovarius.11", "Albimonas.1")
CORES = mcb.CORES


def summary(fl, dna=None):
    sim, per = mcb.simulated(fl, 1.0, dna)
    f = mcb.fit(sim)
    three = [c for c in per.columns if c in THREE]
    f["share_of_three"] = float(per[three].values.sum() / per.values.sum()) if three else 0.0
    f["three_by_model_total"] = {c: float(per[c].sum()) for c in three}
    f["per_model_total"] = {c: float(per[c].sum()) for c in per.columns}
    f["measured_total"] = float(mcb.CH4.sum())
    f["sim_total"] = float(sim.sum())
    return f


def main():
    stored = json.load(open(V.STORED))
    old = json.load(open(V.S / "clauded_insights" / "modelling_consistent_basis_r232.json"))
    rep = json.load(open(V.OUT / "step4c_repair_and_rerun.json"))
    primary = json.load(open(V.OUT / "step4_primary_repairs.json"))
    inv = {V.iterative_id(m): m for m in V.model_ids()}
    out = {"basis": "areal (consistent basis of modelling_consistent_basis_r232.py) and relab_only sensitivity",
           "primary_repairs": primary, "regression": {}, "areal": {}, "relab_only": {}, "alternatives_areal": {}}
    mean_dna = mcb.VW.mean()
    for r in ("0.07", "0.33"):
        before = stored[r]["methane"]
        after = dict(before)
        for iid in THREE:
            mid = inv[iid]
            x = rep[iid]["+".join(primary[iid])]
            after[mid] = {c: x["cell27"][r]["cores"][c]["methane"] for c in CORES}
        b, a = summary(before), summary(after)
        name = f"MePn_{r}_20"
        reg = {k: (b[k], old["areal"][name][k]) for k in ("model_over_measured_total", "spearman_rho", "model_over_measured_min",
                                                           "model_over_measured_max", "pearson_r2")}
        out["regression"][name] = {"before_vs_stored_json": reg,
                                   "max_abs_diff": max(abs(u - v) for u, v in reg.values())}
        out["areal"][r] = {"before": b, "after": a}
        out["relab_only"][r] = {"before": summary(before, mean_dna), "after": summary(after, mean_dna)}
        # every alternative single-reaction repair (robustness of the corrected total)
        # every alternative repair that removes the obligacy, one model varied at a time (others at their primary):
        # 'metabolic' = template additions that are not transporters (step4c), the RAST-supported single reversals and
        # the MSGapfill sets (step4d); 'export_sink' = transporters run as exporters with an export-only exchange
        # (step4d), i.e. the model dumps a C1-carrying metabolite.
        more = json.load(open(V.OUT / "step4d_more_repairs.json"))
        cands = []
        for iid in THREE:
            for repname, x in rep[iid].items():
                nm = x["added"][0]["equation"]
                is_tr = any(m.endswith("[e0]") for m in nm.replace("<=>", "+").replace("-->", "+").replace("<--", "+").split("+"))
                if x["b_MePn_blocked"] > 1e-6 and not is_tr:
                    cands.append((iid, repname, "metabolic", x["cell27"]))
        for x in more:
            if x["kind"] == "rast_additions_only" or not x["b_MePn_blocked"] > 1e-6:
                continue
            cls = "metabolic" if (x["kind"] == "reversal" or "+" in x["label"] or x["label"] in ("rxn00336_c0", "rxn00770_c0")) else "export_sink"
            cands.append((x["iterative_id"], x["label"], cls, x["cell27"]))
        alts = {}
        for iid, repname, cls, c27 in cands:
            alt = dict(after)
            alt[inv[iid]] = {c: c27[r]["cores"][c]["methane"] for c in CORES}
            s = summary(alt)
            alts[f"{iid}:{repname}"] = {"class": cls, "k_methane_per_abundance_R1_A": c27[r]["k_methane_per_abundance"]["R1_A"],
                                        **{k: s[k] for k in ("model_over_measured_total", "spearman_rho", "share_of_three")}}
        out["alternatives_areal"][r] = alts
        out.setdefault("alternatives_range_areal", {})[r] = {
            cls: [min(v["model_over_measured_total"] for v in alts.values() if v["class"] == cls),
                  max(v["model_over_measured_total"] for v in alts.values() if v["class"] == cls),
                  sum(1 for v in alts.values() if v["class"] == cls)]
            for cls in ("metabolic", "export_sink")}
    json.dump(out, open(V.OUT / "step5_consistent_basis.json", "w"), indent=1, default=float)
    for r in ("0.07", "0.33"):
        print(f"r = {r}: regression max |diff| vs modelling_consistent_basis_r232.json = {out['regression'][f'MePn_{r}_20']['max_abs_diff']:.2e}")
        for basis in ("areal", "relab_only"):
            for lab in ("before", "after"):
                s = out[basis][r][lab]
                print(f"   {basis:10s} {lab:6s} model/measured total {s['model_over_measured_total']:.4f} "
                      f"(per core {s['model_over_measured_min']:.3f}-{s['model_over_measured_max']:.3f}) rho {s['spearman_rho']:+.4f} "
                      f"(p {s['spearman_p']:.3f}) R2 {s['pearson_r2']:.3f} share of 3 {s['share_of_three']:.4f} sim_total {s['sim_total']:.1f} "
                      f"by core {s['sim_by_core']}")
        for cls, (lo, hi, n) in out["alternatives_range_areal"][r].items():
            print(f"   alternatives [{cls}, n={n}] (one model's repair varied at a time): model/measured total {lo:.4f}-{hi:.4f}")


if __name__ == "__main__":
    main()
