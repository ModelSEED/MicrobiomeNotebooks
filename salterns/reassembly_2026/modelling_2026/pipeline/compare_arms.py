"""Pre-registered annotation comparison: metrics for every arm and the decision rule of ../ANNOTATION_CHOICE_CRITERIA.md.

Re-runnable; reads the build reports in ../annotation_comparison/runs/<arm>/*.report.json and writes
  ../annotation_comparison/per_model_metrics.tsv   (one row per model x arm, every metric)
  ../annotation_comparison/summary.json            (medians, criteria, decision)
The prose RESULTS.md is written by hand from these two files.

ARITHMETIC, FIXED BEFORE THE FULL RUN WAS INSPECTED (2026-11-20, after 4 + 3 test builds only):
- Unit ("model"): one build of one species representative under one guild protocol. 57 species give 58 models,
  because R2_sediment__metabat2.bin.434_sub (UBA8670) is in both the MePn and the betaine guild and is built under
  both protocols. "The median over the 57 models" is therefore taken over these 58 builds.
- Paired set: medians use the models whose build finished (status ok) in all four arms. Builds that errored are
  listed per arm and excluded from every arm's medians. Models that do not grow after gap-filling stay in the
  medians (gap-filling was attempted; its count is still meaningful) and are counted per arm.
- (a) median(m1, arm X) >= 1.10 x median(m1, arm A). m1 = reactions with a gene association in the draft that are in
      the final model and were not gap-filled (hand-added protocol reactions are not in the draft, so not counted).
- (b) median(m2, X) <= 0.80 x median(m2, A). m2 = reactions added or reversed by the guild-medium gap-fills (all stages);
      the ATP-correction core gap-fill is reported separately (m2_atp_core_gapfill) and not counted.
- (c) for each guild: median of m3 (fraction of the guild pathway's reactions gene-associated in the draft) in arm X
      >= that in arm A. Guild means are reported too, but the median decides.
- (d) no model that passes the free-ATP test (ATPM = 0 on all three carbon-free media) in arm A fails it in arm X,
      AND median ATP-test pass rate (X) >= median (A). ATP tests are built by each arm's own ATP correction, so the
      number of tests is reported next to the pass rates.
- (e) median blocked fraction (X) - median blocked fraction (A) <= 0.02 (blocked = zero FVA range with all exchanges
      open, over internal, non-boundary reactions). Dead-end metabolite fractions are reported, not decided on.
- Decision: the arm with the most added evidence (D > C > B) that meets all of (a)-(e); otherwise A.
- Ratios of medians decide. Medians of per-model paired changes are reported alongside; a disagreement between the
  two is flagged in summary.json, not used to override.
usage: python compare_arms.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
AC = HERE.parent / "annotation_comparison"
ARMS = ["A", "B", "C", "D"]
ORDER_RICHEST_FIRST = ["D", "C", "B"]


def load():
    rows = []
    for arm in ARMS:
        for p in sorted((AC / "runs" / arm).glob("*.report.json")):
            r = json.load(open(p))
            base = {"arm": arm, "mag_id": r["mag_id"], "guild": r["guild"], "status": r["status"],
                    "error": r.get("error", "")}
            if r["status"] != "ok":
                rows.append(base)
                continue
            m, q, d = r["metrics"], r["qc"], r["draft"]
            rows.append({**base, "genus": r.get("genus"), "phylum": r.get("phylum"), "template": r["template_file"],
                         "draft_reactions": d["n_reactions"], "draft_gene_associated": d["n_gene_associated"],
                         "draft_rast_supported": d["n_rast_supported"], "draft_extra_only": d["n_extra_only"],
                         "m1": m["m1_gene_associated_non_gapfilled"], "m2": m["m2_gapfilled"],
                         "m2_atp_core_gapfill": m.get("m2_atp_core_gapfill"),
                         "m2_by_stage": json.dumps(m["m2_gapfilled_by_stage"], sort_keys=True),
                         "m3": m["m3_guild_pathway_fraction"],
                         "m3_detail": json.dumps(r["guild_pathway"]["gene_supported"], sort_keys=True),
                         "free_atp_pass": m["m4_free_atp_pass"], "free_atp": json.dumps(q["free_atp"]),
                         "draft_free_atp_max": max(r["draft_free_atp_before_atp_correction"].values()),
                         "atp_test_pass_rate": m["m4_atp_test_pass_rate"], "n_atp_tests": len(q["atp_tests"]),
                         "blocked_frac": m["m5_blocked_frac"], "dead_end_met_frac": m["m5_dead_end_metabolite_frac"],
                         "n_internal": q["n_internal_reactions"], "m6_weak": m["m6_weak_only_reactions"],
                         "extra_only_final": m["extra_only_reactions_in_final"],
                         "extra_only_by_type": json.dumps(m.get("extra_only_by_template_type", {}), sort_keys=True),
                         "grows": m["grows"], "guild_growth": q["guild_growth"],
                         "final_reactions": r["final"]["n_reactions"], "final_genes": r["final"]["n_genes"],
                         "obligate": (q.get("growth_triplet") or {}).get("obligate"),
                         "wall_s": sum(r["times_s"].values())})
    return pd.DataFrame(rows)


def fmt(v, nd=3):
    if isinstance(v, bool):
        return "yes" if v else "**no**"
    if isinstance(v, float):
        return f"{v:.{nd}f}"
    return str(v)


def write_tables(out, P):
    """Markdown tables for RESULTS.md (written to tables.md)."""
    L = []
    med, mean = out["medians"], out["means"]
    L.append("### Metric medians (means in brackets) over the paired models\n")
    L.append("| Metric | " + " | ".join(f"Arm {a}" for a in ARMS) + " |")
    L.append("|---|" + "---:|" * len(ARMS))
    names = [("draft_gene_associated", "Draft reactions with a gene association", 0),
             ("m1", "1. Gene-associated, non-gap-filled reactions (final model)", 0),
             ("m2", "2. Gap-filled reactions, guild-medium gap-fills", 0),
             ("m2_atp_core_gapfill", "2'. ATP-correction core gap-fill (not in rule)", 1),
             ("m3", "3. Guild-pathway gene support (fraction)", 3),
             ("atp_test_pass_rate", "4. ATP-test pass rate", 3), ("n_atp_tests", "4'. ATP tests per model", 1),
             ("blocked_frac", "5. Blocked reaction fraction", 3), ("dead_end_met_frac", "5'. Dead-end metabolite fraction", 3),
             ("m6_weak", "6. Reactions with weak single-source non-RAST evidence only", 0),
             ("extra_only_final", "Reactions with non-RAST evidence only (final model)", 0),
             ("final_reactions", "Final model reactions", 0), ("guild_growth", "Growth on guild medium (1/h)", 3)]
    for k, lab, nd in names:
        L.append(f"| {lab} | " + " | ".join(f"{med[a][k]:.{nd}f} ({mean[a][k]:.{nd}f})" for a in ARMS) + " |")
    L.append("| 4. Models passing the free-ATP test | " + " | ".join(f"{med[a]['free_atp_pass_n']}/{out['n_paired_models']}" for a in ARMS) + " |")
    L.append("| Models growing on the guild medium | " + " | ".join(f"{med[a]['grows_n']}/{out['n_paired_models']}" for a in ARMS) + " |")
    L.append("| MePn-obligate models | " + " | ".join(f"{med[a]['obligate_n']}" for a in ARMS) + " |")
    L.append("\n### Per guild (medians; m3 mean in brackets)\n")
    L.append("| Guild | n | " + " | ".join(f"m3 {a}" for a in ARMS) + " | " + " | ".join(f"m1 {a}" for a in ARMS) +
             " | " + " | ".join(f"m2 {a}" for a in ARMS) + " | " + " | ".join(f"blocked {a}" for a in ARMS) + " |")
    L.append("|---|---:|" + "---:|" * 16)
    for g in sorted(out["per_guild"]["A"]):
        pg = {a: out["per_guild"][a][g] for a in ARMS}
        L.append(f"| {g} | {pg['A']['n']} | " + " | ".join(f"{pg[a]['m3_median']:.2f} ({pg[a]['m3_mean']:.2f})" for a in ARMS)
                 + " | " + " | ".join(f"{pg[a]['m1_median']:.0f}" for a in ARMS) + " | "
                 + " | ".join(f"{pg[a]['m2_median']:.0f}" for a in ARMS) + " | "
                 + " | ".join(f"{pg[a]['blocked_median']:.3f}" for a in ARMS) + " |")
    L.append("\n### Decision rule, each merged arm against arm A\n")
    L.append("| Criterion | Arm B | Arm C | Arm D |")
    L.append("|---|---|---|---|")
    cr = out["criteria"]
    L.append("| (a) m1 ratio to A >= 1.10 | " + " | ".join(f"{cr[a]['a_m1_ratio']:.3f} {fmt(cr[a]['a_pass'])}" for a in "BCD") + " |")
    L.append("| (b) m2 ratio to A <= 0.80 | " + " | ".join(f"{cr[a]['b_m2_ratio']:.3f} {fmt(cr[a]['b_pass'])}" for a in "BCD") + " |")
    L.append("| (c) guild support not lower, every guild | " + " | ".join(
        ", ".join(f"{g.split('_')[0]} {fmt(v)}" for g, v in cr[a]["c_per_guild"].items()) for a in "BCD") + " |")
    L.append("| (d) no model gains free ATP; ATP-test pass rate not lower | " + " | ".join(
        f"{len(cr[a]['d_models_gaining_free_atp'])} gain; {cr[a]['d_atp_pass_rate'][1]:.3f} vs {cr[a]['d_atp_pass_rate'][0]:.3f} {fmt(cr[a]['d_pass'])}" for a in "BCD") + " |")
    L.append("| (e) blocked fraction change <= +2 pp | " + " | ".join(f"{cr[a]['e_blocked_change_pp']:+.2f} pp {fmt(cr[a]['e_pass'])}" for a in "BCD") + " |")
    L.append("| **All met** | " + " | ".join(fmt(cr[a]["all_pass"]) for a in "BCD") + " |")
    L.append("| Paired view (median per-model change): m1 / m2 / blocked | " + " | ".join(
        f"{out['paired_change_medians'][a]['m1_rel_change']:+.1%} / {out['paired_change_medians'][a]['m2_rel_change']:+.1%} / {100*out['paired_change_medians'][a]['blocked_frac_change']:+.2f} pp" for a in "BCD") + " |")
    L.append(f"\n**Decision: arm {out['decision']}.**\n")
    L.append("### Per model\n")
    L.append("m1 = gene-associated non-gap-filled reactions; m2 = gap-filled; m3 = guild-pathway gene support; "
             "blk = blocked fraction; ATP = free-ATP test passed / ATP tests passed; grow = growth on guild medium.\n")
    hdr = ["Guild", "MAG", "Genus"] + [f"{c} {a}" for c in ("m1", "m2", "m3", "blk") for a in ARMS] + [f"ATP {a}" for a in ARMS] + [f"grow {a}" for a in ARMS]
    L.append("| " + " | ".join(hdr) + " |")
    L.append("|" + "---|" * 3 + "---:|" * (len(hdr) - 3))
    W = P.set_index(["guild", "mag_id", "arm"])
    for (g, mag), t in P.groupby(["guild", "mag_id"]):
        t = t.set_index("arm")
        cells = [g.split("_")[0], f"`{mag}`", str(t.loc["A", "genus"])]
        cells += [f"{int(t.loc[a, 'm1'])}" for a in ARMS] + [f"{int(t.loc[a, 'm2'])}" for a in ARMS]
        cells += [f"{t.loc[a, 'm3']:.2f}" for a in ARMS] + [f"{t.loc[a, 'blocked_frac']:.3f}" for a in ARMS]
        cells += [("ok" if t.loc[a, "free_atp_pass"] else "**LOOP**") + f" {round(t.loc[a, 'atp_test_pass_rate'] * t.loc[a, 'n_atp_tests'])}/{t.loc[a, 'n_atp_tests']}" for a in ARMS]
        cells += [f"{t.loc[a, 'guild_growth']:.3f}" if t.loc[a, "grows"] else "**0**" for a in ARMS]
        L.append("| " + " | ".join(cells) + " |")
    (AC / "tables.md").write_text("\n".join(L) + "\n")


def main():
    df = load()
    df.to_csv(AC / "per_model_metrics.tsv", sep="\t", index=False)
    ok = df[df.status == "ok"]
    key = ["mag_id", "guild"]
    n_ok = ok.groupby(key).arm.nunique()
    paired_keys = set(n_ok[n_ok == len(ARMS)].index)
    P = ok[ok.set_index(key).index.isin(paired_keys)]
    out = {"n_models_expected_per_arm": int(df.groupby("arm").size().max()),
           "n_paired_models": len(paired_keys),
           "errors_per_arm": {a: df[(df.arm == a) & (df.status != "ok")][key + ["error"]].values.tolist() for a in ARMS},
           "non_growing_per_arm": {a: P[(P.arm == a) & (~P.grows.astype(bool))][key].values.tolist() for a in ARMS},
           "medians": {}, "means": {}, "per_guild": {}, "paired_change_medians": {}, "criteria": {}}
    metrics = ["m1", "m2", "m2_atp_core_gapfill", "m3", "atp_test_pass_rate", "n_atp_tests", "blocked_frac",
               "dead_end_met_frac", "m6_weak", "extra_only_final", "draft_gene_associated", "final_reactions",
               "guild_growth"]
    for a in ARMS:
        s = P[P.arm == a]
        out["medians"][a] = {k: float(np.median(s[k].astype(float))) for k in metrics}
        out["means"][a] = {k: float(np.mean(s[k].astype(float))) for k in metrics}
        out["medians"][a]["free_atp_pass_n"] = int(s.free_atp_pass.astype(bool).sum())
        out["medians"][a]["grows_n"] = int(s.grows.astype(bool).sum())
        out["medians"][a]["obligate_n"] = int((s.obligate == True).sum())  # noqa: E712
        out["per_guild"][a] = {g: {"n": int(len(t)), "m3_median": float(np.median(t.m3)), "m3_mean": float(np.mean(t.m3)),
                                   "m1_median": float(np.median(t.m1)), "m2_median": float(np.median(t.m2)),
                                   "blocked_median": float(np.median(t.blocked_frac))}
                               for g, t in s.groupby("guild")}
    A = P[P.arm == "A"].set_index(key)
    for a in ORDER_RICHEST_FIRST:
        X = P[P.arm == a].set_index(key).loc[A.index]
        mA, mX = out["medians"]["A"], out["medians"][a]
        ca = mX["m1"] >= 1.10 * mA["m1"]
        cb = mX["m2"] <= 0.80 * mA["m2"]
        cc_g = {g: out["per_guild"][a][g]["m3_median"] >= out["per_guild"]["A"][g]["m3_median"]
                for g in out["per_guild"]["A"]}
        cc = all(cc_g.values())
        gained_loops = [list(k) for k in A.index if bool(A.loc[k, "free_atp_pass"]) and not bool(X.loc[k, "free_atp_pass"])]
        cd = (len(gained_loops) == 0) and (mX["atp_test_pass_rate"] >= mA["atp_test_pass_rate"])
        ce = (mX["blocked_frac"] - mA["blocked_frac"]) <= 0.02
        rel = lambda c: float(np.median((X[c].astype(float) - A[c].astype(float)) / A[c].astype(float).replace(0, np.nan)))  # noqa: E731
        pc = {"m1_rel_change": rel("m1"), "m2_rel_change": rel("m2"),
              "blocked_frac_change": float(np.median(X.blocked_frac - A.blocked_frac))}
        out["paired_change_medians"][a] = pc
        out["criteria"][a] = {
            "a_m1_ratio": mX["m1"] / mA["m1"], "a_pass": bool(ca),
            "b_m2_ratio": (mX["m2"] / mA["m2"]) if mA["m2"] else None, "b_pass": bool(cb),
            "c_per_guild": cc_g, "c_pass": bool(cc),
            "d_models_gaining_free_atp": gained_loops, "d_atp_pass_rate": [mA["atp_test_pass_rate"], mX["atp_test_pass_rate"]],
            "d_pass": bool(cd),
            "e_blocked_change_pp": 100 * (mX["blocked_frac"] - mA["blocked_frac"]), "e_pass": bool(ce),
            "all_pass": bool(ca and cb and cc and cd and ce),
            "paired_view_disagrees": {"a": bool((pc["m1_rel_change"] >= 0.10) != ca),
                                      "b": bool((pc["m2_rel_change"] <= -0.20) != cb),
                                      "e": bool((pc["blocked_frac_change"] <= 0.02) != ce)}}
    winner = next((a for a in ORDER_RICHEST_FIRST if out["criteria"][a]["all_pass"]), "A")
    out["decision"] = winner
    json.dump(out, open(AC / "summary.json", "w"), indent=1, default=str)
    write_tables(out, P)
    print(json.dumps({k: out[k] for k in ("n_paired_models", "decision", "criteria", "medians")}, indent=1, default=str))


if __name__ == "__main__":
    main()
