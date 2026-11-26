#!/usr/bin/env python
"""Summarise results/*.json from validate_build.py as Markdown tables (stdout) and
results/summary.json."""
import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
rows = {}
for p in sorted(glob.glob(os.path.join(HERE, "results", "*__*.json"))):
    d = json.load(open(p))
    rows[(d["label"].split("__")[0], d["arm"])] = d

ORDER = ["Roseovarius_R2BD2_mb38", "BM714_R1AD1_mb51", "Methanolobus_R1AD2_cc21"]
MEDIA = ["water_H", "mineral_O2", "mineral_anox"]


def fa(stage):
    if not stage or "free_atp" not in stage:
        return "n/a"
    return " / ".join(f"{stage['free_atp'][m]['atp']:g}" for m in MEDIA)


summary = {}
print("| MAG | template | arm | draft rxns (gene-assoc.) | free ATP draft | ATP media w/o gapfill | selected ATP media | core gapfill rxns | expansion-filtered | free ATP after ATP corr. | gapfill new / reversed | growth, complete | ATP tests pass (final) | free ATP gapfilled | final rxns |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for mag in ORDER:
    for arm in ("upstream", "combined"):
        d = rows.get((mag, arm))
        if not d:
            print(f"| {mag} | | {arm} | missing |")
            continue
        s = d["stages"]
        dr, ac, gf = s.get("draft", {}), s.get("atp_correction", {}), s.get("gapfilled", {})
        tests = gf.get("test_results") or []
        rec = {
            "draft_reactions": dr["counts"]["reactions"], "draft_gene_assoc": dr["counts"]["gene_associated"],
            "free_atp_draft": {m: dr["free_atp"][m]["atp"] for m in MEDIA},
            "atp_media_without_gapfill": len(ac.get("media_making_atp_without_gapfill", [])),
            "media_evaluated": ac.get("media_evaluated"),
            "selected_media": ac.get("selected_media"),
            "core_gapfilled": len(ac.get("core_gapfilled_reactions", [])),
            "filtered_noncore": len(ac.get("filtered_noncore", [])),
            "free_atp_after_atp_correction": {m: ac["free_atp"][m]["atp"] for m in MEDIA} if ac else None,
            "gapfill_status": gf.get("status"),
            "gapfill_new": len(gf.get("new", {})), "gapfill_reversed": len(gf.get("reversed", {})),
            "growth_complete": gf.get("growth_complete"),
            "atp_tests_pass": f"{sum(t['pass'] for t in tests)}/{len(tests)}" if tests else None,
            "free_atp_gapfilled": {m: gf["free_atp"][m]["atp"] for m in MEDIA} if gf.get("free_atp") else None,
            "final_reactions": gf.get("counts", {}).get("reactions"),
            "seconds": d.get("timing", {}).get("total"),
        }
        summary[f"{mag}__{arm}"] = rec
        print(f"| {mag} | {d['template_kind']} | {arm} | {rec['draft_reactions']} ({rec['draft_gene_assoc']}) | {fa(dr)} "
              f"| {rec['atp_media_without_gapfill']}/{rec['media_evaluated']} | {len(rec['selected_media'] or [])} "
              f"| {rec['core_gapfilled']} | {rec['filtered_noncore']} | {fa(ac)} | {rec['gapfill_new']} / {rec['gapfill_reversed']} "
              f"| {rec['growth_complete']} | {rec['atp_tests_pass']} | {fa(gf)} | {rec['final_reactions']} |")

# draft reaction differences between arms
print()
for mag in ORDER:
    a, b = rows.get((mag, "upstream")), rows.get((mag, "combined"))
    if not (a and b):
        continue
    ra, rb = set(a["stages"]["draft"]["reaction_ids"]), set(b["stages"]["draft"]["reaction_ids"])
    summary[f"{mag}__draft_diff"] = {"only_upstream": sorted(ra - rb), "only_combined": sorted(rb - ra)}
    print(f"{mag}: draft reactions only upstream {sorted(ra - rb)}; only combined {sorted(rb - ra)}")
    for arm, d in (("upstream", a), ("combined", b)):
        for st in ("draft", "atp_correction", "gapfilled"):
            fa_ = d["stages"].get(st, {}).get("free_atp", {})
            for m, v in fa_.items():
                if v["atp"] > 1e-6:
                    act = list(v.get("active", {}))
                    print(f"   leak {arm}/{st}/{m}: ATP {v['atp']}, {v.get('n_active')} active reactions, top: {act[:25]}")
json.dump(summary, open(os.path.join(HERE, "results", "summary.json"), "w"), indent=1)
