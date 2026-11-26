#!/usr/bin/env python
"""Comparison tables for the 67 grouped GROW samples (run at the end; safe to run any time, missing sides are NA).

Reads the local CheckM collections (collect_checkm_select.py output, checkm_<tag>.tsv: completeness, contamination,
n50, size_bp, pass_50_10 per bin):
  per-sample  ../assembly_binning_2026/work/checkm/<sample>/checkm_<sample>.tsv
  sitecov     work/checkm/<sample>.sitecov/checkm_<sample>.sitecov.tsv
  coassembly  work/checkm/<group>.coassembly/checkm_<group>.coassembly.tsv
Writes comparison/:
  per_sample_vs_sitecov.tsv      one row per grouped sample: MAGs >= 50/<10 per-sample vs sitecov, and their medians
  coassembly_vs_sitecov_union.tsv  one row per co-assembly pilot group: co-assembly MAGs vs the union of members' sitecov
  distributions.tsv              completeness / contamination / N50 quartiles of the >= 50/<10 MAGs per set
MAG sets are never merged or dereplicated together; counts are compared side by side.
Usage: ~/Documents/py_venv/bin/python compare.py
"""
import csv
import os
import statistics
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
PS_CHECKM = os.path.normpath(os.path.join(HERE, "..", "assembly_binning_2026", "work", "checkm"))
CHECKM = os.path.join(HERE, "work", "checkm")
OUT = os.path.join(HERE, "comparison")


def mags(base, tag):
    """The >= 50/<10 bins of one unit, or None if its CheckM table is not there (yet)."""
    p = os.path.join(base, tag, f"checkm_{tag}.tsv")
    if not os.path.exists(p):
        return None
    with open(p, newline="") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if r["pass_50_10"] == "1"]


def med(rows, col):
    vals = [float(r[col]) for r in rows or [] if r.get(col) not in (None, "", "None")]
    return round(statistics.median(vals), 2) if vals else "NA"


def quart(rows, col):
    vals = sorted(float(r[col]) for r in rows if r.get(col) not in (None, "", "None"))
    if not vals:
        return ["NA"] * 5
    q = statistics.quantiles(vals, n=4) if len(vals) > 1 else [vals[0]] * 3
    return [round(x, 2) for x in (vals[0], q[0], q[1], q[2], vals[-1])]


def main():
    with open(os.path.join(HERE, "manifest.tsv"), newline="") as fh:
        man = list(csv.DictReader(fh, delimiter="\t"))
    os.makedirs(OUT, exist_ok=True)
    sets = defaultdict(list)
    rows = []
    for r in man:
        s = r["sample"]
        ps, sc = mags(PS_CHECKM, s), mags(CHECKM, f"{s}.sitecov")
        sets["per-sample (67 grouped samples)"] += ps or []
        sets["sitecov"] += sc or []
        rows.append({"site_group": r["site_group"], "sample": s, "per_sample_mags": len(ps) if ps is not None else "NA",
                     "sitecov_mags": len(sc) if sc is not None else "NA",
                     "delta": (len(sc) - len(ps)) if ps is not None and sc is not None else "NA",
                     "per_sample_median_completeness": med(ps, "completeness"),
                     "sitecov_median_completeness": med(sc, "completeness"),
                     "per_sample_median_contamination": med(ps, "contamination"),
                     "sitecov_median_contamination": med(sc, "contamination"),
                     "per_sample_median_n50": med(ps, "n50"), "sitecov_median_n50": med(sc, "n50")})
    with open(os.path.join(OUT, "per_sample_vs_sitecov.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    groups = defaultdict(list)
    for r in man:
        if r["coassembly_pilot"] == "True":
            groups[r["site_group"]].append(r["sample"])
    grows = []
    for g, members in sorted(groups.items()):
        ca = mags(CHECKM, f"{g}.coassembly")
        msc = [mags(CHECKM, f"{s}.sitecov") for s in members]
        union = [x for m in msc if m for x in m]
        sets["coassembly"] += ca or []
        sets["sitecov union, co-assembly groups"] += union
        grows.append({"site_group": g, "members": len(members), "coassembly_mags": len(ca) if ca is not None else "NA",
                      "sitecov_union_mags": len(union) if all(m is not None for m in msc) else "NA",
                      "coassembly_median_completeness": med(ca, "completeness"),
                      "sitecov_union_median_completeness": med(union, "completeness"),
                      "coassembly_median_contamination": med(ca, "contamination"),
                      "sitecov_union_median_contamination": med(union, "contamination"),
                      "coassembly_median_n50": med(ca, "n50"), "sitecov_union_median_n50": med(union, "n50")})
    with open(os.path.join(OUT, "coassembly_vs_sitecov_union.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(grows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(grows)
    with open(os.path.join(OUT, "distributions.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["set", "n_mags", "metric", "min", "q1", "median", "q3", "max"])
        for name, rs in sets.items():
            for col in ("completeness", "contamination", "n50"):
                w.writerow([name, len(rs), col] + quart(rs, col))
    done_ps = sum(1 for r in rows if r["per_sample_mags"] != "NA")
    done_sc = sum(1 for r in rows if r["sitecov_mags"] != "NA")
    print(f"per-sample CheckM collected for {done_ps}/67, sitecov for {done_sc}/67, co-assembly for "
          f"{sum(1 for r in grows if r['coassembly_mags'] != 'NA')}/8 -> {OUT}")


if __name__ == "__main__":
    main()
