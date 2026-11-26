#!/usr/bin/env python
"""Step 12 (v2): regression check v1 -> v2 for every label (and merged class) present in both surveys.
  results/regression_v1_v2_per_sample.tsv  per sample x label: v1 / v2 reads, ambiguous_best, copies_per_GE, rel. deviation
  results/regression_v1_v2_summary.tsv     per label: median and max |relative deviation| of copies_per_GE over the
                                           samples (all samples with v1 reads > 0, and samples with >= 20 v1 reads),
                                           pooled 39-sample ratio, change in reads and in ambiguous best hits
  results/regression_genome_equivalents.tsv  GE depth v1 vs v2 per sample
  results/regression_mcrA_genus_pooled.tsv   mcrA reads of the 12 unrestored-core samples by genus, v1 vs v2 shares
                                             (v1 genus names harmonized to GTDB R232 with the v2 KEGG-organism map)
  results/regression_lost_reads_<sample>.tsv  (optional, `attribute <sample> ...`) for the v1 reads of a label that v2
                                             does not count for that label: what v2 does with them (other label,
                                             decoy, ambiguous, low identity, no hit), from the v2 DIAMOND table.
Relative deviation = (v2 - v1) / v1 of copies_per_GE.
"""
import os, sys, gzip
import numpy as np
import pandas as pd

W = "/scratch1/afreiburger/salterns_2026/read_markers_v2"
V1 = "/scratch1/afreiburger/salterns_2026/read_markers"
os.chdir(W)


def main():
    a = pd.read_csv(f"{V1}/results/per_sample_KO.tsv", sep="\t")
    b = pd.read_csv("results/per_sample_KO.tsv", sep="\t")
    keep = ["sample", "label", "reads", "ambiguous_best", "low_identity_best", "copies_per_GE"]
    x = a[keep].merge(b[keep], on=["sample", "label"], suffixes=("_v1", "_v2"), how="inner")
    x["rel_dev"] = np.where(x.copies_per_GE_v1 > 0, (x.copies_per_GE_v2 - x.copies_per_GE_v1) / x.copies_per_GE_v1, np.nan)
    x.to_csv("results/regression_v1_v2_per_sample.tsv", sep="\t", index=False, float_format="%.6g")
    rows = []
    for lab, g in x.groupby("label"):
        g1 = g[g.reads_v1 > 0]; g20 = g[g.reads_v1 >= 20]
        rows.append(dict(label=lab, n_samples_v1_reads=len(g1), reads_v1=g.reads_v1.sum(), reads_v2=g.reads_v2.sum(),
                         pooled_ratio_reads=g.reads_v2.sum() / g.reads_v1.sum() if g.reads_v1.sum() else np.nan,
                         median_abs_rel_dev=g1.rel_dev.abs().median(), max_abs_rel_dev=g1.rel_dev.abs().max(),
                         median_abs_rel_dev_ge20=g20.rel_dev.abs().median(), max_abs_rel_dev_ge20=g20.rel_dev.abs().max(),
                         median_rel_dev=g1.rel_dev.median(), n_samples_ge20=len(g20),
                         ambiguous_v1=g.ambiguous_best_v1.sum(), ambiguous_v2=g.ambiguous_best_v2.sum(),
                         low_identity_v1=g.low_identity_best_v1.sum(), low_identity_v2=g.low_identity_best_v2.sum()))
    s = pd.DataFrame(rows).sort_values("max_abs_rel_dev_ge20", ascending=False)
    s.to_csv("results/regression_v1_v2_summary.tsv", sep="\t", index=False, float_format="%.4g")
    ga = pd.read_csv(f"{V1}/results/genome_equivalents.tsv", sep="\t", index_col=0)
    gb = pd.read_csv("results/genome_equivalents.tsv", sep="\t", index_col=0)
    gg = pd.DataFrame({"GE_depth_v1": ga.GE_depth_mean, "GE_depth_v2": gb.GE_depth_mean})
    gg["rel_dev"] = gg.GE_depth_v2 / gg.GE_depth_v1 - 1
    gg.to_csv("results/regression_genome_equivalents.tsv", sep="\t", float_format="%.6g")
    # mcrA genus shares, 12 unrestored cores pooled
    org = pd.read_csv("meta/kegg_org_gtdb.tsv", sep="\t", dtype=str).fillna("")
    kname = org.kegg_name.str.replace(r"^(Candidatus|uncultured) ", "", regex=True).str.split().str[0].str.strip("[]")
    nm = pd.DataFrame({"k": kname, "g": org.gtdb_genus})
    nm = nm[nm.g != ""].groupby("k").g.agg(lambda v: v.value_counts().index[0]).to_dict()   # majority R232 genus per KEGG name
    p1 = pd.read_csv(f"{V1}/results/mcrA_taxonomy_unrestored_pooled.tsv", sep="\t", dtype={"genus": str}).fillna({"genus": ""})
    p2 = pd.read_csv("results/mcrA_taxonomy_unrestored_pooled.tsv", sep="\t", dtype={"genus": str}).fillna({"genus": ""})
    p1["genus_h"] = [nm.get(g, g) for g in p1.genus]
    t1 = p1.groupby("genus_h").reads.sum(); t2 = p2.groupby("genus").reads.sum()
    t = pd.DataFrame({"reads_v1": t1, "reads_v2": t2}).fillna(0)
    t["share_v1"] = t.reads_v1 / t.reads_v1.sum(); t["share_v2"] = t.reads_v2 / t.reads_v2.sum()
    t.index.name = "genus (GTDB R232)"
    t.sort_values("reads_v2", ascending=False).to_csv("results/regression_mcrA_genus_pooled.tsv", sep="\t", float_format="%.4g")
    pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
    print(s.to_string())
    print(gg.describe().to_string())
    print(t.sort_values("reads_v2", ascending=False).head(15).to_string())


def attribute(samples):
    """For v1-counted reads of each label not counted for the same label in v2: the v2 decision."""
    sys.argv = ["x"]
    import importlib.util
    spec = importlib.util.spec_from_file_location("cl", "scripts/08_classify.py"); cl = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cl)
    cut = pd.read_csv("cutoffs.tsv", sep="\t", dtype=str); g = dict(zip(cut.param, cut.value))
    emax, lmin, margin, pdef = float(g["EMAX"]), int(g["LMIN"]), float(g["MARGIN"]), float(g["PMIN_DEFAULT"])
    pmin = {k[5:]: float(v) for k, v in g.items() if k.startswith("PMIN:")}
    for s in samples:
        r1 = pd.read_csv(f"{V1}/results/per_sample/{s}.reads.tsv.gz", sep="\t", dtype={"read": str, "mate": str})
        r2 = pd.read_csv(f"results/per_sample/{s}.reads.tsv.gz", sep="\t", dtype={"read": str, "mate": str})
        x = r1.merge(r2[["read", "mate", "label"]], on=["read", "mate"], how="left", suffixes=("_v1", "_v2"))
        lost = x[x.label_v2 != x.label_v1]
        want = {(r, m) for r, m in zip(lost.read, lost.mate)}
        dec = {}
        for m in ("1", "2"):
            for q, hs in cl.groups([f"dmnd/{s}_{m}.tsv.gz"]):
                if (q, m) in want:
                    L, b = cl.decide(hs, emax, lmin, margin, pmin, pdef)
                    dec[(q, m)] = (L or "no_hit_passing") if L not in ("ambiguous", "low_identity") else \
                        f"{L}(best={cl.LAB[b[4]]})"
        lost = lost.assign(v2_decision=[dec.get((r, m), "no_v2_hit") for r, m in zip(lost.read, lost.mate)])
        o = lost.groupby(["label_v1", "v2_decision"]).size().reset_index(name="reads")
        o.to_csv(f"results/regression_lost_reads_{s}.tsv", sep="\t", index=False)
        gained = x.label_v2.notna() & (x.label_v2 == x.label_v1)
        print(s, "v1 counted", len(r1), "same label in v2", int(gained.sum()), "changed", len(lost))
        print(o.sort_values("reads", ascending=False).head(30).to_string())


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "attribute":
        attribute(sys.argv[2:])
    else:
        main()
