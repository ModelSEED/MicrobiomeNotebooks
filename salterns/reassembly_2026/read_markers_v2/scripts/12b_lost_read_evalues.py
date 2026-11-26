#!/usr/bin/env python
"""Step 12b (v2): for v1-counted reads of selected labels that v2 does not count for any label ('no_hit_passing'),
the E-value / aligned length of the same read's best hit to its v1 representative in v1 and in v2 (both DBs contain that
representative unchanged). Tests whether the losses are E-value inflation from the ~2.5x larger v2 DB.
Usage: 12b_lost_read_evalues.py <sample> <label> [<label> ...] -> results/regression_lost_read_evalues_<sample>.tsv"""
import sys, gzip, os
import numpy as np
import pandas as pd
W = "/scratch1/afreiburger/salterns_2026/read_markers_v2"; V1 = "/scratch1/afreiburger/salterns_2026/read_markers"
os.chdir(W)
s, labs = sys.argv[1], set(sys.argv[2:])
r1 = pd.read_csv(f"{V1}/results/per_sample/{s}.reads.tsv.gz", sep="\t", dtype={"read": str, "mate": str})
r2 = pd.read_csv(f"results/per_sample/{s}.reads.tsv.gz", sep="\t", dtype={"read": str, "mate": str})
r1 = r1[r1.label.isin(labs)]
x = r1.merge(r2[["read", "mate", "label"]], on=["read", "mate"], how="left", suffixes=("", "_v2"))
lost = x[x.label_v2.isna()]
want = {(r, m): rep for r, m, rep in zip(lost.read, lost.mate, lost.rep)}
res = {}
for tag, base in (("v1", V1), ("v2", W)):
    for m in ("1", "2"):
        with gzip.open(f"{base}/dmnd/{s}_{m}.tsv.gz", "rt") as fh:
            for l in fh:
                f = l.split("\t", 12)
                k = (f[0], m)
                if k in want and f[1] == want[k]:
                    res.setdefault(k, {})[tag] = (float(f[10]), int(f[3]), float(f[11]))
rows = []
for k, rep in want.items():
    d = res.get(k, {})
    v1 = d.get("v1", (np.nan, np.nan, np.nan)); v2 = d.get("v2", (np.nan, np.nan, np.nan))
    rows.append(dict(read=k[0], mate=k[1], rep=rep, label=lost.set_index(["read", "mate"]).loc[k, "label"],
                     evalue_v1=v1[0], alen_v1=v1[1], bits_v1=v1[2], evalue_v2=v2[0], alen_v2=v2[1], bits_v2=v2[2]))
o = pd.DataFrame(rows)
o["ratio_v2_v1"] = o.evalue_v2 / o.evalue_v1
o.to_csv(f"results/regression_lost_read_evalues_{s}.tsv", sep="\t", index=False, float_format="%.3g")
print(s, "lost reads", len(o), "same rep found in v2 table:", o.evalue_v2.notna().sum())
print("v1 E of the lost reads: median %.2g, 90th pct %.2g" % (o.evalue_v1.median(), o.evalue_v1.quantile(0.9)))
print("v2 E (same rep):          median %.2g, min %.2g" % (o.evalue_v2.median(), o.evalue_v2.min()))
print("E ratio v2/v1: median %.2f" % o.ratio_v2_v1.median(), "| v2 E in (1e-5, 1e-4]:", int(((o.evalue_v2 > 1e-5) & (o.evalue_v2 <= 1e-4)).sum()),
      "| v1 E in (1e-6, 1e-5]:", int(((o.evalue_v1 > 1e-6) & (o.evalue_v1 <= 1e-5)).sum()), "| same rep missing in v2 (not in top 5):", int(o.evalue_v2.isna().sum()))
print(o.groupby("label").size().to_dict())
