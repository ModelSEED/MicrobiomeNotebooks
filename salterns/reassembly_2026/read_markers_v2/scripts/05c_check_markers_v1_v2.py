#!/usr/bin/env python
"""Step 5c (v2): are the marker-family references of v2 identical to v1's? Compares, per label present in both DBs,
the 99 % cluster representatives and member sets (refs/clu/clusters99.tsv) and the representative sequences.
If the 8 marker families are identical, GE-depth differences v1 -> v2 can only come from DIAMOND hits / ties against
the larger v2 DB. Writes results/check_labels_v1_v2.tsv."""
import os
import pandas as pd

W = "/scratch1/afreiburger/salterns_2026/read_markers_v2"
V1 = "/scratch1/afreiburger/salterns_2026/read_markers"
os.chdir(W)
cols = ["label", "rep", "member"]
a = pd.read_csv(f"{V1}/refs/clu/clusters99.tsv", sep="\t", header=None, names=cols, dtype=str)
b = pd.read_csv("refs/clu/clusters99.tsv", sep="\t", header=None, names=cols, dtype=str)
rows = []
for lab in sorted(set(a.label) & set(b.label)):
    x, y = a[a.label == lab], b[b.label == lab]
    mx, my = set(x.member), set(y.member)
    rx, ry = set(x.rep), set(y.rep)
    px = set(map(tuple, x[["rep", "member"]].values)); py = set(map(tuple, y[["rep", "member"]].values))
    rows.append(dict(label=lab, members_v1=len(mx), members_v2=len(my), members_only_v1=len(mx - my), members_only_v2=len(my - mx),
                     reps_v1=len(rx), reps_v2=len(ry), reps_only_v1=len(rx - ry), reps_only_v2=len(ry - rx),
                     identical_clustering=px == py))
o = pd.DataFrame(rows)
o.to_csv("results/check_labels_v1_v2.tsv", sep="\t", index=False)
pd.set_option("display.width", 200); pd.set_option("display.max_rows", 200)
print(o.to_string())
M = ["uL2", "uL11", "uS3", "uS11", "uS12", "uS9", "uL4", "uL16", "euk_marker"]
print("marker families identical:", bool(o[o.label.isin(M)].identical_clustering.all()))
print("labels only in v1:", sorted(set(a.label) - set(b.label)))
