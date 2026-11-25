#!/usr/bin/env python
"""Step 9c: reduced-reference re-mapping (09) vs the full self-mapping BAM of the binning test (09b), per classified
read (>= 95 % identity in both): same placement status (kept_mag / other_bin / unbinned / unassembled) and same contig.
Writes results/fullbam_check.tsv (per sample x label) and prints the overall agreement."""
import os, sys
import pandas as pd

os.chdir("/scratch1/afreiburger/salterns_2026/read_markers")
st = pd.read_csv("meta/contig_status.tsv.gz", sep="\t", dtype=str)
rows = []
for s in sys.argv[1:]:
    rd = pd.read_csv(f"results/per_sample/{s}.reads.tsv.gz", sep="\t", dtype={"read": str, "mate": str})[["read", "mate", "label"]]
    cols = ["read", "mate", "contig", "pos", "mapq", "ident"]
    a = pd.read_csv(f"map/{s}.hits.tsv.gz", sep="\t", header=None, names=cols, dtype={"read": str, "mate": str})
    b = pd.read_csv(f"map/{s}.fullbam_hits.tsv.gz", sep="\t", header=None, names=cols, dtype={"read": str, "mate": str})
    a, b = a[a.ident >= 0.95], b[b.ident >= 0.95]
    sts = st[st["sample"] == s].set_index("contig").status
    x = rd.merge(a[["read", "mate", "contig"]], on=["read", "mate"], how="left") \
          .merge(b[["read", "mate", "contig"]], on=["read", "mate"], how="left", suffixes=("_reduced", "_full"))
    x["status_reduced"] = x.contig_reduced.map(sts).fillna("unassembled")
    x["status_full"] = x.contig_full.map(sts).fillna("unassembled")
    for lab, g in x.groupby("label"):
        rows.append(dict(sample=s, label=lab, reads=len(g), same_status=(g.status_reduced == g.status_full).mean(),
                         same_contig_when_both_mapped=(g.contig_reduced == g.contig_full)[g.contig_reduced.notna() & g.contig_full.notna()].mean(),
                         unassembled_reduced=(g.status_reduced == "unassembled").mean(),
                         unassembled_full=(g.status_full == "unassembled").mean(),
                         kept_mag_reduced=(g.status_reduced == "kept_mag").mean(), kept_mag_full=(g.status_full == "kept_mag").mean()))
    print(s, "reads", len(x), "same status", round((x.status_reduced == x.status_full).mean(), 4),
          "unassembled reduced/full", round((x.status_reduced == "unassembled").mean(), 4), round((x.status_full == "unassembled").mean(), 4))
pd.DataFrame(rows).to_csv("results/fullbam_check.tsv", sep="\t", index=False, float_format="%.4f")
