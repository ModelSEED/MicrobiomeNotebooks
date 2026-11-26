#!/usr/bin/env python
"""Step 8c: summarise the validation grids (val/grid_<mode>.tsv) per parameter set:
macro-averaged recall and precision over the target labels and over the 8 markers (labels with >= 50 held-out reads),
min precision over targets, background false positives. Writes val/grid_summary.tsv and prints the top rows."""
import os, sys
import pandas as pd

os.chdir("/scratch1/afreiburger/salterns_2026/read_markers")
rows = []
for mode in ["default", "mid-sensitive", "sensitive"]:
    p = f"val/grid_{mode}.tsv"
    if not os.path.exists(p):
        continue
    g = pd.read_csv(p, sep="\t")
    for key, gg in g.groupby(["EMAX", "LMIN", "MARGIN", "PMIN"]):
        t = gg[(gg.role == "target") & (gg.n_true >= 50)]
        m = gg[(gg.role == "marker") & (gg.n_true >= 50)]
        rows.append(dict(mode=mode, EMAX=key[0], LMIN=key[1], MARGIN=key[2], PMIN=key[3],
                         target_recall=t.recall.mean(), target_precision=t.precision.mean(),
                         target_min_precision=t.precision.min(), marker_recall=m.recall.mean(),
                         marker_precision=m.precision.mean(), fp_background=gg.fp_background.sum(),
                         fp_holdout=gg.fp_holdout.sum(), uh_assigned=gg.uh_assigned.sum(),
                         f1_targets=(2 * t.recall * t.precision / (t.recall + t.precision)).mean()))
s = pd.DataFrame(rows)
s.to_csv("val/grid_summary.tsv", sep="\t", index=False, float_format="%.4f")
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
print(s.sort_values("f1_targets", ascending=False).head(40).to_string())
