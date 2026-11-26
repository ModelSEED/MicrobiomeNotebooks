#!/usr/bin/env python
"""Figure S10 for the 2026 guild models: core (ATP) gap-filling per ATP test medium of each arm-A draft
(scan/*.json from ../pipeline/atp_media_scan.py, run on TF), for the same 22 media as the original Figure S10, as
two clustered heatmaps: (a) the MePn-consumer models, (b) the betaine-reducer models. The aerobic L-lactate medium
(LLac.O2) used for the simulations is boxed in red. Cells without a gap-filling solution are grey ('–').
Run: ~/Documents/py_venv/bin/python figureS10_final.py -> figureS10_final.png/.svg, figureS10_counts.tsv, figureS10_summary.json
"""
import glob
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patches
from scipy.cluster.hierarchy import leaves_list, linkage

HERE = Path(__file__).resolve().parent
S = HERE.parents[2]
MEDIA = ["Glc", "Ac.O2", "Pyr.O2", "Glc.O2", "Glyc.O2", "Dlac.O2", "LLac.O2", "Fum.O2", "Fum", "Pyr", "For.O2",
         "Etho.O2", "Akg", "Akg.O2", "Succ.O2", "Glyc", "Llac", "Succ", "Dlac", "Methanol", "Ac", "For"]  # original S10
SIM = "LLac.O2"
IDS = json.load(open(S / "correlational/final585/iterativeIDs_final585.json"))
rows = {}
for f in sorted(glob.glob(str(HERE / "scan" / "*.json"))):
    r = json.load(open(f))
    rows[(r["guild"], r["mag_id"])] = r["core_gapfill_by_medium"]
df = pd.DataFrame(rows).T[MEDIA].astype(float)
df.index.names = ["guild", "mag_id"]
out = df.copy()
out.insert(0, "label", [IDS[m] for _, m in df.index])
out.to_csv(HERE / "figureS10_counts.tsv", sep="\t")


def order(m):
    filled = m.fillna(np.nanmax(m.values) + 2)
    r = leaves_list(linkage(filled.values, "average")) if len(m) > 2 else np.arange(len(m))
    c = leaves_list(linkage(filled.values.T, "average"))
    return m.iloc[r, c]


summ = {}
panels = [("mepn_consumer", "a", "Core gap filling – MePn consumers"), ("betaine_reducer", "b", "Core gap filling – betaine reducers")]
heights = [len(df.loc[g]) for g, _, _ in panels]
fig, axes = plt.subplots(2, 1, figsize=(11, 0.26 * sum(heights) + 4.5), dpi=200,
                         gridspec_kw={"height_ratios": heights, "hspace": 0.32})
vmax = np.nanmax(df.values)
for ax, (g, letter, title) in zip(axes, panels):
    m = df.loc[g].copy()
    m.index = [IDS[i] for i in m.index]
    m = order(m)
    arr = np.ma.masked_invalid(m.values)
    cmap = plt.cm.viridis.copy(); cmap.set_bad("#bdbdbd")
    im = ax.imshow(arr, cmap=cmap, vmin=0, vmax=vmax, aspect="auto")
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            v = m.values[i, j]
            ax.text(j, i, "–" if np.isnan(v) else f"{int(v)}", ha="center", va="center", fontsize=6.5,
                    color="white" if (np.isnan(v) or v < vmax * 0.6) else "black")
    ax.set_xticks(range(m.shape[1])); ax.set_xticklabels(m.columns, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(m.shape[0])); ax.set_yticklabels(m.index, fontsize=7)
    j = list(m.columns).index(SIM)
    ax.add_patch(patches.Rectangle((j - 0.5, -0.5), 1, m.shape[0], fill=False, ec="red", lw=2))
    ax.set_title(title, fontsize=11)
    ax.text(-0.02, 1.02, letter, transform=ax.transAxes, fontsize=16, fontweight="bold", ha="right", va="bottom")
    med = m.median()
    rank = int((med < med[SIM]).sum()) + 1
    summ[g] = {"n_models": int(len(m)), "median_by_medium": med.round(2).to_dict(), "LLac.O2_median": float(med[SIM]),
               "LLac.O2_rank_by_median": rank, "n_media": len(MEDIA),
               "media_with_lower_median": sorted(med[med < med[SIM]].index),
               "no_solution_cells": int(m.isna().sum().sum())}
cb = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.01)
cb.set_label("reactions added or reversed by core gap-filling", fontsize=9)
fig.savefig(HERE / "figureS10_final.png", bbox_inches="tight")
fig.savefig(HERE / "figureS10_final.svg", bbox_inches="tight")
json.dump(summ, open(HERE / "figureS10_summary.json", "w"), indent=1)
for g, s in summ.items():
    print(g, "n", s["n_models"], "LLac.O2 median", s["LLac.O2_median"], "rank", s["LLac.O2_rank_by_median"], "of",
          s["n_media"], "lower:", s["media_with_lower_median"], "no-solution cells", s["no_solution_cells"])
