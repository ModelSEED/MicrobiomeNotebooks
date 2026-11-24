#!/usr/bin/env python
"""Figure S11 for the final re-assembled MAG set (585 species representatives): within-site co-variation of prevalent
taxa (lower-triangle clustered heatmap). Port of ../r232/figureS11_r232.py (rendering of the 2020 cell
`within_condition_triangle_heatmap` in ../analysis_2_marimo.py); inputs are the outputs of rerun_correlational_final585.py,
no statistics are recomputed:
  C_within_final585.csv            Pearson r between within-site CLR residuals
  glasso_direct_edges_final585.csv Graphical Lasso edges; |partial r| > 0.05 are marked with a green asterisk
  guilds_final585.csv              KOfam guilds (any MAG of the species), collapsed to taxa; archaea marked with a dagger
The canvas grows with the number of taxa (34 in for up to 169 taxa, as the R232 figure).
Run from this folder: ~/Documents/py_venv/bin/python figureS11_final585.py -> figureS11_final585.png (+ .order.json)
"""
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import patches
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
SIZE, LABEL_FONT, HIGHLIGHT_FONT, DPI = 34.0, 12.0, 13.0, 200
STAR_COLOR, BOX_LW = "#90EE90", 2.0
LATIN = re.compile(r"^[A-Z][a-z]+$")          # genus names (not placeholder IDs or family names ending -aceae)

C = pd.read_csv(HERE / "C_within_final585.csv", index_col=0)
gl = pd.read_csv(HERE / "glasso_direct_edges_final585.csv")
g = pd.read_csv(HERE / "guilds_final585.csv")
g["taxon"] = g.iterativeID_final585.str.rsplit(".", n=1).str[0]
tx = lambda col: set(g.loc[g[col], "taxon"])
GROUPS = [  # (legend name, colour, taxa); first match wins
    ("Methanogen (mcrA)", "#8E44AD", tx("methanogen")),
    ("Betaine reducers (grdH or grdI)", "#1f77b4", tx("betaine_reducer")),
    ("MePn consumers (phnJ)", "#e8671b", tx("mepn_consumer")),
]
ARCHAEA = tx("archaea")
assert not C.isna().any().any(), "NaN in C_within - clustering would fail"
OVERLAP = {(a, b): sorted(tx(a) & tx(b)) for a, b in (("methanogen", "betaine_reducer"), ("methanogen", "mepn_consumer"), ("betaine_reducer", "mepn_consumer"))}
print("taxa in two guilds (first listed group wins the colour):", {k: v for k, v in OVERLAP.items() if v})
N = len(C)
SIZE = max(SIZE, SIZE * N / 169)


def group_of(t):
    return next(((n, c) for n, c, m in GROUPS if t in m), (None, None))


def shown(t):
    return f"{t} †" if t in ARCHAEA else t


cm = sns.clustermap(C, cmap="coolwarm_r", center=0, vmin=-1, vmax=1, figsize=(SIZE, SIZE),
                    dendrogram_ratio=(0.08, 0.08), xticklabels=True, yticklabels=True)
cm.ax_row_dendrogram.set_visible(False); cm.ax_col_dendrogram.set_visible(False)
ax = cm.ax_heatmap
ax.yaxis.set_ticks_position("left"); ax.yaxis.set_label_position("left")
for lbls, rot, setter in ((ax.get_xticklabels(), 90, ax.set_xticklabels), (ax.get_yticklabels(), 0, ax.set_yticklabels)):
    for lbl in lbls:
        t = lbl.get_text(); _, col = group_of(t)
        lbl.set_rotation(rot)
        lbl.set_fontstyle("italic" if LATIN.match(t) else "normal")
        if col:
            lbl.set_color(col); lbl.set_fontweight("bold"); lbl.set_fontsize(HIGHLIGHT_FONT)
        else:
            lbl.set_fontsize(LABEL_FONT)
        lbl.set_text(shown(t))
    setter(lbls)

ro, co = cm.dendrogram_row.reordered_ind, cm.dendrogram_col.reordered_ind
reo = C.iloc[ro, co]
assert list(reo.index) == list(reo.columns), "row and column orders differ"
mask = np.triu(np.ones_like(reo, dtype=bool), k=1)
mesh = ax.collections[0]
arr = mesh.get_array().reshape(reo.shape); arr[mask] = np.nan; mesh.set_array(arr.ravel())

direct = set()
for r in gl.itertuples():
    if abs(r.partial_r) > 0.05:
        direct |= {(r.taxon_a, r.taxon_b), (r.taxon_b, r.taxon_a)}
order = list(reo.index)
n_star = 0
for i in range(N):
    for j in range(i):
        if (order[i], order[j]) in direct:
            ax.text(j + 0.5, i + 0.5, "*", ha="center", va="center", color=STAR_COLOR,
                    fontsize=14, fontweight="bold", clip_on=False); n_star += 1
for i, t in enumerate(order):
    _, col = group_of(t)
    if col:
        ax.add_patch(patches.Rectangle((0, i), i + 1, 1, fill=False, edgecolor=col, linewidth=BOX_LW, clip_on=False))
        ax.add_patch(patches.Rectangle((i, i), 1, N - i, fill=False, edgecolor=col, linewidth=BOX_LW, clip_on=False))
ax.set_xlabel(""); ax.set_ylabel("")
ax.set_title("Within-site Pearson r of CLR residuals, final MAG set (585 species), GTDB R232 taxa  "
             r"(asterisks: Graphical Lasso direct association, $|r_{\rm partial}| > 0.05$)", fontsize=30, pad=20)
hm = ax.get_position()
cm.cax.set_position([hm.x0 + hm.width * 0.78, hm.y0 + hm.height * 0.62, hm.width * 0.025, hm.height * 0.30])
cm.cax.tick_params(labelsize=14); cm.cax.set_title("Pearson r", fontsize=18, pad=10)
present = [(n, c, m & set(C.index)) for n, c, m in GROUPS]
handles = [Line2D([], [], ls="", marker="s", ms=16, mfc="none", mec=c, mew=2.5, label=f"{n}: {len(m)} taxa")
           for n, c, m in present if m]
handles += [Line2D([], [], ls="", marker="$†$", ms=16, color="black", label=f"Archaeal taxon: {len(ARCHAEA & set(C.index))} taxa"),
            Line2D([], [], ls="", marker="$*$", ms=18, color=STAR_COLOR, mec="0.4", mew=0.3,
                   label=f"Graphical Lasso edge ({n_star} pairs)")]
ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.52, 0.55), fontsize=20, frameon=True,
          title="Labels, boxes and marks", title_fontsize=22)
cm.figure.savefig(HERE / "figureS11_final585.png", dpi=DPI, bbox_inches="tight")
json.dump({"order": order, "asterisks": n_star, "taxa": N,
           "groups": {n: sorted(m) for n, _, m in present}, "archaea": sorted(ARCHAEA & set(C.index))},
          open(HERE / "figureS11_final585.order.json", "w"), indent=1)
print(f"taxa={N} asterisks={n_star}", {n: len(m) for n, _, m in present}, "archaea", len(ARCHAEA & set(C.index)))
