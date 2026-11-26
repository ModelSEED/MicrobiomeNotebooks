#!/usr/bin/env python
"""Figure 10: the hypothesized methane sources of the unrestored salterns, one row per pathway:
(a) the pathway (substrate -> key enzyme -> product; TMA from betaine and from choline/carnitine feeds methylotrophic
    methanogenesis); (b) its abundance per 100 genome equivalents in all reads (open) and in the MAGs (filled), median
    and range over the 12 samples of the six unrestored cores; (c) its most abundant organisms (reads: genus of the best
    hit, share of the pathway's reads; MAGs: the most abundant guild species); (d) the turnover each pathway would need to
    supply the measured CH4 of the six cores on its own: substrate pool turnover (diamonds; measured porewater pools) and
    the doubling time of its organisms (squares; per-gram yields x measured biomass, MAG-based filled, read-based open).
Colours: the dataviz reference palette, slots 1-3 only (blue = methylphosphonate, orange = TMA producers, aqua =
methanogenesis); identity is carried by the row labels, never by colour alone.
Inputs: pathway_summary.json, pathway_abundance.tsv, pathway_turnover.tsv (pathway_summary.py).
Run: ~/Documents/py_venv/bin/python figure10_pathways.py -> figure10_pathways.png / .svg
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = Path(__file__).resolve().parent
INK, INK2, GRID = "#1f1f1f", "#5a5a5a", "#e6e6e6"
COL = {"mpn": "#2a78d6", "tma": "#eb6834", "meth": "#1baf7a"}
TINT = {"mpn": "#e3eefb", "tma": "#fde8df", "meth": "#ddf4eb"}
CORES = ["R1_A", "R1_B", "R1_C", "R2_A", "R2_B", "R2_C"]

J = json.load(open(HERE / "pathway_summary.json"))
A = pd.read_csv(HERE / "pathway_abundance.tsv", sep="\t", index_col=0)
T = pd.read_csv(HERE / "pathway_turnover.tsv", sep="\t", index_col=0)
PW = J["pathways"]
n = len(PW)
SUB = {"mpn": "Methylphosphonate", "betaine": "Glycine betaine", "choline": "Choline, carnitine",
       "methylamines": "TMA, DMA, MMA", "methanol": "Methanol", "hydrogenotrophic": "H$_2$/CO$_2$, formate",
       "acetoclastic": "Acetate"}
ENZ = {"mpn": "C–P lyase (PhnJ)", "betaine": "betaine reductase (GrdI)", "choline": "CutC, CntA",
       "methylamines": "MttB, MtbB, MtmB + Mcr", "methanol": "MtaB + Mcr", "hydrogenotrophic": "Mcr (hydrogenotrophic)",
       "acetoclastic": "Mcr (acetoclastic)"}
PROD = {"mpn": "CH$_4$ + P$_i$", "betaine": "TMA", "choline": "TMA", "methylamines": "CH$_4$", "methanol": "CH$_4$",
        "hydrogenotrophic": "CH$_4$", "acetoclastic": "CH$_4$"}
WHERE = {"mpn": "aerobic bacteria, haloarchaea", "betaine": "anaerobic Bacillota", "choline": "bacteria",
         "methylamines": "methylotrophic methanogens", "methanol": "methylotrophic methanogens",
         "hydrogenotrophic": "hydrogenotrophic methanogens", "acetoclastic": "acetoclastic methanogens"}


def it(g):
    """Italicise Latin genus names (not GTDB placeholders such as M55B127)."""
    return f"$\\it{{{g}}}$" if g[:1].isupper() and g[1:].islower() and g.isalpha() else g


def mag_label(lbl):
    base, num = lbl.rsplit(".", 1)
    return f"{it(base)}.{num}"


fig = plt.figure(figsize=(14.5, 8.2), dpi=300)
gs = fig.add_gridspec(1, 4, width_ratios=[1.4, 0.95, 1.45, 1.3], wspace=0.08, left=0.01, right=0.99, top=0.88,
                      bottom=0.13)
axA, axB, axC, axD = (fig.add_subplot(gs[0, i]) for i in range(4))
ys = np.arange(n)[::-1].astype(float)                       # row 0 at the top
for ax in (axA, axB, axC, axD):
    ax.set_ylim(-0.6, n - 0.4)
    for y in ys:
        ax.axhspan(y - 0.5, y + 0.5, color="#f6f6f5" if int(y) % 2 == 0 else "white", zorder=0, lw=0)

# ---- (a) pathways -------------------------------------------------------------------------------------------------------
axA.set_xlim(0, 10); axA.axis("off")
ypos = {p["key"]: y for p, y in zip(PW, ys)}
for p, y in zip(PW, ys):
    k, c = p["key"], p["cls"]
    for x0, w, txt in ((0.15, 3.3, SUB[k]), (7.55, 2.3, PROD[k])):
        axA.add_patch(FancyBboxPatch((x0, y - 0.2), w, 0.4, boxstyle="round,pad=0.02,rounding_size=0.08",
                                     fc=TINT[c], ec=COL[c], lw=1.1, zorder=2))
        axA.text(x0 + w / 2, y, txt, ha="center", va="center", fontsize=8.1, color=INK, zorder=3)
    axA.add_patch(FancyArrowPatch((3.5, y), (7.5, y), arrowstyle="-|>", mutation_scale=10, color=COL[c], lw=1.4,
                                  zorder=2))
    axA.text(5.5, y + 0.07, ENZ[k], ha="center", va="bottom", fontsize=7.4, color=INK)
    axA.text(5.5, y - 0.07, WHERE[k], ha="center", va="top", fontsize=6.4, color=INK2, style="italic")
ym = ypos["methylamines"] + 0.4                                  # TMA feeds the methylamine row (right-margin route)
for src in ("betaine", "choline"):
    axA.plot([9.85, 9.97, 9.97], [ypos[src], ypos[src], ym], color=COL["tma"], lw=0.9, ls=(0, (3, 2)), zorder=1,
             clip_on=False)
axA.plot([9.97, 1.8], [ym, ym], color=COL["tma"], lw=0.9, ls=(0, (3, 2)), zorder=1, clip_on=False)
axA.annotate("", xy=(1.8, ypos["methylamines"] + 0.21), xytext=(1.8, ym),
             arrowprops=dict(arrowstyle="-|>", color=COL["tma"], lw=0.9, mutation_scale=8))
axA.text(5.9, ym + 0.03, "TMA", fontsize=7, color=COL["tma"], ha="center", va="bottom")
axA.set_title("a  Pathway", loc="left", fontsize=11, color=INK)

# ---- (b) abundance ------------------------------------------------------------------------------------------------------
lo, hi = 3e-4, 40
for p, y in zip(PW, ys):
    k, c = p["key"], COL[p["cls"]]
    for col, dy, filled in ((f"{k}_reads", 0.13, False), (f"{k}_MAGs", -0.13, True)):
        v = A[col]
        if (v <= 0).all():
            axB.text(lo * 1.4, y + dy, "no MAG", fontsize=7.2, color=INK2, va="center", ha="left")
            continue
        v = v[v > 0]
        axB.hlines(y + dy, v.min(), v.max(), color=c, lw=1.6, zorder=2)
        axB.scatter([v.median()], [y + dy], s=42, marker="o", facecolors=c if filled else "white", edgecolors=c,
                    linewidths=1.3, zorder=3)
    nreads = J["mcrA_groups"]["reads"].get(k) if p["cls"] == "meth" and k not in ("methanol", "methylamines") else None
    if nreads is not None:
        axB.text(hi * 0.85, y + 0.13, f"n = {nreads} mcrA reads", fontsize=6.6, color=INK2, ha="right", va="center")
axB.set_xscale("log"); axB.set_xlim(lo, hi); axB.set_yticks([])
axB.set_xlabel("Copies (reads) or genomes (MAGs)\nper 100 genome equivalents", fontsize=8.6, color=INK)
axB.grid(True, axis="x", which="major", color=GRID, lw=0.7); axB.set_axisbelow(True)
axB.tick_params(labelsize=8)
for sp in ("top", "right", "left"):
    axB.spines[sp].set_visible(False)
axB.set_title("b  Abundance", loc="left", fontsize=11, color=INK)

# ---- (c) organisms ------------------------------------------------------------------------------------------------------
axC.set_xlim(0, 1); axC.axis("off")
for p, y in zip(PW, ys):
    top = p["top_reads"] if sum(len(g) for g, _ in p["top_reads"]) < 42 else p["top_reads"][:2]
    tr = ", ".join(f"{it(g)} {s * 100:.0f}%" for g, s in top) or "—"
    tm = ", ".join(mag_label(m) for m in p["top_mags"]) or "none"
    axC.text(0.03, y + 0.13, "Reads: " + tr, fontsize=7.6, color=INK, va="center")
    axC.text(0.03, y - 0.15, "MAGs: " + tm, fontsize=7.6, color=INK, va="center")
axC.set_title("c  Most abundant organisms", loc="left", fontsize=11, color=INK)

# ---- (d) turnover -------------------------------------------------------------------------------------------------------
tlo, thi = 1 / 300, 24 * 365 * 3
for p, y in zip(PW, ys):
    k, c = p["key"], COL[p["cls"]]
    pool = T.loc[k, [f"pool_turnover_d_{cc}" for cc in CORES]].astype(float) * 24
    if pool.notna().all():
        axD.hlines(y + 0.16, pool.min(), pool.max(), color=c, lw=1.4, zorder=2)
        axD.scatter([pool.median()], [y + 0.16], s=40, marker="D", color=c, edgecolor="white", linewidths=0.7, zorder=3)
    else:
        axD.text(tlo * 1.5, y + 0.16, "pool not quantified" if k == "mpn" else "pool not measured", fontsize=7,
                 color=INK2, va="center",
                 bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.9))
    rows = [(kind, off) for kind, off in (("MAG", -0.05), ("reads", -0.2), ("reads_lower_bound", -0.2))
            if f"doubling_h_{kind}_{CORES[0]}" in T.columns and pd.notna(T.loc[k, f"doubling_h_{kind}_{CORES[0]}"])]
    if not rows:
        axD.text(tlo * 1.5, y - 0.12, "no organism model" if k != "methanol" else "same methanogens as the row above",
                 fontsize=7, color=INK2, va="center", bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.9))
    for kind, off in rows:
        d = T.loc[k, [f"doubling_h_{kind}_{cc}" for cc in CORES]].astype(float)
        axD.hlines(y + off, d.min(), d.max(), color=c, lw=1.4, zorder=2, ls="-" if kind != "reads_lower_bound" else ":")
        filled = kind == "MAG"
        axD.scatter([d.median()], [y + off], s=36, marker="s", facecolors=c if filled else "white", edgecolors=c,
                    linewidths=1.2, zorder=3)
        if kind == "reads_lower_bound":
            axD.annotate("", xy=(d.max() * 6, y + off), xytext=(d.max() * 1.3, y + off),
                         arrowprops=dict(arrowstyle="-|>", color=c, lw=0.9, mutation_scale=7))
    if k == "mpn":                                               # MPn as the only P source (read-based)
        d = T.loc[k, [f"doubling_h_reads_onlyP_{cc}" for cc in CORES]].astype(float)
        axD.hlines(y - 0.36, d.min(), d.max(), color=c, lw=1.0, zorder=2)
        axD.scatter([d.median()], [y - 0.36], s=26, marker="s", facecolors="white", edgecolors=c, linewidths=1.0, zorder=3)
        axD.text(d.max() * 1.4, y - 0.36, "MPn the only P source", fontsize=6.6, color=INK2, va="center",
                 bbox=dict(fc="white", ec="none", pad=0.4, alpha=0.9))
TT = ((1 / 60, "1 min"), (1, "1 h"), (24, "1 d"), (24 * 7, "1 wk"), (24 * 365, "1 yr"))
for h_, _ in TT:
    axD.axvline(h_, color=GRID, lw=0.8, zorder=1)
axD.set_xscale("log"); axD.set_xlim(tlo, thi); axD.set_yticks([])
axD.set_xticks([h_ for h_, _ in TT]); axD.set_xticklabels([l_ for _, l_ in TT]); axD.minorticks_off()
axD.set_xlabel("Turnover needed to supply the measured CH$_4$", fontsize=8.6, color=INK)
axD.tick_params(labelsize=8)
for sp in ("top", "right", "left"):
    axD.spines[sp].set_visible(False)
axD.set_title("d  Turnover needed", loc="left", fontsize=11, color=INK)

# ---- legend -------------------------------------------------------------------------------------------------------------
h = [Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor=INK2, markersize=6),
     Line2D([], [], marker="o", ls="", color=INK2, markersize=6),
     Line2D([], [], marker="D", ls="", color=INK2, markersize=5),
     Line2D([], [], marker="s", ls="", color=INK2, markersize=5),
     Line2D([], [], marker="s", ls="", markerfacecolor="white", markeredgecolor=INK2, markersize=5),
     Line2D([], [], ls=":", color=INK2)]
l = ["b: all reads (median, range of 12 samples)", "b: MAGs", "d: substrate pool turnover (6 cores)",
     "d: doubling time, MAG biomass", "d: doubling time, read-based biomass",
     "d: lower bound (no model; methylotroph yield)"]
fig.legend(h, l, loc="lower center", ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(0.62, 0.0))
fig.savefig(HERE / "figure10_pathways.png", bbox_inches="tight")
fig.savefig(HERE / "figure10_pathways.svg", bbox_inches="tight")
print("wrote figure10_pathways.png")
