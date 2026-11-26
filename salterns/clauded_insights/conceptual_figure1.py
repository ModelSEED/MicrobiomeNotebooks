#!/usr/bin/env python
"""
Proposed conceptual Figure 1: the study system and how the evidence converges.
Left  : the restoration / salinity-sulfate gradient across the four site types.
Middle: the two candidate methane pathways in a saltern sediment column.
Right : the five lines of evidence converging on the methylphosphonate conclusion.
Schematic (no data). Run: ~/Documents/py_venv/bin/python conceptual_figure1.py
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

MEPN, BET, ENV, INK = "#D55E00", "#0072B2", "#009E73", "#222222"
fig, ax = plt.subplots(figsize=(15, 8.2), dpi=200)
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")


def box(x, y, w, h, text, fc, ec=None, fs=10, tc="white", bold=True, alpha=1.0):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6,rounding_size=1.6",
                 fc=fc, ec=ec or fc, lw=1.4, alpha=alpha, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color=tc, fontweight="bold" if bold else "normal", zorder=3, wrap=True)


def arrow(x1, y1, x2, y2, color=INK, lw=2.2, style="-|>"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=16,
                 color=color, lw=lw, zorder=1))


ax.text(50, 97, "The hypersaline saltern methane system", ha="center", fontsize=17, fontweight="bold")

# ---------- LEFT: gradient ----------
ax.text(15, 90, "A  Restoration / salinity gradient", ha="center", fontsize=12, fontweight="bold", color=INK)
sites = [("Unrestored salterns (R1, R2)\nhypersaline · high SO$_4$ · HIGH CH$_4$", MEPN, 74),
         ("Restored saltern (SF2)", "#7aa8c4", 60),
         ("Reference wetland (R2A)\nlow salt · LOW CH$_4$", ENV, 46)]
for txt, c, yy in sites:
    box(3, yy, 24, 11, txt, c, fs=8.6, tc="white")
arrow(29, 50, 29, 84, color=INK, lw=3)
ax.text(31.5, 67, "salinity, sulfate,\nCH$_4$ increase", rotation=90, va="center",
        ha="center", fontsize=8.5, color=INK)
ax.text(15, 40, "3 cores × 2 depths per site\n(24 metagenomes; 12 depth-merged cores)",
        ha="center", fontsize=8.2, color="0.35", style="italic")

# ---------- MIDDLE: sediment + pathways ----------
ax.text(52, 90, "B  Two candidate CH$_4$ pathways", ha="center", fontsize=12, fontweight="bold", color=INK)
ax.add_patch(Rectangle((38, 60), 28, 14, fc="#eaf3f7", ec="0.6", zorder=0))       # oxic
ax.add_patch(Rectangle((38, 34), 28, 26, fc="#e9e2d0", ec="0.6", zorder=0))       # anoxic
ax.text(64.5, 72.5, "oxic", ha="right", fontsize=8, color="0.4", style="italic")
ax.text(64.5, 36, "anoxic", ha="right", fontsize=8, color="0.4", style="italic")
box(40, 63.5, 24, 8, "Methylphosphonate  --C–P lyase-->  CH$_4$ + P$_i$\n(Roseovarius, Halomonas · P scavenging)",
    MEPN, fs=8.0)
box(40, 45, 24, 9, "Glycine betaine --> trimethylamine (TMA)\n--> CH$_4$   (methylotrophic methanogen, archaeal)",
    BET, fs=8.0)
ax.text(52, 39.5, "TMA also from betaine reduction by\nHalanaerobiales / Clostridia (bacterial)",
        ha="center", fontsize=7.6, color="0.35", style="italic")
arrow(52, 60, 52, 54.5, color="0.4", lw=1.6, style="-|>")

# ---------- RIGHT: evidence ----------
ax.text(86, 90, "C  Converging lines of evidence", ha="center", fontsize=12, fontweight="bold", color=INK)
ev = [("Substrate-addition incubations (Fig 4)\nMePn-acid --> more CH$_4$", MEPN, 79, "causal"),
      ("Halomonas culturing (Fig 5-6)\nC–P lyase strain makes CH$_4$ on MePn", MEPN, 70, "causal"),
      ("Metabolomics / NMR (Fig 2-3)\nTMA rho=0.86, betaine 0.70 vs CH$_4$", BET, 61, "correlative"),
      ("Metagenomics / MAGs (Fig 1, 8)\nphn genes; guild networks", "0.5", 52, "descriptive"),
      ("Metabolic modelling (Fig 7)\nMePn CAN yield the CH$_4$; Roseovarius 77%", MEPN, 43, "sufficiency")]
for txt, c, yy, tag in ev:
    box(73, yy, 25, 7, txt, c, fs=7.4, tc="white")
    ax.text(99.4, yy + 3.5, tag, ha="right", va="center", fontsize=6.6, color="0.45",
            rotation=0, style="italic")
    arrow(72.8, yy + 3.5, 69, 28.5, color="0.7", lw=1.1)
box(40, 22, 40, 9,
    "Methylphosphonate degradation is the leading\nCH$_4$ source in the unrestored salterns",
    INK, fs=10.5)
ax.text(50, 16.5, "strongest support = the manipulative experiments (incubations, culturing);\n"
        "cross-sample correlations track the salinity/sulfate gradient and are not organism-specific",
        ha="center", fontsize=8.0, color="0.35", style="italic")

plt.tight_layout()
fig.savefig("conceptual_figure1.png", bbox_inches="tight")
fig.savefig("conceptual_figure1.svg", bbox_inches="tight")
print("Saved conceptual_figure1.png/.svg")
