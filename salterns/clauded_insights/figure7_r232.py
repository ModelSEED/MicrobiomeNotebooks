#!/usr/bin/env python
"""
Figure 7 for the R232 revision: simulated versus measured methane in the six unrestored cores, on the dimensionally
consistent areal basis of modelling_consistent_basis_r232.py (run that first).

- Methylphosphonate: the 20 simulated MePn-consumer models at phosphate source ratio 0.07, each MAG's flux per unit of
  simulated abundance rescaled to its biomass in the core (relative abundance x DNA x layer thickness, summed over
  D1 and D2). The earlier 'areal' version (figure7_corrected.py panel a) kept a 1/DNA term from the model input and is
  superseded.
- Betaine: the nine grdH carriers (methane = 0.06 x TMA), converted with the core's volume-weighted DNA, because
  their TMA excretion is mostly set by the uptake bound rather than by biomass.
Panel (b) of the submitted figure (abundance vs methane) is superseded by Figure 9, Figure S8 and Table 5.

Run: ~/Documents/py_venv/bin/python figure7_r232.py   -> figure7_r232.png / .svg / _stats.json
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
from matplotlib.ticker import FuncFormatter, LogFormatterMathtext, LogLocator

HERE = Path(__file__).resolve().parent
C_MEPN, C_BET = "#D55E00", "#0072B2"      # Okabe-Ito, as in figure7_corrected.py
CORES = ["R1_A", "R1_B", "R1_C", "R2_A", "R2_B", "R2_C"]

d = json.load(open(HERE / "modelling_consistent_basis_r232.json"))
x = np.array([d["measured_CH4"][c] for c in CORES])
series = [(d["areal"]["MePn_0.07_20"], C_MEPN, "o", "Methylphosphonate consumers (20 models)", "Methylphosphonate"),
          (d["areal"]["betaine_grdH9"], C_BET, "s", "Betaine reducers (9 grdH carriers)", "Betaine")]

fig, ax = plt.subplots(figsize=(7.2, 6.4), dpi=300)
texts, ys_all = [], []
for s, c, mk, lab, _ in series:
    y = np.array([s["sim_by_core"][k] for k in CORES]); ys_all.append(y)
    ax.scatter(x, y, s=62 if mk == "s" else 48, marker=mk, color=c, edgecolor="white", lw=1.0,
               zorder=4 if mk == "o" else 3, label=lab)
    xs = np.linspace(x.min(), x.max(), 400)
    yl = s["slope"] * xs + s["intercept"]
    ax.plot(xs[yl > 0], yl[yl > 0], "--", color=c, lw=1.6, zorder=2)
    texts += [ax.text(xi, yi, k.replace("_", " "), fontsize=7.5, color="0.25") for xi, yi, k in zip(x, y, CORES)]
ax.set_xscale("log"); ax.set_yscale("log")
for axis in (ax.xaxis, ax.yaxis):
    axis.set_major_locator(LogLocator(base=10, subs=(1.0,)))
    axis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
    axis.set_major_formatter(LogFormatterMathtext(base=10, labelOnlyBase=True))
    axis.set_minor_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}" if round(v) in (200, 300, 500, 2000, 3000, 5000) else ""))
    axis.set_tick_params(which="minor", labelsize=8)
ax.grid(True, which="both", linestyle="-", color="0.9", lw=0.6, zorder=0)
ax.set_axisbelow(True)
lo = min(x.min(), *(y.min() for y in ys_all)) * 0.5
hi = max(x.max(), *(y.max() for y in ys_all)) * 2.0
ax.plot([lo, hi], [lo, hi], ":", color="0.35", lw=1.3, zorder=1, label="1:1")
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
adjust_text(texts, ax=ax, expand=(1.3, 1.6), arrowprops=dict(arrowstyle="-", color="0.6", lw=0.5))
ax.set_xlabel(r"Measured CH$_4$ emission (µmol m$^{-2}$ d$^{-1}$)", fontsize=11)
ax.set_ylabel(r"Simulated CH$_4$ emission (µmol m$^{-2}$ d$^{-1}$)", fontsize=11)
for (s, c, _, _, name), y0 in zip(series, (0.97, 0.80)):
    txt = (f"{name}: y = {s['slope']:.2f}x {'+' if s['intercept'] >= 0 else '−'} {abs(s['intercept']):.0f}\n"
           f"$R^2$ = {s['pearson_r2']:.2f};  Spearman ρ = {s['spearman_rho']:.2f} (n = 6)").replace("-", "−")
    ax.text(0.03, y0, txt, transform=ax.transAxes, fontsize=8.5, color="0.15", va="top",
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=c, lw=1.2))
ax.legend(loc="lower right", fontsize=8.5, framealpha=0.95)
fig.tight_layout()
fig.savefig(HERE / "figure7_r232.png", bbox_inches="tight")
fig.savefig(HERE / "figure7_r232.svg", bbox_inches="tight")
json.dump({"MePn_0.07_20models": series[0][0], "betaine_grdH9": series[1][0], "measured_CH4": d["measured_CH4"],
           "basis": d["basis"]}, open(HERE / "figure7_r232_stats.json", "w"), indent=1)
print("MePn", {k: round(v, 3) for k, v in series[0][0].items() if isinstance(v, float)})
print("betaine", {k: round(v, 3) for k, v in series[1][0].items() if isinstance(v, float)})
