#!/usr/bin/env python
"""Figure 7 for the final MAG set: simulated versus measured methane in the six unrestored cores, from the 2026 guild
models (one per MePn-consumer and betaine-reducer species; ../models/), on the consistent areal basis of
consistent_basis_final.py (run that first).

- Methylphosphonate: the 38 MePn-consumer models at phosphate source ratio 0.07 (filled circles, with fit), and at 0.33
  (open circles; same pattern, 3.8-fold higher). Each model's CH4 per unit of simulated biomass (reference run) times
  its biomass in the core (relative abundance x DNA x layer thickness, summed over D1 and D2).
- Betaine: the 18 betaine-reducer models (methane = 0.06 x TMA), converted with the core's volume-weighted DNA,
  because their TMA excretion is mostly set by the 0.1 betaine-uptake bound rather than by biomass.

Run: ~/Documents/py_venv/bin/python figure7_final.py   -> figure7_final.png / .svg / _stats.json
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
C_MEPN, C_BET = "#D55E00", "#0072B2"      # Okabe-Ito, as in figure7_r232.py
CORES = ["R1_A", "R1_B", "R1_C", "R2_A", "R2_B", "R2_C"]

d = json.load(open(HERE / "consistent_basis_final.json"))
nm, nb = d["n_mepn_models"], d["n_betaine_models"]
x = np.array([d["measured_CH4"][c] for c in CORES])
m07, m33, bet = d["areal"]["MePn_0.07_38"], d["areal"]["MePn_0.33_38"], d["areal"]["betaine_18"]
series = [(m07, C_MEPN, "o", f"Methylphosphonate consumers ({nm} models), ratio 0.07", "Methylphosphonate (0.07)"),
          (bet, C_BET, "s", f"Betaine reducers ({nb} models)", "Betaine")]

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
y33 = np.array([m33["sim_by_core"][k] for k in CORES]); ys_all.append(y33)
ax.scatter(x, y33, s=48, marker="o", facecolor="none", edgecolor=C_MEPN, lw=1.1, zorder=3,
           label=f"Methylphosphonate consumers ({nm} models), ratio 0.33 (no fit)")
for xi, y0, y1 in zip(x, ys_all[0], y33):
    ax.plot([xi, xi], [y0, y1], "-", color=C_MEPN, lw=0.6, alpha=0.5, zorder=1)
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
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=1, fontsize=8.5, frameon=False)
fig.tight_layout()
fig.savefig(HERE / "figure7_final.png", bbox_inches="tight")
fig.savefig(HERE / "figure7_final.svg", bbox_inches="tight")
json.dump({"MePn_0.07_38": m07, "MePn_0.33_38": m33, "betaine_18": bet, "measured_CH4": d["measured_CH4"]},
          open(HERE / "figure7_final_stats.json", "w"), indent=1)
print("MePn 0.07", {k: round(v, 3) for k, v in m07.items() if isinstance(v, float)})
print("betaine", {k: round(v, 3) for k, v in bet.items() if isinstance(v, float)})
