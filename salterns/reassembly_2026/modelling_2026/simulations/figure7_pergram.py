#!/usr/bin/env python
"""Figure 7 (per-gram framing): (a) methane produced per gram of biomass formed by each guild model; (b) the doubling
time each guild's measured biomass would need to supply the measured methane of each unrestored core.
Inputs: per_gram_yields.tsv, specific_rates.tsv / .json (specific_rates.py), ../models/methanogen/curation_k.json.
Colours: the dataviz reference palette in fixed slot order (slot 1 blue = MePn consumers, as an ordinal blue ramp for
the three ordered phosphorus scenarios; slot 2 orange = betaine reducers; slot 3 aqua = methanogen), each series also
with its own marker shape, so identity never rests on colour alone.
Run: ~/Documents/py_venv/bin/python figure7_pergram.py -> figure7_pergram.png / .svg
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
INK, INK2, GRID = "#1f1f1f", "#5a5a5a", "#e6e6e6"
SER = {"MePn, ratio 0.07": ("#6da7ec", "o"), "MePn, ratio 0.33": ("#2a78d6", "^"), "MePn, only P source": ("#184f95", "s"),
       "Betaine, 0.06 CH4 per TMA": ("#eb6834", "D"), "Methanogen, methylamines": ("#1baf7a", "P")}
LAB = {"MePn, ratio 0.07": "MPn consumers, MPn:Pi 0.07", "MePn, ratio 0.33": "MPn consumers, MPn:Pi 0.33",
       "MePn, only P source": "MPn consumers, MPn only P source", "Betaine, 0.06 CH4 per TMA": "Betaine reducers (0.06 CH$_4$ per TMA)",
       "Methanogen, methylamines": "$Methanolobus$ (methylamines)"}
CORES = ["R1_A", "R1_B", "R1_C", "R2_A", "R2_B", "R2_C"]

Y = pd.read_csv(HERE / "per_gram_yields.tsv", sep="\t")
YM, YB = Y[Y.guild == "mepn_consumer"], Y[Y.guild == "betaine_reducer"]
km = json.load(open(HERE.parent / "models/methanogen/curation_k.json"))
sr = pd.read_csv(HERE / "specific_rates.tsv", sep="\t")
sr = sr[sr.scenario.isin(SER)]

fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.9), dpi=300, gridspec_kw={"width_ratios": [1, 1.25], "wspace": 0.32})
rng = np.random.default_rng(0)
cols = [("MePn, ratio 0.07", YM["k_CH4_r0.07"]), ("MePn, ratio 0.33", YM["k_CH4_r0.33"]),
        ("MePn, only P source", YM["k_CH4_MePn_only_P"]), ("Betaine, 0.06 CH4 per TMA", 0.06 * YB["k_TMA_unbounded_min"]),
        ("Methanogen, methylamines", pd.Series([km["methylamines"]]))]
for i, (name, vals) in enumerate(cols):
    c, mk = SER[name]
    v = vals.dropna().values
    x = i + rng.uniform(-0.18, 0.18, len(v)) if len(v) > 1 else np.array([i])
    a.scatter(x, v, s=34 if len(v) > 1 else 90, marker=mk, color=c, edgecolor="white" if len(v) > 1 else c,
              linewidths=0.8, zorder=5 if len(v) == 1 else 3)
    a.hlines(np.median(v), i - 0.3, i + 0.3, color=INK, lw=1.6, zorder=4)
    if name.startswith("Methanogen"):
        lo, hi = km["coupling_range"]
        a.vlines(i, lo, hi, color=INK2, lw=1.2, zorder=2)
    a.text(i + 0.33, np.median(v), f"{np.median(v):.3g}", ha="left", va="center", fontsize=8, color=INK)
a.set_yscale("log"); a.set_ylim(0.01, 2000)
SHORT = ["MPn consumers\nMPn:Pi 0.07", "MPn consumers\nMPn:Pi 0.33", "MPn consumers\nMPn only P", "Betaine\nreducers", "$Methanolobus$"]
a.set_xticks(range(len(cols))); a.set_xticklabels(SHORT, fontsize=7.8, rotation=30, ha="right", rotation_mode="anchor")
a.set_xlim(-0.6, len(cols) - 0.2)
a.set_ylabel("CH$_4$ produced per gram of biomass formed (mmol g$^{-1}$)", fontsize=9.5, color=INK)
a.grid(True, axis="y", which="major", color=GRID, lw=0.7); a.set_axisbelow(True)
a.set_title("a  Methane per gram of biomass formed", loc="left", fontsize=11, color=INK)
for sp in ("top", "right"):
    a.spines[sp].set_visible(False)

offs = np.linspace(-0.28, 0.28, len(SER))
for (name, (c, mk)), off in zip(SER.items(), offs):
    d = sr[sr.scenario == name].set_index("core").reindex(CORES)
    b.scatter(np.arange(6) + off, d.doubling_time_h, s=40, marker=mk, color=c, edgecolor="white", linewidths=0.8,
              zorder=3, label=LAB[name])
for h, lab in ((1 / 6, "10 min"), (1, "1 h"), (24, "1 d"), (24 * 7, "1 wk")):
    b.axhline(h, color=INK2, lw=0.8, ls=":", zorder=1)
    b.text(-0.45, h * 1.08, lab, fontsize=7.5, color=INK2, ha="left", va="bottom")
b.set_yscale("log")
b.set_xticks(range(6)); b.set_xticklabels([c.replace("_", " ") for c in CORES], fontsize=9)
b.set_ylabel("Doubling time needed to supply measured CH$_4$ (h)", fontsize=9.5, color=INK)
b.grid(True, axis="y", which="major", color=GRID, lw=0.7); b.set_axisbelow(True)
b.set_title("b  Growth needed in each unrestored core", loc="left", fontsize=11, color=INK)
for sp in ("top", "right"):
    b.spines[sp].set_visible(False)
b.legend(loc="upper center", bbox_to_anchor=(0.45, -0.12), ncol=2, fontsize=8, frameon=False)
fig.savefig(HERE / "figure7_pergram.png", bbox_inches="tight")
fig.savefig(HERE / "figure7_pergram.svg", bbox_inches="tight")
print("wrote figure7_pergram.png")
