#!/usr/bin/env python
"""
Build-out figure: taxonomic breakdown of the modelled methane.

Turns the abstract Figure-7 R^2 into a taxon-resolved, mechanistic result by
showing WHICH families carry the modelled methylphosphonate-derived methane in each
unrestored core.  This directly visualises the manuscript's claim (Results, para on
metabolic modelling) that Rhodobacteraceae -- Roseovarius / Albimonas -- dominate the
modelled MePn methane and are uniquely abundant in the unrestored ponds.

Per-MAG modelled emission = fluxes_to_emissions(that MAG's own CH4 flux) at the
best-fit MePn:phosphate ratio 0.07, grouped by GTDB family.
Run:  ~/Documents/py_venv/bin/python family_contribution.py
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
DATA, CORR = REPO / "data", REPO / "correlational"
OUT = Path(__file__).resolve().parent

DNA = pd.read_csv(DATA / "SaltPondsDNA.csv").set_index("Sample")["DNA conc (µg/g soil)"]
PREFIXES = sorted({"_".join(i.split("_")[:-1]) for i in DNA.index})
VW_DNA = pd.Series({p: (5 * DNA[p + "_D1"] + 10 * DNA[p + "_D2"]) / 15 for p in PREFIXES})


def to_emission(row):
    area_m2 = math.pi * (5 / 2) ** 2 * 0.01 ** 2
    ug = VW_DNA.reindex(row.index).fillna(2.0)
    gbio = 0.7 * 294.52 * ug / 0.1 / 1e6
    return row / area_m2 * 1000.0 * 24 * gbio


tax = json.load(open(DATA / "disaggregated_taxonomy.json"))
mepn = json.load(open(CORR / "saltern_fluxes.json"))["0.07"]["methane"]   # {MAG:{sample:flux}}

# per-MAG emission per sample -> genus (family as fallback)
rows = []
for mag, byS in mepn.items():
    t = tax.get(mag.replace(".contigs", ""), {})
    lab = t.get("genus") or t.get("family") or "Unclassified"
    em = to_emission(pd.Series(byS).abs())
    rows.append((lab, em))
fam_df = pd.DataFrame({i: r[1] for i, r in enumerate(rows)}).T
fam_df["family"] = [r[0] for r in rows]
by_family = fam_df.groupby("family").sum(numeric_only=True)
by_family = by_family.loc[by_family.sum(axis=1).sort_values(ascending=False).index]  # biggest first
samples = [c for c in by_family.columns if c in PREFIXES]
by_family = by_family[samples]

# keep top families, lump the rest
TOPN = 6
top = by_family.index[:TOPN].tolist()
plot_df = by_family.loc[top].copy()
if len(by_family) > TOPN:
    plot_df.loc["Other"] = by_family.iloc[TOPN:].sum()

# Okabe-Ito-ish qualitative set (colourblind-safe), Rhodobacteraceae -> vermillion
palette = ["#D55E00", "#0072B2", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#999999"]
colors = {f: palette[i % len(palette)] for i, f in enumerate(plot_df.index)}

fig, ax = plt.subplots(figsize=(9, 5.6), dpi=200)
bottom = np.zeros(len(plot_df.columns))
x = np.arange(len(plot_df.columns))
for fam in plot_df.index:
    vals = plot_df.loc[fam].values
    ax.bar(x, vals, bottom=bottom, width=0.72, label=fam, color=colors[fam],
           edgecolor="white", linewidth=0.8)
    bottom += vals
ax.set_xticks(x); ax.set_xticklabels(plot_df.columns, fontsize=11)
ax.set_ylabel(r"Modelled MePn-derived CH$_4$" + "\n" + r"[$\mu$mol m$^{-2}$ d$^{-1}$]", fontsize=12)
ax.set_xlabel("Unrestored saltern core", fontsize=12)
ax.set_title("Modelled methylphosphonate methane by genus\n(MePn:P ratio 0.07, biomass = in-situ abundance)",
             fontsize=13, fontweight="bold")
ax.legend(title="GTDB genus", fontsize=9, title_fontsize=10, bbox_to_anchor=(1.01, 1), loc="upper left")
ax.grid(axis="y", linestyle="--", alpha=0.3)
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
fig.savefig(OUT / "family_contribution.png", bbox_inches="tight")
fig.savefig(OUT / "family_contribution.svg", bbox_inches="tight")
plt.close(fig)

print("Family share of total modelled MePn methane (all unrestored cores):")
share = (by_family.sum(axis=1) / by_family.sum().sum() * 100).round(1)
for fam, pct in share.items():
    print(f"  {pct:5.1f}%  {fam}")
print("\nSaved family_contribution.png/.svg")
