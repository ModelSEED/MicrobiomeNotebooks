#!/usr/bin/env python
"""
Figure 8 (co-occurrence / Graphical-Lasso network) -- re-rendered for legibility.

The submitted Figure 8 (docx image14) is scientifically fine but the per-community
spring layout collided node labels badly (e.g. "BM003KS3-K002", "SIND01UBA5335",
"ActinobacteraceaePWM001").  This re-render keeps the same 216-edge Graphical-Lasso
graph (|partial r| > 0.05) and Louvain communities, but:
  * a single global spring layout with a larger repulsion constant,
  * adjustText label repulsion + small white label halos (no overlaps),
  * nodes coloured by GTDB phylum, sized by degree,
  * betaine/MePn consumer nodes outlined and bold,
  * a caveat annotation that the network is inferred from n=12 (p >> n).

Run:  ~/Documents/py_venv/bin/python figure8_relabeled.py
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
from adjustText import adjust_text

REPO = Path(__file__).resolve().parents[1]
CORR = REPO / "correlational"
OUT = Path(__file__).resolve().parent
THRESH = 0.05                        # same |partial r| filter as the pipeline

edges = pd.read_csv(CORR / "glasso_direct_edges.csv")
edges = edges[edges["partial_r"].abs() > THRESH]
phylum_of = json.load(open(CORR / "iterativeID_phylums.json"))
color_of = {k: tuple(v) for k, v in json.load(open(CORR / "Phylum_color_map.json")).items()}

# manuscript-named consumer taxa to highlight
BET = {"Acidaminobacteraceae", "Halanaerobiales", "UBA8670", "UBA6794", "Mcinerneyibacterium"}
MEPN = {"Phaeodactylibacter", "Salinarimonas", "Roseovarius", "Albimonas"}


def phylum(node: str) -> str:
    if node in phylum_of:
        return phylum_of[node]
    for k, v in phylum_of.items():           # strip trailing .N and match
        if k.rsplit(".", 1)[0] == node:
            return v
    return "Unknown"


G = nx.Graph()
for _, r in edges.iterrows():
    G.add_edge(r["taxon_a"], r["taxon_b"], w=abs(r["partial_r"]), r=r["partial_r"])
comms = nx.community.louvain_communities(G, weight="w", resolution=1.0, seed=42)
node2comm = {n: i for i, c in enumerate(comms) for n in c}
Q = nx.community.modularity(G, comms, weight="w")

pos = nx.spring_layout(G, weight="w", k=2.6 / np.sqrt(G.number_of_nodes()), iterations=400, seed=42)

fig, ax = plt.subplots(figsize=(15, 12), dpi=170)
# edges coloured by sign of partial correlation
for u, v, d in G.edges(data=True):
    ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
            color=("#2166AC" if d["r"] > 0 else "#B2182B"),
            lw=0.4 + 4.0 * d["w"], alpha=0.35, zorder=1)
# nodes
deg = dict(G.degree())
for n in G.nodes():
    ph = phylum(n)
    col = color_of.get(ph, (0.7, 0.7, 0.7, 1.0))
    is_con = n in BET or n in MEPN
    ax.scatter([pos[n][0]], [pos[n][1]], s=60 + 55 * deg[n],
               color=col, edgecolor=("black" if is_con else "white"),
               linewidth=(2.0 if is_con else 0.6), zorder=3)
# labels (repelled)
texts = []
for n in G.nodes():
    con = n in BET or n in MEPN
    texts.append(ax.text(pos[n][0], pos[n][1], n,
                         fontsize=(9 if con else 7.2),
                         fontweight=("bold" if con else "normal"),
                         color=("black" if con else "0.15"),
                         bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.7)))
adjust_text(texts, ax=ax, expand=(1.25, 1.5),
            arrowprops=dict(arrowstyle="-", color="0.6", lw=0.4))

ax.axis("off")
ax.set_title("Figure 8 (re-rendered)  ·  Direct interactions: Graphical-Lasso partial correlations\n"
             f"within-condition residuals · {G.number_of_nodes()} taxa, {G.number_of_edges()} edges · "
             f"Louvain Q = {Q:.2f} · {len(comms)} communities",
             fontsize=14, fontweight="bold")

# legends
phyla_present = sorted({phylum(n) for n in G.nodes()})
ph_handles = [Line2D([0], [0], marker="o", ls="", mfc=color_of.get(p, (.7, .7, .7, 1)),
                     mec="white", ms=10, label=p) for p in phyla_present if p != "Unknown"]
leg1 = ax.legend(handles=ph_handles, loc="upper left", fontsize=8, title="Phylum",
                 ncol=2, framealpha=0.9)
ax.add_artist(leg1)
extra = [Line2D([0], [0], color="#2166AC", lw=3, label="positive partial r (coupling)"),
         Line2D([0], [0], color="#B2182B", lw=3, label="negative partial r (partitioning)"),
         Line2D([0], [0], marker="o", ls="", mfc="0.8", mec="black", mew=2, ms=11,
                label="betaine / MePn consumer")]
ax.legend(handles=extra, loc="lower left", fontsize=9, framealpha=0.9)
ax.text(0.5, -0.02, "Caveat: precision matrix inferred from n = 12 samples over ~140 taxa (p >> n); "
        "edges are unstable — treat as exploratory (bootstrap edge-stability recommended).",
        transform=ax.transAxes, ha="center", fontsize=9, style="italic", color="0.35")

plt.tight_layout()
fig.savefig(OUT / "figure8_relabeled.png", bbox_inches="tight")
fig.savefig(OUT / "figure8_relabeled.svg", bbox_inches="tight")
plt.close(fig)
print(f"nodes={G.number_of_nodes()} edges={G.number_of_edges()} communities={len(comms)} Q={Q:.3f}")
print("community sizes:", sorted((len(c) for c in comms), reverse=True))
print("Saved figure8_relabeled.png/.svg")
