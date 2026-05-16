"""Louvain-modularity co-occurrence network on the same Glasso edges that feed
within_condition_glasso_network.png. Layout strategy is ported verbatim (in
spirit) from the codiffusion_bioreactor data_processing.ipynb cell:
    1. Louvain communities on the positive-partial_r subgraph, weighted by |partial_r|.
    2. Per-module independent spring_layout, placed around a ring whose radius
       is sized so module hulls don't overlap.
    3. SAT-shrunk convex-hull shading per module.
    4. Singleton (negative-only) nodes interleaved on an outer ring.

Inputs (already on disk, identical to the Glasso figure's inputs):
    glasso_direct_edges.csv          - edges from GraphicalLassoCV on CLR residuals
    within_condition_meta.json       - taxa list / condition labels
    averaged_normalized_MAG_abundances.csv  - for mean rel. abundance -> node size
    iterativeIDs.json                - MAG-id -> short taxonomic id
    ../data/Saltern_phylogeny.json   - short id -> Phylum (node color)

Outputs (new naming, parallel to within_condition_glasso_network.png):
    within_condition_louvain_network.png
    within_condition_louvain_modules.json
"""

import math
from json import dump, load
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import cm as mcm
from scipy.spatial import ConvexHull

HERE = Path(__file__).resolve().parent
PARTIAL_THRESHOLD = 0.05          # same edge filter as the Glasso figure
LOUVAIN_RESOLUTION = 1.0          # codiffusion default
LOUVAIN_SEED = 42

# ----- load Glasso edges + meta ----------------------------------------------
edges_df = pd.read_csv(HERE / "glasso_direct_edges.csv")
meta = load(open(HERE / "within_condition_meta.json"))
taxa_in_panel = meta["taxa"]
print(f"loaded {len(edges_df)} Glasso edges across {len(taxa_in_panel)} taxa")

# ----- short-id -> phylum (node color) ---------------------------------------
iterativeIDs = load(open(HERE / "iterativeIDs.json"))
phylogeny    = load(open(HERE.parent / "data" / "Saltern_phylogeny.json"))

def _to_short(full_id):
    bare = full_id
    for suf in (".contigs__", ".contigs"):
        if bare.endswith(suf):
            bare = bare[: -len(suf)]
            break
    iid = iterativeIDs.get(bare, bare)
    return ".".join(iid.split(".")[:-1]) if "." in iid else iid

id_to_phylum = {}
for full_id, taxonomy in phylogeny.items():
    short = _to_short(full_id)
    ph = taxonomy.get("Phylum", "Unknown") if isinstance(taxonomy, dict) else "Unknown"
    ph = ph or "Unknown"
    if id_to_phylum.get(short, "Unknown") in ("", "Unknown", None):
        id_to_phylum[short] = ph

# ----- mean rel. abundance per short-taxon (-> node size) --------------------
mag_abund = pd.read_csv(HERE / "averaged_normalized_MAG_abundances.csv", index_col=0)
mag_short = mag_abund.index.to_series().map(_to_short)
per_short = mag_abund.groupby(mag_short).sum()
per_sample_rel = per_short.div(per_short.sum(axis=0), axis=1)
mean_rel_abund = per_sample_rel.mean(axis=1).to_dict()

# ----- build graph from the same Glasso edges --------------------------------
filtered = edges_df[edges_df["partial_r"].abs() > PARTIAL_THRESHOLD]
print(f"|partial_r| > {PARTIAL_THRESHOLD}: kept {len(filtered)} of {len(edges_df)} edges")

G = nx.Graph()
for _, row in filtered.iterrows():
    G.add_edge(row["taxon_a"], row["taxon_b"],
               partial_r=row["partial_r"], weight=abs(row["partial_r"]))
G.remove_nodes_from(list(nx.isolates(G)))
print(f"graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# ----- Louvain on positive-partial_r subgraph --------------------------------
pos_sub = G.edge_subgraph(
    [(u, v) for u, v, d in G.edges(data=True) if d["partial_r"] > 0]
).copy()
layout_comms = list(nx.community.louvain_communities(
    pos_sub, weight="weight",
    resolution=LOUVAIN_RESOLUTION, seed=LOUVAIN_SEED))

# nodes only connected by negative edges become singleton modules
node_module = {n: i for i, c in enumerate(layout_comms) for n in c}
isolates = [n for n in G.nodes() if n not in node_module]
for iso in isolates:
    node_module[iso] = len(layout_comms)
    layout_comms.append({iso})

multi_comms = [c for c in layout_comms if len(c) >= 2]
singleton_comms = [c for c in layout_comms if len(c) == 1]
Q = nx.community.modularity(pos_sub, multi_comms, weight="weight")
print(f"Louvain: {len(multi_comms)} multi-node modules + "
      f"{len(singleton_comms)} singletons, Q = {Q:.3f}")
print(f"sizes (louvain): {sorted([len(c) for c in multi_comms], reverse=True)}")

# ----- module-membership export ----------------------------------------------
module_export = {}
louv_n = sing_n = 0
for comm in sorted(layout_comms, key=lambda c: -len(c)):
    members = sorted(comm)
    if len(comm) >= 2:
        louv_n += 1
        module_export[f"module_{louv_n}"] = {
            "size": len(comm), "type": "louvain", "members": members}
    else:
        sing_n += 1
        module_export[f"singleton_{sing_n}"] = {
            "size": 1, "type": "isolate_singleton", "members": members}
dump(module_export, open(HERE / "within_condition_louvain_modules.json", "w"), indent=2)
print(f"wrote within_condition_louvain_modules.json "
      f"({louv_n} modules + {sing_n} singletons)")

# ----- per-module spring layout + ring placement -----------------------------
np.random.seed(LOUVAIN_SEED)
scale_factor = 20.0
module_layouts = {}
for mi, comm in enumerate(layout_comms):
    if len(comm) < 2:
        continue
    sub = G.subgraph(comm).copy()
    module_scale = math.sqrt(len(comm)) * scale_factor
    if sub.number_of_edges() == 0:
        module_layouts[mi] = {n: np.random.normal(0, module_scale * 0.3, 2) for n in comm}
    else:
        module_layouts[mi] = nx.spring_layout(
            sub, seed=LOUVAIN_SEED, iterations=3000, k=5.0,
            scale=module_scale, weight="weight")

if module_layouts:
    max_module_radius = max(
        max(np.linalg.norm(p) for p in module_layouts[mi].values())
        for mi in module_layouts)
else:
    max_module_radius = 1.0

hull_buffer = 1.02
gap_factor  = 1.01
n_multi_mods = len(module_layouts)
min_chord = 2 * max_module_radius * hull_buffer * gap_factor
sep_radius = (min_chord / (2 * np.sin(np.pi / n_multi_mods))
              if n_multi_mods > 1 else 0)
print(f"per-module layout: ring radius={sep_radius:.2f}, "
      f"max module radius={max_module_radius:.2f}")

pos = {}
sorted_mis = sorted(module_layouts.keys(), key=lambda mi: -len(layout_comms[mi]))
for i, mi in enumerate(sorted_mis):
    theta = 2 * np.pi * i / n_multi_mods
    center = np.array([np.cos(theta), np.sin(theta)]) * sep_radius
    for n, local in module_layouts[mi].items():
        pos[n] = center + local

singleton_mis = [mi for mi, c in enumerate(layout_comms) if len(c) == 1]
if singleton_mis:
    singleton_radius = sep_radius + max_module_radius * 2.5
    for i, mi in enumerate(singleton_mis):
        theta = 2 * np.pi * (i + 0.5) / max(len(singleton_mis), 1)
        for n in layout_comms[mi]:
            pos[n] = np.array([np.cos(theta), np.sin(theta)]) * singleton_radius

# ----- styling: phylum colors, node sizes, edge colors -----------------------
present_phyla = sorted({id_to_phylum.get(t, "Unknown") for t in G.nodes()})
archaea_keywords = ("archae", "halobacterota", "methanobacteriota",
                    "thermoplasmatota", "micrarchaeota", "asgardarchaeota",
                    "nanohaloarchaeota")
archaea_phyla  = sorted([p for p in present_phyla
                         if p and any(k in p.lower() for k in archaea_keywords)])
bacteria_phyla = sorted([p for p in present_phyla
                         if p and p not in archaea_phyla and p != "Unknown"])
taxa_color_map = {}
for i, ph in enumerate(archaea_phyla):
    taxa_color_map[ph] = plt.cm.turbo(0.05 + 0.20 * (i / max(len(archaea_phyla), 1)))
for i, ph in enumerate(bacteria_phyla):
    taxa_color_map[ph] = plt.cm.turbo(0.30 + 0.65 * (i / max(len(bacteria_phyla), 1)))
DEFAULT_COLOR = "lightgray"

node_list   = list(G.nodes())
NODE_BASE, NODE_SCALE = 800, 28000
node_sizes  = np.array([NODE_BASE + NODE_SCALE * np.sqrt(max(mean_rel_abund.get(n, 0), 0))
                        for n in node_list])
diam_pt     = 2 * np.sqrt(node_sizes / np.pi)
node_colors = [taxa_color_map.get(id_to_phylum.get(n, "Unknown"), DEFAULT_COLOR)
               for n in node_list]

pcorrs   = np.array([d["partial_r"] for _, _, d in G.edges(data=True)])
p_max    = max(np.abs(pcorrs).max(), 1e-6)
widths   = 0.6 + 4.0 * np.abs(pcorrs) / p_max
norm     = mcolors.TwoSlopeNorm(vmin=-p_max, vcenter=0, vmax=p_max)
edge_cmap = mcm.RdBu
ecolors  = [edge_cmap(norm(r)) for r in pcorrs]

# ----- render -----------------------------------------------------------------
fig, ax = plt.subplots(figsize=(34, 26))
nx.draw_networkx_edges(G, pos, width=widths, edge_color=ecolors, alpha=0.7, ax=ax)
nx.draw_networkx_nodes(G, pos, nodelist=node_list,
                       node_size=node_sizes, node_color=node_colors,
                       edgecolors="black", linewidths=0.4, ax=ax)

for n, d_pt in zip(node_list, diam_pt):
    x, y = pos[n]
    ax.text(x, y, str(n), fontsize=max(6, 0.45 * d_pt),
            color="black", fontweight="normal",
            ha="center", va="center", clip_on=False)

# ----- SAT-shrunk convex-hull shading per module -----------------------------
palette_multi = plt.cm.tab10(np.linspace(0, 1, 10))
module_palette = []
multi_idx = 0
for c in layout_comms:
    if len(c) >= 2:
        module_palette.append(palette_multi[multi_idx % 10])
        multi_idx += 1
    else:
        module_palette.append((0.65, 0.65, 0.65, 1.0))

def _convex_polys_intersect(p1, p2):
    for poly, other in ((p1, p2), (p2, p1)):
        for i in range(len(poly)):
            edge = poly[(i + 1) % len(poly)] - poly[i]
            normal = np.array([-edge[1], edge[0]])
            proj_a = poly @ normal
            proj_b = other @ normal
            if proj_a.max() < proj_b.min() or proj_b.max() < proj_a.min():
                return False
    return True

hull_data = {}
initial_scale = 1.30
for mi, comm in enumerate(layout_comms):
    if len(comm) < 2:
        continue
    pts = np.array([pos[n] for n in comm if n in pos])
    centroid = pts.mean(axis=0)
    if len(pts) >= 3:
        h = ConvexHull(pts)
        verts = pts[h.vertices]
    else:
        r = max(0.05, np.linalg.norm(pts - centroid, axis=1).max())
        theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
        verts = centroid + r * np.column_stack([np.cos(theta), np.sin(theta)])
    hull_data[mi] = [centroid, verts, initial_scale]

min_scale = 0.75
for it in range(120):
    overlap = False
    ids = list(hull_data.keys())
    for ia in range(len(ids)):
        mi_a = ids[ia]
        ca, va, sa = hull_data[mi_a]
        poly_a = ca + (va - ca) * sa
        for ib in range(ia + 1, len(ids)):
            mi_b = ids[ib]
            cb, vb, sb = hull_data[mi_b]
            poly_b = cb + (vb - cb) * sb
            if _convex_polys_intersect(poly_a, poly_b):
                overlap = True
                if hull_data[mi_a][2] > min_scale:
                    hull_data[mi_a][2] = max(min_scale, hull_data[mi_a][2] * 0.95)
                if hull_data[mi_b][2] > min_scale:
                    hull_data[mi_b][2] = max(min_scale, hull_data[mi_b][2] * 0.95)
    if not overlap:
        print(f"hull non-overlap converged after {it + 1} iterations")
        break
else:
    print(f"hull non-overlap floor reached at scale {min_scale}")

for mi, (ca, verts, sv) in hull_data.items():
    poly_pts = ca + (verts - ca) * sv
    poly = mpatches.Polygon(poly_pts, closed=True,
                            facecolor=module_palette[mi],
                            edgecolor=module_palette[mi],
                            linewidth=2, alpha=0.22, zorder=0)
    ax.add_patch(poly)

# ----- titles, colorbar, legends ---------------------------------------------
ax.set_axis_off()
ax.set_title("Louvain modularity on the Graphical-Lasso direct-interaction network\n"
             "(same edges as within_condition_glasso_network.png; per-module "
             "spring layout, ring placement)",
             fontsize=22, pad=20)

sm = mcm.ScalarMappable(cmap=edge_cmap, norm=norm); sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, shrink=0.4, pad=0.01, location="right")
cbar.set_label("partial correlation (controlling for all other taxa)", fontsize=20)
cbar.ax.tick_params(labelsize=16)
cbar.ax.yaxis.set_label_position("left")
cbar.ax.set_position([0.78, 0.08, 0.018, 0.30])

def header_patch(title):
    return mpatches.Patch(color="none", label=fr"$\bf{{{title}}}$")
handles = []
if archaea_phyla:
    handles.append(header_patch("Archaea"))
    handles += [mpatches.Patch(color=taxa_color_map[p], label=p) for p in archaea_phyla]
if bacteria_phyla:
    handles.append(header_patch("Bacteria"))
    handles += [mpatches.Patch(color=taxa_color_map[p], label=p) for p in bacteria_phyla]
leg1 = ax.legend(handles=handles, title="Phylum",
                 loc="upper left", bbox_to_anchor=(0.04, 0.55),
                 bbox_transform=fig.transFigure,
                 fontsize=16, title_fontsize=20, frameon=True)
ax.add_artist(leg1)

# horizontal abundance legend (matches Glasso figure)
legend_entries = [(0.001, "0.1%"), (0.01, "1%"), (0.05, "5%")]
sizes_pt2_leg  = [NODE_BASE + NODE_SCALE * np.sqrt(a) for a, _ in legend_entries]
diameters_pt   = [2 * np.sqrt(s / np.pi) for s in sizes_pt2_leg]
breather_pt    = 18
offsets_pt     = [diameters_pt[0] / 2]
for k in range(1, len(legend_entries)):
    offsets_pt.append(offsets_pt[-1]
                      + (diameters_pt[k - 1] + diameters_pt[k]) / 2 + breather_pt)
total_pt   = offsets_pt[-1] + diameters_pt[-1] / 2
fig_w_pts  = fig.get_figwidth() * 72
fig_h_pts  = fig.get_figheight() * 72
x_start    = 0.04
y_circles  = 0.91
title_y    = y_circles + max(diameters_pt) / fig_h_pts / 2 + 0.01
label_y    = y_circles - max(diameters_pt) / fig_h_pts / 2 - 0.018
title_x    = x_start + (total_pt / fig_w_pts) / 2
ax.text(title_x, title_y, "Mean rel. abundance",
        transform=fig.transFigure, ha="center", va="bottom",
        fontweight="bold", fontsize=18, clip_on=False)
for (abund, label), s, off in zip(legend_entries, sizes_pt2_leg, offsets_pt):
    x = x_start + off / fig_w_pts
    ax.scatter([x], [y_circles], s=s, color="slategray", alpha=0.85,
               edgecolor="black", linewidth=0.4,
               transform=fig.transFigure, clip_on=False)
    ax.text(x, label_y, label,
            transform=fig.transFigure, ha="center", va="top",
            fontsize=16, clip_on=False)

# modularity summary box
summary = (r"$\bf{Modularity}$" + "\n"
           + f"Q = {Q:.3f}" + "\n"
           + f"modules = {len(multi_comms)} louvain + "
             f"{len(singleton_comms)} singletons (total = {len(layout_comms)})\n"
           + f"sizes (louvain): {sorted([len(c) for c in multi_comms], reverse=True)}")
ax.text(-0.10, 0.7, summary, transform=ax.transAxes,
        fontsize=18, ha="left", va="bottom",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                  edgecolor="black", alpha=0.85, linewidth=1.5))

fig.savefig(HERE / "within_condition_louvain_network.png",
            dpi=300, bbox_inches="tight")
plt.close(fig)
print("wrote within_condition_louvain_network.png")
