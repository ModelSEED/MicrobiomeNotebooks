#!/usr/bin/env python
"""Salterns correlational studies on the FINAL re-assembled MAG set (2026-11-19): 585 species representatives.

Same statistics and parameters as ../r232/rerun_correlational_r232.py (the 310 original MAGs, GTDB R232 labels); only
the community changes:
  * members   - the 585 species representatives of the final set (935 second-pass MAGs + 520 kept first-pass site
                co-assembly MAGs; skani 95 % ANI / 50 % AF species clusters; ../../reassembly_2026/comparison_rollout/).
  * abundance - CoverM 0.7.0 read counts of the 39 metagenomes against the 585 representatives
                (../../reassembly_2026/abundance_final/final585_relabund_24sediment.tsv: fraction of QC'd reads, 24
                sediment samples, plus unbinned_reads). Same DNA-weighted D1/D2 merge into 12 cores.
  * labels    - iterativeIDs_final585.json (make_ids_final585.py); species with an original MAG keep its R232 ID.
  * guilds    - KOfam 2026-10-26 calls, a species is in a guild when ANY of its MAGs carries the marker
                (../../reassembly_2026/comparison_rollout/final_guilds_species.tsv, final_set_guilds.py):
                  MePn consumer   phnJ                 betaine reducer  grdH or grdI
                  methanogen      mcrA                 TMA producer     cutC+cutD or cntA+cntB (none in the set)
                  TMA candidate   cutC or cntA without its partner subunit
                The representative's own call is kept next to it (<guild>_rep) as a sensitivity column.
                The 2020 curated modelling sets are gone: the 2026 models are built for every MePn and betaine species
                (comparison_rollout/modelling_species.tsv), so "modelled set" == guild set and is not repeated.
Outputs: this folder, suffix _final585. Nothing outside it is modified.
Run (from this folder): ~/Documents/py_venv/bin/python rerun_correlational_final585.py [--parts guilds,corr,heatmap,network,table3]
"""
import argparse
import collections
import json
import subprocess
import sys
from itertools import combinations, permutations, product
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
COR = HERE.parent
S = COR.parent
RA = S / "reassembly_2026"
CR = RA / "comparison_rollout"
TAG = "final585"
ZERO = 1e-5                     # detection level used throughout the 2020 code
N_NULL = 4000                   # random-subset draws (as in reassessment_figure.py)
REST = {"unrestored": ("#D55E00", "o"), "restored": ("#009E73", "s"), "reference": ("#56B4E9", "^")}


def restoration(core):
    return "restored" if core.startswith("SF2") else "reference" if core.startswith("R2A") else "unrestored"


# ---- shared inputs ------------------------------------------------------------------------------
ids = json.load(open(HERE / f"iterativeIDs_{TAG}.json"))                    # representative MAG -> ID
taxa_of = json.load(open(HERE / f"phylogeny_{TAG}.json"))                   # representative MAG -> ranks
rel = pd.read_csv(RA / "abundance_final/final585_relabund_24sediment.tsv", sep="\t", index_col=0)
dna = pd.read_csv(S / "data/SaltPondsDNA.csv").set_index("Sample")["DNA conc (µg/g soil)"]
meta = pd.read_csv(S / "data/Cliff_Sample_Metadata_BGC_NMR.csv").set_index("Sample")
MAGS = [m for m in rel.index if m != "unbinned_reads"]
assert len(MAGS) == 585 and set(MAGS) == set(ids), "relabund rows must be the 585 representatives"
assert set(rel.columns) <= set(dna.index)
PRE = sorted({"_".join(s.split("_")[:-1]) for s in dna.index})        # 12 cores
UNRESTORED = [p for p in PRE if p.startswith(("R1_", "R2_"))]         # R2A_* is the reference site


def merge_depths(df):
    """DNA-weighted D1/D2 merge, as merge_depths_by_DNA in analysis_2_marimo.py."""
    return pd.DataFrame({p: (df[f"{p}_D1"] * dna[f"{p}_D1"] + df[f"{p}_D2"] * dna[f"{p}_D2"])
                         / (dna[f"{p}_D1"] + dna[f"{p}_D2"]) for p in PRE})


ch4 = meta["CH4_umol_m2_d"].astype(float)
AB = {"raw": merge_depths(rel), "DNA-normalised": merge_depths(rel.div(dna[rel.columns], axis=1))}
CH4 = {"raw": pd.Series({p: ch4[[f"{p}_D1", f"{p}_D2"]].mean() for p in PRE}),
       "DNA-normalised": pd.Series({p: (ch4[[f"{p}_D1", f"{p}_D2"]] / dna[[f"{p}_D1", f"{p}_D2"]]).mean()
                                    for p in PRE})}
GCALL = ["mepn_consumer", "betaine_reducer", "methanogen", "methylamine_user", "tma_producer", "tma_candidate_one_sided"]


# ---- 1. guilds ------------------------------------------------------------------------------------
def build_guilds():
    d = pd.read_csv(CR / "mags_three_way.tsv", sep="\t")
    reps = d[d.final_species_rep == True].copy()
    reps["mag"] = reps.genome.str.split("|", n=1).str[1]
    sp = pd.read_csv(CR / "final_guilds_species.tsv", sep="\t", index_col=0)
    fg = pd.read_csv(CR / "final_guilds.tsv", sep="\t")
    for c in GCALL:
        fg[c] = fg[c].astype("boolean")
    has = fg.kofam_source != "missing"
    assert has.all(), "KOfam calls missing for some final-set MAGs"
    per = fg.groupby("species_cluster")
    n_mags = per.size()
    n_call = per[GCALL].sum()
    rep_call = fg[fg.final_species_rep].set_index("species_cluster")[GCALL]
    ms = pd.read_csv(CR / "modelling_species.tsv", sep="\t")
    sc_col = "species_cluster" if "species_cluster" in ms.columns else ms.columns[0]
    modelled = set(ms[sc_col])
    orig = pd.read_csv(HERE / f"iterativeIDs_{TAG}.csv").set_index("mag")
    rows = []
    for r in reps.itertuples():
        t, sc = taxa_of[r.mag], r.species_cluster
        row = {"mag": r.mag, "set": r.set, "species_cluster": sc, f"iterativeID_{TAG}": ids[r.mag],
               "originals_in_cluster_r232": orig.loc[r.mag, "originals_in_cluster_r232"]
               if isinstance(orig.loc[r.mag, "originals_in_cluster_r232"], str) else "",
               "domain": t["Kingdom"], "phylum": t["Phylum"], "class": t["Class"], "order": t["Order"],
               "family": t["Family"], "genus": t["Genus"], "species": t["Species"],
               "completeness": r.completeness, "contamination": r.contamination, "mimag": r.mimag,
               "n_mags_in_species": int(n_mags[sc])}
        for c in GCALL:
            row[c] = bool(sp.loc[sc, c])
            row[f"{c}_rep"] = bool(rep_call.loc[sc, c])
            row[f"{c}_n_mags"] = int(n_call.loc[sc, c])
        row |= {"archaea": t["Kingdom"] == "Archaea", "methanosarcinales": t["Order"] == "Methanosarcinales",
                "modelled_2026": sc in modelled, "mean_rel_abund_12cores": AB["raw"].loc[r.mag].mean()}
        rows.append(row)
    g = pd.DataFrame(rows)
    assert len(g) == 585 and g.mag.is_unique
    g.to_csv(HERE / f"guilds_{TAG}.csv", index=False)
    sizes = {c: int(g[c].sum()) for c in GCALL + ["archaea", "methanosarcinales", "modelled_2026"]}
    print("guild sizes (species):", sizes)
    print("  representative's own call:", {c: int(g[f"{c}_rep"].sum()) for c in GCALL})
    for c in ("mepn_consumer", "betaine_reducer"):
        print(f"  modelled_2026 vs {c}: {int((g[c] & ~g.modelled_2026).sum())} {c} not modelled, "
              f"{int((g.modelled_2026 & ~g[c]).sum())} modelled without the call (other guild)")
    print("methanogen species:", g.loc[g.methanogen, [f"iterativeID_{TAG}", "species", "methanogen_n_mags", "n_mags_in_species"]].to_string(index=False))
    return g


GUILDS = [("mepn_consumer", "MePn consumers (phnJ)"), ("betaine_reducer", "Betaine reducers (grdH or grdI)"),
          ("methanogen", "Methanogens (mcrA)"), ("tma_candidate_one_sided", "TMA-route candidates (cutC or cntA only)"),
          ("methanosarcinales", "Methanosarcinales"), ("archaea", "All archaea")]


def active_guilds(g):
    """Guilds with members; a guild whose membership equals an earlier one is dropped (e.g. methanogen == order)."""
    out, seen = [], []
    for col, lab in GUILDS:
        mem = frozenset(g.loc[g[col], "mag"])
        if mem and mem not in seen:
            out.append((col, lab)); seen.append(mem)
    return out


# ---- 2. guild abundance vs methane ---------------------------------------------------------------
def guild_correlations(g):
    rng = np.random.default_rng(0)
    rows = []
    for col, label in active_guilds(g):
        members = list(g.loc[g[col], "mag"])
        for basis in AB:
            for subset, cores in (("12 cores", PRE), ("6 unrestored", UNRESTORED)):
                y = CH4[basis][cores].values
                x = AB[basis].loc[members, cores].sum().values
                rho, p = stats.spearmanr(x, y)
                null = np.array([stats.spearmanr(AB[basis].loc[list(rng.choice(MAGS, len(members), replace=False)),
                                                                cores].sum().values, y)[0] for _ in range(N_NULL)])
                null = null[np.isfinite(null)]
                raw = AB["raw"].loc[members, cores].sum().values   # detection judged on rel. abundance
                rows.append({"guild": col, "label": label, "n_species": len(members),
                             "n_mags": int(g.loc[g[col], f"{col}_n_mags" if f"{col}_n_mags" in g else "n_mags_in_species"].sum()),
                             "n_mags_marker": f"{col}_n_mags" in g, "basis": basis,
                             "cores": subset, "n_cores": len(cores), "spearman_rho": rho, "p": p,
                             "null_frac_ge_rho": (null >= rho).mean() if np.isfinite(rho) else np.nan,
                             "mean_rel_abund": raw.mean(), "max_rel_abund": raw.max(),
                             "cores_detected": int((raw > ZERO).sum())})
    out = pd.DataFrame(rows)
    out.to_csv(HERE / f"guild_methane_correlations_{TAG}.csv", index=False)
    mlb = g.loc[g.methanogen, "mag"].tolist()
    per_core = pd.DataFrame({"core": PRE, "restoration": [restoration(c) for c in PRE],
                             "CH4_umol_m2_d": CH4["raw"][PRE].values,
                             "methanogen_rel_abund": AB["raw"].loc[mlb, PRE].sum().values,
                             "archaea_rel_abund": AB["raw"].loc[list(g.loc[g.archaea, "mag"]), PRE].sum().values,
                             "unbinned_reads": merge_depths(rel.loc[["unbinned_reads"]]).loc["unbinned_reads", PRE].values})
    per_core.to_csv(HERE / f"methanogen_archaea_by_core_{TAG}.csv", index=False)
    return out


def per_mag_correlations(g):
    rows = []
    for mag in MAGS:
        r = {"mag": mag, f"iterativeID_{TAG}": ids[mag]}
        for basis, tag in (("raw", "raw"), ("DNA-normalised", "dna")):
            x = AB[basis].loc[mag, PRE].values
            rho, p = stats.spearmanr(x, CH4[basis][PRE].values) if np.ptp(x) > 0 else (np.nan, np.nan)
            r[f"rho_{tag}"], r[f"p_{tag}"] = rho, p
        rows.append(r)
    m = pd.DataFrame(rows)
    for tag in ("raw", "dna"):
        ok = m[f"p_{tag}"].notna()
        m.loc[ok, f"q_{tag}"] = multipletests(m.loc[ok, f"p_{tag}"], method="fdr_bh")[1]
    m = m.merge(g[["mag"] + [c for c, _ in GUILDS] + ["mean_rel_abund_12cores"]], on="mag")
    m.sort_values("rho_raw", ascending=False).to_csv(HERE / f"mag_methane_correlations_{TAG}.csv", index=False)
    print(f"per-species: {int((m.q_raw < 0.05).sum())} of {len(m)} correlated with CH4 at q < 0.05 (areal), "
          f"{int((m.q_dna < 0.05).sum())} (DNA-normalised)")
    return m


def plot_guilds(g, corr, basis):
    panels = [c for c, _ in active_guilds(g) if c != "methanosarcinales"]
    ncol = 3 if len(panels) > 4 else 2
    nrow = int(np.ceil(len(panels) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 4.75 * nrow), dpi=200, squeeze=False)
    unit = (r"CH$_4$ flux [$\mu$mol m$^{-2}$ d$^{-1}$, symlog; <0 = net uptake]" if basis == "raw" else
            r"CH$_4$ flux / DNA [symlog; <0 = net uptake]")
    for ax in axes.flat[len(panels):]:
        ax.set_visible(False)
    for ax, col in zip(axes.flat, panels):
        label = dict(GUILDS)[col]
        members = list(g.loc[g[col], "mag"])
        y = AB[basis].loc[members, PRE].sum() * (100 if basis == "raw" else 1)
        floor = (y[y > 0].min() / 3) if (y > 0).any() else 1e-6
        for cat, (colr, mk) in REST.items():
            for c in [c for c in PRE if restoration(c) == cat]:
                det = AB["raw"].loc[members, c].sum() > ZERO
                ax.scatter(CH4[basis][c], y[c] if y[c] > 0 else floor, s=70, marker=mk,
                           color=colr if det else "white", edgecolor=colr if not det else "#1f1f1f",
                           linewidths=1.6 if not det else 0.6, zorder=3)
        ax.set_xscale("symlog", linthresh=10 if basis == "raw" else 1); ax.set_yscale("log")
        ax.axvline(0, color="#9a9a9a", lw=0.8, zorder=1)
        sub = corr[(corr.guild == col) & (corr.basis == basis)].set_index("cores")
        a, b = sub.loc["12 cores"], sub.loc["6 unrestored"]
        txt = (f"12 cores: ρ = {a.spearman_rho:+.2f} (p = {a.p:.3f})\n"
               f"6 unrestored: ρ = {b.spearman_rho:+.2f} (p = {b.p:.3f})\n"
               f"null: {a.null_frac_ge_rho:.0%} of random {len(members)}-species sets ≥ ρ")
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", fontsize=8.5, color="#1f1f1f",
                bbox=dict(boxstyle="round", fc="white", ec="#9a9a9a", alpha=0.9))
        nm = int(a.n_mags)
        what = "with the marker" if a.n_mags_marker else "in these species"
        ax.set_title(f"{label}\n(n = {len(members)} species; {nm} MAG{'s' if nm != 1 else ''} {what})",
                     fontsize=10.5, loc="left", color="#1f1f1f")
        ax.set_xlabel(unit, fontsize=9, color="#3a3a3a")
        ax.set_ylabel("summed rel. abundance [% of reads]" if basis == "raw" else
                      "summed rel. abundance / DNA", fontsize=9, color="#3a3a3a")
        ax.grid(alpha=0.25, which="both", ls="--"); ax.tick_params(labelsize=8, colors="#3a3a3a")
    handles = [plt.Line2D([], [], ls="", marker=mk, ms=8, color=c, mec="#1f1f1f", label=cat)
               for cat, (c, mk) in REST.items()]
    handles.append(plt.Line2D([], [], ls="", marker="o", ms=8, color="white", mec="#5a5a5a", mew=1.6,
                              label="below detection (rel. abund. ≤ 1e-5; zeros at axis floor)"))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=10)
    fig.suptitle(f"Guild abundance vs. methane, final MAG set (585 species; {basis} basis)", fontsize=13, x=0.01,
                 ha="left", color="#1f1f1f")
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    out = HERE / f"guild_methane_{'raw' if basis == 'raw' else 'dnanorm'}_{TAG}.png"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print("wrote", out.name)


# ---- 3. co-occurrence -------------------------------------------------------------------------------
def collapse_to_taxa(df_cols_mag):
    """Columns MAG -> ID -> taxon (suffix dropped), summed; unbinned dropped."""
    df = df_cols_mag.loc[:, [c for c in df_cols_mag.columns if c in ids]]
    df.columns = [ids[c].rsplit(".", 1)[0] for c in df.columns]
    return df.T.groupby(df.columns).sum().T


def taxon_guilds(g):
    tg = collections.defaultdict(set)
    for _, r in g.iterrows():
        t = r[f"iterativeID_{TAG}"].rsplit(".", 1)[0]
        for col in ("mepn_consumer", "betaine_reducer", "methanogen"):
            if r[col]:
                tg[col].add(t)
    return tg


def taxon_phylum_kingdom():
    out = {}
    for mag, i in ids.items():
        t = taxa_of[mag]
        out.setdefault(i.rsplit(".", 1)[0], (t["Phylum"], t["Kingdom"]))
    return out


def heatmap(g):
    """4-site Spearman heatmap, as the 'all v all heatmap' cell (reduced=True, defined=True)."""
    import seaborn as sns
    from matplotlib import patches
    ab = merge_depths(rel.loc[MAGS])
    sites = sorted({p.split("_")[0] for p in PRE})
    ab = pd.DataFrame({s: (ab[f"{s}_A"] + ab[f"{s}_B"] + ab[f"{s}_C"]) / 3 for s in sites}).T
    ab = ab.loc[:, (ab.fillna(0) > ZERO).sum() >= 3]
    df = collapse_to_taxa(ab)
    n = df.shape[1]
    corr = pd.DataFrame(np.ones((n, n)), index=df.columns, columns=df.columns)
    pval = corr.copy()
    for i in range(n):
        for j in range(i + 1, n):
            c, p = stats.spearmanr(df.iloc[:, i], df.iloc[:, j])
            corr.iloc[i, j] = corr.iloc[j, i] = c
            pval.iloc[i, j] = pval.iloc[j, i] = p
    pres = (df > ZERO).astype(int)
    cooc = collections.Counter(pair for row in pres.itertuples(index=False)
                               for pair in combinations(sorted(c for c, v in zip(pres.columns, row) if v), 2))
    pairs = [(a, b, corr.loc[a, b], pval.loc[a, b]) for (a, b) in cooc if np.isfinite(corr.loc[a, b])]
    reject = multipletests([p for *_, p in pairs], alpha=0.05, method="fdr_bh")[0]
    json.dump(sorted({x for (a, b, *_), r in zip(pairs, reject) if r for x in (a, b)}),
              open(HERE / f"significantly_connected_organisms_{TAG}.json", "w"))
    corr.to_csv(HERE / f"heatmap_4site_spearman_{TAG}.csv")
    pk = taxon_phylum_kingdom()
    arch = sorted({pk[t][0] for t in df.columns if pk[t][1] == "Archaea"})
    bact = sorted({pk[t][0] for t in df.columns if pk[t][1] != "Archaea"})
    cmap = {p: plt.cm.turbo(i / max(len(arch), 1) * 0.15) for i, p in enumerate(arch)}
    cmap |= {p: plt.cm.turbo(0.2 + i / max(len(bact), 1) * 0.8) for i, p in enumerate(bact)}
    colors = pd.Series({t: cmap[pk[t][0]] for t in df.columns}, name="Phylum")
    scale = max(1.0, n / 109)                       # the 2020/R232 canvas (90 x 100 in) held 109 taxa
    cm = sns.clustermap(corr, row_colors=colors, col_colors=colors, cmap="coolwarm_r", center=0,
                        figsize=(90 * scale, 100 * scale), dendrogram_ratio=(0.1, 0.2))
    cm.figure.subplots_adjust(bottom=0.15, top=0.95)
    cm.ax_row_dendrogram.set_visible(False); cm.ax_col_dendrogram.set_visible(False)
    hm = cm.ax_heatmap.get_position()
    fw, fh = cm.figure.get_figwidth(), cm.figure.get_figheight()
    cm.ax_row_colors.set_position([hm.x0 - 0.015, hm.y0, 0.015, hm.height])
    cm.ax_col_colors.set_position([hm.x0, hm.y0 - 0.015, hm.width, 0.015])
    cm.ax_heatmap.yaxis.set_ticks_position("left"); cm.ax_heatmap.yaxis.set_label_position("left")
    cm.ax_heatmap.tick_params(axis="y", pad=0.015 * fw * 72 + 15)
    cm.ax_heatmap.tick_params(axis="x", pad=0.015 * fh * 72 + 15)
    cm.ax_cbar.set_ylabel("")
    cm.ax_cbar.tick_params(labelsize=45 * scale)
    cm.ax_cbar.set_xlabel("Spearman Correlation", fontsize=45 * scale, labelpad=30)
    LS = 40
    cm.ax_heatmap.set_yticks(np.arange(n) + 0.5)
    cm.ax_heatmap.set_xticks(np.arange(n) + 0.5)
    rows_o, cols_o = cm.dendrogram_row.reordered_ind, cm.dendrogram_col.reordered_ind
    cm.ax_heatmap.set_yticklabels([corr.index[i] for i in rows_o], fontsize=LS, rotation=0)
    cm.ax_heatmap.set_xticklabels([corr.columns[i] for i in cols_o], fontsize=LS, rotation=60, ha="right",
                                  rotation_mode="anchor")
    tg = taxon_guilds(g)
    gcol = {"methanogen": "#009E73", "betaine_reducer": "#1f77b4", "mepn_consumer": "#e8671b"}
    levels = json.load(open(HERE / f"iterativeID_levels_{TAG}.json"))
    genus_taxa = {i.rsplit(".", 1)[0] for i, lv in levels.items() if lv == "Genus"}
    boxed = set()
    for lab in list(cm.ax_heatmap.get_yticklabels()) + list(cm.ax_heatmap.get_xticklabels()):
        t = lab.get_text()
        for col in ("methanogen", "betaine_reducer", "mepn_consumer"):
            if t in tg[col]:
                lab.set_color(gcol[col]); lab.set_fontweight("bold"); lab.set_fontsize(LS * 1.2); boxed.add(t)
                break
        if t in genus_taxa:
            lab.set_fontstyle("italic")
    reo = corr.iloc[rows_o, cols_o]
    mask = np.triu(np.ones_like(reo, dtype=bool), k=1)
    mesh = cm.ax_heatmap.collections[0]
    arr = mesh.get_array().reshape(reo.shape); arr[mask] = np.nan; mesh.set_array(arr.ravel())
    cm.ax_cbar.set_position([hm.x1 - 0.1, hm.y0 + 0.3, 0.06, hm.height / 5])
    N = len(corr.index)
    for t in boxed:
        rp = rows_o.index(corr.index.get_loc(t)); cp = cols_o.index(corr.columns.get_loc(t))
        cm.ax_heatmap.add_patch(patches.Rectangle((0, rp), rp + 1, 1, lw=6, ec="black", fc="none", clip_on=False))
        cm.ax_heatmap.add_patch(patches.Rectangle((cp, N), 1, -(N - cp), lw=6, ec="black", fc="none", clip_on=False))
    preo = pval.iloc[rows_o, cols_o]
    for i in range(N):
        for j in range(N):
            if preo.iloc[i, j] < 0.05 and not mask[i, j] and i != j:
                cm.ax_heatmap.add_patch(patches.Rectangle((j, i), 1, 1, fill=False, ec="lightgreen", lw=8))
    cm.ax_heatmap.set_xlabel("Member taxa (GTDB R232 iterativeIDs, species summed per taxon)", fontsize=50 * scale)
    cm.ax_heatmap.set_ylabel("Member taxa", fontsize=50 * scale)
    head = lambda s: patches.Patch(color="none", label=rf"$\bf{{{s.replace(' ', chr(92) + ' ')}}}$")
    guild_lab = {"mepn_consumer": "MePn consumer (phnJ)", "betaine_reducer": "Betaine reducer (grdH/grdI)",
                 "methanogen": "Methanogen (mcrA)"}
    shown = [col for col in guild_lab if tg[col] & set(df.columns)]
    handles = ([head("Archaea")] + [patches.Patch(color=cmap[p], label=p) for p in arch] if arch else []) + \
              [head("Bacteria")] + [patches.Patch(color=cmap[p], label=p) for p in bact] + \
              [head("Guild (label colour)")] + [patches.Patch(color=gcol[c], label=guild_lab[c]) for c in shown]
    absent = {col: sorted(tg[col] - set(df.columns)) for col in guild_lab}
    print("guild taxa absent from the heatmap (present at < 3 of 4 sites):", absent)
    cm.ax_heatmap.legend(handles=handles, title="Taxonomic Phylum", title_fontsize=8 * fw / 10, loc="lower left",
                         bbox_to_anchor=(0.65, 0.5), fontsize=7 * fw / 10, frameon=True)
    out = HERE / f"abundance_correlatons_one_triangle_reduced_defined_{TAG}.png"
    cm.figure.savefig(out, dpi=150 if scale > 1.5 else 300, bbox_inches="tight"); plt.close(cm.figure)
    p_floor = pval.values[np.triu_indices(N, 1)]
    print(f"heatmap: {N} taxa over {len(sites)} sites; {int(reject.sum())} FDR-significant pairs; "
          f"pairs with nominal p<0.05: {int((p_floor < 0.05).sum())} (n=4 -> only |rho|=1 reaches it)")
    print("wrote", out.name)
    return {"heatmap_taxa": N, "heatmap_archaeal_taxa": sum(pk[t][1] == "Archaea" for t in df.columns),
            "heatmap_guild_taxa_absent": absent}


def within_condition(EFFECT_THRESHOLD=0.7, FDR_ALPHA=0.05, MIN_COOCC=6):
    """Port of within_condition_correlation() in analysis_2_marimo.py, same parameters."""
    from skbio.stats.composition import clr
    from sklearn.covariance import GraphicalLassoCV
    ab = merge_depths(rel.loc[MAGS])
    sites = sorted({p.split("_")[0] for p in PRE})
    rows, cond = {}, []
    for site in sites:
        for rep in ("A", "B", "C"):
            rows[f"{site}_{rep}"] = ab[f"{site}_{rep}"]; cond.append(site)
    long = pd.DataFrame(rows).T
    cond = np.array(cond)
    long = long.loc[:, (long.fillna(0) > ZERO).sum() >= MIN_COOCC]
    n_mags_kept = long.shape[1]
    long = collapse_to_taxa(long)
    vals = long.values
    pseudo = vals[vals > 0].min() / 2
    X = clr(vals + pseudo)
    taxa = long.columns.to_numpy()
    V = X.shape[1]
    C_marg = np.corrcoef(X.T)
    R = X.copy()
    for c in np.unique(cond):
        R[cond == c] -= R[cond == c].mean(axis=0, keepdims=True)
    C_within = np.corrcoef(R.T)
    C_betw, _ = stats.spearmanr(np.vstack([X[cond == c].mean(0) for c in np.unique(cond)]))
    Cm, Cw, Cb = (pd.DataFrame(a, index=taxa, columns=taxa) for a in (C_marg, C_within, C_betw))
    pres = (long > ZERO).astype(int)
    cooc = collections.Counter(pair for row in pres.itertuples(index=False)
                               for pair in combinations(sorted(c for c, v in zip(pres.columns, row) if v), 2))
    idxs = [np.where(cond == c)[0] for c in np.unique(cond)]
    PERM = []
    for combo in product(*[list(permutations(ix)) for ix in idxs]):
        pi = np.arange(len(cond))
        for orig, new in zip(idxs, combo):
            pi[orig] = list(new)
        PERM.append(pi)
    PERM = np.array(PERM)
    tl = {t: i for i, t in enumerate(taxa)}
    cmean = lambda v: v - np.array([v[cond == c].mean() for c in cond])

    def perm_p(a, b, rw):
        xr, yr = cmean(X[:, tl[a]]), cmean(X[:, tl[b]])
        d = np.sqrt((xr ** 2).sum() * (yr ** 2).sum())
        null = (np.abs(yr[PERM] @ xr) / d if d > 0 else np.zeros(len(PERM))) + 1e-12
        return (np.sum(null >= abs(rw)) + 1) / (len(null) + 1)

    cand = []
    for (a, b), cnt in cooc.items():
        if cnt < MIN_COOCC:
            continue
        rw = Cw.loc[a, b]
        if np.isnan(rw):
            continue
        cand.append((a, b, rw, Cm.loc[a, b], Cb.loc[a, b], cnt, perm_p(a, b, rw)))
    cdf = pd.DataFrame(cand, columns=["taxon_a", "taxon_b", "r_within", "r_marginal", "r_between",
                                      "cooccurrence", "perm_p"])
    rej, q, _, _ = multipletests(cdf.perm_p.values, alpha=FDR_ALPHA, method="fdr_bh")
    cdf["q_BH"], cdf["fdr_pass"] = q, rej
    cdf["effect_pass"] = cdf.r_within.abs() >= EFFECT_THRESHOLD
    cdf["both_pass"] = cdf.effect_pass & cdf.fdr_pass
    cdf.to_csv(HERE / f"all_candidate_pairs_{TAG}.csv", index=False)
    cdf[cdf.both_pass].sort_values(["q_BH", "perm_p"]).to_csv(HERE / f"direct_interactions_{TAG}.csv", index=False)
    rng = np.random.default_rng(42)
    iu = np.triu_indices(V, 1)
    maxT, n_ge = [], 0
    for _ in range(500):
        Rs = R.copy()
        for col in range(V):
            for ix in idxs:
                Rs[ix, col] = Rs[rng.permutation(ix), col]
        maxT.append(np.abs(np.corrcoef(Rs.T)[iu]).max())
    maxT = np.array(maxT)
    q95 = np.percentile(maxT, 95)
    gl = GraphicalLassoCV(alphas=10, max_iter=200).fit(R)
    P = gl.precision_
    edges = [(a, b, -P[i, j] / np.sqrt(P[i, i] * P[j, j])) for i, a in enumerate(taxa) for j, b in enumerate(taxa)
             if j > i and abs(P[i, j]) > 1e-8]
    gdf = pd.DataFrame(edges, columns=["taxon_a", "taxon_b", "partial_r"]).sort_values(
        "partial_r", key=lambda s: s.abs(), ascending=False)
    gdf.to_csv(HERE / f"glasso_direct_edges_{TAG}.csv", index=False)
    for name, m in (("C_marginal", Cm), ("C_within", Cw), ("C_between", Cb)):
        m.to_csv(HERE / f"{name}_{TAG}.csv")
    json.dump({"condition": cond.tolist(), "samples": list(long.index), "taxa": list(taxa)},
              open(HERE / f"within_condition_meta_{TAG}.json", "w"))
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(C_marg[iu], C_within[iu], s=2, alpha=0.3, color="#2a78d6")
    ax.plot([-1, 1], [-1, 1], "k--", lw=1)
    for s in (1, -1):
        ax.axhline(s * EFFECT_THRESHOLD, color="#9a9a9a", ls=":", lw=0.8)
    ax.set_xlabel("marginal r (12 samples, Pearson on CLR)"); ax.set_ylabel("within-condition r (residuals)")
    ax.set_title("Coupling vs. condition response (final set, 585 species)"); ax.set_aspect("equal")
    fig.tight_layout(); fig.savefig(HERE / f"marginal_vs_within_{TAG}.png", dpi=200); plt.close(fig)
    arch = sum(taxon_phylum_kingdom()[t][1] == "Archaea" for t in taxa)
    w = Cw.values[iu]
    print(f"within-condition: params EFFECT={EFFECT_THRESHOLD} FDR={FDR_ALPHA} MIN_COOCC={MIN_COOCC} ZERO={ZERO} "
          f"pseudocount={pseudo:.3e} GraphicalLassoCV(alphas=10,max_iter=200) alpha_={gl.alpha_:.4g}")
    print(f"  species kept {n_mags_kept}; taxa {V} ({arch} archaeal), candidates {len(cdf)}, effect-pass "
          f"{int(cdf.effect_pass.sum())}, FDR-pass {int(cdf.fdr_pass.sum())}, both {int(cdf.both_pass.sum())}; "
          f"maxT 95th pct |r| = {q95:.3f} (pairs above: {int((np.abs(w) >= q95).sum())}); Glasso edges {len(gdf)}")
    return {"within_species_kept": n_mags_kept, "within_taxa": int(V), "within_archaeal_taxa": int(arch),
            "candidates": len(cdf), "effect_pass": int(cdf.effect_pass.sum()), "fdr_pass": int(cdf.fdr_pass.sum()),
            "maxT_q95": float(q95), "pairs_ge_maxT_q95": int((np.abs(w) >= q95).sum()),
            "glasso_alpha": float(gl.alpha_), "glasso_edges_all": len(gdf),
            "glasso_edges_gt005": int((gdf.partial_r.abs() > 0.05).sum()), "pseudocount": float(pseudo)}


def network(g):
    abf = HERE / f"averaged_normalized_MAG_abundances_{TAG}.csv"
    AB["raw"].loc[MAGS].to_csv(abf)
    cmd = [sys.executable, "-I", str(COR / "louvain_cooccurrence_network.py"),
           "--edges", str(HERE / f"glasso_direct_edges_{TAG}.csv"), "--meta", str(HERE / f"within_condition_meta_{TAG}.json"),
           "--ids", str(HERE / f"iterativeIDs_{TAG}.json"), "--phylogeny", str(HERE / f"phylogeny_{TAG}.json"),
           "--abundances", str(abf),
           "--out-png", str(HERE / f"within_condition_louvain_network_{TAG}.png"),
           "--out-modules", str(HERE / f"within_condition_louvain_modules_{TAG}.json"),
           "--subtitle", "final MAG set (585 species), GTDB R232 labels; glasso_direct_edges_final585.csv"]
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    print(res.stdout)
    mods = json.load(open(HERE / f"within_condition_louvain_modules_{TAG}.json"))
    Q = float(next(l for l in res.stdout.splitlines() if l.startswith("Louvain:")).split("Q = ")[1])
    gdf = pd.read_csv(HERE / f"glasso_direct_edges_{TAG}.csv")
    gdf = gdf[gdf.partial_r.abs() > 0.05]
    nodes = set(gdf.taxon_a) | set(gdf.taxon_b)
    pk = taxon_phylum_kingdom()
    tg = taxon_guilds(g)
    arch_nodes = {t for t in nodes if pk[t][1] == "Archaea"}
    mod_of = {m: k for k, v in mods.items() for m in v["members"]}
    same_phylum = int(sum(pk[a][0] == pk[b][0] for a, b in zip(gdf.taxon_a, gdf.taxon_b)))
    meta_taxa = set(json.load(open(HERE / f"within_condition_meta_{TAG}.json"))["taxa"])
    deg = collections.Counter(list(gdf.taxon_a) + list(gdf.taxon_b))
    return {"network_nodes": len(nodes), "network_edges": len(gdf), "louvain_Q": Q,
            "louvain_modules": sum(v["type"] == "louvain" for v in mods.values()),
            "module_sizes": sorted([v["size"] for v in mods.values() if v["type"] == "louvain"], reverse=True),
            "network_archaeal_taxa": len(arch_nodes),
            "archaeal_modules": sorted({mod_of[t] for t in arch_nodes}),
            "same_phylum_edges": same_phylum,
            "guild_taxa_in_network": {c: sorted(tg[c] & nodes) for c in tg},
            "guild_taxa_in_panel_without_edges": {c: sorted((tg[c] & meta_taxa) - nodes) for c in tg},
            "guild_taxa_not_in_panel": {c: sorted(tg[c] - meta_taxa) for c in tg},
            "guild_taxa_degree": {c: {t: deg[t] for t in sorted(tg[c] & nodes)} for c in tg}}


def table3():
    """Family share (%) of binned abundance per site type (GTDB R232 families of the 585 representatives)."""
    types = {"Unrestored (R1, R2)": ("R1_", "R2_"), "Restored (SF2)": ("SF2",), "Reference (R2A)": ("R2A",)}
    fam = [taxa_of[m]["Family"] or "no family assigned" for m in MAGS]
    cols, unb = {}, {}
    for t, pre in types.items():
        cs = [c for c in rel.columns if c.startswith(pre)]
        v = rel.loc[MAGS, cs].mean(axis=1)
        cols[t] = (v / v.sum()).groupby(fam).sum() * 100
        unb[t] = rel.loc["unbinned_reads", cs].mean() * 100
    tab = pd.DataFrame(cols).sort_values("Unrestored (R1, R2)", ascending=False)
    top, other = tab.head(10), tab.iloc[10:].sum()
    arch = {taxa_of[m]["Family"] for m in MAGS if taxa_of[m]["Kingdom"] == "Archaea"}
    rows = [[f, *map(float, r)] for f, r in zip(top.index, top.values)] + [["Other families", *map(float, other.values)]]
    out = {"columns": list(types), "rows": rows, "unbinned": unb, "archaeal_in_top": sorted(set(top.index) & arch)}
    json.dump(out, open(HERE / f"table3_{TAG}.json", "w"), indent=1)
    tab.to_csv(HERE / f"table3_all_families_{TAG}.csv")
    print(pd.DataFrame([r[1:] for r in rows], index=[r[0] for r in rows], columns=list(types)).round(1).to_string())
    print("unbinned %:", {k: round(v, 1) for k, v in unb.items()})
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parts", default="guilds,corr,heatmap,network,table3")
    parts = set(ap.parse_args().parts.split(","))
    summ_f = HERE / f"summary_{TAG}.json"
    summ = json.load(open(summ_f)) if summ_f.exists() else {}
    g = build_guilds() if parts & {"guilds", "corr", "heatmap", "network"} else None
    if g is not None:
        summ["guild_species"] = {c: int(g[c].sum()) for c in GCALL + ["archaea", "methanosarcinales"]}
        summ["guild_mags"] = {c: int(g.loc[g[c], f"{c}_n_mags"].sum()) for c in GCALL}
        summ["n_species"], summ["n_archaeal_species"] = len(g), int(g.archaea.sum())
    if "corr" in parts:
        corr = guild_correlations(g)
        m = per_mag_correlations(g)
        summ["per_species_q05_raw"], summ["per_species_q05_dna"] = int((m.q_raw < 0.05).sum()), int((m.q_dna < 0.05).sum())
        for basis in AB:
            plot_guilds(g, corr, basis)
        show = corr[corr.basis == "raw"][["label", "n_species", "cores", "spearman_rho", "p", "null_frac_ge_rho",
                                           "mean_rel_abund", "cores_detected"]]
        print(show.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    if "heatmap" in parts:
        summ |= heatmap(g)
    if "network" in parts:
        summ |= within_condition()
        summ |= network(g)
    if "table3" in parts:
        table3()
    json.dump(summ, open(summ_f, "w"), indent=1, default=str)
