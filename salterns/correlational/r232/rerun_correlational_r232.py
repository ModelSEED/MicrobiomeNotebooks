#!/usr/bin/env python
"""Rerun of the salterns correlational studies with the GTDB r232 labels.

Changes vs the 2020 analyses in ../analysis_2_marimo.py:
  * labels  - ../iterativeIDs_r232.json (GTDB-Tk 2.8.0 / r232) for all 310 dereplicated MAGs. The
              30 archaea now have IDs, so they enter the co-occurrence analyses (the 2020 "defined
              IDs" filter dropped every MAG without an ID, i.e. all archaea).
  * guilds  - re-derived per MAG from the RAST gene annotation with explicit rules (guilds_r232.csv),
              kept next to the hand-curated 2020 sets:
                MePn consumer    phnJ, the C-P lyase (EC 4.7.1.1)   [2020: any phn gene -> 27 MAGs]
                betaine reducer  grdH, betaine reductase B           [2020: 15 model-curated MAGs]
                TMA candidates   cntA (carnitine monooxygenase) or cutC-cluster genes; no MAG has CutC
                                 itself annotated, only genes "clustered with choline TMA-lyase"
                methanogen       mcrA (methyl-coenzyme M reductase)
                context          all archaea; Methanosarcinales
  * guild abundance vs CH4 on both bases the repo uses (raw areal flux; DNA-normalised), all 12 cores
    and the 6 unrestored cores, with the random-subset specificity null of
    ../../clauded_insights/reassessment_figure.py (subsets drawn from the 310 MAGs, not unbinned reads).
  * co-occurrence - the 4-site Spearman heatmap and the within-condition (4 conditions x 3 cores)
    pipeline + Graphical Lasso, recomputed with the 2020 parameters; the Louvain figure comes from
    ../louvain_cooccurrence_network.py pointed at the r232 files.
Not rerun: the FBA simulations (fits / fluxes) - modelling, not correlation.
Every output goes to this folder; no 2020 file is modified.

Run (from this folder):  ~/Documents/py_venv/bin/python rerun_correlational_r232.py [--parts guilds,corr,heatmap,network]
"""
import argparse
import collections
import json
import re
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
RAST = S / "jplfaria/Cliff MAGs/RASTified_MAG_json"
PAT = re.compile(r"Salt_Pond_MetaG_?(.+?)_MG_DASTool_bins[._](concoct_out|metabat|maxbin)\.(\d+)")
ZERO = 1e-5                     # detection level used throughout the 2020 code
N_NULL = 4000                   # random-subset draws (as in reassessment_figure.py)
REST = {"unrestored": ("#D55E00", "o"), "restored": ("#009E73", "s"), "reference": ("#56B4E9", "^")}
GENES = {   # RAST function-string patterns (lower-case); evidence strings are kept in guilds_r232.csv
    "phnJ": r"c-p lyase|ec 4\.7\.1\.1\b",
    "phn_accessory": r"methylphosphonate 5-triphosphate (synthase|diphosphatase)|phnn|ec 2\.7\.4\.23",
    "grdH": r"betaine reductase component b",
    "cutC_cluster": r"choline trimethylamine.lyase",   # matches only "...clustered with choline TMA-lyase"
    "cntA": r"carnitine monooxygenase",
    "mcrA": r"methyl.coenzyme m reductase",
    "mtbB_mtmB": r"(dimethylamine|monomethylamine)\s*(:corrinoid )?methyltransferase",
}


def key(name):
    m = PAT.search(str(name))
    return f"{m[1]}|{m[2]}|{int(m[3])}" if m else None


def restoration(core):
    return "restored" if core.startswith("SF2") else "reference" if core.startswith("R2A") else "unrestored"


# ---- shared inputs ------------------------------------------------------------------------------
ids = json.load(open(COR / "iterativeIDs_r232.json"))                 # MAG -> r232 ID
ids_2020 = json.load(open(COR / "iterativeIDs.json"))
phylo = json.load(open(S / "data/Saltern_phylogeny_r232.json"))        # MAG.contigs__ -> ranks
rel = pd.read_csv(S / "data/Cliff_310MAG_relabund.tsv", sep="\t", index_col=0)
dna = pd.read_csv(S / "data/SaltPondsDNA.csv").set_index("Sample")["DNA conc (µg/g soil)"]
meta = pd.read_csv(S / "data/Cliff_Sample_Metadata_BGC_NMR.csv").set_index("Sample")
MAGS = [m for m in rel.index if key(m)]                                # 310; drops unbinned_reads
assert len(MAGS) == 310 and all(m in ids for m in MAGS), "relabund rows must all have r232 IDs"
taxa_of = {m.replace(".contigs__", ""): t for m, t in phylo.items()}
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


def marimo_set(name):
    src = (COR / "analysis_2_marimo.py").read_text()
    body = re.search(rf"\n    {name} = \{{(.*?)\}}", src, re.S)[1]
    return {key(x) for x in re.findall(r"'(Salt_Pond[^']+)'", body)}


# ---- 1. guilds ------------------------------------------------------------------------------------
def build_guilds():
    fn_block = re.compile(r'"functions"\s*:\s*\[(.*?)\]', re.S)
    evidence = {g: {} for g in GENES}
    for mag in MAGS:
        path = RAST / f"{mag}.contigs__.RAST.json"
        text = path.read_text(errors="ignore")
        funcs = {s for blk in fn_block.findall(text) for s in re.findall(r'"([^"]*)"', blk)}
        for g, pat in GENES.items():
            hits = sorted(s for s in funcs if re.search(pat, s.lower()))
            if hits:
                evidence[g][mag] = hits
    mepn_2020, bet_2020 = marimo_set("mepe"), marimo_set("betaine")
    checkm = pd.read_csv(S / "MAG_QC/checkm/summary/CheckM_summary_table.tsv", sep="\t")
    checkm.index = checkm["Bin Name"].map(key)
    rows = []
    for mag in MAGS:
        t, k = taxa_of[mag], key(mag)
        rows.append({
            "mag": mag, "iterativeID_r232": ids[mag], "iterativeID_2020": ids_2020.get(mag, ""),
            "domain": t["Kingdom"], "phylum": t["Phylum"], "order": t["Order"], "family": t["Family"],
            "genus": t["Genus"], "species": t["Species"],
            "completeness": checkm.loc[k, "Completeness"], "contamination": checkm.loc[k, "Contamination"],
            "mepn_consumer": mag in evidence["phnJ"], "mepn_2020": k in mepn_2020,
            "betaine_reducer": mag in evidence["grdH"], "betaine_2020": k in bet_2020,
            "tma_producer": mag in evidence["cutC_cluster"] or mag in evidence["cntA"],
            "methanogen": mag in evidence["mcrA"], "archaea": t["Kingdom"] == "Archaea",
            "methanosarcinales": t["Order"] == "Methanosarcinales",
            **{f"ev_{g}": " | ".join(evidence[g].get(mag, [])) for g in GENES},
            "mean_rel_abund_12cores": AB["raw"].loc[mag].mean(),
        })
    g = pd.DataFrame(rows)
    g.to_csv(HERE / "guilds_r232.csv", index=False)
    print("guild sizes:", {c: int(g[c].sum()) for c in ("mepn_consumer", "mepn_2020", "betaine_reducer",
          "betaine_2020", "tma_producer", "methanogen", "archaea", "methanosarcinales")})
    for a, b in (("mepn_2020", "mepn_consumer"), ("betaine_2020", "betaine_reducer")):
        print(f"  {a} not in {b}:", sorted(g.loc[g[a] & ~g[b], "iterativeID_r232"]))
        print(f"  {b} not in {a}:", sorted(g.loc[g[b] & ~g[a], "iterativeID_r232"]))
    print("TMA-candidate evidence (RAST strings):")
    for _, r in g[g.tma_producer].iterrows():
        print(f"  {r.iterativeID_r232}: {(r.ev_cutC_cluster or r.ev_cntA)[:150]}")
    return g


GUILDS = [("mepn_consumer", "MePn consumers (phnJ)"), ("mepn_2020", "MePn consumers, 2020 set"),
          ("betaine_reducer", "Betaine reducers (grdH)"), ("betaine_2020", "Betaine reducers, 2020 set"),
          ("tma_producer", "TMA-route candidates (cntA; cutC-cluster genes)"), ("methanogen", "Methanogens (mcrA)"),
          ("methanosarcinales", "Methanosarcinales"), ("archaea", "All archaea")]


# ---- 2. guild abundance vs methane ---------------------------------------------------------------
def guild_correlations(g):
    rng = np.random.default_rng(0)
    rows = []
    for col, label in GUILDS:
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
                rows.append({"guild": col, "label": label, "n_MAGs": len(members), "basis": basis,
                             "cores": subset, "n_cores": len(cores), "spearman_rho": rho, "p": p,
                             "null_frac_ge_rho": (null >= rho).mean() if np.isfinite(rho) else np.nan,
                             "mean_rel_abund": raw.mean(), "max_rel_abund": raw.max(),
                             "cores_detected": int((raw > ZERO).sum())})
    out = pd.DataFrame(rows)
    out.to_csv(HERE / "guild_methane_correlations_r232.csv", index=False)
    mlb = g.loc[g.methanogen, "mag"].tolist()
    per_core = pd.DataFrame({"core": PRE, "restoration": [restoration(c) for c in PRE],
                             "CH4_umol_m2_d": CH4["raw"][PRE].values,
                             "methanogen_rel_abund": AB["raw"].loc[mlb, PRE].sum().values,
                             "archaea_rel_abund": AB["raw"].loc[list(g.loc[g.archaea, "mag"]), PRE].sum().values})
    per_core.to_csv(HERE / "methanogen_archaea_by_core_r232.csv", index=False)
    return out


def per_mag_correlations(g):
    rows = []
    for mag in MAGS:
        r = {"mag": mag, "iterativeID_r232": ids[mag]}
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
    m.sort_values("rho_raw", ascending=False).to_csv(HERE / "mag_methane_correlations_r232.csv", index=False)
    return m


def plot_guilds(g, corr, basis):
    panels = [c for c, _ in GUILDS if c not in ("tma_producer", "methanosarcinales")]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9.5), dpi=200)
    unit = (r"CH$_4$ flux [$\mu$mol m$^{-2}$ d$^{-1}$, symlog; <0 = net uptake]" if basis == "raw" else
            r"CH$_4$ flux / DNA [symlog; <0 = net uptake]")
    for ax, col in zip(axes.flat, panels):
        label = dict(GUILDS)[col]
        members = list(g.loc[g[col], "mag"])
        y = AB[basis].loc[members, PRE].sum() * (100 if basis == "raw" else 1)
        floor = (y[y > 0].min() / 3) if (y > 0).any() else 1e-6
        for cat, (colr, mk) in REST.items():
            cores = [c for c in PRE if restoration(c) == cat]
            for c in cores:
                det = AB["raw"].loc[members, c].sum() > ZERO
                ax.scatter(CH4[basis][c], y[c] if y[c] > 0 else floor, s=70, marker=mk,
                           color=colr if det else "white", edgecolor=colr if not det else "#1f1f1f",
                           linewidths=1.6 if not det else 0.6, zorder=3)
        # reference cores (R2A) have net CH4 uptake (negative flux) -> symmetric-log x axis
        ax.set_xscale("symlog", linthresh=10 if basis == "raw" else 1); ax.set_yscale("log")
        ax.axvline(0, color="#9a9a9a", lw=0.8, zorder=1)
        sub = corr[(corr.guild == col) & (corr.basis == basis)].set_index("cores")
        a, b = sub.loc["12 cores"], sub.loc["6 unrestored"]
        txt = (f"12 cores: ρ = {a.spearman_rho:+.2f} (p = {a.p:.3f})\n"
               f"6 unrestored: ρ = {b.spearman_rho:+.2f} (p = {b.p:.3f})\n"
               f"null: {a.null_frac_ge_rho:.0%} of random {len(members)}-MAG sets ≥ ρ")
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", fontsize=8.5, color="#1f1f1f",
                bbox=dict(boxstyle="round", fc="white", ec="#9a9a9a", alpha=0.9))
        ax.set_title(f"{label}  (n = {len(members)} MAG{'s' if len(members) != 1 else ''})",
                     fontsize=11, loc="left", color="#1f1f1f")
        ax.set_xlabel(unit, fontsize=9, color="#3a3a3a")
        ax.set_ylabel("summed rel. abundance [% of mapped reads]" if basis == "raw" else
                      "summed rel. abundance / DNA", fontsize=9, color="#3a3a3a")
        ax.grid(alpha=0.25, which="both", ls="--"); ax.tick_params(labelsize=8, colors="#3a3a3a")
    handles = [plt.Line2D([], [], ls="", marker=mk, ms=8, color=c, mec="#1f1f1f", label=cat)
               for cat, (c, mk) in REST.items()]
    handles.append(plt.Line2D([], [], ls="", marker="o", ms=8, color="white", mec="#5a5a5a", mew=1.6,
                              label="below detection (rel. abund. ≤ 1e-5; zeros at axis floor)"))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=10)
    fig.suptitle(f"Guild abundance vs. methane, GTDB r232 guilds ({basis} basis)", fontsize=13, x=0.01,
                 ha="left", color="#1f1f1f")
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    out = HERE / f"guild_methane_{'raw' if basis == 'raw' else 'dnanorm'}_r232.png"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print("wrote", out.name)


# ---- 3. co-occurrence -------------------------------------------------------------------------------
def collapse_to_taxa(df_cols_mag):
    """Columns MAG -> r232 ID -> taxon (suffix dropped), summed; unbinned/undefined dropped."""
    df = df_cols_mag.loc[:, [c for c in df_cols_mag.columns if c in ids]]
    df.columns = [ids[c].rsplit(".", 1)[0] for c in df.columns]
    return df.T.groupby(df.columns).sum().T


def taxon_guilds(g):
    tg = collections.defaultdict(set)
    for _, r in g.iterrows():
        t = r.iterativeID_r232.rsplit(".", 1)[0]
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
    ab = merge_depths(rel)
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
              open(HERE / "significantly_connected_organisms_r232.json", "w"))
    pk = taxon_phylum_kingdom()
    arch = sorted({pk[t][0] for t in df.columns if pk[t][1] == "Archaea"})
    bact = sorted({pk[t][0] for t in df.columns if pk[t][1] != "Archaea"})
    cmap = {p: plt.cm.turbo(i / max(len(arch), 1) * 0.15) for i, p in enumerate(arch)}
    cmap |= {p: plt.cm.turbo(0.2 + i / max(len(bact), 1) * 0.8) for i, p in enumerate(bact)}
    colors = pd.Series({t: cmap[pk[t][0]] for t in df.columns}, name="Phylum")
    cm = sns.clustermap(corr, row_colors=colors, col_colors=colors, cmap="coolwarm_r", center=0,
                        figsize=(90, 100), dendrogram_ratio=(0.1, 0.2))
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
    cm.ax_cbar.tick_params(labelsize=45)
    cm.ax_cbar.set_xlabel("Spearman Correlation", fontsize=45, labelpad=30)
    LS = 40
    cm.ax_heatmap.set_yticklabels(cm.ax_heatmap.get_yticklabels(), fontsize=LS, rotation=0)
    cm.ax_heatmap.set_xticklabels(cm.ax_heatmap.get_xticklabels(), fontsize=LS, rotation=60, ha="right",
                                  rotation_mode="anchor")
    tg = taxon_guilds(g)
    gcol = {"methanogen": "#009E73", "betaine_reducer": "#1f77b4", "mepn_consumer": "#e8671b"}
    levels = json.load(open(COR / "iterativeID_levels_r232.json"))
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
    rows_o, cols_o = cm.dendrogram_row.reordered_ind, cm.dendrogram_col.reordered_ind
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
    cm.ax_heatmap.set_xlabel("Member taxa (GTDB r232 iterativeIDs, MAGs summed per taxon)", fontsize=50)
    cm.ax_heatmap.set_ylabel("Member taxa", fontsize=50)
    head = lambda s: patches.Patch(color="none", label=rf"$\bf{{{s.replace(' ', chr(92) + ' ')}}}$")
    guild_lab = {"mepn_consumer": "MePn consumer (phnJ)", "betaine_reducer": "Betaine reducer (grdH)",
                 "methanogen": "Methanogen (mcrA)"}
    shown = [col for col in guild_lab if tg[col] & set(df.columns)]   # only guilds that made the filter
    handles = ([head("Archaea")] + [patches.Patch(color=cmap[p], label=p) for p in arch] if arch else []) + \
              [head("Bacteria")] + [patches.Patch(color=cmap[p], label=p) for p in bact] + \
              [head("Guild (label colour)")] + [patches.Patch(color=gcol[c], label=guild_lab[c]) for c in shown]
    print("guild taxa absent from the heatmap (present at < 3 of 4 sites):",
          {col: sorted(tg[col] - set(df.columns)) for col in guild_lab})
    cm.ax_heatmap.legend(handles=handles, title="Taxonomic Phylum", title_fontsize=8 * fw / 10, loc="lower left",
                         bbox_to_anchor=(0.65, 0.5), fontsize=7 * fw / 10, frameon=True)
    out = HERE / "abundance_correlatons_one_triangle_reduced_defined_r232.png"
    cm.figure.savefig(out, dpi=300, bbox_inches="tight"); plt.close(cm.figure)
    p_floor = pval.values[np.triu_indices(N, 1)]
    print(f"heatmap: {N} taxa over {len(sites)} sites; {int(reject.sum())} FDR-significant pairs; "
          f"pairs with nominal p<0.05: {int((p_floor < 0.05).sum())} (n=4 -> only |rho|=1 reaches it)")
    print("wrote", out.name)


def within_condition(EFFECT_THRESHOLD=0.7, FDR_ALPHA=0.05, MIN_COOCC=6):
    """Port of within_condition_correlation() in analysis_2_marimo.py, r232 labels, same parameters."""
    from skbio.stats.composition import clr
    from sklearn.covariance import GraphicalLassoCV
    ab = merge_depths(rel)
    sites = sorted({p.split("_")[0] for p in PRE})
    rows, cond = {}, []
    for site in sites:
        for rep in ("A", "B", "C"):
            rows[f"{site}_{rep}"] = ab[f"{site}_{rep}"]; cond.append(site)
    long = pd.DataFrame(rows).T
    cond = np.array(cond)
    long = long.loc[:, (long.fillna(0) > ZERO).sum() >= MIN_COOCC]
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
    # exact within-condition permutation null (all 3!^4 = 1296 permutations), vectorised: permuting inside a
    # condition leaves condition means unchanged, so each permuted residual vector is a re-indexing of yr.
    idxs = [np.where(cond == c)[0] for c in np.unique(cond)]
    PERM = []
    for combo in product(*[list(permutations(ix)) for ix in idxs]):
        pi = np.arange(len(cond))
        for orig, new in zip(idxs, combo):
            pi[orig] = list(new)
        PERM.append(pi)
    PERM = np.array(PERM)
    tl = list(taxa)

    def perm_p(a, b, rw, loop=False):
        x, y = X[:, tl.index(a)], X[:, tl.index(b)]
        xr = x - np.array([x[cond == c].mean() for c in cond])
        yr = y - np.array([y[cond == c].mean() for c in cond])
        if loop:   # the original implementation, kept for the equality check below
            null = []
            for pi in PERM:
                yp = y[pi]; ypr = yp - np.array([yp[cond == c].mean() for c in cond])
                d = np.sqrt((xr ** 2).sum() * (ypr ** 2).sum())
                null.append(abs((xr * ypr).sum() / d) if d > 0 else 0.0)
            null = np.array(null)
        else:
            d = np.sqrt((xr ** 2).sum() * (yr ** 2).sum())
            null = np.abs(yr[PERM] @ xr) / d if d > 0 else np.zeros(len(PERM))
            null = null + 1e-12   # the identity permutation equals |r_within| only up to ~4e-16 rounding;
                                  # count it (and exact ties) deterministically instead of by rounding luck
        return (np.sum(null >= abs(rw)) + 1) / (len(null) + 1)

    cand = []
    for (a, b), cnt in cooc.items():
        if cnt < MIN_COOCC:
            continue
        rw = Cw.loc[a, b]
        if np.isnan(rw):
            continue
        cand.append((a, b, rw, Cm.loc[a, b], Cb.loc[a, b], cnt, perm_p(a, b, rw)))
    step = 1 / (len(PERM) + 1)
    check = [abs(perm_p(a, b, rw, loop=True) - p) <= step + 1e-15 for a, b, rw, *_, p in cand[:150]]
    print(f"vectorised permutation p within one permutation step ({step:.2e}) of the original loop "
          f"on {len(check)} pairs: {all(check)}")
    cdf = pd.DataFrame(cand, columns=["taxon_a", "taxon_b", "r_within", "r_marginal", "r_between",
                                      "cooccurrence", "perm_p"])
    rej, q, _, _ = multipletests(cdf.perm_p.values, alpha=FDR_ALPHA, method="fdr_bh")
    cdf["q_BH"], cdf["fdr_pass"] = q, rej
    cdf["effect_pass"] = cdf.r_within.abs() >= EFFECT_THRESHOLD
    cdf["both_pass"] = cdf.effect_pass & cdf.fdr_pass
    cdf.to_csv(HERE / "all_candidate_pairs_r232.csv", index=False)
    cdf[cdf.both_pass].sort_values(["q_BH", "perm_p"]).to_csv(HERE / "direct_interactions_r232.csv", index=False)
    rng = np.random.default_rng(42)
    iu = np.triu_indices(V, 1)
    maxT = []
    for _ in range(500):
        Rs = R.copy()
        for col in range(V):
            for ix in idxs:
                Rs[ix, col] = Rs[rng.permutation(ix), col]
        maxT.append(np.abs(np.corrcoef(Rs.T)[iu]).max())
    maxT = np.array(maxT)
    gl = GraphicalLassoCV(alphas=10, max_iter=200).fit(R)
    P = gl.precision_
    edges = [(a, b, -P[i, j] / np.sqrt(P[i, i] * P[j, j])) for i, a in enumerate(taxa) for j, b in enumerate(taxa)
             if j > i and abs(P[i, j]) > 1e-8]
    gdf = pd.DataFrame(edges, columns=["taxon_a", "taxon_b", "partial_r"]).sort_values(
        "partial_r", key=lambda s: s.abs(), ascending=False)
    gdf.to_csv(HERE / "glasso_direct_edges_r232.csv", index=False)
    for name, m in (("C_marginal", Cm), ("C_within", Cw), ("C_between", Cb)):
        m.to_csv(HERE / f"{name}_r232.csv")
    json.dump({"condition": cond.tolist(), "samples": list(long.index), "taxa": list(taxa)},
              open(HERE / "within_condition_meta_r232.json", "w"))
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(C_marg[iu], C_within[iu], s=2, alpha=0.3, color="#2a78d6")
    ax.plot([-1, 1], [-1, 1], "k--", lw=1)
    for s in (1, -1):
        ax.axhline(s * EFFECT_THRESHOLD, color="#9a9a9a", ls=":", lw=0.8)
    ax.set_xlabel("marginal r (12 samples, Pearson on CLR)"); ax.set_ylabel("within-condition r (residuals)")
    ax.set_title("Coupling vs. condition response (GTDB r232 taxa)"); ax.set_aspect("equal")
    fig.tight_layout(); fig.savefig(HERE / "marginal_vs_within_r232.png", dpi=200); plt.close(fig)
    arch = sum(taxon_phylum_kingdom()[t][1] == "Archaea" for t in taxa)
    print(f"within-condition: params EFFECT={EFFECT_THRESHOLD} FDR={FDR_ALPHA} MIN_COOCC={MIN_COOCC} ZERO={ZERO} "
          f"pseudocount={pseudo:.3e} GraphicalLassoCV(alphas=10,max_iter=200) alpha_={gl.alpha_:.4g}")
    print(f"  taxa {V} ({arch} archaeal), candidates {len(cdf)}, effect-pass {int(cdf.effect_pass.sum())}, "
          f"FDR-pass {int(cdf.fdr_pass.sum())}, both {int(cdf.both_pass.sum())}; maxT 95th pct |r| = "
          f"{np.percentile(maxT, 95):.3f}; Glasso edges {len(gdf)}")


def network():
    cmd = [sys.executable, "-I", str(COR / "louvain_cooccurrence_network.py"),
           "--edges", str(HERE / "glasso_direct_edges_r232.csv"), "--meta", str(HERE / "within_condition_meta_r232.json"),
           "--ids", str(COR / "iterativeIDs_r232.json"), "--phylogeny", str(S / "data/Saltern_phylogeny_r232.json"),
           "--abundances", str(COR / "averaged_normalized_MAG_abundances.csv"),
           "--out-png", str(HERE / "within_condition_louvain_network_r232.png"),
           "--out-modules", str(HERE / "within_condition_louvain_modules_r232.json"),
           "--subtitle", "GTDB r232 labels, archaea included; glasso_direct_edges_r232.csv"]
    print(subprocess.run(cmd, capture_output=True, text=True, check=True).stdout)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parts", default="guilds,corr,heatmap,network")
    parts = set(ap.parse_args().parts.split(","))
    g = build_guilds() if parts & {"guilds", "corr", "heatmap"} else None
    if "corr" in parts:
        corr = guild_correlations(g)
        per_mag_correlations(g)
        for basis in AB:
            plot_guilds(g, corr, basis)
        show = corr[corr.basis == "raw"][["label", "n_MAGs", "cores", "spearman_rho", "p", "null_frac_ge_rho",
                                           "mean_rel_abund", "cores_detected"]]
        print(show.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    if "heatmap" in parts:
        heatmap(g)
    if "network" in parts:
        within_condition()
        network()
