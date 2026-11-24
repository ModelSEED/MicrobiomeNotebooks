#!/usr/bin/env python
"""Sensitivity checks of the within-site co-occurrence network for the final MAG set (585 species), the checks the
manuscript cites for the R232 network (../r232/prefilter_sensitivity_within.py --null 20 / --variance), recomputed here:
  null      GraphicalLassoCV on N within-site shuffles of the baseline residuals (the first N draws of the maxT loop,
            seed 42): edge counts, |partial r| > 0.05 edges, same-phylum share, archaeal-archaeal edges, Louvain Q
  variance  Spearman(residual variance, Glasso degree), and a refit on standardised residuals with the degrees of the
            guild taxa and named genera
  prefilter the MAG-level prevalence filter relaxed (>= 1 core, >= 0 cores): edges of guild taxa vs the baseline
Statistics as within_condition() in rerun_correlational_final585.py (same residuals; asserted against its outputs).
Inputs (this folder): within_condition_meta_final585.json, glasso_direct_edges_final585.csv, iterativeIDs_final585.json,
phylogeny_final585.json, guilds_final585.csv, averaged_normalized_MAG_abundances_final585.csv (the DNA-merged 12-core
table written by rerun_correlational_final585.py --parts network).
Writes network_sensitivity_final585/{null,variance,prefilter}_*.csv and summary.json.
Run: python network_sensitivity_final585.py [--null 20] [--workers 4] [--parts null,variance,prefilter]
"""
import argparse
import collections
import json
import os
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats
from skbio.stats.composition import clr
from sklearn.covariance import GraphicalLassoCV

HERE = Path(__file__).resolve().parent
TAG = "final585"
OUT = HERE / f"network_sensitivity_{TAG}"
OUT.mkdir(exist_ok=True)
ZERO, PARTIAL_THRESHOLD, MAXT_SEED = 1e-5, 0.05, 42

ids = json.load(open(HERE / f"iterativeIDs_{TAG}.json"))
phylo = json.load(open(HERE / f"phylogeny_{TAG}.json"))
taxon_of = {m: i.rsplit(".", 1)[0] for m, i in ids.items()}
PHYLUM, KINGDOM = {}, {}
for m, t in phylo.items():
    PHYLUM.setdefault(taxon_of[m], t["Phylum"]); KINGDOM.setdefault(taxon_of[m], t["Kingdom"])
g = pd.read_csv(HERE / f"guilds_{TAG}.csv")
g["taxon"] = g[f"iterativeID_{TAG}"].str.rsplit(".", n=1).str[0]
GUILD_TAXA = {c: set(g.loc[g[c], "taxon"]) for c in ("mepn_consumer", "betaine_reducer", "methanogen")}
NAMED = ["Roseovarius", "Albimonas", "Marivita", "Halomonas", "Vreelandella", "Halanaerobium", "Methanolobus"]
KEY_TAXA = sorted(set().union(*GUILD_TAXA.values()) | set(NAMED))
ab12 = pd.read_csv(HERE / f"averaged_normalized_MAG_abundances_{TAG}.csv", index_col=0)   # MAG x 12 cores
sites = sorted({c.split("_")[0] for c in ab12.columns})
ORDER = [f"{s}_{r}" for s in sites for r in ("A", "B", "C")]
LONG = ab12[ORDER].T
COND = np.array([c.split("_")[0] for c in ORDER])
IDXS = [np.where(COND == c)[0] for c in np.unique(COND)]


def residuals(min_prev):
    df = LONG.loc[:, (LONG.fillna(0) > ZERO).sum() >= min_prev]
    df = df.copy(); df.columns = [taxon_of[c] for c in df.columns]
    long = df.T.groupby(df.columns).sum().T
    vals = long.values
    X = clr(vals + vals[vals > 0].min() / 2)
    R = X.copy()
    for c in np.unique(COND):
        R[COND == c] -= R[COND == c].mean(axis=0, keepdims=True)
    return R, long.columns.to_numpy(), long


def glasso_edges(R, taxa):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gl = GraphicalLassoCV(alphas=10, max_iter=200).fit(R)
    P = gl.precision_
    e = [(a, b, -P[i, j] / np.sqrt(P[i, i] * P[j, j])) for i, a in enumerate(taxa)
         for j, b in enumerate(taxa) if j > i and abs(P[i, j]) > 1e-8]
    return pd.DataFrame(e, columns=["taxon_a", "taxon_b", "partial_r"]), float(gl.alpha_)


def louvain(gdf):
    """Same graph and Louvain call as louvain_cooccurrence_network.py (no drawing)."""
    f = gdf[gdf.partial_r.abs() > PARTIAL_THRESHOLD]
    G = nx.Graph()
    for r in f.itertuples():
        G.add_edge(r.taxon_a, r.taxon_b, partial_r=r.partial_r, weight=abs(r.partial_r))
    if not G.number_of_edges():
        return {"nodes": 0, "modules": 0, "Q": np.nan}
    pos = G.edge_subgraph([(u, v) for u, v, d in G.edges(data=True) if d["partial_r"] > 0]).copy()
    comms = list(nx.community.louvain_communities(pos, weight="weight", resolution=1.0, seed=42)) if pos.number_of_edges() else []
    multi = [c for c in comms if len(c) >= 2]
    return {"nodes": G.number_of_nodes(), "modules": len(multi),
            "Q": nx.community.modularity(pos, multi, weight="weight") if multi else np.nan}


def coherence(gdf):
    e = gdf[gdf.partial_r.abs() > PARTIAL_THRESHOLD]
    if not len(e):
        return {"edges_gt0.05": 0, "same_phylum_edges": 0, "same_phylum_frac": np.nan, "arch_arch_edges": 0}
    same = [PHYLUM[a] == PHYLUM[b] for a, b in zip(e.taxon_a, e.taxon_b)]
    aa = [KINGDOM[a] == KINGDOM[b] == "Archaea" for a, b in zip(e.taxon_a, e.taxon_b)]
    return {"edges_gt0.05": len(e), "same_phylum_edges": int(np.sum(same)), "same_phylum_frac": float(np.mean(same)),
            "arch_arch_edges": int(np.sum(aa))}


def shuffled(R, n):
    """The first n within-site shuffles of the maxT loop in within_condition() (seed 42, same draw order)."""
    rng = np.random.default_rng(MAXT_SEED)
    out = []
    for _ in range(n):
        Rs = R.copy()
        for col in range(R.shape[1]):
            for ix in IDXS:
                Rs[ix, col] = Rs[rng.permutation(ix), col]
        out.append(Rs)
    return out


def _null_fit(args):
    it, Rs, taxa, standardise = args
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    t0 = time.time()
    ng, alpha = glasso_edges(Rs / Rs.std(0) if standardise else Rs, taxa)
    return {"shuffle": it, "standardised": standardise, "glasso_alpha": alpha, "glasso_edges": len(ng), **coherence(ng),
            **{f"louvain_{k}": v for k, v in louvain(ng).items()}, "seconds": round(time.time() - t0, 1)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--null", type=int, default=20)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--parts", default="null,variance,prefilter")
    a = ap.parse_args()
    parts = set(a.parts.split(","))
    R, taxa, long = residuals(6)
    meta = json.load(open(HERE / f"within_condition_meta_{TAG}.json"))
    assert list(taxa) == meta["taxa"], "baseline taxa differ from within_condition()"
    obs = pd.read_csv(HERE / f"glasso_direct_edges_{TAG}.csv")
    summ_f = OUT / "summary.json"
    summ = json.load(open(summ_f)) if summ_f.exists() else {}
    o = {"glasso_edges": len(obs), **coherence(obs), **{f"louvain_{k}": v for k, v in louvain(obs).items()}}
    summ["observed"] = o
    print("observed:", o, flush=True)
    if "null" in parts:
        Rs = shuffled(R, a.null)
        jobs = [(i, x, taxa, False) for i, x in enumerate(Rs)] + [(i, x, taxa, True) for i, x in enumerate(Rs)]
        with ProcessPoolExecutor(a.workers) as ex:
            rows = list(ex.map(_null_fit, jobs))
        d = pd.DataFrame(rows)
        d.to_csv(OUT / "null_glasso_shuffles.csv", index=False)
        for std, sub in d.groupby("standardised"):
            k = "null_standardised" if std else "null"
            summ[k] = {"n": len(sub), "edges_min": int(sub.glasso_edges.min()), "edges_max": int(sub.glasso_edges.max()),
                       "edges_median": float(sub.glasso_edges.median()),
                       "edges_gt005_min": int(sub["edges_gt0.05"].min()), "edges_gt005_max": int(sub["edges_gt0.05"].max()),
                       "same_phylum_frac_median": float(sub.same_phylum_frac.median()),
                       "same_phylum_frac_max": float(sub.same_phylum_frac.max()),
                       "arch_arch_max": int(sub.arch_arch_edges.max()),
                       "Q_min": float(sub.louvain_Q.min()), "Q_max": float(sub.louvain_Q.max()),
                       "Q_median": float(sub.louvain_Q.median()),
                       "frac_edges_ge_obs": float((sub.glasso_edges >= o["glasso_edges"]).mean()),
                       "frac_Q_ge_obs": float((sub.louvain_Q >= o["louvain_Q"]).mean())}
            print(k, summ[k], flush=True)
    if "variance" in parts:
        deg = pd.concat([obs.taxon_a, obs.taxon_b]).value_counts().reindex(taxa, fill_value=0)
        var = pd.Series(R.var(0), index=taxa)
        rho, p = stats.spearmanr(var, deg)
        gz, alpha = glasso_edges(R / R.std(0), taxa)
        gz.to_csv(OUT / "variance_standardised_edges.csv", index=False)
        gzf = gz[gz.partial_r.abs() > PARTIAL_THRESHOLD]
        dz = collections.Counter(list(gzf.taxon_a) + list(gzf.taxon_b))
        obsf = obs[obs.partial_r.abs() > PARTIAL_THRESHOLD]
        do = collections.Counter(list(obsf.taxon_a) + list(obsf.taxon_b))
        rows = [{"taxon": t, "in_panel": t in set(taxa), "resid_var": float(var.get(t, np.nan)),
                 "degree_gt005": do.get(t, 0), "degree_gt005_standardised": dz.get(t, 0),
                 "guilds": ",".join(c for c, s in GUILD_TAXA.items() if t in s)} for t in KEY_TAXA]
        pd.DataFrame(rows).to_csv(OUT / "variance_key_taxa.csv", index=False)
        pd.DataFrame({"taxon": taxa, "resid_var": var.values, "glasso_degree": deg.values}).to_csv(
            OUT / "variance_degree.csv", index=False)
        summ["variance"] = {"spearman_resid_var_vs_degree": float(rho), "p": float(p),
                            "standardised_alpha": alpha, "standardised_edges": len(gz),
                            "standardised_edges_gt005": len(gzf),
                            "standardised_louvain": louvain(gz), "standardised_coherence": coherence(gz),
                            "key_taxa_in_panel_without_edges": [r["taxon"] for r in rows if r["in_panel"] and not r["degree_gt005"]],
                            "their_standardised_degree": {r["taxon"]: r["degree_gt005_standardised"] for r in rows
                                                          if r["in_panel"] and not r["degree_gt005"]}}
        print("variance", summ["variance"], flush=True)
    if "prefilter" in parts:
        base = obs[obs.partial_r.abs() > PARTIAL_THRESHOLD]
        bdeg = collections.Counter(list(base.taxon_a) + list(base.taxon_b))
        out = {}
        for mp in (1, 0):
            Rm, tm, _ = residuals(mp)
            gm, alpha = glasso_edges(Rm, tm)
            gm.to_csv(OUT / f"prefilter_min{mp}_edges.csv", index=False)
            gf = gm[gm.partial_r.abs() > PARTIAL_THRESHOLD]
            dm = collections.Counter(list(gf.taxon_a) + list(gf.taxon_b))
            gained = {t: dm[t] for t in KEY_TAXA if dm.get(t, 0) and not bdeg.get(t, 0)}
            new_edges = gf[~gf.apply(lambda r: frozenset((r.taxon_a, r.taxon_b)) in
                                     {frozenset(x) for x in zip(base.taxon_a, base.taxon_b)}, axis=1)]
            out[f"min{mp}"] = {"taxa": len(tm), "alpha": alpha, "edges": len(gm), "edges_gt005": len(gf),
                               "key_taxa_gaining_edges": gained, "new_edges_gt005": len(new_edges),
                               **{f"louvain_{k}": v for k, v in louvain(gm).items()}}
            print(f"prefilter min{mp}", out[f"min{mp}"], flush=True)
        summ["prefilter"] = out
    json.dump(summ, open(summ_f, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
