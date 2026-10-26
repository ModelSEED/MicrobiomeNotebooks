#!/usr/bin/env python
"""Prefilter sensitivity of the within-condition co-occurrence pipeline (GTDB r232 labels).

What is varied (everything else is the 2020 / r232 pipeline, copied from within_condition() in
rerun_correlational_r232.py, which is NOT imported because it loads its inputs at import time):
  * MIN_COOCC in {6 (baseline), 4, 3, 2, 1, 0}.  The one parameter acts twice in the original code:
      (a) MAG-level prevalence filter BEFORE the taxon collapse  (rerun_correlational_r232.py:386,
          analysis_2_marimo.py:1996): keep MAGs with rel. abundance > ZERO in >= MIN_COOCC of 12 cores;
      (b) pair-level co-occurrence filter on the collapsed taxa for the permutation/BH test
          (rerun_correlational_r232.py:435, analysis_2_marimo.py:2058).
    Here both move together (as in the code); MIN_COOCC = 1 is the task's "no prevalence filter", and
    0 additionally keeps the 25 MAGs never above 1e-5 in any core (and tests pairs that never co-occur).
  * the 2020 "defined IDs only" restriction (analysis_2_marimo.py:1997-2001): ON = drop the 30 MAGs
    that had no 2020 iterativeID (all archaea) before the collapse; labels stay r232 so the settings are
    comparable. This is not the 140-taxon 2020 reproduction (that used the 2020 labels).
  * EFFECT gate (|r_within| >= 0.7) on/off: it is a post-hoc label in this pipeline
    (rerun_correlational_r232.py:449), so "off" only changes both_pass -> fdr_pass; reported as columns.
  * a CLR decomposition run: MIN_COOCC = 1 and 0 with the pseudocount pinned to the baseline value, to
    separate "more taxa" from "smaller pseudocount" (the pseudocount = half the smallest non-zero value
    of the FILTERED matrix, rerun_correlational_r232.py:389).
Extra diagnostics: null counts of |r| >= 0.7 from the maxT shuffles, the BH resolution arithmetic for
the 1,296-permutation p-value floor, the share of residual variance carried by sub-detection entries,
Graphical Lasso refits on within-condition-shuffled residuals (the maxT null), edge Jaccard and
Louvain adjusted Rand index against the baseline.

Outputs: prefilter_sensitivity/within_*  (see within_REPORT.md).  No existing file is modified.
Run (from this folder), in this order:
  ~/Documents/py_venv/bin/python prefilter_sensitivity_within.py            # 14 settings (~15 min)
  ~/Documents/py_venv/bin/python prefilter_sensitivity_within.py --null 20  # 20-shuffle Glasso null (~20 min)
  ~/Documents/py_venv/bin/python prefilter_sensitivity_within.py --variance # variance / standardised Glasso
  (--quick: baseline + verification against the committed r232 files only)
"""
import argparse
import json
import math
import re
import sys
import time
import warnings
from itertools import permutations, product
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats
from skbio.stats.composition import clr
from sklearn.covariance import GraphicalLassoCV
from sklearn.metrics import adjusted_rand_score
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent          # correlational/r232
COR = HERE.parent
S = COR.parent
OUT = HERE / "prefilter_sensitivity"
OUT.mkdir(exist_ok=True)
PAT = re.compile(r"Salt_Pond_MetaG_?(.+?)_MG_DASTool_bins[._](concoct_out|metabat|maxbin)\.(\d+)")
ZERO, EFFECT, FDR, PARTIAL_THRESHOLD = 1e-5, 0.7, 0.05, 0.05
N_MAXT, MAXT_SEED, N_GLASSO_NULL = 500, 42, 3


class Tee:
    def __init__(self, path):
        self.f, self.o = open(path, "w"), sys.stdout

    def write(self, s):
        self.f.write(s); self.o.write(s); self.f.flush()

    def flush(self):
        self.f.flush(); self.o.flush()


# ---- inputs (as rerun_correlational_r232.py:73-83) ---------------------------------------------------
ids = json.load(open(COR / "iterativeIDs_r232.json"))
ids_2020 = json.load(open(COR / "iterativeIDs.json"))
rel = pd.read_csv(S / "data/Cliff_310MAG_relabund.tsv", sep="\t", index_col=0)
dna = pd.read_csv(S / "data/SaltPondsDNA.csv").set_index("Sample")["DNA conc (µg/g soil)"]
guilds = pd.read_csv(HERE / "guilds_r232.csv")
PRE = sorted({"_".join(s.split("_")[:-1]) for s in dna.index})
taxon_of = {m: i.rsplit(".", 1)[0] for m, i in ids.items()}


def merge_depths(df):
    """DNA-weighted D1/D2 merge (rerun_correlational_r232.py:86-89)."""
    return pd.DataFrame({p: (df[f"{p}_D1"] * dna[f"{p}_D1"] + df[f"{p}_D2"] * dna[f"{p}_D2"])
                         / (dna[f"{p}_D1"] + dna[f"{p}_D2"]) for p in PRE})


def long_form():
    """12 cores x MAGs (+ unbinned_reads), rows ordered R1, R2, R2A, SF2 x A, B, C (:378-385)."""
    ab = merge_depths(rel)
    sites = sorted({p.split("_")[0] for p in PRE})
    rows, cond = {}, []
    for site in sites:
        for rep in ("A", "B", "C"):
            rows[f"{site}_{rep}"] = ab[f"{site}_{rep}"]; cond.append(site)
    return pd.DataFrame(rows).T, np.array(cond)


LONG, COND = long_form()
MAG_PREV = (LONG.fillna(0) > ZERO).sum()          # cores (of 12) with rel. abundance > ZERO


def filter_collapse(min_prev, defined_2020):
    """(a) MAG prevalence filter, optional 2020-ID restriction, then sum MAGs per r232 taxon (:386-387)."""
    df = LONG.loc[:, (LONG.fillna(0) > ZERO).sum() >= min_prev]
    mags = [c for c in df.columns if c in ids and (not defined_2020 or c in ids_2020)]
    df = df[mags]
    df.columns = [taxon_of[c] for c in mags]
    return df.T.groupby(df.columns).sum().T, mags


# ---- guild roles per taxon (from guilds_r232.csv) ------------------------------------------------------
ROLE_COLS = [("mepn_consumer", "MePn(phnJ)"), ("mepn_2020", "MePn2020"), ("betaine_reducer", "betaine(grdH)"),
             ("betaine_2020", "betaine2020"), ("tma_producer", "TMA-cand"), ("methanogen", "methanogen(mcrA)")]
guilds["taxon"] = guilds.iterativeID_r232.str.rsplit(".", n=1).str[0]
KINGDOM = guilds.groupby("taxon").domain.first().to_dict()
ROLES = {}
for t, sub in guilds.groupby("taxon"):
    r = [f"{lab}:{','.join(sub.loc[sub[c], 'iterativeID_r232'])}" for c, lab in ROLE_COLS if sub[c].any()]
    if KINGDOM[t] == "Archaea":
        r.append("archaea")
    ROLES[t] = "; ".join(r)
NAMED = ["Roseovarius", "Albimonas", "Marivita", "Halomonas", "Vreelandella", "Halanaerobium", "Methanolobus"]
KEY_TAXA = sorted({t for t, r in ROLES.items() if r} | set(NAMED))


# ---- statistics (copied from within_condition(), rerun_correlational_r232.py:388-466) ----------------
def perm_matrix(cond):
    idxs = [np.where(cond == c)[0] for c in np.unique(cond)]
    PERM = []
    for combo in product(*[list(permutations(ix)) for ix in idxs]):
        pi = np.arange(len(cond))
        for orig, new in zip(idxs, combo):
            pi[orig] = list(new)
        PERM.append(pi)
    return np.array(PERM), idxs


PERM, IDXS = perm_matrix(COND)                    # 1296 x 12, includes the identity


def perm_pvalues(R, absCw):
    """Exact within-condition permutation p for every pair, same statistic and tie rule as perm_p()
    (:415-431): column j (taxon_b) permuted, column i (taxon_a) fixed, null + 1e-12 >= |r_within|."""
    V = R.shape[1]
    norms = np.sqrt((R ** 2).sum(0))
    cnt = np.zeros((V, V), dtype=np.int64)
    with np.errstate(divide="ignore", invalid="ignore"):
        for j in range(V):
            num = np.abs(R[PERM, j] @ R)                       # 1296 x V
            d = norms[j] * norms
            null = np.where(d > 0, num / d, 0.0) + 1e-12
            cnt[:, j] = (null >= absCw[:, j]).sum(0)
    return (cnt + 1) / (len(PERM) + 1)


def louvain(gdf):
    """Same graph and Louvain call as louvain_cooccurrence_network.py:95-125 (no drawing)."""
    filt = gdf[gdf.partial_r.abs() > PARTIAL_THRESHOLD]
    G = nx.Graph()
    for _, row in filt.iterrows():
        G.add_edge(row["taxon_a"], row["taxon_b"], partial_r=row["partial_r"], weight=abs(row["partial_r"]))
    G.remove_nodes_from(list(nx.isolates(G)))
    out = {"nodes": G.number_of_nodes(), "edges": G.number_of_edges(), "modules": 0, "singletons": 0,
           "Q": np.nan, "sizes": [], "export": {}, "node_module": {}}
    if G.number_of_edges() == 0:
        return out
    pos_sub = G.edge_subgraph([(u, v) for u, v, d in G.edges(data=True) if d["partial_r"] > 0]).copy()
    comms = (list(nx.community.louvain_communities(pos_sub, weight="weight", resolution=1.0, seed=42))
             if pos_sub.number_of_edges() else [])
    in_mod = {n for c in comms for n in c}
    comms += [{n} for n in G.nodes() if n not in in_mod]
    multi = [c for c in comms if len(c) >= 2]
    out["Q"] = nx.community.modularity(pos_sub, multi, weight="weight") if multi else np.nan
    out["modules"], out["singletons"] = len(multi), len(comms) - len(multi)
    out["sizes"] = sorted([len(c) for c in multi], reverse=True)
    ln = sn = 0
    for comm in sorted(comms, key=lambda c: -len(c)):
        if len(comm) >= 2:
            ln += 1; name = f"module_{ln}"
            out["export"][name] = {"size": len(comm), "type": "louvain", "members": sorted(comm)}
        else:
            sn += 1; name = f"singleton_{sn}"
            out["export"][name] = {"size": 1, "type": "isolate_singleton", "members": sorted(comm)}
        for n in comm:
            out["node_module"][n] = name
    return out


def glasso_edges(R, taxa):
    """GraphicalLassoCV(alphas=10, max_iter=200) on the residuals; edges |P_ij| > 1e-8 (:463-468)."""
    info = {"alpha": np.nan, "n_iter": np.nan, "precision_finite": None, "runtime_warnings": 0, "error": ""}
    t0 = time.time()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            gl = GraphicalLassoCV(alphas=10, max_iter=200).fit(R)
        except Exception as exc:                                  # a failure is a result
            info["error"] = f"{type(exc).__name__}: {exc}"[:300]
            gl = None
    info["runtime_warnings"] = sum(issubclass(x.category, RuntimeWarning) for x in w)
    info["convergence_warnings"] = sum("converge" in str(x.message).lower() for x in w)
    info["glasso_seconds"] = round(time.time() - t0, 1)
    if gl is None:
        return pd.DataFrame(columns=["taxon_a", "taxon_b", "partial_r"]), info
    P = gl.precision_
    info.update(alpha=gl.alpha_, n_iter=gl.n_iter_, precision_finite=bool(np.isfinite(P).all()),
                alpha_grid_min=float(np.min(gl.cv_results_["alphas"])),
                alpha_grid_max=float(np.max(gl.cv_results_["alphas"])))
    edges = [(a, b, -P[i, j] / np.sqrt(P[i, i] * P[j, j])) for i, a in enumerate(taxa)
             for j, b in enumerate(taxa) if j > i and abs(P[i, j]) > 1e-8]
    gdf = pd.DataFrame(edges, columns=["taxon_a", "taxon_b", "partial_r"]).sort_values(
        "partial_r", key=lambda s: s.abs(), ascending=False)
    return gdf, info


def run(label, min_prev, defined_2020=False, pseudo_fixed=None, glasso_null=0):
    t0 = time.time()
    long, mags = filter_collapse(min_prev, defined_2020)
    taxa = long.columns.to_numpy()
    assert list(taxa) == sorted(taxa)
    vals = long.values
    pseudo = pseudo_fixed if pseudo_fixed is not None else vals[vals > 0].min() / 2
    X = clr(vals + pseudo)
    V = X.shape[1]
    C_marg = np.corrcoef(X.T)
    R = X.copy()
    for c in np.unique(COND):
        R[COND == c] -= R[COND == c].mean(axis=0, keepdims=True)
    Cw = np.corrcoef(R.T)
    Cb, _ = stats.spearmanr(np.vstack([X[COND == c].mean(0) for c in np.unique(COND)]))
    pres = (long > ZERO).astype(int).values
    cooc = pres.T @ pres                                         # cores where both taxa > ZERO
    iu = np.triu_indices(V, 1)
    P = perm_pvalues(R, np.abs(Cw))
    allp = pd.DataFrame({"taxon_a": taxa[iu[0]], "taxon_b": taxa[iu[1]], "r_within": Cw[iu],
                         "r_marginal": C_marg[iu], "r_between": Cb[iu], "cooccurrence": cooc[iu],
                         "perm_p": P[iu]})
    allp = allp[np.isfinite(allp.r_within)]
    cdf = allp[allp.cooccurrence >= min_prev].copy()            # (b) pair co-occurrence filter
    rej, q, _, _ = multipletests(cdf.perm_p.values, alpha=FDR, method="fdr_bh")
    cdf["q_BH"], cdf["fdr_pass"] = q, rej
    cdf["effect_pass"] = cdf.r_within.abs() >= EFFECT
    cdf["both_pass"] = cdf.effect_pass & cdf.fdr_pass
    cand6 = allp[allp.cooccurrence >= 6]                         # decoupled: pair filter kept at 6
    rej6 = multipletests(cand6.perm_p.values, alpha=FDR, method="fdr_bh")[0] if len(cand6) else []
    # maxT loop verbatim (:453-462), plus null counts of |r| >= EFFECT; first shuffles kept for the Glasso null
    rng = np.random.default_rng(MAXT_SEED)
    tested = np.zeros((V, V), bool)
    ia = {t: k for k, t in enumerate(taxa)}
    tested[cdf.taxon_a.map(ia).values, cdf.taxon_b.map(ia).values] = True
    tested = tested[iu]
    maxT, null_eff_all, null_eff_tested, shuffled = [], [], [], []
    for it in range(N_MAXT):
        Rs = R.copy()
        for col in range(V):
            for ix in IDXS:
                Rs[ix, col] = Rs[rng.permutation(ix), col]
        cu = np.abs(np.corrcoef(Rs.T)[iu])
        maxT.append(cu.max()); null_eff_all.append((cu >= EFFECT).sum()); null_eff_tested.append((cu[tested] >= EFFECT).sum())
        if it < glasso_null:
            shuffled.append(Rs)
    maxT = np.array(maxT)
    gdf, ginfo = glasso_edges(R, taxa)
    lv = louvain(gdf)
    # p-value resolution vs BH
    m = len(cdf)
    pmin_possible = 2 / (len(PERM) + 1)                          # identity permutation always counted
    steps = np.arange(2, 7) / (len(PERM) + 1)
    obs_steps = [int(np.isclose(cdf.perm_p, s).sum()) for s in steps]
    subdet = vals <= ZERO
    k = taxa.tolist()
    prev_taxon = pd.Series(pres.sum(0), index=taxa)
    big = gdf[gdf.partial_r.abs() > PARTIAL_THRESHOLD]
    low_prev = set(prev_taxon[prev_taxon < 6].index)
    summ = {
        "setting": label, "MIN_COOCC": min_prev, "defined_2020_ids_only": defined_2020,
        "pseudocount": pseudo, "pseudocount_pinned": pseudo_fixed is not None,
        "MAGs_kept": len(mags), "archaeal_MAGs_kept": int(sum(KINGDOM[taxon_of[x]] == "Archaea" for x in mags)),
        "taxa": V, "archaeal_taxa": int(sum(KINGDOM[t] == "Archaea" for t in taxa)),
        "taxa_prev_lt6": len(low_prev),
        "frac_entries_exact_zero": float((vals == 0).mean()), "frac_entries_le_ZERO": float(subdet.mean()),
        "share_residual_SS_from_le_ZERO_entries": float((R[subdet] ** 2).sum() / (R ** 2).sum()),
        "pairs_total": len(iu[0]), "pairs_tested": m,
        "pairs_cooc0": int((allp.cooccurrence == 0).sum()),
        "pairs_cooc0_absr_ge_0.7": int(((allp.cooccurrence == 0) & (allp.r_within.abs() >= EFFECT)).sum()),
        "min_perm_p": float(cdf.perm_p.min()), "perm_p_floor_possible": pmin_possible,
        "n_at_floor": obs_steps[0], "null_expected_at_floor": m / len(PERM),
        "n_at_floor_needed_for_any_BH": math.ceil(m * pmin_possible / FDR),
        "obs_counts_first5_p_steps": " ".join(map(str, obs_steps)),
        "BH_pass": int(cdf.fdr_pass.sum()), "min_q_BH": float(cdf.q_BH.min()),
        "BH_pass_pair_filter_fixed6": int(np.sum(rej6)),
        "effect_pass_obs": int(cdf.effect_pass.sum()),
        "effect_pass_null_mean": float(np.mean(null_eff_tested)),
        "effect_pass_null_95pct": float(np.percentile(null_eff_tested, 95)),
        "both_pass_gate_on": int(cdf.both_pass.sum()), "pass_gate_off_FDR_only": int(cdf.fdr_pass.sum()),
        "maxT_95": float(np.percentile(maxT, 95)), "maxT_99": float(np.percentile(maxT, 99)),
        "glasso_alpha": ginfo["alpha"], "glasso_alpha_grid": f"{ginfo.get('alpha_grid_min', np.nan):.3g}-"
                                                            f"{ginfo.get('alpha_grid_max', np.nan):.3g}",
        "glasso_n_iter": ginfo["n_iter"], "glasso_precision_finite": ginfo["precision_finite"],
        "glasso_runtime_warnings": ginfo["runtime_warnings"], "glasso_convergence_warnings": ginfo["convergence_warnings"],
        "glasso_error": ginfo["error"], "glasso_edges": len(gdf), "glasso_edges_gt0.05": len(big),
        "glasso_edges_gt0.05_touching_prev_lt6": int(sum((a in low_prev) or (b in low_prev)
                                                         for a, b in zip(big.taxon_a, big.taxon_b))),
        "louvain_nodes": lv["nodes"], "louvain_edges": lv["edges"], "louvain_modules": lv["modules"],
        "louvain_singletons": lv["singletons"], "louvain_Q": lv["Q"],
        "louvain_sizes": " ".join(map(str, lv["sizes"])),
    }
    # Glasso on within-condition-shuffled residuals (the maxT null: coupling destroyed, design kept)
    nulls = []
    for s_i, Rs in enumerate(shuffled):
        ng, ni = glasso_edges(Rs, taxa)
        nl = louvain(ng)
        nulls.append({"setting": label, "shuffle": s_i, "taxa": V, "glasso_alpha": ni["alpha"],
                      "glasso_edges": len(ng), "glasso_edges_gt0.05": int((ng.partial_r.abs() > PARTIAL_THRESHOLD).sum()),
                      "max_abs_partial_r": float(ng.partial_r.abs().max()) if len(ng) else np.nan,
                      "louvain_nodes": nl["nodes"], "louvain_modules": nl["modules"], "louvain_Q": nl["Q"],
                      "glasso_error": ni["error"]})
    # key taxa
    key_rows = []
    for t in KEY_TAXA:
        r = {"setting": label, "taxon": t, "kingdom": KINGDOM[t], "roles": ROLES[t],
             "MAGs_in_taxon": " ".join(sorted(ids[x] for x in ids if taxon_of[x] == t)),
             "MAG_prevalence": " ".join(f"{ids[x]}:{MAG_PREV[x]}" for x in sorted(
                 (x for x in ids if taxon_of[x] == t), key=lambda x: ids[x])),
             "in_matrix": t in k}
        if t in k:
            r["MAGs_included"] = " ".join(sorted(ids[x] for x in mags if taxon_of[x] == t))
            r["taxon_prevalence"] = int(prev_taxon[t])
            sub = cdf[(cdf.taxon_a == t) | (cdf.taxon_b == t)]
            r["pairs_tested"] = len(sub)
            r["min_perm_p"] = float(sub.perm_p.min()) if len(sub) else np.nan
            r["min_q_BH"] = float(sub.q_BH.min()) if len(sub) else np.nan
            r["effect_pass_pairs"] = int(sub.effect_pass.sum())
            e = gdf[(gdf.taxon_a == t) | (gdf.taxon_b == t)]
            fmt = lambda s: " ".join(f"{(x.taxon_b if x.taxon_a == t else x.taxon_a)}({x.partial_r:+.3f})"
                                     for x in s.itertuples())
            r["glasso_degree"] = len(e); r["glasso_partners"] = fmt(e)
            e5 = e[e.partial_r.abs() > PARTIAL_THRESHOLD]
            r["degree_gt0.05"] = len(e5)
            r["louvain_module"] = lv["node_module"].get(t, "")
            mod = lv["export"].get(r["louvain_module"], {}).get("members", [])
            r["module_size"] = len(mod)
            r["module_members"] = " ".join(x for x in mod if x != t)
        key_rows.append(r)
    taxon_rows = pd.DataFrame({"setting": label, "taxon": taxa, "taxon_prevalence": prev_taxon.values,
                               "glasso_degree": [int(((gdf.taxon_a == t) | (gdf.taxon_b == t)).sum()) for t in taxa],
                               "louvain_module": [lv["node_module"].get(t, "") for t in taxa]})
    summ["runtime_s"] = round(time.time() - t0, 1)
    print(f"[{label}] MAGs {len(mags)} taxa {V} ({summ['archaeal_taxa']} arch) pseudo {pseudo:.3e} "
          f"tested {m} BH {summ['BH_pass']} floor {obs_steps[0]}/{summ['n_at_floor_needed_for_any_BH']} needed "
          f"effect {summ['effect_pass_obs']} (null {summ['effect_pass_null_mean']:.0f}) maxT95 {summ['maxT_95']:.3f} "
          f"glasso alpha {ginfo['alpha']:.4g} edges {len(gdf)} (>0.05: {len(big)}) louvain {lv['nodes']}n/"
          f"{lv['edges']}e/{lv['modules']}m Q={lv['Q']:.3f}  [{summ['runtime_s']} s]", flush=True)
    return dict(summary=summ, cdf=cdf, gdf=gdf, lv=lv, key=key_rows, taxa_rows=taxon_rows, nulls=nulls)


PHYLUM = guilds.groupby("taxon").phylum.first().to_dict()


def coherence(gdf):
    """Share of |partial r| > 0.05 edges joining taxa of the same phylum / both archaeal."""
    e = gdf[gdf.partial_r.abs() > PARTIAL_THRESHOLD]
    if not len(e):
        return {"edges_gt0.05": 0, "same_phylum_frac": np.nan, "arch_arch_edges": 0, "pos_frac": np.nan}
    same = [PHYLUM[a] == PHYLUM[b] for a, b in zip(e.taxon_a, e.taxon_b)]
    aa = [KINGDOM[a] == KINGDOM[b] == "Archaea" for a, b in zip(e.taxon_a, e.taxon_b)]
    return {"edges_gt0.05": len(e), "same_phylum_frac": float(np.mean(same)), "arch_arch_edges": int(np.sum(aa)),
            "pos_frac": float((e.partial_r > 0).mean())}


def residuals(min_prev, defined_2020=False, pseudo_fixed=None):
    long, _ = filter_collapse(min_prev, defined_2020)
    vals = long.values
    X = clr(vals + (pseudo_fixed if pseudo_fixed is not None else vals[vals > 0].min() / 2))
    R = X.copy()
    for c in np.unique(COND):
        R[COND == c] -= R[COND == c].mean(axis=0, keepdims=True)
    return R, long.columns.to_numpy()


def glasso_null_extended(n):
    """n within-condition shuffles of the residuals (the first n draws of the maxT loop, seed 42) refitted with
    GraphicalLassoCV, for the baseline and the no-filter setting; observed coherence for every setting."""
    rows = []
    for label, mp in (("base_min6", 6), ("min0", 0)):
        R, taxa = residuals(mp)
        V = R.shape[1]
        rng = np.random.default_rng(MAXT_SEED)
        for it in range(n):
            Rs = R.copy()
            for col in range(V):
                for ix in IDXS:
                    Rs[ix, col] = Rs[rng.permutation(ix), col]
            ng, ni = glasso_edges(Rs, taxa)
            nl = louvain(ng)
            rows.append({"setting": label, "shuffle": it, "taxa": V, "glasso_alpha": ni["alpha"],
                         "glasso_edges": len(ng), **coherence(ng),
                         "max_abs_partial_r": float(ng.partial_r.abs().max()) if len(ng) else np.nan,
                         "louvain_nodes": nl["nodes"], "louvain_modules": nl["modules"], "louvain_Q": nl["Q"],
                         "glasso_error": ni["error"]})
            print(f"null {label} shuffle {it}: alpha {ni['alpha']:.3g} edges {len(ng)} "
                  f"(>0.05 {rows[-1]['edges_gt0.05']}, same-phylum {rows[-1]['same_phylum_frac']:.2f}, "
                  f"arch-arch {rows[-1]['arch_arch_edges']}) Q {nl['Q']:.3f}", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "within_glasso_null_extended.csv", index=False)
    obs = pd.read_csv(OUT / "within_glasso_edges_by_setting.csv")
    summ = pd.read_csv(OUT / "within_settings_summary.csv").set_index("setting")
    orows = [{"setting": s, "glasso_edges": len(g), **coherence(g), "louvain_Q": summ.loc[s, "louvain_Q"],
              "louvain_modules": summ.loc[s, "louvain_modules"], "glasso_alpha": summ.loc[s, "glasso_alpha"]}
             for s, g in obs.groupby("setting", sort=False)]
    pd.DataFrame(orows).to_csv(OUT / "within_glasso_observed_coherence.csv", index=False)
    nd = pd.DataFrame(rows)
    for label in ("base_min6", "min0"):
        o = next(r for r in orows if r["setting"] == label)
        d = nd[nd.setting == label]
        print(f"{label}: observed edges {o['glasso_edges']} (>0.05 {o['edges_gt0.05']}), same-phylum "
              f"{o['same_phylum_frac']:.3f}, arch-arch {o['arch_arch_edges']}, Q {o['louvain_Q']:.3f} | null (n={len(d)}) "
              f"edges {d.glasso_edges.min()}-{d.glasso_edges.max()} (median {d.glasso_edges.median():.0f}), "
              f">0.05 median {d['edges_gt0.05'].median():.0f}, same-phylum median {d.same_phylum_frac.median():.3f} "
              f"(max {d.same_phylum_frac.max():.3f}), arch-arch max {d.arch_arch_edges.max()}, Q median "
              f"{d.louvain_Q.median():.3f} (range {d.louvain_Q.min():.3f}-{d.louvain_Q.max():.3f}); "
              f"share of null fits with edges >= observed {(d.glasso_edges >= o['glasso_edges']).mean():.2f}, "
              f"Q >= observed {(d.louvain_Q >= o['louvain_Q']).mean():.2f}", flush=True)


def variance_check():
    """GraphicalLassoCV is fitted on the unstandardised residuals, so its L1 penalty acts on covariances:
    (1) per-taxon residual variance vs Glasso degree in every r232-label setting; (2) refit on standardised
    residuals (correlation scale) for base_min6 and min0, with the key-taxon degrees."""
    obs = pd.read_csv(OUT / "within_glasso_edges_by_setting.csv")
    summ = pd.read_csv(OUT / "within_settings_summary.csv")
    # pairs whose |r_within| exceeds the maxT 95th percentile (FWER 5 %) in each setting
    ex = []
    for s in summ.itertuples():
        R, taxa = residuals(s.MIN_COOCC, bool(s.defined_2020_ids_only),
                            s.pseudocount if s.pseudocount_pinned else None)
        Cw = np.corrcoef(R.T)
        iu = np.triu_indices(len(taxa), 1)
        hit = np.where(np.abs(Cw[iu]) >= s.maxT_95)[0]
        ex.append({"setting": s.setting, "max_abs_r_within": float(np.nanmax(np.abs(Cw[iu]))), "maxT_95": s.maxT_95,
                   "maxT_99": s.maxT_99, "pairs_ge_maxT95": len(hit),
                   "pairs": " ".join(f"{taxa[iu[0][h]]}~{taxa[iu[1][h]]}({Cw[iu][h]:+.4f})" for h in hit)})
        print(f"{s.setting}: max |r_within| {ex[-1]['max_abs_r_within']:.4f} vs maxT95 {s.maxT_95:.4f}: "
              f"{len(hit)} pair(s) {ex[-1]['pairs']}", flush=True)
    pd.DataFrame(ex).to_csv(OUT / "within_maxT_exceedances.csv", index=False)
    rows, stats_rows = [], []
    for label, mp in (("base_min6", 6), ("min4", 4), ("min3", 3), ("min2", 2), ("min1", 1), ("min0", 0)):
        R, taxa = residuals(mp)
        long, _ = filter_collapse(mp, False)
        g = obs[obs.setting == label]
        deg = pd.concat([g.taxon_a, g.taxon_b]).value_counts().reindex(taxa, fill_value=0)
        var = pd.Series(R.var(0), index=taxa)
        prev = (long > ZERO).sum()
        rho_v, p_v = stats.spearmanr(var, deg)
        rho_p, p_p = stats.spearmanr(prev, var)
        stats_rows.append({"setting": label, "taxa": len(taxa), "spearman_degree_vs_resid_var": rho_v, "p": p_v,
                           "spearman_prevalence_vs_resid_var": rho_p, "p_prev": p_p,
                           "median_var_degree0": float(var[deg == 0].median()),
                           "median_var_degree_ge1": float(var[deg > 0].median())})
        rows += [{"setting": label, "taxon": t, "resid_var": var[t], "taxon_prevalence": int(prev[t]),
                  "glasso_degree": int(deg[t]), "key": t in KEY_TAXA} for t in taxa]
        print(f"{label}: Spearman(resid var, Glasso degree) = {rho_v:+.2f} (p = {p_v:.1e}); "
              f"Spearman(prevalence, resid var) = {rho_p:+.2f}; median var degree 0 / >=1: "
              f"{stats_rows[-1]['median_var_degree0']:.3f} / {stats_rows[-1]['median_var_degree_ge1']:.3f}", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "within_glasso_degree_vs_variance.csv", index=False)
    std_rows, std_edges = [], []
    for label, mp in (("base_min6", 6), ("min0", 0)):
        R, taxa = residuals(mp)
        Z = R / R.std(0)
        gz, info = glasso_edges(Z, taxa)
        lv = louvain(gz)
        std_edges.append(gz.assign(setting=f"{label}_standardised"))
        r = {"setting": f"{label}_standardised", "taxa": len(taxa), "glasso_alpha": info["alpha"],
             "glasso_edges": len(gz), **coherence(gz), "louvain_nodes": lv["nodes"], "louvain_modules": lv["modules"],
             "louvain_Q": lv["Q"], "glasso_error": info["error"]}
        for t in KEY_TAXA:
            if t in taxa:
                e = gz[((gz.taxon_a == t) | (gz.taxon_b == t)) & (gz.partial_r.abs() > PARTIAL_THRESHOLD)]
                r[f"deg_gt0.05_{t}"] = len(e)
                r[f"partners_{t}"] = " ".join(f"{(x.taxon_b if x.taxon_a == t else x.taxon_a)}({x.partial_r:+.3f})"
                                              for x in e.itertuples())
        std_rows.append(r)
        print(f"{label} standardised: alpha {info['alpha']:.3g} edges {len(gz)} (>0.05 {r['edges_gt0.05']}) "
              f"louvain {lv['nodes']}n/{lv['modules']}m Q={lv['Q']:.3f}; key taxa with >0.05 edges: "
              f"{ {t: r[f'deg_gt0.05_{t}'] for t in KEY_TAXA if r.get(f'deg_gt0.05_{t}')} }", flush=True)
    pd.DataFrame(std_rows).to_csv(OUT / "within_glasso_standardised.csv", index=False)
    pd.concat(std_edges)[["setting", "taxon_a", "taxon_b", "partial_r"]].to_csv(
        OUT / "within_glasso_standardised_edges.csv", index=False)
    pd.DataFrame(stats_rows).to_csv(OUT / "within_glasso_degree_vs_variance_summary.csv", index=False)


def null_standardised(n):
    """Null for the standardised refit: the first n maxT shuffles of the baseline residuals, standardised."""
    R, taxa = residuals(6)
    V = R.shape[1]
    rng = np.random.default_rng(MAXT_SEED)
    rows = []
    for it in range(n):
        Rs = R.copy()
        for col in range(V):
            for ix in IDXS:
                Rs[ix, col] = Rs[rng.permutation(ix), col]
        ng, ni = glasso_edges(Rs / Rs.std(0), taxa)
        nl = louvain(ng)
        rows.append({"setting": "base_min6_standardised", "shuffle": it, "taxa": V, "glasso_alpha": ni["alpha"],
                     "glasso_edges": len(ng), **coherence(ng), "louvain_nodes": nl["nodes"],
                     "louvain_modules": nl["modules"], "louvain_Q": nl["Q"], "glasso_error": ni["error"]})
        print(f"null standardised shuffle {it}: alpha {ni['alpha']:.3g} edges {len(ng)} (>0.05 "
              f"{rows[-1]['edges_gt0.05']}, same-phylum {rows[-1]['same_phylum_frac']:.2f}, arch-arch "
              f"{rows[-1]['arch_arch_edges']}) Louvain {nl['nodes']}n/{nl['modules']}m Q {nl['Q']:.3f}", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "within_glasso_null_standardised.csv", index=False)


def edge_set(g, thr=0.0):
    return {frozenset((a, b)) for a, b, p in zip(g.taxon_a, g.taxon_b, g.partial_r) if abs(p) > thr}


def jaccard(a, b):
    return len(a & b) / len(a | b) if a | b else np.nan


def verify_baseline(res):
    """File-level comparison with the committed r232 outputs."""
    out = {}
    meta = json.load(open(HERE / "within_condition_meta_r232.json"))
    out["taxa_identical"] = list(res["taxa_rows"].taxon) == meta["taxa"]
    old = pd.read_csv(HERE / "all_candidate_pairs_r232.csv")
    mg = old.merge(res["cdf"], on=["taxon_a", "taxon_b"], suffixes=("_old", "_new"), how="outer", indicator=True)
    out["candidate_pairs_old_new_both"] = mg._merge.value_counts().to_dict()
    both = mg[mg._merge == "both"]
    for c in ("r_within", "r_marginal", "r_between", "perm_p", "q_BH"):
        out[f"max_abs_diff_{c}"] = float((both[f"{c}_old"] - both[f"{c}_new"]).abs().max())
    out["fdr_pass_equal"] = bool((both.fdr_pass_old == both.fdr_pass_new).all())
    out["effect_pass_equal"] = bool((both.effect_pass_old == both.effect_pass_new).all())
    og = pd.read_csv(HERE / "glasso_direct_edges_r232.csv")
    out["glasso_edges_old_new"] = (len(og), len(res["gdf"]))
    out["glasso_edge_sets_equal"] = edge_set(og) == edge_set(res["gdf"])
    m2 = og.merge(res["gdf"], on=["taxon_a", "taxon_b"], suffixes=("_old", "_new"))
    out["max_abs_diff_partial_r"] = float((m2.partial_r_old - m2.partial_r_new).abs().max())
    om = json.load(open(HERE / "within_condition_louvain_modules_r232.json"))
    out["louvain_modules_identical"] = om == res["lv"]["export"]
    out["louvain_on_committed_edge_csv_identical"] = louvain(og)["export"] == om
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="baseline only (verification)")
    ap.add_argument("--null", type=int, default=0, metavar="N",
                    help="only the extended Graphical-Lasso null: N shuffles for base_min6 and min0 "
                         "(needs the outputs of a full run)")
    ap.add_argument("--variance", action="store_true",
                    help="only the residual-variance / standardised-Glasso check (needs a full run)")
    ap.add_argument("--null-std", type=int, default=0, metavar="N",
                    help="only the null for the standardised refit: N shuffles of the baseline residuals")
    args = ap.parse_args()
    if args.null_std:
        sys.stdout = Tee(OUT / "within_run_null_standardised.log")
        null_standardised(args.null_std)
        return
    if args.variance:
        sys.stdout = Tee(OUT / "within_run_variance.log")
        variance_check()
        return
    if args.null:
        sys.stdout = Tee(OUT / "within_run_null.log")
        glasso_null_extended(args.null)
        return
    sys.stdout = Tee(OUT / ("within_run_quick.log" if args.quick else "within_run.log"))
    print(f"inputs: {len(LONG.columns) - 1} MAGs + unbinned_reads, 12 cores; MAG prevalence (cores > {ZERO}) "
          f"counts: {MAG_PREV.drop('unbinned_reads').value_counts().sort_index().to_dict()}")
    base = run("base_min6", 6, glasso_null=N_GLASSO_NULL)
    ver = verify_baseline(base)
    print("baseline verification vs committed r232 outputs:")
    for k_, v in ver.items():
        print(f"  {k_}: {v}")
    if args.quick:
        return
    results = [base]
    for mp in (4, 3, 2, 1, 0):
        results.append(run(f"min{mp}", mp, glasso_null=N_GLASSO_NULL if mp in (1, 0) else 0))
    for mp in (6, 4, 3, 2, 1, 0):
        results.append(run(f"min{mp}_2020ids", mp, defined_2020=True))
    p0 = base["summary"]["pseudocount"]
    for mp in (1, 0):
        results.append(run(f"min{mp}_pseudo_pinned", mp, pseudo_fixed=p0))
    # comparisons with the baseline
    be, be5 = edge_set(base["gdf"]), edge_set(base["gdf"], PARTIAL_THRESHOLD)
    for r in results:
        s = r["summary"]
        g, g5 = edge_set(r["gdf"]), edge_set(r["gdf"], PARTIAL_THRESHOLD)
        common = set(base["taxa_rows"].taxon) & set(r["taxa_rows"].taxon)
        within = lambda es: {e for e in es if e <= common}
        s["edge_jaccard_vs_base"] = jaccard(g, be)
        s["edge_jaccard_gt0.05_vs_base"] = jaccard(g5, be5)
        s["edge_jaccard_gt0.05_vs_base_common_taxa"] = jaccard(within(g5), within(be5))
        s["base_edges_gt0.05_retained"] = len(be5 & g5)
        nodes = sorted(set(base["lv"]["node_module"]) & set(r["lv"]["node_module"]))
        s["louvain_common_nodes"] = len(nodes)
        s["louvain_ARI_vs_base_common_nodes"] = (adjusted_rand_score(
            [base["lv"]["node_module"][n] for n in nodes], [r["lv"]["node_module"][n] for n in nodes])
            if len(nodes) > 1 else np.nan)
    summary = pd.DataFrame([r["summary"] for r in results])
    summary.to_csv(OUT / "within_settings_summary.csv", index=False)
    pd.DataFrame([row for r in results for row in r["key"]]).to_csv(OUT / "within_key_taxa_by_setting.csv", index=False)
    pd.concat([r["gdf"].assign(setting=r["summary"]["setting"]) for r in results])[
        ["setting", "taxon_a", "taxon_b", "partial_r"]].to_csv(OUT / "within_glasso_edges_by_setting.csv", index=False)
    json.dump({r["summary"]["setting"]: r["lv"]["export"] for r in results},
              open(OUT / "within_louvain_modules_by_setting.json", "w"), indent=1)
    pd.DataFrame([n for r in results for n in r["nulls"]]).to_csv(OUT / "within_glasso_null_shuffles.csv", index=False)
    # one row per taxon: prevalence, which settings admit it, Glasso degree per setting (r232-ID settings)
    tr = pd.concat([r["taxa_rows"] for r in results if not r["summary"]["defined_2020_ids_only"]
                    and not r["summary"]["pseudocount_pinned"]])
    wide = tr.pivot(index="taxon", columns="setting", values="glasso_degree")
    cols = [r["summary"]["setting"] for r in results if r["summary"]["setting"] in wide.columns]
    wide = wide[cols].add_prefix("glasso_degree_")
    mods = tr.pivot(index="taxon", columns="setting", values="louvain_module")[cols].add_prefix("module_")
    allt = sorted(set(taxon_of.values()))
    tm = pd.DataFrame(index=pd.Index(allt, name="taxon"))
    tm["kingdom"] = [KINGDOM[t] for t in allt]
    tm["roles"] = [ROLES[t] for t in allt]
    tm["max_MAG_prevalence"] = [int(max(MAG_PREV[x] for x in ids if taxon_of[x] == t)) for t in allt]
    tm["in_base_min6"] = tm.index.isin(base["taxa_rows"].taxon)
    tm = tm.join(wide).join(mods)
    tm.to_csv(OUT / "within_taxa_membership.csv")
    # key MAGs: prevalence per site and the threshold at which each passes the MAG filter
    rows = []
    for _, gr in guilds.iterrows():
        if not (ROLES.get(gr.taxon) or gr.taxon in NAMED):
            continue
        m_ = gr.mag
        pr = LONG[m_] > ZERO
        rows.append({"MAG_r232": gr.iterativeID_r232, "MAG_2020": gr.iterativeID_2020 if isinstance(gr.iterativeID_2020, str) else "",
                     "taxon": gr.taxon, "domain": gr.domain,
                     "flags": " ".join(lab for c, lab in ROLE_COLS if gr[c]) + (" archaea" if gr.archaea else ""),
                     "cores_gt_ZERO": int(pr.sum()),
                     **{f"cores_{s}": int(pr[[f"{s}_{x}" for x in "ABC"]].sum()) for s in ("R1", "R2", "R2A", "SF2")},
                     "max_rel_abund": float(LONG[m_].max()), "mean_rel_abund": float(LONG[m_].mean()),
                     **{f"passes_min{t}": bool(pr.sum() >= t) for t in (6, 4, 3, 2, 1, 0)}})
    pd.DataFrame(rows).sort_values(["taxon", "MAG_r232"]).to_csv(OUT / "within_key_MAG_prevalence.csv", index=False)
    json.dump(ver, open(OUT / "within_baseline_verification.json", "w"), indent=1, default=str)
    show = ["setting", "MAGs_kept", "taxa", "archaeal_taxa", "pseudocount", "frac_entries_le_ZERO",
            "share_residual_SS_from_le_ZERO_entries", "pairs_tested", "n_at_floor", "n_at_floor_needed_for_any_BH",
            "BH_pass", "min_q_BH", "effect_pass_obs", "effect_pass_null_mean", "maxT_95", "glasso_alpha",
            "glasso_edges", "glasso_edges_gt0.05", "louvain_nodes", "louvain_modules", "louvain_Q",
            "edge_jaccard_gt0.05_vs_base", "louvain_ARI_vs_base_common_nodes"]
    print(summary[show].to_string(index=False, float_format=lambda v: f"{v:.3g}"))


if __name__ == "__main__":
    main()
