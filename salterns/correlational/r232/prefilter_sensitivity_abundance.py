#!/usr/bin/env python
"""Sensitivity of the salterns correlational results to the abundance / prevalence prefilters that sit
OUTSIDE the within-condition (CLR-residual permutation + BH + Graphical Lasso) pipeline.

Filters examined (rule, code location, what it removes, what happens without it):
  (a) 2020 "defined iterativeID" filter: MAGs without a 2020 iterativeID are dropped
      (../analysis_2_marimo.py:1348-1350, :1795-1796, :2000-2001; analysis(2).ipynb cells 47/55).
      All 30 archaea had no 2020 ID.
  (b) 4-site "all v all" Spearman heatmap: MAG kept if its site mean (3 cores, DNA-weighted depth merge)
      is > 1e-5 at >= 3 of the 4 sites (../analysis_2_marimo.py:1344-1345; rerun_correlational_r232.py:283),
      applied per MAG BEFORE MAGs are summed into taxa ('reduced'); variants 'reduced' / 'reduced_defined'.
  (c) Figure S11 / within-condition matrix: MAG kept if > 1e-5 in >= 6 of the 12 depth-merged cores
      (_wc_run, ../analysis_2_marimo.py:1996; rerun_correlational_r232.py:386). The older cell
      (../analysis_2_marimo.py:1790-1791; analysis(2).ipynb cell 55) used >= 3 of 12 instead.
      Only the filter and the properties of the matrix it feeds are examined here; the permutation / BH /
      Graphical-Lasso stage belongs to the within-condition sensitivity analysis.
  (d) ZERO = 1e-5 in the guild-vs-CH4 and per-MAG-vs-CH4 correlations (rerun_correlational_r232.py:156-201)
      and the DNA-weighted D1/D2 merge (rerun_correlational_r232.py:86-89; ../analysis_2_marimo.py:65-66).

Outputs (all in ./prefilter_sensitivity/, prefixed abund_):
  abund_mag_prevalence.csv          one row per MAG: prevalence at 24 samples / 12 cores / 4 sites, pass/fail per filter
  abund_key_organisms.csv           the named organisms only (guild members, Methanolobus, TMA candidates, archaea, ...)
  abund_filter_summary.csv          per filter: MAGs and taxa kept/removed, archaea, guild members, binned reads kept
  abund_partial_taxa.csv            taxa that survive (b)/(c) only as a partial sum of their MAGs; guild colouring check
  abund_spearman_n4_exact_null.csv  exact permutation distribution of Spearman rho for n = 4 vs scipy's p
  abund_heatmap_variants.csv        4-site heatmap with / without the filter: m, |rho| = 1 pairs, BH (scipy vs exact p)
  abund_heatmap_rho_values.csv      the rho values that occur in each heatmap variant
  abund_s11_matrix_properties.csv   12-core S11-type matrix with / without the filter: zeros, ties, CLR, pseudocount
  abund_bh_vs_m.csv                 per-MAG vs CH4 BH-FDR on families of different size m (prefilter pools)
  abund_guild_sum_sensitivity.csv   guild-sum vs CH4 rho with members restricted by each filter, censoring, null pools
  abund_depth_merge_sensitivity.csv detection and guild rho under alternative D1/D2 merges
  abund_key_organism_prevalence.png key organisms: cores detected per site and pass/fail of each filter
  abund_n4_spearman_bh.png          n = 4 exact vs scipy p; BH discoveries vs number of tests m

Run (from this folder):  ~/Documents/py_venv/bin/python prefilter_sensitivity_abundance.py
No existing file is modified; everything is written to ./prefilter_sensitivity/.
"""
import collections
import itertools
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rerun_correlational_r232 as rr  # noqa: E402  (module-level code only loads the inputs)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

warnings.filterwarnings("ignore", category=stats.ConstantInputWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

OUT = HERE / "prefilter_sensitivity"
OUT.mkdir(exist_ok=True)
ZERO, ALPHA = rr.ZERO, 0.05
SITES = ["R1", "R2", "SF2", "R2A"]                       # unrestored, unrestored, restored, reference
UNRESTORED_SITES = {"R1", "R2"}
PRE = rr.PRE
CORES_BY_SITE = {s: [p for p in PRE if p.split("_")[0] == s] for s in SITES}
rel = rr.rel                                              # 311 x 24 (310 MAGs + unbinned_reads)
MAGS = rr.MAGS
ids, ids20 = rr.ids, rr.ids_2020
g = pd.read_csv(HERE / "guilds_r232.csv").set_index("mag")
assert list(g.index) == MAGS

tax232 = {m: ids[m].rsplit(".", 1)[0] for m in MAGS}
tax20 = {m: ids20[m].rsplit(".", 1)[0] for m in MAGS if m in ids20}
is_arch = g["archaea"].to_dict()


def save(df, name, **kw):
    df.to_csv(OUT / name, **kw)
    print("wrote", name, df.shape)


# =================================================================================================
# 0. prevalence tables
# =================================================================================================
AB12 = rr.merge_depths(rel)                               # DNA-weighted D1/D2 merge, 311 x 12
SITE4 = pd.DataFrame({s: AB12[CORES_BY_SITE[s]].mean(axis=1) for s in SITES})   # merge_triplicates
DET24, DET12, DET4 = rel > ZERO, AB12 > ZERO, SITE4 > ZERO

ROLE_COLS = [("mepn_consumer", "MePn phnJ"), ("mepn_2020", "MePn 2020 set"), ("betaine_reducer", "betaine grdH"),
             ("betaine_2020", "betaine 2020 set"), ("methanogen", "methanogen mcrA"),
             ("tma_producer", "TMA-route candidate"), ("archaea", "archaeon")]
NAMED = ("Roseovarius", "Albimonas", "Marivita", "Halomonas", "Vreelandella", "Halanaerobium", "Methanolobus")


def roles(m):
    r = [lab for col, lab in ROLE_COLS if g.loc[m, col]]
    r += [f"named: {n}" for n in NAMED if ids[m].startswith(n + ".")]
    return "; ".join(r)


prev = pd.DataFrame(index=MAGS)
prev["iterativeID_r232"] = [ids[m] for m in MAGS]
prev["iterativeID_2020"] = [ids20.get(m, "") for m in MAGS]
prev["taxon_r232"] = [tax232[m] for m in MAGS]
prev["domain"] = g["domain"]
prev["roles"] = [roles(m) for m in MAGS]
prev["n_samples_det_24"] = DET24.loc[MAGS].sum(axis=1)
prev["n_cores_det_12"] = DET12.loc[MAGS].sum(axis=1)
prev["n_unrestored_cores_det_6"] = DET12.loc[MAGS, CORES_BY_SITE["R1"] + CORES_BY_SITE["R2"]].sum(axis=1)
prev["n_sites_det_4"] = DET4.loc[MAGS].sum(axis=1)
for s in SITES:
    prev[f"cores_det_{s}"] = DET12.loc[MAGS, CORES_BY_SITE[s]].sum(axis=1)
for s in SITES:
    prev[f"site_mean_{s}"] = SITE4.loc[MAGS, s]
prev["max_rel_abund_12cores"] = AB12.loc[MAGS].max(axis=1)
prev["mean_rel_abund_12cores"] = AB12.loc[MAGS].mean(axis=1)
prev["only_unrestored_sites"] = [bool(DET12.loc[m, CORES_BY_SITE["R1"] + CORES_BY_SITE["R2"]].any()
                                      and not DET12.loc[m, CORES_BY_SITE["SF2"] + CORES_BY_SITE["R2A"]].any())
                                 for m in MAGS]
prev["n_sites_with_any_core_det"] = sum((prev[f"cores_det_{s}"] > 0).astype(int) for s in SITES)
prev["pass_a_defined2020"] = [m in ids20 for m in MAGS]
prev["pass_b_ge3of4sites"] = prev["n_sites_det_4"] >= 3
prev["pass_c_ge6of12cores"] = prev["n_cores_det_12"] >= 6
prev["pass_c_old_ge3of12cores"] = prev["n_cores_det_12"] >= 3
prev["in_2020_heatmap"] = prev.pass_a_defined2020 & prev.pass_b_ge3of4sites
prev["in_2020_S11"] = prev.pass_a_defined2020 & prev.pass_c_ge6of12cores
prev["pass_d_guild_sums"] = True                          # guild sums use every member (rerun:164)
x = AB12.loc[MAGS, PRE]
prev["pass_d_per_mag_nonconstant"] = (x.max(axis=1) - x.min(axis=1)) > 0   # rerun:192 np.ptp(x) > 0
save(prev.reset_index(names="mag"), "abund_mag_prevalence.csv", index=False)

key = prev[prev.roles != ""].copy()
order = {lab: i for i, (_, lab) in enumerate(ROLE_COLS)}
key["_o"] = [min([order.get(r, 9) for r in s.split("; ")]) for s in key.roles]
key = key.sort_values(["_o", "iterativeID_r232"]).drop(columns="_o")
save(key.reset_index(names="mag"), "abund_key_organisms.csv", index=False)

# =================================================================================================
# 1. filter summary (MAGs, taxa, archaea, guild members) + partial-sum taxa
# =================================================================================================
SETS = {lab: set(g.index[g[col]]) for col, lab in ROLE_COLS}
SETS["Roseovarius/Albimonas/Marivita"] = {m for m in MAGS if ids[m].split(".")[0] in ("Roseovarius", "Albimonas", "Marivita")}
SETS["Halomonas/Vreelandella"] = {m for m in MAGS if ids[m].split(".")[0] in ("Halomonas", "Vreelandella")}
SETS["Halanaerobium"] = {m for m in MAGS if ids[m].startswith("Halanaerobium.")}
binned = AB12.loc[MAGS, PRE]

FILTERS = [
    ("a: 2020 defined iterativeID", "analysis_2_marimo.py:1348-1350 (also :1795-1796, :2000-2001)", "pass_a_defined2020", "2020"),
    ("b: >=3 of 4 sites > 1e-5 (2020 heatmap, 2020 IDs + defined)", "analysis_2_marimo.py:1344-1345", "in_2020_heatmap", "2020"),
    ("b: >=3 of 4 sites > 1e-5 (r232 heatmap)", "rerun_correlational_r232.py:283", "pass_b_ge3of4sites", "r232"),
    ("c: >=6 of 12 cores > 1e-5 (2020 S11 / _wc_run, 2020 IDs + defined)", "analysis_2_marimo.py:1996 (+ defined drop :2000-2001)", "in_2020_S11", "2020"),
    ("c: >=6 of 12 cores > 1e-5 (r232 within-condition)", "rerun_correlational_r232.py:386", "pass_c_ge6of12cores", "r232"),
    ("c-old: >=3 of 12 cores > 1e-5 (older within-condition cell)", "analysis_2_marimo.py:1790-1791; analysis(2).ipynb cell 55", "pass_c_old_ge3of12cores", "r232"),
    ("d: guild sums (no prevalence filter)", "rerun_correlational_r232.py:164", "pass_d_guild_sums", "r232"),
    ("d: per-MAG vs CH4 (non-constant only)", "rerun_correlational_r232.py:192", "pass_d_per_mag_nonconstant", "r232"),
]
rows = []
for name, where, col, lab in FILTERS:
    keep = set(prev.index[prev[col]])
    tx_all = set(tax232.values()) if lab == "r232" else set(tax20.values())
    tmap = tax232 if lab == "r232" else tax20
    tx_keep = {tmap[m] for m in keep if m in tmap}
    r = {"filter": name, "code": where, "labels": lab, "MAGs_kept": len(keep), "MAGs_removed": 310 - len(keep),
         "taxa_kept": len(tx_keep), "taxa_total_in_labelling": len(tx_all),
         "archaeal_MAGs_kept": sum(is_arch[m] for m in keep),
         "binned_reads_kept_pct_mean12": 100 * binned.loc[sorted(keep)].sum().mean() / binned.sum().mean()}
    rem = prev.loc[[m for m in MAGS if m not in keep]]
    r["removed: detected only at R1/R2 (unrestored)"] = int(rem.only_unrestored_sites.sum())
    r["removed: detected in cores of 1 site only"] = int((rem.n_sites_with_any_core_det == 1).sum())
    r["removed: detected in no core"] = int((rem.n_cores_det_12 == 0).sum())
    r["removed: archaea"] = int(rem.domain.eq("Archaea").sum())
    for s, members in SETS.items():
        r[f"{s} kept"] = f"{len(members & keep)}/{len(members)}"
    rows.append(r)
summary = pd.DataFrame(rows)
save(summary, "abund_filter_summary.csv", index=False)

# partial-sum taxa: taxa that survive (b) or (c) with some of their MAGs dropped (filter is applied per MAG
# before collapse_to_taxa), and guild-coloured heatmap taxa whose gene-carrying MAGs were all dropped
part = []
for fcol, fname in (("pass_b_ge3of4sites", "b (r232 heatmap)"), ("pass_c_ge6of12cores", "c (r232 S11/within)")):
    for t, mags in pd.Series(tax232).groupby(pd.Series(tax232)).groups.items():
        mags = list(mags)
        kept = [m for m in mags if prev.loc[m, fcol]]
        if not kept:
            continue
        carriers = {col: [m for m in mags if g.loc[m, col]] for col in ("mepn_consumer", "betaine_reducer", "methanogen")}
        rec = {"filter": fname, "taxon_r232": t, "MAGs_in_taxon": len(mags), "MAGs_kept": len(kept),
               "partial_sum": len(kept) < len(mags),
               "abund_share_kept_pct": 100 * AB12.loc[kept, PRE].sum().sum() / AB12.loc[mags, PRE].sum().sum(),
               "kept_MAGs": " ".join(ids[m] for m in kept),
               "dropped_MAGs": " ".join(ids[m] for m in mags if m not in kept)}
        flags = []
        for col, lab in (("mepn_consumer", "MePn"), ("betaine_reducer", "betaine"), ("methanogen", "methanogen")):
            if carriers[col]:
                kc = [m for m in carriers[col] if m in kept]
                rec[f"{lab}_carriers_kept"] = f"{len(kc)}/{len(carriers[col])}"
                if not kc:
                    flags.append(f"coloured {lab} but no {lab} carrier survives")
            else:
                rec[f"{lab}_carriers_kept"] = ""
        rec["colouring_flag"] = "; ".join(flags)
        part.append(rec)
part = pd.DataFrame(part)
save(part[part.partial_sum | (part.colouring_flag != "")], "abund_partial_taxa.csv", index=False)

# =================================================================================================
# 2. n = 4 Spearman: exact permutation null vs scipy's t-approximation
# =================================================================================================
perm4 = np.array(list(itertools.permutations(range(4))))
base = np.arange(4)
rho_perm = np.array([stats.spearmanr(base, p)[0] for p in perm4])
vals = np.round(np.unique(np.round(rho_perm, 10)), 10)
null_rows = []
for v in vals:
    example = perm4[np.argmin(np.abs(rho_perm - v))]
    null_rows.append({"rho": v, "n_permutations_of_24": int(np.sum(np.isclose(rho_perm, v))),
                      "exact_p_two_sided": float(np.mean(np.abs(rho_perm) >= abs(v) - 1e-12)),
                      "exact_p_one_sided_upper": float(np.mean(rho_perm >= v - 1e-12)),
                      "scipy_spearmanr_p": float(stats.spearmanr(base, example)[1])})
n4 = pd.DataFrame(null_rows)
assert len(n4) == 11 and n4.loc[n4.rho == 1, "exact_p_two_sided"].item() == 2 / 24
assert n4.loc[n4.rho == 1, "scipy_spearmanr_p"].item() == 0.0
save(n4, "abund_spearman_n4_exact_null.csv", index=False)


def exact_perm_p(D):
    """Exact two-sided permutation p for every column pair of a 4 x N matrix (all 24 orders; ties allowed)."""
    Rk = stats.rankdata(D, axis=0)
    Z = Rk - Rk.mean(axis=0)
    nrm = np.sqrt((Z ** 2).sum(axis=0))
    with np.errstate(invalid="ignore", divide="ignore"):
        Z = Z / nrm
    obs = Z.T @ Z
    cnt = np.zeros_like(obs)
    for p in perm4:
        cnt += (np.abs(Z.T @ Z[p]) >= np.abs(obs) - 1e-9)
    P = cnt / len(perm4)
    P[~np.isfinite(obs)] = np.nan
    return obs, P


def spearman_matrix(D):
    """Spearman rho and scipy's two-sided t-approximation p for all column pairs of D (n x N).
    Same numbers as scipy.stats.spearmanr(D) (checked below); written out because spearmanr returns a
    scalar NaN instead of a matrix for some inputs with many constant columns."""
    n = D.shape[0]
    Rk = stats.rankdata(D, axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = np.corrcoef(Rk.T)
        rho = np.clip(rho, -1, 1)
        t = rho * np.sqrt((n - 2) / ((rho + 1.0) * (1.0 - rho)))
        p = 2 * stats.t.sf(np.abs(t), n - 2)
    p[np.isclose(np.abs(rho), 1)] = 0.0                     # scipy returns p = 0 at |rho| = 1
    return rho, p


_chk = np.random.default_rng(1).random((4, 30))
_r, _p = spearman_matrix(_chk)
_rs, _ps = stats.spearmanr(_chk)
assert np.allclose(_r, _rs) and np.allclose(_p, _ps, atol=1e-12)


def heatmap_variant(name, D, filtered_taxa=None, note=""):
    """D: sites x taxa. Reproduces heatmap() of rerun_correlational_r232.py (scipy p, BH over co-present pairs)."""
    N = D.shape[1]
    rho, p = spearman_matrix(D.values)
    _, pex = exact_perm_p(D.values)
    pres = (D > ZERO).values
    co = pres.T.astype(int) @ pres.astype(int)
    iu = np.triu_indices(N, 1)
    fin = np.isfinite(rho[iu])
    copres = (co[iu] >= 1) & fin
    r_t, p_t, pe_t = rho[iu][copres], p[iu][copres], pex[iu][copres]
    m = len(p_t)
    rej_sc = multipletests(p_t, alpha=ALPHA, method="fdr_bh")[0] if m else np.array([], bool)
    rej_ex = multipletests(pe_t, alpha=ALPHA, method="fdr_bh")[0] if m else np.array([], bool)
    pm1 = np.isclose(np.abs(r_t), 1)
    k = int(pm1.sum())
    const = int((np.ptp(D.values, axis=0) == 0).sum())
    ties = int(sum(len(np.unique(D.values[:, j])) < 4 for j in range(N)))
    cols = np.array(D.columns)
    a_idx, b_idx = iu[0][copres], iu[1][copres]
    extra = {}
    if filtered_taxa is not None:
        out_a = ~np.isin(cols[a_idx], list(filtered_taxa))
        out_b = ~np.isin(cols[b_idx], list(filtered_taxa))
        extra["pm1_pairs_with_a_taxon_failing_filter"] = int((pm1 & (out_a | out_b)).sum())
        extra["pm1_pairs_both_taxa_failing_filter"] = int((pm1 & out_a & out_b).sum())
    for t in ("Methanolobus", "Halanaerobium", "Roseovarius", "Haloquadratum", "Natrinema"):
        if t in cols:
            j = list(cols).index(t)
            sel = (a_idx == j) | (b_idx == j)
            extra[f"pm1_partners_{t}"] = int((pm1 & sel).sum())
            extra[f"tested_pairs_{t}"] = int(sel.sum())
        else:
            extra[f"pm1_partners_{t}"] = "not in matrix"
            extra[f"tested_pairs_{t}"] = 0
    rec = {"variant": name, "note": note, "n_taxa": N,
           "n_archaeal_taxa": int(sum(1 for c in cols if c in ARCH_TAXA)),
           "constant_columns": const, "columns_with_tied_site_values": ties,
           "all_pairs_C(N,2)": N * (N - 1) // 2, "pairs_with_finite_rho": int(fin.sum()),
           "m_tested_copresent_pairs": m, "pairs_abs_rho_eq_1": k, "scipy_p_eq_0": int((p_t == 0).sum()),
           "expected_abs_rho_1_if_independent (m/12)": m / 12,
           "BH_rejections_scipy_p": int(rej_sc.sum()),
           "BH_rejections_exact_p": int(rej_ex.sum()),
           "min_exact_p": float(np.nanmin(pe_t)) if m else np.nan,
           "BH_cut_needed_at_k=#pm1 (k*0.05/m)": (k * ALPHA / m) if m and k else np.nan,
           "BH_cut_needed_at_k=m (0.05)": ALPHA, **extra}
    vals_ = pd.Series(np.round(r_t, 6)).value_counts().sort_index()
    return rec, vals_


ARCH_TAXA = {tax232[m] for m in MAGS if is_arch[m]}


def to_taxa(df_sites_mags, tmap, keep_unlabelled=False):
    df = df_sites_mags.copy()
    if keep_unlabelled:   # 2020 'reduced' without 'defined': unlabelled MAG names lose their bin number
        cols = [tmap.get(c, c) if c in tmap else (".".join(c.split(".")[:-1]) or "unbinned") for c in df.columns]
        cols = [c if c else "unbinned" for c in cols]
    else:
        df = df.loc[:, [c for c in df.columns if c in tmap]]
        cols = [tmap[c] for c in df.columns]
    df.columns = cols
    return df.T.groupby(level=0).sum().T


S4 = SITE4.T                                               # 4 sites x 311
f_b = S4.loc[:, (S4 > ZERO).sum() >= 3]                    # the (b) rule on MAG columns
taxa_b_r232 = set(to_taxa(f_b, tax232).columns)
variants = [
    ("2020 heatmap reproduced (2020 IDs, defined, reduced, >=3/4 sites)", to_taxa(f_b, tax20), None,
     "check: 2020 pvals.txt has 4371 tests, 1006 with p = 0, 1006 BH rejections"),
    ("2020 'reduced' without 'defined' (>=3/4 sites)", to_taxa(f_b, tax20, keep_unlabelled=True), None,
     "unlabelled archaea become bin-name stems; unbinned_reads becomes 'unbinned'"),
    ("r232 heatmap as published in r232/ (>=3/4 sites)", to_taxa(f_b, tax232), None,
     "check: heatmap_r232.log reports 109 taxa, 1285 FDR-significant pairs"),
    ("r232, NO prevalence filter", to_taxa(S4, tax232), taxa_b_r232,
     "sub-1e-5 site means kept as measured (mapping-noise level, mostly non-zero)"),
    ("r232, NO prevalence filter, values <= 1e-5 censored to 0", to_taxa(S4.where(S4 > ZERO, 0.0), tax232), taxa_b_r232,
     "censoring creates ties and constant (all-zero) columns"),
    ("2020 IDs (defined), NO prevalence filter", to_taxa(S4, tax20), set(to_taxa(f_b, tax20).columns), ""),
]
hv_rows, hv_vals = [], {}
for name, D, ftaxa, note in variants:
    rec, vc = heatmap_variant(name, D, ftaxa, note)
    hv_rows.append(rec)
    hv_vals[name] = vc
hv = pd.DataFrame(hv_rows)
assert hv.loc[0, "m_tested_copresent_pairs"] == 4371 and hv.loc[0, "BH_rejections_scipy_p"] == 1006
assert hv.loc[2, "n_taxa"] == 109 and hv.loc[2, "BH_rejections_scipy_p"] == 1285
save(hv, "abund_heatmap_variants.csv", index=False)
save(pd.DataFrame(hv_vals).fillna(0).astype(int).rename_axis("rho"), "abund_heatmap_rho_values.csv")
reduced_nodef_cols = list(to_taxa(f_b, tax20, keep_unlabelled=True).columns)
pseudo_taxa = [c for c in reduced_nodef_cols if c.startswith("Salt_Pond") or c == "unbinned"]

# =================================================================================================
# 3. S11-type matrix (12 cores x taxa) with / without the prevalence filter
# =================================================================================================
from skbio.stats.composition import clr  # noqa: E402

COND = np.array([p.split("_")[0] for p in PRE])


def s11_matrix(min_cores, tmap):
    long = AB12.loc[MAGS, PRE].T                           # 12 x 310 (unbinned excluded as in 'defined')
    if min_cores:
        long = long.loc[:, (long > ZERO).sum() >= min_cores]
    return to_taxa(long, tmap)


def residuals(X):
    R = X.copy()
    for c in np.unique(COND):
        R[COND == c] -= R[COND == c].mean(axis=0, keepdims=True)
    return R


def corr_upper(M):
    C = np.corrcoef(M.T)
    return C


s11_rows = []
r_within_sets = {}
for lab, tmap in (("r232", tax232), ("2020 (defined)", tax20)):
    for min_cores, fname in ((6, ">=6 of 12 (S11 caption, _wc_run)"), (3, ">=3 of 12 (older cell)"), (0, "none")):
        L = s11_matrix(min_cores, tmap)
        v = L.values
        nz = v[v > 0]
        pc = nz.min() / 2
        X = clr(v + pc)
        R = residuals(X)
        below = v <= ZERO
        ss = (X ** 2)
        V = L.shape[1]
        pres = (L > ZERO).values.astype(int)
        co = pres.T @ pres
        iu = np.triu_indices(V, 1)
        # censored variant: <= 1e-5 -> 0, same pseudocount rule
        vc = np.where(v > ZERO, v, 0.0)
        pcc = vc[vc > 0].min() / 2
        Xc = clr(vc + pcc)
        Rc = residuals(Xc)
        Cc = np.corrcoef(Rc.T)
        allzero = (vc == 0).all(axis=0)
        pm1c = np.isclose(np.abs(Cc[iu]), 1, atol=1e-9)
        Cw = np.corrcoef(R.T)
        r_within_sets[(lab, fname)] = pd.DataFrame(Cw, index=L.columns, columns=L.columns)
        old = clr(v + 0.5) if min_cores == 3 else None    # the older cell's +0.5 pseudocount
        s11_rows.append({
            "labels": lab, "prevalence_filter": fname, "n_taxa": V,
            "n_archaeal_taxa": int(sum(c in ARCH_TAXA for c in L.columns)) if lab == "r232" else 0,
            "cells": v.size, "frac_cells_le_1e-5": below.mean(), "exact_zero_cells": int((v == 0).sum()),
            "min_nonzero": nz.min(), "pseudocount_half_min_nonzero": pc,
            "CLR_min": X.min(), "CLR_max": X.max(), "CLR_range": X.max() - X.min(),
            "share_of_sum_CLR2_from_cells_le_1e-5": ss[below].sum() / ss.sum(),
            "taxa_with_tied_values": int(sum(len(np.unique(v[:, j])) < 12 for j in range(V))),
            "constant_taxa_raw": int((np.ptp(v, axis=0) == 0).sum()),
            "NaN_r_within_taxa_raw": int(np.isnan(np.diag(Cw)).sum()),
            "pairs_C(V,2)": V * (V - 1) // 2, "pairs_cooccurring_ge6": int((co[iu] >= 6).sum()),
            "censored_taxa_all_below_1e-5": int(allzero.sum()),
            "censored_taxa_with_ties": int(sum(len(np.unique(vc[:, j])) < 12 for j in range(V))),
            "censored_pseudocount": pcc,
            "censored_pairs_abs_r_within_eq_1": int(pm1c.sum()),
            "censored_NaN_r_within_taxa": int(np.isnan(np.diag(Cc)).sum()),
            "older_cell_+0.5_CLR_range": (old.max() - old.min()) if old is not None else np.nan,
        })
s11 = pd.DataFrame(s11_rows)
assert s11.loc[(s11.labels == "r232") & (s11.prevalence_filter.str.startswith(">=6")), "n_taxa"].item() == 169
assert s11.loc[(s11.labels == "2020 (defined)") & (s11.prevalence_filter.str.startswith(">=6")), "n_taxa"].item() == 140

# CLR-reference sensitivity: r_within among the S11 (>=6/12) taxa when CLR is taken over the filtered set
# (as published) vs over all taxa (no prefilter)
clr_cmp = []
for lab in ("r232", "2020 (defined)"):
    A = r_within_sets[(lab, ">=6 of 12 (S11 caption, _wc_run)")]
    B = r_within_sets[(lab, "none")].loc[A.index, A.columns]
    iu = np.triu_indices(len(A), 1)
    a, b = A.values[iu], B.values[iu]
    clr_cmp.append({"labels": lab, "pairs": len(a), "pearson_r_between_versions": np.corrcoef(a, b)[0, 1],
                    "mean_abs_diff": np.abs(a - b).mean(), "max_abs_diff": np.abs(a - b).max(),
                    "pairs_crossing_abs_r_0.7": int(((np.abs(a) >= 0.7) != (np.abs(b) >= 0.7)).sum()),
                    "pairs_abs_r_ge_0.7_filtered": int((np.abs(a) >= 0.7).sum()),
                    "pairs_abs_r_ge_0.7_unfiltered_CLR": int((np.abs(b) >= 0.7).sum())})
clr_cmp = pd.DataFrame(clr_cmp)
s11 = s11.merge(clr_cmp.add_prefix("CLRref_").rename(columns={"CLRref_labels": "labels"}).assign(
    prevalence_filter=">=6 of 12 (S11 caption, _wc_run)"), on=["labels", "prevalence_filter"], how="left")
save(s11, "abund_s11_matrix_properties.csv", index=False)

# =================================================================================================
# 4. BH-FDR vs number of tests m: per-MAG vs CH4 (12 cores), families defined by each prefilter
# =================================================================================================
pm = pd.read_csv(HERE / "mag_methane_correlations_r232.csv").set_index("mag")
# recompute to be sure the stored p-values are what the current inputs give
for basis, tag in (("raw", "raw"), ("DNA-normalised", "dna")):
    for m in MAGS[:25]:
        r_, p_ = stats.spearmanr(rr.AB[basis].loc[m, PRE].values, rr.CH4[basis][PRE].values)
        assert np.isclose(p_, pm.loc[m, f"p_{tag}"])
guild_union = set().union(*[SETS[lab] for _, lab in ROLE_COLS])
POOLS = [("all 310 MAGs (as run)", set(MAGS)),
         ("a: 2020-defined (280 bacteria)", set(prev.index[prev.pass_a_defined2020])),
         ("b: >=3 of 4 sites", set(prev.index[prev.pass_b_ge3of4sites])),
         ("c: >=6 of 12 cores", set(prev.index[prev.pass_c_ge6of12cores])),
         ("c-old: >=3 of 12 cores", set(prev.index[prev.pass_c_old_ge3of12cores])),
         ("guild members only (union of key sets)", guild_union)]


def bh(p):
    p = np.asarray(p)
    m = len(p)
    rej, q, _, _ = multipletests(p, alpha=ALPHA, method="fdr_bh")
    srt = np.sort(p)
    ok = np.nonzero(srt <= (np.arange(1, m + 1) * ALPHA / m))[0]
    cut = srt[ok[-1]] if len(ok) else np.nan
    return rej, q, cut, (ok[-1] + 1) * ALPHA / m if len(ok) else ALPHA / m


bh_rows = []
for tag in ("raw", "dna"):
    for pname, pool in POOLS:
        sub = pm.loc[sorted(pool)]
        sub = sub[sub[f"p_{tag}"].notna()]
        rej, q, cut, kcut = bh(sub[f"p_{tag}"].values)
        sig = set(sub.index[rej])
        r = {"basis": tag, "family": pname, "m_tests": len(sub), "BH_discoveries_q<0.05": int(rej.sum()),
             "nominal_p<0.05": int((sub[f"p_{tag}"] < 0.05).sum()), "largest_rejected_p": cut,
             "BH_line_at_last_rejection (k*0.05/m)": kcut, "BH_line_at_k=1 (0.05/m)": ALPHA / len(sub)}
        for lab in ("MePn phnJ", "betaine grdH", "methanogen mcrA", "archaeon"):
            s = SETS[lab]
            r[f"{lab}: significant/in family/total"] = f"{len(s & sig)}/{len(s & set(sub.index))}/{len(s)}"
        bh_rows.append(r)
bhm = pd.DataFrame(bh_rows)
assert bhm.loc[(bhm.basis == "raw") & (bhm.family.str.startswith("all 310")), "BH_discoveries_q<0.05"].item() == 254

# BH under a pure-null family: how many of the +-1 pairs survive when the same k signals sit in families of size m
synth = []
p_sig = np.full(5, 1e-4)
for m in (5, 10, 50, 100, 310, 1000, 5886, 30000):
    p = np.concatenate([p_sig, np.linspace(1 / m, 1, max(m - 5, 0))]) if m > 5 else p_sig[:m]
    synth.append({"basis": "illustration", "family": f"5 tests at p = 1e-4 + {max(m - 5, 0)} null tests",
                  "m_tests": m, "BH_discoveries_q<0.05": int(multipletests(p, alpha=ALPHA, method='fdr_bh')[0].sum()),
                  "BH_line_at_k=1 (0.05/m)": ALPHA / m, "BH_line_at_last_rejection (k*0.05/m)": 5 * ALPHA / m})
bhm = pd.concat([bhm, pd.DataFrame(synth)], ignore_index=True)
save(bhm, "abund_bh_vs_m.csv", index=False)

# =================================================================================================
# 5. guild sums vs CH4: every member (as run) vs members restricted by each filter, and censoring
# =================================================================================================
N_NULL = 20000                  # 4000 in rerun_correlational_r232.py; 20000 keeps the Monte-Carlo SE <= 0.35 points
rng = np.random.default_rng(0)


def draw(n_pool, k):
    """N_NULL random k-subsets (without replacement) of range(n_pool), vectorised."""
    return np.argsort(rng.random((N_NULL, n_pool)), axis=1)[:, :k]
AB_cens = {b: rr.AB[b].where(rr.AB["raw"] > ZERO, 0.0) for b in rr.AB}


def vec_spearman(S, y):
    """Spearman rho of each row of S (draws x cores) with y."""
    Rs = stats.rankdata(S, axis=1)
    Rs = Rs - Rs.mean(axis=1, keepdims=True)
    ry = stats.rankdata(y) - stats.rankdata(y).mean()
    with np.errstate(invalid="ignore", divide="ignore"):
        return (Rs @ ry) / np.sqrt((Rs ** 2).sum(axis=1) * (ry ** 2).sum())


ref = pd.read_csv(HERE / "guild_methane_correlations_r232.csv")
G_ROWS = []
GUILDS = [("mepn_consumer", "MePn phnJ"), ("mepn_2020", "MePn 2020 set"), ("betaine_reducer", "betaine grdH"),
          ("betaine_2020", "betaine 2020 set"), ("tma_producer", "TMA-route candidate"),
          ("methanogen", "methanogen mcrA"), ("archaea", "archaeon")]
VARIANTS = [("all members (as run)", "all", False), ("a: 2020-defined members", "pass_a_defined2020", False),
            ("b: members >=3 of 4 sites", "pass_b_ge3of4sites", False),
            ("c: members >=6 of 12 cores", "pass_c_ge6of12cores", False),
            ("all members, values <= 1e-5 censored to 0", "all", True)]
for col, lab in GUILDS:
    members_all = list(g.index[g[col]])
    tot = AB12.loc[members_all, PRE].sum().sum()
    for vname, fcol, cens in VARIANTS:
        pool = MAGS if fcol == "all" else list(prev.index[prev[fcol]])
        members = [m for m in members_all if m in pool]
        for basis in rr.AB:
            A = (AB_cens if cens else rr.AB)[basis]
            for subset, cores in (("12 cores", PRE), ("6 unrestored", rr.UNRESTORED)):
                y = rr.CH4[basis][cores].values
                if members:
                    rho, p = stats.spearmanr(A.loc[members, cores].sum().values, y)
                    M = A.loc[pool, cores].values
                    S = M[draw(len(pool), len(members))].sum(axis=1)
                    nul = vec_spearman(S, y)
                    nul = nul[np.isfinite(nul)]
                    nf = (nul >= rho - 1e-12).mean() if np.isfinite(rho) else np.nan
                else:
                    rho = p = nf = np.nan
                rec = {"guild": lab, "variant": vname, "n_members": len(members), "n_members_all": len(members_all),
                       "abund_share_kept_pct": 100 * AB12.loc[members, PRE].sum().sum() / tot if members else 0.0,
                       "basis": basis, "cores": subset, "spearman_rho": rho, "p": p,
                       "null_frac_ge_rho (pool = MAGs passing the same filter)": nf,
                       "dropped_members": " ".join(ids[m] for m in members_all if m not in members)}
                if vname.startswith("all members (as run)"):
                    rr_row = ref[(ref.guild == col) & (ref.basis == basis) & (ref.cores == subset)]
                    rec["rho_in_guild_methane_correlations_r232.csv"] = rr_row.spearman_rho.item()
                    rec["null_in_guild_methane_correlations_r232.csv"] = rr_row.null_frac_ge_rho.item()
                    assert np.isclose(rr_row.spearman_rho.item(), rho)
                G_ROWS.append(rec)
# the random-subset null itself depends on which MAGs it may draw: same as-run guild sums (all members),
# null drawn only from MAGs passing a filter
for col, lab in GUILDS:
    members = list(g.index[g[col]])
    for pname, fcol in (("all 310 MAGs (as run)", None), ("a: 2020-defined MAGs", "pass_a_defined2020"),
                        ("b: MAGs >=3 of 4 sites", "pass_b_ge3of4sites"), ("c: MAGs >=6 of 12 cores", "pass_c_ge6of12cores")):
        pool = MAGS if fcol is None else list(prev.index[prev[fcol]])
        for basis in rr.AB:
            y = rr.CH4[basis][PRE].values
            rho = stats.spearmanr(rr.AB[basis].loc[members, PRE].sum().values, y)[0]
            M = rr.AB[basis].loc[pool, PRE].values
            nul = vec_spearman(M[draw(len(pool), len(members))].sum(axis=1), y)
            nul = nul[np.isfinite(nul)]
            G_ROWS.append({"guild": lab, "variant": f"all members; null pool = {pname}", "n_members": len(members),
                           "n_members_all": len(members), "abund_share_kept_pct": 100.0, "basis": basis,
                           "cores": "12 cores", "spearman_rho": rho, "p": np.nan,
                           "null_frac_ge_rho (pool = MAGs passing the same filter)": (nul >= rho - 1e-12).mean(),
                           "dropped_members": "", "null_pool_size": len(pool)})
gs = pd.DataFrame(G_ROWS)
save(gs, "abund_guild_sum_sensitivity.csv", index=False)

# =================================================================================================
# 6. depth merge: DNA-weighted (as run) vs alternatives
# =================================================================================================
d = rr.dna
assert np.allclose(rr.AB["DNA-normalised"].values,
                   pd.DataFrame({p: (rel[f"{p}_D1"] + rel[f"{p}_D2"]) / (d[f"{p}_D1"] + d[f"{p}_D2"]) for p in PRE}).values)
MERGES = {
    "DNA-weighted (as run)": lambda a, b, d1, d2: (a * d1 + b * d2) / (d1 + d2),
    "plain mean": lambda a, b, d1, d2: (a + b) / 2,
    "thickness 5:10 (D1 = 0-5, D2 = 5-15 cm; inferred)": lambda a, b, d1, d2: (5 * a + 10 * b) / 15,
    "DNA x thickness": lambda a, b, d1, d2: (5 * d1 * a + 10 * d2 * b) / (5 * d1 + 10 * d2),
    "max of the two depths (detection upper bound)": lambda a, b, d1, d2: np.maximum(a, b),
}
w_d1 = pd.Series({p: d[f"{p}_D1"] / (d[f"{p}_D1"] + d[f"{p}_D2"]) for p in PRE})
dm_rows = []
ref_det = None
for mname, f in MERGES.items():
    M = pd.DataFrame({p: f(rel.loc[MAGS, f"{p}_D1"], rel.loc[MAGS, f"{p}_D2"], d[f"{p}_D1"], d[f"{p}_D2"]) for p in PRE})
    det = M > ZERO
    s4 = pd.DataFrame({s: M[CORES_BY_SITE[s]].mean(axis=1) for s in SITES}) > ZERO
    pb, pc = s4.sum(axis=1) >= 3, det.sum(axis=1) >= 6
    if ref_det is None:
        ref_det, ref_pb, ref_pc = det, pb, pc
    flips_b = [ids[m] for m in MAGS if pb[m] != ref_pb[m] and m in guild_union]
    flips_c = [ids[m] for m in MAGS if pc[m] != ref_pc[m] and m in guild_union]
    r = {"merge": mname, "MAG_core_cells_detected_of_3720": int(det.values.sum()),
         "cells_flipping_vs_DNA_weighted": int((det != ref_det).values.sum()),
         "MAGs_pass_b_ge3of4": int(pb.sum()), "MAGs_pass_c_ge6of12": int(pc.sum()),
         "key_organisms_flipping_b": " ".join(flips_b), "key_organisms_flipping_c": " ".join(flips_c)}
    for col, lab in GUILDS:
        mem = list(g.index[g[col]])
        r[f"rho_raw_12cores {lab}"] = stats.spearmanr(M.loc[mem, PRE].sum().values, rr.CH4["raw"][PRE].values)[0]
    dm_rows.append(r)
dm = pd.DataFrame(dm_rows)
# cells where one depth is detected but the DNA-weighted merge falls to <= 1e-5
D1, D2 = rel.loc[MAGS, [f"{p}_D1" for p in PRE]].values, rel.loc[MAGS, [f"{p}_D2" for p in PRE]].values
diluted = ((D1 > ZERO) | (D2 > ZERO)) & ~(AB12.loc[MAGS, PRE].values > ZERO)
dm["note"] = ""
dm.loc[0, "note"] = (f"{int(diluted.sum())} MAG x core cells have one depth > 1e-5 but the DNA-weighted merge <= 1e-5; "
                     f"D1 weight per core ranges {w_d1.min():.3f}-{w_d1.max():.3f}; "
                     f"on the DNA-normalised basis the merge reduces to (rel_D1 + rel_D2)/(DNA_D1 + DNA_D2)")
save(dm, "abund_depth_merge_sensitivity.csv", index=False)

# =================================================================================================
# 7. figures
# =================================================================================================
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SEQ = ["#fcfcfb", "#b7d3f6", "#5598e7", "#184f95"]      # 0..3 cores detected (sequential blue)
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})

# --- figure 1: key organisms x filters
fig_rows = key.copy()
group = []
for r in fig_rows.roles:
    if "MePn" in r:
        group.append("MePn consumers (phnJ and/or 2020 set)")
    elif "betaine" in r:
        group.append("betaine reducers (grdH and/or 2020 set)")
    elif "methanogen" in r:
        group.append("methanogen")
    elif "TMA" in r:
        group.append("TMA-route candidates")
    elif "archaeon" in r:
        group.append("other archaea")
    else:
        group.append("other named taxa")
fig_rows["group"] = group
gorder = ["MePn consumers (phnJ and/or 2020 set)", "betaine reducers (grdH and/or 2020 set)", "methanogen",
          "TMA-route candidates", "other named taxa", "other archaea"]
fig_rows["_g"] = [gorder.index(x) for x in fig_rows.group]
fig_rows = fig_rows.sort_values(["_g", "n_cores_det_12", "iterativeID_r232"], ascending=[True, False, True])
n = len(fig_rows)
fig, ax = plt.subplots(figsize=(9.2, 0.205 * n + 1.9), dpi=170)
ycur, yt, ylab, gpos = 0, [], [], []
for gname in gorder:
    sub = fig_rows[fig_rows.group == gname]
    if sub.empty:
        continue
    gpos.append((ycur - 0.2, gname))
    ycur += 1
    for m, r in sub.iterrows():
        for j, s in enumerate(SITES):
            c = int(r[f"cores_det_{s}"])
            ax.add_patch(plt.Rectangle((j - 0.45, ycur - 0.42), 0.9, 0.84, fc=SEQ[c], ec="#c3c2b7", lw=0.4))
            ax.text(j, ycur, str(c), ha="center", va="center", fontsize=6.5, color=INK if c < 2 else "#fcfcfb")
        ax.text(4.0, ycur, f"{int(r.n_cores_det_12)}", ha="center", va="center", fontsize=7, color=INK)
        for k, fc in enumerate(["pass_a_defined2020", "pass_b_ge3of4sites", "pass_c_ge6of12cores", "pass_c_old_ge3of12cores"]):
            ok = bool(r[fc])
            ax.scatter(5.2 + k, ycur, s=34, marker="o", facecolor=INK if ok else "white", edgecolor=INK if ok else MUTED,
                       lw=1.0, zorder=3)
        lab = r.iterativeID_r232 + ("  (2020: " + r.iterativeID_2020 + ")" if r.iterativeID_2020 and
                                    r.iterativeID_2020 != r.iterativeID_r232 else "")
        if g.loc[m, "mepn_2020"] and not g.loc[m, "mepn_consumer"]:
            lab += "  [2020 set only]"
        if g.loc[m, "betaine_2020"] and not g.loc[m, "betaine_reducer"]:
            lab += "  [2020 set only]"
        yt.append(ycur); ylab.append(lab)
        ycur += 1
    ycur += 0.4
for y0, gname in gpos:
    ax.text(-0.75, y0 + 0.55, gname, ha="right", va="center", fontsize=8, fontweight="bold", color=INK)
ax.set_yticks(yt); ax.set_yticklabels(ylab, fontsize=6.8)
ax.set_ylim(ycur, -1.6)
ax.set_xlim(-0.6, 8.6)
tops = SITES + ["cores\n/12"] + ["(a) 2020\nID", "(b) ≥3/4\nsites", "(c) ≥6/12\ncores", "≥3/12\n(old)"]
ax.set_xticks(list(range(4)) + [4.0] + [5.2, 6.2, 7.2, 8.2]); ax.set_xticklabels(tops, fontsize=7.2)
ax.xaxis.tick_top()
ax.tick_params(length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
ax.text(-0.6, -1.25, "cores with rel. abund. > 1e-5 per site (of 3)", fontsize=7, color=INK2, ha="left", va="bottom")
ax.text(5.2, -1.25, "filled = kept, open = removed", fontsize=7, color=INK2, ha="left", va="bottom")
fig.suptitle("Which key organisms each prefilter keeps (DNA-weighted depth merge; detection = rel. abund. > 1e-5)",
             fontsize=9.5, x=0.01, ha="left", color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.985))
fig.savefig(OUT / "abund_key_organism_prevalence.png", bbox_inches="tight", facecolor="#fcfcfb")
plt.close(fig)
print("wrote abund_key_organism_prevalence.png")

# --- figure 2: n = 4 exact vs scipy p; BH discoveries vs m
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.3), dpi=170, gridspec_kw={"width_ratios": [1, 1.15]})
pos = n4[n4.rho >= 0].sort_values("rho")
xx = np.arange(len(pos))
a1.bar(xx - 0.19, pos.exact_p_two_sided, width=0.36, color=BLUE, label="exact permutation p (24 orders)")
a1.bar(xx + 0.19, pos.scipy_spearmanr_p, width=0.36, color=ORANGE, label="scipy.stats.spearmanr p (t approx.)")
a1.axhline(0.05, color=INK2, lw=1, ls="--")
a1.text(-0.45, 0.07, "p = 0.05", ha="left", fontsize=7.5, color=INK2,
        bbox=dict(boxstyle="square,pad=0.15", fc="#fcfcfb", ec="none"))
a1.set_xticks(xx); a1.set_xticklabels([f"{v:.1f}" for v in pos.rho], fontsize=8)
a1.set_xlabel("|Spearman ρ| (n = 4 sites: the only attainable values)", fontsize=8.5)
a1.set_ylabel("two-sided p", fontsize=8.5)
a1.text(xx[-1] + 0.19, 0.02, "0", ha="center", fontsize=8, color=ORANGE, fontweight="bold")
a1.text(xx[-1] - 0.19, 2 / 24 + 0.02, "1/12", ha="center", fontsize=8, color=BLUE, fontweight="bold")
a1.set_title("a  n = 4: no ρ reaches p < 0.05 exactly; scipy gives p = 0 at |ρ| = 1", fontsize=9, loc="left")
a1.legend(fontsize=7.5, frameon=False, loc="upper right")
a1.grid(axis="y", color=GRID, lw=0.6); a1.set_axisbelow(True)
for sp in ("top", "right"):
    a1.spines[sp].set_visible(False)
# panel b: discoveries vs m
hvp = hv.iloc[[0, 2, 3, 4]].copy()
short = [("2020, ≥3/4 sites: 1006", (-9, 4), "right"), ("r232, ≥3/4 sites: 1285", (9, -3), "left"),
         ("r232, no filter: 3930", (9, 2), "left"), ("r232, no filter, ≤1e-5 set to 0: 2025", (9, -4), "left")]
for (_, r), (s_, off, ha) in zip(hvp.iterrows(), short):
    assert s_.endswith(str(int(r.BH_rejections_scipy_p)))
    a2.scatter(r.m_tested_copresent_pairs, r.BH_rejections_scipy_p, s=70, color=ORANGE, edgecolor="#fcfcfb", lw=1.5, zorder=3)
    a2.scatter(r.m_tested_copresent_pairs, max(r.BH_rejections_exact_p, 0.6), s=60, marker="s", color=BLUE,
               edgecolor="#fcfcfb", lw=1.5, zorder=3)
    a2.annotate(s_, (r.m_tested_copresent_pairs, r.BH_rejections_scipy_p), xytext=off, textcoords="offset points",
                fontsize=6.8, color=INK2, ha=ha, va="center")
for tag, mk, colr, dx in (("raw", "^", AQUA, 0.97), ("dna", "v", INK2, 1.03)):
    sub = bhm[(bhm.basis == tag)]
    a2.scatter(sub.m_tests * dx, sub["BH_discoveries_q<0.05"].clip(lower=0.6), s=55, marker=mk, color=colr, zorder=3,
               edgecolor="#fcfcfb", lw=0.8)
for lab_, mt in (("guilds only", 72), ("(b)", 164), ("(c), ≥3/12, (a), all", 245)):
    yv = bhm[(bhm.basis == "raw") & (bhm.m_tests == mt)]["BH_discoveries_q<0.05"].item()
    a2.annotate(lab_, (mt, yv), xytext=(0, 10), textcoords="offset points", fontsize=6.5, color=INK2,
                ha="center" if mt < 200 else "left")
a2.text(0.02, 0.97, "heatmap labels: number of |ρ| = 1 pairs (= BH discoveries on scipy p)", transform=a2.transAxes, fontsize=7,
        color=INK2, va="top")
a2.set_xscale("log"); a2.set_yscale("log")
a2.set_ylim(0.4, 2e4); a2.set_xlim(40, 2.5e5)
a2.text(0.02, 0.03, "values at 0.6 = zero discoveries", transform=a2.transAxes, fontsize=7, color=MUTED)
hdl = [plt.Line2D([], [], ls="", marker="o", color=ORANGE, ms=7, label="4-site heatmap pairs, BH on scipy p"),
       plt.Line2D([], [], ls="", marker="s", color=BLUE, ms=7, label="4-site heatmap pairs, BH on exact p (0 in every variant)"),
       plt.Line2D([], [], ls="", marker="^", color=AQUA, ms=6, label="per-MAG vs CH₄ (raw), one family per prefilter"),
       plt.Line2D([], [], ls="", marker="v", color=INK2, ms=6, label="per-MAG vs CH₄ (DNA-norm.), same families")]
a2.legend(handles=hdl, fontsize=7, frameon=False, loc="lower left", bbox_to_anchor=(0.0, 0.08))
a2.set_xlabel("number of tests m in the BH family", fontsize=8.5)
a2.set_ylabel("BH discoveries at q < 0.05", fontsize=8.5)
a2.set_title("b  BH discoveries vs m: the heatmap count is the |ρ| = 1 count, whatever m", fontsize=9, loc="left")
a2.grid(color=GRID, lw=0.6, which="major"); a2.set_axisbelow(True)
for sp in ("top", "right"):
    a2.spines[sp].set_visible(False)
fig.tight_layout()
fig.savefig(OUT / "abund_n4_spearman_bh.png", bbox_inches="tight", facecolor="#fcfcfb")
plt.close(fig)
print("wrote abund_n4_spearman_bh.png")

# =================================================================================================
# console summary (numbers quoted in abund_REPORT.md)
# =================================================================================================
pd.set_option("display.width", 250, "display.max_columns", 40, "display.max_colwidth", 80)
print("\n== filter summary ==\n", summary.drop(columns="code").to_string(index=False))
print("\n== heatmap variants ==\n", hv.T.to_string())
print("\n2020 'reduced' (no 'defined') pseudo-taxa:", len(pseudo_taxa), pseudo_taxa[:6])
print("\n== S11 matrix ==\n", s11.T.to_string())
print(s11[["labels", "prevalence_filter", "min_nonzero", "pseudocount_half_min_nonzero", "censored_pseudocount"]].to_string())
print("\n== BH vs m ==\n", bhm.to_string(index=False))
print("\n== depth merge ==\n", dm.T.to_string())
