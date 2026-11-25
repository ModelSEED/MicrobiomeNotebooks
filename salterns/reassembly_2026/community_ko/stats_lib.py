"""Shared statistics for the community gene profiles (Figure 1, Figure S4, Figure S6, Table S1).

site_type(sample)     -> 'unrestored' (R1, R2 cores), 'restored' (SF2), 'reference' (R2A), or None
kw_dunn(values, groups) -> Kruskal-Wallis H and p, Dunn's pairwise z tests (tie-corrected, Holm-adjusted within the
                         family), and compact letter display (groups sharing a letter do not differ at alpha)
bh(p)                 -> Benjamini-Hochberg q values
"""
import itertools
import math

import numpy as np
import pandas as pd
from scipy import stats

SITES = ["R1", "R2", "SF2", "R2A"]
TYPES = ["unrestored", "restored", "reference"]


def site_of(sample):
    s = sample.split("_")[0]
    return s if s in SITES else None


def site_type(sample):
    s = site_of(sample)
    return {"R1": "unrestored", "R2": "unrestored", "SF2": "restored", "R2A": "reference"}.get(s)


def bh(p):
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    n = ok.sum()
    if n == 0:
        return q
    order = np.argsort(p[ok])
    ranked = p[ok][order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1)
    q[ok] = out
    return q


def dunn(values, groups):
    """Dunn's test, tie-corrected; returns {(a, b): p} two-sided, unadjusted."""
    v = pd.Series(np.asarray(values, float))
    g = pd.Series(list(groups))
    r = v.rank()
    N = len(v)
    ties = v.value_counts()
    tie = (ties ** 3 - ties).sum() / (12 * (N - 1)) if N > 1 else 0
    out = {}
    levels = list(dict.fromkeys(g))
    for a, b in itertools.combinations(levels, 2):
        ra, rb = r[g == a].mean(), r[g == b].mean()
        na, nb = (g == a).sum(), (g == b).sum()
        se = math.sqrt(max((N * (N + 1) / 12 - tie) * (1 / na + 1 / nb), 1e-300))
        z = (ra - rb) / se
        out[(a, b)] = 2 * stats.norm.sf(abs(z))
    return out


def holm(pdict):
    items = sorted(pdict.items(), key=lambda kv: kv[1])
    m = len(items)
    adj, run = {}, 0.0
    for i, (k, p) in enumerate(items):
        run = max(run, min(1.0, (m - i) * p))
        adj[k] = run
    return adj


def cld(levels, padj, order_by, alpha=0.05):
    """Compact letter display: 'a' for the highest group in order_by (a dict level -> mean rank or median)."""
    levels = sorted(levels, key=lambda x: -order_by[x])
    diff = {frozenset(k) for k, p in padj.items() if p < alpha}
    letters = []                                   # list of sets of levels
    for lv in levels:
        placed = False
        for s in letters:
            if all(frozenset((lv, o)) not in diff for o in s):
                s.add(lv); placed = True
        if not placed:
            letters.append({lv})
    # drop letter sets that are subsets of others
    letters = [s for s in letters if not any(s < t for t in letters)]
    abc = "abcdefgh"
    return {lv: "".join(abc[i] for i, s in enumerate(letters) if lv in s) for lv in levels}


def kw_dunn(values, groups, alpha=0.05):
    v = np.asarray(values, float)
    g = list(groups)
    levels = list(dict.fromkeys(g))
    samples = [v[[x == lv for x in g]] for lv in levels]
    if len({float(x) for x in v}) < 2:
        return dict(H=0.0, p=1.0, padj={}, letters={lv: "a" for lv in levels}, medians={lv: float(np.median(s)) for lv, s in zip(levels, samples)})
    H, p = stats.kruskal(*samples)
    padj = holm(dunn(v, g))
    ranks = pd.Series(v).rank()
    mean_rank = {lv: ranks[[x == lv for x in g]].mean() for lv in levels}
    return dict(H=float(H), p=float(p), padj=padj, letters=cld(levels, padj, mean_rank, alpha),
                medians={lv: float(np.median(s)) for lv, s in zip(levels, samples)})
