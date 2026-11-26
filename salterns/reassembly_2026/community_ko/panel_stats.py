#!/usr/bin/env python
"""Per-KO statistics of the methylated-compound panel across the three restoration types (Figure 1, Table S1).

For each of the 138 panel KOs detected (>= 1 read, or >= 1 assembled gene) in at least one of the 24 sediment
metagenomes: Kruskal-Wallis across unrestored (R1 + R2, n = 12), restored (SF2, n = 6) and reference (R2A, n = 6), on
copies per genome equivalent; Benjamini-Hochberg FDR over all detected panel KOs (one family, used by Figure 1 and
Table S1); Dunn's test with Holm adjustment and compact letters ('a' = highest mean rank).
source = read (read survey named in ../read_survey_current.json) or assembly (KOfam on the assembled genes).
Writes results/panel_stats_<source>.tsv. Run: panel_stats.py [--source read|assembly] [--results DIR]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import profiles as PR  # noqa: E402
import stats_lib as SL  # noqa: E402

TYPES = ["unrestored", "restored", "reference"]
PAIRS = (("unrestored", "restored"), ("unrestored", "reference"), ("restored", "reference"))


def compute(source="read", results=None):
    P = PR.panel().set_index("ko")
    t = PR.profile(source, results)
    W = PR.wide(t, "copies_per_GE", PR.SED, kos=list(P.index))
    R = PR.wide(t, "reads", PR.SED, kos=list(P.index))
    types = [SL.site_type(s) for s in PR.SED]
    rows = []
    for k in P.index:
        surveyed = bool(W.loc[k].notna().all())
        n = int(R.loc[k].fillna(0).sum()) if surveyed else 0
        row = dict(ko=k, gene=P.loc[k, "gene"], surveyed=surveyed, n_24=n if surveyed else np.nan,
                   n_samples_detected=int((R.loc[k].fillna(0) > 0).sum()) if surveyed else np.nan)
        v = W.loc[k].astype(float).values if surveyed else None
        for ty in TYPES:
            row[f"mean_per100GE_{ty}"] = float(np.mean(v[[x == ty for x in types]]) * 100) if surveyed else np.nan
        if surveyed and n > 0:
            r = SL.kw_dunn(v, types)
            row.update(H=r["H"], p=r["p"], **{f"letter_{ty}": r["letters"][ty] for ty in TYPES},
                       **{f"dunn_p_{a}_vs_{b}": r["padj"].get((a, b), r["padj"].get((b, a), np.nan)) for a, b in PAIRS})
            mr = pd.Series(v).rank().groupby(pd.Series(types)).mean()
            row["highest"], row["lowest"] = mr.idxmax(), mr.idxmin()
        rows.append(row)
    st = pd.DataFrame(rows).set_index("ko")
    det = st.p.notna()
    st["q"] = np.nan
    st.loc[det, "q"] = SL.bh(st.loc[det, "p"].values)
    st["source"] = source
    return st


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", choices=["read", "assembly"], default="read")
    ap.add_argument("--results", default=None)
    a = ap.parse_args()
    st = compute(a.source, a.results)
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    st.to_csv(os.path.join(HERE, "results", f"panel_stats_{a.source}.tsv"), sep="\t")
    d = st[st.p.notna()]
    print(f"{a.source}: {int(st.surveyed.sum())} surveyed, {len(d)} detected in the 24 sediment samples, "
          f"{int((d.q < 0.05).sum())} with KW FDR < 0.05")


if __name__ == "__main__":
    main()
