#!/usr/bin/env python
"""Metabolome vs community composition (Results P63, Discussion P82): Mantel test and Bray-Curtis ranges.

Metabolome: all 50 NMR metabolites of data/Cliff_Sample_Metadata_BGC_NMR.csv in the 19 sediment samples with NMR data
(the whole profile, as for the original Mantel test; metabolite() with min_samples=3 gives the 40 used for correlations). Community: (a) the 585 species of the final MAG
set (abundance_final/final585_relabund.tsv; % of all reads, unbinned reads excluded) and, when present, (b) the read-based
genus profile of the ribosomal-protein marker reads (read_markers_v2/results/marker_genus_profile.tsv).
Both sides: Bray-Curtis on square-root-transformed data (as in the original analysis). Mantel: Pearson r on the
condensed distance vectors, 9,999 permutations of the samples (one-sided, r > 0). Spearman Mantel is also given.
Writes results/metabolome_taxonomy.json.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
R26 = os.path.dirname(HERE)
SAL = os.path.dirname(R26)


def metabolome(min_samples=0):
    m = pd.read_csv(os.path.join(SAL, "data", "Cliff_Sample_Metadata_BGC_NMR.csv"))
    cols = list(m.columns[m.columns.get_loc("2'-Deoxyuridine"):m.columns.get_loc("Valine") + 1])
    m["sample"] = m.Sample_name.str.replace(r"^(Historic|restored|unrestored)_", "", regex=True)
    X = m.set_index("sample")[cols].apply(pd.to_numeric, errors="coerce")
    X = X[X.notna().any(axis=1)].fillna(0.0)
    if min_samples:
        X = X.loc[:, (X > 0).sum() >= min_samples]
    return X


def mantel(D1, D2, n=9999, seed=0):
    a, b = squareform(D1, checks=False), squareform(D2, checks=False)
    r = pearsonr(a, b)[0]
    rs = spearmanr(a, b)[0]
    rng = np.random.default_rng(seed)
    k = D1.shape[0]
    cnt = 0
    for _ in range(n):
        p = rng.permutation(k)
        if pearsonr(a, squareform(D2[np.ix_(p, p)], checks=False))[0] >= r:
            cnt += 1
    return float(r), float(rs), (cnt + 1) / (n + 1)


def bc(X):
    return squareform(pdist(np.sqrt(X.values), metric="braycurtis"))


def main():
    X = metabolome()
    out = {"n_samples": int(X.shape[0]), "n_metabolites": int(X.shape[1]), "samples": list(X.index)}
    Dm = bc(X)
    ra = pd.read_csv(os.path.join(R26, "abundance_final", "final585_relabund.tsv"), sep="\t", index_col=0)
    sp = ra.drop(index="unbinned_reads")[X.index].T * 100
    comms = {"MAG species (585)": sp}
    rp = os.path.join(R26, "read_markers_v2", "results", "marker_genus_profile.tsv")
    if os.path.exists(rp):
        t = pd.read_csv(rp, sep="\t")
        gcol = next(c for c in ("genus_harmonized", "genus", "genus_gtdb") if c in t.columns)
        scol = next(c for c in ("share_of_marker_reads", "marker_share", "share") if c in t.columns)
        comms["read marker genera"] = t.pivot_table(index="sample", columns=gcol, values=scol, aggfunc="sum",
                                                    fill_value=0).reindex(X.index).fillna(0)
    for name, C in comms.items():
        Dt = bc(C)
        r, rs, p = mantel(Dm, Dt)
        iu = np.triu_indices(len(X), 1)
        dm, dt = Dm[iu], Dt[iu]
        top = dt >= np.quantile(dt, 0.75)
        out[name] = dict(mantel_r=r, mantel_spearman=rs, p=p, permutations=9999,
                         metab_bc=[float(dm.min()), float(dm.max())], taxo_bc=[float(dt.min()), float(dt.max())],
                         top_quartile_taxo_bc=[float(dt[top].min()), float(dt[top].max())],
                         metab_bc_in_top_quartile=[float(dm[top].min()), float(dm[top].max())],
                         n_pairs_taxo_ge_083=int((dt >= 0.83).sum()),
                         taxo_bc_ge_083=[float(dt[dt >= 0.83].min()), float(dt[dt >= 0.83].max())] if (dt >= 0.83).any() else None,
                         metab_bc_where_taxo_ge_083=[float(dm[dt >= 0.83].min()), float(dm[dt >= 0.83].max())]
                         if (dt >= 0.83).any() else None)
        print(name, json.dumps(out[name]))
    json.dump(out, open(os.path.join(HERE, "results", "metabolome_taxonomy.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
