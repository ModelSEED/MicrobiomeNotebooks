#!/usr/bin/env python
"""Ceiling check before any community-modelling effort: does any measured input rank the methane flux of the six
unrestored cores? A constraint-based model only turns its inputs (abundance, metabolite concentrations, geochemistry)
into fluxes, so if no input, or simple guild x substrate product, tracks CH4 across the cores, no formulation can.
Candidates are fixed here, before looking at results: guild biomass (areal: r x DNA x layer, as in the models),
guild relative abundance, NMR metabolites that are methanogenic substrates or relevant donors, and P / S / salinity
geochemistry. Spearman with exact permutation p (n = 6: |rho| >= 0.886 for p < 0.05). With ~25 candidates about one
nominal hit is expected by chance, so a hit is a lead, not a result.
Run: ~/Documents/py_venv/bin/python input_signal_check.py -> input_signal_check.tsv
"""
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
S = HERE.parents[2]
CORES = ["R1_A", "R1_B", "R1_C", "R2_A", "R2_B", "R2_C"]
LAYER = {"D1": 5, "D2": 10}
meta = pd.read_csv(S / "data/Cliff_Sample_Metadata_BGC_NMR.csv").set_index("Sample")
dna = pd.read_csv(S / "data/SaltPondsDNA.csv").set_index("Sample")["DNA conc (µg/g soil)"]
rel = pd.read_csv(S / "reassembly_2026/abundance_final/final585_relabund_24sediment.tsv", sep="\t", index_col=0)
g = pd.read_csv(S / "correlational/final585/guilds_final585.csv").set_index("mag")
ch4 = pd.Series({c: meta.loc[[c + "_D1", c + "_D2"], "CH4_umol_m2_d"].mean() for c in CORES})


def core_mean(col):
    return pd.Series({c: meta.loc[[c + "_D1", c + "_D2"], col].astype(float).mean() for c in CORES})


def biomass(mags):
    return pd.Series({c: sum((rel.loc[mags, f"{c}_{d}"] * dna[f"{c}_{d}"] * cm).sum() for d, cm in LAYER.items()) for c in CORES})


def relab(mags):
    return pd.Series({c: rel.loc[mags, [f"{c}_D1", f"{c}_D2"]].sum().mean() for c in CORES})


def exact(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = ~np.isnan(x)
    x, y = x[ok], y[ok]
    r0 = stats.spearmanr(x, y)[0]
    perms = list(itertools.permutations(range(len(x))))
    p = sum(abs(stats.spearmanr(x, y[list(q)])[0]) >= abs(r0) - 1e-12 for q in perms) / len(perms)
    return r0, p, int(ok.sum())


cand = {}
for gname in ("mepn_consumer", "betaine_reducer", "methanogen", "archaea"):
    mags = list(g.index[g[gname]])
    cand[f"{gname} biomass (areal)"] = biomass(mags)
    cand[f"{gname} rel. abundance"] = relab(mags)
cand["all binned species biomass (areal)"] = biomass([m for m in rel.index if m != "unbinned_reads"])
cand["soil DNA (mean of depths)"] = pd.Series({c: dna[[c + "_D1", c + "_D2"]].mean() for c in CORES})
for col in ("Betaine", "Trimethylamine", "Dimethylamine", "Methylamine", "Methanol", "Acetate", "Formate", "Lactate",
            "N,N-Dimethylglycine", "Trehalose", "Olsen_P", "P", "SO4_S", "Salinity", "TOC", "NH4", "NO3"):
    cand[col] = core_mean(col)
cand["MePn biomass x (N/P) [P-limitation proxy]"] = cand["mepn_consumer biomass (areal)"] * (core_mean("N") / core_mean("P"))
cand["betaine-reducer biomass x betaine"] = cand["betaine_reducer biomass (areal)"] * cand["Betaine"]
cand["methanogen biomass x (TMA + methylamines + methanol)"] = cand["methanogen biomass (areal)"] * (
    cand["Trimethylamine"].fillna(0) + cand["Dimethylamine"].fillna(0) + cand["Methylamine"].fillna(0) + cand["Methanol"].fillna(0))
rows = []
for k, v in cand.items():
    r, p, n = exact(v[CORES], ch4[CORES])
    rows.append({"input": k, "n_cores": n, "spearman_rho": round(r, 3), "exact_p": round(p, 3),
                 **{c: v[c] for c in CORES}})
t = pd.DataFrame(rows).sort_values("exact_p")
t.to_csv(HERE / "input_signal_check.tsv", sep="\t", index=False)
print("measured CH4:", ch4.round(0).to_dict())
print(t[["input", "n_cores", "spearman_rho", "exact_p"]].to_string(index=False))
