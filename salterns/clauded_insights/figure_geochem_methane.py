#!/usr/bin/env python
"""
Proposed new figure: methane is governed by the salinity/sulfate/P-limitation gradient,
and MAG-consumer abundance explains almost nothing beyond it.

(a) CH4 vs sulfate (the single dominant correlate).
(b) CH4 vs N:P ratio (supports the phosphate-limitation -> methylphosphonate story).
(c) Variance partitioning (n = 12 cores): of the methane variance, how much is uniquely
    the geochemical gradient, shared, vs. uniquely MePn-consumer abundance.

Run: ~/Documents/py_venv/bin/python figure_geochem_methane.py
"""
from __future__ import annotations
import json
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = "../data"
C_ENV, C_AB, C_SHARE = "#0072B2", "#D55E00", "#999999"

DNA = pd.read_csv(f"{DATA}/SaltPondsDNA.csv").set_index("Sample")["DNA conc (µg/g soil)"]
PRE = sorted({"_".join(i.split("_")[:-1]) for i in DNA.index})
rel = pd.read_csv(f"{DATA}/Cliff_310MAG_relabund.tsv", sep="\t", index_col=0)
ab = pd.DataFrame({p: (rel[p+"_D1"]*DNA[p+"_D1"] + rel[p+"_D2"]*DNA[p+"_D2"]) /
                     (DNA[p+"_D1"]+DNA[p+"_D2"]) for p in PRE})
meta = pd.read_csv(f"{DATA}/Cliff_Sample_Metadata_BGC_NMR.csv").set_index("Sample")
meta["core"] = ["_".join(i.split("_")[:-1]) for i in meta.index]
core = meta.apply(pd.to_numeric, errors="coerce").groupby(meta["core"]).mean(numeric_only=True)
core = core.reindex(list(ab.columns))
mepn = [k.replace(".contigs", "") for k in json.load(open("../correlational/saltern_fluxes.json"))["0.07"]["methane"]]
mepn = [i for i in mepn if i in ab.index]
core["mepn_ab"] = ab.loc[mepn].sum().reindex(core.index).values
core["NP"] = core["N"] / core["P"]
y = core["CH4_umol_m2_d"]


def r2(cols):
    # RANK-based variance partition (consistent with the Spearman correlations):
    # rank-transform response + predictors so the strongly-but-nonlinearly related
    # gradient (sulfate) is represented by its rank association, not a poor linear fit.
    X = core[cols].values
    m = np.isfinite(X).all(1) & np.isfinite(y.values)
    X, yy = X[m], y.values[m]
    Xr = np.column_stack([stats.rankdata(X[:, j]) for j in range(X.shape[1])])
    yr = stats.rankdata(yy)
    Xr = np.c_[np.ones(len(Xr)), (Xr - Xr.mean(0)) / Xr.std(0)]
    beta = np.linalg.lstsq(Xr, yr, rcond=None)[0]
    pred = Xr @ beta
    return 1 - ((yr - pred) ** 2).sum() / ((yr - yr.mean()) ** 2).sum()


R_env, R_ab, R_both = r2(["SO4_S"]), r2(["mepn_ab"]), r2(["SO4_S", "mepn_ab"])
uniq_env = R_both - R_ab
uniq_ab = R_both - R_env
shared = R_env + R_ab - R_both

fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16.5, 5.0), dpi=200)

for ax, col, lab in [(axA, "SO4_S", "Sulfate (SO$_4$-S)"), (axB, "NP", "N : P ratio")]:
    v = core[col]; ok = v.notna() & y.notna()
    r, p = stats.spearmanr(v[ok], y[ok])
    ax.scatter(v[ok], y[ok], s=60, color=C_ENV, edgecolor="white", lw=0.7, zorder=3)
    ax.set_yscale("log")
    b = np.polyfit(v[ok], np.log10(y[ok]), 1)
    xs = np.linspace(v[ok].min(), v[ok].max(), 100)
    ax.plot(xs, 10 ** (b[0] * xs + b[1]), "--", color="#0a3d6b", lw=1.6)
    ax.set_xlabel(lab, fontsize=12)
    ax.set_ylabel(r"CH$_4$ emission [$\mu$mol m$^{-2}$ d$^{-1}$]", fontsize=11)
    ax.set_title(("a" if col == "SO4_S" else "b") + f"  Spearman rho = {r:+.2f}, p = {p:.3f} (n = {ok.sum()})",
                 fontsize=12.5, fontweight="bold", loc="left")
    ax.grid(alpha=0.3, which="both", ls="--")

# variance partition
parts = [max(uniq_env, 0), max(shared, 0), max(uniq_ab, 0)]
labels = [f"gradient only\n(SO4-S)\n{uniq_env:.2f}", f"shared\n{shared:.2f}",
          f"MePn abundance\nUNIQUE\n{uniq_ab:.2f}"]
colors = [C_ENV, C_SHARE, C_AB]
bottom = 0
for pval, lab, c in zip(parts, labels, colors):
    axC.bar(0, pval, bottom=bottom, width=0.6, color=c, edgecolor="white")
    if pval > 0.01:
        axC.text(0, bottom + pval / 2, lab, ha="center", va="center", fontsize=10,
                 color="white" if c != C_SHARE else "black", fontweight="bold")
    bottom += pval
axC.set_xlim(-0.7, 0.7); axC.set_xticks([])
axC.set_ylabel(r"share of CH$_4$ variance (log) explained  [R$^2$]", fontsize=11)
axC.set_title(f"c  Variance partition (n = 12)\nMePn abundance adds only {uniq_ab:.2f} beyond the gradient",
              fontsize=12.5, fontweight="bold", loc="left")
axC.set_ylim(0, max(R_both, 0.6) * 1.1)
axC.grid(axis="y", alpha=0.3, ls="--")

plt.tight_layout()
fig.savefig("figure_geochem_methane.png", bbox_inches="tight")
fig.savefig("figure_geochem_methane.svg", bbox_inches="tight")
print(f"R2(gradient)={R_env:.3f}  R2(MePn abund)={R_ab:.3f}  R2(both)={R_both:.3f}")
print(f"unique gradient={uniq_env:.3f}  shared={shared:.3f}  unique MePn abundance={uniq_ab:.3f}")
print("Saved figure_geochem_methane.png/.svg")
