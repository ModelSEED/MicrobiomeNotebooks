#!/usr/bin/env python
"""
Re-assessment figure: why the cross-sample correlations in Figure 7 / S8 do not,
on their own, establish a MePn- (or betaine-) specific mechanism.

(a) Methane is dominated by ONE environmental axis -- sulfate/salinity/N:P -- more
    strongly than by any biological variable (n = 12 cores).
(b) SPECIFICITY null: the summed abundance of a *random* MAG subset of the same size
    correlates with methane about as well as the MePn / betaine consumer sets, so the
    consumer->methane correlation (Fig 7b / S8) is not consumer-specific.
(c) Measured-metabolite correlations (n = 10 cores): TMA is the strongest, and the
    amine/betaine-route substrates (TMA, MMA, betaine) all correlate with methane.

Run: ~/Documents/py_venv/bin/python reassessment_figure.py
"""
from __future__ import annotations
import json
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
np.seterr(all="ignore")

DATA = "../data"
C1, C2, C3 = "#D55E00", "#0072B2", "#009E73"   # Okabe-Ito

DNA = pd.read_csv(f"{DATA}/SaltPondsDNA.csv").set_index("Sample")["DNA conc (µg/g soil)"]
PRE = sorted({"_".join(i.split("_")[:-1]) for i in DNA.index})
rel = pd.read_csv(f"{DATA}/Cliff_310MAG_relabund.tsv", sep="\t", index_col=0)
ab = pd.DataFrame({p: (rel[p+"_D1"]*DNA[p+"_D1"] + rel[p+"_D2"]*DNA[p+"_D2"]) /
                     (DNA[p+"_D1"]+DNA[p+"_D2"]) for p in PRE})
meta = pd.read_csv(f"{DATA}/Cliff_Sample_Metadata_BGC_NMR.csv").set_index("Sample")
meta["core"] = ["_".join(i.split("_")[:-1]) for i in meta.index]
num = meta.apply(pd.to_numeric, errors="coerce")
core = num.groupby(meta["core"]).mean(numeric_only=True)
order = list(ab.columns)
ch4 = pd.Series({p: meta.loc[[p+"_D1", p+"_D2"], "CH4_umol_m2_d"].astype(float).mean() for p in PRE})
y = ch4.reindex(order).values
yc = core["CH4_umol_m2_d"]

mepn = [k.replace(".contigs", "") for k in json.load(open("../correlational/saltern_fluxes.json"))["0.07"]["methane"]]
bet = [k.replace(".contigs", "") for k in json.load(open("../correlational/saltern_fluxes_betaine.json"))["TMA"]]
mepn = [i for i in mepn if i in ab.index]; bet = [i for i in bet if i in ab.index]

fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16.5, 5.2), dpi=200)

# (a) methane vs the dominant environmental variables
env = [("SO4_S", "Sulfate (SO4-S)"), ("Salinity", "Salinity"), ("Olsen_P", "Olsen-P (inorg. phosphate)")]
cols = [C1, C2, C3]
for (g, lab), c in zip(env, cols):
    v = core[g]; ok = v.notna() & yc.notna()
    r, p = stats.spearmanr(v[ok], yc[ok])
    xn = (v[ok]-v[ok].min())/(v[ok].max()-v[ok].min())     # min-max for shared x
    axA.scatter(xn, yc[ok], color=c, s=55, edgecolor="white", lw=0.6,
                label=f"{lab}: rho={r:+.2f} (p={p:.3f})")
axA.set_yscale("log")
axA.set_xlabel("environmental variable (min-max scaled)", fontsize=11)
axA.set_ylabel(r"CH$_4$ emission [$\mu$mol m$^{-2}$ d$^{-1}$]", fontsize=11)
axA.set_title("a  Methane is set by the salinity/sulfate gradient\n(n = 12 cores)",
              fontsize=12, fontweight="bold", loc="left")
axA.legend(fontsize=8.5, loc="lower right"); axA.grid(alpha=0.3, which="both", ls="--")

# (b) specificity null
rng = np.random.default_rng(0)
def rho_set(ids):
    r, _ = stats.spearmanr(ab.loc[ids].sum().reindex(order).values, y); return r
null20 = np.array([rho_set(list(ab.index[rng.choice(len(ab), 20, replace=False)])) for _ in range(4000)])
null15 = np.array([rho_set(list(ab.index[rng.choice(len(ab), 15, replace=False)])) for _ in range(4000)])
obs_m, obs_b = rho_set(mepn), rho_set(bet)
axB.hist(null20[np.isfinite(null20)], bins=40, color="0.75", alpha=0.8,
         label="random 20-MAG subsets")
axB.axvline(obs_m, color=C1, lw=2.5, label=f"MePn consumers (rho={obs_m:.2f})")
axB.axvline(obs_b, color=C2, lw=2.5, ls="--", label=f"Betaine consumers (rho={obs_b:.2f})")
fm = (null20 >= obs_m).mean() * 100
axB.set_xlabel(r"Spearman rho (subset abundance vs. CH$_4$)", fontsize=11)
axB.set_ylabel("random subsets", fontsize=11)
axB.set_title(f"b  Consumer signal is NOT specific\n{fm:.0f}% of random subsets do as well as MePn",
              fontsize=12, fontweight="bold", loc="left")
axB.legend(fontsize=8.5, loc="upper left"); axB.grid(alpha=0.3, ls="--")

# (c) metabolite correlations
mets = ["Trimethylamine", "Trehalose", "Methylamine", "Acetate", "Betaine", "Methanol", "Dimethylamine"]
rs = []
for m in mets:
    v = core.get(m); ok = v.notna() & yc.notna()
    r, _ = stats.spearmanr(v[ok], yc[ok]); rs.append(r)
colc = [C1 if x > 0 else C2 for x in rs]
axC.barh(range(len(mets)), rs, color=colc, edgecolor="white")
axC.set_yticks(range(len(mets))); axC.set_yticklabels(mets, fontsize=10)
axC.invert_yaxis(); axC.axvline(0, color="0.4", lw=1)
axC.set_xlabel(r"Spearman rho vs. CH$_4$ (n = 10)", fontsize=11)
axC.set_title("c  Measured-metabolite evidence\n(TMA strongest; DMA negative)",
              fontsize=12, fontweight="bold", loc="left")
axC.grid(axis="x", alpha=0.3, ls="--")
for i, r in enumerate(rs):
    axC.text(r + (0.02 if r > 0 else -0.02), i, f"{r:+.2f}", va="center",
             ha="left" if r > 0 else "right", fontsize=9)

plt.tight_layout()
fig.savefig("reassessment_gradient_confound.png", bbox_inches="tight")
fig.savefig("reassessment_gradient_confound.svg", bbox_inches="tight")
print(f"MePn random-subset exceedance: {(null20>=obs_m).mean():.3f}   "
      f"Betaine: {(null15>=obs_b).mean():.3f}")
print("Saved reassessment_gradient_confound.png/.svg")
