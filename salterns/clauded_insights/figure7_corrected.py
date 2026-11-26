#!/usr/bin/env python
"""
Corrected Figure 7  (metabolic-modeling regressions of saltern methane).

Reproduces the manuscript's Figure 7 from the *already-computed* flux/abundance
data in the repo, and fixes the defects identified in review while KEEPING the
two features the author wants: a **linear fitted equation** (intuitive units) and
**log-log axes** (so all points fit compressed in one frame).

What changed vs. the submitted Figure 7 (docx image11)
------------------------------------------------------
1. Trendline rendered correctly.  The old code evaluated y = m*x + b only at the
   6 data points and connected them in raw-x order on log axes, producing a broken
   vertical-then-flat curve.  Here the *same* linear fit is evaluated on a dense
   x-grid spanning the axis and drawn as a smooth curve (clipped where y<=0, which
   log axes cannot show -- itself a flag that a negative-intercept linear fit is
   unphysical for low-abundance samples).
2. Honest statistics printed on every series: n, Pearson R^2, and Spearman rho with
   its p-value.  At n=6 a Spearman rho of 0.6 has p ~ 0.21 (not significant); the
   gap between the high Pearson R^2 and the modest rho is the tell that 1-2
   high-leverage points carry the Pearson fit.
3. Panel (a) placed on a *consistent* unit basis (raw areal emission on both axes,
   umol CH4 m^-2 d^-1) so a genuine 1:1 model-validation line can be drawn.  The
   submitted figure divided the measured axis by DNA density while the modelled
   axis was multiplied by biomass (~DNA), so the two axes were not the same
   quantity and no 1:1 line was possible.  The DNA-normalised numbers are still
   computed and printed for reconciliation with the text.
4. Panel (b) extended to all 12 site-cores (unrestored R1/R2, restored SF2,
   reference R2A) coloured by restoration class -- tripling n and showing the
   restoration gradient that motivates the study -- instead of the 6 unrestored
   points only.

All statistics are printed to stdout and dumped to figure7_stats.json for the
numbers_to_reconcile.md table.  This script is deterministic and self-contained;
run:  ~/Documents/py_venv/bin/python figure7_corrected.py
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatterMathtext, NullFormatter

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"
CORR = REPO / "correlational"
OUT = Path(__file__).resolve().parent

# Okabe-Ito colourblind-safe categorical pair (kept consistent across all figures)
C_MEPN = "#D55E00"   # vermillion  -> methylphosphonate
C_BET  = "#0072B2"   # blue        -> betaine
REST_COLORS = {"unrestored": "#D55E00", "restored": "#009E73", "reference": "#56B4E9"}

# ---------------------------------------------------------------- unit conversion
DNA = pd.read_csv(DATA / "SaltPondsDNA.csv").set_index("Sample")["DNA conc (µg/g soil)"]
PREFIXES = sorted({"_".join(i.split("_")[:-1]) for i in DNA.index})
# volume-weighted DNA density per core (5 cm D1 over 10 cm D2), exactly as the pipeline
VW_DNA = pd.Series({p: (5 * DNA[p + "_D1"] + 10 * DNA[p + "_D2"]) / 15 for p in PREFIXES})


def fluxes_to_emissions(row: pd.Series) -> pd.Series:
    """mmol/hr/g_biomass community flux -> umol CH4 m^-2 d^-1 (verbatim from pipeline)."""
    area_m2 = math.pi * (5 / 2) ** 2 * 0.01 ** 2
    volume_cm3 = 294.52
    ug_dna = VW_DNA.reindex(row.index).fillna(2.0)
    dna_per_biomass = 0.1
    g_soil_cm3 = 0.7
    g_biomass_g_soil = g_soil_cm3 * volume_cm3 * ug_dna / dna_per_biomass / 1e6
    return row / area_m2 * 1000.0 * 24 * g_biomass_g_soil


def merge_depths_by_DNA(frame: pd.DataFrame) -> pd.DataFrame:
    """DNA-weighted average of a per-depth table across the two depths (D1 5cm, D2 10cm)
    -> one column per core.  This is a weighted MEAN (weights in denominator); it does
    NOT divide by DNA, so the result is still on the raw relative-abundance scale."""
    return pd.DataFrame(
        {p: (frame[p + "_D1"] * DNA[p + "_D1"] + frame[p + "_D2"] * DNA[p + "_D2"])
             / (DNA[p + "_D1"] + DNA[p + "_D2"]) for p in PREFIXES},
        index=frame.index)


BETAINE_CH4_PER_TMA = 0.06   # manuscript: methane flux assumed = 0.06 x TMA excretion (Fig 4)


# ---------------------------------------------------------------- measured methane
meta = pd.read_csv(DATA / "Cliff_Sample_Metadata_BGC_NMR.csv").set_index("Sample")
ch4 = meta["CH4_umol_m2_d"]
# raw areal emission per core = mean over the two depth rows
CH4_RAW = pd.Series({p: ch4[[p + "_D1", p + "_D2"]].mean() for p in PREFIXES})
# DNA-normalised (as the submitted figure's x-axis)
ch4_norm = ch4 / DNA
CH4_NORM = pd.Series({p: ch4_norm[[p + "_D1", p + "_D2"]].mean() for p in PREFIXES})

UNRESTORED = [p for p in PREFIXES if p.startswith(("R1_", "R2_")) and "restored" not in p]


def rest_class(prefix: str) -> str:
    if prefix.startswith("R2A"):
        return "reference"
    if prefix.startswith("SF2"):
        return "restored"
    return "unrestored"


# ---------------------------------------------------------------- modelled methane
def summed_flux_to_emissions(flux_by_mag: dict) -> pd.Series:
    summed = pd.DataFrame(flux_by_mag).T.sum().abs()      # sum over member MAGs, per sample
    return fluxes_to_emissions(summed).reindex(UNRESTORED)

_mepn_j = json.load(open(CORR / "saltern_fluxes.json"))["0.07"]
mepn_model = summed_flux_to_emissions(_mepn_j["methane"])           # direct C-P-lyase CH4
_bet_j = json.load(open(CORR / "saltern_fluxes_betaine.json"))
# betaine 'methane' field is empty by construction; CH4 = 0.06 x TMA excretion (Methods)
_bet_ch4 = {mag: {s: BETAINE_CH4_PER_TMA * v for s, v in d.items()} for mag, d in _bet_j["TMA"].items()}
bet_model = summed_flux_to_emissions(_bet_ch4)

# ---------------------------------------------------------------- consumer abundance
# RAW relative abundance (depth-merged, NOT divided by DNA) -> avoids the shared-1/DNA
# spurious-correlation artefact of the submitted figure.  Also keep the DNA-normalised
# version to reproduce/expose the inflated numbers reported in the text.
_raw_rel = pd.read_csv(DATA / "Cliff_310MAG_relabund.tsv", sep="\t", index_col=0)
ab_raw = merge_depths_by_DNA(_raw_rel)                              # 311 x 12, raw scale
ab_dnanorm = pd.read_csv(DATA / "averaged_normalized_MAG_abundances.csv", index_col=0)  # /DNA
mepn_ids = [k.replace(".contigs", "") for k in _mepn_j["methane"] if k.replace(".contigs", "") in ab_raw.index]
bet_ids = [k.replace(".contigs", "") for k in _bet_j["TMA"] if k.replace(".contigs", "") in ab_raw.index]
mepn_abund = ab_raw.loc[mepn_ids].sum()          # per sample (all 12), raw scale
bet_abund = ab_raw.loc[bet_ids].sum()
mepn_abund_dn = ab_dnanorm.loc[mepn_ids].sum()   # DNA-normalised (as submitted)
bet_abund_dn = ab_dnanorm.loc[bet_ids].sum()

# ---------------------------------------------------------------- stats helper
def fit_stats(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    lr = stats.linregress(x, y)                      # linear fit in raw space (intuitive eq)
    rho, p = stats.spearmanr(x, y)
    return dict(n=int(len(x)), slope=lr.slope, intercept=lr.intercept,
                pearson_r2=lr.rvalue ** 2, spearman_rho=rho, spearman_p=p)


def log_axes(ax):
    ax.set_xscale("log"); ax.set_yscale("log")
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(LogLocator(base=10, subs=(1.0,)))
        axis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
        axis.set_major_formatter(LogFormatterMathtext(base=10, labelOnlyBase=True))
        axis.set_minor_formatter(NullFormatter())
    ax.grid(True, which="both", linestyle="--", color="gray", alpha=0.3)


def draw_linear_on_log(ax, xdata, s, color):
    """Draw the linear fit y=slope*x+intercept as a DENSE curve across the x-range
    (smooth on log axes), clipped to y>0 which log can display."""
    xmin, xmax = np.nanmin(xdata), np.nanmax(xdata)
    xs = np.linspace(xmin, xmax, 400)
    ys = s["slope"] * xs + s["intercept"]
    m = ys > 0
    ax.plot(xs[m], ys[m], "--", color=color, lw=1.6, zorder=2)


def annotate(ax, s, color, xy):
    eq = f"y = {s['slope']:.2g}x + {s['intercept']:.2g}"
    txt = (f"{eq}\n$R^2$(Pearson) = {s['pearson_r2']:.2f}\n"
           f"$r_s$ = {s['spearman_rho']:.2f}, p = {s['spearman_p']:.2f}  (n={s['n']})")
    ax.text(*xy, txt, transform=ax.transAxes, fontsize=9, color=color, va="top",
            bbox=dict(boxstyle="round", fc="white", ec=color, alpha=0.9))


def label_points(ax, xs, ys, labels, color, collect):
    for x, y, lab in zip(xs, ys, labels):
        if np.isfinite(x) and np.isfinite(y) and y > 0:
            collect.append(ax.text(x, y, lab, fontsize=7.5, color=color))


# ================================================================ FIGURE
from matplotlib.lines import Line2D
from adjustText import adjust_text
fig, (axa, axb) = plt.subplots(2, 1, figsize=(8.4, 13), dpi=200)

# ---- Panel (a): modelled vs measured areal emission (consistent raw units) ------
xa_meas = CH4_RAW.reindex(UNRESTORED).values
sa_mepn = fit_stats(xa_meas, mepn_model.values)
sa_bet = fit_stats(xa_meas, bet_model.values)
_texts = []
for model, s, c, lab in [(mepn_model, sa_mepn, C_MEPN, "MePn"), (bet_model, sa_bet, C_BET, "Betaine")]:
    axa.scatter(xa_meas, model.values, s=55, color=c, edgecolor="white", lw=0.6, zorder=3, label=lab)
    draw_linear_on_log(axa, xa_meas, s, c)
    label_points(axa, xa_meas, model.values, UNRESTORED, c, _texts)
log_axes(axa)
# 1:1 reference (meaningful now that both axes are raw areal emission)
lim = [min(np.nanmin(xa_meas), mepn_model.min(), bet_model.min()) * 0.6,
       max(np.nanmax(xa_meas), mepn_model.max(), bet_model.max()) * 1.6]
axa.plot(lim, lim, ":", color="0.4", lw=1.3, zorder=1, label="1:1 (perfect model)")
adjust_text(_texts, ax=axa, only_move={"text": "xy"}, arrowprops=dict(arrowstyle="-", color="0.7", lw=0.5))
axa.set_xlabel(r"Measured CH$_4$ emission  [$\mu$mol m$^{-2}$ d$^{-1}$]", fontsize=13)
axa.set_ylabel(r"Modelled CH$_4$ emission  [$\mu$mol m$^{-2}$ d$^{-1}$]", fontsize=13)
axa.set_title("a  Metabolic model vs. measured methane  (n = 6 unrestored)", fontsize=13.5, loc="left", fontweight="bold")
annotate(axa, sa_mepn, C_MEPN, (0.03, 0.97))
annotate(axa, sa_bet, C_BET, (0.03, 0.74))
axa.legend(loc="lower right", fontsize=9, framealpha=0.9)

# ---- Panel (b): measured methane vs consumer abundance, ALL 12 cores ------------
order = list(ab_raw.columns)
xb_mepn = mepn_abund.reindex(order).values
xb_bet = bet_abund.reindex(order).values
yb = CH4_RAW.reindex(order).values
edge = [REST_COLORS[rest_class(p)] for p in order]     # restoration by marker edge colour
sb_mepn = fit_stats(xb_mepn, yb)
sb_bet = fit_stats(xb_bet, yb)
for xb, s, c, lab, mk in [(xb_mepn, sb_mepn, C_MEPN, "MePn consumers", "o"),
                          (xb_bet, sb_bet, C_BET, "Betaine consumers", "s")]:
    axb.scatter(xb, yb, s=70, marker=mk, facecolor=c, edgecolor=edge, lw=1.6, zorder=3, label=lab)
    draw_linear_on_log(axb, xb, s, c)
log_axes(axb)
axb.set_xlabel("Summed consumer relative abundance", fontsize=13)
axb.set_ylabel(r"Measured CH$_4$ emission  [$\mu$mol m$^{-2}$ d$^{-1}$]", fontsize=13)
axb.set_title("b  Methane vs. consumer abundance  (all 12 cores)", fontsize=13.5, loc="left", fontweight="bold")
annotate(axb, sb_mepn, C_MEPN, (0.03, 0.97))
annotate(axb, sb_bet, C_BET, (0.03, 0.74))
series_h = [Line2D([0], [0], marker="o", color="w", mfc=C_MEPN, ms=9, label="MePn consumers"),
            Line2D([0], [0], marker="s", color="w", mfc=C_BET, ms=9, label="Betaine consumers")]
rings = [Line2D([0], [0], marker="o", mfc="0.9", mec=v, mew=1.6, ls="", ms=9, label=k) for k, v in REST_COLORS.items()]
leg1 = axb.legend(handles=series_h, loc="lower right", fontsize=9, framealpha=0.9)
axb.add_artist(leg1)
axb.legend(handles=rings, loc="lower right", bbox_to_anchor=(1.0, 0.16), fontsize=8,
           title="restoration (edge)", framealpha=0.9)

plt.tight_layout()
fig.savefig(OUT / "figure7_corrected.png", bbox_inches="tight")
fig.savefig(OUT / "figure7_corrected.svg", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- reconciliation numbers
def dnanorm_pair():
    """Panel (a) on the submitted figure's DNA-normalised basis, for reconciliation."""
    xmeas = CH4_NORM.reindex(UNRESTORED).values
    return fit_stats(xmeas, mepn_model.values), fit_stats(xmeas, bet_model.values)

dn_mepn, dn_bet = dnanorm_pair()
# panel (b) sensitivity: raw-vs-raw (corrected) vs DNA-normalised-both (submitted, inflated)
sb_mepn_unres = fit_stats(mepn_abund.reindex(UNRESTORED).values, CH4_RAW.reindex(UNRESTORED).values)
sb_bet_unres = fit_stats(bet_abund.reindex(UNRESTORED).values, CH4_RAW.reindex(UNRESTORED).values)
sb_mepn_dn_all = fit_stats(mepn_abund_dn.reindex(order).values, CH4_NORM.reindex(order).values)
sb_bet_dn_all = fit_stats(bet_abund_dn.reindex(order).values, CH4_NORM.reindex(order).values)
sb_mepn_dn_un = fit_stats(mepn_abund_dn.reindex(UNRESTORED).values, CH4_NORM.reindex(UNRESTORED).values)
sb_bet_dn_un = fit_stats(bet_abund_dn.reindex(UNRESTORED).values, CH4_NORM.reindex(UNRESTORED).values)

report = {
    "panel_a_modelled_vs_measured": {
        "raw_areal_basis_corrected": {"MePn": sa_mepn, "Betaine": sa_bet},
        "DNA_normalised_basis_as_submitted": {"MePn": dn_mepn, "Betaine": dn_bet},
    },
    "panel_b_methane_vs_abundance": {
        "raw_all_12_cores_corrected": {"MePn": sb_mepn, "Betaine": sb_bet},
        "raw_unrestored_only_n6": {"MePn": sb_mepn_unres, "Betaine": sb_bet_unres},
        "DNAnorm_both_axes_all12_as_submitted": {"MePn": sb_mepn_dn_all, "Betaine": sb_bet_dn_all},
        "DNAnorm_both_axes_unrestored_n6_as_submitted": {"MePn": sb_mepn_dn_un, "Betaine": sb_bet_dn_un},
    },
}
json.dump(report, open(OUT / "figure7_stats.json", "w"), indent=2, default=float)

print("=" * 74)
print("CORRECTED FIGURE 7  --  computed statistics")
print("=" * 74)
def show(title, d):
    print(f"\n{title}")
    for k, s in d.items():
        print(f"  {k:9s}  y={s['slope']:.3g}x+{s['intercept']:.3g}   "
              f"Pearson R2={s['pearson_r2']:.3f}   rho={s['spearman_rho']:.3f} "
              f"(p={s['spearman_p']:.3f}, n={s['n']})")
show("Panel a - modelled vs measured [CORRECTED raw areal basis]", {"MePn": sa_mepn, "Betaine": sa_bet})
show("Panel a - modelled vs measured [DNA-normalised, as submitted]", {"MePn": dn_mepn, "Betaine": dn_bet})
show("Panel b - methane vs abundance [CORRECTED raw, ALL 12 cores]", {"MePn": sb_mepn, "Betaine": sb_bet})
show("Panel b - methane vs abundance [raw, unrestored only, n=6]", {"MePn": sb_mepn_unres, "Betaine": sb_bet_unres})
show("Panel b - methane vs abundance [DNA-norm BOTH axes, all12  <-artifact]", {"MePn": sb_mepn_dn_all, "Betaine": sb_bet_dn_all})
show("Panel b - methane vs abundance [DNA-norm BOTH axes, n=6  <-text's 0.966/0.798]", {"MePn": sb_mepn_dn_un, "Betaine": sb_bet_dn_un})
print("\nSaved figure7_corrected.png/.svg and figure7_stats.json")
