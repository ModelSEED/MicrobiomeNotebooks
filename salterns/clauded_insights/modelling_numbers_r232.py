#!/usr/bin/env python
"""
Modelling statistics for the R232 review of the saltern manuscript (Results P69-P77, Discussion P95,
Conclusions P99, Abstract P18, Figure 7 caption P129).

Recomputes, from the STORED flux files only (no re-simulation), every number the modelling text quotes:
  * MePn models: correlational/saltern_fluxes.json at the phosphate ratios 0.07 and 0.33 (20 modelled MAGs),
    and the same restricted to the R232 phnJ carriers (17 of the 20);
  * the only stored run WITHOUT the phosphate-ratio constraint: correlational/saltern_fluxes_addition.json
    (analysis(1).ipynb cell 33, the substrate-addition simulation: no ratio constraint, MePn uptake <= 0.001);
  * betaine models: correlational/saltern_fluxes_betaine.json, CH4 = 0.06 x TMA (Methods P36), for the 15
    modelled MAGs and the 9 R232 grdH carriers (the betaine re-run in correlational/betaine_grdH_r232/
    showed the 9 re-simulate to the stored values, max |dTMA| 4.9e-15, so stored = re-simulated);
  * CH4 : MePn stoichiometry and realised MePn : Pi uptake ratios;
  * where the ratio constraint binds: it binds in 17 of the 20 models; Roseovarius.9, Roseovarius.11 and
    Albimonas.1 take up MePn above the minimum (same fluxes at 0.07 and 0.33). Their summed CH4 is reported as
    "inferred_unconstrained_nonbinding3" (an INFERENCE for the draft's unstored unconstrained run: on the old
    both-axes-per-DNA basis it gives R^2 = 0.9705, the draft's 0.970);
  * Soil DNA density vs CH4 (Figure S7: R^2 = 0.22, n = 12) and the shared-1/DNA diagnostic behind R^2 ~ 0.97.

Three comparison bases, n = 6 unrestored cores (R1_A-C, R2_A-C):
  raw_areal            x = measured areal CH4 (umol m-2 d-1); y = modelled areal CH4 (per-core volume-weighted
                       DNA, soil density 0.7 g cm-3), exactly as clauded_insights/figure7_corrected.py panel a.
  DNA_normalised_as_submitted
                       x = measured CH4 / DNA (depth mean); same y. This is the submitted Figure 7a.
  notebook_fits_json   x = measured CH4 / DNA; y = summed flux x a constant (2 ug DNA g-1, 1.4 g cm-3), i.e.
                       analysis(1).ipynb cell 27 -> fits.json / fits_betaine.json. Both axes carry 1/DNA. This is
                       the basis of the draft's R^2 = 0.978 / 0.968 / 0.862 and equations 0.51x+138, 0.95x+224,
                       0.24x+4632.

Asserts: 0.07 MePn and 15-MAG betaine stats equal clauded_insights/figure7_stats.json; notebook-basis R^2 equal
fits.json (0.97828062, 0.96812476), fits_betaine.json (0.86221939) and fits_addition.json (0.91854238).

Run:  ~/Documents/py_venv/bin/python clauded_insights/modelling_numbers_r232.py
Writes clauded_insights/modelling_numbers_r232.json
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
S = HERE.parent
DATA, CORR = S / "data", S / "correlational"
OUT = HERE / "modelling_numbers_r232.json"

# ---------------------------------------------------------------- inputs (as figure7_corrected.py)
DNA = pd.read_csv(DATA / "SaltPondsDNA.csv").set_index("Sample")["DNA conc (µg/g soil)"]
PREFIXES = sorted({"_".join(i.split("_")[:-1]) for i in DNA.index})
VW_DNA = pd.Series({p: (5 * DNA[p + "_D1"] + 10 * DNA[p + "_D2"]) / 15 for p in PREFIXES})
UNRESTORED = [p for p in PREFIXES if p.startswith(("R1_", "R2_"))]
assert UNRESTORED == ["R1_A", "R1_B", "R1_C", "R2_A", "R2_B", "R2_C"]

meta = pd.read_csv(DATA / "Cliff_Sample_Metadata_BGC_NMR.csv").set_index("Sample")
ch4 = meta["CH4_umol_m2_d"]
CH4_RAW = pd.Series({p: ch4[[p + "_D1", p + "_D2"]].mean() for p in PREFIXES})
CH4_NORM = pd.Series({p: (ch4 / DNA)[[p + "_D1", p + "_D2"]].mean() for p in PREFIXES})
METHANE_DIC = pd.Series(json.load(open(CORR / "averaged_normalized_methane_unrestored.json")))
assert np.allclose(METHANE_DIC.reindex(UNRESTORED), CH4_NORM.reindex(UNRESTORED), rtol=1e-12)
MEAN_INV_DNA = pd.Series({p: (1 / DNA[p + "_D1"] + 1 / DNA[p + "_D2"]) / 2 for p in PREFIXES})

AREA_M2 = math.pi * (5 / 2) ** 2 * 0.01 ** 2
VOLUME_CM3 = 294.52


def to_areal(summed: pd.Series) -> pd.Series:
    """figure7_corrected.fluxes_to_emissions: per-core VW DNA, soil density 0.7 (Methods P40)."""
    g_bio = 0.7 * VOLUME_CM3 * VW_DNA.reindex(summed.index).fillna(2.0) / 0.1 / 1e6
    return summed / AREA_M2 * 1000.0 * 24 * g_bio


def to_notebook(summed: pd.Series) -> pd.Series:
    """analysis(1).ipynb cell 3 fluxes_to_emissions: constant 2 ug DNA g-1 and 1.4 g cm-3."""
    return summed / AREA_M2 * 1e3 * 24 * (1.4 * VOLUME_CM3 * 2 / 0.1 / 1e6)


def fit(x, y) -> dict:
    x = np.asarray(x, float); y = np.asarray(y, float)
    lr = stats.linregress(x, y)
    rho, rp = stats.spearmanr(x, y)
    return dict(n=int(len(x)), slope=lr.slope, intercept=lr.intercept, pearson_r=lr.rvalue,
                pearson_r2=lr.rvalue ** 2, pearson_p=lr.pvalue, spearman_rho=rho, spearman_p=rp)


def summed(flux_by_mag: dict, mags=None) -> pd.Series:
    df = pd.DataFrame(flux_by_mag)              # rows = cores, cols = MAGs
    if mags is not None:
        df = df[[m for m in df.columns if m in mags]]
    return df.sum(axis=1).abs().reindex(UNRESTORED)   # missing (infeasible) cells count as 0, as in the notebook


def all_bases(s: pd.Series) -> dict:
    areal = to_areal(s)
    nb = to_notebook(s)
    xr, xn = CH4_RAW.reindex(UNRESTORED), CH4_NORM.reindex(UNRESTORED)
    ratio = areal / xr
    return {
        "raw_areal": fit(xr, areal),
        "DNA_normalised_as_submitted": fit(xn, areal),
        "notebook_fits_json": fit(xn, nb),
        "modelled_areal_CH4_by_core": areal.round(3).to_dict(),
        "modelled_notebook_basis_by_core": nb.round(3).to_dict(),
        "model_over_measured_areal_by_core": ratio.round(4).to_dict(),
        "model_over_measured_areal_summed": float(areal.sum() / xr.sum()),
        "model_over_measured_areal_median": float(ratio.median()),
    }


# ---------------------------------------------------------------- guilds / IDs
IDS = json.load(open(CORR / "iterativeIDs_r232.json"))
G = pd.read_csv(CORR / "r232" / "guilds_r232.csv").set_index("mag")
strip = lambda k: k.replace(".contigs", "")
rid = lambda k: IDS.get(strip(k), strip(k))

MEPN = json.load(open(CORR / "saltern_fluxes.json"))
ADD = json.load(open(CORR / "saltern_fluxes_addition.json"))
BET = json.load(open(CORR / "saltern_fluxes_betaine.json"))
BET_CH4 = {m: {c: 0.06 * v for c, v in d.items()} for m, d in BET["TMA"].items()}   # Methods P36

mepn_mags = list(MEPN["0.07"]["methane"])
phnj_mags = [m for m in mepn_mags if bool(G.loc[strip(m), "mepn_consumer"])]
bet_mags = list(BET["TMA"])
grdh_mags = [m for m in bet_mags if bool(G.loc[strip(m), "betaine_reducer"])]
assert len(mepn_mags) == 20 and len(phnj_mags) == 17 and len(bet_mags) == 15 and len(grdh_mags) == 9

report: dict = {"meta": {
    "n_cores": 6, "cores": UNRESTORED,
    "bases": {
        "raw_areal": "x = measured areal CH4 (mean of D1/D2 rows), y = modelled areal CH4 via per-core volume-weighted "
                     "DNA and 0.7 g cm-3 (figure7_corrected.py panel a)",
        "DNA_normalised_as_submitted": "x = measured CH4/DNA (depth mean), same y (submitted Figure 7a)",
        "notebook_fits_json": "x = measured CH4/DNA, y = summed flux x constant (2 ug DNA g-1, 1.4 g cm-3); "
                              "analysis(1).ipynb cell 27 -> fits.json; both axes carry 1/DNA"},
    "measured_CH4_raw": CH4_RAW.reindex(UNRESTORED).round(2).to_dict(),
    "measured_CH4_per_DNA": CH4_NORM.reindex(UNRESTORED).round(3).to_dict(),
    "VW_DNA_ug_g": VW_DNA.reindex(UNRESTORED).round(4).to_dict(),
    "mean_inverse_DNA_by_core": MEAN_INV_DNA.reindex(UNRESTORED).round(4).to_dict(),
}}

# ---------------------------------------------------------------- MePn: 0.07, 0.33, phnJ subsets, addition run
mepn = {}
for ratio in ("0.07", "0.33"):
    mepn[ratio] = {"n_MAGs": len(mepn_mags), **all_bases(summed(MEPN[ratio]["methane"]))}
    mepn[ratio + "_phnJ17"] = {"n_MAGs": len(phnj_mags), **all_bases(summed(MEPN[ratio]["methane"], phnj_mags))}
mepn["addition_no_ratio_constraint"] = {
    "n_MAGs_simulated": len(ADD["methane"]),
    "MAGs_with_CH4": sorted(rid(m) for m, d in ADD["methane"].items() if any(abs(v) > 1e-12 for v in d.values())),
    "cores_with_zero_modelled_CH4": [c for c, v in summed(ADD["methane"]).items() if v == 0],
    "note": "analysis(1).ipynb cell 33 (substrate-addition simulation): no phosphate-ratio constraint, MePn uptake "
            "<= 0.001; infeasible MAG x core cells were skipped and count as 0. The three CH4-producing models "
            "were infeasible in R2_B and R2_C. fits_addition.json R^2 = 0.918 is this run on the notebook basis "
            "with a NEGATIVE slope.",
    **all_bases(summed(ADD["methane"]))}
report["MePn"] = mepn
report["MePn_modelled_set"] = {
    "modelled_20": sorted(rid(m) for m in mepn_mags),
    "modelled_without_phnJ": sorted(rid(m) for m in mepn_mags if m not in phnj_mags),
    "curated_2020_not_modelled": {
        "no_MePn_model_built (archaea; models/MePn_consumers has no file, base RAST models exist in correlational/models)": ["Halorubrum.3", "Haloquadratum.4", "Natrinema.1"],
        "no_growth_on_lactate_medium (dontGrow)": ["Desulfotignum.1", "Desulfotignum.2", "Gaopeijia.1",
                                                   "Phaeodactylibacter.1"]},
    "phnJ_carriers_not_modelled": sorted(G.loc[G.mepn_consumer, "iterativeID_r232"].pipe(
        lambda s: s[~s.index.isin([strip(m) for m in mepn_mags])]).tolist()),
}

# share of modelled MePn CH4 by R232 genus / family (raw areal basis, summed over the 6 cores)
for ratio in ("0.07", "0.33"):
    df = pd.DataFrame(MEPN[ratio]["methane"]).reindex(UNRESTORED).abs()
    areal = df.apply(to_areal)
    by_mag = areal.sum()
    for rank in ("genus", "family"):
        lab = {m: G.loc[strip(m), rank] for m in by_mag.index}
        share = (by_mag.groupby(lab).sum() / by_mag.sum()).sort_values(ascending=False)
        report["MePn"][ratio][f"CH4_share_by_{rank}_areal"] = share.round(4).to_dict()

# Where the phosphate-ratio constraint binds.  MAGs whose MePn : Pi uptake exceeds the minimum take up MePn
# without being forced to; their fluxes are the same at 0.07 and 0.33.  INFERENCE: an in-situ run without the
# ratio constraint would leave these solutions optimal and give the other MAGs no reason to take up MePn, so
# the summed CH4 of these MAGs approximates the "unconstrained" run of P73 (which is not stored).
_pi = pd.DataFrame(MEPN["0.07"]["phosphate"]); _mp = pd.DataFrame(MEPN["0.07"]["MePn"])
_r = (_mp / _pi)
nonbinding = [m for m in mepn_mags if (_r[m] > 0.07 + 1e-6).all()]
binding = [m for m in mepn_mags if np.allclose(_r[m], 0.07, atol=1e-6)]
assert len(nonbinding) + len(binding) == len(mepn_mags)
_d = (pd.DataFrame(MEPN["0.07"]["methane"])[nonbinding] - pd.DataFrame(MEPN["0.33"]["methane"])[nonbinding]).abs().max().max()
_ab = pd.read_csv(CORR / "betaine_grdH_r232" / "inputs" / "averaged_normalized_abundances_unrestored_22f4629.csv", index_col=0)
_k = {rid(m): pd.Series(MEPN["0.07"]["methane"][m]).reindex(UNRESTORED) / _ab.loc[strip(m)].reindex(UNRESTORED)
      for m in mepn_mags}
report["MePn"]["ratio_constraint"] = {
    "binding_at_0.07_MAGs": len(binding),
    "nonbinding_MAGs": sorted(rid(m) for m in nonbinding),
    "nonbinding_MePn_to_Pi_uptake": {rid(m): round(float(_r[m].mean()), 4) for m in nonbinding},
    "nonbinding_CH4_identical_at_0.07_and_0.33_max_abs_diff": float(_d),
    "CH4_per_unit_abundance_constants": sorted({round(float(v.mean()), 5) for v in _k.values()}),
    "CH4_per_unit_abundance_max_CV_across_cores": float(max(v.std() / v.mean() for v in _k.values())),
    "abundance_basis": "correlational/betaine_grdH_r232/inputs/averaged_normalized_abundances_unrestored_22f4629.csv",
    "note": "Modelled CH4 of every MePn MAG is a constant times its (DNA-normalised) abundance in all 6 cores.",
}
s_nb = summed(MEPN["0.07"]["methane"], nonbinding)
mepn["inferred_unconstrained_nonbinding3"] = {
    "n_MAGs": len(nonbinding), "MAGs": sorted(rid(m) for m in nonbinding),
    "share_of_0.07_CH4_areal": float(to_areal(s_nb).sum() / to_areal(summed(MEPN["0.07"]["methane"])).sum()),
    "note": "INFERENCE, not a stored run: summed CH4 of the MAGs whose ratio constraint is not binding",
    **all_bases(s_nb)}

# ---------------------------------------------------------------- betaine (CH4 = 0.06 x TMA)
report["betaine"] = {
    "orig15": {"n_MAGs": 15, **all_bases(summed(BET_CH4))},
    "grdH9": {"n_MAGs": 9, "MAGs": sorted(rid(m) for m in grdh_mags), **all_bases(summed(BET_CH4, grdh_mags))},
    "note": "stored fluxes; the 9 grdH carriers re-simulate to these values exactly "
            "(correlational/betaine_grdH_r232/stats_grdH_r232.json)",
}
fva = pd.read_csv(CORR / "betaine_grdH_r232" / "fva_tma_by_core_r232.csv")
f9 = fva[fva.grdH]
reg = f9.groupby("iterativeID_r232")["regime"].agg(lambda s: sorted(set(s)))
sat_all = [m for m, r in reg.items() if r == ["saturated_at_0.1_uptake_bound"]]
zero_any = sorted(set(f9.loc[f9.TMA_existing.abs() < 1e-9, "iterativeID_r232"]))
report["betaine"]["grdH9_regimes"] = {
    "regimes_by_MAG": reg.to_dict(),
    "fixed_at_0.1_bound_in_all_6_cores": sat_all,
    "families_of_those": sorted({G.loc[G.iterativeID_r232 == m, "family"].iloc[0] for m in sat_all}),
    "zero_TMA_in_stored_solution_any_core": zero_any,
    "cells_at_bound_stored": int((f9.TMA_existing >= 0.1 - 1e-9).sum()), "cells": int(len(f9)),
    "cells_zero_stored": int((f9.TMA_existing.abs() < 1e-9).sum()),
}
f15 = fva
report["betaine"]["orig15_regimes"] = {
    "fixed_at_0.1_bound_in_all_6_cores": [m for m, r in f15.groupby("iterativeID_r232")["regime"].agg(
        lambda s: sorted(set(s))).items() if r == ["saturated_at_0.1_uptake_bound"]],
    "zero_TMA_in_stored_solution_all_cores": sorted(m for m, s in f15.groupby("iterativeID_r232")["TMA_existing"]
                                                    if (s.abs() < 1e-9).all()),
}
bs = json.load(open(CORR / "betaine_grdH_r232" / "stats_grdH_r232.json"))["betaine"]
for k_mine, k_theirs in (("orig15", "orig15_existing"), ("grdH9", "grdH9_existing")):
    for b_mine, b_theirs in (("raw_areal", "raw_areal"), ("DNA_normalised_as_submitted", "DNA_normalised_as_submitted"),
                             ("notebook_fits_json", "notebook_basis_P74")):
        assert abs(report["betaine"][k_mine][b_mine]["pearson_r2"] - bs[k_theirs][b_theirs]["pearson_r2"]) < 1e-10

# ---------------------------------------------------------------- reproduction checks
f7 = json.load(open(HERE / "figure7_stats.json"))["panel_a_modelled_vs_measured"]
chk = {}
for basis_f7, basis in (("raw_areal_basis_corrected", "raw_areal"), ("DNA_normalised_basis_as_submitted",
                                                                       "DNA_normalised_as_submitted")):
    for ser, mine in (("MePn", mepn["0.07"][basis]), ("Betaine", report["betaine"]["orig15"][basis])):
        d = max(abs(f7[basis_f7][ser][k] - mine[k]) for k in ("slope", "intercept", "pearson_r2", "spearman_rho",
                                                              "spearman_p"))
        chk[f"figure7_stats/{basis_f7}/{ser}"] = d
        assert d < 1e-9, (basis_f7, ser, d)
r2txt = lambda f: float(json.load(open(CORR / f))["R2"].split("=")[1]) if "R2" in json.load(open(CORR / f)) else None
fits = json.load(open(CORR / "fits.json"))
for ratio in ("0.07", "0.33"):
    stored = float(fits[ratio]["R2"].split("=")[1])
    chk[f"fits.json[{ratio}]"] = abs(stored - mepn[ratio]["notebook_fits_json"]["pearson_r2"])
    assert chk[f"fits.json[{ratio}]"] < 1e-8
    eq = fits[ratio]["eq"]; s = mepn[ratio]["notebook_fits_json"]
    assert eq == f"y = {s['slope']:.2f}x + {s['intercept']:.2f}", (eq, s["slope"], s["intercept"])
chk["fits_betaine.json"] = abs(r2txt("fits_betaine.json") - report["betaine"]["orig15"]["notebook_fits_json"]["pearson_r2"])
chk["fits_addition.json"] = abs(r2txt("fits_addition.json") -
                                mepn["addition_no_ratio_constraint"]["notebook_fits_json"]["pearson_r2"])
assert chk["fits_betaine.json"] < 1e-8 and chk["fits_addition.json"] < 1e-8
report["reproduction_max_abs_diff"] = chk

# ---------------------------------------------------------------- where the draft's numbers come from
nb = lambda d: round(d["notebook_fits_json"]["pearson_r2"], 3)
report["draft_values_traced"] = {
    "P69 R2=0.966 (MePn, 0.07)": "matches no stored fit at any commit; fits.json[0.07] = %.3f" % nb(mepn["0.07"]),
    "P95 R2=0.978 (MePn)": "fits.json[0.07] = %.5f (notebook basis)" % mepn["0.07"]["notebook_fits_json"]["pearson_r2"],
    "P73 R2=0.968 (ratio 0.33)": "fits.json[0.33] = %.5f (notebook basis)" % mepn["0.33"]["notebook_fits_json"]["pearson_r2"],
    "P73 R2=0.970 (unconstrained)": "no unconstrained in-situ run is stored at any commit (ratio loop was only "
                                     "[0.07] then [0.07, 0.33]); the only run without the ratio constraint is the "
                                     "substrate-addition run: R2 = %.3f with r = %.2f (inverse). INFERENCE: the three MAGs whose ratio "
                                     "constraint is not binding (Roseovarius.9, Roseovarius.11, Albimonas.1; the "
                                     "draft's 'two Roseovarius and one Albimonas') give R2 = %.4f on this basis" % (
                                         nb(mepn["addition_no_ratio_constraint"]),
                                         mepn["addition_no_ratio_constraint"]["notebook_fits_json"]["pearson_r"],
                                         mepn["inferred_unconstrained_nonbinding3"]["notebook_fits_json"]["pearson_r2"]),
    "P74 R2=0.862 (betaine)": "fits_betaine.json = %.5f (notebook basis, 15 MAGs)" % report["betaine"]["orig15"]["notebook_fits_json"]["pearson_r2"],
    "P70 eq 0.51x+138": "fits.json[0.07] notebook basis", "P72 eq 0.95x+224": "fits.json[0.33] notebook basis",
    "P75 eq 0.24x+4632": "betaine 15 MAGs, notebook basis (fits_betaine.json stores a different equation, "
                         "y = -3.95x - 77202.90, from an earlier cell state)",
    "soil_density": "notebook basis used 1.4 g cm-3 and a constant 2 ug DNA g-1; Methods P40 states 0.7",
}

# ---------------------------------------------------------------- stoichiometry and realised P-source ratios
sto = {}
for ratio in ("0.07", "0.33"):
    me = pd.DataFrame(MEPN[ratio]["methane"]); mp = pd.DataFrame(MEPN[ratio]["MePn"]); pi = pd.DataFrame(MEPN[ratio]["phosphate"])
    per = (me / -mp).replace([np.inf, -np.inf], np.nan)
    sto[ratio] = {"CH4_per_MePn_min": float(np.nanmin(per.values)), "CH4_per_MePn_max": float(np.nanmax(per.values)),
                  "community_CH4_per_MePn_by_core": (me.sum(axis=1) / -mp.sum(axis=1)).reindex(UNRESTORED).round(6).to_dict(),
                  "realised_MePn_to_Pi_uptake_by_core": (mp.sum(axis=1) / pi.sum(axis=1)).reindex(UNRESTORED).round(4).to_dict(),
                  "MePn_share_of_total_P_by_core": (mp.sum(axis=1) / (mp.sum(axis=1) + pi.sum(axis=1))).reindex(UNRESTORED).round(4).to_dict()}
me = pd.DataFrame(ADD["methane"]); mp = pd.DataFrame(ADD["MePn"])
sto["addition"] = {"community_CH4_per_MePn_by_core": {c: (None if pd.isna(v) else v) for c, v in
                   (me.sum(axis=1) / -mp.sum(axis=1)).reindex(UNRESTORED).round(6).items()}}   # None = no MePn used
sto["note"] = ("every stored MePn solution (and every committed version of saltern_fluxes.json) gives exactly 1 CH4 per "
               "MePn; the draft's 1.4 : 1 is not supported. Betaine CH4 = 0.06 x TMA is an assumption (Methods P36), "
               "not a model output.")
report["stoichiometry"] = sto

# ---------------------------------------------------------------- shared-denominator diagnostic and DNA vs CH4
xn = CH4_NORM.reindex(UNRESTORED); inv = MEAN_INV_DNA.reindex(UNRESTORED)
report["diagnostics"] = {
    "measured_CH4_per_DNA_vs_mean_1_over_DNA_unrestored": fit(inv, xn),
    "MePn_0.07_notebook_y_vs_mean_1_over_DNA_unrestored": fit(inv, pd.Series(mepn["0.07"]["modelled_notebook_basis_by_core"]).reindex(UNRESTORED)),
    "DNA_vs_CH4_12_cores_VW_DNA": fit(VW_DNA.reindex(PREFIXES), CH4_RAW.reindex(PREFIXES)),
    "DNA_vs_CH4_24_samples": fit(DNA.reindex(ch4.index), ch4),
}

json.dump(report, open(OUT, "w"), indent=1, default=float)

# ---------------------------------------------------------------- printout
def line(name, d):
    return (f"  {name:28s} R2={d['pearson_r2']:.3f} r={d['pearson_r']:+.2f} (p={d['pearson_p']:.3f})  "
            f"rho={d['spearman_rho']:+.2f} (p={d['spearman_p']:.2f})  y={d['slope']:.3g}x{d['intercept']:+.4g}")
for grp, dd in (("MePn", mepn), ("betaine", {k: v for k, v in report["betaine"].items() if k in ("orig15", "grdH9")})):
    for k, d in dd.items():
        if "raw_areal" not in d:
            continue
        print(f"{grp} {k}  (model/measured areal: summed {d['model_over_measured_areal_summed']:.2f}, "
              f"median {d['model_over_measured_areal_median']:.2f})")
        for b in ("raw_areal", "DNA_normalised_as_submitted", "notebook_fits_json"):
            print(line(b, d[b]))
print("stoichiometry 0.07:", sto["0.07"]["CH4_per_MePn_min"], sto["0.07"]["CH4_per_MePn_max"])
print("reproduction:", {k: f"{v:.1e}" for k, v in chk.items()})
print("wrote", OUT)
