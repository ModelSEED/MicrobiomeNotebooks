#!/usr/bin/env python
"""
Modelled betaine CH4 per unrestored core vs measured CH4 for the original 15-MAG betaine set
and the 9 grdH carriers, computed EXACTLY as clauded_insights/figure7_corrected.py does
(same fluxes_to_emissions, summed_flux_to_emissions, CH4_RAW / CH4_NORM and fit_stats; the
functions are copied verbatim because importing that module draws the figure).

Before anything else it reproduces the 15-MAG betaine and the MePn numbers stored in
clauded_insights/figure7_stats.json and stops if they differ.

Inputs : correlational/saltern_fluxes_betaine.json (existing), correlational/saltern_fluxes.json
         (MePn comparator), fluxes_grdH_r232.json (written by rerun_betaine_grdH_r232.py)
Output : stats_grdH_r232.json
Run    : ~/Documents/py_venv/bin/python betaine_stats_grdH_r232.py
"""
from __future__ import annotations

import json
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
CORR = HERE.parent
REPO = CORR.parent
DATA = REPO / "data"
FIG7_STATS = REPO / "clauded_insights" / "figure7_stats.json"

# ------------------------------------------------- verbatim from clauded_insights/figure7_corrected.py
DNA = pd.read_csv(DATA / "SaltPondsDNA.csv").set_index("Sample")["DNA conc (µg/g soil)"]
PREFIXES = sorted({"_".join(i.split("_")[:-1]) for i in DNA.index})
VW_DNA = pd.Series({p: (5 * DNA[p + "_D1"] + 10 * DNA[p + "_D2"]) / 15 for p in PREFIXES})


def fluxes_to_emissions(row: pd.Series) -> pd.Series:
    area_m2 = math.pi * (5 / 2) ** 2 * 0.01 ** 2
    volume_cm3 = 294.52
    ug_dna = VW_DNA.reindex(row.index).fillna(2.0)
    dna_per_biomass = 0.1
    g_soil_cm3 = 0.7
    g_biomass_g_soil = g_soil_cm3 * volume_cm3 * ug_dna / dna_per_biomass / 1e6
    return row / area_m2 * 1000.0 * 24 * g_biomass_g_soil


BETAINE_CH4_PER_TMA = 0.06
meta = pd.read_csv(DATA / "Cliff_Sample_Metadata_BGC_NMR.csv").set_index("Sample")
ch4 = meta["CH4_umol_m2_d"]
CH4_RAW = pd.Series({p: ch4[[p + "_D1", p + "_D2"]].mean() for p in PREFIXES})
ch4_norm = ch4 / DNA
CH4_NORM = pd.Series({p: ch4_norm[[p + "_D1", p + "_D2"]].mean() for p in PREFIXES})
UNRESTORED = [p for p in PREFIXES if p.startswith(("R1_", "R2_")) and "restored" not in p]


def summed_flux_to_emissions(flux_by_mag: dict) -> pd.Series:
    summed = pd.DataFrame(flux_by_mag).T.sum().abs()
    return fluxes_to_emissions(summed).reindex(UNRESTORED)


def fit_stats(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    lr = stats.linregress(x, y)
    rho, p = stats.spearmanr(x, y)
    return dict(n=int(len(x)), slope=lr.slope, intercept=lr.intercept,
                pearson_r2=lr.rvalue ** 2, spearman_rho=rho, spearman_p=p)
# ------------------------------------------------- end verbatim


def fit_stats_ext(x, y):
    """fit_stats + Pearson r and p; NaN instead of an exception for constant series."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if np.nanstd(np.asarray(y, float)) == 0 or np.nanstd(np.asarray(x, float)) == 0:
            return dict(n=int(np.isfinite(y).sum()), slope=np.nan, intercept=np.nan, pearson_r2=np.nan,
                        pearson_r=np.nan, pearson_p=np.nan, spearman_rho=np.nan, spearman_p=np.nan)
        s = fit_stats(x, y)
        lr = stats.linregress(np.asarray(x, float), np.asarray(y, float))
        s.update(pearson_r=lr.rvalue, pearson_p=lr.pvalue)
        return s


def betaine_ch4(tma_by_mag: dict) -> pd.Series:
    """modelled betaine CH4 per unrestored core = 0.06 x TMA, as figure7_corrected.py"""
    return summed_flux_to_emissions({m: {s: BETAINE_CH4_PER_TMA * v for s, v in d.items()}
                                     for m, d in tma_by_mag.items()})


# P74's R2 = 0.862 (fits_betaine.json) came from analysis(1).ipynb at 22f4629: x = DNA-normalised
# measured CH4 (averaged_normalized_methane_unrestored.json), y = summed betaine flux x a CONSTANT
# unit factor (ugDNA_gSoil = 2, gSoil_cm3 = 1.4).  Pearson R2 is invariant to that constant and
# to the 0.06 factor, so y = summed TMA reproduces it.
METHANE_DIC = pd.Series(json.load(open(CORR / "averaged_normalized_methane_unrestored.json")))


def notebook_basis(tma_by_mag):
    area_m2 = math.pi * (5 / 2) ** 2 * 1e-4
    k = 1 / area_m2 * 1e3 * 24 * (1.4 * 294.52 * 2 / 0.1 / 1e6)
    y = pd.DataFrame(tma_by_mag).T.sum().reindex(UNRESTORED) * k * BETAINE_CH4_PER_TMA
    return fit_stats_ext(METHANE_DIC.reindex(UNRESTORED).values, y.values)


def saturation(tma_by_mag: dict, abund: pd.DataFrame, modelled: pd.Series, bound=0.1, tol=1e-6):
    vals = np.array([v for d in tma_by_mag.values() for v in d.values()])
    summed = pd.DataFrame(tma_by_mag).T.sum().reindex(UNRESTORED)
    ab = abund.loc[[m.replace(".contigs", "") for m in tma_by_mag]].sum().reindex(UNRESTORED)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # round so that float noise between saturated sums (e.g. 0.7 vs 0.7000000000000001) is a tie
        rho_ab, p_ab = stats.spearmanr(np.round(summed.values, 10), ab.values)
        r_dna = stats.linregress(VW_DNA.reindex(UNRESTORED).values, modelled.values)
        rho_dna, _ = stats.spearmanr(VW_DNA.reindex(UNRESTORED).values, modelled.values)
    dna = VW_DNA.reindex(UNRESTORED)
    return dict(
        n_cells=int(vals.size), n_cells_at_0p1_bound=int((np.abs(vals - bound) < tol).sum()),
        frac_cells_at_0p1_bound=float((np.abs(vals - bound) < tol).mean()),
        n_cells_zero=int((np.abs(vals) < tol).sum()),
        summed_TMA_by_core=summed.to_dict(),
        summed_TMA_CV=float(summed.std(ddof=0) / summed.mean()) if summed.mean() else np.nan,
        summed_abundance_by_core_sim_basis=ab.to_dict(),
        summed_abundance_CV=float(ab.std(ddof=0) / ab.mean()),
        spearman_summedTMA_vs_summedAbundance=dict(rho=rho_ab, p=p_ab),
        modelled_CH4_vs_VW_DNA=dict(pearson_r2=r_dna.rvalue ** 2, spearman_rho=rho_dna),
        VW_DNA_CV=float(dna.std(ddof=0) / dna.mean()),
        max_over_min_across_cores=dict(summed_flux=float(summed.max() / summed.min()) if summed.min() > 0 else np.inf,
                                       summed_abundance=float(ab.max() / ab.min()),
                                       VW_DNA=float(dna.max() / dna.min())),
    )


def main():
    existing = json.load(open(CORR / "saltern_fluxes_betaine.json"))
    mepn = json.load(open(CORR / "saltern_fluxes.json"))["0.07"]
    fx = json.load(open(HERE / "fluxes_grdH_r232.json"))
    grdh9 = [m + ".contigs" for m in fx["meta"]["grdH9"]]
    mags15 = list(existing["TMA"])
    non6 = [m for m in mags15 if m not in grdh9]
    ab_run = pd.read_csv(HERE / "inputs" / "averaged_normalized_abundances_unrestored_22f4629.csv", index_col=0)
    ab_cur = pd.read_csv(DATA / "averaged_normalized_abundances_unrestored.csv", index_col=0)

    # ---------------------------------------------- 0. reproduce figure7_stats.json
    f7 = json.load(open(FIG7_STATS))["panel_a_modelled_vs_measured"]
    bet15 = betaine_ch4(existing["TMA"])
    mepn_model = summed_flux_to_emissions(mepn["methane"])
    checks = {
        ("raw_areal_basis_corrected", "Betaine"): fit_stats(CH4_RAW.reindex(UNRESTORED).values, bet15.values),
        ("DNA_normalised_basis_as_submitted", "Betaine"): fit_stats(CH4_NORM.reindex(UNRESTORED).values, bet15.values),
        ("raw_areal_basis_corrected", "MePn"): fit_stats(CH4_RAW.reindex(UNRESTORED).values, mepn_model.values),
        ("DNA_normalised_basis_as_submitted", "MePn"): fit_stats(CH4_NORM.reindex(UNRESTORED).values, mepn_model.values),
    }
    repro = {}
    for (basis, series), s in checks.items():
        ref = f7[basis][series]
        diff = max(abs(s[k] - ref[k]) for k in ("slope", "intercept", "pearson_r2", "spearman_rho", "spearman_p"))
        repro[f"{basis}/{series}"] = diff
        assert diff < 1e-9, (basis, series, diff)
    nb = notebook_basis(existing["TMA"])
    assert abs(nb["pearson_r2"] - 0.86221939) < 1e-8, nb  # fits_betaine.json "$R^2$ = 0.86221939"
    methane_dic_equals_CH4_NORM = bool(np.allclose(METHANE_DIC.reindex(UNRESTORED).values,
                                                   CH4_NORM.reindex(UNRESTORED).values, rtol=1e-12))
    print("reproduced figure7_stats.json (max abs diff):", repro)
    print("reproduced fits_betaine.json R2 =", round(nb["pearson_r2"], 8))

    # ---------------------------------------------- 1. sets
    pick = lambda d, mags: {m: d[m] for m in mags}
    fva_run = fx["fva_asrun_15"]
    sets = {
        "orig15_existing": (existing["TMA"], ab_run),
        "grdH9_existing": (pick(existing["TMA"], grdh9), ab_run),
        "grdH9_resim_same_settings_glpk": (fx["grdH9_resim_glpk_asrun"]["TMA"], ab_run),
        "grdH9_fva_min": (pick(fva_run["TMA_min"], grdh9), ab_run),
        "grdH9_fva_max": (pick(fva_run["TMA_max"], grdh9), ab_run),
        "nonGrdH6_existing": (pick(existing["TMA"], non6), ab_run),
        "orig15_fva_min": (fva_run["TMA_min"], ab_run),
        "orig15_fva_max": (fva_run["TMA_max"], ab_run),
    }
    if fx.get("resim_gurobi_asrun_15"):
        sets["grdH9_resim_gurobi"] = (pick(fx["resim_gurobi_asrun_15"]["TMA"], grdh9), ab_run)
    sets["SENS_grdH9_resim_glpk_current_abundance"] = (fx["grdH9_resim_glpk_current_abundance"]["TMA"], ab_cur)
    sets["SENS_grdH9_fva_min_current_abundance"] = (fx["grdH9_fva_current_abundance"]["TMA_min"], ab_cur)
    sets["SENS_grdH9_fva_max_current_abundance"] = (fx["grdH9_fva_current_abundance"]["TMA_max"], ab_cur)

    x_raw = CH4_RAW.reindex(UNRESTORED).values
    x_dn = CH4_NORM.reindex(UNRESTORED).values
    results = {}
    for name, (tma, ab) in sets.items():
        mod = betaine_ch4(tma)
        results[name] = {
            "n_MAGs": len(tma),
            "modelled_betaine_CH4_umol_m2_d_by_core": mod.to_dict(),
            "raw_areal": fit_stats_ext(x_raw, mod.values),
            "DNA_normalised_as_submitted": fit_stats_ext(x_dn, mod.values),
            "notebook_basis_P74": notebook_basis(tma),
            "saturation": saturation(tma, ab, mod),
        }

    # MePn comparator (existing saltern_fluxes.json, 0.07 ratio)
    mepn_sat = saturation(mepn["methane"], ab_run, mepn_model, bound=np.inf)
    for k in ("n_cells_at_0p1_bound", "frac_cells_at_0p1_bound"):
        mepn_sat.pop(k)
    mepn_sat["summed_CH4_flux_by_core"] = mepn_sat.pop("summed_TMA_by_core")
    mepn_sat["summed_CH4_flux_CV"] = mepn_sat.pop("summed_TMA_CV")
    mepn_sat["spearman_summedCH4flux_vs_summedAbundance"] = mepn_sat.pop("spearman_summedTMA_vs_summedAbundance")
    ratio = [v / ab_run.loc[m.replace(".contigs", ""), c] for m, d in mepn["methane"].items() for c, v in d.items()]
    mepn_out = {
        "n_MAGs": len(mepn["methane"]),
        "modelled_MePn_CH4_umol_m2_d_by_core": mepn_model.to_dict(),
        "raw_areal": fit_stats_ext(x_raw, mepn_model.values),
        "DNA_normalised_as_submitted": fit_stats_ext(x_dn, mepn_model.values),
        "abundance_response": mepn_sat,
        "CH4_flux_per_unit_abundance_distinct_values": sorted({round(r, 5) for r in ratio}),
    }

    out = {
        "meta": {
            "method": "modelled CH4 per core = fluxes_to_emissions(|sum_MAG 0.06 x TMA|) as in clauded_insights/figure7_corrected.py; "
                      "x = measured CH4 of the 6 unrestored cores (raw areal = CH4_RAW; DNA-normalised = CH4_NORM, as submitted). "
                      "notebook_basis_P74 = analysis(1).ipynb fit behind fits_betaine.json / P74 (x = DNA-normalised CH4, y = summed flux x constant).",
            "unrestored_cores": UNRESTORED,
            "measured_CH4_raw_umol_m2_d": CH4_RAW.reindex(UNRESTORED).to_dict(),
            "measured_CH4_DNA_normalised": CH4_NORM.reindex(UNRESTORED).to_dict(),
            "VW_DNA_ug_g": VW_DNA.reindex(UNRESTORED).to_dict(),
            "reproduction_of_figure7_stats_max_abs_diff": repro,
            "reproduction_of_fits_betaine_R2_0.86221939": nb["pearson_r2"],
            "methane_dic_equals_CH4_NORM": methane_dic_equals_CH4_NORM,
            "grdH9": fx["meta"]["grdH9"],
            "set_definitions": {
                "orig15_existing": "15 modelled betaine MAGs, TMA from saltern_fluxes_betaine.json",
                "grdH9_existing": "9 grdH carriers, TMA subset of saltern_fluxes_betaine.json",
                "grdH9_resim_same_settings_glpk": "9 grdH carriers re-simulated with the original settings (GLPK, as-run abundances); no structural model change was annotation-supported",
                "grdH9_fva_min / grdH9_fva_max": "lowest / highest TMA export compatible with the same optimal lactate uptake in every core (alternative optima)",
                "nonGrdH6_existing": "the 6 modelled MAGs without grdH (TMA from the hand-added, gene-less rxn17220_c0)",
                "orig15_fva_min / orig15_fva_max": "FVA envelope for the original 15",
                "grdH9_resim_gurobi": "solver sensitivity",
                "SENS_*_current_abundance": "sensitivity: same models/settings with the current data/averaged_normalized_abundances_unrestored.csv (DNA-weighted depth merge)",
            },
        },
        "betaine": results,
        "MePn_comparator_existing": mepn_out,
    }
    json.dump(out, open(HERE / "stats_grdH_r232.json", "w"), indent=2, default=float)

    print(f"\n{'set':46s} {'basis':28s} {'R2':>6s} {'rho':>6s} {'p_rho':>6s} {'p_r':>6s}  sat  rho(TMA,ab) R2(CH4~DNA)")
    for name, r in results.items():
        for basis in ("raw_areal", "DNA_normalised_as_submitted", "notebook_basis_P74"):
            s = r[basis]; sat = r["saturation"]
            print(f"{name:46s} {basis:28s} {s['pearson_r2']:6.3f} {s['spearman_rho']:6.3f} {s['spearman_p']:6.3f} "
                  f"{s['pearson_p']:6.3f}  {sat['frac_cells_at_0p1_bound']:.2f} "
                  f"{sat['spearman_summedTMA_vs_summedAbundance']['rho']:6.2f}      {sat['modelled_CH4_vs_VW_DNA']['pearson_r2']:.3f}")
    for basis in ("raw_areal", "DNA_normalised_as_submitted"):
        s = mepn_out[basis]
        print(f"{'MePn (existing)':46s} {basis:28s} {s['pearson_r2']:6.3f} {s['spearman_rho']:6.3f} {s['spearman_p']:6.3f} {s['pearson_p']:6.3f}")
    print("MePn R2(CH4~DNA) =", round(mepn_sat["modelled_CH4_vs_VW_DNA"]["pearson_r2"], 3),
          " rho(flux, abundance) =", round(mepn_sat["spearman_summedCH4flux_vs_summedAbundance"]["rho"], 3))
    print("\nwrote stats_grdH_r232.json")


if __name__ == "__main__":
    main()
