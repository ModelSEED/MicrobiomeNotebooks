#!/usr/bin/env python
"""
Simulated vs measured methane on a dimensionally consistent areal basis (supersedes the 'areal' basis of
figure7_corrected.py and modelling_numbers_r232.json for the model-vs-measured comparison).

The stored fluxes (correlational/saltern_fluxes.json, saltern_fluxes_betaine.json; analysis(1).ipynb at commit
22f4629) were simulated with each model's biomass fixed to a_i = mean over depths of (r_id / DNA_d): relative
abundance divided by each depth's DNA. Two depths have almost no DNA (R2_B_D1 and R2_C_D2, 0.16 ug/g soil), so a_i is
inflated in R2_B and R2_C. figure7_corrected.py multiplied the summed flux by the core's volume-weighted DNA, which
does not cancel that 1/DNA term: the factor VW_DNA x mean(1/DNA) alone has the same Spearman rho (0.60) against
measured methane as the 'areal' model, and the draft's R^2 = 0.978 / 0.862 come from the same term.

Consistent basis (used here): each MAG's flux per unit of simulated abundance, k_i = flux_i / a_i (a constant for the
MePn models, CV < 1e-8), times its biomass in the core: the sum over the two layers of r_id x DNA_d x layer volume,
converted with the pipeline's constants (0.7 g soil cm-3, DNA = 10 % of biomass, 5 cm core, 5 + 10 cm layers).
Betaine models: their TMA excretion is mostly set by the 0.1 betaine-uptake bound (33 of 54 grdH cells), not by
the biomass input, so the rescaling does not apply (it would give 0.1-27x the measured emissions). Their primary basis
keeps the pipeline's conversion of the summed flux with the core's volume-weighted DNA ('areal' for betaine); the
per-core rescaling is reported as the sensitivity 'betaine_rescaled'.
Sensitivity basis 'relab_only': the same with a constant DNA (the six-core mean), i.e. biomass-specific emission.

Run: ~/Documents/py_venv/bin/python modelling_consistent_basis_r232.py  -> modelling_consistent_basis_r232.json
"""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
S = HERE.parent
DATA, CORR = S / "data", S / "correlational"
CORES = ["R1_A", "R1_B", "R1_C", "R2_A", "R2_B", "R2_C"]
LAYER = {"D1": 5, "D2": 10}                                   # cm
AREA_M2 = math.pi * (5 / 2) ** 2 * 1e-4
K = 0.7 * AREA_M2 * 1e4 / 0.1 / 1e6 / AREA_M2 * 1000 * 24      # g soil/cm3 * cm2 / DNA-per-biomass ... per cm of layer

DNA = pd.read_csv(DATA / "SaltPondsDNA.csv").set_index("Sample")["DNA conc (µg/g soil)"]
REL = pd.read_csv(DATA / "Cliff_310MAG_relabund.tsv", sep="\t", index_col=0)
OLD = pd.read_csv(CORR / "betaine_grdH_r232" / "inputs" / "averaged_normalized_abundances_unrestored_22f4629.csv", index_col=0)
meta = pd.read_csv(DATA / "Cliff_Sample_Metadata_BGC_NMR.csv").set_index("Sample")["CH4_umol_m2_d"]
CH4 = pd.Series({c: meta[[c + "_D1", c + "_D2"]].mean() for c in CORES})
VW = pd.Series({c: (5 * DNA[c + "_D1"] + 10 * DNA[c + "_D2"]) / 15 for c in CORES})
GUILDS = pd.read_csv(CORR / "r232" / "guilds_r232.csv")
IDS = json.load(open(CORR / "iterativeIDs_r232.json"))


def key(name):
    return name.replace(".contigs", "").replace("bins.", "bins_")


RELK, OLDK = {key(m): m for m in REL.index}, {key(m): m for m in OLD.index}
IDK = {key(m): v for m, v in IDS.items()}


def biomass(mag, core, dna=None):
    """sum over layers of r x DNA x layer thickness, in the pipeline's areal units (per unit specific flux)."""
    m = RELK[key(mag)]
    return sum(REL.loc[m, f"{core}_{d}"] * (dna if dna is not None else DNA[f"{core}_{d}"]) * cm for d, cm in LAYER.items()) * K


def simulated(flux_by_mag, scale=1.0, dna=None):
    """Areal simulated emission per core (umol m-2 d-1) and per-MAG contributions."""
    per = {}
    for mag, by_core in flux_by_mag.items():
        a = OLD.loc[OLDK[key(mag)], CORES]
        f = pd.Series(by_core)[CORES].abs() * scale
        k = (f / a).where(a > 0, 0.0)
        per[IDK[key(mag)]] = pd.Series({c: k[c] * biomass(mag, c, dna) for c in CORES})
    per = pd.DataFrame(per)
    return per.sum(axis=1), per


def fit(sim):
    lr = stats.linregress(CH4[CORES], sim[CORES])
    rho, p = stats.spearmanr(CH4[CORES], sim[CORES])
    ratio = sim[CORES] / CH4[CORES]
    return dict(n=6, slope=lr.slope, intercept=lr.intercept, pearson_r2=lr.rvalue ** 2, pearson_p=lr.pvalue,
                spearman_rho=rho, spearman_p=p, model_over_measured_total=sim.sum() / CH4.sum(),
                model_over_measured_min=ratio.min(), model_over_measured_max=ratio.max(),
                sim_by_core=sim[CORES].round(3).to_dict())


def main():
    mepn = json.load(open(CORR / "saltern_fluxes.json"))
    bet = json.load(open(CORR / "saltern_fluxes_betaine.json"))["TMA"]
    phnj = {key(m) for m in GUILDS.loc[GUILDS.mepn_consumer, "mag"]}
    grdh = {key(m) for m in GUILDS.loc[GUILDS.betaine_reducer, "mag"]}
    sets = {
        "MePn_0.07_20": (mepn["0.07"]["methane"], 1.0),
        "MePn_0.33_20": (mepn["0.33"]["methane"], 1.0),
        "MePn_0.07_phnJ17": ({m: v for m, v in mepn["0.07"]["methane"].items() if key(m) in phnj}, 1.0),
        "MePn_0.07_nonbinding3": ({m: v for m, v in mepn["0.07"]["methane"].items()
                                   if IDK[key(m)] in ("Roseovarius.9", "Roseovarius.11", "Albimonas.1")}, 1.0),
        "betaine_grdH9": ({m: v for m, v in bet.items() if key(m) in grdh}, 0.06),
        "betaine_orig15": (bet, 0.06),
    }
    assert len(sets["MePn_0.07_phnJ17"][0]) == 17 and len(sets["MePn_0.07_nonbinding3"][0]) == 3
    assert len(sets["betaine_grdH9"][0]) == 9
    out = {"basis": __doc__.split("Consistent basis (used here):")[1].split("Run:")[0].strip(), "measured_CH4": CH4.round(2).to_dict(),
           "VW_DNA": VW.round(4).to_dict(), "areal": {}, "relab_only": {}, "stored_x_VW_superseded": {}}
    mean_dna = VW.mean()
    out["betaine_rescaled"] = {}
    for name, (fl, sc) in sets.items():
        sim, per = simulated(fl, sc)
        stored = pd.DataFrame(fl).T.abs().sum()[CORES] * sc
        if name.startswith("betaine"):
            out["betaine_rescaled"][name] = fit(sim)
            sim = stored * VW * 15 * K                     # pipeline conversion: summed flux x core biomass
        out["areal"][name] = fit(sim)
        out["relab_only"][name] = fit(simulated(fl, sc, dna=mean_dna)[0])
        out["stored_x_VW_superseded"][name] = {k: v for k, v in fit(stored * VW).items() if k in ("pearson_r2", "spearman_rho", "spearman_p")}
        if name == "MePn_0.07_20":
            fam = {i: GUILDS.set_index("iterativeID_r232").loc[i, "family"] for i in per.columns}
            out["MePn_0.07_20_family_share"] = (per.sum().groupby(fam).sum() / per.values.sum()).round(4).to_dict()
            out["MePn_0.07_20_genus_share"] = (per.sum().groupby(lambda i: i.rsplit(".", 1)[0]).sum() / per.values.sum()).round(4).to_dict()
    # the 1/DNA factor alone (what the superseded 'areal' basis carried)
    fac = VW * pd.Series({c: (1 / DNA[c + "_D1"] + 1 / DNA[c + "_D2"]) / 2 for c in CORES})
    rho, p = stats.spearmanr(CH4[CORES], fac[CORES])
    out["factor_VW_x_mean_inverse_DNA"] = {"spearman_rho": rho, "spearman_p": p, "by_core": fac.round(3).to_dict()}
    json.dump(out, open(HERE / "modelling_consistent_basis_r232.json", "w"), indent=1)
    for basis in ("areal", "relab_only", "betaine_rescaled", "stored_x_VW_superseded"):
        print(basis)
        for name, s in out[basis].items():
            extra = f"  model/measured {s['model_over_measured_total']:.2f} ({s['model_over_measured_min']:.2f}-{s['model_over_measured_max']:.2f})" if "model_over_measured_total" in s else ""
            print(f"  {name:24s} R2={s['pearson_r2']:.3f} rho={s['spearman_rho']:+.2f} p={s['spearman_p']:.2f}{extra}")
    print("family share 0.07:", out["MePn_0.07_20_family_share"])
    print("genus share 0.07:", {k: v for k, v in out["MePn_0.07_20_genus_share"].items() if v > 0.01})
    print("1/DNA factor alone:", out["factor_VW_x_mean_inverse_DNA"])


if __name__ == "__main__":
    main()
