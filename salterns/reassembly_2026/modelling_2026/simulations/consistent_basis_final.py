#!/usr/bin/env python
"""Simulated vs measured methane in the six unrestored cores for the 2026 guild models (final MAG set), on the
dimensionally consistent areal basis of clauded_insights/modelling_consistent_basis_r232.py:

- MePn consumers: each model's CH4 per unit of simulated abundance, k = flux / a (a = mean over depths of relative
  abundance / DNA, the simulation input), taken from the reference run at a fixed biomass of 1e-3 (reference_k.json;
  k is a per-model constant in this protocol, asserted against each model's highest-abundance core; at core abundances
  of about 1e-4 and below the fluxes fall under GLPK's feasibility tolerance, so one cell is 'infeasible' and others
  return 0), times its biomass in the core, i.e. the sum over the two layers of
  r x DNA x layer thickness, converted with the pipeline's constants (0.7 g soil cm-3, DNA = 10 % of biomass, 5 cm
  core, 5 + 10 cm layers).
- Betaine reducers: methane = 0.06 x TMA (the original assumption). TMA export is mostly set by the 0.1 betaine
  uptake bound, so the primary basis keeps the pipeline's conversion of the summed flux with the core's
  volume-weighted DNA ('betaine_18'); the per-core rescaling is the sensitivity 'betaine_18_rescaled', and the FVA
  minimum / maximum of TMA give 'betaine_18_fva_min' / '_fva_max'.
Sensitivities: 'relab_only' (constant DNA = the six-core mean of the volume-weighted DNA), MePn models with the RAST
methane monooxygenase closed ('_mmo_blocked'), and the MePn models whose C-P lyase is fully gene-supported ('_cplyase29').
Inputs: fluxes_final.json (simulate_guilds.py, run on TF), ../../abundance_final/final585_relabund_24sediment.tsv,
../models/models_summary.tsv, ../../../correlational/final585/{iterativeIDs_final585.json,guilds_final585.csv},
../../../data/{SaltPondsDNA.csv,Cliff_Sample_Metadata_BGC_NMR.csv}.
Run: ~/Documents/py_venv/bin/python consistent_basis_final.py -> consistent_basis_final.json
"""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
S = HERE.parents[2]
CORES = ["R1_A", "R1_B", "R1_C", "R2_A", "R2_B", "R2_C"]
LAYER = {"D1": 5, "D2": 10}
AREA_M2 = math.pi * (5 / 2) ** 2 * 1e-4
K = 0.7 * AREA_M2 * 1e4 / 0.1 / 1e6 / AREA_M2 * 1000 * 24      # as modelling_consistent_basis_r232.py
DNA = pd.read_csv(S / "data/SaltPondsDNA.csv").set_index("Sample")["DNA conc (µg/g soil)"]
REL = pd.read_csv(S / "reassembly_2026/abundance_final/final585_relabund_24sediment.tsv", sep="\t", index_col=0)
meta = pd.read_csv(S / "data/Cliff_Sample_Metadata_BGC_NMR.csv").set_index("Sample")["CH4_umol_m2_d"]
CH4 = pd.Series({c: meta[[c + "_D1", c + "_D2"]].mean() for c in CORES})
VW = pd.Series({c: (5 * DNA[c + "_D1"] + 10 * DNA[c + "_D2"]) / 15 for c in CORES})
IDS = json.load(open(S / "correlational/final585/iterativeIDs_final585.json"))
G = pd.read_csv(S / "correlational/final585/guilds_final585.csv").set_index("mag")
MS = pd.read_csv(HERE.parent / "models/models_summary.tsv", sep="\t")


def biomass(mag, core, dna=None):
    return sum(REL.loc[mag, f"{core}_{d}"] * (dna if dna is not None else DNA[f"{core}_{d}"]) * cm for d, cm in LAYER.items()) * K


A_MIN = 1e-4
KCHECK = {}


REF = json.load(open(HERE / "reference_k.json"))


def k_model(by_core, key, scale):
    if key == "methane":                                  # reference-biomass k (MePn); checked in main()
        rec = by_core["_ref"]
        return abs(rec["methane"]) * scale / REF["reference_biomass"], 0.0, 1
    a = pd.Series({c: by_core[c]["abundance"] for c in CORES})
    f = pd.Series({c: by_core[c][key] for c in CORES}).abs() * scale
    ok = (a >= A_MIN) & f.notna()
    kk = f[ok] / a[ok]
    if not len(kk):
        kk = (f / a)[(a > 0) & f.notna()]
    k = float(kk.median()) if len(kk) else 0.0
    cv = float(kk.std() / kk.mean()) if len(kk) > 1 and kk.mean() > 0 else 0.0
    return k, cv, int(ok.sum())


def simulated(flux_by_mag, scale=1.0, dna=None, key="methane"):
    """flux_by_mag: mag -> core -> record; returns areal emission per core and per-model contributions, plus k."""
    per, ks = {}, {}
    for mag, by_core in flux_by_mag.items():
        kv, cv, n = k_model(by_core, key, scale)
        KCHECK[(IDS[mag], key)] = (kv, cv, n)
        k = pd.Series({c: kv for c in CORES})
        ks[IDS[mag]] = k
        per[IDS[mag]] = pd.Series({c: k[c] * biomass(mag, c, dna) for c in CORES})
    return pd.DataFrame(per).sum(axis=1), pd.DataFrame(per), pd.DataFrame(ks)


def fit(sim):
    lr = stats.linregress(CH4[CORES], sim[CORES])
    rho, p = stats.spearmanr(CH4[CORES], sim[CORES])
    ratio = sim[CORES] / CH4[CORES]
    return dict(n=6, slope=lr.slope, intercept=lr.intercept, pearson_r2=lr.rvalue ** 2, pearson_p=lr.pvalue,
                spearman_rho=rho, spearman_p=p, model_over_measured_total=sim.sum() / CH4.sum(),
                model_over_measured_min=ratio.min(), model_over_measured_max=ratio.max(),
                sim_by_core=sim[CORES].round(3).to_dict())


def main():
    fx = json.load(open(HERE / "fluxes_final.json"))
    mepn, bet = fx["mepn"], fx["betaine"]
    for g in (mepn, bet):
        for mag, d in g.items():
            runs = [d[r] for r in ("0.07", "0.33")] if g is mepn else [d]
            for run in runs:
                for c, x in run.items():
                    assert x["status"] == "optimal" or x["abundance"] < A_MIN, (mag, c, x["status"])
    cpl = set(MS[(MS.guild == "mepn_consumer") & (MS.kept == True) & (MS.cp_lyase_gene_supported.astype(str) == "True")].mag_id) \
        if "cp_lyase_gene_supported" in MS else set()
    out = {"basis": __doc__.split("on the\ndimensionally")[0], "measured_CH4": CH4.round(2).to_dict(),
           "VW_DNA": VW.round(4).to_dict(), "n_mepn_models": len(mepn), "n_betaine_models": len(bet),
           "areal": {}, "relab_only": {}, "k_methane_per_abundance": {}, "betaine_tma": {}}
    mean_dna = VW.mean()
    sets = {}
    ref = REF["mepn"]
    assert set(ref) == set(mepn)
    kcheck = []
    for m, d in mepn.items():                              # reference k == k of the highest-abundance core
        for r in ("0.07", "0.33"):
            a = {c: x["abundance"] for c, x in d[r].items()}
            top = max(a, key=a.get)
            if a[top] >= 5e-4:
                kc, kr = abs(d[r][top]["methane"]) / a[top], abs(ref[m][r]["methane"]) / REF["reference_biomass"]
                kcheck.append(abs(kc - kr) / kr if kr > 0 else abs(kc))
    assert max(kcheck) < 1e-4, max(kcheck)
    out["k_reference_vs_top_core_max_rel_diff"] = max(kcheck)

    def with_ref(m, d, r, blocked=False):
        x = dict(d["mmo_blocked"][r] if blocked and "mmo_blocked" in d else d[r])
        x["_ref"] = ref[m]["mmo_blocked"][r] if blocked and "mmo_blocked" in ref[m] else ref[m][r]
        return x
    for r in ("0.07", "0.33"):
        sets[f"MePn_{r}_38"] = ({m: with_ref(m, d, r) for m, d in mepn.items()}, 1.0, "methane")
        sets[f"MePn_{r}_38_mmo_blocked"] = ({m: with_ref(m, d, r, True) for m, d in mepn.items()}, 1.0, "methane")
        if cpl:
            sets[f"MePn_{r}_cplyase{len(cpl)}"] = ({m: with_ref(m, d, r) for m, d in mepn.items() if m in cpl}, 1.0, "methane")
    sets["betaine_18"] = (bet, 0.06, "TMA")
    for name, (fl, sc, key) in sets.items():
        sim, per, ks = simulated(fl, sc, key=key)
        if name.startswith("MePn_") and name.endswith("_38"):
            kk = {i: KCHECK[(i, "methane")] for i in ks.columns}
            out["k_methane_per_abundance"][name] = {"by_model": {i: round(v[0], 6) for i, v in kk.items()},
                                                    "max_cv_across_cores_a_ge_1e-4": max(v[1] for v in kk.values()),
                                                    "models_without_core_a_ge_1e-4": sorted(i for i, v in kk.items() if v[2] == 0),
                                                    "min": min(v[0] for v in kk.values()), "max": max(v[0] for v in kk.values())}
            fam = {i: G.loc[m, "family"] for m, i in IDS.items() if i in per.columns}
            gen = {i: i.rsplit(".", 1)[0] for i in per.columns}
            out[f"{name}_family_share"] = (per.sum().groupby(fam).sum() / per.values.sum()).sort_values(ascending=False).round(4).to_dict()
            out[f"{name}_genus_share"] = (per.sum().groupby(gen).sum() / per.values.sum()).sort_values(ascending=False).round(4).to_dict()
            out[f"{name}_model_share"] = (per.sum() / per.values.sum()).sort_values(ascending=False).round(4).to_dict()
        if name == "betaine_18":
            out["areal"]["betaine_18_rescaled"] = fit(sim)
            stored = pd.DataFrame({m: {c: abs(d[c]["TMA"]) for c in CORES} for m, d in fl.items()}).T.sum()[CORES] * sc
            sim = stored * VW * 15 * K
            for b in ("fva_min", "fva_max"):
                st = pd.DataFrame({m: {c: abs(d[c].get(f"TMA_{b}", 0.0)) for c in CORES} for m, d in fl.items()}).T.sum()[CORES] * sc
                out["areal"][f"betaine_18_{b}"] = fit(st * VW * 15 * K)
            tma = pd.DataFrame({m: {c: d[c]["TMA"] for c in CORES} for m, d in fl.items()}).T
            lo = pd.DataFrame({m: {c: d[c].get("TMA_fva_min", np.nan) for c in CORES} for m, d in fl.items()}).T
            hi = pd.DataFrame({m: {c: d[c].get("TMA_fva_max", np.nan) for c in CORES} for m, d in fl.items()}).T
            out["betaine_tma"] = {"cells": int(tma.size), "at_bound_0.1": int((tma.abs() >= 0.1 - 1e-9).sum().sum()),
                                  "zero": int((tma.abs() < 1e-9).sum().sum()),
                                  "not_unique": int(((hi - lo).abs() > 1e-6).sum().sum()),
                                  "summed_TMA_fold_range": float(tma.abs().sum().max() / tma.abs().sum().min())}
        out["areal"][name] = fit(sim)
        if not name.startswith("betaine"):
            out["relab_only"][name] = fit(simulated(fl, sc, dna=mean_dna, key=key)[0])
    json.dump(out, open(HERE / "consistent_basis_final.json", "w"), indent=1, default=float)
    for name, s in out["areal"].items():
        print(f"{name:28s} model/measured {s['model_over_measured_total']:.3f} ({s['model_over_measured_min']:.3f}-"
              f"{s['model_over_measured_max']:.3f}) rho {s['spearman_rho']:+.3f} (p {s['spearman_p']:.2f}) R2 {s['pearson_r2']:.3f}")
    for name, s in out["relab_only"].items():
        print(f"relab_only {name:17s} model/measured {s['model_over_measured_total']:.3f} rho {s['spearman_rho']:+.3f}")
    for name, k in out["k_methane_per_abundance"].items():
        print(name, "k range", round(k["min"], 6), "-", round(k["max"], 6), "max CV (a>=1e-4)",
              f"{k['max_cv_across_cores_a_ge_1e-4']:.1e}", "no core >= 1e-4:", k["models_without_core_a_ge_1e-4"])
    print("genus share 0.07:", {k: v for k, v in out["MePn_0.07_38_genus_share"].items() if v > 0.02})
    print("family share 0.07:", {k: v for k, v in out["MePn_0.07_38_family_share"].items() if v > 0.02})
    print("betaine TMA:", out["betaine_tma"])


if __name__ == "__main__":
    main()
