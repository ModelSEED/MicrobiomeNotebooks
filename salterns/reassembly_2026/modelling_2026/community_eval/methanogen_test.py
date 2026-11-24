#!/usr/bin/env python
"""Can abundance-weighted methane production by the methanogens, given their substrates, regenerate the measured
methane of the six unrestored cores? Compared on the same footing with the MePn consumers. Everything is per core,
areal (per m2 over the 0-15 cm sampled: D1 0-5 cm, D2 5-15 cm).

Biomass (g dry biomass m-2): relative abundance x soil DNA x layer thickness x 0.7 g soil cm-3 / 0.1 (DNA = 10 % of
  biomass), the constants of the manuscript's conversion (clauded_insights/modelling_consistent_basis_r232.py).
Yields (CH4 per g biomass formed, k):
  MePn consumers: the 2026 models' reference run (simulations/reference_k.json; phosphate-source ratio 0.07 / 0.33).
  Methanogen (Methanolobus.1, the only mcrA species): the 2026 build is acetogenic (arm A, no CH4); the KOfam-arm
  builds (annotation_comparison/runs/{B,C,D}) are methanogenic on methanol. k = CH4 / bio1 at a fixed biomass of 1e-3 on
  methanol as the only substrate, lactate-free, no external acceptor, minimising methanol uptake (as the protocol
  minimises lactate). Their maximum specific methanogenesis rate (CH4 at maximum growth) is recorded too.
The manuscript's conversion multiplies k by standing biomass per hour, i.e. it assumes every member grows at 1 h-1.
Here the rate is explicit instead: for each guild, the specific growth rate (and, for the methanogen, the specific
methanogenesis rate q) its measured biomass would need to supply the measured CH4.
Substrates: porewater NMR concentrations (uM) of methanol, mono-/di-/trimethylamine, x porewater volume (porosity
  1 - 0.7/2.65 = 0.736, from the same bulk density) x the methanogenesis stoichiometry (CH4 per substrate: methanol
  0.75, MMA 0.75, DMA 1.5, TMA 2.25) gives the standing CH4 potential per m2, and the turnover time at which it would
  supply the measured CH4. Betaine is reported separately (its methanogenic use by Methanolobus is not established).
Ranking: Spearman of measured CH4 with methanogen biomass, and with methanogen biomass x CH4 potential (substrate-
  weighted, the linear Michaelis-Menten regime), exact permutation p, n = 6.
Run: ~/Documents/py_venv/bin/python methanogen_test.py -> methanogen_test.json, methanogen_test_by_core.tsv
"""
import itertools
import json
from pathlib import Path

import cobra
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
M = HERE.parent
S = HERE.parents[2]
CORES = ["R1_A", "R1_B", "R1_C", "R2_A", "R2_B", "R2_C"]
LAYER = {"D1": 5, "D2": 10}
BULK, PARTICLE = 0.7, 2.65
POROSITY = 1 - BULK / PARTICLE
meta = pd.read_csv(S / "data/Cliff_Sample_Metadata_BGC_NMR.csv").set_index("Sample")
dna = pd.read_csv(S / "data/SaltPondsDNA.csv").set_index("Sample")["DNA conc (µg/g soil)"]
rel = pd.read_csv(S / "reassembly_2026/abundance_final/final585_relabund_24sediment.tsv", sep="\t", index_col=0)
g = pd.read_csv(S / "correlational/final585/guilds_final585.csv").set_index("mag")
ref = json.load(open(M / "simulations/reference_k.json"))
ch4 = pd.Series({c: meta.loc[[c + "_D1", c + "_D2"], "CH4_umol_m2_d"].mean() for c in CORES})
YIELD = {"Methanol": 0.75, "Methylamine": 0.75, "Dimethylamine": 1.5, "Trimethylamine": 2.25}


def biomass_g_m2(mag, core):
    return sum(rel.loc[mag, f"{core}_{d}"] * dna[f"{core}_{d}"] * cm * BULK * 1e4 / 0.1 / 1e6 for d, cm in LAYER.items())


def pool_umol_m2(col, core):
    return sum(float(meta.loc[f"{core}_{d}", col]) * POROSITY * cm * 1e4 / 1000 for d, cm in LAYER.items())


def exact(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    r0 = stats.spearmanr(x, y)[0]
    p = sum(abs(stats.spearmanr(x, y[list(q)])[0]) >= abs(r0) - 1e-12 for q in itertools.permutations(range(6))) / 720
    return float(r0), float(p)


def methanogen_k():
    from importlib.util import spec_from_file_location, module_from_spec
    sp = spec_from_file_location("bm", str(M / "pipeline/build_model.py"))
    bm = module_from_spec(sp); sp.loader.exec_module(bm)
    med = {k: v for k, v in bm.METHYL_MEDIUM.items() if k not in ("EX_cpd00187_e0", "EX_cpd00425_e0", "EX_cpd00441_e0")}
    out = {}
    for arm in "ABCD":
        f = M / f"annotation_comparison/runs/{arm}/R1_A_D1__concoct.bin.007__methanogen__{arm}.json"
        if arm == "A":
            f = M / "models/methanogen/R1_A_D1__concoct.bin.007.json"
        m = cobra.io.load_json_model(str(f))
        for k in med:
            if k not in m.reactions:
                m.add_boundary(m.metabolites.get_by_id(k[3:] .replace("_e0", "_e0")), type="exchange") if k[3:] in m.metabolites else None
        rec = {}
        with m:
            m.medium = {k: v for k, v in med.items() if k in m.reactions}
            m.objective = "bio1"
            sol = m.optimize()
            rec["mu_max_methanol"] = float(sol.objective_value) if sol.status == "optimal" else 0.0
            rec["q_ch4_at_mu_max"] = float(sol.fluxes.get("EX_cpd01024_e0", 0.0)) if sol.status == "optimal" else 0.0
        with m:
            m.medium = {k: v for k, v in med.items() if k in m.reactions}
            if rec["mu_max_methanol"] > 1e-3:
                m.reactions.bio1.bounds = (1e-3, 1e-3)
                m.objective = m.reactions.EX_cpd00116_e0          # maximise (least negative) = least methanol uptake
                m.reactions.EX_cpd00116_e0.upper_bound = 0
                s2 = m.optimize()
                rec["k_ch4_per_g"] = float(s2.fluxes["EX_cpd01024_e0"] / 1e-3) if s2.status == "optimal" else None
                rec["methanol_per_g"] = float(-s2.fluxes["EX_cpd00116_e0"] / 1e-3) if s2.status == "optimal" else None
            else:
                rec["k_ch4_per_g"] = None
        out[arm] = rec
    return out


def main():
    mk = methanogen_k()
    meth = list(g.index[g.methanogen])
    mepn = [m for m in ref["mepn"]]
    rows = []
    for c in CORES:
        Bm = sum(biomass_g_m2(m, c) for m in meth)
        cap07 = sum(abs(ref["mepn"][m]["0.07"]["methane"]) / ref["reference_biomass"] * biomass_g_m2(m, c) for m in mepn)
        cap33 = sum(abs(ref["mepn"][m]["0.33"]["methane"]) / ref["reference_biomass"] * biomass_g_m2(m, c) for m in mepn)
        pools = {s: pool_umol_m2(s, c) for s in list(YIELD) + ["Betaine"]}
        pot = sum(pools[s] * y for s, y in YIELD.items())
        meas = ch4[c]                                                   # umol m-2 d-1
        rows.append({"core": c, "CH4_measured_umol_m2_d": meas,
                     "methanogen_biomass_g_m2": Bm, "mepn_biomass_g_m2": sum(biomass_g_m2(m, c) for m in mepn),
                     # MePn: CH4 per hour of growth at 1 h-1 = cap (mmol m-2 per unit mu); required mu (h-1)
                     "mepn_mu_required_h_r0.07": meas / 1000 / 24 / cap07,
                     "mepn_mu_required_h_r0.33": meas / 1000 / 24 / cap33,
                     # methanogen: required specific methanogenesis rate q (mmol CH4 g-1 h-1)
                     "methanogen_q_required_mmol_g_h": meas / 1000 / 24 / Bm if Bm > 0 else np.inf,
                     "methyl_substrate_CH4_potential_umol_m2": pot,
                     "methyl_substrate_turnover_days_required": pot / meas,
                     "betaine_pool_umol_m2": pools["Betaine"],
                     **{f"pool_{s}_umol_m2": v for s, v in pools.items()}})
    t = pd.DataFrame(rows).set_index("core")
    for arm, rec in mk.items():
        if rec.get("k_ch4_per_g"):
            t[f"methanogen_mu_required_h_arm{arm}"] = t.methanogen_q_required_mmol_g_h / rec["k_ch4_per_g"]
    t.to_csv(HERE / "methanogen_test_by_core.tsv", sep="\t")
    rank = {"methanogen biomass": exact(t.methanogen_biomass_g_m2, t.CH4_measured_umol_m2_d),
            "methanogen biomass x methyl-substrate CH4 potential": exact(
                t.methanogen_biomass_g_m2 * t.methyl_substrate_CH4_potential_umol_m2, t.CH4_measured_umol_m2_d),
            "methyl-substrate CH4 potential alone": exact(t.methyl_substrate_CH4_potential_umol_m2, t.CH4_measured_umol_m2_d),
            "MePn biomass (= MePn model)": exact(t.mepn_biomass_g_m2, t.CH4_measured_umol_m2_d)}
    out = {"assumptions": {"bulk_density_g_cm3": BULK, "porosity": POROSITY, "DNA_fraction_of_biomass": 0.1,
                           "depth_cm": 15, "CH4_per_substrate": YIELD},
           "methanogen_models": mk, "ranking_spearman_exact": rank,
           "by_core": t.round(6).reset_index().to_dict(orient="records")}
    json.dump(out, open(HERE / "methanogen_test.json", "w"), indent=1, default=float)
    pd.set_option("display.width", 250)
    print("methanogen models (methanol only):", json.dumps(mk, indent=0, default=float))
    cols = ["CH4_measured_umol_m2_d", "methanogen_biomass_g_m2", "mepn_biomass_g_m2", "mepn_mu_required_h_r0.07",
            "mepn_mu_required_h_r0.33", "methanogen_q_required_mmol_g_h"] + [c for c in t.columns if c.startswith("methanogen_mu_required")] + [
            "methyl_substrate_CH4_potential_umol_m2", "methyl_substrate_turnover_days_required", "betaine_pool_umol_m2"]
    print(t[cols].round(4).T.to_string())
    print("ranking:", json.dumps(rank, indent=0))


if __name__ == "__main__":
    main()
