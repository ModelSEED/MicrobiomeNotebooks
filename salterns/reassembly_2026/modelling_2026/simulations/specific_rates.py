#!/usr/bin/env python
"""Per-gram framing of the model-vs-measured comparison (replaces the areal conversion that assumed every member grows
at 1 h-1). For each guild and each unrestored core:
  yield k        = CH4 produced per g of biomass formed (per_gram_yields.tsv; methanogen from the curated model's
                   CURATION results when present, else the methanol-growing KOfam-arm builds, flagged);
  biomass B      = measured: relative abundance x soil DNA x layer thickness x 0.7 g soil cm-3 / 0.1 (DNA = 10 % of
                   dry biomass), summed over D1 (0-5 cm) and D2 (5-15 cm), per species, as measured (no adjustment);
  required mu    = measured CH4 / (sum over species of k x B), i.e. the specific growth rate at which the guild's measured
                   biomass would supply the measured flux by growth-coupled production; doubling time = ln 2 / mu.
Scenarios: MePn consumers at phosphate-source ratios 0.07 and 0.33 and with MePn as the only P source; betaine
reducers with TMA at the protocol optimum, methane = 0.06 x TMA (the manuscript's community conversion) or 2.25 x TMA
(complete methylotrophic conversion); the methanogen (Methanolobus.1, curated model) on methylamines.
Each scenario is repeated as "<name> (read-based abundance)": the guild's biomass in each sample is multiplied by the
read/MAG fold of its marker gene (mcrA, phnJ, grdI copies per genome equivalent in all reads / the guild MAGs' genomes
per genome equivalent; ../../read_markers/results/key_unrestored_cores.tsv), to count members not recovered as MAGs.
For the methanogen this assumes every mcrA carrier has Methanolobus's k (the required activity CH4 / B does not).
Ranking is unaffected by mu: Spearman of measured CH4 with each guild's k x B (exact permutation p, n = 6).
Run after per_gram_yields.py: ~/Documents/py_venv/bin/python specific_rates.py -> specific_rates.tsv, specific_rates.json
"""
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
M = HERE.parent
S = HERE.parents[2]
CORES = ["R1_A", "R1_B", "R1_C", "R2_A", "R2_B", "R2_C"]
LAYER = {"D1": 5, "D2": 10}
meta = pd.read_csv(S / "data/Cliff_Sample_Metadata_BGC_NMR.csv").set_index("Sample")
dna = pd.read_csv(S / "data/SaltPondsDNA.csv").set_index("Sample")["DNA conc (µg/g soil)"]
rel = pd.read_csv(S / "reassembly_2026/abundance_final/final585_relabund_24sediment.tsv", sep="\t", index_col=0)
g = pd.read_csv(S / "correlational/final585/guilds_final585.csv").set_index("mag")
ch4 = pd.Series({c: meta.loc[[c + "_D1", c + "_D2"], "CH4_umol_m2_d"].mean() for c in CORES})
Y0 = pd.read_csv(HERE / "per_gram_yields.tsv", sep="\t")
YM = Y0[Y0.guild == "mepn_consumer"].set_index("mag_id")      # UBA8670 has a model in both guilds
YB = Y0[Y0.guild == "betaine_reducer"].set_index("mag_id")


def biomass(mag, core, fold=None):                         # g dry biomass m-2; fold: per-sample read/MAG multiplier
    return sum(rel.loc[mag, f"{core}_{d}"] * dna[f"{core}_{d}"] * cm * 0.7 * 1e4 / 0.1 / 1e6
               * (fold[f"{core}_{d}"] if fold is not None else 1.0) for d, cm in LAYER.items())


# read-based abundance (../../read_markers/results/key_unrestored_cores.tsv): per sample, the guild's marker-gene copies
# per genome equivalent in all reads / its MAG genomes per genome equivalent (mcrA, phnJ, grdI)
RM = pd.read_csv(S / "reassembly_2026/read_markers/results/key_unrestored_cores.tsv", sep="\t")
FOLD = {gname: RM[RM.guild == gname].set_index("sample").fold_read_vs_MAG.astype(float).to_dict()
        for gname in ("methanogen", "mepn_consumer", "betaine_reducer")}


def methanogen_k():
    cur = M / "models/methanogen/curation_k.json"         # written by the curation (CURATION.md)
    if cur.exists():
        k = json.load(open(cur))
        return {"methylamines": float(k["methylamines"])}, "curated model (MMA, DMA and TMA give the same k)"
    mt = json.load(open(HERE.parent / "community_eval/methanogen_test.json"))["methanogen_models"]
    ks = [r["k_ch4_per_g"] for r in mt.values() if r.get("k_ch4_per_g")]
    return {"methanol": float(np.median(ks))}, "KOfam-arm builds on methanol (interim)"


def exact(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    r0 = stats.spearmanr(x, y)[0]
    p = sum(abs(stats.spearmanr(x, y[list(q)])[0]) >= abs(r0) - 1e-12 for q in itertools.permutations(range(6))) / 720
    return float(r0), float(p)


def main():
    mep, bet = list(YM.index), list(YB.index)
    meth = list(g.index[g.methanogen])
    km, ksrc = methanogen_k()
    scen = {
        "MePn, ratio 0.07": {m: YM.loc[m, "k_CH4_r0.07"] for m in mep},
        "MePn, ratio 0.33": {m: YM.loc[m, "k_CH4_r0.33"] for m in mep},
        "MePn, only P source": {m: YM.loc[m, "k_CH4_MePn_only_P"] for m in mep},
        "Betaine, 0.06 CH4 per TMA": {m: 0.06 * YB.loc[m, "k_TMA_unbounded_min"] for m in bet},
        "Betaine, 2.25 CH4 per TMA": {m: 2.25 * YB.loc[m, "k_TMA_unbounded_min"] for m in bet},
        **{f"Methanogen, {s}": {mm: v for mm in meth} for s, v in km.items()},
    }
    fold_of = {"MePn": "mepn_consumer", "Betaine": "betaine_reducer", "Methanogen": "methanogen"}
    for name in list(scen):
        scen[name + " (read-based abundance)"] = scen[name]
    rows, ranks = [], {}
    for name, kk in scen.items():
        fold = FOLD[fold_of[name.split(",")[0]]] if name.endswith("(read-based abundance)") else None
        cap = pd.Series({c: sum(k * biomass(m, c, fold) for m, k in kk.items()) for c in CORES})   # mmol CH4 m-2 per unit mu (h)
        mu = ch4 / 1000 / 24 / cap                                                            # h-1
        for c in CORES:
            rows.append({"scenario": name, "core": c, "biomass_g_m2": sum(biomass(m, c, fold) for m in kk), "CH4_measured_umol_m2_d": ch4[c], "k_x_B_mmol_per_m2": cap[c],
                         "mu_required_h": mu[c], "doubling_time_h": math.log(2) / mu[c]})
        ranks[name] = exact(cap[CORES], ch4[CORES])
    t = pd.DataFrame(rows)
    t.to_csv(HERE / "specific_rates.tsv", sep="\t", index=False)
    ky = {"MePn ratio 0.07 (median)": float(YM.loc[mep, "k_CH4_r0.07"].median()),
          "MePn ratio 0.33 (median)": float(YM.loc[mep, "k_CH4_r0.33"].median()),
          "MePn only P source (median)": float(YM.loc[mep, "k_CH4_MePn_only_P"].median()),
          "Betaine TMA per g (range)": [float(YB.loc[bet, "k_TMA_unbounded_min"].min()), float(YB.loc[bet, "k_TMA_unbounded_min"].max())],
          **{f"Methanogen {s}": v for s, v in km.items()}, "methanogen_k_source": ksrc}
    summ = {"yields_mmol_CH4_per_g": ky, "ranking_spearman_exact": ranks,
            "doubling_time_h_range": t.groupby("scenario").doubling_time_h.agg(["min", "max"]).round(3).to_dict("index"),
            "mu_required_range": t.groupby("scenario").mu_required_h.agg(["min", "max"]).round(4).to_dict("index")}
    json.dump(summ, open(HERE / "specific_rates.json", "w"), indent=1)
    pd.set_option("display.width", 200)
    print(json.dumps(ky, indent=0))
    print(t.pivot(index="scenario", columns="core", values="doubling_time_h").round(2).to_string())
    print("ranking:", {k: (round(v[0], 2), round(v[1], 2)) for k, v in ranks.items()})


if __name__ == "__main__":
    main()
