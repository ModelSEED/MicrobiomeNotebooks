#!/usr/bin/env python
"""Methane (or TMA) produced per gram of biomass formed, k, for every 2026 guild model, at a fixed biomass flux of 1e-3
h-1 with the community protocol's objective (least lactate uptake):
  MePn consumers: ratio 0.07 and 0.33 (from reference_k.json), and MePn as the only phosphorus source (Pi exchange
    closed; the upper limit, k = biomass P content);
  betaine reducers: TMA per g with the 0.1 betaine-uptake bound of the protocol, and with betaine uptake unbounded
    (how much betaine reduction the model's energy budget favours); methane = 0.06 x TMA (the manuscript's community
    conversion) or 2.25 x TMA (full methylotrophic conversion, 4 TMA -> 9 CH4).
Run: ~/Documents/py_venv/bin/python per_gram_yields.py -> per_gram_yields.tsv
"""
import json
from pathlib import Path
import cobra
import pandas as pd
import simulate_guilds as sg

HERE = Path(__file__).resolve().parent
M = HERE.parent
B = 1e-3
ref = json.load(open(HERE / "reference_k.json"))["mepn"]
rows = []
for f in sorted((M / "models/mepn_consumer").glob("*.json")):
    m = cobra.io.load_json_model(str(f))
    med = dict(sg.MEDIA); med[sg.PI_EX] = 0.0                     # MePn the only P source
    with m:
        sg.set_medium_strict(m, med)
        m.reactions.bio1.bounds = (B, B)
        m.objective = m.reactions.get_by_id(sg.LAC_EX); m.reactions.get_by_id(sg.LAC_EX).bounds = (-1000, 0)
        s = m.optimize()
        kP = s.fluxes[sg.CH4_EX] / B if s.status == "optimal" else float("nan")
    rows.append({"guild": "mepn_consumer", "mag_id": f.stem, "k_CH4_r0.07": ref[f.stem]["0.07"]["methane"] / B,
                 "k_CH4_r0.33": ref[f.stem]["0.33"]["methane"] / B, "k_CH4_MePn_only_P": kP})
for f in sorted((M / "models/betaine_reducer").glob("*.json")):
    m = cobra.io.load_json_model(str(f))
    rec = {"guild": "betaine_reducer", "mag_id": f.stem}
    for lab, cap in (("bound0.1", 0.1), ("unbounded", 1000.0)):
        med = {k: v for k, v in sg.MEDIA.items() if k != sg.MEPN_EX}; med[sg.BET_EX] = cap
        with m:
            sg.set_medium_strict(m, med)
            m.reactions.bio1.bounds = (B, B)
            m.objective = m.reactions.get_by_id(sg.LAC_EX); m.reactions.get_by_id(sg.LAC_EX).bounds = (-1000, 0)
            s = m.optimize()
            if s.status == "optimal":
                lac = s.objective_value
                fva = cobra.flux_analysis.flux_variability_analysis(m, [sg.TMA_EX], fraction_of_optimum=1.0)
                rec[f"k_TMA_{lab}_min"] = fva.loc[sg.TMA_EX, "minimum"] / B
                rec[f"k_TMA_{lab}_max"] = fva.loc[sg.TMA_EX, "maximum"] / B
    rows.append(rec)
t = pd.DataFrame(rows)
t.to_csv(HERE / "per_gram_yields.tsv", sep="\t", index=False)
pd.set_option("display.width", 220)
print(t.groupby("guild").median(numeric_only=True).T.round(4))
print(t[t.guild == "betaine_reducer"][["mag_id", "k_TMA_bound0.1_min", "k_TMA_bound0.1_max", "k_TMA_unbounded_min", "k_TMA_unbounded_max"]].round(2).to_string(index=False))
