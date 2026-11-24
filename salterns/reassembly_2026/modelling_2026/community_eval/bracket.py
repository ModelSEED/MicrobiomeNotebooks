#!/usr/bin/env python
"""Bracket for methanogen_test.py: simulated areal CH4 / measured CH4 (total over the six unrestored cores) for each
guild, with the rate made explicit: CH4 = k x biomass x f x mu (x 24 h d-1), where
  k   CH4 per g biomass formed: MePn consumers at a phosphate-source ratio r of 0.07, 0.33, or 'all P from MePn'
      (k = biomass P content = k_0.07 x 1.07 / 0.07 per model); methanogen (Methanolobus.1) from the methanogenic
      KOfam-arm builds on methanol (arms B-D; arm A, the 2026 build, makes no CH4);
  f   biomass factor: 1 = the pipeline's constants (DNA = 10 % of dry biomass, complete extraction); 10 = e.g. DNA at
      about 3 % of dry biomass with about 30 % extraction efficiency;
  mu  specific growth rate, 0.001 / 0.01 / 0.05 h-1 (sediment doubling times of about four weeks / three days / 14 h).
The manuscript's conversion corresponds to f = 1, mu = 1 h-1. Growth-coupled CH4 only; methanogen maintenance
methanogenesis would add to the methanogen column. Run after methanogen_test.py -> bracket.tsv
"""
import json
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
t = pd.read_csv(HERE / "methanogen_test_by_core.tsv", sep="\t", index_col=0)
mt = json.load(open(HERE / "methanogen_test.json"))
ref = json.load(open(HERE.parent / "simulations/reference_k.json"))
meas = t.CH4_measured_umol_m2_d.sum()
# MePn capacity per unit mu (umol m-2 d-1 per h-1) from methanogen_test (meas / mu_required)
cap = {"0.07": (t.CH4_measured_umol_m2_d / t["mepn_mu_required_h_r0.07"]).sum(),
       "0.33": (t.CH4_measured_umol_m2_d / t["mepn_mu_required_h_r0.33"]).sum()}
cap["all P from MePn"] = cap["0.07"] * 1.07 / 0.07
kmeth = {a: r["k_ch4_per_g"] for a, r in mt["methanogen_models"].items() if r.get("k_ch4_per_g")}
Bm_total = (t.methanogen_biomass_g_m2).sum()
rows = []
for f in (1, 10):
    for mu in (0.001, 0.01, 0.05, 1.0):
        row = {"biomass_factor": f, "mu_h": mu}
        for r, c in cap.items():
            row[f"MePn r={r}"] = c * f * mu / meas
        for a, k in kmeth.items():
            row[f"methanogen arm {a}"] = k * Bm_total * f * mu * 1000 * 24 / meas
        rows.append(row)
b = pd.DataFrame(rows)
b.to_csv(HERE / "bracket.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print("simulated / measured CH4 (total of six cores)")
print(b.round(4).to_string(index=False))
