#!/usr/bin/env python
"""Synthesis of the hypothesized methane sources of the six unrestored cores (R1_A/B/C, R2_A/B/C; 12 samples D1, D2):
for each pathway, its read- and MAG-based abundance, its most abundant organisms, the turnover each would need to supply
the measured CH4 on its own, and an attribution of the measured CH4 among the pathways.

Pathways (the manuscript's hypotheses, Conclusions): methylphosphonate (C-P lyase); betaine -> TMA (betaine reductase)
and choline / carnitine -> TMA (CutC, CntA), which feed methylotrophic methanogenesis; methylamines -> CH4; methanol
-> CH4; hydrogenotrophic (H2/CO2, formate) and acetoclastic methanogenesis.

Abundance (per 100 genome equivalents, GE; ../../read_markers/):
  reads: marker-gene copies per GE in all reads: phnJ; grdI; cutC + cntA; mcrA x the pooled share of the methanogen
         group among the mcrA reads (12 samples; ANME-1-related Methanospirare and unassigned reads excluded from the
         groups); methanol: mtaB x its pooled archaeal share (10d_methanol_taxonomy.py: most mtaB reads are bacterial).
  MAGs:  the guild MAGs' genomes per GE (key_unrestored_cores.tsv); 0 where no MAG carries the pathway.
Top organisms: reads, genus of the best hit (read_taxonomy.tsv, methanol_taxonomy.tsv, mcrA_taxonomy_by_sample.tsv);
  MAGs, the guild species with the highest mean relative abundance in the 12 samples.
Turnover needed for the pathway alone to supply the measured CH4 of each core:
  substrate pool: pool CH4 potential / measured flux (d), from porewater NMR (D1 0-5 cm, D2 5-15 cm; porosity 0.736 as
     community_eval/methanogen_test.py): methylamines (CH4 per TMA 2.25, DMA 1.5, MMA 0.75), betaine via TMA (2.25),
     methanol (0.75), formate (0.25; H2 not measured), acetate (1). MPn, choline, carnitine and H2 were not measured.
  organism doubling time (h): per-gram yields and measured biomass (../simulations/specific_rates.py) for the MPn
     consumers, the betaine reducers (TMA for complete conversion, 2.25 CH4 per TMA) and the methylotrophic
     methanogens (curated Methanolobus model, k = 126 mmol CH4 g-1). Read-based methanogen biomass = the Methanolobus
     MAG biomass x the mcrA read/MAG fold (all mcrA), x the group's share. Hydrogenotrophic and acetoclastic
     methanogens have no model; their doubling time at k = 126 is a lower bound (a lower growth yield, as expected for
     these pathways, means more CH4 per gram formed and allows slower growth).
Attribution: at equal specific growth rates, each guild's CH4 is proportional to k x B (read-based); methanogen CH4
  is split among the groups by their mcrA share. Sensitivities: MPn as the only P source of the MPn consumers; MPn
  consumers growing 10 times faster than the methanogens. Coverage: the fraction of the measured CH4 that growth-coupled
  production would supply at a doubling time of one day.
Run: ~/Documents/py_venv/bin/python pathway_summary.py -> pathway_summary.json, pathway_abundance.tsv,
     pathway_turnover.tsv, attribution.tsv
"""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
M = HERE.parent
S = HERE.parents[2]
RMK = S / "reassembly_2026" / json.load(open(S / "reassembly_2026/read_survey_current.json"))["results"]
CORES = ["R1_A", "R1_B", "R1_C", "R2_A", "R2_B", "R2_C"]
SAMPLES = [f"{c}_{d}" for c in CORES for d in ("D1", "D2")]
LAYER = {"D1": 5, "D2": 10}
POROSITY = 1 - 0.7 / 2.65
LN2 = math.log(2)

meta = pd.read_csv(S / "data/Cliff_Sample_Metadata_BGC_NMR.csv").set_index("Sample")
ch4 = pd.Series({c: meta.loc[[c + "_D1", c + "_D2"], "CH4_umol_m2_d"].mean() for c in CORES})   # umol m-2 d-1
ko = pd.read_csv(RMK / "per_sample_KO.tsv", sep="\t")
ko = ko[ko["sample"].isin(SAMPLES)].pivot(index="sample", columns="label", values="copies_per_GE")
key = pd.read_csv(RMK / "key_unrestored_cores.tsv", sep="\t").pivot(index="sample", columns="guild",
                                                                       values="mag_genomes_per_GE")
rtax = pd.read_csv(RMK / "read_taxonomy.tsv", sep="\t")
rtax = rtax[rtax["sample"].isin(SAMPLES)]
mtax = pd.read_csv(RMK / "mcrA_taxonomy_by_sample.tsv", sep="\t")
mtax = mtax[mtax["sample"].isin(SAMPLES)]
otax = pd.read_csv(RMK / "methanol_taxonomy.tsv", sep="\t")
otax = otax[otax["sample"].isin(SAMPLES)]
sr = pd.read_csv(M / "simulations/specific_rates.tsv", sep="\t")
km = json.load(open(M / "models/methanogen/curation_k.json"))["methylamines"]
Y = pd.read_csv(M / "simulations/per_gram_yields.tsv", sep="\t")
g = pd.read_csv(S / "correlational/final585/guilds_final585.csv").set_index("mag")
ids = json.load(open(S / "correlational/final585/iterativeIDs_final585.json"))

# ---- methanogen groups of the mcrA reads (pooled over the 12 samples) ----------------------------------------------------
METHYLO = {"Methanolobus", "Methanohalophilus", "Methanococcoides", "Methanohalobium", "Methanomethylovorans",
           "Methanosalsum", "Methanosarcina",
           # Methanomassiliicoccales: H2-dependent methylotrophs
           "Methanomethylophilus", "Methanomassiliicoccus", "Methanogranum", "Candidatus Methanoplasma",
           "Methanomassiliicoccales", "Thermoplasmatales", "Methanogenic"}
HYDRO = {"Methanogenium", "Methanoplanus", "Methanoculleus", "Methanofollis", "Methanobacterium", "Methanomicrobiales",
         "Methanocalculus", "Methanospirillum", "Methanobrevibacter", "Methanococcus"}
ACETO = {"Methanothrix", "Methanocrinis"}
NONSOURCE = {"Methanospirare"}                                   # ANME-1-related (Methanophagales): alkyl/methane oxidation


def mgroup(genus):
    if genus in METHYLO:
        return "methylotrophic"
    if genus in HYDRO:
        return "hydrogenotrophic"
    if genus in ACETO:
        return "acetoclastic"
    if genus in NONSOURCE:
        return "anme"
    return "unassigned"


mtax = mtax.assign(group=mtax.genus.fillna("").map(mgroup))
mg = mtax.groupby("group").reads.sum()
assert set(mg.index) <= {"methylotrophic", "hydrogenotrophic", "acetoclastic", "anme", "unassigned"}, mg
n_meth = int(mg.drop(["anme", "unassigned"], errors="ignore").sum())
share = {k: float(mg.get(k, 0) / n_meth) for k in ("methylotrophic", "hydrogenotrophic", "acetoclastic")}
methanogen_frac_of_mcrA = n_meth / mg.sum() + mg.get("unassigned", 0) / mg.sum()   # all but the ANME-related reads
arch_mtaB = otax[otax.label == "mtaB"]
mtaB_arch_share = float(arch_mtaB[arch_mtaB.domain == "Archaea"].reads.sum() / arch_mtaB.reads.sum())


def top_genera(df, n=3, min_share=0.04):
    x = df.groupby("genus").reads.sum().sort_values(ascending=False)
    x = x / x.sum()
    return [(k, round(float(v), 3)) for k, v in x.head(n).items() if v >= min_share]


def top_mags(guild, n=2):
    x = g[g[guild] == True].sort_values("mean_rel_abund_12cores", ascending=False)
    return [ids.get(m, m) for m in x.index[:n]]


# ---- pathways -----------------------------------------------------------------------------------------------------------
P = [
    dict(key="mpn", label="Methylphosphonate", enzyme="C–P lyase (PhnJ)", product="CH4 + Pi", cls="mpn",
         reads=ko["phnJ"], mags=key["mepn_consumer"], top_reads=top_genera(rtax[rtax.label == "phnJ"]),
         top_mags=top_mags("mepn_consumer"), pool=None, pool_note="detected (GC-MS), not quantified"),
    dict(key="betaine", label="Glycine betaine", enzyme="betaine reductase (GrdI)", product="TMA", cls="tma",
         reads=ko["grdI"], mags=key["betaine_reducer"], top_reads=top_genera(rtax[rtax.label == "grdI"]),
         top_mags=top_mags("betaine_reducer"), pool={"Betaine": 2.25}),
    dict(key="choline", label="Choline, carnitine", enzyme="CutC, CntA", product="TMA", cls="tma",
         reads=ko["cutC"] + ko["cntA"], mags=key["tma_candidate_one_sided"],
         top_reads=top_genera(rtax[rtax.label.isin(["cutC", "cntA"])]), top_mags=top_mags("tma_candidate_one_sided"),
         pool=None, pool_note="not measured"),
    dict(key="methylamines", label="TMA, DMA, MMA", enzyme="MttB/MtbB/MtmB + Mcr", product="CH4", cls="meth",
         reads=ko["mcrA"] * share["methylotrophic"], mags=key["methanogen"],
         top_reads=top_genera(mtax[mtax.group == "methylotrophic"]), top_mags=top_mags("methanogen", 1),
         pool={"Trimethylamine": 2.25, "Dimethylamine": 1.5, "Methylamine": 0.75}),
    dict(key="methanol", label="Methanol", enzyme="MtaB + Mcr", product="CH4", cls="meth",
         reads=ko["mtaB"] * mtaB_arch_share, mags=pd.Series(0.0, index=SAMPLES),
         top_reads=top_genera(arch_mtaB[arch_mtaB.domain == "Archaea"]), top_mags=[], pool={"Methanol": 0.75}),
    dict(key="hydrogenotrophic", label="H2/CO2, formate", enzyme="hydrogenotrophic Mcr", product="CH4", cls="meth",
         reads=ko["mcrA"] * share["hydrogenotrophic"], mags=pd.Series(0.0, index=SAMPLES),
         top_reads=top_genera(mtax[mtax.group == "hydrogenotrophic"]), top_mags=[], pool={"Formate": 0.25},
         pool_note="formate only; H2 not measured"),
    dict(key="acetoclastic", label="Acetate", enzyme="acetoclastic Mcr", product="CH4", cls="meth",
         reads=ko["mcrA"] * share["acetoclastic"], mags=pd.Series(0.0, index=SAMPLES),
         top_reads=top_genera(mtax[mtax.group == "acetoclastic"], min_share=0), top_mags=[], pool={"Acetate": 1.0}),
]


def pool_ch4(sub, c):
    return sum(float(meta.loc[f"{c}_{d}", s]) * y * POROSITY * cm * 10 for s, y in sub.items() for d, cm in LAYER.items())


def dbl(name):
    return sr[sr.scenario == name].set_index("core").reindex(CORES)


B = {name: dbl(name).biomass_g_m2 for name in ("MePn, ratio 0.07", "Methanogen, methylamines",
                                                "Methanogen, methylamines (read-based abundance)")}
Bmeth_read = B["Methanogen, methylamines (read-based abundance)"] * methanogen_frac_of_mcrA   # without ANME-related
Bgroup = {k: Bmeth_read * v for k, v in share.items()}
ch4_mmol_h = ch4 / 1000 / 24


def doubling_from(cap):                                          # cap: mmol CH4 m-2 per unit mu
    return LN2 / (ch4_mmol_h / cap)


turn = []
for p in P:
    row = {"pathway": p["key"]}
    for c in CORES:
        row[f"pool_turnover_d_{c}"] = pool_ch4(p["pool"], c) / ch4[c] if p["pool"] else np.nan
    turn.append(row)
dt = {
    "mpn": {"MAG": dbl("MePn, ratio 0.07").doubling_time_h, "reads": dbl("MePn, ratio 0.07 (read-based abundance)").doubling_time_h,
            "MAG_onlyP": dbl("MePn, only P source").doubling_time_h,
            "reads_onlyP": dbl("MePn, only P source (read-based abundance)").doubling_time_h},
    "betaine": {"MAG": dbl("Betaine, 2.25 CH4 per TMA").doubling_time_h,
                "reads": dbl("Betaine, 2.25 CH4 per TMA (read-based abundance)").doubling_time_h},
    "methylamines": {"MAG": dbl("Methanogen, methylamines").doubling_time_h,
                     "reads": doubling_from(km * Bgroup["methylotrophic"])},
    "hydrogenotrophic": {"reads_lower_bound": doubling_from(km * Bgroup["hydrogenotrophic"])},
    "acetoclastic": {"reads_lower_bound": doubling_from(km * Bgroup["acetoclastic"])},
}
for row in turn:
    for kind, ser in dt.get(row["pathway"], {}).items():
        for c in CORES:
            row[f"doubling_h_{kind}_{c}"] = float(ser[c])
T = pd.DataFrame(turn).set_index("pathway")
T.to_csv(HERE / "pathway_turnover.tsv", sep="\t", float_format="%.5g")

A = pd.DataFrame({p["key"]: (p["reads"] * 100) for p in P}).add_suffix("_reads").join(
    pd.DataFrame({p["key"]: (p["mags"] * 100) for p in P}).add_suffix("_MAGs"))
A.index.name = "sample"
A.to_csv(HERE / "pathway_abundance.tsv", sep="\t", float_format="%.5g")

# ---- attribution at equal specific growth rates (read-based biomass) -------------------------------------------------
Bmep = dbl("MePn, ratio 0.07 (read-based abundance)").biomass_g_m2
k07, k33, kP = (Y[Y.guild == "mepn_consumer"][c].median() for c in ("k_CH4_r0.07", "k_CH4_r0.33", "k_CH4_MePn_only_P"))
# capacity of the MPn consumers from the per-model k (as in specific_rates): cap = ln2 / doubling x CH4 rate
cap = {"mpn_0.07": ch4_mmol_h * dbl("MePn, ratio 0.07 (read-based abundance)").doubling_time_h / LN2,
       "mpn_0.33": ch4_mmol_h * dbl("MePn, ratio 0.33 (read-based abundance)").doubling_time_h / LN2,
       "mpn_onlyP": ch4_mmol_h * dbl("MePn, only P source (read-based abundance)").doubling_time_h / LN2,
       "methanogens": km * Bmeth_read}
att = []
for scen, mult in (("mpn_0.07", 1), ("mpn_0.33", 1), ("mpn_onlyP", 1), ("mpn_0.07", 10), ("mpn_onlyP", 10)):
    for c in CORES:
        m_ = cap[scen][c] * mult
        tot = m_ + cap["methanogens"][c]
        att.append({"scenario": f"{scen} x{mult}", "core": c, "MPn": m_ / tot,
                    **{k: cap["methanogens"][c] / tot * v for k, v in share.items()}})
AT = pd.DataFrame(att)
AT.to_csv(HERE / "attribution.tsv", sep="\t", index=False, float_format="%.5g")
# growth-coupled coverage at a doubling time of one day (all methanogens, read-based; MPn at 0.07 and only P)
cover = {"methanogens": (cap["methanogens"] * LN2 / 24 / ch4_mmol_h),
         "mpn_0.07": cap["mpn_0.07"] * LN2 / 24 / ch4_mmol_h, "mpn_onlyP": cap["mpn_onlyP"] * LN2 / 24 / ch4_mmol_h}
q_meth = ch4_mmol_h / Bmeth_read                                  # mmol CH4 g-1 h-1, all methanogens (read-based)

R1, R2 = CORES[:3], CORES[3:]


def rng(s, cores=CORES):
    return [float(s[cores].min()), float(s[cores].max())]


out = {
    "pathways": [{k: (v if not isinstance(v, pd.Series) else None) for k, v in p.items() if k not in ("reads", "mags")}
                 for p in P],
    "abundance_per_100GE": {p["key"]: {"reads": rng(p["reads"] * 100, SAMPLES), "MAGs": rng(p["mags"] * 100, SAMPLES),
                                       "reads_median": float((p["reads"] * 100).median()),
                                       "MAGs_median": float((p["mags"] * 100).median())} for p in P},
    "mcrA_groups": {"reads": {k: int(v) for k, v in mg.items()}, "share_of_methanogen_reads": share,
                    "mtaB_archaeal_share": mtaB_arch_share},
    "pool_turnover_d": {r: rng(T.loc[r, [f"pool_turnover_d_{c}" for c in CORES]].set_axis(CORES))
                        for r in T.index if T.loc[r, [f"pool_turnover_d_{c}" for c in CORES]].notna().all()},
    "doubling_h": {p_: {kind: {"R1": rng(s, R1), "R2": rng(s, R2), "all": rng(s)} for kind, s in d.items()}
                   for p_, d in dt.items()},
    "attribution_equal_growth": {s: {col: rng(AT[AT.scenario == s].set_index("core")[col])
                                     for col in ("MPn", "methylotrophic", "hydrogenotrophic", "acetoclastic")}
                                 for s in AT.scenario.unique()},
    "coverage_at_1d_doubling": {k: {"R1": rng(v, R1), "R2": rng(v, R2)} for k, v in cover.items()},
    "q_required_methanogens_read": {"R1": rng(q_meth, R1), "R2": rng(q_meth, R2)},
    "k_mmol_per_g": {"mpn_0.07": float(k07), "mpn_0.33": float(k33), "mpn_onlyP": float(kP), "methanogen": km},
}
json.dump(out, open(HERE / "pathway_summary.json", "w"), indent=1, default=float)
pd.set_option("display.width", 220)
print(json.dumps({k: out[k] for k in ("abundance_per_100GE", "mcrA_groups", "pool_turnover_d")}, indent=0, default=float))
print(json.dumps({k: out[k] for k in ("doubling_h", "attribution_equal_growth", "coverage_at_1d_doubling",
                                      "q_required_methanogens_read")}, indent=0, default=float))
for p in P:
    print(p["key"], p["top_reads"], p["top_mags"])
