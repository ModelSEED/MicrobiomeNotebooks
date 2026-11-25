#!/usr/bin/env python
"""Per-gram framing of the modelling text (2026-11-24): writes r232_review_edits_2026-11-24.json from
r232_review_edits_2026-11-23.json (run make_modelling_edits_2026-11-23.py first). The 11-23 modelling text compared
'simulated' with measured areal emissions through a conversion that multiplied CH4 per g biomass formed by standing
biomass per hour, i.e. it assumed every member grows at 1 h-1. Here the model output is the methane produced per gram of
biomass formed, and the comparison is the growth (or, for the methanogen, the specific methanogenic activity) each
guild's measured biomass would need to supply the measured emissions. Sources (reassembly_2026/modelling_2026/):
  simulations/per_gram_yields.tsv, specific_rates.tsv/.json (per_gram_yields.py, specific_rates.py; PrmA blocked)
  models/methanogen/curation_k.json, CURATION.md (curated Methanolobus model)
  community_eval/methanogen_test_by_core.tsv, methanogen_test_curated.tsv (specific activity, methylamine pools)
  ../read_markers/results/ (read-based guild abundance: key_unrestored_cores.tsv, mcrA_taxonomy_unrestored_pooled.tsv,
  binned_unbinned_split.tsv); specific_rates.py's "(read-based abundance)" scenarios scale each guild's biomass by it
The wording's preconditions are asserted. Run: ~/Documents/py_venv/bin/python make_modelling_edits_2026-11-24.py
"""
import copy
import json
import math
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
S = HERE.parent
M = S / "reassembly_2026" / "modelling_2026"
SIM = M / "simulations"
CE = M / "community_eval"
DATE = "2026-11-24T00:00:00Z"
CORES = ["R1_A", "R1_B", "R1_C", "R2_A", "R2_B", "R2_C"]
SRC = ("reassembly_2026/modelling_2026/simulations/{per_gram_yields,specific_rates}.tsv (per_gram_yields.py, "
       "specific_rates.py); models/methanogen/CURATION.md; community_eval/methanogen_test_curated.tsv; "
       "reassembly_2026/read_markers/README.md (read-based guild abundance)")

base = json.load(open(HERE / "r232_review_edits_2026-11-23.json"))
Y = pd.read_csv(SIM / "per_gram_yields.tsv", sep="\t")
YM, YB = Y[Y.guild == "mepn_consumer"], Y[Y.guild == "betaine_reducer"]
SR = pd.read_csv(SIM / "specific_rates.tsv", sep="\t")
SRJ = json.load(open(SIM / "specific_rates.json"))
KM = json.load(open(M / "models/methanogen/curation_k.json"))
MT = pd.read_csv(CE / "methanogen_test_curated.tsv", sep="\t", index_col=0)
MS = pd.read_csv(M / "models/models_summary.tsv", sep="\t")
IDS = json.load(open(S / "correlational/final585/iterativeIDs_final585.json"))
G = pd.read_csv(S / "correlational/final585/guilds_final585.csv").set_index("mag")
rel = pd.read_csv(S / "reassembly_2026/abundance_final/final585_relabund_24sediment.tsv", sep="\t", index_col=0)
dna = pd.read_csv(S / "data/SaltPondsDNA.csv").set_index("Sample")["DNA conc (µg/g soil)"]
RMK = S / "reassembly_2026/read_markers/results"
RM = pd.read_csv(RMK / "key_unrestored_cores.tsv", sep="\t")
MTAX = pd.read_csv(RMK / "mcrA_taxonomy_unrestored_pooled.tsv", sep="\t")
SPLIT = pd.read_csv(RMK / "binned_unbinned_split.tsv", sep="\t")


def dt(name):
    d = SR[SR.scenario == name].set_index("core").doubling_time_h
    return d


def fmt_t(h):
    """hours -> compact text (min below 1 h, h below 2 d, d above)."""
    if h < 0.95:
        return f"{h * 60:.0f} min"
    if h < 48:
        return f"{h:.1f} h" if h < 10 else f"{h:.0f} h"
    return f"{h / 24:.0f} d"


def rng_t(name, cores=CORES):
    d = dt(name)[cores]
    return f"{fmt_t(d.min())} to {fmt_t(d.max())}"


def num(v, spec):
    return format(v, spec).replace("-", "−")


# ---- facts and the wording's preconditions ------------------------------------------------------------------------------
k07, k33, kP = YM["k_CH4_r0.07"].median(), YM["k_CH4_r0.33"].median(), YM["k_CH4_MePn_only_P"].median()
kpos07 = YM["k_CH4_r0.07"][YM["k_CH4_r0.07"] > 1e-9]
km = KM["methylamines"]
assert abs(km - SRJ["yields_mmol_CH4_per_g"]["Methanogen methylamines"]) < 1e-9, "methanogen k differs from curation_k.json"
assert (YM["k_CH4_r0.07"] > 1e-9).all(), "an MPn model releases no methane (PrmA block missing?)"
fold = round(km / k07, -2)                                   # "about 2,900"
tma = YB["k_TMA_unbounded_min"]
d07, d33, dP = dt("MePn, ratio 0.07"), dt("MePn, ratio 0.33"), dt("MePn, only P source")
db06, dch = dt("Betaine, 0.06 CH4 per TMA"), dt("Betaine, 2.25 CH4 per TMA")
dm = dt("Methanogen, methylamines")
assert d07.max() < 1 and d33.max() < 1.5 and dP.max() < 24, "MPn doubling times no longer minutes-to-hours"
q = MT["methanogen_q_required_mmol_g_h"]
turn = MT["methylamine_turnover_days_required"]
R1, R2 = ["R1_A", "R1_B", "R1_C"], ["R2_A", "R2_B", "R2_C"]
assert dm[R2].max() < 2 and dm[R1].min() > 3, "methanogen R1/R2 contrast changed"
RB = " (read-based abundance)"
rho_all = SRJ["ranking_spearman_exact"]
rho = {k: v for k, v in rho_all.items() if not k.endswith(RB)}              # MAG-based abundance (primary)
rho_rb = {k[:-len(RB)]: v for k, v in rho_all.items() if k.endswith(RB)}     # read-based abundance (sensitivity)
assert abs(rho["MePn, ratio 0.07"][0] - rho["MePn, only P source"][0]) < 1e-9 and abs(rho["MePn, ratio 0.07"][0]) < 0.3
assert all(rho[k][1] > 0.2 for k in rho), "a guild now ranks the cores (MAG-based)"
assert all(v[1] > 0.15 and v[0] < 0.3 for v in rho_rb.values()), "a guild now ranks the cores (read-based)"
# ---- read-based guild abundance (../read_markers): per sample, marker copies per genome equivalent (GE) in all reads /
# the guild MAGs' genomes per GE; the methanogen doubling times assume every mcrA-carrying organism has Methanolobus's k
FOLD = {gname: RM[RM.guild == gname].set_index("sample").fold_read_vs_MAG.astype(float)
        for gname in ("methanogen", "mepn_consumer", "betaine_reducer")}
rmm = RM[RM.guild == "methanogen"]
mcrA_pooled = float(rmm.mcrA_copies_per_GE.sum() / rmm.mag_genomes_per_GE.sum())
big = rmm[rmm.fold_read_vs_MAG > 3]
assert set(big["sample"]) == {"R2_A_D1", "R2_A_D2", "R2_C_D2"}, "the samples with a several-fold mcrA excess changed"
mcrA_big_reads = [int(big.mcrA_reads.min()), int(big.mcrA_reads.max())]
other = MTAX[MTAX.genus != "Methanolobus"]
other_share = float(other.reads.sum() / MTAX.reads.sum())
other_unasm = float(other.unassembled.sum() / other.reads.sum())
top_other = other.groupby("genus").reads.sum().sort_values(ascending=False)
assert list(top_other.index[:3]) == ["Methanohalophilus", "Methanogenium", "Methanococcoides"], top_other.head()
dm_rb, d07_rb, dP_rb = dt("Methanogen, methylamines" + RB), dt("MePn, ratio 0.07" + RB), dt("MePn, only P source" + RB)
d33_rb, db06_rb = dt("MePn, ratio 0.33" + RB), dt("Betaine, 0.06 CH4 per TMA" + RB)
q_rb = SR[SR.scenario == "Methanogen, methylamines" + RB].set_index("core").mu_required_h * km
assert dm_rb[R2].max() < 6 and dm_rb[R2].min() > 1.5 and q_rb[R2].min() > 10, "read-based R2 methanogen wording"
assert (dm_rb >= dm - 1e-9).all() and (d07_rb >= d07 - 1e-9).all(), "read-based abundance should only lengthen doubling"
fmep, fbet = FOLD["mepn_consumer"], FOLD["betaine_reducer"]
mk = SPLIT[SPLIT.label.isin(["uL11", "uL16", "uL2", "uL4", "uS11", "uS12", "uS3", "uS9"])
           & SPLIT["sample"].str.contains("_D[12]$")].copy()
mk["site"] = mk["sample"].str.split("_").str[0]
mks = mk.groupby("site")[["reads", "unassembled"]].sum()
unasm = (mks.unassembled / mks.reads * 100).round(0)
assert unasm[["R2A", "SF2"]].min() > 70 and unasm[["R1", "R2"]].max() < 50, "unassembled marker shares changed"
DNA3 = 0.10 / 0.03                                            # biomass multiplier if DNA were 3 % of dry biomass
n_mep = int(((MS.guild == "mepn_consumer") & (MS.kept == True)).sum())
n_bet = int(((MS.guild == "betaine_reducer") & (MS.kept == True)).sum())
drop_id = IDS[MS[(MS.guild == "mepn_consumer") & (MS.kept != True)].mag_id.iloc[0]]
# shares of the MPn consumers' methane per unit of growth (k x measured biomass, summed over the six cores)


def bio(m):
    return sum(rel.loc[m, f"{c}_{d}"] * dna[f"{c}_{d}"] * cm for c in CORES for d, cm in (("D1", 5), ("D2", 10)))


contrib = pd.Series({m: YM.set_index("mag_id").loc[m, "k_CH4_r0.07"] * bio(m) for m in YM.mag_id})
fam = contrib.groupby(lambda m: G.loc[m, "family"]).sum() / contrib.sum()
gen = contrib.groupby(lambda m: G.loc[m, "genus"]).sum() / contrib.sum()
FACTS = dict(k07=k07, k33=k33, kP=kP, k_meth=km, k_meth_range=KM["coupling_range"], fold=fold,
             tma_range=[float(tma.min()), float(tma.max())], mepn_doubling_h={"0.07": [d07.min(), d07.max()],
             "0.33": [d33.min(), d33.max()], "only_P": [dP.min(), dP.max()]}, betaine_doubling_h_006=[db06.min(), db06.max()],
             betaine_chain_doubling_h=[dch.min(), dch.max()], methanogen_doubling_h=dm.round(2).to_dict(),
             methanogen_q_required=q.round(2).to_dict(), methylamine_turnover_d=turn.round(2).to_dict(),
             rho={k: v[0] for k, v in rho.items()}, rhodobacteraceae_share=float(fam.get("Rhodobacteraceae", 0)),
             roseovarius_share=float(gen.get("Roseovarius", 0)),
             read_based=dict(rho={k: v for k, v in rho_rb.items()}, mcrA_fold_pooled=mcrA_pooled,
                             fold_range={g_: [float(v.min()), float(v.max())] for g_, v in FOLD.items()},
                             mcrA_other_methanogens_share=other_share, mcrA_other_unassembled=other_unasm,
                             methanogen_q_required=q_rb.round(2).to_dict(), methanogen_doubling_h=dm_rb.round(2).to_dict(),
                             mepn_doubling_h={"0.07": [d07_rb.min(), d07_rb.max()], "0.33": [d33_rb.min(), d33_rb.max()],
                                              "only_P": [dP_rb.min(), dP_rb.max()]},
                             betaine_doubling_h_006=[db06_rb.min(), db06_rb.max()],
                             marker_reads_unassembled_pct_by_site=unasm.to_dict()))

NEW = {}


def put(key, paragraph, kind, old, new, rationale, mirror_index=None):
    e = {"paragraph": paragraph, "kind": kind, "old": old, "new": new, "rationale": rationale, "sources": SRC,
         "origin": "modelling_2026_per_gram", "date": DATE}
    if mirror_index is not None:
        e["mirror_index"] = mirror_index
    NEW[key] = e


put("abstract", 18, "replace",
    "This was corroborated by metabolic modeling, where trace consumption of methylphosphonic acid was sufficient to reproduce the experimental methane emissions.",
    f"In metabolic models, methylphosphonate-consuming bacteria released {k07:.2f}–{kP:.2f} mmol of methane per gram of "
    f"biomass formed, set by the phosphorus they assimilate, whereas the methylotrophic methanogen detected in the same "
    f"sediments released {km:.0f} mmol; to supply the measured emissions, the methylphosphonate consumers would have to "
    f"double every {fmt_t(d07.min())} to {fmt_t(dP.max())}, and the methanogen would need a specific methanogenic activity "
    f"of {q.min():.0f}–{q.max():.0f} mmol g⁻¹ h⁻¹.",
    "per-gram framing; the Abstract states yields and required rates, the inference is in the P18 and P100 comments",
    mirror_index=19)
put("p18_comment", 18, "comment", "a primary methanogenic substrate",
    f"Review (2026-11-24): per gram of biomass formed, methylphosphonate assimilation releases about {fold:,.0f}-fold less "
    "methane than methylotrophic methanogenesis (Figure 7), and the measured emissions would need the MPn consumers to "
    "double every minutes to hours. The models therefore do not support methylphosphonate as the primary substrate "
    "unless MPn is degraded well beyond assimilatory phosphorus demand, which they do not represent. Consider "
    "'a' instead of 'a primary' unless the in vitro data alone carry 'primary'.", "Abstract claim", mirror_index=19)
p35 = next(e for e in base["edits"] if e["paragraph"] == 35 and e["kind"] == "replace_paragraph")
put("p35", 35, "replace_paragraph", p35["old"],
    p35["new"] + " The model of the only *mcrA*-carrying species among the MAGs (*Methanolobus*.1) grew by converting "
    "methyl groups to acetate instead of methane, because it could not reduce methanophenazine; its F420H2 "
    "dehydrogenase (Fpo), membrane-bound heterodisulfide reductase (HdrDE), Rnf complex and methyl-coenzyme M reductase "
    "were corrected to the gene-supported forms, after which it grew by methylotrophic methanogenesis on mono-, di- and "
    "trimethylamine, with the textbook methane yields, but not on methanol (its genome lacks *mtaB* and *mtaC*). A gene "
    "annotated by RAST as methane monooxygenase in two *Roseovarius* models is a propane 2-monooxygenase (*prmA*), and "
    "its reaction was disabled.",
    "curated methanogen (models/methanogen/CURATION.md) and PrmA block (pipeline/apply_prma_block.py)")
put("p36_biomass", 36, "replace",
    "constrained to its relative abundance in each simulated sample (based on read mapping), which scales the member’s metabolism to match its contribution in the sample’s consortium and specifically allows the simulated methane emission to quantitatively approximate the organism’s in situ production",
    "fixed, and the methane (or trimethylamine, TMA) each model released was expressed per gram of biomass formed; in "
    "these models this yield does not depend on the growth rate imposed, and it was combined with each member’s "
    "measured biomass in each core (see below)",
    "per-gram framing replaces 'biomass = relative abundance / DNA'")
put("p36_ratio_binds", 36, "replace", "to create a realistic consumption and excretion profile.",
    "to create a realistic consumption and excretion profile. The ratio constraint bound in every model. A C–P lyase "
    "releases one methane per methylphosphonate it cleaves, and the models cleave methylphosphonate only for the "
    "phosphorus they assimilate, so each MPn-consumer model released a fixed amount of methane per gram of biomass "
    f"formed ({kpos07.min():.3f}–{kpos07.max():.3f} mmol g⁻¹ at a ratio of 0.07, depending on the template’s biomass), "
    f"and at most the phosphorus content of biomass (about {kP:.2f} mmol g⁻¹) if all phosphorus came from "
    "methylphosphonate.", "ratio binds in all models; biochemistry of the yield")
put("p36_betaine_pergram", 36, "insert_before", "One TMA was excreted for each consumed betaine",
    "For the per-gram comparison the uptake bound was lifted, and the least TMA excreted per gram of biomass formed at "
    "the optimum was used (more TMA was possible without changing the optimum). ", "betaine per-gram yield")
put("p36_methanogen", 36, "insert_before", "The metabolic models were finally minimized",
    "The methanogen model was simulated on methylamines without an external electron acceptor, minimizing methylamine "
    "uptake; methane is the end product of its energy metabolism, so its yield per gram of biomass formed is set by the "
    "ATP the cell needs to grow. ", "methanogen protocol")
put("p38", 38, "replace_paragraph", "Linear regression between the experimental and simulated",
    "The modeled yields were compared with the measured emissions through the measured biomass of each member in each "
    "core.", "the regression and the unit conversion of the original comparison are replaced (P39, P40)")
put("p40", 40, "replace_paragraph", "for direct comparison and consistent flux interpretation",
    "The biomass of each member in each core was its relative abundance times the soil DNA concentration of each depth, "
    "summed over the two layers (0–5 and 5–15 cm) with a soil bulk density of 0.7 g cm⁻³ estimated from other studies "
    "of this region, and with DNA taken as 10% of dry biomass. The product of a guild’s methane yield and its biomass is "
    "the methane it produces per unit of growth; dividing the measured emission by this product gives the specific "
    "growth rate µ at which the guild’s measured biomass would supply the measured emission (doubling time ln 2 / µ), "
    "and, for the methanogen, dividing by its biomass alone gives the specific methanogenic activity required. If DNA "
    f"were 3% of dry biomass, every biomass would be {DNA3:.1f} times larger and every required rate {DNA3:.1f} times "
    "smaller. Spearman rank correlations (exact permutation p values) between each guild’s methane per unit of growth "
    "and the measured emissions of the six unrestored cores test whether the models follow the differences among "
    "cores. To include guild members not recovered as MAGs, each guild’s biomass was also scaled, in each sample, by "
    "the copies of its marker gene per genome equivalent in all sequenced reads (*mcrA* for methanogens, *phnJ* for MPn "
    "consumers, *grdI* for betaine reducers; see Metagenomics) divided by its MAGs’ genomes per genome equivalent. The "
    "required specific methanogenic activity does not depend on the methane yield; the doubling times computed from "
    "the read-based methanogen biomass assume that every *mcrA*-carrying organism has the yield of the *Methanolobus* "
    "model.", "explicit-rate comparison; the 10 % DNA constant and its sensitivity are stated; read-based scaling")
put("p69_dna", 69, "replace",
    "Soil DNA density was significantly positively correlated with methane flux (Figure S7), therefore, we only used methane emissions and fluxes that were normalized by DNA density for all modeling efforts.",
    "Soil DNA concentration was weakly and not significantly correlated with methane flux (Figure S7; R^2^ = 0.22, "
    "Spearman ρ = 0.52, p = 0.08, n = 12 cores).", "DNA sentence; no DNA normalization in the per-gram comparison")
put("p69_mepn", 69, "replace",
    "The simulated methane production of the MPn-consuming community members correlated remarkably well (R2 = 0.966) with the experimental methane emissions (Figure S8), and was most correlated when the phosphate source ratio was 0.07 (at least 1/14 of the total phosphate consumed in methylphosphonate versus inorganic phosphate).",
    "At a phosphate source ratio of 0.07 (at least 0.07 mol of methylphosphonate per mol of inorganic phosphate, about "
    f"6.5% of the phosphate taken up), the {n_mep} MPn-consumer models released {k07:.3f} mmol of methane per gram of "
    f"biomass formed (median; Figure 7a). Given their measured biomass, supplying the measured emissions of the six "
    f"unrestored cores would require these organisms to double every {rng_t('MePn, ratio 0.07')} (Figure 7b).",
    "per-gram MPn result")
put("p69_leadin", 69, "replace", "The line-of-best-fit equation from this regression", "",
    "the regression equation (P70) is removed")
put("p71", 71, "replace_paragraph", "reveals that generally the models predict half",
    f"At a ratio of 0.33 (at least 0.33 mol of methylphosphonate per mol of inorganic phosphate, about 25% of the "
    f"phosphate taken up) the yield rose to {k33:.3f} mmol g⁻¹, and the required doubling times to "
    f"{rng_t('MePn, ratio 0.33')}; even with methylphosphonate as the only phosphorus source the yield reached only "
    f"{kP:.2f} mmol g⁻¹, and the doubling times {rng_t('MePn, only P source')}. C–P lyase (*phnJ*) genes were "
    f"{fmep.min():.1f}–{fmep.max():.1f} times as abundant in all sequenced reads as in the MPn-consumer MAGs, so "
    f"counting MPn consumers not recovered as MAGs lengthens these doubling times only to {rng_t('MePn, ratio 0.07' + RB)} "
    f"at a ratio of 0.07 and {rng_t('MePn, only P source' + RB)} with methylphosphonate as the only phosphorus source "
    "(Figure 7b, open symbols). Such rates are implausibly fast for "
    "sediment microorganisms, so assimilation of methylphosphonate by the modeled consumers cannot supply the measured "
    "emissions unless methylphosphonate is degraded well beyond the phosphorus they assimilate (for example, with the "
    "phosphate excreted or stored as polyphosphate), which the models do not represent.",
    "per-gram MPn result at 0.33 and with MPn the only P source; the loophole the models cannot close")
put("p73_track", 73, "replace", "and yield an almost identical fit to the best fitting regression (R2 = 0.968).", "",
    "regression removed")
put("p73_obligacy", 73, "replace",
    "It is worth noting that unconstraining the phosphate source ratio produced almost an identical fit (R2 = 0.970) and still resulted in methylphosphonate consumption by several marine organisms in the Rhodobacteraceae family: two species in the Roseovarius genera and one species in the Albimonas genera.",
    f"None of the {n_mep} models needed methylphosphonate to grow: each also grew with inorganic phosphate, or with "
    "methylphosphonate, as its only phosphorus source, so the methylphosphonate taken up, and with it the methane, was "
    "set by the assumed ratio. Their methane per unit of growth did not follow the differences in emissions among the "
    f"six cores (Spearman ρ = {num(rho['MePn, ratio 0.07'][0], '.2f')}, n = 6), because it is proportional to their biomass, and "
    "most of it came from the Rhodobacteraceae family.", "obligacy, ranking, family")
put("p73_family", 73, "replace", "which supports that this taxa is uniquely responsible for the methane emissions.",
    f"which is consistent with this family contributing to the methane emissions (Rhodobacteraceae accounted for "
    f"{fam.get('Rhodobacteraceae', 0) * 100:.0f}% of the MPn consumers’ methane per unit of growth, *Roseovarius* alone "
    f"{gen.get('Roseovarius', 0) * 100:.0f}%).", "family shares from yield x measured biomass (PrmA blocked)")
put("p73_eq_comment", 73, "comment", "This family is by far the most abundant",
    "Review (2026-11-24): the regression equations (P70, P72, P75), the unit-conversion equation (P39) and their "
    "lead-ins are tracked deletions. Simulated areal emissions need an assumed growth rate: the original pipeline set "
    "each member's growth to its relative abundance divided by DNA, and the 2026-10-26 to 11-23 revisions in effect "
    "multiplied the methane per gram of biomass formed by the standing biomass per hour (every member growing at "
    "1 h⁻¹); any fit only rescales with that choice. The model output is now the methane per gram of biomass formed, and the "
    "comparison is the growth each guild would need (Figure 7). For the implausibility of minute-to-hour doubling "
    "times in sediments, consider citing Hoehler and Jørgensen 2013 (Nat Rev Microbiol 11:83–94). LibreOffice ignores "
    "revisions on equations and still shows them.", "equation removal")
put("p74", 74, "replace",
    "The simulated methane production of the betaine-consuming community members correlated slightly less well (R2 = 0.862) with the experimental methane emissions (Figure 7),",
    f"The {n_bet} betaine-reducer models excreted {tma.min():.0f}–{tma.max():.0f} mmol of TMA per gram of biomass "
    f"formed; at the assumed 0.06 methane per TMA (Methods), their measured biomass would supply the measured emissions "
    f"at doubling times of {rng_t('Betaine, 0.06 CH4 per TMA')}, and their methane per unit of growth was not correlated "
    f"with the emissions either (Spearman ρ = {num(rho['Betaine, 0.06 CH4 per TMA'][0], '.2f')}, n = 6; Figure 7).",
    "per-gram betaine result")
put("p74_leadin", 74, "replace", "with the line-of-best-fit", "", "the regression equation (P75) is removed")
put("p76", 76, "replace_paragraph", "requiring substantially more adjustments in the slope and th",
    "Methane from betaine depends on a second organism that converts the TMA, and only one species among the MAGs "
    f"carries *mcrA*: *Methanolobus*.1, whose curated model released {km:.0f} mmol of methane per gram of biomass formed "
    f"on methylamines ({KM['coupling_range'][0]:.0f}–{KM['coupling_range'][1]:.0f} mmol g⁻¹ across published ion "
    f"couplings), about {fold:,.0f} times the MPn consumers’ yield at a ratio of 0.07. Because a methanogen also "
    "spends methane on maintenance, its required specific methanogenic activity is the more direct measure: given its "
    f"measured biomass, the measured emissions require {q[R1].min():.0f}–{q[R1].max():.0f} mmol g⁻¹ h⁻¹ in the R1 cores "
    f"and {q[R2].min():.0f}–{q[R2].max():.0f} mmol g⁻¹ h⁻¹ in the R2 cores, equivalent to growth-coupled doubling "
    f"times of {rng_t('Methanogen, methylamines', R1)} and {rng_t('Methanogen, methylamines', R2)}. Methanogens not "
    "recovered as MAGs lower these requirements but do not remove them. In all sequenced reads of the six cores, *mcrA* "
    f"was {mcrA_pooled:.1f} times as abundant as the *Methanolobus* MAG (several-fold only in three R2 samples where the "
    f"MAG was scarce and *mcrA* reads few, {mcrA_big_reads[0]}–{mcrA_big_reads[1]} per sample), and "
    f"{other_share * 100:.0f}% of the *mcrA* reads came from other methanogens, mostly *Methanohalophilus*, "
    f"*Methanogenium* and *Methanococcoides*, {other_unasm * 100:.0f}% of them in no assembled contig. With the "
    "methanogens’ biomass scaled to the *mcrA* reads, the measured emissions require "
    f"{q_rb[R1].min():.0f}–{q_rb[R1].max():.0f} mmol g⁻¹ h⁻¹ in the R1 cores and {q_rb[R2].min():.0f}–{q_rb[R2].max():.0f} "
    f"mmol g⁻¹ h⁻¹ in the R2 cores (doubling times of {rng_t('Methanogen, methylamines' + RB, R1)} and "
    f"{rng_t('Methanogen, methylamines' + RB, R2)}). Doubling every few hours, as in the R2 cores, still seems unlikely "
    "for sediment methanogens unless non-growth methanogenesis or a lower DNA content of biomass add to the methane. The "
    "standing porewater pools of mono-, di- and trimethylamine could supply the measured "
    f"emissions if turned over every {turn.min():.1f}–{turn.max():.1f} days, and supplying the TMA for complete "
    f"conversion would need the betaine reducers to double only every {rng_t('Betaine, 2.25 CH4 per TMA')}. Neither "
    f"the methanogens’ biomass (Spearman ρ = {num(rho['Methanogen, methylamines'][0], '.2f')} for the MAG and "
    f"{num(rho_rb['Methanogen, methylamines'][0], '.2f')} counting all *mcrA* reads, p ≥ "
    f"{min(rho['Methanogen, methylamines'][1], rho_rb['Methanogen, methylamines'][1]):.2f}) nor that of any other "
    "guild followed the differences in emissions among the six cores.",
    "betaine -> TMA -> methanogen chain; required specific activity first; R1/R2 contrast; methanogens in the reads")
assert dm[R2].min() > 0.75 and max(dm[R2].max(), dm_rb[R2].max()) < 5.5, "'one to five hours' in the R2 cores"
put("p95_models", 95, "replace",
    "This was mechanistically corroborated by metabolic models that reproduced the experimental methane emissions with high accuracy (R2 = 0.978), and demonstrate that methane is generated stoichiometrically from methylphosphonate.",
    "Metabolic models show why: a C–P lyase releases one methane per methylphosphonate cleaved for phosphorus, so the "
    f"modeled MPn consumers released only {k07:.2f}–{kP:.2f} mmol of methane per gram of biomass formed, about "
    f"{fold:,.0f} times less than the methylotrophic methanogen, for which methane is the product of its energy "
    "metabolism. Their measured biomass could supply the measured emissions only at doubling times of minutes to "
    "hours, whereas the methanogens’ would need doubling times of hours to days in the R1 cores and one to five hours "
    "in the R2 cores, also when the methanogens found only in the reads are counted; neither guild’s biomass followed "
    "the differences in emissions among the six cores.",
    "per-gram Discussion")
put("p95_param", 95, "replace",
    "The metabolic models only parameterized the relative abundances of each member and trace consumption of methylphosphonate, which encourages that these are natural outcomes of the metabolism and not forced behaviors during optimization.",
    "The MPn-consumer models only parameterized a minimum ratio of methylphosphonate to inorganic phosphate uptake: "
    "every model released a fixed amount of methane per gram of biomass formed, set by the ratio and the phosphorus "
    "content of its biomass, so the models translate MPn-consumer biomass and the assumed ratio into methane rather "
    "than testing differences in their metabolism.", "parameterization")
put("p99", 99, "replace",
    "Thirdly, metabolic models reproduced methane emissions from trace amounts of methylphosphonate in near stoichiometric equivalence, supporting the experimental findings. The metabolic models systematically underestimated the quantity of methane emitted, which suggests",
    f"Thirdly, metabolic models showed that assimilating methylphosphonate releases only {k07:.2f}–{kP:.2f} mmol of "
    "methane per gram of biomass formed, so the measured emissions would require the MPn consumers to double every "
    "minutes to hours, whereas methylotrophic methanogenesis releases about "
    f"{fold:,.0f} times more methane per gram. That the MPn models cannot reach the measured emissions at plausible "
    "growth suggests", "per-gram Conclusions")
put("p100_comment", 100, "comment", "We therefore propose methylphosphonate as the primary source of methane",
    "Consider revising this conclusion. Per gram of biomass formed, methylphosphonate assimilation releases about "
    f"{fold:,.0f}-fold less methane than methylotrophic methanogenesis, and the measured emissions would require the MPn "
    "consumers to double every minutes to hours; the methanogens would need doubling times of hours to days in the R1 "
    "cores and one to five hours in the R2 cores, counting the Methanolobus MAG alone or all mcrA reads (Figure 7). In the "
    "re-assembled MAG set no guild's abundance tracks "
    "methane more closely than random sets of species (Table 5). The support for methylphosphonate rests on the "
    "substrate-addition (Figure 4) and culturing (Figure 6) results; the models fit methylphosphonate only if it is "
    "degraded far beyond assimilatory phosphorus demand. For example: 'We propose that methylphosphonate degradation "
    "contributes to the methane …, alongside methylotrophic methanogenesis'.", "hedge")
put("fig7", 129, "replace_paragraph", "Figure 7. Assessing methylphosphonate and betaine metabolism",
    "Figure 7. Methane released per gram of biomass formed by the metabolic models (a), and the growth each guild’s "
    "measured biomass would need to supply the measured methane emissions of the six unrestored cores (b). (a) Each "
    f"point is a model: the {n_mep} methylphosphonate (MPn)-consuming models at an MPn:phosphate uptake ratio of 0.07 "
    "and 0.33 and with MPn as the only phosphorus source (blue, light to dark), the "
    f"{n_bet} betaine-reducing models with methane assumed to be 0.06 of the TMA excreted (orange), and the curated "
    f"*Methanolobus* model on methylamines (green; bar, {KM['coupling_range'][0]:.0f}–{KM['coupling_range'][1]:.0f} "
    "mmol g⁻¹ across published ion couplings); black lines are medians. (b) Doubling time (ln 2 / µ) at which each "
    "guild’s measured biomass in each core (relative abundance × soil DNA × layer thickness, DNA taken as 10% of dry "
    "biomass) would supply the measured emission by growth-coupled methane release (filled symbols); open symbols "
    "scale each guild’s biomass by its marker gene’s abundance in all sequenced reads relative to its MAGs (*phnJ*, "
    "*grdI*, *mcrA*); dotted lines mark 10 min, 1 h, "
    f"1 day and 1 week. {drop_id} did not grow and is not modeled. No guild’s yield × biomass is correlated with the "
    f"measured emissions across the six cores (Spearman ρ from {num(min(v[0] for v in rho.values()), '.2f')} to "
    f"{num(max(v[0] for v in rho.values()), '.2f')}; read-based, {num(min(v[0] for v in rho_rb.values()), '.2f')} to "
    f"{num(max(v[0] for v in rho_rb.values()), '.2f')}).",
    "Figure 7 replaced by figure7_pergram.png")

put("p96_unassembled", 96, "insert_before", "The unrestored systems are relatively well defined",
    "Much of this fraction was never assembled: of the reads of eight universal single-copy ribosomal-protein genes, "
    f"{unasm['R1']:.0f}% (R1), {unasm['R2']:.0f}% (R2), {unasm['R2A']:.0f}% (R2A) and {unasm['SF2']:.0f}% (SF2) lay in "
    "no contig of at least 1 kb, so the unbinned fraction mainly reflects organisms too rare or too diverse to assemble "
    "at this sequencing depth rather than contigs left out of bins. Counted in all reads, the guilds’ marker genes were "
    f"{fmep.min():.1f}–{fmep.max():.1f} (*phnJ*), {fbet.min():.1f}–{fbet.max():.1f} (*grdI*) and, pooled, "
    f"{mcrA_pooled:.1f} (*mcrA*) times as abundant as the guilds’ MAGs in the unrestored cores, which lengthens the "
    "doubling times required of each guild by about these factors but does not change the conclusions of the models "
    "(Figure 7b). ",
    "why the unbinned fraction is large (assembly-limited), and the read-based guild abundance")
put("p36_comment", 36, "comment", "The MAG-derived models were subjected",
    "2026 models (reassembly_2026/modelling_2026/models/README.md). Built with ModelSEEDpy (fork sync-cshenry-main "
    "311f720) on the curated templates of modelling_2026/templates/; RAST (arm A) won the pre-registered annotation "
    f"comparison. {drop_id} was dropped (cannot make ribose-5-phosphate; sedoheptulose-1,7-bisphosphatase, rxn01345, "
    "rescues it but fails two ATP tests). The Methanolobus model was curated by hand (models/methanogen/CURATION.md: "
    "19 gene-justified changes to Fpo, HdrDE, Rnf, Mcr and related reactions); its genome lacks mtaB/mtaC, so it cannot "
    "use methanol, and it has no betaine methyltransferase (mtgB). The Roseovarius 'methane monooxygenase' is PrmA "
    "(KOfam PrmA 1016 vs MmoX 228) and was disabled. Template curation queue: betaine and glycine reductase directions "
    "(rxn11678, rxn11676), Fpo (rxn24614) and Rnf, reversible N-cycle reactions, and the missing betaine reductase "
    "complex.", "model decisions and caveats")

DROP = [(36, "replace", "constrained to its relative abundance"), (38, "replace", "Linear regression between"),
        (40, "replace", "The values were further normalized"), (69, "comment", "The line-of-best-fit"),
        (71, "replace", ""), (71, "comment", ""), (74, "comment", "with the line-of-best-fit")]
REPL = {(36, "comment", "The MAG-derived models"): "p36_comment", (18, "replace", "This was corroborated"): "abstract", (18, "comment", "a primary methanogenic"): "p18_comment",
        (35, "replace_paragraph", ""): "p35", (36, "replace", "to create a realistic"): "p36_ratio_binds",
        (69, "replace", "Soil DNA density"): "p69_dna", (69, "replace", "The simulated methane production"): "p69_mepn",
        (73, "replace", "and yield an almost"): "p73_track", (73, "replace", "It is worth noting"): "p73_obligacy",
        (73, "replace", "which supports that"): "p73_family", (74, "replace", "The simulated methane production"): "p74",
        (76, "replace_paragraph", ""): "p76", (95, "replace", "This was mechanistically"): "p95_models",
        (95, "replace", "The metabolic models only"): "p95_param", (99, "replace", "Thirdly, metabolic models"): "p99",
        (100, "comment", "We therefore propose"): "p100_comment", (129, "replace_paragraph", ""): "fig7"}
out = copy.deepcopy(base)
edits, used, dropped = [], set(), 0
for e in base["edits"]:
    if any(e["paragraph"] == p and e["kind"] == k and e["old"].startswith(pre) for p, k, pre in DROP):
        dropped += 1
        continue
    key = next((v for (p, k, pre), v in REPL.items() if e["paragraph"] == p and e["kind"] == k and e["old"].startswith(pre)), None)
    if key is None:
        edits.append(e)
        continue
    assert key not in used, key
    used.add(key)
    edits.append(NEW[key])
assert used == set(REPL.values()), set(REPL.values()) - used
edits += [NEW[k] for k in NEW if k not in used]
out["edits"] = edits
out["note"] = ("r232_review_edits_2026-11-23.json with the modelling text rewritten in the per-gram framing by "
               "make_modelling_edits_2026-11-24.py (dated 2026-11-24).")
out["facts_2026_11_24"] = json.loads(json.dumps(FACTS, default=float))
json.dump(out, open(HERE / "r232_review_edits_2026-11-24.json", "w"), indent=1, ensure_ascii=False)
print(f"wrote r232_review_edits_2026-11-24.json: {len(edits)} edits ({len(used)} replaced, {dropped} dropped, "
      f"{len(NEW) - len(used)} new)")
print(json.dumps(out["facts_2026_11_24"], indent=0)[:1500])
