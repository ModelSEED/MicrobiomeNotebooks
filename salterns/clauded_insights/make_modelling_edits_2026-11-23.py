#!/usr/bin/env python
"""Modelling edits for the 2026 guild models of the final MAG set: writes r232_review_edits_2026-11-23.json from
r232_review_edits_2026-11-19.json, replacing every modelling entry (Abstract P18/19, Methods P35-P36, Results P69-P76,
Discussion P95, Conclusions P99-P100, Figure 7 caption P129) with text generated from
  reassembly_2026/modelling_2026/simulations/{consistent_basis_final.json, reference_k.json, fluxes_final.json}
  reassembly_2026/modelling_2026/models/models_summary.tsv, annotation_comparison/summary.json
and adding: the P35 Methods rewrite (replace_paragraph), a P73 comment recommending deletion of the Albimonas
sentence, and a P37 comment on Figure S10. Non-modelling entries (Halomonas isolates, P37/P38/P40 wording, P69 DNA sentence, P71 ratio
wording, P77) are kept unchanged. The wording's preconditions are asserted, so a model refresh that changes them fails.
Run: ~/Documents/py_venv/bin/python make_modelling_edits_2026-11-23.py
"""
import copy
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
S = HERE.parent
M = S / "reassembly_2026" / "modelling_2026"
SIM = M / "simulations"
DATE = "2026-11-23T00:00:00Z"
CORES = ["R1_A", "R1_B", "R1_C", "R2_A", "R2_B", "R2_C"]
SRC = ("reassembly_2026/modelling_2026/simulations/consistent_basis_final.json (consistent_basis_final.py; fluxes from "
       "simulate_guilds.py on TF); models: reassembly_2026/modelling_2026/models/README.md, models_summary.tsv")

CB = json.load(open(SIM / "consistent_basis_final.json"))
REF = json.load(open(SIM / "reference_k.json"))
FX = json.load(open(SIM / "fluxes_final.json"))
MS = pd.read_csv(M / "models" / "models_summary.tsv", sep="\t")
AC = json.load(open(M / "annotation_comparison" / "summary.json"))
S10 = json.load(open(M / "figureS10" / "figureS10_summary.json"))
IDS = json.load(open(S / "correlational" / "final585" / "iterativeIDs_final585.json"))
G = pd.read_csv(S / "correlational" / "final585" / "guilds_final585.csv")
A = CB["areal"]
x_meas = np.array([CB["measured_CH4"][c] for c in CORES])


def exact_p(key):
    y = np.array([A[key]["sim_by_core"][c] for c in CORES])
    r0 = stats.spearmanr(x_meas, y)[0]
    hits = sum(abs(stats.spearmanr(x_meas, y[list(p)])[0]) >= abs(r0) - 1e-12 for p in itertools.permutations(range(6)))
    return r0, hits / 720


def pct(v):
    return f"{v * 100:.0f}"


FRACS = {"one half": 1 / 2, "one third": 1 / 3, "one quarter": 1 / 4, "one fifth": 1 / 5, "one sixth": 1 / 6,
         "one seventh": 1 / 7, "one eighth": 1 / 8, "two thirds": 2 / 3, "three quarters": 3 / 4}


def words(v):
    name, f = min(FRACS.items(), key=lambda kv: abs(kv[1] - v) / kv[1])
    assert abs(f - v) / f < 0.12, (v, name)
    return name


def num(v, spec):
    return format(v, spec).replace("-", "−")


def eq(key):
    s = A[key]
    return f"{num(s['slope'], '.2f')} × experimental {'+' if s['intercept'] >= 0 else '−'} {abs(s['intercept']):.0f}"


# ---- facts ----------------------------------------------------------------------------------------------------------
kept = MS[MS.kept == True]
mep_k, bet_k = kept[kept.guild == "mepn_consumer"], kept[kept.guild == "betaine_reducer"]
n_mep, n_bet = len(mep_k), len(bet_k)
n_mep_sp, n_bet_sp = int(G.mepn_consumer.sum()), int(G.betaine_reducer.sum())
dropped = MS[(MS.kept != True) & (MS.guild == "mepn_consumer")]
assert len(dropped) == 1 and n_mep == CB["n_mepn_models"] and n_bet == CB["n_betaine_models"] == n_bet_sp
drop_id = IDS[dropped.mag_id.iloc[0]]
assert not mep_k.mepn_obligate.astype(str).eq("True").any(), "an MePn model is MePn-obligate"
tpl = kept[kept.guild != "methanogen"].template.value_counts().to_dict()
n_gn, n_gp, n_ar = (tpl.get(t, 0) for t in ("GramNegModelTemplateV7.json", "GramPosModelTemplateV7.json", "ArchaeaTemplateV6.json"))
binds = [abs(abs(d[r]["MePn"]) / abs(d[r]["phosphate"]) - float(r)) < 1e-6 for d in REF["mepn"].values() for r in ("0.07", "0.33")]
assert all(binds), "the ratio constraint does not bind in every model"
k07 = {m: d["0.07"]["methane"] / REF["reference_biomass"] for m, d in REF["mepn"].items()}
kpos = [v for v in k07.values() if v > 1e-9]
mmo = [m for m, d in REF["mepn"].items() if "mmo_blocked" in d and d["0.07"]["methane"] < 1e-9]
assert len(mmo) == 1
mmo_id = IDS[mmo[0]]
p07, p33 = A["MePn_0.07_38"], A["MePn_0.33_38"]
m07b, m33b = A["MePn_0.07_38_mmo_blocked"], A["MePn_0.33_38_mmo_blocked"]
bet, bmin, bmax = A["betaine_18"], A["betaine_18_fva_min"], A["betaine_18_fva_max"]
assert 0.10 < p07["model_over_measured_total"] < 0.20 and abs(m07b["model_over_measured_total"] / p07["model_over_measured_total"] - 1) < 0.10
assert bet["model_over_measured_min"] > 1, "betaine no longer exceeds the measured emissions in every core"
r07, pp07 = exact_p("MePn_0.07_38")
r33, pp33 = exact_p("MePn_0.33_38")
rb, ppb = exact_p("betaine_18")
assert pp07 > 0.5 and pp33 > 0.5 and ppb > 0.05
fold = p33["model_over_measured_total"] / p07["model_over_measured_total"]
fam = CB["MePn_0.07_38_family_share"]
gen = CB["MePn_0.07_38_genus_share"]
bt = CB["betaine_tma"]
ab = pd.DataFrame({m: {c: FX["betaine"][m][c]["abundance"] for c in CORES} for m in FX["betaine"]}).T
ab_fold = ab.sum().max() / ab.sum().min()
crit = AC["criteria"]
blk = sorted(crit[a]["e_blocked_change_pp"] for a in ("B", "C", "D"))
gain = sorted({IDS[m] for a in ("C", "D") for m, _ in crit[a]["d_models_gaining_free_atp"]})
assert AC["decision"] == "A" and all(not crit[a]["e_pass"] for a in ("B", "C", "D"))
w07, w33 = words(p07["model_over_measured_total"]), words(p33["model_over_measured_total"])
t07, t33 = pct(p07["model_over_measured_total"]), pct(p33["model_over_measured_total"])
FACTS = dict(n_mep=n_mep, n_bet=n_bet, drop_id=drop_id, mmo_id=mmo_id, t07=t07, t33=t33, w07=w07, w33=w33,
             rho07=r07, p07=pp07, rho_b=rb, p_b=ppb, fold=fold, k=(min(kpos), max(kpos)))

# ---- new entries --------------------------------------------------------------------------------------------------------
NEW = {}


def put(key, paragraph, kind, old, new, rationale, mirror_index=None):
    e = {"paragraph": paragraph, "kind": kind, "old": old, "new": new, "rationale": rationale, "sources": SRC,
         "origin": "modelling_2026", "date": DATE}
    if mirror_index is not None:
        e["mirror_index"] = mirror_index
    NEW[key] = e


put("abstract", 18, "replace",
    "This was corroborated by metabolic modeling, where trace consumption of methylphosphonic acid was sufficient to reproduce the experimental methane emissions.",
    f"Metabolic modeling partly corroborated this: trace consumption of methylphosphonic acid by the "
    f"methylphosphonate-consuming community members could supply about {w07} to {w33} of the measured methane "
    f"emissions, depending on the assumed share of methylphosphonate in phosphate uptake, although the simulated emissions "
    f"did not track the differences among cores.",
    f"2026 models of the final MAG set ({n_mep} MePn-consumer models): {t07}% of the measured methane at ratio 0.07 and "
    f"{t33}% at 0.33; Spearman ρ = {r07:.2f} (n = 6).", mirror_index=19)
put("p35", 35, "replace_paragraph", "Metabolic models were developed using dereplicated medium and high quality",
    f"Metabolic models were built for the representative MAG of each MPn-consumer and betaine-reducer species of the "
    f"re-assembled MAG set ({n_mep_sp} and {n_bet_sp} species; Table 4), one model per species and guild. Prodigal genes "
    "were annotated with the RAST service (53), and draft models were built with ModelSEEDpy, which implements the KBase "
    "modeling pipeline (52), from the current ModelSEED templates with our curations (Gram-negative v7, Gram-positive v7 "
    f"and Archaea v6 templates for {n_gn}, {n_gp} and {n_ar} models, respectively). RAST alone was chosen over drafts "
    "that also used KOfam, eggNOG-mapper v2.1.13 and Bakta v1.12.0 annotations by criteria fixed in advance "
    "(gene-associated and gapfilled reactions, guild-pathway support, ATP production without a carbon source, and the "
    "fraction of blocked reactions): the added annotations raised the number of gene-associated reactions, but they "
    f"also raised the fraction of blocked reactions by {blk[0]:.1f}–{blk[-1]:.1f} percentage points, beyond the "
    "pre-set limit of 2, and the two that included eggNOG-mapper let one model make ATP from ammonia and oxygen. ATP "
    "production was corrected on 54 test media, and the models were gapfilled as in the original protocol: MPn "
    "consumers first on aerobic lactate with phosphate and then on lactate with MPn as the only phosphorus source, after "
    "a C–P lyase pathway and an MPn transporter were added where absent; betaine reducers on aerobic lactate, after which "
    "betaine reductase, glycine reductase and TMA export were added. Reactions that let a model oxidize ammonia "
    "(reversed nitrite and hydroxylamine reductases) or consume methane without annotation support were restricted to "
    "their reductive or methane-producing directions, and models that lost growth were gapfilled again. One "
    f"MPn-consumer model ({drop_id}) could not grow and was excluded, leaving {n_mep} MPn-consumer and {n_bet} "
    "betaine-reducer models; none of the MPn-consumer models required MPn to grow.",
    "Methods for the 2026 models (models/README.md, annotation_comparison/RESULTS.md). The draft's KBase MS2 build of the "
    "2020 MAGs and the methanogen medium sentence no longer apply: the methanogen model is not simulated (comment on P36).")
put("p36_ratio_binds", 36, "replace", "to create a realistic consumption and excretion profile.",
    "to create a realistic consumption and excretion profile. The ratio constraint bound in every model, so each "
    "MPn-consumer model produced a fixed amount of methane per unit of biomass, set by the ratio and the phosphorus "
    f"content of its template biomass ({FACTS['k'][0]:.3f}–{FACTS['k'][1]:.3f} mmol CH~4~ per g of biomass formed at a "
    "ratio of 0.07, depending on the template). A methane monooxygenase annotated by RAST in one *Roseovarius* model was "
    "kept, because it is supported by the annotation (see Results).",
    "ratio binds in 76/76 reference runs (reference_k.json); k clusters by template biomass P content.")
put("p36_betaine", 36, "replace", "The betaine-utilizing organisms were simulated with a betaine consumption rate of",
    f"The betaine-reducing organisms (the {n_bet} betaine-reducer models; Table 4) were simulated with a maximum betaine "
    "uptake rate of", "18 betaine-reducer species, one model each.")
put("p36_comment", 36, "comment", "The MAG-derived models were subjected",
    f"2026 models (added 2026-11-23; reassembly_2026/modelling_2026/models/README.md). Built on TF with ModelSEEDpy "
    "(fork sync-cshenry-main 311f720) and the curated templates of modelling_2026/templates/; arm A (RAST) won the "
    "pre-registered annotation comparison (annotation_comparison/RESULTS.md). Decisions for the authors: (1) "
    f"{drop_id} was dropped because it cannot make ribose-5-phosphate; one reaction (sedoheptulose-1,7-bisphosphatase, "
    "rxn01345) rescues it but fails two ATP tests. (2) The methanogen model (Methanolobus.1) grows as a methyl-to-acetate "
    "acetogen without making methane, so it was not simulated; the manuscript's 'dimethylamine and hydrogen' medium has no "
    "surviving code (the model grows on it, 1.32 h−1, but makes methane only apart from growth). (3) "
    f"{mmo_id} carries a RAST methane monooxygenase (rxn00844) that oxidizes its MePn methane at the lactate optimum; "
    f"closing it gives {pct(m07b['model_over_measured_total'])}% instead of {t07}% of the measured methane at 0.07. "
    "(4) Template curation queue: betaine reductase rxn11678 and glycine reductase rxn11676 are written in the synthesis "
    "direction; Fpo rxn24614 is '<' only; four N-cycle reactions are reversible; the betaine reductase roles have no "
    "complex. The earlier repair of three MePn-obligate 2020 models does not apply: no 2026 model needs MePn.",
    "replaces the 2026-11-08 repair comment")
s10m, s10b = S10["mepn_consumer"], S10["betaine_reducer"]
assert s10m["LLac.O2_rank_by_median"] == 1 and s10b["LLac.O2_rank_by_median"] > 5, "S10 wording assumptions changed"
assert s10b["median_by_medium"]["Pyr"] == 0 and s10b["median_by_medium"]["Fum"] == 0, "S10: pyruvate/fumarate wording"
anaer_lower = [m for m in s10b["media_with_lower_median"] if not m.endswith(".O2")]
put("p37_lactate", 37, "replace", "for all methylphosphonate- and betaine- consuming members",
    "for the methylphosphonate-consuming members",
    "Figure S10 regenerated for the 2026 models (reassembly_2026/modelling_2026/figureS10/)")
put("p37_betaine", 37, "insert_before", "Glycerol and acetate are other candidate carbon sources",
    f"For the betaine-reducing members, however, aerobic L-lactate required more reactions than "
    f"{len(s10b['media_with_lower_median'])} of the {s10b['n_media']} media tested, including anaerobic pyruvate and "
    f"fumarate (median {s10b['LLac.O2_median']:g} against 0; Figure S10). ",
    "Figure S10 regenerated: lactate is not among the cheapest media for the betaine reducers")
put("p37_s10", 37, "comment", "Lactate was selected as the simulated carbon source",
    "Figure S10 is regenerated for the 2026 models (figureS10_final.py: MSATPCorrection core gap-filling of each arm-A "
    "draft on the 22 ATP test media of the original figure, run with the build code on TF). Aerobic L-lactate remains "
    f"the cheapest medium for the MePn consumers, but for the betaine reducers it ranks {s10b['LLac.O2_rank_by_median']} "
    f"of {s10b['n_media']} (anaerobic pyruvate and fumarate need none), and no betaine model grows without O2 although these "
    "organisms (Bacillota such as Halanaerobium, M55B127) are anaerobic Stickland fermenters. Simulating the betaine "
    "reducers on an anaerobic medium (e.g. pyruvate) would be more realistic; the protocol was kept to match the "
    "original analysis.", "S10 regenerated")
put("s10_caption", 157, "replace_paragraph", "Figure S10.",
    f"Figure S10. Core gap filling of the draft metabolic models: the number of reactions that ATP correction had to add "
    f"or reverse for each draft to produce ATP on each of {s10m['n_media']} test media (a single carbon source, with or "
    f"without O~2~), for a) the {s10m['n_models']} methylphosphonate-consuming and b) the {s10b['n_models']} "
    "betaine-reducing models (species labels as in Table 4). Rows and columns are ordered by hierarchical clustering; "
    f"grey cells had no gap-filling solution. The red box marks the aerobic L-lactate medium used for the simulations, "
    f"which needed no reactions for most methylphosphonate consumers (median {s10m['LLac.O2_median']:g}, tied for the "
    f"lowest of the {s10m['n_media']} media) but a median of {s10b['LLac.O2_median']:g} for the betaine reducers, more "
    f"than {len(s10b['media_with_lower_median'])} other media, including anaerobic pyruvate and fumarate (median 0).",
    "Figure S10 replaced by figureS10_final.png")
put("p69_mepn", 69, "replace",
    "The simulated methane production of the MPn-consuming community members correlated remarkably well (R2 = 0.966) with the experimental methane emissions (Figure S8), and was most correlated when the phosphate source ratio was 0.07 (at least 1/14 of the total phosphate consumed in methylphosphonate versus inorganic phosphate).",
    "With a phosphate source ratio of 0.07 (at least 0.07 mol of methylphosphonate per mol of inorganic phosphate, about "
    f"6.5% of the phosphate taken up), the simulated methane production of the {n_mep} MPn-consumer models amounted to "
    f"{t07}% of the experimental methane emissions in total and did not track their differences among the six unrestored "
    f"cores (Figure 7; Spearman ρ = {r07:.2f}, exact p = {pp07:.1f}, n = 6; R^2^ = {p07['pearson_r2']:.2f}).",
    "consistent areal basis, 2026 models")
put("p69_comment", 69, "comment", "The line-of-best-fit equation from this regression",
    f"The equation below is replaced (tracked) by the fit of the 2026 models on the consistent areal basis at ratio 0.07: "
    f"simulated = {eq('MePn_0.07_38')}, in µmol CH4 m−2 d−1 (R² = {p07['pearson_r2']:.2f}, Spearman ρ = {r07:.2f}, n = 6). "
    "The draft's fits (R² 0.978/0.968/0.862) had both axes divided by DNA, and two depths have almost no DNA (R2_B_D1, "
    "R2_C_D2, 0.16 µg g−1 soil), so they were driven by a 1/DNA term. Here each model's methane per unit of simulated "
    "biomass (a constant per model) is rescaled to its biomass in the core (relative abundance × DNA × layer thickness). "
    "Word shows the equation swap as a tracked change; LibreOffice and some other viewers display both equations.",
    "equation swap explanation")
put("p71_total", 71, "replace", "reveals that generally the models predict half of the experimental methane flux.",
    f"reveals that the models predict about {w07} of the experimental methane flux in total ({t07}%), from "
    f"{pct(p07['model_over_measured_min'])}% to {pct(p07['model_over_measured_max'])}% among cores.", "consistent basis")
put("p71_033", 71, "replace", "the models almost perfectly recapitulate the experimental methane emissions",
    f"the simulated emissions reach {t33}% of the measured ones in total but range from "
    f"{pct(p33['model_over_measured_min'])}% to {pct(p33['model_over_measured_max'])}% among cores", "consistent basis")
put("p71_comment", 71, "comment", "If the phosphate ratio is",
    f"The equation below is replaced (tracked) by the 2026 models' consistent-basis fit at ratio 0.33: simulated = "
    f"{eq('MePn_0.33_38')}, in µmol CH4 m−2 d−1 (R² = {p33['pearson_r2']:.2f}, Spearman ρ = {r33:.2f}). The previous "
    "0.95x + 224 was the fit with both axes divided by DNA (R² = 0.968). Word shows the swap as a tracked change; "
    "LibreOffice shows both equations.", "equation swap explanation")
put("p73_track", 73, "replace", "and yield an almost identical fit to the best fitting regression (R2 = 0.968).",
    f"and likewise do not track the measured emissions among cores (Spearman ρ = {r33:.2f}, n = 6; R^2^ = "
    f"{p33['pearson_r2']:.2f}): the ratio rescales the simulated emissions ({fold:.1f}-fold from 0.07 to 0.33) but does "
    "not change their pattern.", "consistent basis")
put("p73_obligacy", 73, "replace",
    "It is worth noting that unconstraining the phosphate source ratio produced almost an identical fit (R2 = 0.970) and still resulted in methylphosphonate consumption by several marine organisms in the Rhodobacteraceae family: two species in the Roseovarius genera and one species in the Albimonas genera.",
    f"It is worth noting that none of the {n_mep} models needed methylphosphonate to grow: each also grew with inorganic "
    "phosphate, or with methylphosphonate, as its only phosphorus source, so the methylphosphonate taken up, and with it "
    f"the methane, was set by the assumed ratio. In one *Roseovarius* model ({mmo_id}), a RAST-annotated methane "
    "monooxygenase oxidized the methane released by the C–P lyase (all of it at a ratio of 0.07); with that reaction "
    f"closed, the simulated total rises to {pct(m07b['model_over_measured_total'])}% (0.07) and "
    f"{pct(m33b['model_over_measured_total'])}% (0.33) of the measured emissions. Most of the simulated methane came from "
    "models of the Rhodobacteraceae family.", "growth_triplet; mmo sensitivity")
put("p73_family", 73, "replace", "which supports that this taxa is uniquely responsible for the methane emissions.",
    "which is consistent with this family contributing to the methane emissions (Rhodobacteraceae models produced "
    f"{pct(fam['Rhodobacteraceae'])}% of the simulated methane at a ratio of 0.07, *Roseovarius* alone "
    f"{pct(gen['Roseovarius'])}%).", "family/genus shares")
put("p73_albimonas_comment", 73, "comment", "Roseovarius has been experimentally observed",
    "Consider deleting the closing Albimonas sentence: KOfam found no phnJ in any Albimonas MAG of the re-assembled set "
    "(the C–P lyase genes of the original Albimonas.1 are frameshifted), so Albimonas has no MPn-consumer model and no "
    "longer appears in these results. Reference 78 would then be uncited. (Not deleted here: the sentence ends inside "
    "an existing comment range and a citation link.)", "Albimonas")
put("p74", 74, "replace",
    "The simulated methane production of the betaine-consuming community members correlated slightly less well (R2 = 0.862) with the experimental methane emissions (Figure 7),",
    f"The simulated methane production of the betaine-reducing community members (the {n_bet} betaine-reducer models; "
    f"Table 4) was not significantly correlated with the experimental methane emissions either (Figure 7; Spearman ρ = "
    f"{rb:.2f}, exact p = {ppb:.2f}, n = 6; R^2^ = {bet['pearson_r2']:.2f}),", "consistent basis (betaine: pipeline conversion)")
put("p74_comment", 74, "comment", "with the line-of-best-fit",
    f"The equation below is replaced (tracked) by the fit of the {n_bet} 2026 betaine models: simulated = "
    f"{eq('betaine_18')}, in µmol CH4 m−2 d−1. The previous 0.24x + 4632 was the 15-MAG fit with both axes divided by DNA "
    "(R² = 0.862). Betaine-model caveats: the betaine reductase reaction rxn17220_c0 is added by hand (as in the original "
    "protocol) and carries the RAST betaine-reductase genes as its gene rule, because the templates have no complex for "
    "those roles; the models are aerobic (no betaine model grows without O2), so an anaerobic enzyme is simulated in "
    f"aerobic models; TMA has no unique optimum in {bt['not_unique']} of the {bt['cells']} model×core cells (with every "
    f"value at its FVA minimum or maximum, the total is {bmin['model_over_measured_total']:.1f}- or "
    f"{bmax['model_over_measured_total']:.1f}-fold the measured emissions). The magnitude depends on the assumed uptake "
    "bound, the 0.06 methane per TMA, the soil bulk density and the DNA content of biomass. Word shows the equation swap "
    "as a tracked change; LibreOffice shows both equations.", "betaine caveats")
put("p76", 76, "replace_paragraph", "requiring substantially more adjustments in the slope and th",
    f"which has a {'positive' if bet['slope'] > 0 else 'negative'} slope, but the simulated betaine-derived methane "
    f"exceeded the measured emissions in all six cores ({bet['model_over_measured_min']:.1f}- to "
    f"{bet['model_over_measured_max']:.1f}-fold; {bet['model_over_measured_total']:.1f}-fold in total). Its magnitude is "
    "bounded by two assumptions, the betaine uptake bound of 0.1 and the conversion of 0.06 methane per TMA, and within "
    "those bounds it was often not determined by the genomes: TMA excretion reached the uptake bound in "
    f"{bt['at_bound_0.1']} of the {bt['cells']} model × core simulations and was not unique in {bt['not_unique']} (with "
    "every value at its minimum or maximum, the total is "
    f"{bmin['model_over_measured_total']:.1f}- to {bmax['model_over_measured_total']:.1f}-fold the measured emissions). "
    f"Across cores, the summed TMA excretion varied {bt['summed_TMA_fold_range']:.1f}-fold while the DNA-normalized "
    f"abundance of the betaine reducers varied {ab_fold:.0f}-fold, so this comparison is not a strong test of the betaine "
    "hypothesis.", "betaine regime")
put("p95_models", 95, "replace",
    "This was mechanistically corroborated by metabolic models that reproduced the experimental methane emissions with high accuracy (R2 = 0.978), and demonstrate that methane is generated stoichiometrically from methylphosphonate.",
    f"Metabolic models of the MPn consumers produced methane equivalent to about {w07} of the experimental emissions in "
    f"the unrestored cores at a phosphate source ratio of 0.07 (about {w33} at 0.33), with one methane generated per "
    "methylphosphonate consumed, but neither these models nor the betaine models tracked the differences in emissions "
    f"among the six cores (Spearman ρ = {r07:.2f} and {rb:.2f}, p ≥ {min(pp07, ppb):.2f}), and the betaine models exceeded "
    "the measured emissions in every core.", "consistent basis")
put("p95_param", 95, "replace",
    "The metabolic models only parameterized the relative abundances of each member and trace consumption of methylphosphonate, which encourages that these are natural outcomes of the metabolism and not forced behaviors during optimization.",
    "The metabolic models only parameterized the relative abundances of each member and a minimum ratio of "
    "methylphosphonate to inorganic phosphate uptake. Every MPn-consumer model produced a fixed amount of methane per unit "
    "of biomass, set by the ratio and the phosphorus content of its template biomass, so the simulated emissions translate "
    "the biomass of MPn consumers and the assumed ratio into methane rather than testing differences in their metabolism.",
    "ratio binds in every model")
put("p99", 99, "replace",
    "Thirdly, metabolic models reproduced methane emissions from trace amounts of methylphosphonate in near stoichiometric equivalence, supporting the experimental findings. The metabolic models systematically underestimated the quantity of methane emitted, which suggests",
    "Thirdly, metabolic models showed that trace consumption of methylphosphonate, with one methane per methylphosphonate "
    f"consumed, can account for part of the methane emitted from the unrestored salterns ({t07}% at a phosphate source "
    f"ratio of 0.07 and {t33}% at 0.33), although the simulated emissions did not track the differences among cores. That "
    "the models underestimated the total methane emitted suggests", "consistent basis")
put("p100_comment", 100, "comment", "We therefore propose methylphosphonate as the primary source of methane",
    f"Consider hedging this conclusion. The 2026 models show that methylphosphonate consumption can supply part of the "
    f"measured methane (about {w07} at a ratio of 0.07, about {w33} at 0.33) but not the pattern among cores (Figure 7), "
    "and in the re-assembled MAG set no guild's abundance tracks methane more closely than random sets of species (Table "
    "5). The strongest support is the substrate-addition (Figure 4) and culturing (Figure 6) evidence. For example: 'We "
    "therefore propose methylphosphonate as the most consistent candidate source of methane…'.", "hedge")
put("fig7", 129, "replace_paragraph", "Figure 7. Assessing methylphosphonate and betaine metabolism",
    "Figure 7. Simulated versus measured methane emissions in the six unrestored cores. Orange: summed methane flux of "
    f"the methylphosphonate-consuming models ({n_mep} of the {n_mep_sp} MPn-consumer species; {drop_id} did not grow), "
    "each model's methane per unit of simulated biomass rescaled to its biomass in the core (relative abundance × soil "
    "DNA concentration of each depth, weighted by depth); filled circles at a phosphate source ratio of 0.07 and open "
    f"circles at 0.33. One *Roseovarius* model ({mmo_id}) releases no methane at 0.07 because its annotated methane "
    "monooxygenase oxidizes it. Blue: summed methane of the "
    f"{n_bet} betaine-reducing models (methane = 0.06 × TMA excretion), converted with the core’s volume-weighted DNA "
    "concentration because their TMA excretion is partly set by the uptake bound. Both axes are in µmol CH~4~ m^−2^ "
    "d^−1^ on log scales; the dotted line marks 1:1 agreement and the dashed lines are linear fits (methylphosphonate "
    f"models at 0.07, and betaine models). Methylphosphonate-derived methane is {t07}% of the measured emissions in total "
    f"at 0.07 ({pct(p07['model_over_measured_min'])}–{pct(p07['model_over_measured_max'])}% among cores; {t33}% at 0.33) "
    f"and does not track them (Spearman ρ = {r07:.2f}; R^2^ = {p07['pearson_r2']:.2f}); betaine-derived methane exceeds "
    f"them in every core ({bet['model_over_measured_min']:.1f}- to {bet['model_over_measured_max']:.1f}-fold) and is not "
    f"significantly correlated with them (ρ = {rb:.2f}, exact p = {ppb:.2f}; R^2^ = {bet['pearson_r2']:.2f}). Guild "
    "abundance versus methane is shown in Figure 9 and Table 5.", "Figure 7 replaced by figure7_final.png")
put("p18_comment", 18, "comment", "a primary methanogenic substrate",
    f"Review: the 2026 models supply about {w07} to {w33} of the measured methane and do not track the differences among "
    "cores (Figure 7), so the modeling supports methylphosphonate as a source of methane but not as the primary "
    "substrate. Consider 'an important' or 'a' instead of 'a primary', unless the in vitro data alone carry 'primary'.",
    "Abstract claim", mirror_index=19)

# ---- merge with the 2026-11-19 file ------------------------------------------------------------------------------------
REPLACES = {  # (paragraph, kind, old prefix) of the 11-19 entry -> key of its replacement (None = dropped)
    (18, "replace", "This was corroborated"): "abstract",
    (36, "replace", "to create a realistic"): "p36_ratio_binds",
    (36, "replace", "The betaine-utilizing"): "p36_betaine",
    (36, "comment", "The MAG-derived models"): "p36_comment",
    (69, "replace", "The simulated methane production"): "p69_mepn",
    (69, "comment", "The line-of-best-fit"): "p69_comment",
    (71, "replace", "reveals that generally"): "p71_total",
    (71, "replace", "the models almost perfectly"): "p71_033",
    (71, "comment", "If the phosphate ratio"): "p71_comment",
    (73, "replace", "and yield an almost"): "p73_track",
    (73, "replace", "It is worth noting"): "p73_obligacy",
    (73, "replace", "which supports that"): "p73_family",
    (74, "replace", "The simulated methane production"): "p74",
    (74, "comment", "with the line-of-best-fit"): "p74_comment",
    (76, "replace_paragraph", "requiring substantially"): "p76",
    (95, "replace", "This was mechanistically"): "p95_models",
    (95, "replace", "The metabolic models only"): "p95_param",
    (99, "replace", "Thirdly, metabolic models"): "p99",
    (100, "comment", "We therefore propose"): "p100_comment",
    (129, "replace_paragraph", "Figure 7. Assessing"): "fig7",
    (18, "comment", "a primary methanogenic"): "p18_comment",
}
old = json.load(open(HERE / "r232_review_edits_2026-11-19.json"))
out = copy.deepcopy(old)
used, edits = set(), []
for e in old["edits"]:
    key = next((v for (pg, kd, pre), v in REPLACES.items() if e["paragraph"] == pg and e["kind"] == kd
                and e["old"].startswith(pre)), "keep")
    if key == "keep":
        edits.append(e)
        continue
    assert key not in used, key
    used.add(key)
    n = NEW[key]
    assert n["old"] == e["old"] or n["kind"] == "replace_paragraph", (key, e["old"][:60])
    edits.append(n)
assert used == set(REPLACES.values()), set(REPLACES.values()) - used
edits += [NEW[k] for k in NEW if k not in used]           # p35, p37_s10, p73_albimonas(+comment)
out["edits"] = edits
out["note"] = ("Copy of r232_review_edits_2026-11-19.json with every modelling entry regenerated for the 2026 guild models of "
               "the final MAG set by make_modelling_edits_2026-11-23.py (dated 2026-11-23), plus the P35 Methods rewrite, the "
               "P73 Albimonas comment, and a P37 comment on Figure S10. 'paragraph' stays in July-export "
               "numbering; 'mirror_index' is the absolute paragraph index in the 2026-11-17/11-20 mirrors.")
out["facts"] = {k: (list(v) if isinstance(v, tuple) else v) for k, v in FACTS.items()}
json.dump(out, open(HERE / "r232_review_edits_2026-11-23.json", "w"), indent=1, ensure_ascii=False)
print(f"wrote r232_review_edits_2026-11-23.json: {len(edits)} edits ({len(used)} replaced, {len(NEW) - len(used)} new)")
print(json.dumps(out["facts"], indent=0, default=float))
