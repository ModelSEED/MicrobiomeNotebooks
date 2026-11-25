#!/usr/bin/env python
"""Synthesis of the methane sources (2026-11-25): writes r232_review_edits_2026-11-25.json from
r232_review_edits_2026-11-24.json (run make_modelling_edits_2026-11-24.py first). Adds, as tracked changes:
  - a Methods paragraph after P40 (pathway abundance, substrate pools, attribution at equal growth rates);
  - a Results paragraph after P78 (the end of the Results) presenting Figure 10;
  - a Conclusions paragraph after P99 attributing the measured methane among the hypothesized sources;
  - the P98 and P100 Conclusions sentences rewritten to agree with it (the P100 'Consider revising' comment is dropped
    because the revision is now in the text);
  - Figure 10 (reassembly_2026/modelling_2026/synthesis/figure10_pathways.png) and its caption, after Figure 9.
Numbers come from reassembly_2026/modelling_2026/synthesis/pathway_summary.json (pathway_summary.py); the wording's
preconditions are asserted. Run: ~/Documents/py_venv/bin/python make_synthesis_edits_2026-11-25.py
"""
import copy
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
S = HERE.parent
SYN = S / "reassembly_2026/modelling_2026/synthesis"
DATE = "2026-11-25T00:00:00Z"
SRC = ("reassembly_2026/modelling_2026/synthesis/pathway_summary.{py,json}, figure10_pathways.py; "
       "reassembly_2026/read_markers/README.md")

base = json.load(open(HERE / "r232_review_edits_2026-11-24.json"))
J = json.load(open(SYN / "pathway_summary.json"))
ab, pt, dbl = J["abundance_per_100GE"], J["pool_turnover_d"], J["doubling_h"]
att, cov, q = J["attribution_equal_growth"], J["coverage_at_1d_doubling"], J["q_required_methanogens_read"]
sh, mg = J["mcrA_groups"]["share_of_methanogen_reads"], J["mcrA_groups"]["reads"]
mtaB_arch = J["mcrA_groups"]["mtaB_archaeal_share"]
top = {p["key"]: p for p in J["pathways"]}


def pct(x, nd=0):
    return f"{x * 100:.{nd}f}"


def r(lohi, spec=".0f"):
    return f"{format(lohi[0], spec)}–{format(lohi[1], spec)}"


def rp(lohi):                                         # a share range in percent, at least one significant figure
    lo, hi = lohi
    return f"{max(round(lo * 100), 1)}–{round(hi * 100)}"


# ---- preconditions of the wording ---------------------------------------------------------------------------------------
assert sh["methylotrophic"] > 0.85 and 0.05 < sh["hydrogenotrophic"] < 0.12 and sh["acetoclastic"] < 0.01, sh
assert mg["acetoclastic"] <= 5 and 0.15 < mtaB_arch < 0.3
assert ab["mpn"]["reads_median"] > ab["betaine"]["reads_median"] > ab["methylamines"]["reads_median"] * 10
assert pt["methylamines"][1] < 4 and pt["betaine"][0] > 30 and pt["acetoclastic"][0] > 5, pt
a07, a33, aP, a07x10 = att["mpn_0.07 x1"], att["mpn_0.33 x1"], att["mpn_onlyP x1"], att["mpn_0.07 x10"]
assert a07["MPn"][1] < 0.1 and a07["methylotrophic"][0] > 0.8 and aP["MPn"][1] < 0.6, att
assert cov["methanogens"]["R2"][1] < 0.25 and cov["methanogens"]["R1"][1] > 1, cov
assert dbl["betaine"]["MAG"]["all"][0] > 24 and dbl["acetoclastic"]["reads_lower_bound"]["all"][1] < 1
def _and(xs):
    xs = list(xs)
    return xs[0] if len(xs) == 1 else ", ".join(xs[:-1]) + " and " + xs[-1]


tr = {k: _and(f"*{g}*" if g[:1].isupper() and g[1:].islower() and g.isalpha() else g for g, _ in v["top_reads"][:3])
      for k, v in top.items()}
assert top["mpn"]["top_reads"][0][0] == "Roseovarius" and top["hydrogenotrophic"]["top_reads"][0][0] == "Methanogenium"
assert [g for g, _ in top["methanol"]["top_reads"][:2]] == ["Methanohalophilus", "Methanococcoides"]

FACTS = dict(shares=sh, mcrA_reads=mg, mtaB_archaeal=mtaB_arch, pool_turnover_d=pt, attribution=att, coverage=cov,
             q_required=q)
NEW = {}


def put(key, paragraph, kind, old, new, rationale):
    NEW[key] = {"paragraph": paragraph, "kind": kind, "old": old, "new": new, "rationale": rationale, "sources": SRC,
                "origin": "synthesis_2026_11_25", "date": DATE}


put("p40_methods_synthesis", 40, "new_paragraph_after", "",
    "To compare the hypothesized methane sources (Figure 10), the abundance of each pathway was taken from its marker "
    "gene in all reads and from its MAGs (see Metagenomics); the *mcrA* reads were assigned to methylotrophic, "
    "hydrogenotrophic and acetoclastic methanogens by the genus of their best hit, pooled over the 12 samples of the "
    "unrestored cores (reads of the ANME-1-related Methanophagales were excluded), and only the archaeal reads of the "
    "methanol methyltransferase gene (*mtaB*) were counted. The standing porewater pool of each measured substrate was "
    "its NMR concentration times the porewater volume of each layer (porosity 0.74, from the bulk density of 0.7 g cm⁻³ "
    "and a particle density of 2.65 g cm⁻³), converted to methane with the stoichiometry of methanogenesis (methane per "
    "molecule: TMA 2.25, DMA 1.5, MMA 0.75, methanol 0.75, formate 0.25, acetate 1, and glycine betaine 2.25 through "
    "TMA); the pool divided by the measured emission is the turnover time at which it would supply the emission. To "
    "apportion the emission among the pathways, each guild’s methane per unit of growth (yield × read-based biomass) was "
    "compared at equal specific growth rates, and the methanogens’ share was divided among the methanogen groups by "
    "their *mcrA* reads.", "Methods for Figure 10 and the attribution")

put("p78_results_synthesis", 78, "new_paragraph_after", "",
    "Figure 10 sets the hypothesized methane sources side by side. Among the pathway genes in all reads, C–P lyase "
    f"(*phnJ*) was the most abundant ({r(ab['mpn']['reads'])} copies per 100 genome equivalents; "
    f"{r(ab['mpn']['MAGs'])} genomes per 100 in the MAGs), mostly of {tr['mpn']}, followed by the TMA-producing "
    f"betaine reductase (*grdI*; {r(ab['betaine']['reads'])}; mostly {tr['betaine']}) and choline and carnitine "
    f"genes (*cutC* + *cntA*; {r(ab['choline']['reads'])}), whereas methanogens were rare (*mcrA* "
    f"{r([ab['methylamines']['reads'][0] / sh['methylotrophic'], ab['methylamines']['reads'][1] / sh['methylotrophic']], '.1f')}). "
    f"Of the *mcrA* reads, {pct(sh['methylotrophic'])}% belonged to methylotrophic methanogens (*Methanolobus*, "
    "*Methanohalophilus*, *Methanococcoides*, and hydrogen-dependent Methanomassiliicoccales), "
    f"{pct(sh['hydrogenotrophic'])}% to hydrogenotrophic methanogens, mostly *Methanogenium*, and "
    f"{mg['acetoclastic']} reads to the acetoclastic *Methanothrix* and *Methanocrinis*. Only {pct(mtaB_arch)}% of the "
    "reads of the methanol methyltransferase gene *mtaB* were archaeal (mostly *Methanohalophilus* and "
    "*Methanococcoides*), and no MAG carried methanol-, hydrogen- or acetate-based methanogenesis. The measured porewater "
    "pools would supply the measured emissions if turned over every "
    f"{r(pt['methylamines'], '.1f')} days (mono-, di- and trimethylamine), {r(pt['methanol'], '.1f')} days (methanol), "
    f"{r([pt['hydrogenotrophic'][0] * 24, pt['hydrogenotrophic'][1] * 24])} hours (formate), "
    f"{r(pt['acetoclastic'])} days (acetate) and {r(pt['betaine'])} days (glycine betaine, if fully converted through "
    "TMA). The betaine reducers could supply the TMA for complete conversion at doubling times of "
    f"{dbl['betaine']['MAG']['all'][0]:.0f} hours or longer, whereas the acetoclastic methanogens, at their read-based "
    "abundance, would need doubling times of minutes.", "Results for Figure 10")

put("p98_conclusion", 98, "replace",
    "although methylphosphonate degradation by halophilic aerobic bacteria such as Halomonas appears to be the most "
    "corroborated conclusion from multiple perspectives.",
    "and methylphosphonate degradation by halophilic aerobic bacteria. Combining the read, MAG, substrate and model "
    "evidence attributes most of the methane to methylotrophic methanogenesis, with methylphosphonate degradation a "
    "secondary source (Figure 10).", "Conclusions opening agrees with the synthesis")

put("p99_synthesis", 99, "new_paragraph_after", "",
    "Taken together, these results attribute most of the measured methane to methylotrophic methanogenesis (Figure 10). "
    f"Methanogens that use methyl groups carry {pct(sh['methylotrophic'])}% of the methanogen *mcrA* reads, "
    f"hydrogenotrophic methanogens {pct(sh['hydrogenotrophic'])}%, and acetoclastic methanogens fewer than 1%, although "
    "acetate was abundant. Their substrates are mainly methylated amines (the dominant *Methanolobus* cannot use "
    "methanol, which only the less abundant *Methanohalophilus* and *Methanococcoides* may use): the porewater mono-, "
    "di- and trimethylamine "
    f"would supply the measured emissions if turned over every {r(pt['methylamines'], '.1f')} days, and the large "
    f"glycine betaine pool (enough for {r(pt['betaine'])} days of emissions if fully converted through TMA), degraded by "
    "abundant betaine-reducing Bacillota, together with choline and carnitine degradation, can resupply TMA without fast "
    "growth. Methylphosphonate degradation by C–P lyase bacteria (mostly *Roseovarius* and *Marivita* in situ, and "
    "*Halomonas* in culture) does contribute, but it is limited by their phosphorus demand: at equal specific growth "
    f"rates the MPn consumers would supply {rp(a07['MPn'])}% of the methane at an MPn:Pi uptake ratio of 0.07, "
    f"{rp(a33['MPn'])}% at 0.33, and {rp(aP['MPn'])}% only if methylphosphonate were their sole phosphorus source; even "
    f"growing ten times faster than the methanogens, they would supply {rp(a07x10['MPn'])}% at a ratio of 0.07. "
    "Growth-coupled production does not close the budget of the R2 cores: at a doubling time of one day, the "
    f"methanogens there would supply only {rp(cov['methanogens']['R2'])}% of the measured methane, so their methane "
    "must come largely from activity beyond growth "
    f"({r(q['R2'])} mmol CH₄ per gram of methanogen biomass per hour) or from more methanogen biomass than the DNA "
    "suggests. In the microcosms, only the TMA addition produced as much methane as the controls under anoxic "
    "conditions (Figure 4b), and methylphosphonic acid did not raise methane significantly above the controls under "
    "oxic conditions (Figure 4a). We therefore attribute the measured methane mainly to methylotrophic "
    "methanogenesis fed by the turnover of compatible solutes, with methylphosphonate degradation as a secondary, "
    "phosphorus-limited source and hydrogenotrophic methanogenesis as a minor one.",
    "the synthesis paragraph the user asked for (attribution among the hypotheses)")

put("p100_propose", 100, "replace",
    "We therefore propose methylphosphonate as the primary source of methane in the hypersaline unrestored salterns. "
    "The C-P lyase pathway in aerobic halotolerant bacteria appears to be the precise mechanism.",
    "We therefore propose that the methane of the hypersaline unrestored salterns comes mainly from methylotrophic "
    "methanogens that use methylated amines derived from compatible solutes, and secondarily from the C–P lyase "
    "pathway of aerobic halotolerant bacteria degrading methylphosphonate.",
    "the proposal agrees with the synthesis paragraph (replaces the P100 'Consider revising' comment)")
put("p100_comment", 100, "comment", "This project improves our understanding",
    "Review (2026-11-25): the proposal is rewritten to match the new synthesis paragraph above and Figure 10. The "
    "attribution rests on the read survey (all 2.08 billion reads), the MAGs, the porewater pools and the per-gram "
    "models; the shares assume equal specific growth rates of the guilds, with sensitivities given in the text. If "
    "the authors prefer to keep methylphosphonate as the primary source, the Abstract, the Introduction (P24) and the "
    "models section would need evidence that methylphosphonate is degraded far beyond assimilatory phosphorus demand.",
    "rationale for the rewritten proposal")

FIG10 = ("Figure 10. The hypothesized methane sources of the unrestored salterns. (a) Pathways: cleavage of "
         "methylphosphonate (MPn) by C–P lyase (blue); degradation of glycine betaine, choline and carnitine to TMA "
         "(orange), which feeds methylotrophic methanogenesis (dashed); and methanogenesis from methylamines, methanol, "
         "H₂/CO₂ or formate, and acetate (green). (b) Abundance per 100 genome equivalents in all reads (open: copies "
         "of *phnJ*, *grdI*, *cutC* + *cntA*, *mcrA* split by the methanogen group of the reads’ best hits, and the "
         "archaeal share of *mtaB*) and in the MAGs (filled: genomes of the guild’s species), median and range over the "
         "12 samples of the six unrestored cores. (c) The most abundant organisms: genus of the reads’ best hits (share "
         "of the pathway’s reads) and the most abundant MAG species. (d) Turnover each pathway alone would need to "
         "supply the measured methane of the six cores: turnover of the measured porewater pool (diamonds; methane per "
         "molecule: TMA 2.25, DMA 1.5, MMA 0.75, methanol 0.75, formate 0.25, acetate 1, glycine betaine 2.25 through "
         "TMA), and doubling time of the pathway’s organisms (squares; per-gram yields × measured biomass as in "
         "Figure 7, MAG-based filled, read-based open; MPn consumers at an MPn:Pi uptake ratio of 0.07 and with MPn as "
         "the only phosphorus source; betaine reducers supplying TMA for complete conversion). Hydrogenotrophic and "
         "acetoclastic methanogens have no model; their doubling times at the methylotrophic model’s yield are lower "
         "bounds (dotted; a lower growth yield allows slower growth). Methylphosphonic acid was detected but not "
         "quantified; choline, carnitine and H₂ were not measured.")

out = copy.deepcopy(base)
edits = []
dropped = 0
for e in out["edits"]:
    if e["paragraph"] == 100 and e["kind"] == "comment" and e["old"].startswith("We therefore propose"):
        dropped += 1                                    # anchor text is now a tracked deletion; see p100_comment
        continue
    if e["kind"] == "comment" and e["paragraph"] in (18, 24) and e.get("origin") == "modelling_2026_per_gram":
        e = dict(e, new=e["new"] + " See also the synthesis paragraph in the Conclusions and Figure 10 (2026-11-25).")
    edits.append(e)
assert dropped == 1
edits += list(NEW.values())
out["edits"] = edits
out["figure10"] = {"image": "reassembly_2026/modelling_2026/synthesis/figure10_pathways.png", "caption": FIG10,
                   "after": "Figure 9 caption"}
out["note"] = ("r232_review_edits_2026-11-24.json plus the 2026-11-25 synthesis of the methane sources "
               "(make_synthesis_edits_2026-11-25.py).")
out["facts_2026_11_25"] = json.loads(json.dumps(FACTS, default=float))
json.dump(out, open(HERE / "r232_review_edits_2026-11-25.json", "w"), indent=1, ensure_ascii=False)
print(f"wrote r232_review_edits_2026-11-25.json: {len(edits)} edits ({len(NEW)} new, {dropped} dropped)")
for k in ("p78_results_synthesis", "p99_synthesis"):
    print(f"\n{k}:\n{NEW[k]['new']}")
