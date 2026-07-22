#!/usr/bin/env python
"""
Assemble manuscript_corrected.md from the extracted manuscript text, applying the
Figure-7 number corrections (see numbers_to_reconcile.md).  Heading markers are
converted to Markdown; only the metabolic-modelling claims are edited -- every other
section is passed through verbatim.  Each edit is asserted so a silent text drift
fails loudly.
"""
from pathlib import Path

SCRATCH = Path("/tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns-salterns"
               "/c7120780-1bee-460b-a0c9-f3faaafbee99/scratchpad/manuscript_text.txt")
OUT = Path(__file__).resolve().parent / "manuscript_corrected.md"

raw = SCRATCH.read_text().split("\n")

# ---- convert heading markers to markdown ----------------------------------------
md = []
for ln in raw:
    if ln.startswith("[[FIGURE/IMAGE HERE]]"):
        continue
    for mark, pre in (("<Heading1> ", "# "), ("<Heading3> ", "## "), ("<Heading5> ", "**CAP**")):
        if ln.startswith(mark):
            body = ln[len(mark):]
            ln = (f"**{body}**" if pre == "**CAP**" else pre + body)
            break
    md.append(ln)
text = "\n".join(md)

# ---- targeted corrections (old -> new); assert each applies exactly once ---------
EDITS = [
 # Para 67 -- MePn modelling vs measured
 ("The simulated methane production of the MPn-consuming community members correlated "
  "remarkably well (R2 = 0.966) with the experimental methane emissions (Figure S8), and "
  "was most correlated when the phosphate source ratio was 0.07 (at least 1/14 of the total "
  "phosphate consumed in methylphosphonate versus inorganic phosphate).",
  "The simulated methane production of the MPn-consuming community members tracked the rank "
  "order of the experimental methane emissions (Figure 7a; Spearman rho = 0.60), and was most "
  "correlated when the phosphate source ratio was 0.07 (at least 1/14 of the total phosphate "
  "consumed in methylphosphonate versus inorganic phosphate). Because this comparison rests on "
  "only the six unrestored cores, the association is suggestive rather than statistically "
  "significant (p = 0.21), and the coefficient of determination is sensitive to normalization: "
  "R2 = 0.37 when both axes are expressed as raw areal emission and R2 = 0.89 when both are "
  "divided by DNA density as in earlier drafts."),
 # Para 69 -- ratio sensitivity + Rhodobacteraceae
 ("and yield an almost identical fit to the best fitting regression (R2 = 0.968). It is worth "
  "noting that unconstraining the phosphate source ratio produced almost an identical fit "
  "(R2 = 0.970) and still resulted in methylphosphonate consumption by several marine organisms "
  "in the Rhodobacteraceae family: two species in the Roseovarius genera and one species in the "
  "Albimonas genera.",
  "and yield a comparable fit. Importantly, unconstraining the phosphate source ratio entirely -- "
  "introducing no tuned parameter -- produced a comparable fit and still resulted in "
  "methylphosphonate consumption dominated by marine organisms in the Rhodobacteraceae family: "
  "of the modelled methylphosphonate-derived methane, 77% originated from Roseovarius and 11% "
  "from Albimonas, with the remainder from other Rhodobacteraceae genera (Figure 7c)."),
 # Para 70 -- betaine modelling
 ("The simulated methane production of the betaine-consuming community members correlated "
  "slightly less well (R2 = 0.862) with the experimental methane emissions (Figure 7), with the "
  "line-of-best-fit",
  "By contrast, the simulated methane production of the betaine-consuming community members "
  "showed essentially no relationship with the experimental methane emissions (Figure 7a; "
  "Spearman rho = -0.09, R2 = 0.10 on the areal basis), with the line-of-best-fit"),
 # Para 73 -- abundance correlations
 ("The total experimental abundances of methylphosphonate-consuming organisms was significantly "
  "positively correlated with methane emissions (R2 = 0.974 for all samples and R2 = 0.966 for "
  "unrestored samples), and more correlated than the betaine-reducing organisms (R2 = 0.852 for "
  "all samples and R2 = 0.798 for unrestored samples). This experimentally supports the evidence "
  "from metabolic modeling that methylphosphonate is much more strongly associated with methane "
  "emissions than metabolism of the compatible solute betaine.",
  "Across all twelve sediment cores, the summed relative abundance of methylphosphonate-consuming "
  "organisms was significantly and positively correlated with methane emissions (Figure 7b; "
  "Spearman rho = 0.73, p = 0.007). The abundance of betaine-reducing organisms was similarly "
  "correlated (rho = 0.70, p = 0.011); the two associations are of comparable strength once all "
  "cores are considered, although methylphosphonate consumers explain more of the variance "
  "(R2 = 0.75 versus 0.38). This supports the metabolic-modeling evidence that methylphosphonate "
  "metabolism is at least as strongly -- and by variance explained, more strongly -- associated "
  "with methane emissions as metabolism of the compatible solute betaine. Restricting the "
  "analysis to the six unrestored cores alone renders neither association significant, "
  "underscoring the value of the full restoration gradient."),
 # Para 90 -- discussion
 ("MPn-consuming organisms were much more associated with methane emissions than betaine-consuming "
  "organisms (R2 = 0.798). This was mechanistically corroborated by metabolic models that "
  "reproduced the experimental methane emissions with high accuracy (R2 = 0.978), and demonstrate "
  "that methane is generated stoichiometrically from methylphosphonate.",
  "MPn-consuming organisms were at least as strongly associated with methane emissions as "
  "betaine-consuming organisms, and explained more of its variance across the full set of twelve "
  "cores (R2 = 0.75 versus 0.38; Spearman rho = 0.73 versus 0.70). This was mechanistically "
  "corroborated by metabolic models that reproduced the rank order of the experimental methane "
  "emissions (Spearman rho = 0.60), and demonstrate that methane is generated stoichiometrically "
  "from methylphosphonate. Because the modelling comparison rests on six unrestored cores, we "
  "treat it as consistency-checking evidence that complements the independent culturing and "
  "substrate-addition experiments rather than as a standalone statistical validation."),
]

for old, new in EDITS:
    assert text.count(old) == 1, f"edit did not match uniquely ({text.count(old)}x): {old[:70]!r}"
    text = text.replace(old, new)

# ---- add a normalization/statistics note into the Analyses methods --------------
anchor = ("The metabolomics concentrations were also correlated with the relative abundances to "
          "empirically inspect how the organisms interacted with the substrates and syntrophic compounds.")
assert text.count(anchor) == 1
text = text.replace(anchor, anchor + " Correlations between modelled and measured methane, and "
    "between consumer abundance and methane, were quantified with Pearson's R2 and Spearman's rho "
    "(with the exact p-value reported given the small sample size); modelled and measured methane "
    "were compared on a consistent areal-emission basis (umol CH4 m-2 d-1) to permit a 1:1 "
    "reference line, and abundance-methane correlations used all twelve depth-merged cores rather "
    "than the unrestored subset alone.")

OUT.write_text(text.rstrip() + "\n")
print(f"wrote {OUT}  ({len(text.split())} words, {text.count(chr(10))+1} lines)")
print("applied", len(EDITS) + 1, "edits; all matched uniquely.")
