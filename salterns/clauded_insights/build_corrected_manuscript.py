#!/usr/bin/env python
"""
Assemble manuscript_corrected.md from the text of "Metabolic Manuscript.docx" (2 July 2026),
applying the Figure-7 number corrections (see numbers_to_reconcile.md) and the three draft
errors found in the R232 review (Table 3, P35 MAG quality, P69 figure citation).  Heading markers
are converted to Markdown; only the metabolic-modelling claims, the P35 MAG-quality phrase and
Table 3 (caption + values) are edited -- every other section is passed through verbatim.  Each
edit is asserted so a silent text drift fails loudly.

Input text.  The July build read a one-off text dump from a session scratchpad that no longer
exists.  extract_text() regenerates that dump from the DOCX (accepted-changes view: text inside
w:del is skipped, w:ins text kept; tabs and line breaks dropped; one line per non-empty paragraph,
with <HeadingN> markers for Heading 1/3/5 and a [[FIGURE/IMAGE HERE]] line before any paragraph
holding a drawing; table cells one per line, merged cells once, empty cells skipped).  With the
July edits only, it reproduces the July manuscript_corrected.md byte for byte (checked 2026-10-26).

P69 (Figure S8 -> Figure 7a) is already handled by the July "Para 67" edit below.
Table 3 values come from "Table 3 - corrected values (2020 GTDB families).csv"
(table3_2020_corrected.py; 2020 GTDB families, all 310 MAGs incl. the 30 archaea).

Usage: ~/Documents/py_venv/bin/python build_corrected_manuscript.py
"""
import csv
from pathlib import Path

from docx import Document
from docx.oxml.ns import qn
from docx.table import Table
from docx.text.paragraph import Paragraph

S = Path(__file__).resolve().parents[1]
SRC = S / "Metabolic Manuscript.docx"
TABLE3_CSV = S / "Table 3 - corrected values (2020 GTDB families).csv"
OUT = Path(__file__).resolve().parent / "manuscript_corrected.md"

_MARK = {"Heading 1": "<Heading1> ", "Heading 3": "<Heading3> ", "Heading 5": "<Heading5> "}
_W_T, _W_TAB, _W_BR, _W_CR, _W_DEL = (qn(t) for t in ("w:t", "w:tab", "w:br", "w:cr", "w:del"))


def _accepted_text(el):
    """Text of `el` with tracked deletions skipped (insertions kept); tabs and breaks dropped."""
    out = []
    for e in el.iter(_W_T, _W_TAB, _W_BR, _W_CR):
        if e.tag != _W_T:
            continue
        a = e.getparent()
        while a is not None and a is not el and a.tag != _W_DEL:
            a = a.getparent()
        if a is None or a is el:
            out.append(e.text or "")
    return "".join(out)


def extract_text(path):
    doc = Document(str(path))
    lines = []
    for child in doc.element.body.iterchildren():
        if child.tag == qn("w:p"):
            p = Paragraph(child, doc._body)
            if child.findall(".//" + qn("w:drawing")) or child.findall(".//" + qn("w:pict")):
                lines.append("[[FIGURE/IMAGE HERE]]")
            t = _accepted_text(child)
            if t.strip():
                lines.append(_MARK.get(p.style.name if p.style is not None else None, "") + t)
        elif child.tag == qn("w:tbl"):
            for row in Table(child, doc._body).rows:
                seen = set()
                for c in row.cells:
                    if id(c._tc) in seen:
                        continue
                    seen.add(id(c._tc))
                    t = "\n".join(_accepted_text(pp) for pp in c._tc.iter(qn("w:p")))
                    if t.strip():
                        lines.append(t)
    return "\n".join(lines)


raw = extract_text(SRC).split("\n")

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

# ---- the three draft errors from the R232 review (2026-10-21) ---------------------
# P69 "(Figure S8)" -> "(Figure 7a)" is already part of the Para 67 edit above.
# P35: 35 of the 310 MAGs exceed 10% CheckM contamination, so they are not all "medium and
# high quality".  (The R232 review cites Table S5 here; these July copies have no Table S5.)
EDITS.append(
 ("dereplicated medium and high quality metagenome-assembled genomes (MAGs)",
  "310 dereplicated metagenome-assembled genomes (MAGs; 275 meet the MIMAG medium-quality criteria)"))

# Table 3: the old values came from correlational/abundances.ipynb cell 1, which pooled the R2A
# cores into "unrestored", dropped the 30 archaea with the unbinned reads, and was copied with
# the Restored and Reference columns swapped.  Corrected values: 2020 GTDB families, all 310 MAGs.
TABLE3_HEADER = ["GTDB-Tk Family", "Unrestored total", "Restored total", " Reference total"]
TABLE3_CAPTION = ("**Table 3. The percent of total biomass that each family comprises in the unrestored, "
                  "restored, and reference (natural) wetlands.  These percentages were computed after "
                  "removing the unbinned reads, which constituted {}, {}, and {} of the total abundance in "
                  "each sample type, respectively. **")
TABLE3_OLD = [
    ["Rhodobacteraceae", "0.205579", "0.023704", "0.008131"],
    ["Balneolaceae", "0.112063", "0.001370", "0.085078"],
    ["Wenzhouxiangellaceae", "0.096556", "0.001906", "0.000975"],
    ["undetermined", "0.081655", "0.169926", "0.165546"],
    ["Saprospiraceae", "0.075393", "0.000535", "0.000266"],
    ["Woeseiaceae", "0.032407", "0.019029", "0.005789"],
    ["UBA12077", "0.031803", "0.000287", "0.000523"],
    ["Sediminispirochaetaceae", "0.027220", "0.000319", "0.000519"],
    ["Flavobacteriaceae", "0.025992", "0.063411", "0.019261"],
    ["Oleiphilaceae", "0.023174", "0.001592", "0.000722"],
]


def table3_corrected():
    """(10 rows [family, unrestored, restored, reference], 3 unbinned shares) from the CSV."""
    with open(TABLE3_CSV, newline="") as fh:
        rows = list(csv.DictReader(fh))
    cols = ("Unrestored (R1, R2)", "Restored (SF2)", "Reference (R2A)")
    top = [r for r in rows if r["In corrected top 10"] == "yes"]
    assert len(top) == 10, len(top)
    unb = [r for r in rows if r["Family (2020 GTDB)"].startswith("Unbinned share")]
    assert len(unb) == 1
    return ([[r["Family (2020 GTDB)"], *(f"{float(r[c]):.6f}" for c in cols)] for r in top],
            [f"{float(unb[0][c]) * 100:.0f}%" for c in cols])


T3_ROWS, T3_UNBINNED = table3_corrected()
EDITS.append(
 ("\n".join([TABLE3_CAPTION.format("68%", "89%", "89%"), *TABLE3_HEADER, *sum(TABLE3_OLD, [])]),
  "\n".join([TABLE3_CAPTION.format(*T3_UNBINNED), *TABLE3_HEADER, *sum(T3_ROWS, [])])))

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
