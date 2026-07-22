#!/usr/bin/env python
"""
Produce a copy of the manuscript DOCX with REAL Word tracked changes (w:ins / w:del)
implementing the review edits:
  * 5 reworded metabolic-modelling paragraphs (numbers corrected; text taken verbatim
    from manuscript_corrected.md) -> old paragraph content marked deleted, corrected
    text marked inserted;
  * new inserted paragraphs (all wrapped as insertions): geochem/gradient result,
    specificity caveat, modelling-as-attribution reframe, and a proposed-Figure-1 note.
Open in Word with Review > All Markup to see/accept the changes.
"""
import copy
from pathlib import Path
from docx import Document
from docx.text.paragraph import Paragraph
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

SRC = Path(__file__).resolve().parents[1] / "Metabolic Manuscript.docx"
CORR = Path(__file__).resolve().parent / "manuscript_corrected.md"
OUT = Path(__file__).resolve().parent / "Metabolic Manuscript (tracked changes).docx"
AUTHOR, DATE = "Claude (review)", "2026-07-16T00:00:00Z"
_id = iter(range(1000, 9000))


def _mkrun(text, italic=False):
    r = OxmlElement("w:r")
    if italic:
        rpr = OxmlElement("w:rPr"); rpr.append(OxmlElement("w:i")); r.append(rpr)
    t = OxmlElement("w:t"); t.set(qn("xml:space"), "preserve"); t.text = text; r.append(t)
    return r


def _wrap_ins(*runs):
    ins = OxmlElement("w:ins")
    ins.set(qn("w:id"), str(next(_id))); ins.set(qn("w:author"), AUTHOR); ins.set(qn("w:date"), DATE)
    for r in runs:
        ins.append(r)
    return ins


def replace_paragraph(paragraph, new_text):
    """Mark ALL of the paragraph's content (runs + hyperlinked citations) deleted and
    insert new_text (tracked). Paragraphs here contain no pre-existing revisions."""
    p = paragraph._p
    pPr = p.find(qn("w:pPr"))
    content = [c for c in list(p) if c is not pPr]
    if content:
        wdel = OxmlElement("w:del")
        wdel.set(qn("w:id"), str(next(_id))); wdel.set(qn("w:author"), AUTHOR); wdel.set(qn("w:date"), DATE)
        content[0].addprevious(wdel)
        for c in content:
            for t in c.iter(qn("w:t")):        # recurse into hyperlinks etc.
                t.tag = qn("w:delText")
            wdel.append(c)                     # move the whole child into the deletion
    p.append(_wrap_ins(_mkrun(new_text)))


def insert_after(paragraph, text, italic=False):
    """Insert a fully-tracked new paragraph after `paragraph`; returns the new Paragraph."""
    new_p = OxmlElement("w:p")
    src = paragraph._p.find(qn("w:pPr"))
    pPr = copy.deepcopy(src) if src is not None else OxmlElement("w:pPr")
    rPr = pPr.find(qn("w:rPr"))
    if rPr is None:
        rPr = OxmlElement("w:rPr"); pPr.append(rPr)
    ins_mark = OxmlElement("w:ins")
    ins_mark.set(qn("w:id"), str(next(_id))); ins_mark.set(qn("w:author"), AUTHOR); ins_mark.set(qn("w:date"), DATE)
    rPr.insert(0, ins_mark)                    # mark the paragraph-end as inserted
    new_p.append(pPr)
    new_p.append(_wrap_ins(_mkrun(text, italic)))
    paragraph._p.addnext(new_p)
    return Paragraph(new_p, paragraph._parent)


# ---- corrected paragraph text (verbatim from manuscript_corrected.md) --------------
corr_lines = [ln for ln in CORR.read_text().split("\n") if ln.strip()]
def corrected(anchor):
    hits = [ln for ln in corr_lines if anchor in ln]
    assert len(hits) == 1, f"corrected md anchor not unique ({len(hits)}): {anchor!r}"
    return hits[0].strip()

REPLACE = [   # anchor present in BOTH the original DOCX paragraph and its corrected md line
    "Soil DNA density was significantly positively correlated",
    "Roseovarius has been experimentally observed to consume",
    "betaine-consuming community members",
    "methylphosphonate-consuming organisms was significantly",
    "MPn-consuming organisms were",
]

# ---- new inserted paragraphs (a, c, b, Fig1) --------------------------------------
GEOCHEM = ("We further asked whether the environmental gradient itself, rather than any "
    "specific taxon, structures these emissions. Across the twelve cores, sediment sulfate "
    "(Spearman rho = 0.97, p < 0.001), chloride and salinity (rho = 0.94 and 0.79), and the "
    "inorganic N:P ratio (rho = 0.73, p = 0.007) were each more strongly correlated with "
    "methane emission than any biological variable, and a rank-based variance partition "
    "attributed only about 2% of the methane variance uniquely to methylphosphonate-consumer "
    "abundance beyond the sulfate gradient (42% unique to the gradient, 51% shared). This is "
    "consistent with a system organised by a single salinity/sulfate axis along which methane, "
    "metabolites, and community composition co-vary; the elevated N:P ratios (severe phosphate "
    "limitation) are themselves consistent with the methylphosphonate-scavenging hypothesis "
    "(proposed geochemistry figure).")
SPECIFICITY = ("Because the summed abundance of even a random set of MAGs of the same size "
    "correlates with methane about as well as the methylphosphonate- or betaine-consumer sets "
    "(37% and 45% of random subsets, respectively, matched or exceeded the observed "
    "coefficients), and because the association weakens after controlling for total DNA density, "
    "these cross-sample abundance correlations should be read as consistent with, rather than as "
    "specific evidence for, a particular methanogenic guild. The manipulative substrate-addition "
    "and culturing experiments provide the organism-specificity that correlation across the "
    "salinity gradient cannot.")
ATTRIBUTION = ("We therefore interpret the metabolic modelling as a demonstration of mechanistic "
    "sufficiency and taxonomic attribution -- methylphosphonate consumption can generate methane "
    "of the observed magnitude, and this capacity is concentrated in a few Rhodobacteraceae "
    "(Roseovarius, 77%; Albimonas, 11%; proposed Figure 7c) -- rather than as an independent "
    "statistical validation, which the six-core sample and the shared environmental gradient "
    "cannot support.")
FIG1 = ("[Proposed Figure 1: a conceptual overview of the study system -- the restoration / "
    "salinity-sulfate gradient, the two candidate methane pathways (aerobic methylphosphonate "
    "degradation and anaerobic betaine->TMA methylotrophic methanogenesis), and the five "
    "converging lines of evidence -- would orient the reader; a draft is provided in the "
    "clauded_insights folder.]")

doc = Document(str(SRC))

def find(anchor):
    hits = [p for p in doc.paragraphs if anchor in p.text]
    assert len(hits) == 1, f"DOCX anchor not unique ({len(hits)}): {anchor!r}"
    return hits[0]

# resolve ALL anchor paragraphs first (paragraph.text goes blank once wrapped in ins/del)
P = {a: find(a) for a in REPLACE}
P["_fig1"] = find("phnGHILMJ genes for the carbon-phosphorus lyase")

# 1) replace the 5 modelling paragraphs
for anchor in REPLACE:
    replace_paragraph(P[anchor], corrected(anchor))

# 2) insertions (using the stored paragraph objects)
a = insert_after(P["methylphosphonate-consuming organisms was significantly"], GEOCHEM)
insert_after(a, SPECIFICITY)
insert_after(P["MPn-consuming organisms were"], ATTRIBUTION)
insert_after(P["_fig1"], FIG1, italic=True)

doc.save(str(OUT))
print("saved", OUT.name)
print("tracked edits: 5 paragraph rewrites + 4 inserted paragraphs")
