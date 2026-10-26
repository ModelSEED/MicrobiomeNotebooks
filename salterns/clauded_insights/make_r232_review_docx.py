#!/usr/bin/env python
"""
R232 review of the saltern manuscript: a copy of ../Metabolic Manuscript.docx with real Word tracked
changes (w:ins / w:del) and review comments that bring in the 2026 MAG re-analysis:
CheckM re-run, GTDB-Tk 2.8.0 / R232 taxonomy, gene-based functional guilds, guild-methane correlations
with a random-set comparison, the methanogen-methane correlation, and the R232 co-occurrence network
(../correlational/r232/). "Reject all" in Word restores the original text, tables and figures.

Input is the manuscript after the docx skill's merge_runs.py (adjacent identically formatted runs
merged, text unchanged), so that sentences can be located inside single runs:
    python <docx-skill>/scripts/merge_runs.py "../Metabolic Manuscript.docx" -o merged.docx
    ~/Documents/py_venv/bin/python make_r232_review_docx.py merged.docx
Output: ../Metabolic Manuscript - R232 review (tracked changes).docx

The 2026-10-26 revision adds the tracked edits listed in r232_review_edits_2026-10-26.json:
- the modelling text and numbers on the areal basis (modelling_numbers_r232.py);
- betaine models restricted to the nine grdH carriers (../correlational/betaine_grdH_r232/);
- the R232 reclassification of the two Halomonas isolates (../MAG_QC/gtdbtk_isolates_r232/);
- smaller corrections.

It also replaces Figure 7 with figure7_r232.png (figure7_r232.py).
"""
import copy
import json
import re
import sys
from pathlib import Path

import pandas as pd
from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches
from docx.text.run import Run

HERE = Path(__file__).resolve().parent
S = HERE.parent
R = S / "correlational" / "r232"
OUT = S / "Metabolic Manuscript - R232 review (tracked changes).docx"
AUTHOR, INITIALS, DATE = "Claude (R232 review)", "CR", "2026-10-26T00:00:00Z"
EDITS = HERE / "r232_review_edits_2026-10-26.json"
_ids = iter(range(7000, 999999))
CONTENT = {qn(t) for t in ("w:t", "w:tab", "w:br", "w:drawing", "w:pict", "w:sym", "w:footnoteReference")}
MARKUP = re.compile(r"(\*[^*]+\*|\^[^^]+\^|~[^~]+~)")


# ---- low-level XML helpers --------------------------------------------------------------------------
def el(tag, **attrs):
    e = OxmlElement(tag)
    for k, v in attrs.items():
        e.set(qn(k), str(v))
    return e


def stamp(e):
    e.set(qn("w:id"), str(next(_ids))); e.set(qn("w:author"), AUTHOR); e.set(qn("w:date"), DATE)
    return e


def make_rpr(tpl, bold=False, italic=False, valign=None, size=None):
    """Schema-ordered rPr: rFonts, b, bCs, i, iCs, color, sz, szCs, vertAlign."""
    rpr = el("w:rPr")
    if tpl is not None and tpl.find(qn("w:rFonts")) is not None:
        rpr.append(copy.deepcopy(tpl.find(qn("w:rFonts"))))
    if bold:
        rpr.append(el("w:b")); rpr.append(el("w:bCs"))
    if italic:
        rpr.append(el("w:i")); rpr.append(el("w:iCs"))
    if tpl is not None and tpl.find(qn("w:color")) is not None:
        rpr.append(copy.deepcopy(tpl.find(qn("w:color"))))
    sz = size or (tpl.find(qn("w:sz")).get(qn("w:val")) if tpl is not None and tpl.find(qn("w:sz")) is not None else None)
    if sz:
        rpr.append(el("w:sz", **{"w:val": sz})); rpr.append(el("w:szCs", **{"w:val": sz}))
    if valign:
        rpr.append(el("w:vertAlign", **{"w:val": valign}))
    return rpr


def runs_from(markup, tpl, bold=False, size=None, lead_tab=False):
    """*italic*, ^superscript^, ~subscript~ -> list of w:r (lead_tab: start with a tab, like the body text)."""
    out = []
    if lead_tab:
        r = el("w:r"); r.append(make_rpr(tpl, bold=bold, size=size)); r.append(el("w:tab")); out.append(r)
    for tok in MARKUP.split(markup):
        if not tok:
            continue
        kw = {}
        if tok[0] == tok[-1] == "*":
            tok, kw = tok[1:-1], {"italic": True}
        elif tok[0] == tok[-1] == "^":
            tok, kw = tok[1:-1], {"valign": "superscript"}
        elif tok[0] == tok[-1] == "~":
            tok, kw = tok[1:-1], {"valign": "subscript"}
        r = el("w:r")
        r.append(make_rpr(tpl, bold=bold, size=size, **kw))
        t = el("w:t"); t.set(qn("xml:space"), "preserve"); t.text = tok
        r.append(t); out.append(r)
    return out


def ins(runs):
    w = stamp(el("w:ins"))
    for r in runs:
        w.append(r)
    return w


def mark_pmark(p, kind):
    """Track the paragraph mark itself (w:ins / w:del go FIRST in pPr/rPr)."""
    pPr = p.find(qn("w:pPr"))
    if pPr is None:
        pPr = el("w:pPr"); p.insert(0, pPr)
    rPr = pPr.find(qn("w:rPr"))
    if rPr is None:
        rPr = el("w:rPr")
        tail = [c for c in pPr if c.tag in (qn("w:sectPr"), qn("w:pPrChange"))]
        (tail[0].addprevious(rPr) if tail else pPr.append(rPr))
    rPr.insert(0, stamp(el(kind)))


def clean_pPr(pPr):
    pPr = copy.deepcopy(pPr)
    for e in list(pPr.iter(qn("w:ins"))) + list(pPr.iter(qn("w:del"))):
        e.getparent().remove(e)
    return pPr


def is_content(r):
    return r.tag == qn("w:r") and any(c.tag in CONTENT for c in r)


def wrap_del(r):
    d = stamp(el("w:del")); r.addprevious(d); d.append(r)
    for t in r.iter(qn("w:t")):
        t.tag = qn("w:delText")
    return d


def delete_contents(p):
    for child in list(p):
        if is_content(child):
            wrap_del(child)
        elif child.tag == qn("w:hyperlink"):
            for r in [r for r in child if is_content(r)]:
                wrap_del(r)


def first_rpr(p):
    for r in p.iter(qn("w:r")):
        if r.find(qn("w:t")) is not None or r.find(qn("w:delText")) is not None:
            return r.find(qn("w:rPr"))
    return None


# ---- paragraph-level edits ------------------------------------------------------------------------------
def replace_paragraph(p, markup, tpl=None, lead_tab=False):
    tpl = tpl if tpl is not None else first_rpr(p)
    delete_contents(p)
    w = ins(runs_from(markup, tpl, lead_tab=lead_tab))
    p.append(w)
    return w


def delete_paragraph(p):
    delete_contents(p); mark_pmark(p, "w:del")


def new_par_after(anchor, pPr_src, markup, tpl, bold=False, size=None, lead_tab=False):
    p = el("w:p")
    if pPr_src is not None:
        p.append(clean_pPr(pPr_src))
    mark_pmark(p, "w:ins")
    w = ins(runs_from(markup, tpl, bold=bold, size=size, lead_tab=lead_tab)); p.append(w)
    anchor.addnext(p)
    return p, w


def _text_runs(p):
    return [r for r in p if r.tag == qn("w:r") and r.find(qn("w:t")) is not None]


def _split(r, k):
    """Split a single-w:t run at character k; returns the right-hand run."""
    t = r.find(qn("w:t")); text = t.text
    right = copy.deepcopy(r)
    t.text = text[:k]; t.set(qn("xml:space"), "preserve")
    rt = right.find(qn("w:t")); rt.text = text[k:]; rt.set(qn("xml:space"), "preserve")
    # other content (e.g. a leading w:tab) stays on its own side of the w:t, not in both halves
    kids, rkids = list(r), list(right)
    ti = kids.index(t)
    for c in kids[ti + 1:]:
        r.remove(c)
    for c in rkids[:ti]:
        if c.tag != qn("w:rPr"):
            right.remove(c)
    r.addnext(right)
    return right


def _locate(p, s):
    runs = _text_runs(p)
    assert all(len(r.findall(qn("w:t"))) == 1 for r in runs)
    full = "".join(r.find(qn("w:t")).text or "" for r in runs)
    assert full.count(s) == 1, f"anchor occurs {full.count(s)}x: {s[:60]!r}"
    return runs, full.index(s)


def _cut_at(p, offset):
    """Ensure a run boundary at character offset (over direct text runs); returns run starting there."""
    pos = 0
    for r in _text_runs(p):
        n = len(r.find(qn("w:t")).text or "")
        if pos == offset:
            return r
        if pos < offset < pos + n:
            return _split(r, offset - pos)
        pos += n
    return None


def tracked_replace(p, old, new_markup):
    _, start = _locate(p, old)
    first = _cut_at(p, start)
    _cut_at(p, start + len(old))
    covered, pos = [], 0
    for r in _text_runs(p):
        n = len(r.find(qn("w:t")).text or "")
        if start <= pos and pos + n <= start + len(old) and n:
            covered.append(r)
        pos += n
    tpl = first.find(qn("w:rPr"))
    for r in covered:                       # a leading w:tab is the paragraph indent: keep it, undeleted, in front
        kids = list(r); t = r.find(qn("w:t"))
        tabs = [c for c in kids[:kids.index(t)] if c.tag == qn("w:tab")]
        if tabs:
            keep = el("w:r")
            if r.find(qn("w:rPr")) is not None:
                keep.append(copy.deepcopy(r.find(qn("w:rPr"))))
            for c in tabs:
                r.remove(c); keep.append(c)
            r.addprevious(keep)
    dels = [wrap_del(r) for r in covered]
    w = ins(runs_from(new_markup, tpl)) if new_markup else None
    if w is not None:
        dels[-1].addnext(w)
    return w


def tracked_insert_before(p, anchor_text, markup):
    _, start = _locate(p, anchor_text)
    r = _cut_at(p, start)
    w = ins(runs_from(markup, r.find(qn("w:rPr"))))
    r.addprevious(w)
    return w


# ---- figures and tables ---------------------------------------------------------------------------------
def picture_after(doc, anchor, path, width_in, pPr_src):
    para = doc.add_paragraph()
    run = para.add_run(); run.add_picture(str(path), width=Inches(width_in))
    p = para._p
    old = p.find(qn("w:pPr"))
    if old is not None:
        p.remove(old)
    p.insert(0, clean_pPr(pPr_src))
    mark_pmark(p, "w:ins")
    w = stamp(el("w:ins")); run._r.addprevious(w); w.append(run._r)
    anchor.addnext(p)
    return p


def build_table(tblPr_src, tcMar_src, header, rows, widths, tpl, size=20, center_from=1, align=None):
    tbl = el("w:tbl")
    tblPr = copy.deepcopy(tblPr_src)
    tblPr.find(qn("w:tblW")).set(qn("w:w"), str(sum(widths)))
    tbl.append(tblPr)
    grid = el("w:tblGrid")
    for w in widths:
        grid.append(el("w:gridCol", **{"w:w": w}))
    tbl.append(grid)
    for ri, cells in enumerate([header] + rows):
        tr = el("w:tr"); trPr = el("w:trPr")
        if ri == 0:
            trPr.append(el("w:tblHeader"))
        trPr.append(stamp(el("w:ins")))            # w:ins goes LAST in trPr
        tr.append(trPr)
        for ci, (txt, w) in enumerate(zip(cells, widths)):
            tc = el("w:tc"); tcPr = el("w:tcPr")
            tcPr.append(el("w:tcW", **{"w:w": w, "w:type": "dxa"}))
            if tcMar_src is not None:
                tcPr.append(copy.deepcopy(tcMar_src))
            tc.append(tcPr)
            p = el("w:p"); pPr = el("w:pPr")
            pPr.append(el("w:spacing", **{"w:line": 240, "w:lineRule": "auto"}))
            jc = align[ci] if (align and ri > 0) else ("center" if (ri == 0 or ci >= center_from) else "left")
            pPr.append(el("w:jc", **{"w:val": jc}))
            p.append(pPr); mark_pmark(p, "w:ins")
            p.append(ins(runs_from(str(txt), tpl, bold=(ri == 0), size=size)))
            tc.append(p); tr.append(tc)
        tbl.append(tr)
    return tbl


def delete_table(tbl):
    for tr in tbl.findall(qn("w:tr")):
        trPr = tr.find(qn("w:trPr"))
        if trPr is None:
            trPr = el("w:trPr"); tr.insert(0, trPr)
        trPr.append(stamp(el("w:del")))
        for p in tr.iter(qn("w:p")):
            delete_contents(p)


def comment(doc, para_obj, run_el, text):
    doc.add_comment(Run(run_el, para_obj), text=text, author=AUTHOR, initials=INITIALS)


M_NS = "{http://schemas.openxmlformats.org/officeDocument/2006/math}"


def tracked_math_swap(p, replacements, drop_fraction=None):
    """Tracked replacement of a paragraph's whole Office Math object: the old m:oMath goes into a w:del and an edited
    copy into a w:ins (both allowed as run-level content). Swapping the whole object, rather than wrapping single math
    runs, keeps the equation intact in viewers that ignore revisions inside math. `replacements` maps exact m:t texts
    to new texts; `drop_fraction` = (numerator, denominator) text of an m:f to remove from the new copy."""
    maths = p.findall(M_NS + "oMath")
    assert len(maths) == 1, f"{len(maths)} m:oMath in paragraph"
    old = maths[0]
    new = copy.deepcopy(old)
    for t in new.iter(M_NS + "t"):
        if t.text in replacements:
            t.text = replacements.pop(t.text)
    assert not replacements, f"math texts not found: {list(replacements)}"
    if drop_fraction:
        fr = [f for f in new.iter(M_NS + "f")
              if "".join(t.text for t in f.find(M_NS + "num").iter(M_NS + "t")) == drop_fraction[0]
              and "".join(t.text for t in f.find(M_NS + "den").iter(M_NS + "t")) == drop_fraction[1]]
        assert len(fr) == 1, f"fraction {drop_fraction} occurs {len(fr)}x"
        fr[0].getparent().remove(fr[0])
    d = stamp(el("w:del")); old.addprevious(d); d.append(old)
    w = stamp(el("w:ins")); w.append(new); d.addnext(w)


def apply_edits(P, pending_comments, path):
    """Tracked edits from a JSON list; 'paragraph' is the paragraph index in the original draft.
    Text edits go first, then comments, whose anchors must lie in text that is still undeleted."""
    for e in sorted(json.load(open(path))["edits"], key=lambda e: e["kind"] == "comment"):
        p = P[e["paragraph"]]
        if e["kind"] == "replace":
            tracked_replace(p._p, e["old"], e["new"])
        elif e["kind"] == "insert_before":
            tracked_insert_before(p._p, e["old"], e["new"])
        elif e["kind"] == "replace_paragraph":
            replace_paragraph(p._p, e["new"], lead_tab=p.text.startswith("\t"))
        else:
            _, start = _locate(p._p, e["old"])
            pending_comments.append((p, _cut_at(p._p, start), e["new"]))


# ---- content -------------------------------------------------------------------------------------------------
LATIN = re.compile(r"^[A-Z][a-z]+$")


def it_id(i):
    """Italicise Latin genus names in iterativeIDs (e.g. *Roseovarius*.3), not placeholder names."""
    base, n = i.rsplit(".", 1)
    return f"*{base}*.{n}" if LATIN.match(base) else i


def compress(ids):
    by = {}
    for i in ids:
        b, n = i.rsplit(".", 1); by.setdefault(b, []).append(int(n))
    parts = []
    for b, ns in by.items():
        ns.sort(); spans, s = [], ns[0]
        for a, c in zip(ns, ns[1:] + [None]):
            if c != a + 1:
                spans.append(f"{s}" if s == a else f"{s}–{a}"); s = c
        name = f"*{b}*" if LATIN.match(b) else b
        parts.append(f"{name}.{', '.join(spans)}")
    return "; ".join(parts)


def num(x, spec):
    return format(x, spec).replace("-", "−")


def fmt_p(p):
    return "<0.001" if p < 0.001 else f"{p:.3f}" if p < 0.1 else f"{p:.2f}"


def table3():
    """Family share (%) of binned abundance per site type, GTDB R232 families, all 310 MAGs."""
    ph = json.load(open(S / "data" / "Saltern_phylogeny_r232.json"))
    rel = pd.read_csv(S / "data" / "Cliff_310MAG_relabund.tsv", sep="\t", index_col=0)
    mags = [m for m in rel.index if m + ".contigs__" in ph]
    assert len(mags) == 310
    types = {"Unrestored (R1, R2)": ("R1_", "R2_"), "Restored (SF2)": ("SF2",), "Reference (R2A)": ("R2A",)}
    fam = [ph[m + ".contigs__"]["Family"] or "no family assigned" for m in mags]
    cols, unb = {}, {}
    for t, pre in types.items():
        cs = [c for c in rel.columns if c.startswith(pre)]
        v = rel.loc[mags, cs].mean(axis=1)
        cols[t] = (v / v.sum()).groupby(fam).sum() * 100
        unb[t] = rel.loc["unbinned_reads", cs].mean() * 100
    tab = pd.DataFrame(cols).sort_values("Unrestored (R1, R2)", ascending=False)
    top, other = tab.head(10), tab.iloc[10:].sum()
    arch = {ph[m + ".contigs__"]["Family"] for m in mags if ph[m + ".contigs__"]["Kingdom"] == "Archaea"}
    rest = list(tab.index[10:])
    n_named, n_unassigned = sum(f != "no family assigned" for f in rest), fam.count("no family assigned")
    label = f"Other families ({n_named})" + (f" and MAGs without a family ({n_unassigned})" if "no family assigned" in rest else "")
    rows = [[f, *top.loc[f].tolist()] for f in top.index] + [[label, *other.tolist()]]
    return {"rows": rows, "unbinned": unb, "archaeal_in_top": [f for f in top.index if f in arch]}


def main(merged):
    doc = Document(merged)
    P = doc.paragraphs
    body_pPr, body_rpr = P[30]._p.find(qn("w:pPr")), first_rpr(P[30]._p)
    h3_pPr, h3_rpr = P[29]._p.find(qn("w:pPr")), first_rpr(P[29]._p)
    cap_pPr, cap_rpr = P[134]._p.find(qn("w:pPr")), first_rpr(P[134]._p)
    tcap_pPr, tcap_rpr = P[113]._p.find(qn("w:pPr")), first_rpr(P[113]._p)
    scap_pPr, scap_rpr = P[175]._p.find(qn("w:pPr")), first_rpr(P[175]._p)
    img_pPr = P[133]._p.find(qn("w:pPr"))
    t3 = doc.tables[2]._tbl
    tblPr_src = t3.find(qn("w:tblPr"))
    tcMar_src = next(t3.iter(qn("w:tcMar")), None)
    data = {"table3": table3()}
    corr = pd.read_csv(R / "guild_methane_correlations_r232.csv")
    percore = pd.read_csv(R / "methanogen_archaea_by_core_r232.csv")
    g = pd.read_csv(R / "guilds_r232.csv")
    guild = {c: sorted(g.loc[g[c], "iterativeID_r232"], key=lambda s: (s.rsplit(".", 1)[0], int(s.rsplit(".", 1)[1])))
             for c in ("mepn_consumer", "mepn_2020", "betaine_reducer", "betaine_2020", "methanogen", "tma_producer")}
    pending_comments = []   # (paragraph object, run element, text) -- added at the end

    # Keep handles to the paragraphs we edit (indices shift once we insert)
    p = {i: P[i] for i in (18, 30, 35, 47, 52, 53, 57, 69, 76, 77, 78, 83, 86, 95, 99, 102, 113, 129, 132, 133, 134, 148,
                           149, 151, 163, 175)}
    first_par = next(x for x in P if x.text.strip())

    # ---- Methods ---------------------------------------------------------------------------------------------------
    h, _ = new_par_after(p[30]._p, h3_pPr, "*MAG quality, taxonomy, and functional guilds*", h3_rpr)
    new_par_after(h, body_pPr,
        "The 310 dereplicated MAGs from Zhou et al. (2021) were re-assessed for this study. Completeness and "
        "contamination were estimated with CheckM v1.0.18 (67) (lineage workflow, reduced reference tree) in KBase (51), "
        "and the MAGs were classified with GTDB-Tk v2.8.0 (68) against GTDB release R232 (Table S5). Each MAG was "
        "labeled with its deepest named GTDB rank down to genus and numbered within that taxon (e.g., *Roseovarius*.3). "
        "Functional guilds were assigned from each MAG’s RAST annotation (53): methylphosphonate (MPn) consumers carry "
        "*phnJ*, the catalytic subunit of the C–P lyase (EC 4.7.1.1); betaine reducers carry *grdH*, the betaine-specific "
        "subunit of betaine reductase; and methanogens carry *mcrA*. MAGs with carnitine monooxygenase (*cntA*) or with "
        "genes of the choline TMA-lyase gene cluster were treated as candidate TMA producers; *cutC* itself was not "
        "annotated in any MAG.", body_rpr, lead_tab=True)
    tracked_replace(p[35]._p, "dereplicated medium and high quality metagenome-assembled genomes (MAGs)",
                    "310 dereplicated metagenome-assembled genomes (MAGs; 275 meet the MIMAG medium-quality criteria, Table S5)")
    replace_paragraph(p[52]._p,
        "To test which organisms track the observed methane emissions, the relative abundance of each MAG, and the summed "
        "abundance of each guild, were correlated with methane flux (Spearman’s ρ) across the 12 cores, with the D1 and D2 "
        "abundances merged by weighting with their DNA concentrations, and separately across the six unrestored cores. "
        "Correlations were computed on two bases: relative abundance versus areal CH~4~ flux, and both quantities divided "
        "by soil DNA concentration (the normalization of the model inputs). Because methane differs strongly among site types, each guild’s "
        "correlation was compared with those of 4,000 random sets of the same number of MAGs, and the fraction of random "
        "sets with an equal or higher ρ is reported. P values for the 310 individual MAGs were corrected with the "
        "Benjamini–Hochberg false discovery rate (FDR).", lead_tab=True)
    new_par_after(p[53]._p, body_pPr,
        "Co-occurrence among taxa was analyzed on centered log-ratio (CLR) transformed abundances of the 12 cores: only "
        "MAGs present (relative abundance > 10^-5^) in at least six cores were retained and then summed per taxon, and "
        "pairs were tested only if both taxa were present together in at least six cores. Site "
        "means were subtracted so that the within-site residuals reflect coupling rather than a shared response to site "
        "conditions. Residual correlations (Pearson) were tested against all 1,296 within-site permutations with "
        "Benjamini–Hochberg FDR correction, a sparse network of direct (conditional) associations was estimated from the "
        "residuals with the graphical lasso (scikit-learn GraphicalLassoCV), and modules were identified by Louvain "
        "modularity on the positive partial-correlation edges (|partial r| > 0.05).", body_rpr, lead_tab=True)

    # ---- Results -------------------------------------------------------------------------------------------------------
    h, _ = new_par_after(p[57]._p, h3_pPr, "*MAG quality, taxonomy, and functional guilds*", h3_rpr)
    new_par_after(h, body_pPr,
        "All 310 dereplicated MAGs were at least 50% complete, and the re-run of CheckM reproduced the original "
        "completeness and contamination values exactly. Of the 310 MAGs, 275 met the MIMAG medium-quality criteria, while "
        "35 exceeded 10% contamination, including two of the *Roseovarius* MAGs modeled as MPn consumers (*Roseovarius*.11 "
        "and .12, 20–21% contamination; Table S5). GTDB R232 classified 280 MAGs as Bacteria and 30 as Archaea, and "
        "resolved 290 MAGs to genus and 176 to species, compared with 201 and 19 in the original classification. Most "
        "other differences were renamed groups (e.g., Proteobacteria are now Pseudomonadota); the genera central to the "
        "MPn hypothesis (*Roseovarius*, *Albimonas*, and *Marivita*) were unchanged, whereas one of the two *Halomonas* "
        "MAGs is now placed in *Vreelandella*. Twenty-two MAGs carried *phnJ* (Table 4): 11 of the 12 *Roseovarius* MAGs, "
        "two each of *Marivita* and *Desulfotignum*, single MAGs of five other bacterial genera, and two haloarchaea "
        "(*Natrinema* and *Haloquadratum*). Nine MAGs carried *grdH* (three *Halanaerobium*, two BM714, two M55B127, "
        "DATFJG01, and UBA8670), and a single MAG carried *mcrA*: *Methanolobus* sp020055655 (Methanosarcinaceae; 98% "
        "complete). These gene-based guilds are subsets of the curated sets used for metabolic modeling: five of the 27 "
        "MAGs previously treated as MPn consumers carry only accessory *phn* genes, and six of the 15 modeled betaine "
        "reducers lack *grdH*.", body_rpr, lead_tab=True)
    # P69 (Figure S8 citation and R² values), P74-P77 (betaine models): see the edits JSON
    hg, _ = new_par_after(p[77]._p, h3_pPr, "*Guild abundances and methane flux*", h3_rpr)
    c = corr.set_index(["guild", "basis", "cores"])
    v = lambda gname, basis, cores, k: c.loc[(gname, basis, cores), k]
    replace_paragraph(p[78]._p,
        f"Across the 12 cores, the summed abundance of MPn consumers was positively correlated with methane flux "
        f"(Spearman ρ = {v('mepn_consumer', 'raw', '12 cores', 'spearman_rho'):.2f}, p = "
        f"{fmt_p(v('mepn_consumer', 'raw', '12 cores', 'p'))}; Figure 9, Table 5), as were betaine reducers "
        f"(ρ = {v('betaine_reducer', 'raw', '12 cores', 'spearman_rho'):.2f}, p = {fmt_p(v('betaine_reducer', 'raw', '12 cores', 'p'))}), "
        f"the single methanogen MAG (ρ = {v('methanogen', 'raw', '12 cores', 'spearman_rho'):.2f}, p = "
        f"{fmt_p(v('methanogen', 'raw', '12 cores', 'p'))}), and the archaea as a whole (ρ = "
        f"{v('archaea', 'raw', '12 cores', 'spearman_rho'):.2f}, p = {fmt_p(v('archaea', 'raw', '12 cores', 'p'))}). These "
        "correlations largely reflect the separation of the high-methane unrestored salterns from the other sites: 254 of "
        "the 310 MAGs were individually correlated with methane (FDR q < 0.05), and on this areal basis no guild correlated "
        "more strongly than random sets of the same number of MAGs (24–46% of random sets did as well). When abundances and "
        f"fluxes were both normalized by soil DNA (the normalization of the model inputs), MPn consumers stood out (ρ = "
        f"{v('mepn_consumer', 'DNA-normalised', '12 cores', 'spearman_rho'):.2f}; higher than 99% of random 22-MAG sets), "
        f"whereas betaine reducers did not (ρ = {v('betaine_reducer', 'DNA-normalised', '12 cores', 'spearman_rho'):.2f}; "
        f"{v('betaine_reducer', 'DNA-normalised', '12 cores', 'null_frac_ge_rho'):.0%} of random sets did as well; Figure S8). "
        "Within the six unrestored cores, no guild was significantly correlated with methane on either basis (MPn "
        f"consumers: ρ = {v('mepn_consumer', 'raw', '6 unrestored', 'spearman_rho'):.2f} areal and "
        f"{v('mepn_consumer', 'DNA-normalised', '6 unrestored', 'spearman_rho'):.2f} DNA-normalized, p ≥ 0.07).", lead_tab=True)
    pc = percore.set_index("core")
    unr = pc[pc.restoration == "unrestored"].methanogen_rel_abund * 100
    p_meth, w_meth = new_par_after(p[78]._p, body_pPr,
        "The only methanogen MAG, *Methanolobus* sp020055655, was detected (relative abundance > 10^-5^) only in the six "
        f"unrestored cores, where it made up {unr.min():.2f}–{unr.max():.2f}% of mapped reads (Table S6). Its correlation "
        "with methane across all 12 cores therefore reflects presence versus absence; within the unrestored cores its "
        f"abundance did not track methane (ρ = {num(v('methanogen', 'raw', '6 unrestored', 'spearman_rho'), '.2f')}, p = "
        f"{fmt_p(v('methanogen', 'raw', '6 unrestored', 'p'))}, n = 6).", body_rpr, lead_tab=True)
    hc, _ = new_par_after(p_meth, h3_pPr, "*Co-occurrence*", h3_rpr)
    new_par_after(hc, body_pPr,
        "With the R232 labels, the 30 archaeal MAGs entered the co-occurrence analysis (169 taxa, 13 of them archaeal). No "
        "within-site correlation passed FDR correction, and at most two pairs reached the 95th percentile of the "
        "permutation maximum |r| (0.993), so the network edges are conditional associations rather than individually significant pairs. "
        "The graphical-lasso network (82 taxa, 186 edges) resolved six modules (Q = 0.747; Figure 8), and all nine archaeal "
        "taxa in the network formed a single haloarchaeal module. Graphical-lasso fits to 20 within-site permutations of "
        "the residuals gave networks of similar size (167–534 nonzero edges, against 326 before the |partial r| > 0.05 "
        "threshold) and modularity (Q = 0.60–0.80); what exceeded the permutations was that edges joined related taxa "
        "(27% of edges joined taxa of the same phylum, against at most 15%; 18 of these 50 edges joined haloarchaea). *Roseovarius*, *Methanolobus*, *Halanaerobium*, and "
        "*Halomonas*/*Vreelandella* had no direct edges, but this mainly reflects their low residual variance, because the "
        "lasso penalty acts on unstandardized residuals; with standardized residuals each had 3–11 weak edges. Relaxing "
        "the prevalence filter gave no edge to any guild taxon that lacked one; the new edges mostly joined rarely "
        "detected taxa.", body_rpr, lead_tab=True)

    # ---- Discussion and Conclusions ----------------------------------------------------------------------------------------
    w = tracked_insert_before(p[83]._p, "Concentrations of formate",
        "In the dereplicated MAG set analyzed here, the only *mcrA*-containing MAG (classified by GTDB R232 as "
        "*Methanolobus* sp020055655, Methanosarcinaceae) was restricted to the unrestored salterns, but its abundance did "
        "not track methane within them (ρ = −0.54, p = 0.27, n = 6). ")
    pending_comments.append((p[83], w[0],
        "Naming: is this MAG the Ca. Methanosalis sp. SBSPR1A genome of ref. 58? GTDB R232 places it in Methanolobus "
        "(sp020055655); please reconcile the name used in the text."))
    tracked_insert_before(p[86]._p, "Here we propose that it also likely occurs",
        "However, the two *Nitrosopumilus* MAGs recovered here occurred almost only in the reference wetland (about 0.1% "
        "of reads each) and not in the unrestored salterns, and no *mpnS* was annotated by RAST in any MAG, so the source "
        "of methylphosphonate in the unrestored salterns remains unidentified. ")
    tracked_replace(p[95]._p,
        "MPn-consuming organisms were much more associated with methane emissions than betaine-consuming organisms (R2 = 0.798).",
        "MPn-consuming organisms were more specifically associated with methane emissions than betaine-consuming organisms "
        "only when abundances and emissions were normalized by soil DNA; on an areal basis, both guilds, the single "
        "methanogen MAG, and most of the community followed the same unrestored-versus-restored gradient (Figure 9, Table 5).")
    tracked_replace(p[99]._p,
        "Secondly, the observed methane emissions were most correlated with the abundances of methylphosphonate-consuming "
        "microbes, more so than with the abundances of betaine-consuming microbes, which was our second hypothesized "
        "methanogenic substrate.",
        "Secondly, MPn consumers (MAGs carrying the C–P lyase gene *phnJ*) were the only guild whose abundance tracked "
        "methane emissions more closely than random sets of MAGs, although only when abundances and emissions were "
        "normalized by soil DNA; on an areal basis, betaine reducers, the single methanogen MAG, and most of the community "
        "followed the same site gradient.")
    p[102]._p.append(ins(runs_from(
        "MAG quality assessments, GTDB R232 classifications, guild assignments, and the correlation and co-occurrence "
        "analyses are in the MicrobiomeNotebooks repository (https://github.com/ModelSEED/MicrobiomeNotebooks; folders "
        "salterns/MAG_QC and salterns/correlational/r232).", first_rpr(p[102]._p))))
    # P47 (Halomonas isolates, GTDB-Tk v1.7.0 -> v2.8.0 / R232): see the edits JSON
    pending_comments.append((p[18], p[18]._p.findall(qn("w:r"))[0],
        "The modeling clause is now qualified (rank order, not magnitude). If the abstract keeps a correlation claim, also "
        "note that guild abundance tracks methane more specifically for MPn consumers only on the DNA-normalized basis, "
        "and that the single methanogen MAG does not track methane within the unrestored salterns (Table 5)."))

    # ---- Table 3 (replaced), Tables 4 and 5 (new) -------------------------------------------------------------------------
    unb = data["table3"]["unbinned"] if data else None
    w = replace_paragraph(p[113]._p,
        "Table 3. Percent of binned (MAG) relative abundance contributed by each GTDB R232 family in the unrestored (R1, R2), "
        "restored (SF2), and reference (R2A) wetlands, for the ten families most abundant in the unrestored salterns "
        "(† archaeal family). Percentages are of the 310 MAGs, including the 30 archaeal MAGs, and exclude unbinned reads, "
        f"which made up {unb['Unrestored (R1, R2)']:.0f}%, {unb['Restored (SF2)']:.0f}%, and {unb['Reference (R2A)']:.0f}% of "
        "reads in each site type, respectively. Because the table is sorted by the unrestored salterns, most of the restored "
        "and reference abundance falls in “Other families”.")
    pending_comments.append((p[113], w[0],
        "The previous Table 3 came from correlational/abundances.ipynb, which had three faults (fixed there on 2026-10-21): "
        "(1) its selection 'R1' in col or 'R2' in col also took the six R2A cores, so the Unrestored column pooled R1, R2 "
        "and R2A; (2) the 30 archaea had no 2020 taxonomy entry and were dropped with the unbinned reads, so the caption's "
        "68% is unbinned plus archaeal reads over those 18 cores (the data give 52%: 47% R1, 58% R2); (3) the notebook's "
        "column order was unrestored, reference, restored, which was copied under the headers Unrestored, Restored, "
        "Reference, swapping the last two. Recomputed here with GTDB R232 families and all 310 MAGs; the corrected "
        "2020-family values are in 'Table 3 - corrected values (2020 GTDB families).csv'."))
    delete_table(t3)
    arch = set(data["table3"]["archaeal_in_top"])
    rows3 = [[f + ("†" if f in arch else ""), *[f"{x:.1f}" for x in vals]] for f, *vals in data["table3"]["rows"]]
    new_t3 = build_table(tblPr_src, tcMar_src, ["GTDB R232 family", "Unrestored (R1, R2)", "Restored (SF2)", "Reference (R2A)"],
                         rows3, [2880, 1890, 1860, 1890], body_rpr, size=22)
    sep, _ = new_par_after(t3, body_pPr, "", body_rpr)       # keeps Word from merging the deleted and new tables
    sep.addnext(new_t3)
    cap4, _ = new_par_after(new_t3, tcap_pPr,
        "Table 4. Functional guilds assigned from the RAST gene annotations of the 310 MAGs (GTDB R232 labels), compared "
        "with the curated sets used for metabolic modeling. Members are listed by iterativeID (taxon.number; Table S5).",
        tcap_rpr)
    rows4 = [
        ["MPn consumers", "*phnJ* (C–P lyase, EC 4.7.1.1)", str(len(guild["mepn_consumer"])), compress(guild["mepn_consumer"]),
         "5 with only accessory *phn* genes: CAZRZN01.1, DAINOF01.1, *Thiolinea*.1 (simulated); *Gaopeijia*.1 (did not "
         "grow); *Halorubrum*.3 (not simulated)"],
        ["Betaine reducers", "*grdH* (betaine reductase, betaine-specific subunit)", str(len(guild["betaine_reducer"])),
         compress(guild["betaine_reducer"]),
         "6 without *grdH* (all simulated in the original set): " + ", ".join(it_id(i) for i in sorted(set(guild["betaine_2020"]) - set(guild["betaine_reducer"])))],
        ["Methanogens", "*mcrA*", str(len(guild["methanogen"])), compress(guild["methanogen"]) + " (sp020055655)", "—"],
        ["TMA-route candidates", "*cntA*; genes of the choline TMA-lyase cluster (no *cutC* itself)", str(len(guild["tma_producer"])),
         "J153.1 (*cntA*); CAZFJG01.2, T1Sed10-7.1, UBA11579.1 (cluster genes only)", "Not modeled; weak evidence"],
    ]
    t4 = build_table(tblPr_src, tcMar_src, ["Guild", "Marker gene", "n", "Members", "In the curated modeling set but lacking the marker"],
                     rows4, [1350, 1650, 760, 2700, 2060], body_rpr, size=20,
                     align=["left", "left", "center", "left", "left"])
    cap4.addnext(t4)
    cap5, _ = new_par_after(t4, tcap_pPr,
        "Table 5. Spearman correlations (p value) between guild abundance and methane flux across the 12 cores and across "
        "the six unrestored cores, on an areal basis and with both quantities normalized by soil DNA. ‘Random ≥ ρ’ is the "
        "fraction of 4,000 random sets of the same number of MAGs with an equal or higher 12-core ρ.", tcap_rpr)
    labels = [("mepn_consumer", "MPn consumers (*phnJ*), 22"), ("mepn_2020", "MPn modeling set, 27"),
              ("betaine_reducer", "Betaine reducers (*grdH*), 9"), ("betaine_2020", "Betaine modeling set, 15"),
              ("methanogen", "Methanogen (*Methanolobus*.1)"), ("archaea", "All archaea, 30")]
    cell = lambda gname, basis, cores: f"{num(v(gname, basis, cores, 'spearman_rho'), '+.2f')} ({fmt_p(v(gname, basis, cores, 'p'))})"
    rows5 = [[lab, cell(gn, "raw", "12 cores"), cell(gn, "raw", "6 unrestored"), f"{v(gn, 'raw', '12 cores', 'null_frac_ge_rho'):.0%}",
              cell(gn, "DNA-normalised", "12 cores"), cell(gn, "DNA-normalised", "6 unrestored"),
              f"{v(gn, 'DNA-normalised', '12 cores', 'null_frac_ge_rho'):.1%}"] for gn, lab in labels]
    t5 = build_table(tblPr_src, tcMar_src, ["Guild, MAGs", "Areal: 12 cores", "Areal: 6 unrestored", "Areal: random ≥ ρ",
                                            "DNA-norm.: 12 cores", "DNA-norm.: 6 unrestored", "DNA-norm.: random ≥ ρ"],
                     rows5, [2040, 1080, 1080, 1080, 1080, 1080, 1080], body_rpr, size=18)
    cap5.addnext(t5)

    # ---- Figure 8 (replaced), Figure 9 (new), Figure S8 (replaced) ----------------------------------------------------------
    # Figure 7: single panel on the areal basis, betaine = nine grdH carriers (caption P129 in the edits JSON)
    delete_paragraph(P[128]._p)
    picture_after(doc, P[128]._p, HERE / "figure7_r232.png", 5.6, img_pPr)
    for i in (132, 133):
        delete_paragraph(p[i]._p)
    picture_after(doc, p[133]._p, R / "within_condition_louvain_network_r232.png", 6.0, img_pPr)
    w = replace_paragraph(p[134]._p,
        "Figure 8. Co-occurrence network of the community members (GTDB R232 labels; archaea included). Edges are partial "
        "correlations from a graphical lasso (L1-regularized precision matrix, regularization chosen by cross-validation) "
        "fit on within-site residuals of the CLR-transformed abundances of the 12 cores, for the 169 taxa formed from "
        "MAGs present (relative abundance > 10^-5^) in at least six cores; edges with |partial r| > 0.05 are shown (82 taxa, 186 edges). "
        "Node size encodes mean relative abundance, node color the phylum, and edge color the partial correlation (red = "
        "negative, blue = positive). Louvain community detection on the positive edges resolves six modules (Q = 0.747), "
        "shown as shaded hulls; all nine archaeal taxa in the network fall in one haloarchaeal module. No individual "
        "within-site correlation passed FDR correction, so edges indicate conditional associations rather than demonstrated "
        "interactions; graphical-lasso fits to permuted residuals give networks of similar size and modularity, and only the "
        "excess of edges between related taxa exceeds them. Because the lasso penalty acts on unstandardized residuals, taxa "
        "with low residual variance, including *Roseovarius*, *Methanolobus*, *Halanaerobium*, and *Halomonas*/*Vreelandella*, "
        "have no edges above the threshold; the MPn consumers retained in the network are *Haloquadratum* and "
        "*Phaeodactylibacter*, and the betaine reducers are DATFJG01 and UBA8670 (in the same module) and M55B127.")
    pending_comments.append((p[134], w[0],
        "Figure 8 is now a single panel: the previous panel a (plain graphical-lasso layout) was not re-rendered; this "
        "Louvain view shows the same edges. Code: correlational/r232/rerun_correlational_r232.py."))
    p9 = picture_after(doc, p[134]._p, R / "guild_methane_raw_r232.png", 6.5, img_pPr)
    new_par_after(p9, cap_pPr,
        "Figure 9. Guild abundance versus methane flux across the 12 cores, with guilds assigned from RAST gene annotations "
        "of the GTDB R232-labeled MAGs: MPn consumers (*phnJ*), betaine reducers (*grdH*), the curated sets used for "
        "metabolic modeling, the single methanogen MAG (*mcrA*, *Methanolobus* sp020055655), and all archaea. Abundance is "
        "the summed relative abundance (D1 and D2 merged with DNA weighting); methane is the areal CH~4~ flux on a "
        "symmetric-log axis because the reference wetland (R2A) showed net uptake. Each panel gives Spearman’s ρ for all "
        "12 cores and for the six unrestored cores, and the fraction of 4,000 random sets of the same number of MAGs with "
        "an equal or higher 12-core ρ. Hollow markers are below detection (relative abundance ≤ 10^-5^).", cap_rpr)
    delete_paragraph(p[148]._p)
    picture_after(doc, p[148]._p, R / "guild_methane_dnanorm_r232.png", 6.5, img_pPr)
    replace_paragraph(p[149]._p,
        "Figure S8. As Figure 9, with guild abundance and methane flux both normalized by soil DNA concentration (the "
        "normalization of the model inputs). On this basis the MPn consumers (*phnJ*) correlate with methane more "
        "strongly than 99% of random 22-MAG sets (ρ = 0.92), whereas the other guilds do not stand out (Table 5).")
    delete_paragraph(p[151]._p)
    pending_comments.append((p[163], p[163]._p.findall(qn("w:r"))[0],
        "Not regenerated with R232 labels. Regeneration needs no new statistics: the figure is drawn from the within-site "
        "correlation matrix and the graphical-lasso edges, whose R232 versions exist (correlational/r232/C_within_r232.csv, "
        "169 taxa including 13 archaeal; glasso_direct_edges_r232.csv). The 2020 plotting code, copied out of "
        "analysis_2_marimo.py, reproduces this image exactly and renders an R232 version. Still needed: a choice of guild "
        "sets to highlight (gene-based or curated), new colors for the R232 taxa, a larger canvas for 169 rows, porting the "
        "renderer into rerun_correlational_r232.py, and the caption fixes noted in the comment below. Figure S11 is also "
        "not cited anywhere in the text: cite it in the Co-occurrence results or drop it. Details: "
        "correlational/r232/FigureS11_regeneration.md."))

    # ---- Tables S5 (external) and S6 -----------------------------------------------------------------------------------------
    s5, _ = new_par_after(p[175]._p, scap_pPr,
        "Table S5. GTDB R232 classification, CheckM completeness and contamination, functional guild assignment with the "
        "supporting RAST annotation, mean relative abundance by site type, and correlation with methane flux for each of "
        "the 310 dereplicated MAGs (separate file: Table S5 - MAG taxonomy, quality and guilds (GTDB R232).xlsx).", scap_rpr)
    s6, _ = new_par_after(s5, scap_pPr,
        "Table S6. Methane flux and the relative abundance (% of mapped reads; D1 and D2 merged with DNA weighting) of the "
        "methanogen MAG (*Methanolobus* sp020055655) and of all 30 archaeal MAGs in each core.", scap_rpr)
    rows6 = [[r.core, r.restoration, num(r.CH4_umol_m2_d, ".1f"), f"{r.methanogen_rel_abund * 100:.4f}", f"{r.archaea_rel_abund * 100:.2f}"]
             for r in percore.sort_values(["restoration", "core"], key=lambda s: s.map({"unrestored": 0, "restored": 1, "reference": 2}) if s.name == "restoration" else s).itertuples()]
    t6 = build_table(tblPr_src, tcMar_src, ["Core", "Site type", "CH~4~ flux (µmol m^-2^ d^-1^)", "*Methanolobus*.1 (%)", "All archaea (%)"],
                     rows6, [1100, 1400, 2500, 1820, 1700], body_rpr, size=20)
    s6.addnext(t6)

    # ---- 2026-10-26 edits: modelling text and numbers, betaine set, Halomonas isolates, smaller corrections ------------------
    apply_edits(P, pending_comments, EDITS)
    # fitted-line equations (Office Math): consistent areal-basis fits (modelling_consistent_basis_r232.json), units without 1/DNA
    for i, slope, icpt in ((70, ("=0.51*methan", "=−0.20*methan"), ("+138 [", "+907 [")),
                           (72, ("=0.95*methan", "=−0.32*methan"), ("+224 [", "+1563 [")),
                           (75, ("=0.24*methan", "=−0.71*methan"), ("+4632 [", "+3192 ["))):
        tracked_math_swap(P[i]._p, dict([slope, icpt, ("*day*", "*day")]), drop_fraction=("ug DNA", "g soil"))

    # ---- comments, then save --------------------------------------------------------------------------------------------------
    pending_comments.insert(0, (first_par, first_par._p.findall(qn("w:r"))[0],
        "R232 review (2026-10-21, revised 2026-10-26): tracked changes add MAG quality, taxonomy (GTDB R232), and gene-based "
        "guild methods and results; replace the abundance–methane results, Table 3, Figure 8, and Figure S8; and add "
        "Figure 9, Tables 4–5, S5 (Excel), and S6. The 2026-10-26 revision corrects the modeling text and numbers: on a "
        "dimensionally consistent areal basis the simulated methane is of the measured magnitude but does not track the "
        "differences among the six unrestored cores (the earlier high fits came from a 1/DNA term in the model input); it also restricts the betaine models to the nine grdH carriers, replaces "
        "Figure 7 with a single panel, and reclassifies the two Halomonas isolates with GTDB R232. Items flagged but not "
        "changed are in these comments; see the accompanying review notes."))
    for para, run_el, text in pending_comments:
        comment(doc, para, run_el, text)
    doc.save(OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main(sys.argv[1])
