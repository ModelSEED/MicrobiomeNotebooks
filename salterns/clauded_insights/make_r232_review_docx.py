#!/usr/bin/env python
"""
R232 review of the saltern manuscript: a copy of the Google Doc as of 2026-11-03 (the local mirror
../Metabolic Manuscript - Google Doc mirror (2026-11-03).docx, from make_gdoc_mirror_docx.py) with real Word tracked
changes (w:ins / w:del) and review comments that bring in the 2026 MAG re-analysis:
CheckM re-run, GTDB-Tk 2.8.0 / R232 taxonomy, gene-based functional guilds, guild-methane correlations
with a random-set comparison, the methanogen-methane correlation, and the R232 co-occurrence network
(../correlational/r232/). Rejecting the changes by "Claude (R232 review)" restores the Google Doc text, tables and
figures; Susannah Tringe's open suggestions are her own tracked changes in the mirror.

The mirror is built from the 2 July export after the docx skill's merge_runs.py (adjacent identically formatted runs
merged, text unchanged), so that sentences can be located inside single runs:
    python <docx-skill>/scripts/merge_runs.py "../Metabolic Manuscript.docx" -o merged.docx
    ~/Documents/py_venv/bin/python make_gdoc_mirror_docx.py merged.docx
    ~/Documents/py_venv/bin/python make_r232_review_docx.py "../Metabolic Manuscript - Google Doc mirror (2026-11-03).docx"
Output: ../Metabolic Manuscript - R232 review (tracked changes).docx

The 2026-10-26 revision adds the tracked edits listed in r232_review_edits_2026-10-26.json:
- the modelling text and numbers on the areal basis (modelling_numbers_r232.py);
- betaine models restricted to the nine grdH carriers (../correlational/betaine_grdH_r232/);
- the R232 reclassification of the two Halomonas isolates (../MAG_QC/gtdbtk_isolates_r232/);
- smaller corrections.

It also replaces Figure 7 with figure7_r232.png (figure7_r232.py); since 2026-11-24 with figure7_pergram.png
(reassembly_2026/modelling_2026/simulations/figure7_pergram.py), the per-gram framing of the 2026 guild models, and the
modelling text comes from r232_review_edits_2026-11-24.json (make_modelling_edits_2026-11-24.py, layered on the frozen
r232_review_edits_2026-11-23.json).
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
from docx.text.paragraph import Paragraph
from docx.text.run import Run

HERE = Path(__file__).resolve().parent
S = HERE.parent
R = S / "correlational" / "r232"            # the 310 original MAGs, GTDB R232 labels
RF = S / "correlational" / "final585"       # the re-assembled MAG set (585 species): abundance, guilds, co-occurrence
DATE_FINAL585 = "2026-11-20T00:00:00Z"      # correlational analyses moved to the re-assembled MAG set
DATE_MODELS = "2026-11-24T00:00:00Z"        # per-gram modelling text, equation deletions and Figure 7 (2026 guild models)
SIMF = S / "reassembly_2026" / "modelling_2026" / "simulations"
OUT = S / "Metabolic Manuscript - R232 review (tracked changes).docx"
AUTHOR, INITIALS, DATE = "Claude (R232 review)", "CR", "2026-10-26T00:00:00Z"
DATE_MEPN_REPAIR = "2026-11-08T00:00:00Z"
DATE_REASSEMBLY = "2026-11-19T00:00:00Z"   # re-assembly / re-binning / annotation additions (reassembly_2026/)   # edits revised for the MePn-obligate model repair (mepn_obligacy_check/)
EDITS = HERE / "r232_review_edits_2026-10-26.json"
EDITS_1117 = HERE / "r232_review_edits_2026-11-24.json"   # for the 2026-11-17 mirror (Abstract split; see that file's note)
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
    """Track-delete everything visible in p: plain runs, runs in hyperlinks, and the runs of other authors' open
    insertions (deleting an insertion nests w:del inside that w:ins, as Word does)."""
    for child in list(p):
        if is_content(child):
            wrap_del(child)
        elif child.tag in (qn("w:hyperlink"), qn("w:ins")):
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
    global DATE
    default = DATE
    for e in sorted(json.load(open(path))["edits"], key=lambda e: e["kind"] == "comment"):
        DATE = e.get("date", default)                    # tracked changes carry the edit's own date when it has one
        p = P.raw(e["mirror_index"]) if "mirror_index" in e and P.shift else P[e["paragraph"]]
        if e["kind"] == "replace":
            tracked_replace(p._p, e["old"], e["new"])
        elif e["kind"] == "insert_before":
            tracked_insert_before(p._p, e["old"], e["new"])
        elif e["kind"] == "replace_paragraph":
            replace_paragraph(p._p, e["new"], lead_tab=p.text.startswith("\t"))
        else:
            _, start = _locate(p._p, e["old"])
            pending_comments.append((p, _cut_at(p._p, start), e["new"]))
    DATE = default


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
    """Family share (%) of binned abundance per site type, GTDB R232 families of the 585 species representatives of the
    re-assembled MAG set (CoverM read counts, reassembly_2026/abundance_final/; labels correlational/final585/)."""
    ph = json.load(open(RF / "phylogeny_final585.json"))
    rel = pd.read_csv(S / "reassembly_2026" / "abundance_final" / "final585_relabund_24sediment.tsv", sep="\t", index_col=0)
    mags = [m for m in rel.index if m in ph]
    assert len(mags) == 585
    types = {"Unrestored (R1, R2)": ("R1_", "R2_"), "Restored (SF2)": ("SF2",), "Reference (R2A)": ("R2A",)}
    fam = [ph[m]["Family"] or "no family assigned" for m in mags]
    cols, unb = {}, {}
    for t, pre in types.items():
        cs = [c for c in rel.columns if c.startswith(pre)]
        v = rel.loc[mags, cs].mean(axis=1)
        cols[t] = (v / v.sum()).groupby(fam).sum() * 100
        unb[t] = rel.loc["unbinned_reads", cs].mean() * 100
    tab = pd.DataFrame(cols).sort_values("Unrestored (R1, R2)", ascending=False)
    top, other = tab.head(10), tab.iloc[10:].sum()
    arch = {ph[m]["Family"] for m in mags if ph[m]["Kingdom"] == "Archaea"}
    rest = list(tab.index[10:])
    n_named, n_unassigned = sum(f != "no family assigned" for f in rest), fam.count("no family assigned")
    label = f"Other families ({n_named})" + (f" and species without a family ({n_unassigned})" if "no family assigned" in rest else "")
    rows = [[f, *top.loc[f].tolist()] for f in top.index] + [[label, *other.tolist()]]
    return {"rows": rows, "unbinned": unb, "archaeal_in_top": [f for f in top.index if f in arch]}


REASM = S / "reassembly_2026"


def latest_numbers():
    fs = sorted(REASM.glob("manuscript_numbers_*.json"))
    assert fs, "run reassembly_2026/manuscript_numbers.py first"
    return json.load(open(fs[-1])), fs[-1].name


def kb(x):
    return "n/a" if x is None else f"{x / 1000:.1f}"


def add_reassembly_2026(mq_meth, res_anchor, t6, pending_comments, h3_pPr, h3_rpr, body_pPr, body_rpr, scap_pPr, scap_rpr,
                        tblPr_src, tcMar_src, wrap):
    """2026 re-assembly: methods (after the MAG-quality methods), results (after the MAG-quality results, ahead of the
    guild, abundance and co-occurrence results that use the re-assembled MAG set), Tables S7 and S8 (after S6).
    Static numbers come from reassembly_2026/{comparison,binning_experiments,annotation_runs}/ (reviewed); the second-pass
    rollout numbers come from the latest reassembly_2026/manuscript_numbers_<date>.json (partial while it runs)."""
    N, nfile = latest_numbers()
    T = N["totals"]
    sp_, fg_ = REASM / "comparison_rollout" / "summary.json", REASM / "comparison_rollout" / "final_guilds_summary.json"
    S2top = json.load(open(sp_)) if sp_.exists() else None
    FG = json.load(open(fg_)) if fg_.exists() else None
    FS = bool(S2top and FG and "final" in S2top)          # the final MAG set (second pass + kept co-assembly MAGs)
    W = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine"}
    def nw(n):
        return W.get(n, f"{n}")
    if FS:
        CO, F = S2top["coassembly"], S2top["final"]
        co_fp = FG["kofam_sources"].get("kofam (first pass)", 0)       # kept co-assembly MAGs the first pass had kept
        co_new = CO["kept"] - co_fp
    snap = N["snapshot_utc"][:10]
    # ---- Methods
    hm, _ = new_par_after(mq_meth, h3_pPr, "*Re-assembly, re-binning, and functional annotation of the metagenomes*", h3_rpr)
    m1, w_m1 = new_par_after(hm, body_pPr,
        "To test whether the original assemblies and bins limited these analyses, the 39 metagenomes were re-assembled and "
        "re-binned from the raw reads. In a first pass, run in KBase, reads were trimmed with Trimmomatic v0.39 "
        "(ILLUMINACLIP TruSeq3-PE-2:2:30:10, LEADING:3, TRAILING:3, SLIDINGWINDOW:4:15, MINLEN:50) and assembled per "
        "sample with metaSPAdes v3.15.3 (k = 33–127; contigs ≥ 1 kb), and the sediment cores of each site were "
        "co-assembled with MEGAHIT v1.2.9. Contigs were binned with MetaBAT2, MaxBin2 v2.2.4 and CONCOCT v1.1 using the "
        "coverage of every library from the same site (9–12 libraries per site); the three bin sets were refined with "
        "DAS Tool v1.1.7; bins were assessed with CheckM v1.1.2 (lineage workflow, reduced tree); and bins ≥ 50% complete "
        "and < 10% contaminated were classified with GTDB-Tk v2.8.0 against GTDB R232.", body_rpr, lead_tab=True)
    m2, _ = new_par_after(m1, body_pPr,
        "A second pass reproduced the read preparation of the original assemblies, which were built from merged, "
        "Tadpole-corrected reads (Zhou et al. 2022). Adapters, sequencing artifacts and PhiX were removed with BBDuk "
        "(BBTools v39.81); read pairs were merged with BBMerge with k-mer extension (rem, extend2 = 50, k = 62, ecct); "
        "merged and unmerged reads were error-corrected with Tadpole (k = 62) and quality-trimmed (trimq = 10, minimum "
        "length 45 nt); and each sample was assembled with SPAdes v3.15.3 in its standard (non-metagenomic) mode "
        "(--only-assembler, k = 33–127). The BBDuk-trimmed reads of every library from the same site were mapped to each "
        "assembly (contigs ≥ 1 kb) with Bowtie2 v2.5.5, and contigs were binned with MetaBAT2 v2.12.1, MaxBin2 v2.2.7 and "
        "CONCOCT v1.1.0 with the first-pass parameters, followed by DAS Tool v1.1.7. Bins were assessed with CheckM2 v1.1.0, "
        "CheckM v1.1.3 and GUNC v1.1.1 (proGenomes 2.1), and bins ≥ 50% complete and < 10% contaminated by CheckM2 that "
        "passed GUNC were kept. Before all samples were processed, the two passes were compared on two samples (R1_A_D1 "
        "and SF2_C_D1) against pre-specified criteria: the numbers of medium- and high-quality MAGs, N50, contamination, "
        "GUNC pass rate, and recovery of the original MAGs. The second pass met every criterion in R1_A_D1 but exceeded "
        "the contamination criterion in the more diverse SF2_C_D1, so it was applied to all 39 metagenomes with the "
        "CheckM2 and GUNC filter described above. Kept MAGs were classified with GTDB-Tk v2.8.0 against GTDB R232. For "
        "comparison, the first-pass bins and the original MAGs were "
        "re-scored with the same CheckM2 and GUNC settings, and genomes were matched across sets with skani (ANI ≥ 95% "
        "over ≥ 50% of the shorter genome). Species clusters were formed jointly for the three MAG sets by greedy "
        "clustering at the same thresholds, starting from the genome with the highest completeness − 5 × contamination."
        + (f" The {CO['bins']} bins of the four first-pass site co-assemblies were assessed and classified in the same "
           f"way, and the {CO['kept']} that met the same criteria were combined with the second-pass MAGs into the final "
           f"MAG set (Table S8)." if FS else ""),
        body_rpr, lead_tab=True)
    m3, _ = new_par_after(m2, body_pPr,
        (f"Genes of the 1,324 first-pass MAGs, the {T['new_kept']} second-pass MAGs and the {co_new} co-assembly MAGs "
         "kept only under the second-pass criteria" if FS else "Genes of the 1,324 first-pass MAGs") +
        " were predicted with Prodigal v2.6.3 (metagenomic mode) and assigned KEGG "
        "orthologs with hmmsearch against KOfam (release 2026-10-26), using each ortholog's score threshold and keeping "
        "the best-scoring passing ortholog per gene, as DRAM does. Because gene callers stop the selenoproteins of glycine "
        "and betaine reductase at their UGA (selenocysteine) codon, TGA-split gene pairs of GrdA, GrdB, GrdF and GrdH were "
        "joined, translated through the selenocysteine codon and re-scored. DRAM v1.5.0 (KOfam, Pfam, dbCAN and MEROPS; "
        "contigs ≥ 1 kb) was run on one MAG per " + ("first-pass " if FS else "") + "species cluster (558 MAGs) with a custom distillate of methane-related "
        "genes. Guilds were assigned from the KOfam calls: MPn consumers carry *phnJ*; betaine reducers carry *grdH* or "
        "*grdI*; methanogens carry *mcrA*; methylamine users carry *mttB*, *mtbB* or *mtmB* together with *mcrA*; and TMA "
        "producers carry *cutC* with *cutD*, or *cntA* with *cntB*.", body_rpr, lead_tab=True)
    if FS:
        new_par_after(m3, body_pPr,
            f"The relative abundance of each of the {F['species']} species of the final set in each metagenome was estimated "
            "by mapping the BBDuk-trimmed read pairs to the species representatives with CoverM v0.7.0 (genome mode; "
            "minimap2 v2.30; reads aligned with ≥ 95% identity over ≥ 75% of their length) and dividing each species’ read "
            "count by the total number of trimmed reads; reads not assigned to any representative are reported as "
            "unbinned. Each "
            "representative was labeled like the original MAGs (deepest named GTDB rank down to genus, numbered within that "
            "taxon), and a species that contains an original MAG keeps that MAG’s label. A species was assigned to a guild "
            "when any of its MAGs carried the marker gene.", body_rpr, lead_tab=True)
    pending_comments.append((wrap(m1), w_m1.find(qn("w:r")),
        "Added 2026-11-16 (re-assembly). Citations to add to the reference list: Trimmomatic (Bolger et al. 2014), SPAdes "
        "and metaSPAdes (Bankevich et al. 2012; Nurk et al. 2017), MEGAHIT (Li et al. 2015), MetaBAT2 (Kang et al. 2019), "
        "MaxBin2 (Wu et al. 2016), CONCOCT (Alneberg et al. 2014), DAS Tool (Sieber et al. 2018), BBTools/BBMerge "
        "(Bushnell et al. 2017), Bowtie2 (Langmead and Salzberg 2012), CheckM2 (Chklovski et al. 2023), GUNC (Orakov et al. "
        "2021), skani (Shaw and Yu 2023), KOfam/KofamKOALA (Aramaki et al. 2020), DRAM (Shaffer et al. 2020), Prodigal "
        "(Hyatt et al. 2010), CoverM (Aroney et al. 2025), minimap2 (Li 2018). Run records: reassembly_2026/README.md, "
        "WORKFLOW.md and annotation_runs/*/README.md."))
    # ---- Results
    hr, _ = new_par_after(res_anchor, h3_pPr, "*Re-assembly and re-binning of the saltern metagenomes*", h3_rpr)
    r1, _ = new_par_after(hr, body_pPr,
        "The first re-assembly pass yielded 1,324 MAGs ≥ 50% complete and < 10% contaminated (323 with ≥ 90% completeness "
        "and < 5% contamination), which formed 558 species-level clusters, against 296 clusters among the 310 original "
        "MAGs. The two sets shared 250 species; 308 species were found only in the re-assembly and 46 only in the original "
        "set, and archaeal species increased from 30 to 47. The re-assembled MAGs were, however, about twice as fragmented "
        "(median N50 of the best MAG per species 8.4 kb, against 19.1 kb for the original MAGs). The original MAGs were "
        "binned not from the public JGI assemblies but from per-sample re-assemblies, whose length distribution we "
        "reconstructed from the original contig names: summed over the 39 samples, the first-pass assemblies held 3,691 Mbp "
        "in contigs ≥ 2.5 kb and 815 Mbp in contigs ≥ 10 kb, against about 6,300 and 2,740 Mbp in the original assemblies "
        "(median per-sample ratios 0.57 and 0.27).", body_rpr, lead_tab=True)
    r2, _ = new_par_after(r1, body_pPr,
        "Reproducing the original read preparation restored the original contiguity. For R1_A_D1, merged reads assembled "
        "in standard mode gave 250.5 Mbp in contigs ≥ 2.5 kb and 133.8 Mbp in contigs ≥ 10 kb, 0.93 and 1.02 times the "
        "original assembly (first pass: 162.8 and 32.6 Mbp), whereas Tadpole correction without merging gave no "
        "improvement (164.1 and 32.5 Mbp). Binned in the same way, this assembly gave 62 MAGs passing CheckM2 (≥ 50% "
        "complete, < 10% contaminated) against 47 from the first pass, 27 high-quality MAGs against 11, and a median N50 "
        "of 15.4 kb against 4.7 kb, with similar contamination (3.3% against 3.5%), a higher GUNC pass rate (0.87 against "
        "0.79), and 10 of the 11 original dereplicated MAGs from that sample recovered (against 8). In the more diverse "
        "SF2_C_D1 sample, it gave 7 MAGs against 4 and 2 high-quality MAGs against none, but with higher median "
        "contamination (4.7% against 1.7%).", body_rpr, lead_tab=True)
    if T.get("orig2020_kept") is not None:
        orig_txt = (f"and {T['orig2020_kept']} ({T['orig2020_hq']} high-quality; median N50 "
                    f"{kb(N['pooled_orig2020_med_n50'])} kb) of the {T['orig2020_mags']} original MAGs from the same samples "
                    "passed the same criteria")
    else:
        orig_txt = f"while the same samples had yielded {T['orig2020_mags']} original MAGs"
    rec_txt = ""
    if T.get("orig2020_recovered") is not None:
        rec_txt = (f" Of the {T['orig2020_mags']} original MAGs, {T['orig2020_recovered']} were recovered by a second-pass "
                   f"MAG of the same sample (ANI ≥ 95%), including {T['orig2020_lost_recovered']} of the "
                   f"{T['orig2020_lost_in_first_pass']} medium-quality original MAGs whose species the first pass had not "
                   "recovered.")
    final = N["n_complete"] == N["n_samples"]
    lead = ("Applied to all 39 metagenomes (Table S7), the second pass kept" if final else
            f"The second pass is being applied to all 39 metagenomes; {N['n_complete']} samples had been completed when "
            f"this draft was prepared (Table S7). In these samples it kept")
    r3, w_r3 = new_par_after(r2, body_pPr,
        f"{lead} {T['new_kept']} MAGs ({T['new_hq']} high-quality; "
        f"median N50 {kb(N['pooled_new_med_n50'])} kb), against {T['kbase2026_kept']} ({T['kbase2026_hq']}; "
        f"{kb(N['pooled_kbase2026_med_n50'])} kb) from the first pass scored with the same CheckM2 and GUNC criteria, "
        f"{orig_txt}.{rec_txt}", body_rpr, lead_tab=True)
    S2 = json.load(open(REASM / "comparison_rollout" / "summary.json")) if final else None
    if final:
        pending_comments.append((wrap(r3), w_r3.find(qn("w:r")),
            f"Final second-pass numbers ({N['n_complete']} of {N['n_samples']} samples; run finished 2026-11-19; sources "
            f"reassembly_2026/{nfile}, comparison_rollout/summary.json and README.md; samples that fell back to "
            f"metaSPAdes: {', '.join(N['meta_fallback']) or 'none'}). "
            + ((f"Final MAG set and guilds (2026-11-19): comparison_rollout/summary.json and final_guilds_summary.json "
                f"(final_set_guilds.py), Table S8 (make_table_s8.py). KOfam calls cover {FG['n_mags_with_calls']:,} of "
                f"{FG['n_mags']:,} MAGs" + (f"; {FG['n_mags_missing_calls']} co-assembly MAGs are still being annotated "
                "(annotation_runs/kofam_coasm38), so the guild counts will be refreshed." if FG["n_mags_missing_calls"]
                else ".") + " Tables 3–6, Figures 8, 9, S8 and S11 and the abundance, guild and co-occurrence text use "
                "this final set (correlational/final585/), and the metabolic models were rebuilt for it "
                "(reassembly_2026/modelling_2026/: one model per MPn-consumer and betaine-reducer species; Figure 7).")
               if FS else "Pending: KOfam/guild calls of the second-pass MAGs and the co-assembly MAGs.")))
    else:
      pending_comments.append((wrap(r3), w_r3.find(qn("w:r")),
        f"PRELIMINARY: second-pass numbers are a snapshot of {snap} ({N['n_complete']} of {N['n_samples']} samples; "
        f"source reassembly_2026/{nfile}; samples that fell back to metaSPAdes: {', '.join(N['meta_fallback']) or 'none'}). "
        "Regenerate with reassembly_2026/manuscript_numbers.py and make_r232_review_docx.py when the poplar run finishes. "
        "Dereplication, GTDB R232 classification and guild calls of the second-pass MAGs are pending; the guild results "
        "in the next paragraph are for the first-pass MAGs."))
    r_last = r3
    if final:
        r_last, _ = new_par_after(r3, body_pPr,
            f"The {T['new_kept']} second-pass MAGs ({S2['rollout_bacteria']} bacterial and {S2['rollout_archaea']} "
            f"archaeal) formed {S2['n_species']['rollout']} species clusters, {S2['archaeal_species']['rollout']} of them "
            f"archaeal (original MAGs {S2['archaeal_species']['original']}; first pass "
            f"{S2['archaeal_species']['first_pass']}). From any sample, they recovered "
            f"{S2['originals_recovered_setwide']['rollout']} of the {S2['n_mags']['original']} original dereplicated MAGs "
            f"(first pass: {S2['originals_recovered_setwide']['first_pass']}), including "
            f"{S2['lost30_recovered_setwide_rollout']} of the 30 original MAGs that had no first-pass MAG above 80% ANI. "
            f"For the {S2['shared_species_rollout_original']} species shared with the original MAGs, the best second-pass "
            f"MAG scored higher (completeness − 5 × contamination) in {S2['rollout_better_than_original']}, with a median "
            f"N50 ratio of {S2['rollout_n50_ratio_vs_original']:.2f}; against the first pass "
            f"({S2['shared_species_rollout_first_pass']} shared species) it scored higher in "
            f"{S2['rollout_better_than_first_pass']}, with {S2['rollout_n50_ratio_vs_first_pass']:.2f} times the N50. "
            f"The first pass spanned more species ({S2['n_species']['first_pass']}), but "
            f"{S2['first_pass_only_species_coassembly_only']} of the {S2['first_pass_only_species']} species it alone "
            f"recovered came only from the site co-assemblies, which the second pass did not repeat; of the "
            f"{S2['first_pass_per_sample_keep_mags']} per-sample first-pass MAGs that passed the same CheckM2 and GUNC "
            f"criteria, {S2['first_pass_per_sample_keep_mags_with_rollout_species']} belonged to a species that the "
            "second pass also recovered.", body_rpr, lead_tab=True)
    if FS:
        r_last, _ = new_par_after(r_last, body_pPr,
            f"Together with the {CO['kept']} co-assembly MAGs that met the same criteria ({CO['kept_hq']} high-quality; "
            f"{co_fp} of them among the 548 that the first pass had kept with CheckM), the final set comprised {F['mags']:,} "
            f"MAGs in {F['species']} species clusters, {F['hq_representatives']} of them with a high-quality representative "
            f"and {F['archaeal_species']} archaeal (Table S8). Of these species, {F['species_coassembly_only']} came only "
            f"from the co-assemblies and {F['species_rollout_only']} only from the second pass. The final set contained "
            f"{F['first_pass_species_in_final']} of the {F['first_pass_species']} first-pass species and recovered "
            f"{F['originals_recovered_setwide']} of the {S2top['n_mags']['original']} original dereplicated MAGs from any "
            f"sample, including {F['lost30_recovered_setwide']} of the 30 that had no first-pass MAG above 80% ANI.",
            body_rpr, lead_tab=True)
        g, rm = FG["guilds"], FG["rast_mepn"]
        bet = " and ".join(f"{nw(n)} {ph}" for ph, n in FG["betaine_species_by_phylum"].items())
        rh = ("all Rhodobacteraceae" if FG["rhodobacteraceae_with_mttB_family"] == FG["rhodobacteraceae_mags"] else
              f"{FG['rhodobacteraceae_with_mttB_family']} of the {FG['rhodobacteraceae_mags']} Rhodobacteraceae MAGs")
        meth = ", ".join(f"*{k}*" for k in FG["methanogen_genera"])
        idt = pd.read_csv(RF / "iterativeIDs_final585.csv")
        bm4 = idt[idt.originals_in_cluster_r232.fillna("").str.split().apply(lambda xs: "BM714.4" in xs)].iterativeID
        bm4_txt = f" (now in species {bm4.iloc[0]})" if len(bm4) and bm4.iloc[0] != "BM714.4" else ""
        alb_bm = (" No assembly yielded a *phnJ* passing the KOfam threshold for *Albimonas*.1, and the betaine reductase "
                  f"annotation of the original MAG BM714.4{bm4_txt} was not supported: KOfam assigned its genes to the "
                  "glycine reductase subunits GrdE and GrdB." if FG["albimonas"]["mepn_mags"] == 0 and
                  FG["bm714_4_cluster"]["betaine_call"] is False else "")
        new_par_after(r_last, body_pPr,
            f"In the final MAG set (Table 4), KOfam assigned *phnJ* to {g['mepn_consumer']['mags']} MAGs in "
            f"{g['mepn_consumer']['species']} species, {nw(FG['mepn_species_rhodobacteraceae'])} of them Rhodobacteraceae "
            f"and {nw(FG['mepn_species_haloarchaea'])} haloarchaea. Betaine reductase (*grdH* or *grdI*) was found in "
            f"{g['betaine_reducer']['mags']} MAGs of {g['betaine_reducer']['species']} species ({bet}); reading through the "
            f"selenocysteine codon raised the number of *grdH* genes from {FG['grdH_genes_prodigal_only']} to "
            f"{FG['grdH_genes_after_rescue']}. *mcrA* was found only in {meth} ({nw(g['methanogen']['mags'])} MAGs, "
            f"{nw(g['methanogen']['species'])} species), which also carried the methylamine methyltransferases. No MAG "
            f"carried a complete TMA-producing gene pair: *cutC* occurred without *cutD* ({FG['tma_cutC_without_cutD_mags']} "
            f"MAGs) and *cntA* without *cntB* ({FG['tma_cntA_without_cntB_mags']} MAGs). Genes of the pyrrolysine-free MttB "
            f"family (K14083) occurred in {FG['mttB_family_without_mcrA_mags']} MAGs without *mcrA*, including {rh}; they "
            f"are not evidence of methylamine-dependent methanogenesis. Of the {rm['rast_mags']} MPn-consumer MAGs "
            f"identified by RAST among the original MAGs ({rm['rast_species']} species), the species of "
            f"{rm['recovered_with_call']} were recovered with *phnJ*, and {rm['further_final_species_without_original']} "
            f"further MPn-consumer species without an original MAG were found.{alb_bm}", body_rpr, lead_tab=True)
    else:
      new_par_after(r_last, body_pPr,
        "In the 1,324 first-pass MAGs, KOfam assigned *phnJ* to 87 MAGs in 35 species, 19 of them Rhodobacteraceae and "
        "three haloarchaea. Betaine reductase (*grdH* or *grdI*) was found in 53 MAGs of 16 species (15 Bacillota and one "
        "Cloacimonadota); reading through the selenocysteine codon raised the number of *grdH* genes from 9 to 54. *mcrA* "
        "was found only in *Methanolobus* (three MAGs, one species), which also carried the methylamine "
        "methyltransferases. No MAG carried a complete TMA-producing gene pair: *cutC* occurred without *cutD* (13 MAGs) "
        "and *cntA* without *cntB* (49 MAGs). Genes of the pyrrolysine-free MttB family (K14083) occurred in 328 MAGs "
        "without *mcrA*, including all Rhodobacteraceae; they are not evidence of methylamine-dependent methanogenesis. "
        "Of the 22 MPn-consumer MAGs identified by RAST among the original MAGs, the species of 13 were recovered with "
        "*phnJ*, and 19 further MPn-consumer species were found. No assembly yielded a *phnJ* passing the KOfam threshold "
        "for *Albimonas*.1, and the betaine reductase annotation of BM714.4 was not supported: KOfam assigned its genes to "
        "the glycine reductase subunits GrdE and GrdB.", body_rpr, lead_tab=True)
    # ---- Table S7
    cap, _ = new_par_after(t6, scap_pPr,
        f"Table S7. Second-pass re-assembly (merged reads, standard-mode SPAdes) compared with the first pass and with the "
        + ("original MAGs, for all 39 samples. Kept:" if final else
           f"original MAGs, for the {N['n_complete']} samples completed by {snap} (of 39). Kept:") + f" CheckM2 ≥ 50% complete, "
        f"< 10% contaminated and GUNC pass; HQ: also ≥ 90% complete and < 5% contaminated; all three sets scored with the "
        f"same CheckM2 v1.1.0 and GUNC v1.1.1. N50: median of the kept MAGs. Recovered: original MAGs matched by a "
        f"second-pass kept MAG of the same sample (skani ANI ≥ 95%, aligned fraction ≥ 50%).", scap_rpr)
    def cell(a, b=None):
        return "n/a" if a is None else (f"{a} ({b})" if b is not None else f"{a}")
    rows = [[r["sample"], cell(r["new_kept"], r["new_hq"]), kb(r["new_med_n50"]),
             cell(r["kbase2026_kept"], r["kbase2026_hq"]), kb(r["kbase2026_med_n50"]),
             f"{r['orig2020_mags']}", cell(r["orig2020_kept"], r["orig2020_hq"]), kb(r["orig2020_med_n50"]),
             cell(r["orig2020_recovered"])] for r in N["rows"]]
    rows.append(["Total", cell(T["new_kept"], T["new_hq"]), kb(N["pooled_new_med_n50"]), cell(T["kbase2026_kept"], T["kbase2026_hq"]),
                 kb(N["pooled_kbase2026_med_n50"]), f"{T['orig2020_mags']}", cell(T["orig2020_kept"], T["orig2020_hq"]),
                 kb(N["pooled_orig2020_med_n50"]), cell(T["orig2020_recovered"])])
    t7 = build_table(tblPr_src, tcMar_src,
                     ["Sample", "Second pass: kept (HQ)", "N50 (kb)", "First pass: kept (HQ)", "N50 (kb)",
                      "Original MAGs", "Original: kept (HQ)", "N50 (kb)", "Originals recovered"],
                     rows, [1500, 1050, 700, 1050, 700, 900, 1050, 700, 1000], body_rpr, size=16)
    cap.addnext(t7)
    if FS:
        new_par_after(t7, scap_pPr,
            f"Table S8. The final re-assembled MAG set: the {T['new_kept']} second-pass MAGs and the {CO['kept']} kept MAGs "
            "of the four first-pass site co-assemblies, with CheckM2 and CheckM completeness and contamination, GUNC result, "
            "GTDB R232 classification, species cluster (skani ANI ≥ 95%, aligned fraction ≥ 50%), its representative and "
            "label, KOfam guild calls, and, for each species representative, the mean relative abundance by site type and "
            "the correlation of its abundance with methane flux (separate file: Table S8 - re-assembled MAG set taxonomy, "
            "quality and guilds (GTDB R232).xlsx).", scap_rpr)


WORDS = {0: "no", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine",
         10: "ten", 11: "eleven", 12: "twelve"}


def _w(n):
    return WORDS.get(n, str(n))


def _tx(t):
    return f"*{t}*" if LATIN.match(t) else t


def _and(xs):
    xs = list(xs)
    return xs[0] if len(xs) == 1 else ", ".join(xs[:-1]) + " and " + xs[-1]


def final585_network_facts(g, SUM, NS):
    """Derived numbers and taxon lists for the final-set co-occurrence, Discussion and Figure S11 text; every list is
    computed from the correlational/final585 outputs, and the wording assumptions are asserted."""
    g = g.copy()
    g["taxon"] = g.iterativeID_final585.str.rsplit(".", n=1).str[0]
    kingdom = g.groupby("taxon").domain.first()
    meta = set(json.load(open(RF / "within_condition_meta_final585.json"))["taxa"])
    mods = json.load(open(RF / "within_condition_louvain_modules_final585.json"))
    vk = pd.read_csv(RF / "network_sensitivity_final585" / "variance_key_taxa.csv")
    vk = vk[vk.guilds.notna()]
    deg = vk.set_index("taxon").degree_gt005
    gt = {c: set(g.loc[g[c], "taxon"]) for c in ("mepn_consumer", "betaine_reducer", "methanogen")}
    panel = (gt["mepn_consumer"] | gt["betaine_reducer"] | gt["methanogen"]) & meta
    with_e = {t for t in panel if deg.get(t, 0) > 0}
    noedge = panel - with_e
    order = lambda ts: sorted(ts, key=lambda t: (-deg.get(t, 0), t))
    mepn_e = order(with_e & gt["mepn_consumer"] - gt["betaine_reducer"])
    bet_e = order(with_e & gt["betaine_reducer"])
    named = [t for t in ("Roseovarius", "Marivita", "Methanolobus", "Halanaerobium", "M55B127") if t in noedge]
    assert len(named) >= 3, named
    std = vk[vk.taxon.isin(noedge)].degree_gt005_standardised
    assert (std > 0).all(), "a guild taxon without edges has none after standardisation too"
    pf = NS["prefilter"]["min1"]
    gained = pf["key_taxa_gaining_edges"]
    assert not set(gained) & panel, "relaxing the filter gave an edge to a guild taxon in the panel"
    e = pd.read_csv(RF / "glasso_direct_edges_final585.csv"); e = e[e.partial_r.abs() > 0.05]
    base = {frozenset(x) for x in zip(e.taxon_a, e.taxon_b)}
    m1 = pd.read_csv(RF / "network_sensitivity_final585" / "prefilter_min1_edges.csv"); m1 = m1[m1.partial_r.abs() > 0.05]
    new = [(a, b) for a, b in zip(m1.taxon_a, m1.taxon_b) if frozenset((a, b)) not in base]
    below = sum((a not in meta) or (b not in meta) for a, b in new)
    gained_mepn = [t for t in gained if t in gt["mepn_consumer"]]
    assert set(gained_mepn) == set(gained), gained
    arch_net = {t for t in set(e.taxon_a) | set(e.taxon_b) if kingdom.get(t) == "Archaea"}
    top_mod = max(mods.values(), key=lambda m: sum(t in arch_net for t in m["members"]))
    arch_in_top = sum(t in arch_net for t in top_mod["members"])
    arch_mod = [arch_in_top]
    assert len(arch_net) == SUM["network_archaeal_taxa"]
    ab12 = pd.read_csv(RF / "averaged_normalized_MAG_abundances_final585.csv", index_col=0)
    nit = g[g.genus == "Nitrosopumilus"]
    det = ab12.loc[nit.mag] > 1e-5
    assert len(nit) and det.any().index[det.any()].str.startswith("R2A").all(), "Nitrosopumilus outside the reference site"
    ref = ab12.loc[nit.mag, [c for c in ab12.columns if c.startswith("R2A")]].mean(axis=1) * 100
    rel = pd.read_csv(S / "reassembly_2026" / "abundance_final" / "final585_relabund_24sediment.tsv", sep="\t", index_col=0)
    u = rel.loc["unbinned_reads"]
    both = g[g.mepn_consumer & g.betaine_reducer]
    s11_both = sorted(set(both.taxon) & meta)
    npair = SUM["pairs_ge_maxT_q95"]
    return {
        "n_maxT_words": "no pair" if npair == 0 else ("only one pair" if npair == 1 else f"only {_w(npair)} pairs"),
        "modules_word": _w(SUM["louvain_modules"]), "arch_net_word": _w(len(arch_net)), "arch_top_word": _w(arch_mod[0]),
        "arch_module_desc": (f"one module ({top_mod['size']} taxa, {_w(top_mod['size'] - arch_in_top)} of them bacterial)"
                             if top_mod["size"] > arch_in_top else "one haloarchaeal module"),
        "n_guild_panel": len(panel), "n_guild_edges_word": _w(len(with_e)), "n_guild_noedge": len(noedge),
        "mepn_edges": _and([_tx(t) for t in mepn_e]),
        "mepn_edges_note": " (all haloarchaea)" if all(kingdom[t] == "Archaea" for t in mepn_e) else "",
        "bet_edges": _and([_tx(t) for t in bet_e]), "noedge_named": _and([_tx(t) for t in named]),
        "std_deg_min": int(std.min()), "std_deg_max": int(std.max()),
        "pf_new": len(new), "pf_new_below": below,
        "pf_gained": f"edges of {_w(len(gained_mepn))} MPn-consumer taxa that the filter had excluded "
                     f"({_and([_tx(t) for t in sorted(gained_mepn)])})",
        "nitro_n_word": _w(len(nit)), "nitro_mags": int(nit.n_mags_in_species.sum()),
        "nitro_ref_pct": f"{ref.mean():.2f}%",
        "unbinned_by_site": {st: u[[c for c in u.index if c.split("_")[0] == st]].mean() * 100 for st in ("SF2", "R2A", "R2", "R1")},
        "s11_overlap_note": ("" if not s11_both else
                             f"; {_and([_tx(t) for t in s11_both])}, which includes a species carrying both markers, is "
                             "marked as a betaine reducer"),
    }


class Paras:
    """doc.paragraphs addressed by July-export index. The Google Doc of 2026-11-17 split the Abstract (P18) in two,
    so every index >= 19 sits one paragraph later in that mirror; raw(i) is the absolute index."""
    def __init__(self, ps):
        self.ps = ps
        self.shift = len(ps) - 281
        assert self.shift in (0, 1), f"unexpected paragraph count {len(ps)}"
        if self.shift:
            assert ps[19].text.startswith("Here we explore"), "expected the 2026-11-17 Abstract split at P19"

    def __getitem__(self, i):
        return self.ps[i if i <= 18 else i + self.shift]

    def raw(self, i):
        return self.ps[i]

    def __iter__(self):
        return iter(self.ps)

    def __len__(self):
        return len(self.ps)


def main(merged):
    doc = Document(merged)
    P = Paras(doc.paragraphs)
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
    corr = pd.read_csv(RF / "guild_methane_correlations_final585.csv")
    percore = pd.read_csv(RF / "methanogen_archaea_by_core_final585.csv")
    g = pd.read_csv(RF / "guilds_final585.csv")
    GC = ("mepn_consumer", "betaine_reducer", "methanogen", "tma_candidate_one_sided")
    guild = {c: sorted(g.loc[g[c], "iterativeID_final585"], key=lambda s: (s.rsplit(".", 1)[0], int(s.rsplit(".", 1)[1])))
             for c in GC}
    nsp = {c: int(g[c].sum()) for c in GC + ("archaea",)}
    nmg = {c: int(g.loc[g[c], f"{c}_n_mags"].sum()) for c in GC}
    SUM = json.load(open(RF / "summary_final585.json"))
    NS = json.load(open(RF / "network_sensitivity_final585" / "summary.json"))
    FG = json.load(open(S / "reassembly_2026" / "comparison_rollout" / "final_guilds_summary.json"))
    FN = final585_network_facts(g, SUM, NS)
    pending_comments = []   # (paragraph object, run element, text) -- added at the end

    # Keep handles to the paragraphs we edit (indices shift once we insert)
    p = {i: P[i] for i in (18, 30, 35, 47, 52, 53, 57, 69, 76, 77, 78, 83, 86, 95, 96, 99, 102, 113, 129, 132, 133, 134,
                           148, 149, 151, 163, 175)}
    first_par = next(x for x in P if x.text.strip())

    # ---- Methods ---------------------------------------------------------------------------------------------------
    h, _ = new_par_after(p[30]._p, h3_pPr, "*MAG quality, taxonomy, and functional guilds*", h3_rpr)
    mq_meth, _ = new_par_after(h, body_pPr,
        "The 310 dereplicated MAGs from Zhou et al. (2022) were re-assessed for comparison with the re-assembled MAG set "
        "described below. Completeness and "
        "contamination were estimated with CheckM v1.0.18 (67) (lineage workflow, reduced reference tree) in KBase (51), "
        "and the MAGs were classified with GTDB-Tk v2.8.0 (68) against GTDB release R232 (Table S5). Each MAG was "
        "labeled with its deepest named GTDB rank down to genus and numbered within that taxon (e.g., *Roseovarius*.3). "
        "Functional guilds were assigned from each MAG’s RAST annotation (53): methylphosphonate (MPn) consumers carry "
        "*phnJ*, the catalytic subunit of the C–P lyase (EC 4.7.1.1); betaine reducers carry *grdH*, the betaine-specific "
        "subunit of betaine reductase; and methanogens carry *mcrA*. MAGs with carnitine monooxygenase (*cntA*) or with "
        "genes of the choline TMA-lyase gene cluster were treated as candidate TMA producers; *cutC* itself was not "
        "annotated in any MAG.", body_rpr, lead_tab=True)
    # P35 (model building) is rewritten by the edits JSON (replace_paragraph, 2026-11-23)
    replace_paragraph(p[52]._p,
        "To test which organisms track the observed methane emissions, the relative abundance of each species of the "
        "re-assembled MAG set, and the summed abundance of each guild, were correlated with methane flux (Spearman’s ρ) "
        "across the 12 cores, with the D1 and D2 "
        "abundances merged by weighting with their DNA concentrations, and separately across the six unrestored cores. "
        "Correlations were computed on two bases: relative abundance versus areal CH~4~ flux, and both quantities divided "
        "by soil DNA concentration (the normalization of the model inputs). Because methane differs strongly among site types, each guild’s "
        "correlation was compared with those of 4,000 random sets of the same number of species, and the fraction of "
        f"random sets with an equal or higher ρ is reported. P values for the {SUM['n_species']} individual species were "
        "corrected with the Benjamini–Hochberg false discovery rate (FDR).", lead_tab=True)
    new_par_after(p[53]._p, body_pPr,
        "Co-occurrence among taxa was analyzed on centered log-ratio (CLR) transformed abundances of the 12 cores: only "
        "species present (relative abundance > 10^-5^) in at least six cores were retained and then summed per taxon "
        "(the taxon of their label), and "
        "pairs were tested only if both taxa were present together in at least six cores. Site "
        "means were subtracted so that the within-site residuals reflect coupling rather than a shared response to site "
        "conditions. Residual correlations (Pearson) were tested against all 1,296 within-site permutations with "
        "Benjamini–Hochberg FDR correction, a sparse network of direct (conditional) associations was estimated from the "
        "residuals with the graphical lasso (scikit-learn GraphicalLassoCV), and modules were identified by Louvain "
        "modularity on the positive partial-correlation edges (|partial r| > 0.05).", body_rpr, lead_tab=True)

    # ---- Results -------------------------------------------------------------------------------------------------------
    h, _ = new_par_after(p[57]._p, h3_pPr, "*MAG quality, taxonomy, and functional guilds*", h3_rpr)
    mq_res, _ = new_par_after(h, body_pPr,
        "All 310 dereplicated MAGs were at least 50% complete, and the re-run of CheckM reproduced the original "
        "completeness and contamination values exactly. Of the 310 MAGs, 275 met the MIMAG medium-quality criteria, while "
        "35 exceeded 10% contamination, including two of the *Roseovarius* MAGs (*Roseovarius*.11 "
        "and .12, 20–21% contamination; Table S5). GTDB R232 classified 280 MAGs as Bacteria and 30 as Archaea, and "
        "resolved 290 MAGs to genus and 176 to species, compared with 201 and 19 in the original classification. Most "
        "other differences were renamed groups (e.g., Proteobacteria are now Pseudomonadota); the genera central to the "
        "MPn hypothesis (*Roseovarius*, *Albimonas*, and *Marivita*) were unchanged, whereas one of the two *Halomonas* "
        "MAGs is now placed in *Vreelandella*. RAST assigned *phnJ* to 22 MAGs (Table S5): 11 of the 12 *Roseovarius* MAGs, "
        "two each of *Marivita* and *Desulfotignum*, single MAGs of five other bacterial genera, and two haloarchaea "
        "(*Natrinema* and *Haloquadratum*). Nine MAGs carried *grdH* (three *Halanaerobium*, two BM714, two M55B127, "
        "DATFJG01, and UBA8670), and a single MAG carried *mcrA*: *Methanolobus* sp020055655 (Methanosarcinaceae; 98% "
        "complete). The analyses below use the re-assembled MAG set and its KOfam-based guilds.",
        body_rpr, lead_tab=True)
    # P69 (Figure S8 citation and R² values), P74-P77 (betaine models): see the edits JSON
    hg, _ = new_par_after(p[77]._p, h3_pPr, "*Guild abundances and methane flux*", h3_rpr)
    c = corr.set_index(["guild", "basis", "cores"])
    v = lambda gname, basis, cores, k: c.loc[(gname, basis, cores), k]
    G5 = ["mepn_consumer", "betaine_reducer", "methanogen", "tma_candidate_one_sided", "archaea"]
    nf12 = [v(x, "raw", "12 cores", "null_frac_ge_rho") for x in G5]
    nfd = {x: v(x, "DNA-normalised", "12 cores", "null_frac_ge_rho") for x in G5}
    p6 = {(x, b): v(x, b, "6 unrestored", "p") for x in G5 for b in ("raw", "DNA-normalised")}
    # the wording below is conditional on these results (re-check the text if an assertion fails)
    assert min(nf12) > 0.05 and min(nfd.values()) > 0.05, "a guild now beats the random-set null"
    assert p6[("mepn_consumer", "DNA-normalised")] < 0.05 and all(pv >= 0.05 for k, pv in p6.items()
                                                                  if k != ("mepn_consumer", "DNA-normalised"))
    replace_paragraph(p[78]._p,
        f"Across the 12 cores, the summed abundance of MPn consumers was positively correlated with methane flux "
        f"(Spearman ρ = {v('mepn_consumer', 'raw', '12 cores', 'spearman_rho'):.2f}, p = "
        f"{fmt_p(v('mepn_consumer', 'raw', '12 cores', 'p'))}; Figure 9, Table 5), as were betaine reducers "
        f"(ρ = {v('betaine_reducer', 'raw', '12 cores', 'spearman_rho'):.2f}, p = {fmt_p(v('betaine_reducer', 'raw', '12 cores', 'p'))}), "
        f"the single methanogen species (ρ = {v('methanogen', 'raw', '12 cores', 'spearman_rho'):.2f}, p = "
        f"{fmt_p(v('methanogen', 'raw', '12 cores', 'p'))}), and the archaea as a whole (ρ = "
        f"{v('archaea', 'raw', '12 cores', 'spearman_rho'):.2f}, p = {fmt_p(v('archaea', 'raw', '12 cores', 'p'))}). These "
        f"correlations largely reflect the separation of the high-methane unrestored salterns from the other sites: "
        f"{SUM['per_species_q05_raw']} of the {SUM['n_species']} species were individually correlated with methane (FDR "
        f"q < 0.05), and no guild correlated more strongly than random sets of the same number of species "
        f"({min(nf12):.0%}–{max(nf12):.0%} of random sets did as well). The same held when abundances and fluxes were both "
        f"normalized by soil DNA (the normalization of the model inputs): {nfd['mepn_consumer']:.0%} of random sets matched "
        f"the correlation of MPn consumers (ρ = {v('mepn_consumer', 'DNA-normalised', '12 cores', 'spearman_rho'):.2f}) and "
        f"{nfd['betaine_reducer']:.0%} that of betaine reducers (ρ = "
        f"{v('betaine_reducer', 'DNA-normalised', '12 cores', 'spearman_rho'):.2f}; Figure S8). Within the six unrestored "
        f"cores, only the DNA-normalized abundance of MPn consumers was nominally correlated with methane (ρ = "
        f"{v('mepn_consumer', 'DNA-normalised', '6 unrestored', 'spearman_rho'):.2f}, p = "
        f"{fmt_p(p6[('mepn_consumer', 'DNA-normalised')])}, n = 6; "
        f"{v('mepn_consumer', 'DNA-normalised', '6 unrestored', 'null_frac_ge_rho'):.0%} of random sets did as well), "
        f"and on the areal basis no guild was (MPn consumers: ρ = {v('mepn_consumer', 'raw', '6 unrestored', 'spearman_rho'):.2f}, "
        f"p = {fmt_p(p6[('mepn_consumer', 'raw')])}).", lead_tab=True)
    pc = percore.set_index("core")
    unr = pc[pc.restoration == "unrestored"].methanogen_rel_abund * 100
    assert v("methanogen", "raw", "12 cores", "cores_detected") == 6 and (pc[pc.restoration != "unrestored"].methanogen_rel_abund <= 1e-5).all()
    p_meth, w_meth = new_par_after(p[78]._p, body_pPr,
        f"The only methanogen species, *Methanolobus* sp020055655 ({_w(nmg['methanogen'])} MAGs), was detected (relative "
        "abundance > 10^-5^) only in the six "
        f"unrestored cores, where it made up {unr.min():.2f}–{unr.max():.2f}% of reads (Table S6). Its correlation "
        "with methane across all 12 cores therefore reflects presence versus absence; within the unrestored cores its "
        f"abundance did not track methane (ρ = {num(v('methanogen', 'raw', '6 unrestored', 'spearman_rho'), '.2f')}, p = "
        f"{fmt_p(v('methanogen', 'raw', '6 unrestored', 'p'))}, n = 6).", body_rpr, lead_tab=True)
    hc, _ = new_par_after(p_meth, h3_pPr, "*Co-occurrence*", h3_rpr)
    cooc_res, _ = new_par_after(hc, body_pPr,
        f"With the re-assembled MAG set, the {SUM['within_species_kept']} species present in at least six cores formed "
        f"{SUM['within_taxa']} taxa ({SUM['within_archaeal_taxa']} of them archaeal) for the co-occurrence analysis. No "
        f"within-site correlation passed FDR correction, and {FN['n_maxT_words']} reached the 95th percentile of the "
        f"permutation maximum |r| ({SUM['maxT_q95']:.3f}), so the network edges are conditional associations rather than "
        f"individually significant pairs. The graphical-lasso network ({SUM['network_nodes']} taxa, {SUM['network_edges']} "
        f"edges) resolved {FN['modules_word']} modules (Q = {SUM['louvain_Q']:.3f}; Figure 8; Figure S11), and "
        f"{FN['arch_top_word']} of the {FN['arch_net_word']} archaeal taxa in the network fell in {FN['arch_module_desc']}. "
        f"Graphical-lasso fits to {NS['null']['n']} within-site permutations of the residuals gave denser networks "
        f"({NS['null']['edges_min']}–{NS['null']['edges_max']} nonzero edges, against {SUM['glasso_edges_all']} before the "
        f"|partial r| > 0.05 threshold) of similar modularity (Q = {NS['null']['Q_min']:.2f}–{NS['null']['Q_max']:.2f}); what "
        f"exceeded the permutations was that edges joined related taxa ({NS['observed']['same_phylum_frac']:.0%} of edges "
        f"joined taxa of the same phylum, against at most {NS['null']['same_phylum_frac_max']:.0%}). Of the "
        f"{FN['n_guild_panel']} guild taxa in the analysis, {FN['n_guild_edges_word']} had edges: the MPn consumers "
        f"{FN['mepn_edges']}{FN['mepn_edges_note']} and the betaine reducers {FN['bet_edges']}. The other "
        f"{FN['n_guild_noedge']}, including {FN['noedge_named']}, had none, mainly because of their low residual variance: "
        f"the lasso penalty acts on unstandardized residuals, whose variance tracked each taxon's number of edges (Spearman "
        f"ρ = {NS['variance']['spearman_resid_var_vs_degree']:.2f}), and with standardized residuals each of them had "
        f"{FN['std_deg_min']}–{FN['std_deg_max']} weak edges. Relaxing the prevalence filter gave none of them an edge; "
        f"{FN['pf_new_below']} of the {FN['pf_new']} edges it added involved a taxon below the filter, including "
        f"{FN['pf_gained']}.", body_rpr, lead_tab=True)

    # ---- Discussion and Conclusions ----------------------------------------------------------------------------------------
    w = tracked_insert_before(p[83]._p, "Concentrations of formate",
        f"In the re-assembled MAG set, the only *mcrA*-containing species (classified by GTDB R232 as *Methanolobus* "
        f"sp020055655, Methanosarcinaceae; {_w(nmg['methanogen'])} MAGs) was restricted to the unrestored salterns, but its "
        f"abundance did not track methane within them (ρ = {num(v('methanogen', 'raw', '6 unrestored', 'spearman_rho'), '.2f')}, "
        f"p = {fmt_p(v('methanogen', 'raw', '6 unrestored', 'p'))}, n = 6). ")
    pending_comments.append((p[83], w[0],
        "Naming: is this MAG the Ca. Methanosalis sp. SBSPR1A genome of ref. 58? GTDB R232 places it in Methanolobus "
        "(sp020055655); please reconcile the name used in the text."))
    tracked_insert_before(p[86]._p, "Here we propose that it also likely occurs",
        f"However, the {FN['nitro_n_word']} *Nitrosopumilus* species recovered here ({_w(FN['nitro_mags'])} MAGs) occurred only "
        f"in the reference wetland (about {FN['nitro_ref_pct']} of reads each) and not in the unrestored salterns, and no "
        "MAG of the re-assembled set had a KOfam hit to methylphosphonate synthase (*mpnS*; E ≤ 10^-5^), so the source "
        "of methylphosphonate in the unrestored salterns remains unidentified. ")
    tracked_replace(p[95]._p,
        "MPn-consuming organisms were much more associated with methane emissions than betaine-consuming organisms (R2 = 0.798).",
        "Neither MPn-consuming nor betaine-consuming organisms were more closely associated with methane emissions than "
        "random sets of community members: both guilds, the single methanogen species, and most of the community followed "
        "the same unrestored-versus-restored gradient, on an areal basis and with abundances and emissions normalized by "
        "soil DNA (Figure 9, Table 5).")
    tracked_replace(p[99]._p,
        "Secondly, the observed methane emissions were most correlated with the abundances of methylphosphonate-consuming "
        "microbes, more so than with the abundances of betaine-consuming microbes, which was our second hypothesized "
        "methanogenic substrate.",
        "Secondly, the abundances of MPn consumers (organisms carrying the C–P lyase gene *phnJ*) and of betaine reducers "
        "both followed methane emissions across sites, but no more closely than random sets of community members, because "
        "most of the community follows the same site gradient; within the unrestored salterns, only the DNA-normalized "
        "abundance of MPn consumers was nominally correlated with methane (n = 6).")
    p[102]._p.append(ins(runs_from(
        "MAG quality assessments, GTDB R232 classifications, guild assignments, and the correlation and co-occurrence "
        "analyses are in the MicrobiomeNotebooks repository (https://github.com/ModelSEED/MicrobiomeNotebooks; folders "
        "salterns/MAG_QC, salterns/reassembly_2026, salterns/correlational/r232 and salterns/correlational/final585).",
        first_rpr(p[102]._p))))
    unb_site = FN["unbinned_by_site"]
    tracked_replace(p[96]._p, "88% in SF2, 89% in R2A, 58% in R2, and 47% in R1",
                    f"{unb_site['SF2']:.0f}% in SF2, {unb_site['R2A']:.0f}% in R2A, {unb_site['R2']:.0f}% in R2, and "
                    f"{unb_site['R1']:.0f}% in R1 for the re-assembled MAG set")
    # P47 (Halomonas isolates, GTDB-Tk v1.7.0 -> v2.8.0 / R232): see the edits JSON
    pm = P.raw(19) if P.shift else p[18]            # the paragraph that holds the modeling clause
    pending_comments.append((pm, pm._p.findall(qn("w:r"))[0],
        "The modeling clause is now qualified (part of the measured methane, and no tracking of the differences among cores). If the abstract keeps a correlation claim, also "
        "note that in the re-assembled MAG set no guild's abundance tracks methane more closely than random sets of species, "
        "on the areal or the DNA-normalized basis, and that the single methanogen species does not track methane within "
        "the unrestored salterns (Table 5)."))

    # ---- Table 3 (replaced), Tables 4 and 5 (new) -------------------------------------------------------------------------
    unb = data["table3"]["unbinned"] if data else None
    w = replace_paragraph(p[113]._p,
        "Table 3. Percent of binned (MAG) relative abundance contributed by each GTDB R232 family in the unrestored (R1, R2), "
        "restored (SF2), and reference (R2A) wetlands, for the ten families most abundant in the unrestored salterns "
        f"(† archaeal family). Percentages are of the {SUM['n_species']} species of the re-assembled MAG set, including the "
        f"{SUM['n_archaeal_species']} archaeal species, and exclude unbinned reads (reads not assigned to any species), "
        f"which made up {unb['Unrestored (R1, R2)']:.0f}%, {unb['Restored (SF2)']:.0f}%, and {unb['Reference (R2A)']:.0f}% of "
        "reads in each site type, respectively. Because the table is sorted by the unrestored salterns, most of the restored "
        "and reference abundance falls in “Other families”.")
    pending_comments.append((p[113], w[0],
        "The previous Table 3 came from correlational/abundances.ipynb, which had three faults (fixed there on 2026-10-21): "
        "(1) its selection 'R1' in col or 'R2' in col also took the six R2A cores, so the Unrestored column pooled R1, R2 "
        "and R2A; (2) the 30 archaea had no 2020 taxonomy entry and were dropped with the unbinned reads, so the caption's "
        "68% is unbinned plus archaeal reads over those 18 cores (the data give 52%: 47% R1, 58% R2); (3) the notebook's "
        "column order was unrestored, reference, restored, which was copied under the headers Unrestored, Restored, "
        "Reference, swapping the last two. Recomputed here with GTDB R232 families for the 585 species of the re-assembled "
        "MAG set (CoverM read counts; correlational/final585/table3_final585.json); the corrected 2020-family values for the "
        "310 original MAGs are in 'Table 3 - corrected values (2020 GTDB families).csv'."))
    delete_table(t3)
    arch = set(data["table3"]["archaeal_in_top"])
    rows3 = [[f + ("†" if f in arch else ""), *[f"{x:.1f}" for x in vals]] for f, *vals in data["table3"]["rows"]]
    new_t3 = build_table(tblPr_src, tcMar_src, ["GTDB R232 family", "Unrestored (R1, R2)", "Restored (SF2)", "Reference (R2A)"],
                         rows3, [2880, 1890, 1860, 1890], body_rpr, size=22)
    sep, _ = new_par_after(t3, body_pPr, "", body_rpr)       # keeps Word from merging the deleted and new tables
    sep.addnext(new_t3)
    cap4, _ = new_par_after(new_t3, tcap_pPr,
        "Table 4. Functional guilds of the re-assembled MAG set, assigned from KOfam annotations; a species is in a guild "
        "when any of its MAGs carries the marker gene. Members are listed by species label (taxon.number; Table S8); a "
        "species that contains one of the 310 original MAGs keeps that MAG's label (Table S5).", tcap_rpr)
    orig_n = {c: int((g[c] & g.originals_in_cluster_r232.notna()).sum()) for c in GC}
    msum = pd.read_csv(S / "reassembly_2026" / "modelling_2026" / "models" / "models_summary.tsv", sep="\t")
    idf = json.load(open(RF / "iterativeIDs_final585.json"))
    MODN = {"mepn": int(((msum.guild == "mepn_consumer") & (msum.kept == True)).sum()),
            "bet": int(((msum.guild == "betaine_reducer") & (msum.kept == True)).sum()),
            "mepn_drop": ", ".join(idf[m] for m in msum[(msum.guild == "mepn_consumer") & (msum.kept != True)].mag_id)}
    assert MODN["bet"] == nsp["betaine_reducer"]
    rm, rb = FG["rast_mepn"], FG["rast_betaine"]
    rows4 = [
        ["MPn consumers", "*phnJ* (C–P lyase)", f"{nsp['mepn_consumer']} ({nmg['mepn_consumer']})",
         compress(guild["mepn_consumer"]),
         f"{orig_n['mepn_consumer']} contain an original MAG; {rm['recovered_with_call']} of the {rm['rast_species']} species "
         f"with a RAST *phnJ* call among the original MAGs are in this guild; {MODN['mepn']} modeled "
         f"({MODN['mepn_drop']} did not grow)"],
        ["Betaine reducers", "*grdH* or *grdI* (betaine reductase)", f"{nsp['betaine_reducer']} ({nmg['betaine_reducer']})",
         compress(guild["betaine_reducer"]),
         f"{orig_n['betaine_reducer']} contain an original MAG; {rb['recovered_with_call']} of the {rb['rast_species']} species "
         f"with a RAST *grdH* call among the original MAGs are in this guild; *grdH* counted after reading through the "
         f"selenocysteine codon; all {MODN['bet']} modeled"],
        ["Methanogens", "*mcrA*", f"{nsp['methanogen']} ({nmg['methanogen']})", compress(guild["methanogen"]) + " (sp020055655)",
         "Also carries the methylamine methyltransferases (*mttB*, *mtbB*, *mtmB*); no *mtaB*/*mtaC* (methanol). "
         "Modeled after curation (methanogenic on methylamines)"],
        ["TMA-route candidates", "*cutC* or *cntA* without its partner subunit", f"{nsp['tma_candidate_one_sided']} "
         f"({nmg['tma_candidate_one_sided']})", compress(guild["tma_candidate_one_sided"]),
         "No MAG carries *cutC* with *cutD* or *cntA* with *cntB*; not modeled"],
    ]
    t4 = build_table(tblPr_src, tcMar_src, ["Guild", "Marker gene", "Species (MAGs)", "Members", "Notes"],
                     rows4, [1350, 1650, 900, 2860, 1760], body_rpr, size=18,
                     align=["left", "left", "center", "left", "left"])
    cap4.addnext(t4)
    cap5, _ = new_par_after(t4, tcap_pPr,
        "Table 5. Spearman correlations (p value) between guild abundance and methane flux across the 12 cores and across "
        "the six unrestored cores, on an areal basis and with both quantities normalized by soil DNA, for the guilds of "
        "the re-assembled MAG set (Table 4). ‘Random ≥ ρ’ is the fraction of 4,000 random sets of the same number of "
        "species with an equal or higher 12-core ρ.", tcap_rpr)
    labels = [("mepn_consumer", f"MPn consumers (*phnJ*), {nsp['mepn_consumer']}"),
              ("betaine_reducer", f"Betaine reducers (*grdH*/*grdI*), {nsp['betaine_reducer']}"),
              ("methanogen", "Methanogen (*Methanolobus*.1)"),
              ("tma_candidate_one_sided", f"TMA-route candidates, {nsp['tma_candidate_one_sided']}"),
              ("archaea", f"All archaea, {nsp['archaea']}")]
    cell = lambda gname, basis, cores: f"{num(v(gname, basis, cores, 'spearman_rho'), '+.2f')} ({fmt_p(v(gname, basis, cores, 'p'))})"
    rows5 = [[lab, cell(gn, "raw", "12 cores"), cell(gn, "raw", "6 unrestored"), f"{v(gn, 'raw', '12 cores', 'null_frac_ge_rho'):.0%}",
              cell(gn, "DNA-normalised", "12 cores"), cell(gn, "DNA-normalised", "6 unrestored"),
              f"{v(gn, 'DNA-normalised', '12 cores', 'null_frac_ge_rho'):.1%}"] for gn, lab in labels]
    t5 = build_table(tblPr_src, tcMar_src, ["Guild, species", "Areal: 12 cores", "Areal: 6 unrestored", "Areal: random ≥ ρ",
                                            "DNA-norm.: 12 cores", "DNA-norm.: 6 unrestored", "DNA-norm.: random ≥ ρ"],
                     rows5, [2040, 1080, 1080, 1080, 1080, 1080, 1080], body_rpr, size=18)
    cap5.addnext(t5)

    # ---- Figure 8 (replaced), Figure 9 (new), Figure S8 (replaced) ----------------------------------------------------------
    # Figure 7: single panel on the areal basis, betaine = nine grdH carriers (caption P129 in the edits JSON)
    delete_paragraph(P[128]._p)
    picture_after(doc, P[128]._p, SIMF / "figure7_pergram.png", 6.5, img_pPr)
    # Figure S10: core gap-filling of the 2026 drafts (reassembly_2026/modelling_2026/figureS10/figureS10_final.py);
    # its caption (P157) is rewritten by the edits JSON
    delete_paragraph(P[156]._p)
    picture_after(doc, P[156]._p, S / "reassembly_2026" / "modelling_2026" / "figureS10" / "figureS10_final.png", 6.5, img_pPr)
    for i in (132, 133):
        delete_paragraph(p[i]._p)
    picture_after(doc, p[133]._p, RF / "within_condition_louvain_network_final585.png", 6.0, img_pPr)
    w = replace_paragraph(p[134]._p,
        "Figure 8. Co-occurrence network of the community members (re-assembled MAG set; GTDB R232 labels). Edges are "
        "partial correlations from a graphical lasso (L1-regularized precision matrix, regularization chosen by "
        "cross-validation) fit on within-site residuals of the CLR-transformed abundances of the 12 cores, for the "
        f"{SUM['within_taxa']} taxa formed from the {SUM['within_species_kept']} species present (relative abundance > "
        f"10^-5^) in at least six cores; edges with |partial r| > 0.05 are shown ({SUM['network_nodes']} taxa, "
        f"{SUM['network_edges']} edges). Node size encodes mean relative abundance, node color the phylum, and edge color "
        "the partial correlation (red = negative, blue = positive). Louvain community detection on the positive edges "
        f"resolves {FN['modules_word']} modules (Q = {SUM['louvain_Q']:.3f}), shown as shaded hulls; "
        f"{FN['arch_top_word']} of the {FN['arch_net_word']} archaeal taxa in the network fall in {FN['arch_module_desc']}. "
        "No individual within-site correlation passed FDR correction, so edges indicate conditional associations rather "
        "than demonstrated interactions; graphical-lasso fits to permuted residuals give denser networks of similar "
        "modularity, and only the excess of edges between related taxa exceeds them. Because the lasso penalty acts on "
        f"unstandardized residuals, taxa with low residual variance, including {FN['noedge_named']}, have no edges above "
        f"the threshold; the MPn consumers in the network are {FN['mepn_edges']}, and the betaine reducers are "
        f"{FN['bet_edges']}.")
    pending_comments.append((p[134], w[0],
        "Figure 8 is now a single panel: the previous panel a (plain graphical-lasso layout) was not re-rendered; this "
        "Louvain view shows the same edges. Code: correlational/final585/rerun_correlational_final585.py (the same "
        "statistics as correlational/r232/rerun_correlational_r232.py, applied to the re-assembled MAG set); permutation, "
        "variance and prefilter checks: correlational/final585/network_sensitivity_final585.py."))
    p9 = picture_after(doc, p[134]._p, RF / "guild_methane_raw_final585.png", 6.5, img_pPr)
    new_par_after(p9, cap_pPr,
        "Figure 9. Guild abundance versus methane flux across the 12 cores for the species of the re-assembled MAG set, "
        "with guilds assigned from KOfam annotations (a species is in a guild when any of its MAGs carries the marker): MPn "
        "consumers (*phnJ*), betaine reducers (*grdH* or *grdI*), the single methanogen species (*mcrA*, *Methanolobus* "
        "sp020055655), TMA-route candidates (*cutC* or *cntA* without its partner subunit), and all archaea. Abundance is "
        "the summed relative abundance (D1 and D2 merged with DNA weighting); methane is the areal CH~4~ flux on a "
        "symmetric-log axis because the reference wetland (R2A) showed net uptake. Each panel gives Spearman’s ρ for all "
        "12 cores and for the six unrestored cores, and the fraction of 4,000 random sets of the same number of species with "
        "an equal or higher 12-core ρ. Hollow markers are below detection (relative abundance ≤ 10^-5^).", cap_rpr)
    delete_paragraph(p[148]._p)
    picture_after(doc, p[148]._p, RF / "guild_methane_dnanorm_final585.png", 6.5, img_pPr)
    replace_paragraph(p[149]._p,
        "Figure S8. As Figure 9, with guild abundance and methane flux both normalized by soil DNA concentration (the "
        "normalization of the model inputs). On this basis, too, no guild correlates with methane more strongly than "
        "random sets of the same number of species (MPn consumers: ρ = "
        f"{v('mepn_consumer', 'DNA-normalised', '12 cores', 'spearman_rho'):.2f}, "
        f"{nfd['mepn_consumer']:.0%} of random sets did as well; Table 5).")
    delete_paragraph(p[151]._p)
    # Figure S11: regenerated with R232 labels (../correlational/r232/figureS11_r232.py)
    delete_paragraph(P[162]._p)
    picture_after(doc, P[162]._p, RF / "figureS11_final585.png", 6.5, img_pPr)
    for m in p[163]._p.findall(M_NS + "oMath"):            # the old caption's 'ρ=0.7' equation object
        d = stamp(el("w:del")); m.addprevious(d); d.append(m)
    w = replace_paragraph(p[163]._p,
        f"Figure S11. Within-site co-variation of prevalent taxa. Relative abundances of the {SUM['n_species']} species of "
        "the re-assembled MAG set were merged across the two depths of each core, weighted by DNA concentration (12 "
        f"cores). Species above 10^-5^ relative abundance in at least 6 of the 12 cores ({SUM['within_species_kept']} "
        f"species) were summed by GTDB R232 taxon ({SUM['within_taxa']} taxa, {SUM['within_archaeal_taxa']} archaeal, "
        "marked †), CLR-transformed, and the mean of each site (R1, R2, SF2, R2A; three cores each) was subtracted, which "
        "removes co-variation caused by shared responses to site. Colors show Pearson correlations between these "
        "within-site residuals, with rows and columns ordered by hierarchical clustering (average linkage). Green "
        f"asterisks mark the {SUM['network_edges']} pairs joined by an edge of the graphical-lasso network (|partial r| > "
        "0.05; Figure 8). No pair passed the permutation-based false discovery rate threshold, so the heatmap is "
        "descriptive. Labels and boxed rows/columns mark methylphosphonate consumers (*phnJ*; orange), betaine reducers "
        f"(*grdH* or *grdI*; blue{FN['s11_overlap_note']}), and the methanogen *Methanolobus* (*mcrA*; purple).")
    pending_comments.append((p[163], w[0],
        "Regenerated for the re-assembled MAG set with the KOfam-based guilds (correlational/final585/figureS11_final585.py, "
        "the renderer of correlational/r232/figureS11_r232.py), from the final585 within-site matrix and graphical-lasso "
        "edges; no statistics were recomputed. The previous caption described "
        "the asterisks as correlations above ρ = 0.7; in the code they marked graphical-lasso edges, as stated now. "
        "Methanolobus is detected only in the six unrestored cores, so much of its row reflects the cores where it is "
        "absent (CLR of below-detection values). The old caption's ρ = 0.7 equation object is a tracked deletion; "
        "LibreOffice ignores revisions on equations and still shows it."))

    # ---- Tables S5 (external) and S6 -----------------------------------------------------------------------------------------
    s5, _ = new_par_after(p[175]._p, scap_pPr,
        "Table S5. The 310 original dereplicated MAGs (Zhou et al. 2022): GTDB R232 classification, CheckM completeness and "
        "contamination, functional guild assignment with the supporting RAST annotation, and mean relative abundance by "
        "site type and correlation with methane flux in the original abundance table (separate file: Table S5 - MAG "
        "taxonomy, quality and guilds (GTDB R232).xlsx).", scap_rpr)
    s6, _ = new_par_after(s5, scap_pPr,
        "Table S6. Methane flux and the relative abundance (% of reads; D1 and D2 merged with DNA weighting) of the "
        f"methanogen species (*Methanolobus* sp020055655; {_w(nmg['methanogen'])} MAGs) and of all {nsp['archaea']} archaeal "
        "species of the re-assembled MAG set in each core.", scap_rpr)
    rows6 = [[r.core, r.restoration, num(r.CH4_umol_m2_d, ".1f"), f"{r.methanogen_rel_abund * 100:.4f}", f"{r.archaea_rel_abund * 100:.2f}"]
             for r in percore.sort_values(["restoration", "core"], key=lambda s: s.map({"unrestored": 0, "restored": 1, "reference": 2}) if s.name == "restoration" else s).itertuples()]
    t6 = build_table(tblPr_src, tcMar_src, ["Core", "Site type", "CH~4~ flux (µmol m^-2^ d^-1^)", "*Methanolobus*.1 (%)", "All archaea (%)"],
                     rows6, [1100, 1400, 2500, 1820, 1700], body_rpr, size=20)
    s6.addnext(t6)

    # ---- 2026-10-26 edits: modelling text and numbers, betaine set, Halomonas isolates, smaller corrections ------------------
    apply_edits(P, pending_comments, EDITS_1117 if P.shift else EDITS)
    # equations (Office Math): the unit conversion (P39) and the three regression fits (P70, P72, P75) are tracked
    # deletions (per-gram framing, 2026-11-24; see the P73 comment in r232_review_edits_2026-11-24.json)
    global DATE
    DATE, keep = DATE_MODELS, DATE
    for i in (39, 70, 72, 75):
        maths = P[i]._p.findall(".//" + M_NS + "oMath")
        assert maths, f"no equation in P{i}"
        for m in maths:
            d = stamp(el("w:del")); m.addprevious(d); d.append(m)
        delete_paragraph(P[i]._p)
    DATE = keep

    # ---- 2026-11-16: re-assembly, re-binning and annotation of the metagenomes (reassembly_2026/) ---------------------------
    DATE, keep_date = DATE_REASSEMBLY, DATE
    add_reassembly_2026(mq_meth, mq_res, t6, pending_comments, h3_pPr, h3_rpr, body_pPr, body_rpr, scap_pPr, scap_rpr,
                        tblPr_src, tcMar_src, lambda el: Paragraph(el, first_par._parent))
    DATE = keep_date

    # ---- comments, then save --------------------------------------------------------------------------------------------------
    pending_comments.insert(0, (first_par, first_par._p.findall(qn("w:r"))[0],
        "R232 review (2026-10-21, revised 2026-10-26, rebased 2026-11-04 on the Google Doc of 2026-11-03, whose open "
        "suggestions by Susannah Tringe appear as her tracked changes): tracked changes add MAG quality, taxonomy (GTDB R232), and gene-based "
        "guild methods and results; replace the abundance–methane results, Table 3, Figure 8, and Figure S8; and add "
        "Figure 9, Tables 4–5, S5 (Excel), and S6. The 2026-10-26 revision corrects the modeling text and numbers: on a "
        "dimensionally consistent areal basis the simulated methane does not track the "
        "differences among the six unrestored cores (the earlier high fits came from a 1/DNA term in the model input); it also restricts the betaine models to the nine grdH carriers, replaces "
        "Figure 7 with a single panel, and reclassifies the two Halomonas isolates with GTDB R232. A 2026-11-08 revision "
        "repaired three models of the 2020 set that needed methylphosphonate only because of gapfilling gaps; it is "
        "superseded by the 2026-11-23 models below. The 2026-11-19 "
        "revision adds the re-assembly, re-binning and KOfam annotation of the metagenomes (Methods, Results, Tables S7–S8). "
        "The 2026-11-20 revision (rebased on the Google Doc of 2026-11-20) moves Tables 3–6, Figures 8, 9, S8 and S11 and the "
        "abundance, guild and co-occurrence text to the re-assembled MAG set (585 species); on that set no guild tracks "
        "methane more closely than random sets of species on either basis, so the DNA-normalized MPn claim of the earlier "
        "revision is withdrawn. The 2026-11-23 revision replaces the metabolic models with new ones built for every "
        "MPn-consumer and betaine-reducer species of that set (curated ModelSEED templates, RAST annotation chosen by a "
        "pre-registered comparison); Figure 7, the modeling Methods and Results, and the modeling sentences of the Abstract, "
        "Discussion and Conclusions describe them, and the earlier repair of three MPn-obligate models no longer applies "
        "(no new model needs MPn). The 2026-11-24 revision replaces the model-vs-measured comparison, which needed an assumed "
        "growth rate (in effect 1 h⁻¹), with methane per gram of biomass formed and the growth each guild would need "
        "(Figure 7), adds the curated Methanolobus model, and removes the regression equations. Items flagged but not "
        "changed are in these comments; see the accompanying review notes."))
    for para, run_el, text in pending_comments:
        comment(doc, para, run_el, text)
    doc.save(OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main(sys.argv[1])
