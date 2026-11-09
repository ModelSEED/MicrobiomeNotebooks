#!/usr/bin/env python
"""
Local Word mirror of the "Metabolic Manuscript" Google Doc as of 2026-11-03, the base for the R232 review.

The Google Doc could not be exported as .docx through the Drive connector (the export is over the size limit), so the
mirror is rebuilt from the last export, ../Metabolic Manuscript.docx (2026-07-02), plus every change the Google Docs API
reports since then (gdoc_sync/gdoc_sync_2026-11-03.json):
  - 8 of Susannah Tringe's suggestions that A. Freiburger accepted on 2026-08-27 are accepted here;
  - her 37 suggestions of 2026-11-03 (22 edits here) are added as Word tracked changes under her name and date;
  - the comment thread resolved on 2026-08-27 is removed, and the 13 threads and 1 reply added since July are inserted
    with their authors and dates.
  - the 4 threads added on 2026-11-05 are inserted too (the text was unchanged on that date; checked 2026-11-09).
The result is verified paragraph by paragraph: with all tracked changes rejected it must equal the Google Doc text
with its open suggestions rejected, and with all accepted it must equal the text with them accepted
(gdoc_sync/gdoc_views_2026-11-03.json, from the Docs API). Figures were checked separately (same 20 images, same sizes).

Usage (input = the July export after the docx skill's merge_runs.py):
    ~/Documents/py_venv/bin/python make_gdoc_mirror_docx.py merged.docx
Output: ../Metabolic Manuscript - Google Doc mirror (2026-11-03).docx
"""
import copy
import json
import re
import sys
from pathlib import Path

from docx import Document
from docx.oxml.ns import qn
from docx.text.run import Run

import make_r232_review_docx as mr

HERE = Path(__file__).resolve().parent
SPEC = json.load(open(HERE / "gdoc_sync" / "gdoc_sync_2026-11-03.json"))
VIEWS = json.load(open(HERE / "gdoc_sync" / "gdoc_views_2026-11-03.json"))
OUT = HERE.parent / "Metabolic Manuscript - Google Doc mirror (2026-11-03).docx"


def accept_revisions(P, by_paragraph):
    """Accept specific tracked changes (by paragraph index and w:id): unwrap w:ins, drop w:del."""
    for pi, ids in by_paragraph.items():
        p = P[int(pi)]._p
        hits = [e for e in p if e.tag in (qn("w:ins"), qn("w:del")) and e.get(qn("w:id")) in ids]
        assert {e.get(qn("w:id")) for e in hits} == set(ids), (pi, ids)
        for e in hits:
            if e.tag == qn("w:ins"):
                for c in list(e):
                    e.addprevious(c)
            e.getparent().remove(e)


def tracked_superscript(p, context, spans):
    """Tracked formatting change: make the characters at `spans` (offsets within `context`) superscript."""
    _, start = mr._locate(p, context)
    for a, b in spans:
        mr._cut_at(p, start + a); mr._cut_at(p, start + b)
    pos, targets = 0, []
    for r in mr._text_runs(p):
        n = len(r.find(qn("w:t")).text or "")
        if any(start + a <= pos and pos + n <= start + b for a, b in spans) and n:
            targets.append(r)
        pos += n
    assert "".join(r.find(qn("w:t")).text for r in targets) == "".join(context[a:b] for a, b in spans)
    for r in targets:
        rpr = r.find(qn("w:rPr"))
        if rpr is None:
            rpr = mr.el("w:rPr"); r.insert(0, rpr)
        old = copy.deepcopy(rpr)
        va = mr.el("w:vertAlign", **{"w:val": "superscript"})
        anchor = next((c for c in rpr if c.tag in (qn("w:rtl"), qn("w:cs"), qn("w:em"), qn("w:lang"))), None)
        anchor.addprevious(va) if anchor is not None else rpr.append(va)
        chg = mr.stamp(mr.el("w:rPrChange"))
        chg.append(old); old.tag = qn("w:rPr")
        rpr.append(chg)


def remove_comments(doc, ids):
    body = doc.element.body
    for tag in ("w:commentRangeStart", "w:commentRangeEnd"):
        for e in list(body.iter(qn(tag))):
            if e.get(qn("w:id")) in ids:
                e.getparent().remove(e)
    for ref in list(body.iter(qn("w:commentReference"))):
        if ref.get(qn("w:id")) in ids:
            run = ref.getparent(); run.getparent().remove(run)
    ce = doc.comments._comments_elm
    for c in list(ce):
        if c.get(qn("w:id")) in ids:
            ce.remove(c)


def comment_on(doc, para, anchor, text, author, initials, date, context=None):
    p = para._p
    if context:                          # a short anchor that recurs in the paragraph: place it inside its (unique) context
        _, cstart = mr._locate(p, context)
        start = cstart + context.index(anchor)
    else:
        _, start = mr._locate(p, anchor)
    first = mr._cut_at(p, start); mr._cut_at(p, start + len(anchor))
    pos, covered = 0, []
    for r in mr._text_runs(p):
        n = len(r.find(qn("w:t")).text or "")
        if start <= pos and pos + n <= start + len(anchor) and n:
            covered.append(r)
        pos += n
    c = doc.add_comment([Run(covered[0], para), Run(covered[-1], para)], text=text, author=author, initials=initials)
    c._comment_elm.set(qn("w:date"), date)


def reply_to(doc, parent_id, text, author, initials, date):
    """Google exports replies as separate comments over the parent's range; do the same."""
    body = doc.element.body
    s = next(e for e in body.iter(qn("w:commentRangeStart")) if e.get(qn("w:id")) == parent_id)
    end = next(e for e in body.iter(qn("w:commentRangeEnd")) if e.get(qn("w:id")) == parent_id)
    runs, on = [], False
    for e in s.getparent().iter():
        if e is s: on = True
        elif e is end: break
        elif on and e.tag == qn("w:r") and e.find(qn("w:t")) is not None:
            runs.append(e)
    para = next(p for p in doc.paragraphs if p._p is s.getparent())
    c = doc.add_comment([Run(runs[0], para), Run(runs[-1], para)], text=text, author=author, initials=initials)
    c._comment_elm.set(qn("w:date"), date)


def views(doc):
    out = []
    for p in doc.element.body.iter(qn("w:p")):
        rej, acc = [], []
        for e in p.iter(qn("w:t"), qn("w:delText"), qn("w:tab")):
            anc = {a.tag for a in e.iterancestors()}
            t = "\t" if e.tag == qn("w:tab") else (e.text or "")
            if qn("w:ins") not in anc: rej.append(t)
            if qn("w:del") not in anc: acc.append(t)
        out.append(("".join(rej), "".join(acc)))
    return out


def norm(s):
    s = s.replace("[IMG]", "").replace("[EQ]", "").replace(" ", " ").replace("\u000b", " ")
    return re.sub(r"\s+", " ", s).strip()


def main(merged):
    doc = Document(merged)
    P = doc.paragraphs
    accept_revisions(P, SPEC["accept_revisions"]["paragraph_ids"])
    mr.AUTHOR, mr.DATE = SPEC["suggestions"]["author"], SPEC["suggestions"]["date"]
    mr._ids = iter(range(5000, 6999))
    for e in SPEC["suggestions"]["edits"]:
        p = P[e["p"]]._p
        if e["kind"] == "replace":
            mr.tracked_replace(p, e["old"], e["new"])
        elif e["kind"] == "insert_before":
            mr.tracked_insert_before(p, e["old"], e["new"])
        elif e["kind"] == "superscript":
            tracked_superscript(p, e["old"], e["spans"])
    remove_comments(doc, set(SPEC["remove_comments"]["word_ids"]))
    for c in SPEC["add_comments"]:
        hits = [q for q in P if c["anchor"] in "".join((r.find(qn("w:t")).text or "") for r in mr._text_runs(q._p))]
        if "p" in c:                      # anchor text occurs in several paragraphs: use the one the Google Doc names
            hits = [q for q in hits if q._p is P[c["p"]]._p]
        assert len(hits) == 1, (c["id"], len(hits))
        if "context" in c:                # text around the comment in the Drive export, so a wrong "p" fails loudly
            text = "".join((r.find(qn("w:t")).text or "") for r in mr._text_runs(hits[0]._p))
            assert c["context"] in text and c["anchor"] in c["context"], (c["id"], c["context"])
        comment_on(doc, hits[0], c["anchor"], c["text"], c["author"], c["initials"], c["date"],
                   context=c.get("context") if c.get("locate_by_context") else None)
    for r in SPEC["add_replies"]:
        reply_to(doc, r["parent_word_id"], r["text"], r["author"], r["initials"], r["date"])
    # verify against the Google Docs API views
    V = views(doc)
    assert len(V) == len(VIEWS), (len(V), len(VIEWS))
    bad = [(g["i"], k) for g, v in zip(VIEWS, V) for k, gt, mt in (("base", g["base"], v[0]), ("proposed", g["proposed"], v[1]))
           if norm(gt) != norm(mt)]
    assert not bad, f"mirror differs from the Google Doc in {len(bad)} paragraph views: {bad[:10]}"
    doc.save(OUT)
    print(f"wrote {OUT.name}: {len(V)} paragraphs match the Google Doc with suggestions rejected and accepted; "
          f"{len(list(doc.comments))} comments")


if __name__ == "__main__":
    main(sys.argv[1])
