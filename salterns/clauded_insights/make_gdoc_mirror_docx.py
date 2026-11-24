#!/usr/bin/env python
"""
Local Word mirror of the "Metabolic Manuscript" Google Doc, the base for the R232 review. One mirror per sync date:
  2026-11-03 (default): gdoc_sync/gdoc_sync_2026-11-03.json + gdoc_views_2026-11-03.json (described below);
  2026-11-17: gdoc_sync/gdoc_sync_2026-11-17.json + gdoc_views_2026-11-17.json. On top of the 11-03 changes, A. Freiburger's
    edits of 2026-11-16/17: 27 of Susannah Tringe's suggestions accepted (24 from 2026-05-07 in July P23-P24, 3 from
    2026-11-03 at P8, P23, P24), direct rewrites of the Abstract (split into two paragraphs, so 492 paragraphs) and the
    first four Introduction paragraphs (P20-P23; the compatible-solute paragraph P22 now comes before P21), 3 threads
    resolved and 1 orphaned (Word ids 2, 3, 6, 8 removed). The new spec sections (accept_all_revisions,
    accepted_suggestions, direct_edits) are optional, so the 11-03 spec builds exactly as before.

The Google Doc could not be exported as .docx through the Drive connector (the export is over the size limit), so the
mirror is rebuilt from the last export, ../Metabolic Manuscript.docx (2026-07-02), plus every change the Google Docs API
reports since then (for 2026-11-03, gdoc_sync/gdoc_sync_2026-11-03.json):
  - 8 of Susannah Tringe's suggestions that A. Freiburger accepted on 2026-08-27 are accepted here;
  - her 37 suggestions of 2026-11-03 (22 edits here) are added as Word tracked changes under her name and date;
  - the comment thread resolved on 2026-08-27 is removed, and the 13 threads and 1 reply added since July are inserted
    with their authors and dates.
  - the 4 threads added on 2026-11-05 are inserted too (the text was unchanged on that date; checked 2026-11-09).
The result is verified paragraph by paragraph: with all tracked changes rejected it must equal the Google Doc text
with its open suggestions rejected, and with all accepted it must equal the text with them accepted
(gdoc_sync/gdoc_views_2026-11-03.json, from the Docs API). Figures were checked separately (same 20 images, same sizes).

Usage (input = the July export after the docx skill's merge_runs.py):
    ~/Documents/py_venv/bin/python make_gdoc_mirror_docx.py merged.docx [--date 2026-11-17] [--out path.docx]
Output: ../Metabolic Manuscript - Google Doc mirror (<date>).docx (date default 2026-11-03)
"""
import argparse
import copy
import json
import re
from pathlib import Path

from docx import Document
from docx.opc.constants import RELATIONSHIP_TYPE as RT
from docx.oxml.ns import qn
from docx.text.run import Run

import make_r232_review_docx as mr

HERE = Path(__file__).resolve().parent
DEFAULT_DATE = "2026-11-03"


def sync_paths(date):
    return (HERE / "gdoc_sync" / f"gdoc_sync_{date}.json", HERE / "gdoc_sync" / f"gdoc_views_{date}.json",
            HERE.parent / f"Metabolic Manuscript - Google Doc mirror ({date}).docx")


def _load(date):
    spec, views, out = sync_paths(date)
    return json.load(open(spec)), json.load(open(views)), out


SPEC, VIEWS, OUT = _load(DEFAULT_DATE)


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


COMMENT_MARKS = (qn("w:commentRangeStart"), qn("w:commentRangeEnd"), qn("w:commentReference"))


def accept_all_revisions(P, spec, removed_ids):
    """Accept every tracked change left in whole paragraphs (all suggestions there were accepted in the Google Doc).
    spec: {"<p>": {"author": ..., "elements": n}}. Comments anchored only inside deleted text must be ones removed."""
    for pi, want in spec.items():
        p = P[int(pi)]._p
        revs = [e for e in p if e.tag in (qn("w:ins"), qn("w:del"))]
        assert len(revs) == want["elements"], (pi, len(revs), want["elements"])
        assert {e.get(qn("w:author")) for e in revs} == {want["author"]}, (pi, {e.get(qn("w:author")) for e in revs})
        assert not [e for e in p.iter(qn("w:ins"), qn("w:del")) if e.getparent() is not p], pi
        inside = {m.get(qn("w:id")) for e in revs if e.tag == qn("w:del") for m in e.iter(*COMMENT_MARKS)}
        assert inside <= removed_ids, (pi, inside - removed_ids)
        for e in revs:
            if e.tag == qn("w:ins"):
                for c in list(e):
                    e.addprevious(c)
            e.getparent().remove(e)


def plain_edit(p, e):
    """An accepted suggestion: apply it as a tracked change under a placeholder author, then accept it."""
    keep = mr.AUTHOR
    mr.AUTHOR = "__accepted__"
    if e["kind"] == "replace":
        mr.tracked_replace(p, e["old"], e["new"])
    elif e["kind"] == "insert_before":
        mr.tracked_insert_before(p, e["old"], e["new"])
    else:
        raise ValueError(e["kind"])
    mr.AUTHOR = keep
    mine = [x for x in p if x.tag in (qn("w:ins"), qn("w:del")) and x.get(qn("w:author")) == "__accepted__"]
    assert mine, e
    for x in mine:
        if x.tag == qn("w:ins"):
            for c in list(x):
                x.addprevious(c)
        x.getparent().remove(x)


def _seg_rpr(tpl, seg):
    """Run properties for a Docs API segment: the paragraph's own run format plus italic/bold/sub-/superscript."""
    rpr = copy.deepcopy(tpl) if tpl is not None else mr.el("w:rPr")
    for tag in ("w:b", "w:bCs", "w:i", "w:iCs", "w:vertAlign"):
        for x in rpr.findall(qn(tag)):
            rpr.remove(x)
    head = [x for x in rpr if x.tag in (qn("w:rStyle"), qn("w:rFonts"))]
    at = rpr.index(head[-1]) + 1 if head else 0
    flags = []
    if seg.get("b"):
        flags += [mr.el("w:b", **{"w:val": "1"}), mr.el("w:bCs", **{"w:val": "1"})]
    if seg.get("i"):
        flags += [mr.el("w:i", **{"w:val": "1"}), mr.el("w:iCs", **{"w:val": "1"})]
    for k, x in enumerate(flags):
        rpr.insert(at + k, x)
    valign = "superscript" if seg.get("sup") else "subscript" if seg.get("sub") else "baseline" if seg.get("link") else None
    if valign:                            # the July export marks citation links vertAlign=baseline, before w:rtl
        va = mr.el("w:vertAlign", **{"w:val": valign})
        tail = next((x for x in rpr if x.tag in (qn("w:rtl"), qn("w:cs"), qn("w:em"), qn("w:lang"))), None)
        tail.addprevious(va) if tail is not None else rpr.append(va)
    return rpr


def segment_content(doc, segments, tpl):
    """Docs API segments (plain str, or {"t", "i"/"b"/"sup"/"sub"/"link"}) -> w:r / w:hyperlink elements.
    A tab starts its own run (w:tab before the w:t), so each run keeps a single w:t."""
    out = []
    for seg in segments:
        seg = {"t": seg} if isinstance(seg, str) else seg
        assert "\u000b" not in seg["t"] and "\n" not in seg["t"], seg
        runs = []
        for k, piece in enumerate(seg["t"].split("\t")):
            if k == 0 and not piece:
                continue
            r = mr.el("w:r"); r.append(_seg_rpr(tpl, seg))
            if k:
                r.append(mr.el("w:tab"))
            if piece:
                t = mr.el("w:t"); t.set(qn("xml:space"), "preserve"); t.text = piece; r.append(t)
            runs.append(r)
        if seg.get("link"):
            h = mr.el("w:hyperlink", **{"r:id": doc.part.relate_to(seg["link"], RT.HYPERLINK, is_external=True)})
            for r in runs:
                h.append(r)
            out.append(h)
        else:
            out.extend(runs)
    return out


def rewrite_paragraph(doc, para, segments, reanchor=()):
    """Direct (untracked) edit: replace the paragraph's content by the Docs API runs, keeping its paragraph format.
    Comments in the paragraph must all be re-anchored from the export markers (others were removed before)."""
    p = para._p
    have = {m.get(qn("w:id")) for m in p.iter(*COMMENT_MARKS)}
    want = {w for r in reanchor for w in r["word_ids"]}
    assert have == want, ("comments in rewritten paragraph", have, want)
    refs = {m.get(qn("w:id")): copy.deepcopy(m.getparent()) for m in p.iter(qn("w:commentReference"))}
    tpl = copy.deepcopy(mr.first_rpr(p))
    for c in list(p):
        if c.tag != qn("w:pPr"):
            p.remove(c)
    for x in segment_content(doc, segments, tpl):
        p.append(x)
    for r in reanchor:
        _, cstart = mr._locate(p, r["context"])
        start = cstart + r["context"].index(r["anchor"])
        mr._cut_at(p, start); mr._cut_at(p, start + len(r["anchor"]))
        pos, covered = 0, []
        for run in mr._text_runs(p):
            n = len(run.find(qn("w:t")).text or "")
            if start <= pos and pos + n <= start + len(r["anchor"]) and n:
                covered.append(run)
            pos += n
        assert "".join(x.find(qn("w:t")).text for x in covered) == r["anchor"], r
        after = covered[-1]
        for wid in r["word_ids"]:
            covered[0].addprevious(mr.el("w:commentRangeStart", **{"w:id": wid}))
            end = mr.el("w:commentRangeEnd", **{"w:id": wid})
            after.addnext(end); end.addnext(refs[wid]); after = refs[wid]
    return tpl


def apply_direct_edits(doc, P, edits):
    """Index-based rewrites first (July-export numbering), then paragraph inserts and moves."""
    tpls = {}
    for e in edits:
        if e["kind"] == "rewrite":
            tpls[e["p"]] = rewrite_paragraph(doc, P[e["p"]], e["segments"], e.get("reanchor_comments", ()))
    for e in edits:
        if e["kind"] == "insert_after":
            src = P[e["p"]]._p
            q = mr.el("w:p")
            q.append(copy.deepcopy(src.find(qn("w:pPr"))))
            for x in segment_content(doc, e["segments"], tpls.get(e["p"], mr.first_rpr(src))):
                q.append(x)
            src.addnext(q)
        elif e["kind"] == "move_before":
            P[e["before"]]._p.addprevious(P[e["p"]]._p)
        else:
            assert e["kind"] == "rewrite", e["kind"]


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


def main(merged, date=DEFAULT_DATE, out=None):
    global SPEC, VIEWS, OUT
    SPEC, VIEWS, OUT = _load(date)
    OUT = Path(out) if out else OUT
    doc = Document(merged)
    P = doc.paragraphs
    removed = set(SPEC["remove_comments"]["word_ids"])
    accept_revisions(P, SPEC["accept_revisions"]["paragraph_ids"])
    accept_all_revisions(P, SPEC.get("accept_all_revisions", {}).get("paragraphs", {}), removed)
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
    for e in SPEC.get("accepted_suggestions", {}).get("edits", []):
        plain_edit(P[e["p"]]._p, e)
    remove_comments(doc, removed)
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
    apply_direct_edits(doc, P, SPEC.get("direct_edits", {}).get("edits", []))
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
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("merged", help="the 2026-07-02 export after merge_runs.py")
    ap.add_argument("--date", default=DEFAULT_DATE, help="sync date: selects gdoc_sync/gdoc_{sync,views}_<date>.json")
    ap.add_argument("--out", help="output path (default ../Metabolic Manuscript - Google Doc mirror (<date>).docx)")
    a = ap.parse_args()
    main(a.merged, a.date, a.out)
