"""Turn the R232 review's tracked changes (author "Claude (R232 review)") into Google Docs API requests that make the
same changes as SUGGESTIONS in the Google Doc, mapped onto a fresh documents.get read.

    python review_to_suggestions.py <review.docx> <read_doc JSON> <out_dir> [--edits <edits JSON>]

Writes <out_dir>/plan.json (every operation with its API indices, the paragraph it belongs to and its section),
<out_dir>/hand_list.md (what cannot or must not go through the API: pictures, tables, comments, edits that touch an
open suggestion of someone else, paragraphs whose text does not match the Doc), and nothing else; requests for one
chunk are made by `chunk_requests(plan, lo, hi)` (all operations with lo <= index < hi, from the end backward).

Mapping. The review = the Doc mirror + Claude's changes; rejecting Claude's changes gives the mirror (validated), and
the mirror equals the Doc (validated per paragraph). So the review's paragraphs whose paragraph mark is not a Claude
insertion, in document order (table cells included, Claude's new tables excluded), are the Doc's paragraphs in order.
In each, the ORIGINAL text (unchanged + Claude-deleted + others' suggested-deleted characters) must equal the Doc's
base text (text runs without suggestedInsertionIds), character by character; each original character then has an API
index (UTF-16 units). A Claude deletion [s, e) becomes deleteContentRange over those indices; a Claude insertion before
original character k becomes insertText at the index of k (at the end: before the paragraph's newline), followed by
updateTextStyle for its italic / bold / superscript / subscript runs. Claude's new paragraphs become one insertText of
"\n" + their texts at the end of the preceding existing paragraph (headings get the named style of the existing heading
with the same Word style). In SUGGEST mode a deletion leaves the text in place, so only insertions shift indices; the
requests of a chunk are therefore sorted from the highest index down.
"""
import argparse
import json
import re
from pathlib import Path

from docx import Document

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
M = "{http://schemas.openxmlformats.org/officeDocument/2006/math}"
A_BLIP = "{http://schemas.openxmlformats.org/drawingml/2006/main}blip"
CL = "Claude (R232 review)"


def u16(s):
    return len(s.encode("utf-16-le")) // 2


def pmark(p):
    rpr = p.find(f"{W}pPr/{W}rPr")
    if rpr is None:
        return None, None
    for tag, kind in ((f"{W}ins", "ins"), (f"{W}del", "del")):
        e = rpr.find(tag)
        if e is not None:
            return kind, e.get(f"{W}author")
    return None, None


def fmt_of(r):
    rpr = r.find(f"{W}rPr")
    f = {"i": False, "b": False, "va": None}
    if rpr is not None:
        for tag, key in ((f"{W}i", "i"), (f"{W}b", "b")):
            e = rpr.find(tag)
            if e is not None and e.get(f"{W}val") not in ("0", "false"):
                f[key] = True
        v = rpr.find(f"{W}vertAlign")
        if v is not None and v.get(f"{W}val") in ("superscript", "subscript"):
            f["va"] = v.get(f"{W}val")
    return f


def tokens(p):
    """(char, state, fmt) in document order; state: base, cdel, cins, odel, oins. Also flags for pictures and maths."""
    out, flags = [], set()

    def walk(el, state):
        for c in el:
            tag = c.tag
            if tag in (f"{W}ins", f"{W}del"):
                mine = c.get(f"{W}author") == CL
                st = ("c" if mine else "o") + ("ins" if tag == f"{W}ins" else "del")
                if state in ("oins", "oins_c"):              # Claude's change inside another author's insertion:
                    st = "oins_c"                            # not in the Doc's base text; left to the authors
                walk(c, st)
            elif tag == f"{W}r":
                f = fmt_of(c)
                for x in c:
                    if x.tag in (f"{W}t", f"{W}delText"):
                        out.extend((ch, state, f) for ch in (x.text or ""))
                    elif x.tag == f"{W}tab":
                        out.append(("\t", state, f))
                    elif x.tag == f"{W}br":
                        if x.get(f"{W}type") != "page":      # a page break is an element, not text, in the Docs API
                            out.append(("\u000b", state, f))
                    elif x.tag == f"{W}drawing" or x.find(f".//{A_BLIP}") is not None:
                        flags.add("picture")
            elif tag in (f"{W}hyperlink", f"{W}smartTag", f"{W}customXml", f"{W}fldSimple", f"{W}sdt", f"{W}sdtContent"):
                walk(c, state)
            elif tag in (f"{M}oMath", f"{M}oMathPara"):
                flags.add("math")
    walk(p, "base")
    return out, flags


def api_paragraphs(doc):
    body = doc["tabs"][0]["documentTab"]["body"]["content"]
    out = []
    for el in body:
        if "paragraph" in el:
            out.append(el)
        elif "table" in el:
            for row in el["table"]["tableRows"]:
                for cell in row["tableCells"]:
                    out.extend(dict(ce, _table_end=el["endIndex"]) for ce in cell["content"] if "paragraph" in ce)
    return out


def api_chars(el):
    """Base characters of an API paragraph with their indices; and the set of indices that are suggested insertions."""
    base, sugg_ins, sugg_del = [], set(), set()
    for e in el["paragraph"]["elements"]:
        if "textRun" not in e:
            continue
        tr = e["textRun"]; t = tr["content"]; i = e["startIndex"]
        ins = bool(tr.get("suggestedInsertionIds")); dele = bool(tr.get("suggestedDeletionIds"))
        for ch in t:
            if ins:
                sugg_ins.add(i)
            else:
                base.append((ch, i))
                if dele:
                    sugg_del.add(i)
            i += u16(ch)
    if base and base[-1][0] == "\n":
        base = base[:-1]
    return base, sugg_ins, sugg_del


def norm(s):
    return s.replace(" ", " ")


def plan(review_path, api_path, edits_path=None):
    api = json.load(open(api_path)); api = api.get("content", api)
    AP = api_paragraphs(api)
    doc = Document(review_path)
    body = doc.element.body
    # review paragraphs in document order, with the table they sit in
    rev = []
    for el in body:
        if el.tag == f"{W}p":
            rev.append((el, None))
        elif el.tag == f"{W}tbl":
            rows = el.findall(f"{W}tr")
            new_tbl = rows and all(r.find(f"{W}trPr/{W}ins") is not None and
                                   r.find(f"{W}trPr/{W}ins").get(f"{W}author") == CL for r in rows)
            for p in el.iter(f"{W}p"):
                rev.append((p, "new_table" if new_tbl else "table"))
    existing = [(p, t) for p, t in rev if not (pmark(p) == ("ins", CL) or t == "new_table")]
    assert len(existing) == len(AP), (len(existing), len(AP))
    style_of = {}                                       # Word pStyle -> API namedStyleType (from existing paragraphs)
    for (p, _), a in zip(existing, AP):
        ps = p.find(f"{W}pPr/{W}pStyle")
        if ps is not None:
            style_of.setdefault(ps.get(f"{W}val"), a["paragraph"]["paragraphStyle"].get("namedStyleType"))
    ref_style = {}
    for a in AP:
        st_ = a["paragraph"]["paragraphStyle"].get("namedStyleType", "NORMAL_TEXT")
        for e in a["paragraph"]["elements"]:
            ts = e.get("textRun", {}).get("textStyle", {})
            if "weightedFontFamily" in ts and "fontSize" in ts and st_ not in ref_style and \
                    not e["textRun"].get("suggestedInsertionIds") and e["textRun"]["content"].strip():
                ref_style[st_] = {k: ts[k] for k in ("weightedFontFamily", "fontSize") if k in ts}
                ref_style[st_]["foregroundColor"] = ts.get("foregroundColor", {"color": {"rgbColor": {}}})
                break
    ops, hand, mismatched = [], [], []
    ai = {id(p): k for k, (p, _) in enumerate(existing)}
    section = "front"
    pending_new = []                                    # Claude's new paragraphs since the last existing one
    last_existing = None

    def flush_new():
        nonlocal pending_new
        if not pending_new:
            return
        k = last_existing
        a = AP[k]
        texts, styles, fmts = [], [], []
        for p, t in pending_new:
            tk, fl = tokens(p)
            txt = "".join(ch for ch, st, f in tk if st in ("cins", "base"))
            if t == "new_table":
                continue
            if "picture" in fl:
                hand.append({"what": "insert picture (new figure)", "after_api_paragraph": k, "section": section,
                             "context": api_text(a)[-120:]})
                continue
            if not txt.strip():
                continue
            ps = p.find(f"{W}pPr/{W}pStyle")
            styles.append(style_of.get(ps.get(f"{W}val")) if ps is not None else None)
            texts.append(txt)
            fmts.append([f for ch, st, f in tk if st in ("cins", "base")])
        anchor_style = a["paragraph"]["paragraphStyle"].get("namedStyleType", "NORMAL_TEXT")
        at, prefix = a["endIndex"] - 1, False
        if "_table_end" in a:                           # after a table: in front of the paragraph that follows it
            at, prefix = a["_table_end"], True
            nxt = next(p for p in AP if p["startIndex"] == at)
            anchor_style = nxt["paragraph"]["paragraphStyle"].get("namedStyleType", "NORMAL_TEXT")
        for n_, (t_, f_, sty) in enumerate(zip(texts, fmts, styles)):
            target = sty or "NORMAL_TEXT"
            ops.append({"kind": "new_paragraph", "index": at, "prefix": prefix, "seq": n_, "text": t_, "fmt": f_,
                        # a named style is applied only when the inherited (anchor) style differs: in SUGGEST mode it
                        # resets the run styles, so it goes before the text styles
                        "style": target if target != anchor_style else None,
                        "base_style": ref_style.get(target) or ref_style.get("NORMAL_TEXT"),
                        "api_paragraph": k, "section": section, "anchor_tail": api_text(a)[-80:]})
        pending_new = []

    tables_new = sum(1 for p, t in rev if t == "new_table")
    if tables_new:
        hand.append({"what": f"new or replaced tables ({tables_new} cell paragraphs): place by hand from the review"})
    for p, t in rev:
        if (p, t) not in existing and id(p) not in ai:
            if t != "new_table":
                pending_new.append((p, t))
            continue
        flush_new()
        k = ai[id(p)]; last_existing = k
        a = AP[k]
        txt_api = api_text(a)
        if a["paragraph"]["paragraphStyle"].get("namedStyleType", "").startswith("HEADING"):
            section = txt_api.strip()[:40] or section
        tk, fl = tokens(p)
        if not any(st in ("cins", "cdel") for _, st, _ in tk) and pmark(p) != ("del", CL):
            continue
        base, s_ins, s_del = api_chars(a)
        orig = [(ch, st, f) for ch, st, f in tk if st in ("base", "cdel", "odel")]
        o_txt = "".join(ch for ch, _, _ in orig); b_txt = "".join(ch for ch, _ in base)
        if norm(o_txt) != norm(b_txt):
            mismatched.append({"api_paragraph": k, "section": section, "review": o_txt[:160], "doc": b_txt[:160]})
            hand.append({"what": "text differs from the Doc: apply this paragraph's changes by hand", "api_paragraph": k,
                         "section": section, "context": b_txt[:120]})
            continue
        if t == "table":
            if not hand or hand[-1].get("what") != "existing table changed (replaced or deleted): by hand" or \
                    hand[-1].get("section") != section:
                hand.append({"what": "existing table changed (replaced or deleted): by hand", "api_paragraph": k,
                             "section": section, "cells": 0})
            hand[-1]["cells"] += 1
            continue
        if pmark(p) == ("del", CL):                     # whole paragraph deleted
            if "picture" in fl:
                hand.append({"what": "old figure picture: replace by hand", "api_paragraph": k, "section": section,
                             "context": b_txt[:80]})
                continue
            if s_ins or s_del:
                hand.append({"what": "delete paragraph (holds another suggestion): by hand", "api_paragraph": k,
                             "section": section, "context": b_txt[:120]})
                continue
            ops.append({"kind": "delete_paragraph", "index": a["startIndex"], "end": a["endIndex"], "api_paragraph": k,
                        "section": section, "text": b_txt[:200], "math": "math" in fl})
            continue
        # character-level operations
        idx = [i for _, i in base]
        end_index = a["endIndex"] - 1
        kpos = 0; cur_del = None; cur_ins = None
        para_ops = []
        for ch, st, f in tk:
            if st in ("base", "odel"):
                if cur_del: para_ops.append(cur_del); cur_del = None
                if cur_ins: para_ops.append(cur_ins); cur_ins = None
                kpos += 1
            elif st == "cdel":
                if cur_ins: para_ops.append(cur_ins); cur_ins = None
                if cur_del is None:
                    cur_del = {"kind": "delete", "s": kpos, "e": kpos}
                cur_del["e"] = kpos + 1; kpos += 1
            elif st == "cins":
                if cur_del: para_ops.append(cur_del); cur_del = None
                if cur_ins is None:
                    cur_ins = {"kind": "insert", "k": kpos, "text": "", "fmt": []}
                cur_ins["text"] += ch; cur_ins["fmt"].append(f)
            elif st == "oins":
                if cur_del: para_ops.append(cur_del); cur_del = None
                if cur_ins: para_ops.append(cur_ins); cur_ins = None
        for x in (cur_del, cur_ins):
            if x:
                para_ops.append(x)
        final = []
        for o in para_ops:
            o.update(api_paragraph=k, section=section)
            if o["kind"] == "insert":
                o.update(index=idx[o["k"]] if o["k"] < len(idx) else end_index)
                final.append(o)
                continue
            # a deletion: split it around characters that are someone else's suggested deletion, and around their
            # suggested insertions lying inside the range (those stay for the authors to resolve; noted in the hand list)
            lo_, hi_ = idx[o["s"]], idx[o["e"] - 1] + u16(base[o["e"] - 1][0])
            inside_ins = sorted(i for i in s_ins if lo_ <= i < hi_)
            pieces, start = [], None
            for kk in range(o["s"], o["e"]):
                i = idx[kk]
                ok = i not in s_del and not any(lo_ <= j < i for j in inside_ins if (start is not None and idx[start] < j < i))
                if ok and start is None:
                    start = kk
                if (not ok) and start is not None:
                    pieces.append((start, kk)); start = None
                if ok and start is not None and kk + 1 < o["e"] and idx[kk + 1] != i + u16(base[kk][0]):
                    pieces.append((start, kk + 1)); start = None     # an index gap (another run/element) splits it
            if start is not None:
                pieces.append((start, o["e"]))
            if len(pieces) != 1 or inside_ins or any(idx[kk] in s_del for kk in range(o["s"], o["e"])):
                hand.append({"what": "note: a replaced passage contains another author's open suggestion; it is left in "
                             "place (reject it, or resolve it with the replacement)", "api_paragraph": k,
                             "section": section, "context": b_txt[o["s"]:o["s"] + 100]})
            for a_, b_ in pieces:
                final.append(dict(o, s=a_, e=b_, index=idx[a_], end=idx[b_ - 1] + u16(base[b_ - 1][0]),
                                  text=b_txt[a_:b_]))
        ops.extend(final)
    flush_new()
    # comments of the review
    cm = []
    for c in doc.comments if hasattr(doc, "comments") else []:
        if c.author == CL:
            cm.append(c.text)
    if cm:
        hand.append({"what": f"{len(cm)} review comments: add by hand (the API cannot anchor comments)", "texts": cm})
    return {"revisionId": api["revisionId"], "ops": ops, "hand": hand, "mismatched": mismatched}


def api_text(el):
    return "".join(e["textRun"]["content"] for e in el["paragraph"]["elements"] if "textRun" in e)


def _style_requests(start, fmts):
    """One reset of the inserted range to plain text (it would otherwise inherit the style of the character before it),
    then one request per italic / bold / superscript / subscript run."""
    plain = {"i": False, "b": False, "va": None}
    reqs = [{"updateTextStyle": {"range": {"startIndex": start, "endIndex": start + len(fmts)},
                                 "textStyle": {"italic": False, "bold": False, "baselineOffset": "NONE"},
                                 "fields": "italic,bold,baselineOffset"}}] if fmts else []
    i, n = 0, len(fmts)
    while i < n:
        f = fmts[i]; j = i
        while j < n and fmts[j] == f:
            j += 1
        if f != plain:
            ts, fields = {}, []
            if f["i"]:
                ts["italic"] = True; fields.append("italic")
            if f["b"]:
                ts["bold"] = True; fields.append("bold")
            if f["va"]:
                ts["baselineOffset"] = {"superscript": "SUPERSCRIPT", "subscript": "SUBSCRIPT"}[f["va"]]
                fields.append("baselineOffset")
            reqs.append({"updateTextStyle": {"range": {"startIndex": start + i, "endIndex": start + j}, "textStyle": ts,
                                             "fields": ",".join(fields)}})
        i = j
    return reqs


def op_requests(o):
    if o["kind"] in ("delete", "delete_paragraph"):
        return [{"deleteContentRange": {"range": {"startIndex": o["index"], "endIndex": o["end"]}}}]
    if o["kind"] == "insert":
        assert u16(o["text"]) == len(o["text"]), "non-BMP characters"
        return [{"insertText": {"location": {"index": o["index"]}, "text": o["text"]}}] + _style_requests(o["index"], o["fmt"])
    if o["kind"] == "new_paragraph":                   # "\n" + text at the end of the anchor paragraph
        text = o["text"] + "\n" if o.get("prefix") else "\n" + o["text"]
        assert u16(text) == len(text)
        lo = o["index"] + (0 if o.get("prefix") else 1); hi = lo + len(o["text"])
        reqs = [{"insertText": {"location": {"index": o["index"]}, "text": text}}]
        if o["style"]:
            reqs.append({"updateParagraphStyle": {"range": {"startIndex": lo, "endIndex": hi},
                                                  "paragraphStyle": {"namedStyleType": o["style"]},
                                                  "fields": "namedStyleType"}})
        base = dict(o.get("base_style") or {}, italic=False, bold=False, baselineOffset="NONE")
        reqs.append({"updateTextStyle": {"range": {"startIndex": o["index"], "endIndex": hi + (1 if o.get("prefix") else 0)},
                                         "textStyle": base,
                                         "fields": ",".join(base)}})
        reqs += _style_requests(lo, o["fmt"])[1:]
        return reqs
    raise ValueError(o["kind"])


def chunk_requests(P, lo, hi):
    sel = [o for o in P["ops"] if lo <= o["index"] < hi]
    # highest index first; at equal index an insertion (it goes after a suggested deletion ending there) comes first
    # at one index: deletions first (they do not shift text), then new paragraphs, then the inline insertion, which
    # therefore ends up before the new paragraphs (each insertion at an index goes in front of what is already there)
    # new paragraphs at one index go in last-first (each lands in front of the ones already inserted)
    rank = {"delete": 3, "delete_paragraph": 3, "new_paragraph": 2, "insert": 1}
    sel.sort(key=lambda o: (o["index"], rank[o["kind"]], o.get("seq", 0)), reverse=True)
    reqs = []
    for o in sel:
        reqs += op_requests(o)
    return sel, reqs


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("review"); ap.add_argument("api"); ap.add_argument("out")
    a = ap.parse_args()
    P = plan(a.review, a.api)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    json.dump(P, open(out / "plan.json", "w"), indent=1, ensure_ascii=False)
    from collections import Counter
    print("ops:", Counter(o["kind"] for o in P["ops"]), "hand:", len(P["hand"]), "mismatched:", len(P["mismatched"]))
    print("inserted characters:", sum(len(o["text"]) + (o["kind"] == "new_paragraph") for o in P["ops"]
                                       if o["kind"] in ("insert", "new_paragraph")))
    print("sections:", Counter(o["section"] for o in P["ops"]))
