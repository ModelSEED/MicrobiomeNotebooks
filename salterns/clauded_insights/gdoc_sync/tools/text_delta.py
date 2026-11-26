#!/usr/bin/env python
"""Text-anchored suggestions for the Google Doc: apply a list of edits, written against the Doc's PROPOSED text (base +
suggested insertions, without suggested deletions), as SUGGEST-mode requests computed on a fresh read.

  text_delta.py plan   <read_doc.json> <ops.json> [--only key,key] > requests.json
  text_delta.py verify <read_doc.json> <ops.json> [--only ...]      (each op's 'expect' text is in its paragraph)

ops.json: a list of
  {"key", "anchor": text that occurs in exactly one paragraph's proposed text,
   "kind": "replace"     -- "old" (unique in that paragraph's proposed text) -> "new" (markup);
           "insert_before" -- insert "new" in front of "old";
           "replace_all" -- suggest deleting every base character of the paragraph (others' open suggestions are left
                            alone) and insert "new" at the end of the paragraph;
           "append_paragraph" -- a new paragraph "new" after the anchor paragraph (same paragraph style);
   "new": markup (*italic*, ^superscript^, ~subscript~), "expect": text to find afterwards (optional)}
Requests run from the end of the Doc backward, so every index is valid when its request runs; inserted text is reset to
Times New Roman 12 pt black before its italics and offsets are applied (inserted text otherwise inherits the style of the
character before it). Deleting text that is itself a suggested insertion (of the same account) removes it from the
index space and merges the new text into that suggestion (tested 2026-11-26, P96 "80%" -> "79%"); deleting base text
leaves it in place as a suggested deletion. Either way, ops sorted from the highest index down stay valid.
"""
import json
import re
import sys

TEXT = {"weightedFontFamily": {"fontFamily": "Times New Roman", "weight": 400}, "fontSize": {"magnitude": 12, "unit": "PT"},
        "foregroundColor": {"color": {"rgbColor": {}}}, "italic": False, "bold": False, "baselineOffset": "NONE"}
FIELDS = "weightedFontFamily,fontSize,foregroundColor,italic,bold,baselineOffset"
MARK = re.compile(r"(\*[^*]+\*|\^[^^]+\^|~[^~]+~)")


def u16(s):
    return len(s.encode("utf-16-le")) // 2


def load(p):
    d = json.load(open(p))
    return d.get("content", d)


def paragraphs(doc):
    out = []
    for e in doc["tabs"][0]["documentTab"]["body"]["content"]:
        if "paragraph" not in e:
            continue
        chars = []                                    # (index, ch, state) ; state base / ins / del / obj
        for x in e["paragraph"]["elements"]:
            if "textRun" in x:
                tr = x["textRun"]
                st = "ins" if tr.get("suggestedInsertionIds") else ("del" if tr.get("suggestedDeletionIds") else "base")
                i = x["startIndex"]
                for ch in tr["content"]:
                    chars.append((i, ch, st)); i += u16(ch)
            else:
                chars.append((x["startIndex"], "\x00", "obj"))
        out.append({"start": e["startIndex"], "end": e["endIndex"], "chars": chars})
    return out


def proposed(p):
    keep = [c for c in p["chars"] if c[2] in ("base", "ins")]
    return "".join(c[1] for c in keep), keep


def parse(markup):
    text, spans = "", []
    for tok in MARK.split(markup):
        if not tok:
            continue
        kind = None
        if tok[0] == tok[-1] == "*" and len(tok) > 2:
            tok, kind = tok[1:-1], "italic"
        elif tok[0] == tok[-1] == "^" and len(tok) > 2:
            tok, kind = tok[1:-1], "SUPERSCRIPT"
        elif tok[0] == tok[-1] == "~" and len(tok) > 2:
            tok, kind = tok[1:-1], "SUBSCRIPT"
        if kind:
            spans.append((u16(text), u16(text) + u16(tok), kind))
        text += tok
    return text, spans


def insert_requests(at, markup):
    text, spans = parse(markup)
    reqs = [{"insertText": {"text": text, "location": {"index": at}}},
            {"updateTextStyle": {"range": {"startIndex": at, "endIndex": at + u16(text)}, "textStyle": TEXT,
                                 "fields": FIELDS}}]
    for a, b, kind in spans:
        st = {"italic": True} if kind == "italic" else {"baselineOffset": kind}
        reqs.append({"updateTextStyle": {"range": {"startIndex": at + a, "endIndex": at + b}, "textStyle": st,
                                         "fields": "italic" if kind == "italic" else "baselineOffset"}})
    return reqs


def text_end(p):
    """Index just after the paragraph's last text character (so inserted text goes before a trailing page break)."""
    last = [c for c in p["chars"] if c[2] != "obj" and c[1] != "\n"][-1]
    return last[0] + u16(last[1])


def find_paragraph(paras, anchor):
    hits = [p for p in paras if anchor in proposed(p)[0]]
    assert len(hits) == 1, f"anchor in {len(hits)} paragraphs: {anchor[:70]!r}"
    return hits[0]


def op_requests(paras, op):
    p = find_paragraph(paras, op["anchor"])
    text, keep = proposed(p)
    kind = op["kind"]
    if kind == "replace":
        old = op["old"]
        assert text.count(old) == 1, f"{op['key']}: old occurs {text.count(old)}x"
        k = text.index(old)
        seg = keep[k:k + len(old)]
        a, b = seg[0][0], seg[-1][0] + u16(seg[-1][1])
        inner = [c for c in p["chars"] if a <= c[0] < b and c[2] in ("del", "obj")]
        assert not inner, f"{op['key']}: the old text spans a suggested deletion or an object"
        reqs = (insert_requests(b, op["new"]) if op["new"] else []) + \
               [{"deleteContentRange": {"range": {"startIndex": a, "endIndex": b}}}]
        return a, reqs
    if kind == "insert_before":
        old = op["old"]
        assert text.count(old) == 1, f"{op['key']}: old occurs {text.count(old)}x"
        at = keep[text.index(old)][0]
        return at, insert_requests(at, op["new"])
    if kind == "replace_all":
        base = [c for c in p["chars"] if c[2] == "base" and c[1] not in ("\n",)]
        runs, cur = [], None                          # contiguous base ranges
        for c in base:
            if cur and c[0] == cur[1]:
                cur[1] = c[0] + u16(c[1])
            else:
                cur = [c[0], c[0] + u16(c[1])]; runs.append(cur)
        end = text_end(p)                            # after the last text character (before a page break)
        reqs = insert_requests(end, op["new"])
        for a, b in sorted(runs, reverse=True):
            reqs.append({"deleteContentRange": {"range": {"startIndex": a, "endIndex": b}}})
        return p["start"], reqs
    if kind == "append_paragraph":
        end = text_end(p)
        return end, insert_requests(end, "\n" + op["new"])
    raise ValueError(kind)


def plan(read_p, ops_p, only=None):
    doc = load(read_p)
    paras = paragraphs(doc)
    ops = [o for o in json.load(open(ops_p)) if only is None or o["key"] in only]
    built = [op_requests(paras, o) for o in ops]
    anchors = [a for a, _ in built]
    assert len(anchors) == len(set(anchors)), "two ops at one index"
    built.sort(key=lambda x: -x[0])
    return {"revisionId": doc["revisionId"], "requests": [r for _, rs in built for r in rs], "n_ops": len(ops)}


def verify(read_p, ops_p, only=None):
    paras = paragraphs(load(read_p))
    allp = "\n".join(proposed(p)[0] for p in paras)
    ok = True
    for o in json.load(open(ops_p)):
        if only is not None and o["key"] not in only:
            continue
        exp = o.get("expect") or parse(o["new"])[0].strip()
        good = exp in allp
        ok &= good
        print(f"{o['key']:28s}", "OK" if good else "MISSING")
    print("OK" if ok else "MISMATCH")


if __name__ == "__main__":
    cmd, rp, op = sys.argv[1:4]
    only = sys.argv[sys.argv.index("--only") + 1].split(",") if "--only" in sys.argv else None
    if cmd == "plan":
        json.dump(plan(rp, op, only), sys.stdout, ensure_ascii=False)
    else:
        verify(rp, op, only)
