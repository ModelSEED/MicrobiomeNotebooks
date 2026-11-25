#!/usr/bin/env python
"""Swap figure images in the Google Doc as suggestions (tracked changes), and the Figure S6 caption.

  figure_swap.py plan <read_doc.json> <commit sha> [--only "Figure S11,Figure 7"] > requests.json
      Locates each figure by its caption (the paragraph whose proposed text starts with "Figure N." or "Figure N:"):
        - replace: the image-only paragraph(s) right above the caption hold the old image(s). The new image is
          suggested at the first old image's position and the old image(s) are suggested for deletion (Figure 8's two
          panels become one image);
        - insert: no image above the caption (Figures 9 and 10, new in the review). A new paragraph holding the image is
          suggested in front of the caption;
        - inline: the old image sits at the start of the caption paragraph (Figures 1 and S4). The new image is
          suggested in front of it and the old one is suggested for deletion.
      Images come from doc_images.py (images.json). Google fetches them from raw.githubusercontent.com at <commit sha>,
      so that commit must be pushed and contain them. Width 468 pt (6.5 in), or less so that the height stays within
      600 pt. Requests are ordered from the end of the Doc backward, so each one's indices are still valid when it runs.
      Prints {"requests": [...], "revisionId": ...} for update_doc with writeControl {writeMode: SUGGEST,
      requiredRevisionId}.
  figure_swap.py verify <read_doc.json> <commit sha> [--only ...]
      Checks that each figure now has the new image as a suggested insertion (by sourceUri) and its old image(s) as
      suggested deletions.
"""
import json
import re
import sys
from pathlib import Path

S = Path(__file__).resolve().parents[3]
IMAGES = json.load(open(S / "clauded_insights/gdoc_sync/doc_images/images.json"))
REPO = "https://raw.githubusercontent.com/ModelSEED/MicrobiomeNotebooks/{sha}/salterns/{path}"
MODE = {"Figure 1": "inline", "Figure S4": "inline", "Figure 9": "insert", "Figure 10": "insert"}
TEXT = {"weightedFontFamily": {"fontFamily": "Times New Roman", "weight": 400}, "fontSize": {"magnitude": 12, "unit": "PT"},
        "foregroundColor": {"color": {"rgbColor": {}}}, "italic": False, "bold": False, "baselineOffset": "NONE"}
TEXT_FIELDS = "weightedFontFamily,fontSize,foregroundColor,italic,bold,baselineOffset"

# Figure S6 caption (2026 rows): (old text, new text) in reading order; '' old = insert after the anchor
S6_EDITS = [
    (" and metagenomic data (bottom row)",
     ", metagenomic reads (middle row) and metagenome-assembled genomes (MAGs; bottom row) of the 24 sediment samples"),
    ("Metagenomic abundances were computed by the IMG/M pipeline by using lastal to map contigs to the IMG reference "
     "database. ", ""),
]
S6_APPEND_AFTER = "reference database (102)."
S6_APPEND = (" Metagenomic abundance is the percentage of the reads of eight single-copy ribosomal-protein genes whose best "
             "DIAMOND hit is a reference genome of the genus (GTDB R232); Halomonas includes Vreelandella, which GTDB "
             "separates from it, and Clostridium includes the GTDB genera with Clostridium suffixes (e.g. Clostridium_B). "
             "Only 8% of the reads assigned to Clostridium had a best hit at least 95% identical, so that assignment is "
             "uncertain. MAG abundance is the percentage of all reads that mapped to the genus’s species in the "
             "re-assembled MAG set; no species of Idiomarina, Clostridium or Marinilabilia was recovered. In the "
             "metagenomic and MAG rows, letters show Dunn’s tests (Holm-adjusted) after a significant Kruskal–Wallis "
             "test across the four sites (p < 0.05); sites sharing a letter do not differ.")
ITALIC = ["Halomonas", "Vreelandella", "Clostridium", "Idiomarina", "Marinilabilia"]


def u16(s):
    return len(s.encode("utf-16-le")) // 2


def load(p):
    d = json.load(open(p))
    return d.get("content", d)


def body_paragraphs(doc):
    return [e for e in doc["tabs"][0]["documentTab"]["body"]["content"] if "paragraph" in e]


def proposed_text(p):
    t = ""
    for x in p["paragraph"]["elements"]:
        if "textRun" in x and not x["textRun"].get("suggestedDeletionIds"):
            t += x["textRun"]["content"]
    return t


def images_in(p):
    return [x for x in p["paragraph"]["elements"] if "inlineObjectElement" in x]


def image_only(p):
    els = p["paragraph"]["elements"]
    txt = "".join(x["textRun"]["content"] for x in els if "textRun" in x)
    return bool(images_in(p)) and txt.strip() == ""


def caption_index(paras, fig):
    pat = re.compile(r"^\s*" + re.escape(fig) + r"[.:]")
    hits = [i for i, p in enumerate(paras) if pat.match(proposed_text(p))]
    assert len(hits) == 1, (fig, hits)
    return hits[0]


def size(fig):
    w, h = IMAGES[fig]["width_px"], IMAGES[fig]["height_px"]
    wpt = min(468.0, 600.0 * w / h)
    return {"width": {"magnitude": round(wpt, 1), "unit": "PT"}, "height": {"magnitude": round(wpt * h / w, 1), "unit": "PT"}}


def url(fig, sha):
    return REPO.format(sha=sha, path=IMAGES[fig]["path"])


def image_request(fig, sha, index):
    return {"insertInlineImage": {"uri": url(fig, sha), "location": {"index": index}, "objectSize": size(fig)}}


def figure_ops(paras, fig, sha):
    """[(anchor index, [requests])] for one figure."""
    ci = caption_index(paras, fig)
    cap = paras[ci]
    mode = MODE.get(fig, "replace")
    if mode == "insert":
        assert not image_only(paras[ci - 1]), (fig, "has an image above it; use replace")
        at = cap["startIndex"]
        return [(at, [{"insertText": {"text": "\n", "location": {"index": at}}}, image_request(fig, sha, at)])]
    if mode == "inline":
        olds = images_in(cap)
        assert len(olds) == 1 and olds[0]["startIndex"] == cap["startIndex"], (fig, "inline image not at the start")
    else:
        olds, j = [], ci - 1
        while image_only(paras[j]):
            olds = images_in(paras[j]) + olds
            j -= 1
        assert olds, (fig, "no image paragraph above the caption")
    for o in olds:
        assert not o["inlineObjectElement"].get("suggestedDeletionIds"), (fig, "old image already suggested for deletion")
    s, e = olds[0]["startIndex"], olds[-1]["endIndex"]
    # new image first, at the old one's position; then the old image(s), now one index later, are deleted
    return [(s, [image_request(fig, sha, s), {"deleteContentRange": {"range": {"startIndex": s + 1, "endIndex": e + 1}}}])]


def s6_caption_ops(paras):
    cap = paras[caption_index(paras, "Figure S6")]
    runs = [x for x in cap["paragraph"]["elements"] if "textRun" in x]
    assert not any(x["textRun"].get("suggestedInsertionIds") or x["textRun"].get("suggestedDeletionIds") for x in runs), \
        "the S6 caption already holds suggestions"
    # character offset -> Doc index over the caption's text runs
    idx, text = [], ""
    for x in cap["paragraph"]["elements"]:
        if "textRun" in x:
            i = x["startIndex"]
            for ch in x["textRun"]["content"]:
                idx.append(i); text += ch; i += u16(ch)
        else:
            idx.append(x["startIndex"]); text += "\x00"
    ops = []
    for old, new in S6_EDITS:
        assert text.count(old) == 1, old
        k = text.index(old)
        a, b = idx[k], idx[k + len(old) - 1] + u16(old[-1])
        reqs = []
        if new:
            reqs += [{"insertText": {"text": new, "location": {"index": b}}},
                     {"updateTextStyle": {"range": {"startIndex": b, "endIndex": b + u16(new)}, "textStyle": TEXT,
                                          "fields": TEXT_FIELDS}}]
        reqs.append({"deleteContentRange": {"range": {"startIndex": a, "endIndex": b}}})
        ops.append((a, reqs))
    assert text.count(S6_APPEND_AFTER) == 1
    k = text.index(S6_APPEND_AFTER) + len(S6_APPEND_AFTER)
    at = idx[k]
    reqs = [{"insertText": {"text": S6_APPEND, "location": {"index": at}}},
            {"updateTextStyle": {"range": {"startIndex": at, "endIndex": at + u16(S6_APPEND)}, "textStyle": TEXT,
                                 "fields": TEXT_FIELDS}}]
    for w in ITALIC:
        for m in re.finditer(r"\b" + w + r"\b", S6_APPEND):
            st = at + u16(S6_APPEND[:m.start()])
            reqs.append({"updateTextStyle": {"range": {"startIndex": st, "endIndex": st + u16(w)},
                                             "textStyle": {"italic": True}, "fields": "italic"}})
    ops.append((at, reqs))
    return ops


def plan(read_p, sha, only=None):
    doc = load(read_p)
    paras = body_paragraphs(doc)
    figs = [f for f in IMAGES if only is None or f in only]
    ops = []
    for f in figs:
        ops += figure_ops(paras, f, sha)
    if only is None or "Figure S6 caption" in only:
        ops += s6_caption_ops(paras)
    starts = sorted(a for a, _ in ops)
    assert len(starts) == len(set(starts)), "two operations at one index"
    ops.sort(key=lambda o: -o[0])
    return {"revisionId": doc["revisionId"], "requests": [r for _, rs in ops for r in rs]}


def verify(read_p, sha, only=None):
    doc = load(read_p)
    paras = body_paragraphs(doc)
    io = doc["tabs"][0]["documentTab"]["inlineObjects"]
    ok = True
    for f in [f for f in IMAGES if only is None or f in only]:
        ci = caption_index(paras, f)
        j = ci - 1                                     # the image-only paragraphs right above the caption, and the caption
        while j >= 0 and image_only(paras[j]):
            j -= 1
        span = paras[j + 1:ci + 1]
        objs = [x for p in span for x in images_in(p)]
        new = [x for x in objs if io[x["inlineObjectElement"]["inlineObjectId"]]["inlineObjectProperties"]
               ["embeddedObject"].get("imageProperties", {}).get("sourceUri") == url(f, sha)]
        ins = bool(new) and all(x["inlineObjectElement"].get("suggestedInsertionIds") for x in new)
        dels = [x for x in objs if x["inlineObjectElement"].get("suggestedDeletionIds")]
        want_del = MODE.get(f, "replace") != "insert"
        good = ins and len(new) == 1 and (bool(dels) == want_del)
        ok &= good
        print(f"{f:11s} new image suggested: {ins} ({len(new)}), old suggested for deletion: {len(dels)}",
              "OK" if good else "CHECK")
    print("OK" if ok else "MISMATCH")


if __name__ == "__main__":
    cmd, rp, sha = sys.argv[1:4]
    only = None
    if "--only" in sys.argv:
        only = [x.strip() for x in sys.argv[sys.argv.index("--only") + 1].split(",")]
    if cmd == "plan":
        json.dump(plan(rp, sha, only), sys.stdout, ensure_ascii=False)
    else:
        verify(rp, sha, only)
