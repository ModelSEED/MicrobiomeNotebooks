"""Independent verifier for the R232 tracked-change DOCX edits (docs lens).

For each (current, backup) pair:
  1. python-docx opens the current file.
  2. Audit every w:ins / w:del (and rPr-level markers, *PrChange) by author: count, date format, id uniqueness.
  3. Reject the NEW author's changes on the current document.xml and compare it with the backup's document.xml as a
     per-character stream (text, run formatting minus rsids, tracked-change wrapper chain of other authors) plus the
     paragraph / table / drawing structure.  Equality means: every textual or formatting edit is tracked, and rejecting
     all new changes restores the original.
  4. Accept the NEW author's changes and dump every paragraph whose visible text differs from the backup.
  5. Comments: new comments by the new author, anchored, and old comments untouched.
Usage: docx_verify.py <current.docx> <backup.docx> <new author> <out_prefix>
"""
import sys, zipfile, re, json, copy, difflib
from lxml import etree

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
A_NS = "http://schemas.openxmlformats.org/drawingml/2006/main"
WP = "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"
def q(t): return "{%s}%s" % (W, t)

cur_path, bak_path, NEW, outp = sys.argv[1:5]

import docx  # python-docx
d = docx.Document(cur_path)
print(f"python-docx load OK: {len(d.paragraphs)} body paragraphs, {len(d.tables)} tables, "
      f"{len(d.inline_shapes)} inline shapes")
db = docx.Document(bak_path)
print(f"backup: {len(db.paragraphs)} paragraphs, {len(db.tables)} tables, {len(db.inline_shapes)} inline shapes")

zc, zb = zipfile.ZipFile(cur_path), zipfile.ZipFile(bak_path)
parser = etree.XMLParser(remove_blank_text=False, huge_tree=True)
cur = etree.fromstring(zc.read("word/document.xml"), parser)
bak = etree.fromstring(zb.read("word/document.xml"), parser)

REV_TAGS = ["ins", "del", "moveFrom", "moveTo", "rPrChange", "pPrChange", "sectPrChange", "tblPrChange",
            "trPrChange", "tcPrChange", "numberingChange", "cellIns", "cellDel", "cellMerge", "tblGridChange"]

def audit(root, label):
    out = {}
    ids = []
    bad_date = []
    for t in REV_TAGS:
        for e in root.iter(q(t)):
            au = e.get(q("author"))
            parent_tag = etree.QName(e.getparent()).localname
            kind = t + ("@rPr" if parent_tag == "rPr" and t in ("ins", "del") else "")
            out.setdefault(au, {}).setdefault(kind, 0)
            out[au][kind] += 1
            ids.append(e.get(q("id")))
            dt = e.get(q("date"))
            if au == NEW and (dt is None or not re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ", dt)):
                bad_date.append((t, dt))
    dup = {i for i in ids if ids.count(i) > 1} if len(ids) < 20000 else None
    print(f"[{label}] revision elements by author: {json.dumps(out, sort_keys=True)}")
    print(f"[{label}] revision ids: {len(ids)} total, duplicates: {sorted(dup)[:20] if dup else 'none'}")
    if bad_date: print(f"[{label}] NEW-author bad dates: {bad_date[:5]}")
    return out

aud_c = audit(cur, "current")
aud_b = audit(bak, "backup")
if NEW in aud_b: print("!! backup already has new-author revisions")
for au, kinds in aud_b.items():
    print(f"   old author {au!r}: backup {kinds} | current {aud_c.get(au)}")

# ---------------------------------------------------------------- comments
def comments(z):
    if "word/comments.xml" not in z.namelist(): return {}
    c = etree.fromstring(z.read("word/comments.xml"))
    res = {}
    for e in c.iter(q("comment")):
        txt = "".join(e.itertext())
        res[e.get(q("id"))] = (e.get(q("author")), e.get(q("initials")), e.get(q("date")), txt,
                               etree.tostring(e))
    return res
cc, cb = comments(zc), comments(zb)
new_c = {i: v for i, v in cc.items() if i not in cb}
changed_old = [i for i in cb if i not in cc or cc[i][4] != cb[i][4]]
print(f"comments: backup {len(cb)}, current {len(cc)}, new {len(new_c)}, old changed/removed {changed_old}")
for i, v in new_c.items():
    anchors = (len(cur.xpath(f'//w:commentRangeStart[@w:id="{i}"]', namespaces={"w": W})),
               len(cur.xpath(f'//w:commentRangeEnd[@w:id="{i}"]', namespaces={"w": W})),
               len(cur.xpath(f'//w:commentReference[@w:id="{i}"]', namespaces={"w": W})))
    print(f"  new comment {i}: author={v[0]!r} initials={v[1]!r} date={v[2]} anchors(start,end,ref)={anchors} "
          f"text={v[3][:160]!r}{'...' if len(v[3])>160 else ''}")
    if v[0] != NEW: print("  !! new comment by another author")

# ---------------------------------------------------------------- reject / accept
def unwrap(e):
    p = e.getparent(); i = p.index(e)
    for k, ch in enumerate(list(e)):
        p.insert(i + k, ch)
    # tail handling irrelevant in OOXML (no mixed content)
    p.remove(e)

def reject_new(root):
    root = copy.deepcopy(root)
    n = {"ins": 0, "del": 0, "rpr": 0, "chg": 0, "cm": 0}
    for e in list(root.iter(q("ins"))):
        if e.get(q("author")) == NEW and etree.QName(e.getparent()).localname != "rPr":
            e.getparent().remove(e); n["ins"] += 1
    for e in list(root.iter(q("del"))):
        if e.get(q("author")) == NEW and etree.QName(e.getparent()).localname != "rPr":
            for t in e.iter(q("delText")): t.tag = q("t")
            for t in e.iter(q("delInstrText")): t.tag = q("instrText")
            unwrap(e); n["del"] += 1
    for e in list(root.iter(q("ins"), q("del"))):
        if e.get(q("author")) == NEW:  # rPr-level marker
            n["rpr"] += 1
            if e.tag == q("del"): e.getparent().remove(e)
            else: raise SystemExit("inserted paragraph mark by NEW: not handled")
    for t in ["rPrChange", "pPrChange", "sectPrChange", "tblPrChange", "trPrChange", "tcPrChange"]:
        for e in list(root.iter(q(t))):
            if e.get(q("author")) == NEW:
                n["chg"] += 1; raise SystemExit(f"{t} by NEW: not handled")
    newids = set(new_c)
    for tag in ["commentRangeStart", "commentRangeEnd"]:
        for e in list(root.iter(q(tag))):
            if e.get(q("id")) in newids: e.getparent().remove(e); n["cm"] += 1
    for e in list(root.iter(q("commentReference"))):
        if e.get(q("id")) in newids:
            r = e.getparent()
            r.remove(e)
            # remove the carrier run if it is now empty apart from rPr
            if all(etree.QName(ch).localname == "rPr" for ch in r):
                r.getparent().remove(r)
            n["cm"] += 1
    return root, n

def accept_new(root):
    root = copy.deepcopy(root)
    for e in list(root.iter(q("del"))):
        if e.get(q("author")) == NEW and etree.QName(e.getparent()).localname != "rPr":
            e.getparent().remove(e)
    for e in list(root.iter(q("ins"))):
        if e.get(q("author")) == NEW and etree.QName(e.getparent()).localname != "rPr":
            unwrap(e)
    return root

RSID = re.compile(r'\s+w:rsid\w*="[^"]*"')
def canon(e):
    if e is None: return ""
    s = etree.tostring(e, encoding="unicode")
    s = RSID.sub("", s)
    s = re.sub(r'\sxmlns(:\w+)?="[^"]*"', "", s)
    return s

WRAP = {"ins", "del", "hyperlink", "smartTag", "sdt", "sdtContent", "fldSimple", "customXml", "moveFrom", "moveTo"}

def stream(root):
    """Per-character / per-object stream of the body: (kind, value, rPr, wrapper chain)."""
    body = root.find(q("body"))
    out = []
    def chain_of(el, stop):
        ch = []
        x = el.getparent()
        while x is not None and x is not stop:
            ln = etree.QName(x).localname
            if ln in WRAP:
                ch.append((ln, x.get(q("author")), x.get(q("date"))))
            x = x.getparent()
        return tuple(reversed(ch))
    for el in body.iter():
        if not isinstance(el.tag, str): continue
        ln = etree.QName(el).localname
        if ln == "p":
            ppr = el.find(q("pPr"))
            out.append(("P", canon(ppr)))
        elif ln == "tbl": out.append(("TBL", canon(el.find(q("tblPr")))))
        elif ln == "tr": out.append(("TR", canon(el.find(q("trPr")))))
        elif ln == "tc": out.append(("TC", canon(el.find(q("tcPr")))))
        elif ln == "sectPr" and etree.QName(el.getparent()).localname in ("body", "pPr"):
            out.append(("SECT", canon(el)))
        elif ln == "bookmarkStart": out.append(("BM", el.get(q("name"))))
        elif ln in ("commentRangeStart", "commentRangeEnd"): out.append((ln, el.get(q("id"))))
        elif ln == "r":
            rpr = canon(el.find(q("rPr")))
            p = el
            while p is not None and etree.QName(p).localname != "p": p = p.getparent()
            ch = chain_of(el, p)
            for c in el:
                cl = etree.QName(c).localname
                if cl in ("t", "delText"):
                    for chh in (c.text or ""):
                        out.append(("C", chh, rpr, ch, cl))
                elif cl == "rPr": continue
                elif cl == "drawing":
                    emb = c.xpath(".//a:blip/@r:embed", namespaces={"a": A_NS, "r": R})
                    ext = c.xpath(".//wp:extent", namespaces={"wp": WP})
                    out.append(("DRAW", tuple(emb), tuple((x.get("cx"), x.get("cy")) for x in ext), rpr, ch))
                else:
                    out.append(("X", cl, canon(c), rpr, ch))
    return out

rej, n = reject_new(cur)
print(f"reject NEW: removed {n['ins']} w:ins, unwrapped {n['del']} w:del, {n['rpr']} rPr markers, "
      f"{n['cm']} comment markers")
s_rej, s_bak = stream(rej), stream(bak)
print(f"stream lengths: rejected {len(s_rej)}, backup {len(s_bak)}")
if s_rej == s_bak:
    print("REJECT CHECK: PASS (rejecting the new author's changes restores the backup exactly: text, run formatting, "
          "other authors' tracked changes, paragraph/table/section properties, drawings)")
else:
    sm = difflib.SequenceMatcher(None, s_rej, s_bak, autojunk=False)
    k = 0
    for op in sm.get_opcodes():
        if op[0] == "equal": continue
        k += 1
        if k <= 15:
            print("REJECT DIFF", op[0], "rej:", s_rej[op[1]:min(op[2], op[1] + 6)], "| bak:", s_bak[op[3]:min(op[4], op[3] + 6)])
    print(f"REJECT CHECK: FAIL ({k} differing blocks)")

# ---------------------------------------------------------------- accepted text per paragraph
def para_texts(root, accepted):
    res = []
    for p in root.iter(q("p")):
        parts = []
        for e in p.iter(q("t"), q("delText"), q("drawing")):
            # skip text in nested paragraphs (textboxes) handled by their own p
            pp = e.getparent()
            while pp is not None and pp.tag != q("p"): pp = pp.getparent()
            if pp is not p: continue
            if e.tag == q("delText"): continue  # visible text only (others' deletions hidden)
            if e.tag == q("drawing"):
                emb = e.xpath(".//a:blip/@r:embed", namespaces={"a": A_NS, "r": R})
                parts.append(f"[IMG {','.join(emb)}]")
            else:
                parts.append(e.text or "")
        res.append("".join(parts))
    return res
acc = accept_new(cur)
pa, pb = para_texts(acc, True), para_texts(bak, True)
print(f"paragraphs (all, incl. tables): accepted {len(pa)}, backup {len(pb)}")
changes = []
if len(pa) == len(pb):
    for i, (a, b) in enumerate(zip(pa, pb)):
        if a != b: changes.append({"p_all": i, "before": b, "after": a})
else:
    print("!! paragraph count differs")
json.dump(changes, open(outp + "_changed_paragraphs.json", "w"), indent=1, ensure_ascii=False)
print(f"changed paragraphs (accepted view vs backup): {len(changes)} -> {outp}_changed_paragraphs.json")
for c in changes:
    sm = difflib.SequenceMatcher(None, c["before"], c["after"], autojunk=False)
    for op in sm.get_opcodes():
        if op[0] != "equal":
            print(f"  p{c['p_all']}: {op[0]} {c['before'][max(0,op[1]-25):op[2]+5]!r} -> {c['after'][max(0,op[3]-25):op[4]+5]!r}")

# new-author revisions: list each with text
recs = []
for e in cur.iter(q("ins"), q("del")):
    if e.get(q("author")) != NEW or etree.QName(e.getparent()).localname == "rPr": continue
    txt = "".join(x.text or "" for x in e.iter(q("t"), q("delText")))
    img = e.xpath(".//a:blip/@r:embed", namespaces={"a": A_NS, "r": R})
    recs.append((etree.QName(e).localname, e.get(q("id")), e.get(q("date")), txt, img))
json.dump(recs, open(outp + "_new_revisions.json", "w"), indent=1, ensure_ascii=False)
print(f"new-author revision elements: {len(recs)} ({sum(1 for r in recs if r[0]=='ins')} ins, "
      f"{sum(1 for r in recs if r[0]=='del')} del); dates: {sorted({r[2] for r in recs})}")
etree.ElementTree(acc).write(outp + "_accepted_document.xml", encoding="UTF-8", xml_declaration=True, standalone=True)
