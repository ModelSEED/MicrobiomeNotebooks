"""Map Word comments in a mirror to Google comment threads (Docs API) by author + createTime (to the second)."""
import json, re, sys
from docx import Document
from docx.oxml.ns import qn

def api_posts(api_path):
    c = json.load(open(api_path))["content"]
    posts = {}
    for t in c["comments"]:
        for kind, p in [("head", t["headPost"])] + [("reply", r) for r in t.get("replies", [])]:
            if kind == "reply" and not p.get("content"):   # resolve/reopen actions
                continue
            key = (p["author"].get("displayName", "?"), p["createTime"][:19] + "Z")
            posts.setdefault(key, []).append((t["commentId"], kind, t["status"], t.get("plainTextQuote", "")))
    return c["comments"], posts

def export_anchors(export_path):
    fc = json.load(open(export_path))["fileContent"]
    # anchored text per comment id (all segments), with markdown links reduced to their text
    out = {}
    for cid in set(re.findall(r"<comment_start id=([^>]+)>", fc)):
        segs = re.findall(rf"<comment_start id={re.escape(cid)}>(.*?)<comment_end id={re.escape(cid)}>", fc, flags=re.S)
        out[cid] = segs
    return out

def mdnorm(s):
    s = re.sub(r"<comment_(start|end) id=[^>]+>", "", s)
    s = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", s)
    s = s.replace("\\", "").replace("*", "").replace("#", "").replace("&#181;", "µ")
    return re.sub(r"\s+", "", s)

def mirror_comments(docx_path):
    doc = Document(docx_path)
    meta = {c.get(qn("w:id")): (c.get(qn("w:author")), c.get(qn("w:date")), "".join(t.text or "" for t in c.iter(qn("w:t"))))
            for c in doc.comments._comments_elm}
    body_ps = list(doc.element.body.iter(qn("w:p")))
    top = {id(p._p): i for i, p in enumerate(doc.paragraphs)}
    res = {}
    for pi, p in enumerate(body_ps):
        open_, base, acc = {}, [], []
        for e in p.iter():
            if e.tag == qn("w:commentRangeStart"): open_[e.get(qn("w:id"))] = (len("".join(base)), len("".join(acc)))
            elif e.tag == qn("w:commentRangeEnd"):
                cid = e.get(qn("w:id")); b0, a0 = open_.pop(cid, (None, None))
                if b0 is None: continue
                res.setdefault(cid, []).append({"bodyp": pi, "base": "".join(base)[b0:], "acc": "".join(acc)[a0:], "P": top.get(id(p))})
            elif e.tag in (qn("w:t"), qn("w:delText"), qn("w:tab")):
                anc = {a.tag for a in e.iterancestors()}
                t = "\t" if e.tag == qn("w:tab") else (e.text or "")
                if qn("w:ins") not in anc: base.append(t)
                if qn("w:del") not in anc: acc.append(t)
        # fix: slices above taken at end-marker time
    return doc, meta, res

if __name__ == "__main__":
    mirror, api, export = sys.argv[1:4]
    comments, posts = api_posts(api)
    anchors = export_anchors(export)
    doc, meta, res = mirror_comments(mirror)
    for wid, (author, date, text) in sorted(meta.items(), key=lambda x: int(x[0])):
        hits = posts.get((author, date), [])
        r = res.get(wid, [])
        where = ",".join(f"b{x['bodyp']}" for x in r) or "UNANCHORED"
        h = hits[0] if hits else ("?", "?", "?", "")
        print(f"w{wid}\t{where}\t{h[0]}\t{h[1]}\t{h[2]}\tin_export={h[0] in anchors}\t{author}\t{date}\t{text[:50]!r}")
