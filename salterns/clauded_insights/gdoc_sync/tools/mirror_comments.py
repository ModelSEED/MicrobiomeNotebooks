import sys
from docx import Document
from docx.oxml.ns import qn
doc = Document(sys.argv[1]); lo, hi = int(sys.argv[2]), int(sys.argv[3])
cm = {c.get(qn("w:id")): c for c in doc.comments._comments_elm}
for i, p in enumerate(doc.paragraphs[lo:hi], lo):
    # walk text, record comment ranges
    pos, txt, starts, ranges = 0, [], {}, []
    for e in p._p.iter():
        if e.tag == qn("w:commentRangeStart"): starts[e.get(qn("w:id"))] = pos
        elif e.tag == qn("w:commentRangeEnd"):
            cid = e.get(qn("w:id")); ranges.append((cid, starts.get(cid), pos))
        elif e.tag in (qn("w:t"), qn("w:delText")):
            anc = {a.tag for a in e.iterancestors()}
            t = e.text or ""
            txt.append(t); pos += len(t)
        elif e.tag == qn("w:tab"): txt.append("\t"); pos += 1
    s = "".join(txt)
    for cid, a, b in ranges:
        c = cm.get(cid)
        body = "".join(t.text or "" for t in c.iter(qn("w:t"))) if c is not None else "?"
        print(f"P{i} cid={cid} author={c.get(qn('w:author')) if c is not None else '?'} date={c.get(qn('w:date')) if c is not None else '?'} anchor={s[a:b]!r} :: {body[:90]!r}")
