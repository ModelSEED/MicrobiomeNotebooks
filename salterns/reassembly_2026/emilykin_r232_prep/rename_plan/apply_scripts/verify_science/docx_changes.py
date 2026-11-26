"""Print old/new text of every paragraph touched by a given tracked-change author, and its comments (read-only)."""
import sys, zipfile, re
from lxml import etree
W="http://schemas.openxmlformats.org/wordprocessingml/2006/main"; ns={"w":W}
path, author = sys.argv[1], sys.argv[2]
z=zipfile.ZipFile(path)
doc=etree.fromstring(z.read("word/document.xml"))
com={}
if "word/comments.xml" in z.namelist():
    c=etree.fromstring(z.read("word/comments.xml"))
    for e in c.findall("w:comment",ns):
        if e.get(f"{{{W}}}author")==author:
            com[e.get(f"{{{W}}}id")]="".join(e.itertext())
def text(p, mode):
    out=[]
    for el in p.iter():
        tag=etree.QName(el).localname
        if tag not in ("t","delText"): continue
        anc=[etree.QName(a).localname+"|"+(a.get(f"{{{W}}}author") or "") for a in el.iterancestors() if etree.QName(a).localname in("ins","del")]
        mine_ins=any(a==f"ins|{author}" for a in anc); mine_del=any(a==f"del|{author}" for a in anc)
        if mode=="old" and mine_ins: continue
        if mode=="new" and mine_del: continue
        out.append(el.text or "")
    return "".join(out)
body=doc.find("w:body",ns)
for i,p in enumerate(body.iter(f"{{{W}}}p")):
    touched=[a for a in p.iter() if etree.QName(a).localname in ("ins","del") and a.get(f"{{{W}}}author")==author]
    cids=[e.get(f"{{{W}}}id") for e in p.iter(f"{{{W}}}commentRangeStart")]
    if not touched and not any(c in com for c in cids): continue
    o,n=text(p,"old"),text(p,"new")
    print(f"\n### para {i}  ({len(touched)} change elements)")
    if o!=n:
        print("OLD:",o[:3000]); print("NEW:",n[:3000])
    else:
        print("TEXT (unchanged / image-only):",n[:300])
    for c in cids:
        if c in com: print(f"  COMMENT {c}:",com[c][:2500])
