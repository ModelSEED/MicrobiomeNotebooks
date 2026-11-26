import sys, re, zipfile
from lxml import etree
W="{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
def final_text(p):
    out=[]
    for el in p.iter(W+"t", W+"delText", W+"tab", W+"br"):
        # skip deleted
        anc=[a.tag for a in el.iterancestors()]
        if W+"del" in anc or el.tag==W+"delText": continue
        if el.tag==W+"t": out.append(el.text or "")
        elif el.tag==W+"tab": out.append("\t")
    return "".join(out)
path=sys.argv[1]; pats=sys.argv[2:]
root=etree.fromstring(zipfile.ZipFile(path).read("word/document.xml"))
body=root.find(W+"body")
for i,p in enumerate(body.iter(W+"p"),1):
    t=final_text(p)
    if any(re.search(x,t) for x in pats):
        ins=len(p.findall(".//"+W+"ins")); dl=len(p.findall(".//"+W+"del"))
        intbl=any(a.tag==W+"tc" for a in p.iterancestors())
        print(f"## para {i} ins={ins} del={dl} tbl={intbl} runs={len(p.findall('.//'+W+'r'))}\n{t}\n")
