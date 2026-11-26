"""Scan manuscript .md and accepted DOCX text for R220-only GTDB names (no longer used by any MAG in R232)
and for the 51 old MAG IDs, printing each hit with context, for manual classification."""
import csv, re, sys, zipfile, json
from lxml import etree
REPO="/home/freiburger/Documents/EmilyKin"; PLAN="/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan"
names=set()
for r in csv.DictReader(open(f"{PLAN}/data/taxon_name_changes_r220_r232.tsv"),delimiter="\t"):
    if r["r220_name_still_used_in_r232_mag_taxonomy"]=="False" and r["r220_name"]:
        names.add(r["r220_name"])
mp=json.load(open(f"{REPO}/mag_iterativeID_r220_to_r232.json"))
oldids={v["old_id"] for v in mp.values() if v["old_id"]!=v["new_id"]}
W="http://schemas.openxmlformats.org/wordprocessingml/2006/main"
def accepted_paras(path):
    z=zipfile.ZipFile(path); doc=etree.fromstring(z.read("word/document.xml"))
    out=[]
    for p in doc.iter(f"{{{W}}}p"):
        t=[]
        for el in p.iter(f"{{{W}}}t"):
            if any(etree.QName(a).localname=="del" for a in el.iterancestors()): continue
            t.append(el.text or "")
        out.append("".join(t))
    return out
def scan(label, paras):
    hits=0
    for i,p in enumerate(paras):
        for n in sorted(names,key=len,reverse=True):
            for m in re.finditer(r"(?<![A-Za-z0-9_])"+re.escape(n)+r"(?![A-Za-z0-9])",p):
                hits+=1; s=max(0,m.start()-90); print(f"[{label} {i}] NAME {n!r}: ...{p[s:m.end()+60]}...")
        for o in oldids:
            for m in re.finditer(r"(?<![A-Za-z0-9_.])"+re.escape(o)+r"(?![0-9]|_m)",p):
                hits+=1; s=max(0,m.start()-90); print(f"[{label} {i}] OLDID {o!r}: ...{p[s:m.end()+60]}...")
    print(f"== {label}: {hits} hits")
scan("md", open(f"{REPO}/MANUSCRIPT_methods_results_discussion.md").read().split("\n"))
for f in ["papert3_070626 (2)_correction.docx","SI-paper3 (3)_correction.docx"]:
    scan(f[:12], accepted_paras(f"{REPO}/{f}"))
