"""Quantify the skani -n 3 truncation: re-bucket the 66 low_AF_unresolved MAGs of Table S6 using an unrestricted
skani search (same DB, --min-af 0, no -n) with the Table S6 rule. Input: apply_scratch/verify_science/skani_lowaf_non/all_hits.tsv."""
import csv, os, re, gzip
from collections import Counter
V="/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan/apply_scratch/verify_science"
REPO="/home/freiburger/Documents/EmilyKin"; P="/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep"
tax={}
for l in gzip.open(f"{P}/gtdb_downloads/r232/bac120_taxonomy_r232.tsv.gz","rt"):
    a,t=l.rstrip("\n").split("\t"); tax[a[3:]]=t
T={r["MAG"]:r for r in csv.DictReader(open(f"{REPO}/Table_S6_skani_novelty.tsv",newline=""),delimiter="\t")}
hits={}
for r in csv.DictReader(open(f"{V}/skani_lowaf_non/all_hits.tsv"),delimiter="\t"):
    m=os.path.basename(r["Query_file"])[:-3]; acc=re.search(r"(GC[AF]_\d+\.\d+)",r["Ref_file"]).group(1)
    hits.setdefault(m,[]).append((float(r["ANI"]),float(r["Align_fraction_query"]),float(r["Align_fraction_ref"]),acc))
low=sorted(m for m in T if T[m]["ANI_class"]=="low_AF_unresolved")
moves=Counter(); out=[]
for m in low:
    el=[h for h in hits.get(m,[]) if h[1]>=15]
    if el:
        b=max(el,key=lambda h:h[0]); c=">=95" if b[0]>=95 else ("85-95" if b[0]>=85 else "<85")
        out.append((m,T[m]["skani_ANI"],T[m]["skani_AF_query"],c,f"{b[0]:.2f}",f"{b[1]:.2f}",b[3],tax.get(b[3],"").split(";")[-1]))
    else: c="low_AF_unresolved"
    moves[c]+=1
print("low_AF MAGs:",len(low),"| re-bucketed without -n 3:",dict(moves))
cur=Counter(r["ANI_class"] for r in T.values())
alt=Counter(cur); alt["low_AF_unresolved"]-=len(out)
for r in out: alt[r[3]]+=1
print("Table S6 as built (-n 3):",dict(cur),"| without the -n 3 truncation:",dict(alt))
print("MAG\tS6_ANI(-n3)\tS6_AFq(-n3)\tclass_unrestricted\tANI\tAF_query\tref\tref_species")
for r in out: print("\t".join(r))
