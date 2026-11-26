"""Boundary-aware scan of text files for: 51 old MAG IDs (with bin), R220-only GTDB names, release strings."""
import csv, re, sys, json
P="/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan"
idm=list(csv.DictReader(open(f"{P}/docs/id_change_map.tsv"),delimiter="\t"))
tx=list(csv.DictReader(open(f"{P}/data/taxon_name_changes_r220_r232.tsv"),delimiter="\t"))
r220only=sorted({r["r220_name"] for r in tx if r["r220_name_still_used_in_r232_mag_taxonomy"]=="False" and r["r220_name"]})
pats=[]
for r in idm:
    pats.append(("oldID", r["old_id"], r["MAG"]+"->"+r["new_id"]+(" [ASV-shared]" if r["old_is_asv_id"]=="True" else "")))
for n in r220only:
    pats.append(("r220name", n, ""))
for s in ["r220","R220","v2.7.2","release220","155 ASV","62 GTDB","74 GTDB","Patescibacteria"]:
    pats.append(("misc", s, ""))
def rx(s):
    return re.compile(r"(?<![A-Za-z0-9_])"+re.escape(s)+r"(?![0-9A-Za-z_])")
cp=[(k,s,n,rx(s)) for k,s,n in pats]
for f in sys.argv[1:]:
    try: lines=open(f,encoding="utf-8").read().split("\n")
    except Exception as e: print("ERR",f,e); continue
    hits=0
    for i,l in enumerate(lines,1):
        found=[(k,s,n) for k,s,n,r in cp if r.search(l)]
        if found:
            hits+=1
            print(f"{f}:{i}: "+"; ".join(f"{k}:{s}{(' '+n) if n else ''}" for k,s,n in found))
            print("    "+l.strip()[:400])
    print(f"== {f}: {hits} lines")
