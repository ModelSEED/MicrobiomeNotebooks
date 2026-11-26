"""D1 sweep: every modified/new table with a bin column + MAG-ID column equals the R232 map by bin (unique per file
where the file has one row per bin); plus CRLF line counts HEAD vs working tree for every modified text file."""
import csv, io, json, subprocess, os
from collections import defaultdict
csv.field_size_limit(10**9)
R="/home/freiburger/Documents/EmilyKin"
ID={b:v["new_id"] for b,v in json.load(open(f"{R}/mag_iterativeID_r220_to_r232.json")).items()}
st=subprocess.run(["git","-C",R,"status","--porcelain"],capture_output=True,text=True).stdout.splitlines()
mod=[l[3:].strip('"') for l in st if l.startswith(" M ")]
new=[l[3:].strip('"') for l in st if l.startswith("?? ")]
BINC=("MAG","mag","bin","tip","MAG_bin","source_MAG_bin","target_MAG_bin","guild_node_bin","source_bin","target_bin")
IDC=("iterativeID","mag_iterativeID","iterative_ID","MAG_iteration_ID","new_mag_iterativeID","MAG_iterativeID","new_id",
     "source_MAG_iterativeID","target_MAG_iterativeID","guild_node_MAG_iterativeID")
PAIRS={"source_MAG_iterativeID":["source_MAG_bin","source_bin"],"target_MAG_iterativeID":["target_MAG_bin","target_bin"],
       "guild_node_MAG_iterativeID":["guild_node_bin"],"MAG_iterativeID":["MAG_bin","MAG"]}
for f in sorted(mod+[x for x in new if x.endswith((".tsv",".csv"))]):
    if not f.endswith((".tsv",".csv")) or not os.path.isfile(f"{R}/{f}"): continue
    sep="\t" if f.endswith(".tsv") else ","
    rows=list(csv.DictReader(open(f"{R}/{f}",newline="",encoding="utf-8",errors="replace"),delimiter=sep))
    if not rows: continue
    cols=list(rows[0])
    for ic in [c for c in cols if c in IDC]:
        bcs=PAIRS.get(ic,[c for c in BINC if c in cols])
        bc=[c for c in bcs if c in cols]
        if not bc: continue
        bc=bc[0]
        n=bad=blank=0; ex=[]; per_bin=defaultdict(set)
        for r in rows:
            b=r[bc]
            if b not in ID: continue
            n+=1; v=r[ic]
            if v=="": blank+=1; continue
            per_bin[v].add(b)
            if v!=ID[b]: bad+=1; ex.append((b,v,ID[b]))
        dup=[k for k,s in per_bin.items() if len(s)>1]
        if n: print(f"{'OK ' if not bad and not dup else 'BAD'} {f} [{bc}->{ic}] rows={n} wrong={bad} blank={blank} ids_shared_by_bins={len(dup)} {ex[:3]}")
# JSON: mag_abundance_summary.json
d=json.load(open(f"{R}/mag_abundance_summary.json"))
bad=[b for b in d if d[b].get("mag_iterativeID")!=ID[b]]
print("mag_abundance_summary.json mag_iterativeID wrong:",bad)
# EOL
print("\nEOL (CRLF line counts HEAD -> now) for modified text files with a change:")
for f in mod:
    try: c=open(f"{R}/{f}","rb").read()
    except IsADirectoryError: continue
    if b"\0" in c[:8000]: continue
    h=subprocess.run(["git","-C",R,"show",f"HEAD:{f}"],capture_output=True).stdout
    hc,hl=h.count(b"\r\n"),h.count(b"\n"); cc,cl=c.count(b"\r\n"),c.count(b"\n")
    hk="CRLF" if hc==hl and hl else ("LF" if hc==0 else "MIXED"); ck="CRLF" if cc==cl and cl else ("LF" if cc==0 else "MIXED")
    if hk!=ck: print("  EOL CHANGED",f,hk,"->",ck,(hc,hl),(cc,cl))
print("EOL check done over",len(mod),"modified files")
