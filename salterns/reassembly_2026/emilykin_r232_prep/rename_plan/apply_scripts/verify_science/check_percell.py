import csv, io, re, subprocess
from collections import defaultdict
R="/home/freiburger/Documents/EmilyKin"
cur=list(csv.DictReader(open(f"{R}/mag_percell_abundance.csv",newline="")))
hd=list(csv.DictReader(io.StringIO(subprocess.run(["git","-C",R,"show","HEAD:mag_percell_abundance.csv"],capture_output=True,text=True).stdout)))
assert [r["MAG"] for r in cur]==[r["MAG"] for r in hd]
ch=defaultdict(int)
for a,b in zip(hd,cur):
    for k in a:
        if a[k]!=b[k]: ch[k]+=1
print("percell changed columns:",dict(ch))
pc=[c for c in cur[0] if c.startswith("per_cell_abun_CAN_")]
print(pc)
def key(g): return "Unclassified" if not g.strip() else re.sub(r"_[A-Z]+$","",g.strip())
for label,rows in [("HEAD/R220",hd),("now/R232",cur)]:
    g=defaultdict(lambda:[0.0]*5)
    for r in rows:
        for i,c in enumerate(pc): g[key(r["Genus"])][i]+=float(r[c])
    order=sorted(g,key=lambda k:-sum(g[k])/5)
    print(label,"top20:",[ (k,round(sum(g[k])/5,2)) for k in order[:20]])
    for k in ["Competibacter","Accumulibacter","Thauera","Azonexus","Unclassified","Lysobacter","Proximibacter","Simplicispira","JJ008"]:
        print("   ",k,[round(x,2) for x in g.get(k,[0]*5)],"mean",round(sum(g.get(k,[0]*5))/5,2), "MAGs:",[r["MAG"] for r in rows if key(r["Genus"])==k][:8])
