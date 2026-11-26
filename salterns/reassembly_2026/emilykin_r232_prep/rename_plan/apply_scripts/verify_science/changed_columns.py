"""For every modified tracked .tsv/.csv: which columns changed vs HEAD (row-aligned), row count equality."""
import csv, io, subprocess, sys
from collections import Counter
csv.field_size_limit(10**9)
R="/home/freiburger/Documents/EmilyKin"
files=[l[3:] for l in subprocess.run(["git","-C",R,"status","--porcelain"],capture_output=True,text=True).stdout.splitlines() if l.startswith(" M ") and l[3:].strip('"').endswith((".tsv",".csv"))]
for f in files:
    f=f.strip('"')
    sep="\t" if f.endswith(".tsv") else ","
    h=subprocess.run(["git","-C",R,"show",f"HEAD:{f}"],capture_output=True).stdout.decode("utf-8","replace").replace("\r\n","\n")
    c=open(f"{R}/{f}",encoding="utf-8",errors="replace",newline="").read().replace("\r\n","\n")
    H=list(csv.reader(io.StringIO(h),delimiter=sep)); C=list(csv.reader(io.StringIO(c),delimiter=sep))
    if not H: print("EMPTY HEAD",f); continue
    hh,ch=H[0],C[0]
    note=""
    if hh!=ch: note=f" HEADER added={[x for x in ch if x not in hh]} removed={[x for x in hh if x not in ch]}"
    cnt=Counter()
    if len(H)!=len(C): print(f"ROWS {f}: {len(H)} -> {len(C)}{note}"); continue
    idx={x:i for i,x in enumerate(ch)}
    for a,b in zip(H[1:],C[1:]):
        for i,x in enumerate(hh):
            if x in idx and i<len(a) and idx[x]<len(b) and a[i]!=b[idx[x]]: cnt[x]+=1
    print(f"{f}: rows={len(C)-1} changed={dict(cnt)}{note}")
