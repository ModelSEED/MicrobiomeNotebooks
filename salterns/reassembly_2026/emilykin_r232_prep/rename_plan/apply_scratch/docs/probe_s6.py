import openpyxl, csv, io, subprocess, sys
R="/home/freiburger/Documents/EmilyKin"
B="/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan/apply_scratch/data/backups/MAG_target_gene_filled.HEAD.xlsx"
sub=openpyxl.load_workbook(f"{R}/MAG_target_gene_filled_UPDATED (1).xlsx")["Table S6_ANI"]
bh=openpyxl.load_workbook(B)["Table S6_ANI"]
bn=openpyxl.load_workbook(f"{R}/MAG_target_gene_filled.xlsx")["Table S6_ANI"]
def rows(ws,nc):
    hdr=[ws.cell(1,c).value for c in range(1,nc+1)]
    d={}
    for r in range(2,ws.max_row+1):
        d[ws.cell(r,1).value]=[ws.cell(r,c).value for c in range(1,nc+1)]
    return hdr,d
hs,ds=rows(sub,27); hb,db=rows(bh,20); hn,dn=rows(bn,20)
print("hdr A-T equal (except A):", hs[1:20]==hb[1:20], hs[0], hb[0])
from collections import Counter
diff=Counter()
for b in ds:
    for i in range(1,20):
        if ds[b][i]!=db[b][i]: diff[hs[i]]+=1
print("submission vs builder HEAD, A-T diffs:", dict(diff))
t=list(csv.DictReader(io.StringIO(subprocess.check_output(["git","-C",R,"show","HEAD:Table_S6_skani_novelty.tsv"]).decode()),delimiter="\t"))
th={r["MAG"]:r for r in t}
cols=["skani_best_ref","skani_best_taxonomy","skani_ANI","skani_AF_query","skani_AF_ref","ANI_class","Novelty_interpretation"]
diff2=Counter(); ex=[]
for b in ds:
    for j,c in enumerate(cols):
        v=ds[b][20+j]; w=th[b][c]
        vv="" if v is None else (str(v) if not isinstance(v,float) else repr(v))
        if vv!=w:
            try:
                if float(vv)==float(w): continue
            except: pass
            diff2[c]+=1; ex.append((b,c,v,w))
print("submission U-AA vs HEAD Table_S6 tsv:", dict(diff2), ex[:5])
print("types U-AA:", Counter(type(ds[b][22]).__name__ for b in ds))
