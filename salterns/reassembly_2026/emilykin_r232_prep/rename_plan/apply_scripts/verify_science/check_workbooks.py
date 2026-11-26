"""Workbook taxonomy / ID cells vs R232 by bin (read-only, openpyxl)."""
import csv, json, re, sys
import openpyxl
REPO="/home/freiburger/Documents/EmilyKin"; PREP="/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep"
S={r["user_genome"]:r for r in csv.DictReader(open(f"{PREP}/results/gtdbtk.r232.summary.tsv"),delimiter="\t")}
L={b:{t[:3]:t[3:] for t in r["classification"].split(";")} for b,r in S.items()}
ID={b:v["new_id"] for b,v in json.load(open(f"{REPO}/mag_iterativeID_r220_to_r232.json")).items()}
RK={"phylum":"p__","class":"c__","order":"o__","family":"f__","genus":"g__","species":"s__"}
def norm(h): return re.sub(r"\s+"," ",str(h or "")).strip().lower()
for wbn in ["MAG_target_gene_filled.xlsx","MAG_target_gene_filled_UPDATED (1).xlsx","cladeexcel.xlsx","abundance and coverage (1).xlsx"]:
    wb=openpyxl.load_workbook(f"{REPO}/{wbn}",read_only=True,data_only=False)
    for ws in wb.worksheets:
        rows=list(ws.iter_rows(values_only=True))
        # find header row: a row with a cell that looks like a bin column header and a following row with a bin value
        hdr=None
        for i,r in enumerate(rows[:8]):
            if any(norm(c) in ("mag bin name","mags","mag","bin","mag bin","mag_bin","mags bin name") for c in r):
                hdr=i;break
        if hdr is None: print(wbn,"|",ws.title,": no bin header"); continue
        H=[norm(c) for c in rows[hdr]]
        bcol=[j for j,h in enumerate(H) if h in ("mag bin name","mags","mag","bin","mag bin","mag_bin","mags bin name")][0]
        checks={}
        for j,h in enumerate(H):
            if h in RK: checks[j]=("rank",RK[h])
            elif h in ("taxonomy","gtdb taxonomy","classification","gtdb classification","mag_gtdb_taxonomy"): checks[j]=("lin",None)
            elif "iteration id" in h or h in ("iterative id","iterative_id","iterativeid","mag iterative id","mag_iterativeid"): checks[j]=("id",None)
        res={}
        nb=0
        for r in rows[hdr+1:]:
            b=str(r[bcol]).strip() if r[bcol] is not None else ""
            if b not in S: continue
            nb+=1
            for j,(kind,rk) in checks.items():
                v=r[j] if j<len(r) else None; v="" if v is None else str(v).strip()
                if kind=="rank":
                    exp=L[b].get(rk,"")
                    ok = v==exp or v==rk+exp or (exp=="" and v in("","N/A","NA","-"))
                elif kind=="lin": ok = v==S[b]["classification"]
                else: ok = v==ID[b]
                res.setdefault(H[j],[0,0,[]]); res[H[j]][0]+=1
                if not ok: res[H[j]][1]+=1; res[H[j]][2].append((b,v[:60],(L[b].get(rk,"") if kind=="rank" else (ID[b] if kind=="id" else "lin"))))
        print(wbn,"|",ws.title,"| bins",nb,"|",{k:(v[0],v[1],v[2][:3]) for k,v in res.items()})
