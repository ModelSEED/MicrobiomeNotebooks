"""Independent check of Table S7 / guild CSV: tier rules on R220 vs R232 names; HEAD vs working tree."""
import csv, io, re, subprocess, json
REPO="/home/freiburger/Documents/EmilyKin"; PREP="/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep"
def lin(s):
    d={}
    for t in s.split(";"):
        d[t[:3]]=t[3:]
    return d
r232={r["user_genome"]:lin(r["classification"]) for r in csv.DictReader(open(f"{PREP}/results/gtdbtk.r232.summary.tsv"),delimiter="\t")}
h=subprocess.run(["git","-C",REPO,"show","HEAD:gtdbtk.bac120.summary.tsv"],capture_output=True,text=True).stdout.replace("\r\n","\n")
r220={r["user_genome"]:lin(r["classification"]) for r in csv.DictReader(io.StringIO(h),delimiter="\t")}
rules={"canPAO":("g__","Accumulibacter"),"canGAO_g":("g__","Competibacter|Contendobacter|Contendibacter"),"canGAO_f":("f__","Competibacteraceae"),
 "suppPAO":("g__","|".join(["Tetrasphaera","Tetrasphaera_A","Phosphoribacter","Lutibacillus","Dermatophostum","Dechloromonas","Microlunatus","Azonexus"])),
 "suppGAO":("g__","|".join(["Defluviicoccus","Propionivibrio","Micropruina"]))}
rules_head=dict(rules); rules_head["canGAO_g"]=("g__","Competibacter|Contendobacter")
diff=[]
for b in sorted(r232):
    for k,(rk,pat) in rules.items():
        a=bool(re.search(rules_head[k][1],r220[b].get(rk,""),re.I)); c=bool(re.search(pat,r232[b].get(rk,""),re.I))
        c_noD6=bool(re.search(rules_head[k][1],r232[b].get(rk,""),re.I))
        if a!=c or a!=c_noD6: diff.append((b,k,r220[b].get(rk),r232[b].get(rk),a,c,c_noD6))
print("tier flag differences R220(HEAD regex) vs R232(D6 regex / HEAD regex):",diff)
for f,sep,keep in [("Table_S7_MAG_functional_classification.csv",",",{"iterative_ID","Genus","Family","Phylum"}),("mag_guild_classification.csv",",",{"iterative_ID","Genus"})]:
    cur=list(csv.DictReader(open(f"{REPO}/{f}",newline="")))
    hd=list(csv.DictReader(io.StringIO(subprocess.run(["git","-C",REPO,"show",f"HEAD:{f}"],capture_output=True,text=True).stdout)))
    assert [r["MAG"] for r in cur]==[r["MAG"] for r in hd]
    assert list(cur[0])==list(hd[0])
    ch={}
    for a,b in zip(hd,cur):
        for k in a:
            if a[k]!=b[k]: ch.setdefault(k,0); ch[k]+=1
    print(f, "changed columns:", ch, "unexpected:", {k:v for k,v in ch.items() if k not in keep})
    mp=json.load(open(f"{REPO}/mag_iterativeID_r220_to_r232.json"))
    m={x["bin"]:x["new_id"] for x in mp} if isinstance(mp,list) else ({k:(v["new_id"] if isinstance(v,dict) else v) for k,v in mp.items()})
    badid=[r["MAG"] for r in cur if r["iterative_ID"]!=m[r["MAG"]]]
    print("  iterative_ID != R232 map:",len(badid),badid[:5], "| unique:", len({r['iterative_ID'] for r in cur})==len(cur))
