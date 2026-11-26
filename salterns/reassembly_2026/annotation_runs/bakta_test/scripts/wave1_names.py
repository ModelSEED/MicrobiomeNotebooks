import json, sys, collections
import pandas as pd
SP=sys.argv[1]; MARK=sys.argv[2]; OUT=sys.argv[3]
pd.set_option("display.width",250); pd.set_option("display.max_colwidth",75)
mk=pd.read_csv(MARK,sep="\t",index_col=0)
kol=pd.read_csv(f"{SP}/grdscan/ko_list",sep="\t",index_col=0); thr=pd.to_numeric(kol.threshold,errors="coerce")
rows=[]
for l in open(f"{SP}/interim/markers.domtblout"):
    if l.startswith("#"): continue
    f=l.split(); rows.append((f[0],f[3],float(f[6]),float(f[7]),float(f[13])))
h=pd.DataFrame(rows,columns=["gene","ko","E","full","dom"])
h["score"]=h.full.where(h.ko.map(kol.score_type)=="full",h.dom); h["thr"]=h.ko.map(thr); h["sig"]=h.score>h.thr
best=h.sort_values("E").drop_duplicates("gene").set_index("gene")
sig=h[h.sig].sort_values("E").drop_duplicates("gene").set_index("gene")
summ=[]; det=[]
LABS=sys.argv[4].split(",")
for lab in LABS:
    d=json.load(open(f"{SP}/wave1/bakta/{lab}/{lab}.json"))
    cds=[f for f in d["features"] if f["type"]=="cds"]
    src=collections.Counter("IPS/UPS" if (f.get("ups") or f.get("ips")) else "PSC" if (f.get("psc") or {}).get("valid") else "PSCC" if f.get("pscc") else "PSC<50%" if f.get("psc") else "none" for f in cds)
    n=len(cds)
    kx=[f for f in cds if any(x.startswith("KEGG:") for x in f.get("db_xrefs",[]))]
    summ.append(dict(label=lab,cds=n,with_KEGG=len(kx),frac_KEGG=round(len(kx)/n,3),
        frac_EC=round(sum(any(x.startswith("EC:") for x in f.get("db_xrefs",[])) for f in cds)/n,3),
        frac_nonhyp=round(sum(not f.get("hypothetical") for f in cds)/n,3),
        **{f"frac_{k}":round(v/n,3) for k,v in src.items()},
        with_U=sum("U" in f.get("aa","") for f in cds), pseudo=sum(1 for f in d["features"] if f["type"]=="cds" and f.get("pseudo"))+sum(1 for f in d["features"] if f.get("type")=="pseudogene")))
    import re
    tt=open(f"{SP}/wave1/logs/bakta_{lab}.time").read()
    w=[float(x) for x in re.search(r"Elapsed.*: (\S+)",tt)[1].split(":")]
    summ[-1].update(wall_min=round((w[0]*60+w[1]) if len(w)==2 else (w[0]*3600+w[1]*60+w[2])/60,1) if len(w)==3 else round((w[0]*60+w[1])/60,1),
        cpu_min=round((float(re.search(r"User time.*: (\S+)",tt)[1])+float(re.search(r"System time.*: (\S+)",tt)[1]))/60,1),
        rss_gb=round(int(re.search(r"Maximum resident.*: (\d+)",tt)[1])/1e6,1),
        sec_cds=[(f["locus"],len(f["aa"]),f["product"]) for f in cds if "U" in f.get("aa","")])
    idx={(f["sequence"],f["start"],f["stop"],f["strand"]):f for f in cds}
    seqmap={s["id"]:s["orig_id"] for s in d["sequences"]}
    for g in best.index:
        p=g.split("|")
        if p[0]!="bakta" or p[1]!=lab: continue
        c,se,st=p[2].split(":"); s0,s1=map(int,se.split("-")); st=st or "-"
        f=idx.get((c,s0,s1,st))
        b=best.loc[g]
        if b.ko not in mk.index or (b.score<100 and g not in sig.index): continue
        det.append(dict(label=lab,contig=seqmap[c],start=s0,stop=s1,strand=st,locus=f["locus"] if f else "?",aa=len(f["aa"]) if f else None,
            kofam_call=mk.gene.get(sig.ko[g],sig.ko[g]) if g in sig.index else "-",best_marker=f"{mk.gene[b.ko]}:{b.score:.0f}/{b.thr:.0f}",
            bakta_gene=(f or {}).get("gene") or "",bakta_product=(f or {}).get("product",""),
            bakta_KEGG=",".join(x[5:] for x in (f or {}).get("db_xrefs",[]) if x.startswith("KEGG:")),
            src="(not in final json)" if f is None else "IPS/UPS" if (f.get("ups") or f.get("ips")) else "PSC" if (f.get("psc") or {}).get("valid") else "PSCC" if f.get("pscc") else "none"))
S=pd.DataFrame(summ).set_index("label"); D=pd.DataFrame(det).sort_values(["label","contig","start"])
print(S.T.to_string()); 
sel=D[D.best_marker.str.match(r"(phn|mtt|mtg|cnt|grd|dmd|tor)")]
print(sel[["label","contig","start","aa","kofam_call","best_marker","bakta_gene","bakta_product","bakta_KEGG","src"]].to_string(index=False))
S.to_csv(f"{OUT}/wave1_bakta_summary.tsv",sep="\t"); D.to_csv(f"{OUT}/wave1_bakta_marker_names.tsv",sep="\t",index=False)
