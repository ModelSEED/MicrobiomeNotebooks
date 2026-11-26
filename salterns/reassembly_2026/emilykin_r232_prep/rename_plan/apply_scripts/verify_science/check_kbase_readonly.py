"""READ-ONLY KBase check: ws 262319 Assemblies carry the R232 name and GTDB lineage metadata (get_object_info3 only)."""
import json, os, urllib.request, csv
REPO="/home/freiburger/Documents/EmilyKin"; PREP="/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep"
WS="https://kbase.us/services/ws"
tok=open(os.path.expanduser("~/.kbase/token")).read().strip()
def rpc(method, params):
    body=json.dumps({"version":"1.1","method":method,"params":params,"id":"1"}).encode()
    req=urllib.request.Request(WS,data=body,headers={"Authorization":tok,"Content-Type":"application/json"})
    return json.load(urllib.request.urlopen(req,timeout=120))["result"]
log=json.load(open(f"{REPO}/kbase_mag_rename_log.json"))
recs=log if isinstance(log,list) else (log.get("records") or log.get("mags") or list(log.values()))
print(type(log).__name__, (list(log)[:3] if isinstance(log,dict) else recs[0]))
S={r["user_genome"]:r for r in csv.DictReader(open(f"{PREP}/results/gtdbtk.r232.summary.tsv"),delimiter="\t")}
ID={b:v["new_id"] for b,v in json.load(open(f"{REPO}/mag_iterativeID_r220_to_r232.json")).items()}
bins=sorted(log)
# latest version by objid (no version) to confirm v2 is current
objs=[{"ref":"262319/"+log[b]["upa"].split("/")[1]} for b in bins]
infos=[]
for i in range(0,len(objs),100):
    infos+=rpc("Workspace.get_object_info3",[{"objects":objs[i:i+100],"includeMetadata":1}])[0]["infos"]
bad=[]; mk=set()
for b,inf in zip(bins,infos):
    name,ver,meta=inf[1],inf[4],inf[10] or {}
    mk|=set(meta)
    lin=meta.get("GTDB_lineage") or meta.get("gtdb_lineage")
    exp=S[b]["classification"]
    if name!=ID[b] or ver!=2 or lin is None or lin!=exp: bad.append((b,name,ID[b],ver,(lin or "")[:60]))
print("assemblies checked",len(infos),"problems",len(bad),bad[:5])
print("metadata keys seen:",sorted(k for k in mk if "GTDB" in k.upper()))
sample=infos[bins.index("CAN_4_bin.208")]; print("CAN_4_bin.208:",sample[1],sample[4],{k:v for k,v in (sample[10] or {}).items() if "GTDB" in k.upper()})
# AssemblySet latest
si=rpc("Workspace.get_object_info3",[{"objects":[{"ref":"262319/585"}],"includeMetadata":0}])[0]["infos"][0]
print("AssemblySet 585:",si[1],"v",si[4])
# Genome renames (read-only): every renames.tsv object has its new name; Genomes keep their version
ren=list(csv.DictReader(open(f"{PREP}/rename_plan/kbase/renames.tsv"),delimiter="\t"))
inf=rpc("Workspace.get_object_info3",[{"objects":[{"ref":f"262319/{r['objid']}"} for r in ren],"includeMetadata":0}])[0]["infos"]
bad=[(r["objid"],r["kind"],i[1],r["new_name"],i[4],r["ver"]) for r,i in zip(ren,inf) if i[1]!=r["new_name"] or (r["kind"]=="Genome" and str(i[4])!=r["ver"]) or (r["kind"]=="Assembly" and i[4]!=2)]
print("renamed objects checked",len(ren),"by kind",{k:sum(1 for r in ren if r['kind']==k) for k in {r['kind'] for r in ren}},"problems",bad[:5],len(bad))
wi=rpc("Workspace.get_workspace_info",[{"id":262319}])[0]
print("workspace max objid:",wi[4])
