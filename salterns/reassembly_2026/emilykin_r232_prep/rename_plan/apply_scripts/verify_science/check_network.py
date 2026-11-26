"""D2: network results unchanged vs HEAD after removing the fields added on purpose; added MAG-ID fields = R232 map by bin."""
import csv, io, json, subprocess, glob, os
R="/home/freiburger/Documents/EmilyKin"
def head(p): return subprocess.run(["git","-C",R,"show",f"HEAD:{p}"],capture_output=True,text=True).stdout
mp=json.load(open(f"{R}/mag_iterativeID_r220_to_r232.json")); new={b:v["new_id"] for b,v in mp.items()}
ADDED={"source_MAG_iterativeID":"source_MAG_bin","target_MAG_iterativeID":"target_MAG_bin","MAG_iterativeID":"MAG_bin",
       "guild_node_MAG_iterativeID":"guild_node_bin"}
files=["network/compositional_networks/out/asv/spieceasi_glasso_edgelist.csv","network/compositional_networks/out/asv/spieceasi_glasso_edgelist_MAGonly.csv",
 "network/compositional_networks/out/asv/spieceasi_mb_edgelist.csv","network/compositional_networks/out/asv/spieceasi_mb_edgelist_MAGonly.csv",
 "network/compositional_networks/out/asv/spieceasi_mag_nodes.csv","network/compositional_networks/out/asv/spieceasi_module_membership.csv",
 "tablesforpaper/modules_membership_89node.csv","tablesforpaper/guild_flanking_all.csv","tablesforpaper/guild_flanking_nonMAG.csv","tablesforpaper/pao_gao_direct_edges.csv"]
for f in files:
    h=list(csv.DictReader(io.StringIO(head(f))))
    c=list(csv.DictReader(open(f"{R}/{f}",newline="")))
    added=[k for k in c[0] if k not in h[0]]
    stripped=[{k:v for k,v in r.items() if k in h[0]} for r in c]
    same = stripped==h and [k for k in c[0] if k in h[0]]==list(h[0])
    # added MAG-ID fields: bin -> R232
    binmap={"source_MAG_iterativeID":["source_MAG_bin","source_bin"],"target_MAG_iterativeID":["target_MAG_bin","target_bin"],
            "MAG_iterativeID":["MAG_bin"],"guild_node_MAG_iterativeID":["guild_node_bin"]}
    bad=0; n=0
    for r in c:
        for a in added:
            bc=[x for x in binmap.get(a,[]) if x in r][0]
            b=r[bc]
            exp=new.get(b,"") if b else ""
            n+=1
            if r[a]!=exp: bad+=1
    nodes=set()
    for r in h:
        for k in ("source","target","node","name","asv_iterativeID","guild_node","partner"):
            if k in r: nodes.add(r[k])
    print(f"{'OK ' if same and not bad else 'BAD'} {f}: rows HEAD={len(h)} now={len(c)} added={added} rest_identical={same} added_cells={n} wrong={bad} distinct_nodes={len(nodes)}")
# louvain json
hj=json.loads(head("network/louvain_module_mag_membership.json")); cj=json.load(open(f"{R}/network/louvain_module_mag_membership.json"))
def strip(o):
    if isinstance(o,dict): return {k:strip(v) for k,v in o.items() if k not in("mag_iterativeID",)}
    if isinstance(o,list): return [strip(x) for x in o]
    return o
diffkeys=[]
def cmp(a,b,path=""):
    if type(a)!=type(b): diffkeys.append(path); return
    if isinstance(a,dict):
        for k in set(a)|set(b):
            if k not in a or k not in b: diffkeys.append(path+"/"+k); continue
            cmp(a[k],b[k],path+"/"+k)
    elif isinstance(a,list):
        if len(a)!=len(b): diffkeys.append(path+"[len]"); return
        for i,(x,y) in enumerate(zip(a,b)): cmp(x,y,f"{path}[{i}]")
    elif a!=b: diffkeys.append(path)
cmp(strip(hj),strip(cj))
print("louvain json diffs after removing mag_iterativeID:",diffkeys)
# membership pairs and mag_iterativeID correctness
cnt=0;bad=0;pairs_h=[];pairs_c=[]
def walk(o,acc):
    if isinstance(o,dict):
        if "MAG" in o and "iterativeID" in o: acc.append(o)
        for v in o.values(): walk(v,acc)
    elif isinstance(o,list):
        for v in o: walk(v,acc)
acc=[];walk(cj,acc);acch=[];walk(hj,acch)
for o in acc:
    if "mag_iterativeID" in o:
        cnt+=1; bad+= o["mag_iterativeID"]!=new[o["MAG"]]
print("louvain member records now",len(acc),"HEAD",len(acch),"with mag_iterativeID",cnt,"wrong",bad)
print("distinct (node,MAG) now==HEAD:",sorted({(o['iterativeID'],o['MAG']) for o in acc})==sorted({(o['iterativeID'],o['MAG']) for o in acch}))
# all-phases / aggregate / per-phase network files untouched
st=subprocess.run(["git","-C",R,"status","--porcelain","--","network"],capture_output=True,text=True).stdout
print("network/ status lines:\n"+st)
for p in ["network/network_summary_p_value_FDR_allphases.json","network/network_module_membership_p_value_FDR_allphases.json"]:
    if os.path.exists(f"{R}/{p}"):
        d=json.load(open(f"{R}/{p}")); print(p, {k:d[k] for k in list(d)[:8] if not isinstance(d[k],(dict,list))})
