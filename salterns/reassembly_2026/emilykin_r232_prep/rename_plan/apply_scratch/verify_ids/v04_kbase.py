"""Read-only KBase check of ws 262319: names joined on bin. Token never printed."""
import json, os, re, collections, urllib.request
T=open(os.path.expanduser("~/.kbase/token")).read().strip()
ALLOWED=("list_objects","get_objects2","get_object_info3")
def rpc(m,p):
    assert m.split(".")[1] in ALLOWED
    r=urllib.request.Request("https://kbase.us/services/ws",data=json.dumps({"version":"1.1","id":"1","method":m,"params":p}).encode(),headers={"Authorization":T,"Content-Type":"application/json"})
    return json.load(urllib.request.urlopen(r,timeout=600))["result"]
m=json.load(open('ref_map.json')); NEW=m['new']; OLD=m['old']
objs=[]; mn=0
while True:
    res=rpc("Workspace.list_objects",[{"ids":[262319],"showHidden":1,"includeMetadata":1,"minObjectID":mn,"maxObjectID":mn+4999}])[0]
    objs+=res
    if mn>2000: break
    mn+=5000
print('objects',len(objs),collections.Counter(o[2].split('-')[0] for o in objs))
names=collections.Counter(o[1] for o in objs); print('dup names',[n for n,c in names.items() if c>1])
asm=[o for o in objs if o[2].startswith('KBaseGenomeAnnotations.Assembly')]
gen=[o for o in objs if o[2].startswith('KBaseGenomes.Genome-')]
sets=[o for o in objs if 'AssemblySet' in o[2]]
# fetch assembly_id for assemblies
refs=[{"ref":f"262319/{o[0]}","included":["assembly_id","std_lineages"]} for o in asm]
data=rpc("Workspace.get_objects2",[{"objects":refs}])[0]["data"]
objid2bin={}; bad=[]; lin_bad=0
import csv
lin={}
for r in csv.DictReader(open('/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/results/gtdbtk.r232.summary.tsv'),delimiter='\t'): lin[r['user_genome']]=r['classification']
for d in data:
    info=d['info']; b=d['data'].get('assembly_id')
    objid2bin[info[0]]=b
    if b not in NEW: bad.append(('asm bin?',info[0],info[1],b)); continue
    if info[1]!=NEW[b]: bad.append(('asm name',info[0],info[1],b,NEW[b]))
    if info[4]!=2: bad.append(('asm ver',info[0],info[4]))
    g=(d['data'].get('std_lineages') or {}).get('gtdb')
    gl=g.get('lineage') if isinstance(g,dict) else g
    if gl is None or (lin.get(b) and lin[b] not in json.dumps(g)): lin_bad+=1
print('assemblies',len(asm),'bins',len(set(objid2bin.values())),'bad',bad[:10],len(bad),'lineage mismatches',lin_bad)
print('example std_lineages',json.dumps(data[0]['data'].get('std_lineages'))[:300])
refs=[{"ref":f"262319/{o[0]}","included":["assembly_ref","scientific_name"]} for o in gen]
gd=[]
for i in range(0,len(refs),15):
    gd+=rpc("Workspace.get_objects2",[{"objects":refs[i:i+15]}])[0]["data"]
gbad=[]; gbins=set()
for d in gd:
    info=d['info']; ar=d['data'].get('assembly_ref','')
    try: aid=int(ar.split('/')[1])
    except: aid=None
    b=objid2bin.get(aid)
    if b is None:
        gbad.append(('genome no asm',info[0],info[1],ar)); continue
    gbins.add(b)
    mm=re.match(r'^emilykin_MAGs_RAST_code\d+(?:_b\d+)?_(.+)$',info[1])
    if not mm or mm.group(1)!=NEW[b]: gbad.append(('genome name',info[0],info[1],b,NEW[b]))
print('genomes',len(gen),'bins',len(gbins),'bad',gbad[:10],len(gbad))
for s in sets:
    print('set',s[0],s[1],'ver',s[4],s[2])
sd=rpc("Workspace.get_objects2",[{"objects":[{"ref":"262319/585"}]}])[0]["data"][0]
items=sd['data']['items']; ibad=[]
for it in items:
    oid=int(it['ref'].split('/')[1]); ver=int(it['ref'].split('/')[2]); b=objid2bin.get(oid)
    if b is None or it.get('label')!=NEW[b] or ver!=2: ibad.append((it,b))
print('set585 ver',sd['info'][4],'items',len(items),'bad',ibad[:5],len(ibad))
# names in ws equal to any old-only id (old ids that are not new ids of any bin)
newset=set(NEW.values()); oldonly={OLD[b] for b in OLD if OLD[b]!=NEW[b] and OLD[b] not in newset}
hits=[(o[0],o[1],o[2]) for o in objs if o[1] in oldonly or any(o[1].endswith('_'+x) for x in oldonly)]
print('objects still named with an old-only ID',len(hits),hits[:10])
json.dump({'n':len(objs)},open('kbase_check.json','w'))
