import json,re,os,subprocess,collections,csv
R='/home/freiburger/Documents/EmilyKin'
m=json.load(open('ref_map.json'));NEW=m['new'];OLD=m['old']
BIN=re.compile(r'^(CAN_[1-5]_bin\.\d+|coasm_bin\.\d+)$')
files=[f for f in subprocess.check_output(['git','-C',R,'ls-files'],text=True).split('\n')+subprocess.check_output(['git','-C',R,'ls-files','--others','--exclude-standard'],text=True).split('\n') if f.endswith('.json')]
res=collections.defaultdict(collections.Counter)
def walk(o,f):
    if isinstance(o,dict):
        bins={k:v for k,v in o.items() if isinstance(v,str) and BIN.match(v)}
        if len(bins)>=2:
            for k,v in o.items():
                if not isinstance(v,str) or k in bins: continue
                mm=re.match(r'^(.*?)(MAG_iterativeID|mag_iterativeID|iterativeID|MAG_ID)$',k)
                if not mm: continue
                pre=mm.group(1)
                cands=[bk for bk in bins if bk.startswith(pre) and pre]
                if len(cands)!=1: res[f]['ambiguous:'+k]+=1; continue
                b=bins[cands[0]]
                res[f]['ok' if v==NEW.get(b) else ('old' if v==OLD.get(b) else 'other')+':'+k]+=1
        for v in o.values(): walk(v,f)
    elif isinstance(o,list):
        for v in o: walk(v,f)
for f in files:
    p=os.path.join(R,f)
    if not os.path.isfile(p) or os.path.getsize(p)>80_000_000: continue
    raw=open(p,'rb').read()
    if not re.search(rb'(CAN_[1-5]_bin\.\d|coasm_bin\.\d)',raw): continue
    try: o=json.loads(raw)
    except Exception: continue
    walk(o,f)
for f,c in res.items(): print(f,dict(c))
