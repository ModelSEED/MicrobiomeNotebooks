import json, csv, re
from collections import defaultdict, Counter
EK='/home/freiburger/Documents/EmilyKin/'
asv=set(json.load(open(EK+'iterativeIDs.json')))
H=json.load(open('id_history.json'))
chg={r['MAG']:r for r in csv.DictReader(open('id_change_map.tsv'),delimiter='\t')}
holders=defaultdict(set)
for m,h in H['hist'].items():
    for c,i in h.items():
        if i: holders[i].add((m,c))
for m,i in H['head'].items(): holders[i].add((m,'HEAD'))
NEW={m:i for m,i in H['new'].items()}
# aliases of changing MAGs: kind old/hist/new
alias_kind={}
for m,r in chg.items():
    alias_kind.setdefault(r['old_id'],set()).add('old')
    alias_kind.setdefault(r['new_id'],set()).add('new')
    for c,i in H['hist'].get(m,{}).items():
        if i and i!=r['old_id']: alias_kind.setdefault(i,set()).add('hist')
def pat(a):
    if a.startswith('Ca_'):
        core=re.escape(a[3:])
        return r'(?:(?<![A-Za-z0-9_])(?:[*_]?(?:Ca|Candidatus)[._]?[*_]?\s*))?(?<![A-Za-z0-9])'+core+r'(?![0-9A-Za-z_])'
    return r'(?<![A-Za-z0-9_])'+re.escape(a)+r'(?![0-9A-Za-z_])'
# longest first so 'X.1_m' not shadowed
pats=[(a,re.compile(pat(a))) for a in sorted(alias_kind,key=len,reverse=True)]
BIN=re.compile(r'(?<![A-Za-z0-9_])(?:CAN_\d_bin\.\d+|coasm_bin\.\d+)(?![0-9])')
ASVCUE=re.compile(r'\bASVs?\b|amplicon|16S|MiDAS|SILVA|iterativeIDs\b',re.I)
MAGCUE=re.compile(r'\bMAGs?\b|genome|\bbin\b|GTDB|CheckM|completeness',re.I)
out=[]
for line in open('segments.jsonl'):
    s=json.loads(line)
    t=s['text']
    if not t: continue
    found=[]
    for a,p in pats:
        for mo in p.finditer(t):
            # skip if overlapped by an already-found longer alias
            if any(mo.start()<e and mo.end()>b for b,e,_ in found): continue
            found.append((mo.start(),mo.end(),a))
    if not found: continue
    ctxbins=set(BIN.findall(t))|set(BIN.findall(s.get('row','') or ''))
    for b,e,a in sorted(found):
        hold={m for m,_ in holders.get(a,())}
        newhold={m for m,i in NEW.items() if i==a}
        kinds=alias_kind[a]
        isasv=a in asv
        hb=sorted(ctxbins & (hold|newhold))
        if hb: ctx='MAG:bin-confirmed'
        elif not isasv and ctxbins: ctx='MAG:non-ASV id, other bins in row ('+','.join(sorted(ctxbins))[:60]+')'
        elif not isasv: ctx='MAG:non-ASV id'
        else:
            sent=t
            if ASVCUE.search(sent) and not MAGCUE.search(sent): ctx='ASV?:asv cue'
            elif MAGCUE.search(sent) and not ASVCUE.search(sent): ctx='MAG?:mag cue'
            elif ctxbins: ctx='?:other bins in row ('+','.join(sorted(ctxbins))[:60]+')'
            else: ctx='ambiguous'
        q=t[max(0,b-120):e+120].replace('\t',' ').replace('\n',' ')
        out.append(dict(file=s['file'],loc=json.dumps(s['loc'],ensure_ascii=False),state=s['state'],alias=a,match=t[b:e],
            alias_kind='/'.join(sorted(kinds)),alias_is_asv_id=isasv,
            alias_holders=';'.join(sorted(f'{m}@{c}' for m,c in holders.get(a,()))) + (';'+';'.join(f'{m}@R232' for m in sorted(newhold)) if newhold else ''),
            ctx_bins=','.join(sorted(ctxbins))[:200],auto_context=ctx,quote=q))
with open('id_occurrences_raw.tsv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(out[0]),delimiter='\t'); w.writeheader(); w.writerows(out)
c=Counter((o['file']) for o in out)
for k,v in sorted(c.items(),key=lambda x:-x[1]): print(v,k)
