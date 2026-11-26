"""Generic bin-keyed MAG-ID column scanner (read-only on the repo)."""
import csv, json, os, re, sys, subprocess, collections, io
csv.field_size_limit(10**9)
R='/home/freiburger/Documents/EmilyKin'
PLAN='/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan'
m=json.load(open('ref_map.json')); NEW=m['new']; OLD=m['old']
BIN_RE=re.compile(r'^(CAN_[1-5]_bin\.\d+|coasm_bin\.\d+)$')
# all-era MAG IDs
ALL=set(NEW.values())|set(OLD.values())
for r in csv.DictReader(open(f'{PLAN}/docs/alias_holders.tsv'),delimiter='\t'): ALL.add(r['alias'])
for r in csv.DictReader(open(f'{PLAN}/docs/id_change_map.tsv'),delimiter='\t'):
    for k in ('id_924fa18','id_6a04588','id_ffe668c'):
        if r[k]: ALL.add(r[k])
hc=subprocess.check_output(['git','-C',R,'show','HEAD:mag_iterativeID_old_to_new.csv']).decode().splitlines()
for r in csv.DictReader(hc):
    for k in ('old_mag_iterativeID','old_intersection_id'):
        if r[k]: ALL.add(r[k])
ni=json.load(open(f'{R}/scripts/_new_ids.json'))
ASV=set(ni['asv_old2new'].values())
ALL.discard('')
def cat(v,b):
    if b not in NEW: return 'bin_not_in_291'
    if v=='' or v is None: return 'blank'
    if v==NEW[b]: return 'eq_new'
    if v==OLD[b]: return 'eq_old_changed'
    if v in ALL: return 'other_mag_id'
    return 'non_mag'
def analyse_table(rows, fname, sheet=''):
    out=[]
    if not rows: return out
    ncol=max(len(r) for r in rows)
    cols=[[ (r[i] if i<len(r) else '') for r in rows] for i in range(ncol)]
    cols=[[('' if v is None else str(v).strip()) for v in c] for c in cols]
    bincols=[i for i,c in enumerate(cols) if sum(bool(BIN_RE.match(v)) for v in c)>=1 and sum(bool(BIN_RE.match(v)) for v in c)>=0.5*max(1,sum(v!='' for v in c[1:]))]
    if not bincols: return out
    for i,c in enumerate(cols):
        if i in bincols: continue
        nid=sum(v in ALL for v in c)
        if nid==0: continue
        best=None
        for bi in bincols:
            cnt=collections.Counter()
            pairs=collections.defaultdict(set); rev=collections.defaultdict(set)
            for v,b in zip(c,cols[bi]):
                if not BIN_RE.match(b): continue
                k=cat(v,b); cnt[k]+=1
                if k!='blank': pairs[b].add(v); rev[v].add(b)
            score=cnt['eq_new']+cnt['eq_old_changed']
            if best is None or score>best[0]: best=(score,bi,cnt,pairs,rev)
        score,bi,cnt,pairs,rev=best
        multi_id=sum(len(s)>1 for s in pairs.values()); shared=[v for v,s in rev.items() if len(s)>1 and v in ALL]
        out.append(dict(file=fname,sheet=sheet,col=i,header=cols[i][0],bincol=cols[bi][0],n_ids=nid,
            eq_new=cnt['eq_new'],eq_old_changed=cnt['eq_old_changed'],other_mag_id=cnt['other_mag_id'],non_mag=cnt['non_mag'],blank=cnt['blank'],
            bin_not_in_291=cnt['bin_not_in_291'],bins_with_multi_ids=multi_id,ids_shared_by_bins=len(shared),
            ex_old=';'.join(sorted({f'{b}:{v}' for v,b in zip(c,cols[bi]) if BIN_RE.match(b) and cat(v,b) in ('eq_old_changed','other_mag_id')})[:6]),
            ex_shared=';'.join(shared[:4])))
    return out
def json_walk(o, path, ctx, recs):
    if isinstance(o,dict):
        keys=list(o.keys())
        binkeyed= keys and sum(bool(BIN_RE.match(str(k))) for k in keys)>=0.5*len(keys)
        local=ctx
        if not binkeyed:
            bvals=[v for v in o.values() if isinstance(v,str) and BIN_RE.match(v)]
            if len(set(bvals))==1: local=bvals[0]
        for k,v in o.items():
            p=path+'/'+('{bin}' if binkeyed and BIN_RE.match(str(k)) else str(k))
            c2=str(k) if binkeyed and BIN_RE.match(str(k)) else local
            if isinstance(v,str):
                if v in ALL and c2: recs.append((p,v,c2))
            else: json_walk(v,p,c2,recs)
            if binkeyed and BIN_RE.match(str(k)) and isinstance(v,str): pass
    elif isinstance(o,list):
        for v in o:
            if isinstance(v,str):
                if v in ALL and ctx: recs.append((path+'[]',v,ctx))
            else: json_walk(v,path+'[]',ctx,recs)
def analyse_json(fname):
    try: o=json.load(open(os.path.join(R,fname)))
    except Exception as e: return []
    recs=[]; json_walk(o,'',None,recs)
    agg=collections.defaultdict(list)
    for p,v,b in recs: agg[p].append((v,b))
    out=[]
    for p,l in agg.items():
        cnt=collections.Counter(cat(v,b) for v,b in l)
        pairs=collections.defaultdict(set); rev=collections.defaultdict(set)
        for v,b in l: pairs[b].add(v); rev[v].add(b)
        out.append(dict(file=fname,sheet='',col=p,header=p,bincol='ctx',n_ids=len(l),eq_new=cnt['eq_new'],eq_old_changed=cnt['eq_old_changed'],
          other_mag_id=cnt['other_mag_id'],non_mag=cnt['non_mag'],blank=cnt['blank'],bin_not_in_291=cnt['bin_not_in_291'],
          bins_with_multi_ids=sum(len(s)>1 for s in pairs.values()),ids_shared_by_bins=sum(len(s)>1 for s in rev.values()),
          ex_old=';'.join(sorted({f'{b}:{v}' for v,b in l if cat(v,b) in ('eq_old_changed','other_mag_id')})[:6]),ex_shared=''))
    return out
def main(listfile, outfile):
    files=[l.rstrip('\n') for l in open(listfile) if l.strip()]
    res=[]; n=0
    for f in files:
        p=os.path.join(R,f)
        if not os.path.isfile(p): continue
        sz=os.path.getsize(p)
        if sz>80_000_000: continue
        ext=f.rsplit('.',1)[-1].lower() if '.' in f else ''
        try:
            if ext in ('csv','tsv','txt','tab'):
                raw=open(p,'rb').read()
                if b'\0' in raw[:4096]: continue
                txt=raw.decode('utf-8','replace')
                if not re.search(r'(CAN_[1-5]_bin\.\d|coasm_bin\.\d)',txt): continue
                delim=',' if ext=='csv' else '\t'
                if ext=='txt' and '\t' not in txt[:2000]: continue
                rows=list(csv.reader(io.StringIO(txt),delimiter=delim))
                res+=analyse_table(rows,f); n+=1
            elif ext=='json':
                raw=open(p,'rb').read()
                if not re.search(rb'(CAN_[1-5]_bin\.\d|coasm_bin\.\d)',raw): continue
                res+=analyse_json(f); n+=1
            elif ext in ('xlsx','xlsm'):
                import openpyxl
                wb=openpyxl.load_workbook(p,read_only=True,data_only=False)
                for ws in wb.worksheets:
                    rows=[list(r) for r in ws.iter_rows(values_only=True)]
                    res+=analyse_table(rows,f,ws.title)
                n+=1
        except Exception as e:
            print('ERR',f,e,file=sys.stderr)
    keys=['file','sheet','col','header','bincol','n_ids','eq_new','eq_old_changed','other_mag_id','non_mag','blank','bin_not_in_291','bins_with_multi_ids','ids_shared_by_bins','ex_old','ex_shared']
    with open(outfile,'w',newline='') as fh:
        w=csv.DictWriter(fh,keys,delimiter='\t'); w.writeheader(); w.writerows(res)
    print('files scanned',n,'columns',len(res))
if __name__=='__main__': main(sys.argv[1],sys.argv[2])
