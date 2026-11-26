"""Find 'ID near bin' pairs in text files where the ID is the bin's pre-R232 ID (or another bin's MAG ID)."""
import re,json,subprocess,collections,sys,os
R='/home/freiburger/Documents/EmilyKin'
m=json.load(open('ref_map.json'));NEW=m['new'];OLD=m['old']
ch={b for b in NEW if NEW[b]!=OLD[b]}
ids=sorted({OLD[b] for b in ch}|set(NEW.values()),key=len,reverse=True)
idrx='(?:'+'|'.join(map(re.escape,ids))+')'
binrx=r'(?:CAN_[1-5]_bin\.\d+|coasm_bin\.\d+)'
B=r'(?<![A-Za-z0-9_.])'; E=r'(?![0-9A-Za-z_])'
pat=re.compile(rf'{B}(?P<id>{idrx}){E}[^\n]{{0,25}}?(?<![A-Za-z0-9_])(?P<bin>{binrx})(?!\d)|(?<![A-Za-z0-9_])(?P<bin2>{binrx})(?!\d)[^\n]{{0,25}}?{B}(?P<id2>{idrx}){E}')
files=subprocess.check_output(['git','-C',R,'ls-files'],text=True).split('\n')+subprocess.check_output(['git','-C',R,'ls-files','--others','--exclude-standard'],text=True).split('\n')
changed=set(l.strip() for l in open('changed_tracked.txt'))
res=collections.defaultdict(list)
for f in files:
    if not f or not re.search(r'\.(csv|tsv|json|svg|html|htm|nwk|tree|newick|treefile|xml)$',f): continue
    p=os.path.join(R,f)
    if not os.path.isfile(p) or os.path.getsize(p)>50_000_000: continue
    t=open(p,encoding='utf-8',errors='replace').read()
    if not re.search(binrx,t): continue
    for ln,line in enumerate(t.split('\n'),1):
        for mm in pat.finditer(line):
            i=mm.group('id') or mm.group('id2'); b=mm.group('bin') or mm.group('bin2')
            if b not in NEW: continue
            if i==NEW[b]: continue
            if i==OLD[b]: kind='old_of_bin'
            else: kind='other_bins_id'
            res[f].append((ln,kind,i,b,line.strip()[:220]))
out=open('textpairs_tab.tsv','w')
for f in sorted(res):
    k=collections.Counter(x[1] for x in res[f])
    print(('M ' if f in changed else '- ')+f, dict(k))
    for x in res[f]: out.write('\t'.join([f,'M' if f in changed else '-',str(x[0]),x[1],x[2],x[3],x[4]])+'\n')
