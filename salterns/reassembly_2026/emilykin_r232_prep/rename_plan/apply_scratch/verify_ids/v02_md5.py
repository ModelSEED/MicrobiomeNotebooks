import csv, hashlib, os, collections
R='/home/freiburger/Documents/EmilyKin'
A='/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan/apply_log'
rows=list(csv.DictReader(open(f'{A}/asv_md5_before.tsv'),delimiter='\t'))
print(collections.Counter(r['asv_keyed'] for r in rows))
print(collections.Counter(r['action'] for r in rows).most_common(20))
bad=[];missing=[]
for r in rows:
    p=os.path.join(R,r['file'])
    if not os.path.exists(p): missing.append(r['file']); continue
    h=hashlib.md5(open(p,'rb').read()).hexdigest()
    if h!=r['md5']: bad.append((r['file'],r['asv_keyed'],r['action'],r['note']))
print('missing',missing)
print('changed',len(bad))
for b in bad: print(b)
