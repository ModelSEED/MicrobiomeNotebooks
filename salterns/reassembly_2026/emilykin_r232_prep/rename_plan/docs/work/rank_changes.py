import csv
from collections import defaultdict, Counter
P='/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/results/r220_vs_r232_per_mag.tsv'
ranks=['d','p','c','o','f','g','s']
ch=defaultdict(lambda: defaultdict(set))
for r in csv.DictReader(open(P),delimiter='\t'):
    a=r['r220_classification'].split(';'); b=r['r232_classification'].split(';')
    for i,(x,y) in enumerate(zip(a,b)):
        if x!=y:
            ch[(ranks[i],x)][y].add(r['mag'])
out=[]
for (rk,x),d in sorted(ch.items()):
    for y,ms in d.items():
        out.append((rk,x,y,len(ms),','.join(sorted(ms))))
with open('rank_name_changes.tsv','w') as fh:
    fh.write('rank\tr220\tr232\tn_mags\tmags\n')
    for o in out: fh.write('\t'.join(map(str,o))+'\n')
for o in out:
    if o[0]!='s': print(o[0],o[1],'->',o[2],o[3])
print('species changes:',sum(1 for o in out if o[0]=='s'))
for o in out:
    if o[0]=='s': print(o[1],'->',o[2],o[4])
