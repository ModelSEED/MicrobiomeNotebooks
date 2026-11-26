import csv
from collections import defaultdict, Counter
P='/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/results/r220_vs_r232_per_mag.tsv'
ranks='dpcofgs'
tot=Counter(); chg=defaultdict(Counter); r232names=Counter()
rows=list(csv.DictReader(open(P),delimiter='\t'))
for r in rows:
    a=r['r220_classification'].split(';'); b=r['r232_classification'].split(';')
    for i,(x,y) in enumerate(zip(a,b)):
        tot[x]+=1; r232names[y]+=1
        if x!=y: chg[x][y]+=1
out=[]
for x,d in chg.items():
    if x.endswith('__'): continue
    n=tot[x]; nch=sum(d.values())
    out.append((x[0],x[3:],n,nch,'; '.join(f'{y[3:] or "(blank)"}x{c}' for y,c in d.most_common()), r232names.get(x,0)))
out.sort()
with open('r220_name_stats.tsv','w') as fh:
    fh.write('rank\tr220_name\tn_mags_r220\tn_changed\tr232_targets\tn_mags_with_this_name_in_r232\n')
    for o in out: fh.write('\t'.join(map(str,o))+'\n')
for o in out:
    if o[0]!='s': print(*o, sep=' | ')
# also check summary consistency
import os
EK='/home/freiburger/Documents/EmilyKin/'
s=list(csv.DictReader(open(EK+'mag_abundance_summary.tsv'),delimiter='\t'))
cl={r['mag']:r['r220_classification'] for r in rows}
bad=0
for r in s:
    lin=';'.join(f'{k[0].lower()}__{r[k]}' for k in ['Domain','Phylum','Class','Order','Family','Genus','Species'])
    if lin!=cl[r['MAG']]: bad+=1; 
print('summary vs gtdbtk r220 lineage mismatches:',bad)
