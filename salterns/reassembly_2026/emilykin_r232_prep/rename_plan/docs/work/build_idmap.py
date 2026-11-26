import csv, json, sys
from pathlib import Path
W = Path(__file__).parent
EK = Path.home()/'Documents/EmilyKin'
SIM = W.parent.parent.parent/'results/iterativeid_simulation_r232.tsv'
asv = set(json.load(open(EK/'iterativeIDs.json')))
hist = {}
for c in ['924fa18','6a04588','ffe668c','408cb39']:
    with open(W/f'mag_summary_{c}.tsv') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            hist.setdefault(r['MAG'],{})[c]=r['mag_iterativeID']
# HEAD
with open(EK/'mag_abundance_summary.tsv') as fh:
    head = {r['MAG']:r for r in csv.DictReader(fh, delimiter='\t')}
sim = list(csv.DictReader(open(SIM), delimiter='\t'))
rows=[]
allold = {}
for r in sim:
    m=r['MAG']; h=hist.get(m,{})
    assert head[m]['mag_iterativeID']==r['committed_id'], (m, head[m]['mag_iterativeID'], r['committed_id'])
    allold[m]=r['committed_id']
out=[]
newset = {r['simulated_r232_id']:r['MAG'] for r in sim}
oldset = {r['committed_id']:r['MAG'] for r in sim}
for r in sim:
    if r['changed']!='True': continue
    m=r['MAG']; h=hist.get(m,{})
    old=r['committed_id']; new=r['simulated_r232_id']
    out.append(dict(MAG=m, old_id=old, new_id=new, case_old=r['case_r220'], case_new=r['case_r232'],
        old_is_asv_id=old in asv, new_is_asv_id=new in asv,
        old_id_held_by_other_mag_after=newset.get(old,'') if newset.get(old,'')!=m else '',
        new_id_held_by_other_mag_before=oldset.get(new,'') if oldset.get(new,'')!=m else '',
        id_924fa18=h.get('924fa18',''), id_6a04588=h.get('6a04588',''), id_ffe668c=h.get('ffe668c',''),
        r220_genus=r['r220_genus'], r232_genus=r['r232_genus']))
with open(W/'id_change_map.tsv','w',newline='') as fh:
    w=csv.DictWriter(fh, fieldnames=list(out[0]), delimiter='\t'); w.writeheader(); w.writerows(out)
for o in out: print('\t'.join(str(v) for v in o.values()))
# also unchanged MAGs whose 6a04588 id differs from HEAD (stale-in-docs ids)
stale=[(m,h.get('6a04588'),allold[m]) for m,h in hist.items() if m in allold and h.get('6a04588') and h.get('6a04588')!=allold[m]]
print('n MAGs where 6a04588 id != HEAD id:', len(stale))
json.dump({'hist':hist,'head':allold,'new':{r['MAG']:r['simulated_r232_id'] for r in sim}}, open(W/'id_history.json','w'), indent=0)
