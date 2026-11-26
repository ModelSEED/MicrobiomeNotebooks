import json, csv
from collections import defaultdict
H=json.load(open('id_history.json'))
EK='/home/freiburger/Documents/EmilyKin/'
asv=set(json.load(open(EK+'iterativeIDs.json')))
chg={r['MAG']:r for r in csv.DictReader(open('id_change_map.tsv'),delimiter='\t')}
holders=defaultdict(set)   # alias -> {(mag, era)}
for m,h in H['hist'].items():
    for c,i in h.items():
        if i: holders[i].add((m,c))
for m,i in H['head'].items(): holders[i].add((m,'HEAD'))
for m,i in H['new'].items(): holders[i].add((m,'R232'))
# also 6a04588 map stored in mag_guild_classification? check separately
rows=[]
aliases=set()
for m in chg:
    for c,i in H['hist'].get(m,{}).items():
        if i: aliases.add(i)
    aliases.add(chg[m]['old_id']); aliases.add(chg[m]['new_id'])
for a in sorted(aliases):
    hs=sorted(holders[a])
    mags=sorted({x for x,_ in hs})
    rows.append(dict(alias=a, is_asv_id=a in asv, n_mags=len(mags), holders=';'.join(f'{x}@{c}' for x,c in hs)))
with open('alias_holders.tsv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
for r in rows:
    if r['n_mags']>1: print(r)
