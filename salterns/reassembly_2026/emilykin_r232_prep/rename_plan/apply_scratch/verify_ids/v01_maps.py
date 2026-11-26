import csv, json, sys, subprocess
R='/home/freiburger/Documents/EmilyKin'
P='/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep'
sim={}
for r in csv.DictReader(open(f'{P}/results/iterativeid_simulation_r232.tsv'),delimiter='\t'):
    sim[r['MAG']]=r
print('sim rows',len(sim))
new={b:r['simulated_r232_id'] for b,r in sim.items()}
old={b:r['committed_id'] for b,r in sim.items()}
print('changed',sum(new[b]!=old[b] for b in sim), 'root changed', sum(new[b].split('.')[0]!=old[b].split('.')[0] for b in sim))
print('dup new',len(new)-len(set(new.values())),'dup old',len(old)-len(set(old.values())))
from collections import Counter
print('cases r232',Counter(r['case_r232'] for r in sim.values()),'r220',Counter(r['case_r220'] for r in sim.values()))
# old committed at HEAD
head_json=json.loads(subprocess.check_output(['git','-C',R,'show','HEAD:mag_iterativeID_old_to_new.json']))
print('HEAD old_to_new == committed_id', sum(head_json.get(b)==old[b] for b in sim), len(head_json))
# repo map files
j=json.load(open(f'{R}/mag_iterativeID_r220_to_r232.json'))
bad=[b for b in sim if j.get(b,{}).get('new_id')!=new[b] or j[b]['old_id']!=old[b] or j[b]['case_new']!=sim[b]['case_r232'] or j[b]['case_old']!=sim[b]['case_r220']]
print('r220_to_r232.json rows',len(j),'bad',bad)
c={r['bin']:r for r in csv.DictReader(open(f'{R}/mag_iterativeID_r220_to_r232.csv'))}
bad=[b for b in sim if c.get(b,{}).get('new_id')!=new[b] or c[b]['old_id']!=old[b] or c[b]['case_new']!=sim[b]['case_r232'] or c[b]['case_old']!=sim[b]['case_r220']]
print('r220_to_r232.csv rows',len(c),'bad',bad)
o=json.load(open(f'{R}/mag_iterativeID_old_to_new.json'))
print('old_to_new.json rows',len(o),'bad',[b for b in sim if o.get(b)!=new[b]], 'extra',set(o)-set(sim))
oc={r['MAG']:r for r in csv.DictReader(open(f'{R}/mag_iterativeID_old_to_new.csv'))}
print('old_to_new.csv rows',len(oc),'bad',[b for b in sim if oc[b]['new_mag_iterativeID']!=new[b] or oc[b]['assignment_case']!=sim[b]['case_r232']])
# compare old_to_new.csv other columns vs HEAD
hc={r['MAG']:r for r in csv.DictReader(subprocess.check_output(['git','-C',R,'show','HEAD:mag_iterativeID_old_to_new.csv']).decode().splitlines())}
print('old cols same as HEAD', sum(oc[b]['old_mag_iterativeID']==hc[b]['old_mag_iterativeID'] and oc[b]['old_intersection_id']==hc[b]['old_intersection_id'] for b in sim))
ni=json.load(open(f'{R}/scripts/_new_ids.json'))
print('_new_ids keys',list(ni.keys()))
for k,v in ni.items():
    if isinstance(v,dict): print(k,len(v))
mn=ni.get('mag_new')
if mn: print('mag_new bad',[b for b in sim if mn.get(b)!=new[b]])
hni=json.loads(subprocess.check_output(['git','-C',R,'show','HEAD:scripts/_new_ids.json']))
for k in ni:
    if k!='mag_new': print(k,'same as HEAD',ni[k]==hni[k])
pm=json.load(open(f'{P}/rename_plan/data/mag_id_map_r220_r232.json'))
print('plan map == sim',all(pm[b]==new[b] for b in sim),len(pm))
json.dump({'new':new,'old':old},open('ref_map.json','w'))
