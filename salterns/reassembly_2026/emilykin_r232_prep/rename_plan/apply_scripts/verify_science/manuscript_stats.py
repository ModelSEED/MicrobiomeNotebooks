"""Recompute manuscript taxonomy numbers from the R232 GTDB-Tk summary and repo tables (independent)."""
import csv, json, io, subprocess, re
from collections import Counter
REPO="/home/freiburger/Documents/EmilyKin"; PREP="/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep"
S={r["user_genome"]:r for r in csv.DictReader(open(f"{PREP}/results/gtdbtk.r232.summary.tsv"),delimiter="\t")}
def lin(s): return {t[:3]:t[3:] for t in s.split(";")}
L={b:lin(r["classification"]) for b,r in S.items()}
print("methods:",Counter(r["classification_method"] for r in S.values()))
print("domains:",Counter(l["d__"] for l in L.values()))
ph=Counter(l["p__"] for l in L.values()); print("phyla n=",len(ph),ph.most_common())
for rk,name in [("c__","class"),("o__","order"),("f__","family"),("g__","genus"),("s__","species")]:
    n=sum(1 for l in L.values() if l.get(rk)); print(f"resolved at {name}: {n}; distinct {len({l[rk] for l in L.values() if l.get(rk)})}; unnamed {291-n}")
red=[b for b,r in S.items() if "RED" in r["classification_method"]]; print("RED novelty:",len(red))
print("warnings:",Counter(r["warnings"] for r in S.values()))
# cases
mp=json.load(open(f"{REPO}/mag_iterativeID_r220_to_r232.json"))
print("cases new:",Counter(v["case_new"] for v in mp.values()),"old:",Counter(v["case_old"] for v in mp.values()))
print("changed IDs:",sum(v["old_id"]!=v["new_id"] for v in mp.values()),"root changed:",sum(v["old_id"].split(".")[0]!=v["new_id"].split(".")[0] for v in mp.values()))
sim={}
for r in csv.DictReader(open(f"{PREP}/results/iterativeid_simulation_r232.tsv"),delimiter="\t"):
    sim[r.get("bin") or r.get("MAG") or r.get("mag")]=r
k=list(sim.values())[0]; print("sim cols:",list(k)[:12])
idcol=[c for c in k if "new" in c.lower() or "r232" in c.lower()]
print("sim id-ish cols:",idcol)
# Accumulibacter etc
for b in ["coasm_bin.185","coasm_bin.250","coasm_bin.347","coasm_bin.476","CAN_4_bin.64","CAN_1_bin.98","coasm_bin.535","CAN_2_bin.85","CAN_4_bin.192","CAN_5_bin.106","CAN_4_bin.236","CAN_3_bin.203","CAN_4_bin.116","CAN_4_bin.208","CAN_5_bin.40","coasm_bin.55","coasm_bin.260","coasm_bin.481","CAN_5_bin.70","coasm_bin.365","CAN_3_bin.221"]:
    r=S[b]; print(b, L[b]["p__"], L[b]["g__"], "|", L[b]["s__"], "| ref", r["closest_genome_reference"], r["closest_genome_ani"], r["closest_genome_af"], "|", r["classification_method"][:40], "| id", mp[b]["old_id"],"->",mp[b]["new_id"])
gen=Counter(l["g__"] for l in L.values())
print("Accumulibacter MAGs:",gen["Accumulibacter"],"Competibacter+_A:",gen["Competibacter"]+gen["Competibacter_A"], {g:c for g,c in gen.items() if "ompetib" in g or "ontend" in g or "Azonexus" in g or "Thauera" in g})
