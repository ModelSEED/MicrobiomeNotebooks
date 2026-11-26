"""B5. Structure of the original models (models/MePn_consumers/*.json, read-only) and the local GramNeg V6 template
for the reactions the biological argument turns on: C1 outlets, Thr/Gly routes, PRPP routes, 5'-dA fate, EMC."""
import sys, json
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_bio import *
import cobra
cobra.Configuration().solver = "glpk"

KEY = {
 "rxn00692_c0": "GlyA SHMT", "rxn00541_c0": "LtaE Thr aldolase", "rxn01300_c0": "ThrB homoserine kinase",
 "rxn01069_c0": "ThrC Thr synthase", "rxn00424_c0": "SGAT Ser-glyoxylate AT", "rxn00272_c0": "AGT Ala-glyoxylate AT",
 "rxn00330_c0": "malate synthase", "rxn00336_c0": "isocitrate lyase", "rxn01068_c0": "Thr dehydrogenase",
 "rxn00274_c0": "Kbl", "rxn06377_c0": "GcvP", "rxn06600_c0": "GcvT", "rxn00908_c0": "GCS overall",
 "rxn00690_c0": "Fhs formate-THF ligase", "rxn00691_c0": "PurU", "rxn00907_c0": "FolD dehydrogenase",
 "rxn01211_c0": "FolD cyclohydrolase", "rxn00371_c0": "FDH (NAD)", "rxn00375_c0": "HutG N-formylGlu deformylase",
 "rxn00770_c0": "Prs PRPP synthetase", "rxn00131_c0": "Amn AMP nucleosidase", "rxn00139_c0": "APRT",
 "rxn04651_c0": "PhnN", "rxn26447_c0": "PhnP", "rxn26456_c0": "PhnGHIL", "rxn26457_c0": "PhnM",
 "rxs43120_c0": "PhnJ (hand-built)", "rxs00001_c0": "MePn ABC (hand-built)", "SK_cpd03091_c0": "5'-dA sink",
 "EX_cpd00047_e0": "formate exchange", "EX_cpd00033_e0": "glycine exchange",
}
inv = {iid(m): m for m in model_ids()}
out = {}
for t in ["Roseovarius.9", "Roseovarius.11", "Albimonas.1", "Roseovarius.8", "Roseovarius.1"]:
    m = cobra.io.load_json_model(str(MEPN_MODELS / f"{inv[t]}.json"))
    print(f"===== {t} ({len(m.reactions)} rxns)")
    rec = {}
    for rid, lab in KEY.items():
        if rid in m.reactions:
            r = m.reactions.get_by_id(rid)
            rec[rid] = {"lb": r.lower_bound, "ub": r.upper_bound, "gpr": r.gene_reaction_rule}
            gpr = r.gene_reaction_rule
            gpr = (gpr[:60] + "...") if len(gpr) > 60 else gpr
            print(f"  {rid:16s} {lab:28s} bounds ({r.lower_bound:g},{r.upper_bound:g})  GPR: {gpr.replace('Salt_Pond_','') or '(none)'}")
        else:
            rec[rid] = None
            print(f"  {rid:16s} {lab:28s} ABSENT")
    # every reaction touching 5'-deoxyadenosine
    if "cpd03091_c0" in m.metabolites:
        print("  reactions with 5'-deoxyadenosine (cpd03091_c0):")
        for r in m.metabolites.cpd03091_c0.reactions:
            print(f"     {r.id:16s} ({r.lower_bound:g},{r.upper_bound:g}) GPR[{bool(r.gene_reaction_rule)}] {r.build_reaction_string(True)[:120]}")
    # N-formyl-glutamate producers
    for met in m.metabolites:
        if met.name.lower().startswith("n-formyl-l-glutamate"):
            print(f"  reactions with {met.id} {met.name}:")
            for r in met.reactions:
                print(f"     {r.id:16s} ({r.lower_bound:g},{r.upper_bound:g}) GPR[{bool(r.gene_reaction_rule)}] {r.build_reaction_string(True)[:120]}")
    emc = [r for r in m.reactions if any(k in r.name.lower() or k in r.build_reaction_string(True).lower() for k in
           ("ethylmalonyl", "methylsuccinyl", "mesaconyl", "malyl-coa", "crotonyl-coa", "methylmalyl", "citramalyl-coa"))]
    print(f"  EMC-type reactions in model ({len(emc)}):")
    for r in emc:
        print(f"     {r.id:16s} ({r.lower_bound:g},{r.upper_bound:g}) GPR[{bool(r.gene_reaction_rule)}] {r.build_reaction_string(True)[:120]}")
    out[t] = rec
json.dump(out, open(VB / "b5_model_structure.json", "w"), indent=1)

# Template (local GramNeg V6, the family fliu used) directions/roles for the same reactions
tpl = json.load(open("/home/freiburger/Documents/ModelSEEDTemplates/templates/v6.0/GramNegModelTemplateV6.json"))
roles = {r["id"]: r["name"] for r in tpl["roles"]}
cpx = {c["id"]: [roles.get(cr["templaterole_ref"].split("/")[-1], "?") for cr in c["complexroles"]] for c in tpl["complexes"]}
trx = {r["id"]: r for r in tpl["reactions"]}
print("\n===== local template GramNegModelTemplateV6.json")
for rid in ["rxn01300_c", "rxn00690_c", "rxn00691_c", "rxn00131_c", "rxn00770_c", "rxn06377_c", "rxn06600_c",
            "rxn00908_c", "rxn00330_c", "rxn00139_c", "rxn00375_c"]:
    r = trx.get(rid)
    if not r:
        print(f"  {rid}: not in template"); continue
    rl = sorted({x for ref in r["templatecomplex_refs"] for x in cpx.get(ref.split('/')[-1], [])})
    print(f"  {rid}: type {r['type']}, direction {r['direction']}, GapfillDirection {r.get('GapfillDirection')}, roles {rl}")
# template reactions consuming 5'-dA other than sinks
print("  template reactions containing cpd03091_c (5'-dA):")
for r in tpl["reactions"]:
    if any(rr["templatecompcompound_ref"].split("/")[-1] == "cpd03091_c" for rr in r["templateReactionReagents"]):
        rl = sorted({x for ref in r["templatecomplex_refs"] for x in cpx.get(ref.split('/')[-1], [])})
        print(f"     {r['id']} dir {r['direction']} type {r['type']} name {r['name'][:60]} roles {rl[:3]}")
print("  template drains:", [d for d in tpl.get("drain_list", {})][:40] if isinstance(tpl.get("drain_list"), dict) else tpl.get("drain_list"))
