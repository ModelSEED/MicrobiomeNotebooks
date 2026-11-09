"""B13. (i) ModelSEED biochemistry DB (local ModelSEEDDatabase, read-only) reversibility of the reactions the argument
uses; (ii) local GramNeg V6 template: EMC roles that have no template reaction attached (orphan roles)."""
import json, glob, sys
sys.dont_write_bytecode = True
want = ["rxn00131", "rxn00690", "rxn00691", "rxn00770", "rxn01300", "rxn00908", "rxn06377", "rxn06600", "rxn00330",
        "rxn00139", "rxn00375", "rxn00371", "rxn04651", "rxn26447", "rxn26456", "rxn26457", "rxn13649"]
print("== ModelSEED DB reversibility ('>' forward only, '<' reverse only, '=' reversible, '?' unknown)")
for f in sorted(glob.glob("/home/freiburger/Documents/ModelSEEDDatabase/Biochemistry/reaction_*.json")):
    for r in json.load(open(f)):
        if r["id"] in want:
            print(f"  {r['id']} {r.get('reversibility')}  {r.get('name')[:55]:55s} {r.get('definition')[:110]}")
tpl = json.load(open("/home/freiburger/Documents/ModelSEEDTemplates/templates/v6.0/GramNegModelTemplateV6.json"))
roles = {r["id"]: r["name"] for r in tpl["roles"]}
used = set()
cpx_roles = {c["id"]: {cr["templaterole_ref"].split("/")[-1] for cr in c["complexroles"]} for c in tpl["complexes"]}
for r in tpl["reactions"]:
    for ref in r["templatecomplex_refs"]:
        used |= cpx_roles.get(ref.split("/")[-1], set())
print("\n== template direction of rxn00131_c / rxn00690_c / rxn13649_c")
for r in tpl["reactions"]:
    if r["id"] in ("rxn00131_c", "rxn00690_c", "rxn13649_c"):
        print(f"  {r['id']} direction {r['direction']} GapfillDirection {r.get('GapfillDirection')} type {r['type']} name {r['name']}")
print("\n== EMC-core roles in the template and whether any template reaction uses them")
for rid, name in roles.items():
    if any(k in name for k in ("Crotonyl-CoA carboxylase/reductase", "Ethylmalonyl-CoA mutase", "Methylsuccinyl-CoA dehydrogenase",
                               "crotonyl-CoA carboxylase/reductase", "Mesaconyl-CoA hydratase", "Ethylmalonyl-CoA epimerase")):
        print(f"  {rid}: '{name}'  -> attached to a template reaction: {rid in used}")
