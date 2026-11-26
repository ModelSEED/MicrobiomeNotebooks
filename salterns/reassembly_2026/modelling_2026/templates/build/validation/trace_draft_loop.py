#!/usr/bin/env python
"""Trace the draft-level ATP loop left in the Roseovarius draft: rebuild the draft on both arms,
print the loop reactions with their template directions, and test whether each arm still
leaks once the loop's reactions are blocked one at a time."""
import json, logging, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import validate_build as V
for n in ("modelseedpy", "cobra", "optlang"):
    logging.getLogger(n).setLevel(logging.ERROR)
import gurobipy; gurobipy.setParam("OutputFlag", 0); gurobipy.setParam("Threads", 2)
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.core.msbuilder import MSBuilder
G = sys.argv[1]; kind = sys.argv[2]
out = {}
for arm in ("upstream", "combined"):
    tj = V.read_template(arm, kind); tmpl = MSTemplateBuilder.from_dict(tj).build()
    core = MSTemplateBuilder.from_dict(V.read_template(arm, "core")).build()
    g, _ = V.genome_from_rast_json(G)
    m = MSBuilder(g, tmpl, ontology_term="RAST", template_core=core).build_base_model("x", "0", allow_all_non_grp_reactions=False, annotate_with_rast=False)
    MSBuilder.add_atpm(m)
    base = V.free_atp(m, V.FREE_MEDIA["water_H"])
    print(f"== {arm}: free ATP {base['atp']}; active {base.get('active')}")
    loop = [r for r in base.get("active", {}) if r != "ATPM_c0"]
    for rid in loop:
        r = m.reactions.get_by_id(rid); tr = {x["id"]: x for x in tj["reactions"]}[rid[:-1]]
        print(f"   {rid} [{tr["direction"]}/{tr["GapfillDirection"]}] gpr={bool(r.genes)} {r.build_reaction_string(use_metabolite_names=True)}")
    # knock each candidate loop reaction (upstream GAPN set + combined lactate set) and retest
    for ko in ["rxn23850_c0", "rxn00145_c0", "rxn00499_c0", "rxn13689_c0", "rxn10122_c0"]:
        if ko in m.reactions:
            with m:
                m.reactions.get_by_id(ko).knock_out()
                v = V.free_atp(m, V.FREE_MEDIA["water_H"])
            print(f"   KO {ko}: free ATP {v['atp']} active {list(v.get('active', {}))[:12]}")
