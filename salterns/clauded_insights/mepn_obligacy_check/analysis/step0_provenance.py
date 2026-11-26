"""Gapfill provenance: diff each final model (salterns/models/MePn_consumers/<id>.json, byte-identical to
cliff_mags/mac/salterns/genome_scale/<id>.json) against fliu's pre-gapfill draft
(cliff_mags/mac/salterns/models/<id>__.RAST.json; gapfill.ipynb cell 4 loads './models/{i}__.RAST.json').
Known non-gapfill additions (gapfill.ipynb): ATPM_c0 (cell 3), rxn05602_c0 (cell 74 manual core gapfill),
the MePn module rxs00001_c0/rxn26456_c0/rxn26457_c0/rxs43120_c0/rxn26447_c0/rxn04651_c0 (cell 76),
SK_ drains (cell 92 _add_dm), EX_ exchanges. Everything else that is new = cell-93 gapfill solution.
Plain-JSON diff (no solver).  Output: step0_provenance.json"""
import json
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, "/tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns/a24fbb2e-9d33-5e05-aeea-a07d777ccde7/scratchpad/mepn_verify/analysis")
import vcommon as V  # noqa: E402

MODULE = {"rxs00001_c0", "rxn26456_c0", "rxn26457_c0", "rxs43120_c0", "rxn26447_c0", "rxn04651_c0"}


def rx(path):
    d = json.load(open(path))
    return {r["id"]: r for r in d["reactions"]}, {m["id"]: m for m in d["metabolites"]}


out = {}
for mid in V.model_ids():
    fin, finm = rx(V.MEPN_MODELS / f"{mid}.json")
    dra, dram = rx(V.DRAFTS / f"{mid}__.RAST.json")
    new = sorted(set(fin) - set(dra))
    gone = sorted(set(dra) - set(fin))
    cls = {"ATPM": [], "manual_core_rxn05602": [], "mepn_module": [], "EX": [], "SK": [], "gapfilled": [], "other": []}
    for r in new:
        if r == "ATPM_c0":
            cls["ATPM"].append(r)
        elif r == "rxn05602_c0":
            cls["manual_core_rxn05602"].append(r)
        elif r in MODULE:
            cls["mepn_module"].append(r)
        elif r.startswith("EX_"):
            cls["EX"].append(r)
        elif r.startswith("SK_") or r.startswith("DM_"):
            cls["SK"].append(r)
        elif r.startswith("rxn") or r.startswith("rxs") or r.startswith("sul"):
            cls["gapfilled"].append(r)
        else:
            cls["other"].append(r)
    gf = {r: {"bounds": [fin[r]["lower_bound"], fin[r]["upper_bound"]], "GPR": fin[r]["gene_reaction_rule"],
              "name": fin[r]["name"]} for r in cls["gapfilled"]}
    bchg = {r: {"draft": [dra[r]["lower_bound"], dra[r]["upper_bound"]], "final": [fin[r]["lower_bound"], fin[r]["upper_bound"]]}
            for r in set(fin) & set(dra)
            if (dra[r]["lower_bound"], dra[r]["upper_bound"]) != (fin[r]["lower_bound"], fin[r]["upper_bound"])}
    gchg = sorted(r for r in set(fin) & set(dra) if dra[r]["gene_reaction_rule"] != fin[r]["gene_reaction_rule"])
    schg = sorted(r for r in set(fin) & set(dra) if dra[r]["metabolites"] != fin[r]["metabolites"])
    module_in_draft = sorted(MODULE & set(dra))
    out[mid] = {"iterative_id": V.iterative_id(mid), "n_final": len(fin), "n_draft": len(dra),
                "n_new": len(new), "removed_from_draft": gone,
                "classes": {k: len(v) for k, v in cls.items()}, "other_new": cls["other"],
                "mepn_module_added": cls["mepn_module"], "mepn_module_in_draft": module_in_draft,
                "gapfilled": gf, "bounds_changed_vs_draft": bchg, "gpr_changed": gchg, "stoich_changed": schg,
                "new_metabolites": sorted(set(finm) - set(dram))}
json.dump(out, open(V.OUT / "step0_provenance.json", "w"), indent=1)
for mid, d in out.items():
    print(f"{d['iterative_id']:22s} final {d['n_final']:5d} draft {d['n_draft']:5d} new {d['n_new']:4d} {d['classes']} "
          f"module_in_draft={d['mepn_module_in_draft']} bounds_changed={len(d['bounds_changed_vs_draft'])} "
          f"gpr_changed={len(d['gpr_changed'])} stoich_changed={d['stoich_changed']} removed={len(d['removed_from_draft'])} other={d['other_new']}")
