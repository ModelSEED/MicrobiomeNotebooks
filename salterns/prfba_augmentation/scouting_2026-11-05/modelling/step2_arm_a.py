"""Step 2: arm (a) -- a local model of Roseovarius.9 from its own RAST annotation.  Build -> ATPM -> MePn
module -> ATP correction -> gapfill on the cell-5 lactate medium (Pi present) -> MePn protocol (GLPK).
Compares the draft with the original base model (models/<MAG>.contigs__.RAST.json) and the gapfilled model
with the original MePn model (models/MePn_consumers/).  Usage: step2_arm_a.py v7|v6"""
import json
import os
import resource
import sys
import time

os.chdir(os.path.dirname(os.path.abspath(__file__)))
import cobra  # noqa: E402

import common as C  # noqa: E402
import build_pipeline as B  # noqa: E402

version = sys.argv[1] if len(sys.argv) > 1 else "v7"
cobra.Configuration().solver = "glpk"
import modelseedpy  # noqa: E402

rep = {"modelseedpy_file": modelseedpy.__file__, "modelseedpy_version": modelseedpy.__version__,
       "cobra": cobra.__version__, "solver": "glpk"}
t0 = time.perf_counter()
genome, ginfo = C.msgenome_from_rast(C.RAST_DIR / f"{C.MAG}.contigs__.RAST.json")
rep["genome"] = ginfo
rep["t_genome_s"] = round(time.perf_counter() - t0, 2)

# template class from the MAG's own functions (assign once per MAG; share between arms)
try:
    from modelseedpy.helpers import get_classifier
    rep["classifier_class"] = str(get_classifier("knn_ACNP_RAST_filter_01_17_2023").classify(genome))
except Exception as e:
    rep["classifier_class"] = f"error: {str(e)[:120]}"

mid = C.MAG + ".contigs__.RAST"
second = C.M_GENOME_SCALE if (len(sys.argv) > 2 and sys.argv[2] == "pifree") else None
tag = version + ("_pifree" if second else "")
draft, model, brep = B.build_and_gapfill(genome, mid, version=version, second_gapfill_media=second)
rep["build"] = brep
rep["draft_mepn"] = B.mepn_presence(draft)
rep["final_mepn"] = B.mepn_presence(model)
rep["growth"] = B.growth_report(model)

# compare with the original base model and the original MePn model
orig_base = cobra.io.load_json_model(C.BASE_MODELS / f"{C.MAG}.contigs__.RAST.json")
orig_mepn = cobra.io.load_json_model(C.MEPN_MODELS / f"{C.MAG}.contigs.json")


def ids(m, internal=True):
    return {r.id for r in m.reactions if not internal or not r.id.startswith(("EX_", "SK_", "DM_", "bio", "ATPM"))}


def gpr_genes(m):
    return {r.id: frozenset(g.id for g in r.genes) for r in m.reactions}


d0 = cobra.io.load_json_model  # noqa
A, Bset = ids(draft), ids(orig_base)
ga, gb = gpr_genes(draft), gpr_genes(orig_base)
both = A & Bset
rep["draft_vs_original_base"] = {
    "n_internal_draft": len(A), "n_internal_orig_base": len(Bset), "shared": len(both),
    "only_local_draft": sorted(A - Bset), "only_orig_base": sorted(Bset - A),
    "shared_with_identical_gene_set": sum(1 for r in both if ga[r] == gb[r]),
    "n_genes_draft": len(draft.genes), "n_genes_orig_base": len(orig_base.genes)}
F, O = ids(model), ids(orig_mepn)
rep["final_vs_original_mepn_model"] = {
    "n_internal_final": len(F), "n_internal_orig_mepn": len(O), "shared": len(F & O),
    "only_local_final": sorted(F - O), "only_orig_mepn": sorted(O - F),
    "n_reactions_total_final": len(model.reactions), "n_reactions_total_orig_mepn": len(orig_mepn.reactions)}

# MePn protocol on the rebuilt model
ab = C.abundances()
rep["mepn_protocol"] = {}
for mult in (0.07, 0.33):
    t1 = time.perf_counter()
    sim = C.simulate_mepn(model, ab, mult)
    t_sim = time.perf_counter() - t1
    try:
        fva = C.fva_mepn(model, ab, mult)
    except Exception as e:
        fva = {"error": str(e)[:200]}
    rep["mepn_protocol"][str(mult)] = {
        "time_sim_s": round(t_sim, 2),
        "cores": {c: {**sim[c], "CH4_per_abundance": (sim[c]["methane"] / sim[c]["abundance"]) if sim[c]["methane"] is not None else None,
                      "FVA": fva.get(c) if isinstance(fva, dict) else None} for c in C.CORES}}
cobra.io.save_json_model(model, f"arm_a_{tag}_gapfilled.json")
cobra.io.save_json_model(draft, f"arm_a_{tag}_draft.json")
rep["peak_rss_GB"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2, 2)
rep["t_total_s"] = round(time.perf_counter() - t0, 2)
json.dump(rep, open(f"step2_arm_a_{tag}.json", "w"), indent=1, default=str)

b = rep["build"]
print("modelseedpy:", rep["modelseedpy_file"], rep["modelseedpy_version"], "| class:", rep["classifier_class"])
print("genome:", rep["genome"], "| times:", b["times_s"], "| template load", b["t_load_s"], "| total", rep["t_total_s"],
      "s | peak RSS", rep["peak_rss_GB"], "GB")
print("draft rxns", b["draft_n_reactions"], "genes", b["draft_n_genes"], "| after ATP corr", b["n_reactions_after_atp_correction"],
      "| final", b["final_n_reactions"], "| ATP tests", b["atp_tests"], "| core GF", b["atp_core_gapfilling_n"])
for gl in ("gapfill_1", "gapfill_2"):
    g1 = b.get(gl, {}).get("solution") if gl in b else "not run"
    print(gl, "new", g1 if not isinstance(g1, dict) else len(g1["new"]), "reversed", None if not isinstance(g1, dict) else len(g1["reversed"]),
          "| new non-EX:", sorted(k for k in g1["new"] if not k.startswith("EX_")) if (isinstance(g1, dict) and gl == "gapfill_2") else "")
print("MePn module:", b["mepn_module"])
print("draft MePn rxns:", {k: bool(v) for k, v in rep["draft_mepn"].items()})
print("final MePn rxns:", rep["final_mepn"])
print("growth:", rep["growth"])
dv = rep["draft_vs_original_base"]
print("draft vs orig base: local", dv["n_internal_draft"], "orig", dv["n_internal_orig_base"], "shared", dv["shared"],
      "identical GPR gene sets", dv["shared_with_identical_gene_set"], "| only local", len(dv["only_local_draft"]),
      "| only orig", len(dv["only_orig_base"]), dv["only_orig_base"][:30])
fv = rep["final_vs_original_mepn_model"]
print("final vs orig MePn model: local", fv["n_internal_final"], "orig", fv["n_internal_orig_mepn"], "shared", fv["shared"],
      "| only local", len(fv["only_local_final"]), "| only orig", len(fv["only_orig_mepn"]))
for mult, d in rep["mepn_protocol"].items():
    print("MePn protocol multiple", mult, "sim s", d["time_sim_s"])
    for c, v in d["cores"].items():
        fva = v["FVA"]
        print(f"  {c} {v['status']} CH4/a {v['CH4_per_abundance']} MePn {v['MePn']} Pi {v['phosphate']} lact {v['lactate']} "
              f"| FVA CH4 {fva['EX_cpd01024_e0'] if fva else None}")
