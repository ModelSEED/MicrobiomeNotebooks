"""Step 3a: arm (b) CODE-PATH TEST with MOCK reference features (synthetic ids 'mockref.N', invented
probabilities; NOT a prFBA augmentation, no existing augmented genome is used).  The same pipeline as arm (a):
  b1 'inclusion'  : MAG features (p = 1) + mock features (p < 1) -> one draft -> ATP correction -> gapfill(s);
                    reaction.notes['probability'] = max p over GPR genes (MAG genes 1, gene-less 0)
                    (codiffusion_bioreactor/build_prfba_community_models.py lines 356-362)
  b2 'gapfill'    : MAG-only draft; mock features enter only as MSGapfill reaction_scores {rxn: {gene: p}}
                    (penalty x (2 - p), modelseedpy/fbapkg/gapfillingpkg.py compute_gapfilling_penalties)
Then the cell-27 MePn protocol (plain FBA, GLPK).  Output: step3a_arm_b_mock.json, arm_b1/b2 model JSONs."""
import json
import os
import time

os.chdir(os.path.dirname(os.path.abspath(__file__)))
import cobra  # noqa: E402

import common as C  # noqa: E402
import build_pipeline as B  # noqa: E402

cobra.Configuration().solver = "glpk"
MOCK = [  # (id, functions, p) -- invented for the code-path test
    ("mockref.1", ["Phosphoribosyl 1,2-cyclic phosphodiesterase (EC 3.1.4.55)"], 0.8),   # PhnP, template name
    ("mockref.2", ["L-lactate permease"], 0.5),                                           # gapfilled in arm (a)
    ("mockref.3", ["Methane monooxygenase component A alpha chain (EC 1.14.13.25)"], 0.1),
    ("mockref.4", ["Methane monooxygenase component A beta chain (EC 1.14.13.25)"], 0.1),
    ("mockref.5", ["Methane monooxygenase component C (EC 1.14.13.25)"], 0.1),
]
P = {fid: p for fid, _, p in MOCK}
rast = C.RAST_DIR / f"{C.MAG}.contigs__.RAST.json"
mid = C.MAG + ".contigs__.RAST"


def set_probabilities(model, gene_p):
    n = 0
    for r in model.reactions:
        if r.genes:
            r.notes["probability"] = str(max(gene_p.get(g.id, 1.0) for g in r.genes))  # MAG genes p = 1
            n += 1
        else:
            r.notes["probability"] = "0.0"
    return n


def mock_reaction_scores():
    """{rxn_base: {gene: p}} from the template's role->reaction map, for the mock features only."""
    from modelseedpy.core.msbuilder import MSBuilder
    from modelseedpy.core.msgenome import MSGenome, MSFeature
    g = MSGenome()
    feats = []
    for fid, fns, _ in MOCK:
        f = MSFeature(fid, "")
        for fn in fns:
            f.add_ontology_term("RAST", fn)
        feats.append(f)
    g.add_features(feats)
    gs_t, _, _ = B.load_templates("v7")
    sets = MSBuilder(g, gs_t).generate_reaction_complex_sets()
    scores = {}
    for rid, cpx in sets.items():
        genes = {gid for roles in cpx.values() for gl in roles.values() for gid in gl}
        scores[rid.rsplit("_", 1)[0]] = {gid: {"probability": P[gid]} for gid in genes}
    return scores


rep = {"mock_features": MOCK}
t0 = time.perf_counter()
genome_b1, gi = C.msgenome_from_rast(rast, extra_features=[(fid, fns) for fid, fns, _ in MOCK])
draft_b1, model_b1, brep_b1 = B.build_and_gapfill(genome_b1, mid + ".b1", "v7", second_gapfill_media=C.M_GENOME_SCALE)
rep["b1"] = {"genome": gi, "build": brep_b1, "time_s": round(time.perf_counter() - t0, 1),
             "n_rxn_with_genes_probability": set_probabilities(model_b1, P),
             "mock_reactions_in_draft": {r.id: r.gene_reaction_rule for r in draft_b1.reactions
                                         if any(g.id.startswith("mockref") for g in r.genes)}}
t0 = time.perf_counter()
scores = mock_reaction_scores()
rep["b2_reaction_scores"] = scores
genome_a, _ = C.msgenome_from_rast(rast)
draft_b2, model_b2, brep_b2 = B.build_and_gapfill(genome_a, mid + ".b2", "v7", reaction_scores=scores,
                                                  second_gapfill_media=C.M_GENOME_SCALE)
rep["b2"] = {"build": brep_b2, "time_s": round(time.perf_counter() - t0, 1),
             "n_rxn_with_genes_probability": set_probabilities(model_b2, P),
             "gapfilled_reactions_given_mock_genes": {r.id: r.gene_reaction_rule for r in model_b2.reactions
                                                      if any(g.id.startswith("mockref") for g in r.genes)}}
ab = C.abundances()
arm_a = cobra.io.load_json_model("arm_a_v7_pifree_gapfilled.json")
rep["mepn_protocol_R1_A"] = {}
for lab, m in (("arm_a", arm_a), ("b1_inclusion", model_b1), ("b2_gapfill_scores", model_b2)):
    for mult in (0.07, 0.33):
        s = C.simulate_mepn(m, ab[["R1_A"]], mult)["R1_A"]
        rep["mepn_protocol_R1_A"][f"{lab}|{mult}"] = {
            "status": s["status"], "CH4_per_biomass": s["methane"] / s["abundance"] if s["methane"] is not None else None,
            "lactate_per_biomass": -s["lactate"] / s["abundance"] if s["lactate"] is not None else None}
rep["growth_b1"] = B.growth_report(model_b1)
rep["mepn_b1"] = B.mepn_presence(model_b1)
cobra.io.save_json_model(model_b1, "arm_b1_mock_v7_gapfilled.json")
cobra.io.save_json_model(model_b2, "arm_b2_mock_v7_gapfilled.json")
json.dump(rep, open("step3a_arm_b_mock.json", "w"), indent=1, default=str)
for lab in ("b1", "b2"):
    b = rep[lab]["build"]
    print(lab, "time", rep[lab]["time_s"], "s | times", b["times_s"], "| draft", b["draft_n_reactions"], "final", b["final_n_reactions"])
    for gl in ("gapfill_1", "gapfill_2"):
        sol = b.get(gl, {}).get("solution")
        print("  ", gl, None if not sol else sorted(k for k in sol["new"] if not k.startswith("EX_")))
print("b1 mock reactions in draft:", rep["b1"]["mock_reactions_in_draft"])
print("b2 reaction_scores:", scores)
print("b2 gapfilled reactions that received mock genes:", rep["b2"]["gapfilled_reactions_given_mock_genes"])
print("b1 growth:", rep["growth_b1"])
for k, v in rep["mepn_protocol_R1_A"].items():
    print(k, v)
