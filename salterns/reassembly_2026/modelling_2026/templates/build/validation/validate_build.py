#!/usr/bin/env python
"""Build, ATP-correct and gapfill one salterns MAG on one template arm, and test it for ATP
made from nothing at each stage.

Arms
  combined  the four files exported to modelling_2026/templates/ (branch salterns-2026-combined)
  upstream  the same four files at ModelSEED/ModelSEEDTemplates main 2294299, read with git show

Pipeline (the steps of MSBuilder.build_and_gapfill, run one at a time so each stage is measured)
  1. MSGenome from the KBase RAST genome JSON: every feature's 'functions' split with
     modelseedpy.core.rast_client.split_annotation into RAST ontology terms
  2. MSBuilder(genome, template, template_core=Core-V6).build_base_model(annotate_with_rast=False),
     then MSBuilder.add_atpm
  3. MSATPCorrection(core template, load_default_medias()): evaluate_growth_media,
     determine_growth_media, apply_growth_media_gapfilling, expand_model_to_genome_scale,
     build_tests
  4. MSGapfill(template, test_conditions=ATP tests, target bio1).run_gapfilling(complete medium),
     integrated into a copy with MSBuilder.integrate_gapfill_solution + apply_gapfill_reversals

Free-ATP tests (maximise ATPM_c0; every boundary reaction closed for uptake, then only the listed
compounds may enter; pass = ATPM < 1e-6)
  water_H        H2O and H+ only (the prFBA C12 gate and MSATPCorrection's 'empty' test)
  mineral_O2     + O2, NH3, phosphate, sulfate, K+, Na+, Mg2+, Ca2+, Cl-   (no carbon, no donor)
  mineral_anox   the same without O2, with nitrate
Run on the draft, after ATP correction and after gapfilling. Biomass on the same media is also
recorded. When a test is not 0, the reactions carrying the pFBA solution are written out.

    validate_build.py --genome <RAST json> --template-kind gramneg|grampos|archaea \
                      --arm combined|upstream --label NAME --out results/NAME.json
"""
import argparse
import json
import logging
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATES_DIR = os.path.abspath(os.path.join(HERE, "..", ".."))
TEMPLATES_REPO = "/home/freiburger/Documents/ModelSEEDTemplates"
UPSTREAM_REV = "2294299"
KB_COMPLETE = "/home/freiburger/Documents/ModelSEEDDatabase/Media/KBaseMedia.cpd"
FILES = {"gramneg": "GramNegModelTemplateV7.json", "grampos": "GramPosModelTemplateV7.json",
         "archaea": "ArchaeaTemplateV6.json", "core": "Core-V6.json"}
FREE_MEDIA = {
    "water_H": ["cpd00001", "cpd00067"],
    "mineral_O2": ["cpd00001", "cpd00067", "cpd00007", "cpd00013", "cpd00009", "cpd00048",
                   "cpd00205", "cpd00971", "cpd00254", "cpd00063", "cpd00099"],
    "mineral_anox": ["cpd00001", "cpd00067", "cpd00209", "cpd00013", "cpd00009", "cpd00048",
                     "cpd00205", "cpd00971", "cpd00254", "cpd00063", "cpd00099"],
}


def read_template(arm, kind):
    if arm == "combined":
        raw = open(os.path.join(TEMPLATES_DIR, FILES[kind]), "rb").read()
    else:
        raw = subprocess.run(["git", "-C", TEMPLATES_REPO, "show",
                              f"{UPSTREAM_REV}:templates/v7.0/{FILES[kind]}"],
                             capture_output=True, check=True).stdout
    return json.loads(raw)


def genome_from_rast_json(path):
    from modelseedpy.core.msgenome import MSGenome, MSFeature
    from modelseedpy.core.rast_client import split_annotation
    d = json.load(open(path))
    g = MSGenome()
    g.id = d.get("id")
    g.scientific_name = d.get("scientific_name")
    feats, n_roles = [], 0
    for f in d["features"]:
        mf = MSFeature(f["id"], f.get("protein_translation") or "")
        fns = f.get("functions") or ([f["function"]] if f.get("function") else [])
        for fn in fns:
            for role in split_annotation(fn):
                role = role.strip()
                if role:
                    mf.add_ontology_term("RAST", role)
                    n_roles += 1
        feats.append(mf)
    g.add_features(feats)
    return g, {"features": len(feats), "role_terms": n_roles,
               "features_with_roles": sum(1 for f in feats if f.ontology_terms.get("RAST"))}


def model_counts(model):
    ex = [r for r in model.reactions if r.id.startswith("EX_")]
    other_bnd = [r for r in model.boundary if not r.id.startswith("EX_")]
    internal = [r for r in model.reactions if r not in model.boundary and not r.id.startswith("bio")
                and r.id != "ATPM_c0"]
    return {"reactions": len(model.reactions), "metabolites": len(model.metabolites),
            "genes": len(model.genes), "internal_reactions": len(internal),
            "gene_associated": sum(1 for r in internal if r.genes),
            "exchanges": len(ex), "other_boundary": len(other_bnd)}


def ensure_exchange(model, cpd):
    from cobra import Metabolite, Reaction
    mid = f"{cpd}_e0"
    if mid not in model.metabolites:
        model.add_metabolites([Metabolite(mid, None, mid, 0, "e0")])
    ex = f"EX_{mid}"
    if ex not in model.reactions:
        r = Reaction(ex, "Exchange for " + mid, "exchanges", -1000, 1000)
        r.add_metabolites({model.metabolites.get_by_id(mid): -1})
        model.add_reactions([r])
    return model.reactions.get_by_id(ex)


def free_atp(model, cpds, trace=True):
    """Max ATPM_c0 (and bio1) with every boundary closed for uptake except cpds."""
    import cobra.flux_analysis
    out = {}
    with model:
        for r in model.boundary:
            if r.lower_bound < 0:
                r.lower_bound = 0
        for c in cpds:
            ensure_exchange(model, c).bounds = (-1000, 1000)
        model.reactions.ATPM_c0.bounds = (0, 1000)
        model.objective = "ATPM_c0"
        v = model.slim_optimize()
        out["atp"] = 0.0 if v is None or v != v else round(float(v), 6)
        if trace and out["atp"] > 1e-6:
            sol = cobra.flux_analysis.pfba(model)
            fl = {k: round(x, 4) for k, x in sol.fluxes.items() if abs(x) > 1e-6}
            out["n_active"] = len(fl)
            out["active"] = dict(sorted(fl.items(), key=lambda kv: -abs(kv[1]))[:80])
        bio = [r.id for r in model.reactions if r.id.startswith("bio1")]
        if bio:
            model.objective = bio[0]
            b = model.slim_optimize()
            out["biomass"] = 0.0 if b is None or b != b else round(float(b), 6)
    out["pass"] = out["atp"] < 1e-6 and out.get("biomass", 0.0) < 1e-6
    return out


def free_atp_all(model):
    return {k: free_atp(model, v) for k, v in FREE_MEDIA.items()}


def complete_media():
    from modelseedpy import MSMedia
    cpds = [l.strip() for l in open(KB_COMPLETE) if l.strip()]
    m = MSMedia.from_dict({c: 100.0 for c in cpds})
    m.id = m.name = "complete"
    return m


def run_tests(model, tests):
    from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
    res = []
    for t in tests:
        with model:
            MSPackageManager.get_pkg_mgr(model).getpkg("KBaseMediaPkg").build_package(t["media"])
            model.objective = t["objective"]
            v = model.slim_optimize()
            v = 0.0 if v is None or v != v else float(v)
        ok = v <= t["threshold"] + 1e-9 if t["is_max_threshold"] else v >= t["threshold"]
        res.append({"media": t["media"].id, "threshold": round(t["threshold"], 6),
                    "value": round(v, 6), "pass": bool(ok)})
    return res


def library_provenance():
    import cobra
    import modelseedpy
    path = os.path.dirname(os.path.dirname(modelseedpy.__file__))
    head = subprocess.run(["git", "-C", path, "rev-parse", "HEAD"], capture_output=True,
                          text=True).stdout.strip()
    branch = subprocess.run(["git", "-C", path, "rev-parse", "--abbrev-ref", "HEAD"],
                            capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", path, "status", "--porcelain", "--untracked-files=no"],
                           capture_output=True, text=True).stdout.strip() != ""
    out = {"modelseedpy_path": path, "modelseedpy_head": head, "modelseedpy_branch": branch,
           "modelseedpy_dirty": dirty, "cobra": cobra.__version__, "python": sys.executable}
    try:
        import gurobipy
        out["gurobi"] = ".".join(map(str, gurobipy.gurobi.version()))
    except Exception:  # noqa: BLE001
        pass
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--template-kind", required=True, choices=("gramneg", "grampos", "archaea"))
    ap.add_argument("--arm", required=True, choices=("combined", "upstream"))
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=2)
    args = ap.parse_args()
    for name in ("modelseedpy", "cobra", "optlang"):
        logging.getLogger(name).setLevel(logging.ERROR)
    try:
        import gurobipy
        gurobipy.setParam("Threads", args.threads)
        gurobipy.setParam("OutputFlag", 0)
    except Exception:  # noqa: BLE001
        pass
    from modelseedpy.core.mstemplate import MSTemplateBuilder
    from modelseedpy.core.msbuilder import MSBuilder
    from modelseedpy.core.msatpcorrection import MSATPCorrection, load_default_medias
    from modelseedpy import MSGapfill
    from modelseedpy.core.msmodel import get_reaction_constraints_from_direction

    res = {"label": args.label, "arm": args.arm, "template_kind": args.template_kind,
           "genome": args.genome, "libraries": library_provenance(), "stages": {}}
    t0 = time.time()
    tmpl = MSTemplateBuilder.from_dict(read_template(args.arm, args.template_kind)).build()
    core = MSTemplateBuilder.from_dict(read_template(args.arm, "core")).build()
    res["template"] = {"id": tmpl.id, "reactions": len(tmpl.reactions),
                       "core_reactions": len(core.reactions)}
    genome, gstats = genome_from_rast_json(args.genome)
    res["genome_stats"] = gstats

    # 1. draft
    b = MSBuilder(genome, tmpl, ontology_term="RAST", template_core=core)
    model = b.build_base_model(args.label, "0", allow_all_non_grp_reactions=False,
                               annotate_with_rast=False)
    MSBuilder.add_atpm(model)
    res["stages"]["draft"] = {"counts": model_counts(model), "free_atp": free_atp_all(model),
                              "reaction_ids": sorted(r.id for r in model.reactions)}
    res["timing"] = {"draft": round(time.time() - t0, 1)}

    # 2. ATP correction
    t1 = time.time()
    atpc = MSATPCorrection(model, core, load_default_medias(), compartment="c0",
                           atp_hydrolysis_id="ATPM_c0", load_default_medias=False)
    pre = atpc.evaluate_growth_media()
    atpc.determine_growth_media()
    atpc.apply_growth_media_gapfilling()
    atpc.expand_model_to_genome_scale()
    tests = atpc.build_tests()
    gf_stats = {}
    for media, st in atpc.media_gapfill_stats.items():
        gf_stats[media.id] = None if st is None else {
            "new": sorted(st.get("new", {})), "reversed": sorted(st.get("reversed", {}))}
    res["stages"]["atp_correction"] = {
        "atp_before_gapfill": {k: round(float(v or 0), 4) for k, v in pre.items()},
        "media_evaluated": len(pre),
        "media_making_atp_without_gapfill": sorted(k for k, v in pre.items() if (v or 0) > 1e-6),
        "core_gapfill_per_media": gf_stats,
        "selected_media": sorted(m.id for m in atpc.selected_media),
        "core_gapfilled_reactions": [str(x) for x in (atpc.cumulative_core_gapfilling or [])],
        "filtered_noncore": sorted(f"{it[0].id}{it[1]}" for it in atpc.filtered_noncore),
        "tests": [{"media": t["media"].id, "threshold": round(t["threshold"], 6),
                   "is_max_threshold": t["is_max_threshold"]} for t in tests],
        "test_results": run_tests(model, tests),
        "counts": model_counts(model),
        "free_atp": free_atp_all(model),
    }
    res["timing"]["atp_correction"] = round(time.time() - t1, 1)

    # 3. gapfill on complete medium
    t2 = time.time()
    gf = MSGapfill(model, default_gapfill_templates=[tmpl], test_conditions=tests,
                   default_target="bio1")
    cm = complete_media()
    gres = gf.run_gapfilling(cm)
    stage = {}
    if gres is None:
        stage["status"] = "no gapfill solution"
        res["stages"]["gapfilled"] = stage
    else:
        gap_sol = {}
        for rxn_id, d in gres["new"].items():
            if rxn_id[:-1] in tmpl.reactions:
                gap_sol[rxn_id[:-1]] = get_reaction_constraints_from_direction(d)
        mg = model.copy()
        MSBuilder.integrate_gapfill_solution(tmpl, mg, gap_sol)
        rev = MSBuilder.apply_gapfill_reversals(mg, gres.get("reversed"))
        from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
        with mg:
            MSPackageManager.get_pkg_mgr(mg).getpkg("KBaseMediaPkg").build_package(cm)
            mg.objective = "bio1"
            g = mg.slim_optimize()
        stage = {"status": "ok", "new": {k: v for k, v in sorted(gres["new"].items())},
                 "reversed": {k: v for k, v in sorted((gres.get("reversed") or {}).items())},
                 "reversals_applied": rev, "growth_complete": round(float(g or 0), 6),
                 "test_results": run_tests(mg, tests), "counts": model_counts(mg),
                 "free_atp": free_atp_all(mg)}
        res["stages"]["gapfilled"] = stage
    res["timing"]["gapfill"] = round(time.time() - t2, 1)
    res["timing"]["total"] = round(time.time() - t0, 1)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(res, fh, indent=1, default=str)
    print(f"[{args.label}] done in {res['timing']['total']} s")


if __name__ == "__main__":
    main()
