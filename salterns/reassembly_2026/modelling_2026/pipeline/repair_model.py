"""Final-set repairs of one winning-arm model (post-processing; the comparison builds are not changed).

Applied identically to every model of the final set, after the pre-registered comparison was decided:
 1. Ammonia-oxidation block. QC of the comparison builds found models that make ATP from NH3 + O2 (or NH3 +
    nitrate) with no carbon source: a nitrite reductase run in the oxidising direction (rxn13647, ubiquinol:nitrite,
    used NH3 -> nitrite; rxn03978, cytochrome c:nitrite, written NH3 -> nitrite) turns NH3 into a respiratory electron
    donor. A first repair pass (v1) showed the re-gap-fill then finds the next such route: hydroxylamine reductase
    rxn00109 run NH3 -> NH2OH with NAD and rxn00570 run NH2OH + O2 -> nitrite. None of these organisms is a
    nitrifier, and ModelSEEDpy's ATP-test media contain no NH3, so ATP correction cannot see it. So (a) the four
    reactions are restricted to their reductive direction in the model and in the in-memory gap-filling template
    (NH3_OX_BLOCK), and (b) every repair gap-fill also uses the three carbon-free free-ATP media as ATP tests
    (ATPM <= 1e-5), so that no candidate may make ATP from NH3 + O2 / nitrate. The template files are not changed.
 1b. No gap-filled methane consumption. In the energy-limited Pi-free MePn stage (lactate 1) the gap-filler can add
    a way to oxidise the CH4 that the C-P lyase (rxn26459) releases in the cytosol, so that CH4 export is 0: a methane
    monooxygenase (rxn00843; repair v2, 3 models) or reversed methanogenesis (rxn03127 methyl-CoM reductase run
    CH4 -> methyl-CoM; repair v3, 1 model). No annotation supports either, and the original protocol's hand-built lyase
    rxs43120 releases CH4 outside the cell. So every template reaction with methane (except the transport rxn10471)
    is restricted, in the in-memory gap-filling template, to its methane-producing direction, and in the model every
    gap-filled reaction loses its methane-consuming direction (removed if nothing is left). Methane-consuming
    reactions that come from the annotation (draft, gene-associated) are left as annotated and reported.
 2. Growth restoration. Two arm-A builds ended with zero growth although each gap-fill reported a solution (the
    integrated solution lost a reaction; a second gap-fill pass on the same medium restores growth, e.g. one
    reaction, rxn00333, for R2_B_D2__metabat2.bin.056). For every protocol stage, if the model is below the stage's
    minimum (biomass 0.01 / 0.02, TMA or CH4 export 0.01) after step 1, the stage's gap-fill is re-run (same medium,
    target, ATP tests and template; at most 2 passes).
 3. Betaine reductase GPR (betaine reducers). The template has no complex for betaine reductase (rxn11678) although
    it carries the roles; the hand-added rxn17220_c0 is given the genes RAST annotates with those roles, labelled
    'annotation-supported, not template-mapped' (no effect on flux).
 4. Methanogen (report only): whether methanogenic growth becomes possible if the template's Fpo reaction rxn24614
    (direction '<') is freed. The model is not changed by this test.
Then the full QC of build_model.qc_model is re-run.

usage: repair_model.py --root R --runs out/A --tag TAG --out final --msdb PATH
"""
import argparse
import gzip
import json
import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_model as B  # noqa: E402

# bound set to 0 = the direction that oxidises N (NH3 -> NH2OH -> nitrite); each reaction keeps its reductive direction
NH3_OX_BLOCK = {"rxn13647": "lower_bound",   # 2 H+ + nitrite + 3 UQH2 <=> 2 H2O + NH3 + 3 UQ
                "rxn03978": "upper_bound",   # 2 H2O + NH3 + 6 cyt c3+ <=> 8 H+ + nitrite + 6 cyt c2+
                "rxn00109": "upper_bound",   # H2O + NAD + NH3 <=> NADH + 2 H+ + NH2OH
                "rxn00570": "lower_bound"}   # H2O + H+ + nitrite <=> O2 + NH2OH
METHANE = ("cpd01024_c", "cpd01024_e")
METHANE_TRANSPORT = "rxn10471"
BETAINE_REDUCTASE_ROLES = ["betaine reductase (EC 1.21.4.4)",
                           "Betaine reductase component B alpha subunit (EC 1.21.4.4)",
                           "Betaine reductase component B beta subunit (EC 1.21.4.4)"]


def stages(guild):
    """(label, medium, target, minimum) checked in order; gap-fill re-run if below minimum."""
    if guild == "mepn_consumer":
        return [("gapfill_1", B.MEDIA_NO_MEPN, "bio1", 0.01), ("gapfill_2", B.M_GENOME_SCALE, "bio1", 0.02)]
    if guild == "betaine_reducer":
        return [("gapfill_1", B.MEDIA_NO_MEPN, "bio1", 0.01), ("gapfill_2", B.BETAINE_MEDIUM, "EX_cpd00441_e0", 0.01)]
    return [("gapfill_1", B.METHYL_MEDIUM, "bio1", 0.01), ("gapfill_2", B.METHYL_MEDIUM, "EX_cpd01024_e0", 0.01)]


def value(model, medium, target):
    with model:
        B.set_medium(model, medium)
        model.reactions.bio1.lower_bound = 0
        model.objective = target
        v = model.slim_optimize(error_value=float("nan"))
    return 0.0 if v != v else v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--msdb", required=True)
    ap.add_argument("--solver", default="gurobi")
    a = ap.parse_args()
    import cobra
    cobra.Configuration().solver = a.solver
    root, runs, out = Path(a.root), Path(a.runs), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    src_rep = json.load(open(runs / f"{a.tag}.report.json"))
    rep = {"tag": a.tag, "mag_id": src_rep["mag_id"], "guild": src_rep["guild"], "arm": src_rep["arm"],
           "source_model_md5": B.md5(runs / f"{a.tag}.json"), "source_report_md5": B.md5(runs / f"{a.tag}.report.json"),
           "code_md5": {"repair_model.py": B.md5(__file__), "build_model.py": B.md5(B.__file__)}, "status": "started"}
    T = B.Timer()
    try:
        guild = src_rep["guild"]
        model = cobra.io.load_json_model(str(runs / f"{a.tag}.json"))
        from modelseedpy.biochem import modelseed_biochem as MB
        MB.ModelSEEDBiochem.default_biochemistry = MB.from_local(a.msdb)
        from modelseedpy.core.msatpcorrection import load_default_medias
        meds = {x[0].id: x[0] for x in load_default_medias()}
        tests = [{"media": meds[t["media"]], "is_max_threshold": True, "threshold": t["threshold"],
                  "objective": "ATPM_c0"} for t in src_rep["atp_correction"]["tests"]]
        free_tests = [{"media": B.to_msmedia({f"EX_{c}_e0": 1000 for c in cpds}, f"freeATP_{lab}"),
                       "is_max_threshold": True, "threshold": 1e-5, "objective": "ATPM_c0"}
                      for lab, cpds in B.FREE_ATP_MEDIA.items()]
        gs_t = B.load_template(root / "templates" / src_rep["template_file"])
        if guild == "mepn_consumer":
            B.mepn_template_additions(gs_t, {})
        with gzip.open(root / "evidence" / f"{src_rep['mag_id']}.json.gz", "rt") as fh:
            ev = json.load(fh)
        rp = {"free_atp_before": B.free_atp(model)}
        # 1. ammonia-oxidation block
        blk = {}
        for rid, bound in NH3_OX_BLOCK.items():
            if rid + "_c" in gs_t.reactions:
                setattr(gs_t.reactions.get_by_id(rid + "_c"), bound, 0)
            if rid + "_c0" in model.reactions:
                r = model.reactions.get_by_id(rid + "_c0")
                blk[rid + "_c0"] = {"bounds_before": list(r.bounds), "gapfilled": src_rep["gapfilled"].get(rid + "_c0"),
                                    "gene_rule": r.gene_reaction_rule}
                setattr(r, bound, 0)
                blk[rid + "_c0"]["bounds_after"] = list(r.bounds)
        rp["nh3_oxidation_block"] = blk
        n_tpl = 0
        for tr in gs_t.reactions:
            if tr.id.startswith(METHANE_TRANSPORT):
                continue
            for met, c in tr.metabolites.items():
                if met.id in METHANE:
                    if c < 0:
                        tr.upper_bound = 0
                    else:
                        tr.lower_bound = 0
                    n_tpl += 1
        mmo = {}
        for r in list(model.reactions):
            if r.id.startswith((METHANE_TRANSPORT, "EX_")):
                continue
            for met, c in r.metabolites.items():
                if met.id in ("cpd01024_c0", "cpd01024_e0") and ((c < 0 and r.upper_bound > 0) or (c > 0 and r.lower_bound < 0)):
                    g = src_rep["gapfilled"].get(r.id)
                    d = {"gapfilled": g, "gene_rule": r.gene_reaction_rule, "bounds_before": list(r.bounds)}
                    if g:
                        if c < 0:
                            r.upper_bound = 0
                        else:
                            r.lower_bound = 0
                        if r.bounds == (0, 0):
                            model.remove_reactions([r])
                            d["action"] = "removed (gap-filled, methane-consuming only)"
                        else:
                            d["action"] = f"methane-consuming direction blocked (gap-filled), bounds {list(r.bounds)}"
                    else:
                        d["action"] = "kept: annotation-supported (draft) reaction that can consume methane"
                    mmo[r.id] = d
                    break
        rp["methane_consumption"] = {"template_reactions_restricted": n_tpl, "model": mmo}
        rp["free_atp_after_block"] = B.free_atp(model)
        # 2. growth restoration, stage by stage
        regf = {}
        gapfilled = {k: v for k, v in src_rep["gapfilled"].items() if k in model.reactions}
        for lab, med, target, mo in stages(guild):
            if guild == "betaine_reducer" and lab == "gapfill_2":
                B.add_missing_media_exchanges(model, B.BETAINE_MEDIUM)
            v0 = value(model, med, target)
            regf[lab] = {"value_before": v0, "minimum": mo, "passes": []}
            for p in (1, 2):
                if v0 >= mo - 1e-9:
                    break
                r = {}
                ad, rv = B.gapfill_stage(model, gs_t, tests + free_tests, med, lab, target, mo,
                                         f"repair_{lab}_pass{p}", r, T)
                gapfilled.update({x: f"repair_{lab}" for x in ad + rv if x not in gapfilled})
                v0 = value(model, med, target)
                regf[lab]["passes"].append({"added": ad, "reversed": rv, "value_after": v0,
                                            "detail": r[f"repair_{lab}_pass{p}"]})
            regf[lab]["value_final"] = v0
        rp["growth_restoration"] = regf
        # 3. betaine reductase GPR from RAST role annotations
        if guild == "betaine_reducer" and "rxn17220_c0" in model.reactions:
            r = model.reactions.rxn17220_c0
            genes = B.role_genes(ev["rast"], BETAINE_REDUCTASE_ROLES)
            rp["rxn17220_gpr_before"] = r.gene_reaction_rule
            if not r.gene_reaction_rule and genes:
                r.gene_reaction_rule = " or ".join(genes)
                r.notes["gpr_source"] = "annotation-supported, not template-mapped (RAST betaine reductase roles)"
            rp["rxn17220_gpr_after"] = r.gene_reaction_rule
            rp["betaine_reductase_rast_genes"] = genes
        # 4. methanogen: Fpo (rxn24614, template '<') freed, report only
        if guild == "methanogen":
            fpo = {}
            for lab, free in (("as_built", False), ("rxn24614_freed", True)):
                with model:
                    if free:
                        if "rxn24614_c0" in model.reactions:
                            model.reactions.rxn24614_c0.bounds = (-1000, 1000)
                        else:
                            model.add_reactions([gs_t.reactions.rxn24614_c.to_reaction(model)])
                            model.reactions.rxn24614_c0.bounds = (-1000, 1000)
                    B.set_medium(model, B.METHYL_MEDIUM)
                    model.objective = "bio1"
                    s = model.optimize()
                    d = {"growth": s.objective_value if s.status == "optimal" else None,
                         "CH4_at_max_growth": float(s.fluxes.get("EX_cpd01024_e0", 0)) if s.status == "optimal" else None,
                         "acetate_at_max_growth": float(s.fluxes.get("EX_cpd00029_e0", 0)) if s.status == "optimal" else None}
                    if "EX_cpd00029_e0" in model.reactions:
                        model.reactions.EX_cpd00029_e0.upper_bound = 0
                        d["growth_acetate_export_blocked"] = model.slim_optimize(error_value=float("nan"))
                    fpo[lab] = d
            rp["methanogen_fpo_test"] = fpo
        B.add_missing_media_exchanges(model, {"BETAINE": B.BETAINE_MEDIUM, "MEPN": B.MEDIA,
                                              "METH": B.METHYL_MEDIUM}["BETAINE" if guild == "betaine_reducer" else
                                                                      "MEPN" if guild == "mepn_consumer" else "METH"])
        if guild == "mepn_consumer":
            B.add_missing_media_exchanges(model, B.M_GENOME_SCALE)
        if guild == "methanogen":
            B.add_missing_media_exchanges(model, B.MANUSCRIPT_METHANOGEN_MEDIUM)
        rp["exchange_bounds_reset"] = B.normalize_exchanges(model)
        with T("qc"):
            qc = B.qc_model(model, guild, tests)
        rep["repair"] = rp
        rep["qc"] = qc
        rep["gapfilled"] = gapfilled
        rep["hand_added"] = src_rep.get("hand_added", [])
        rep["guild_pathway"] = src_rep["guild_pathway"]
        rep["template_file"] = src_rep["template_file"]
        rep["template_rule"] = src_rep["template_rule"]
        rep["final"] = {"n_reactions": len(model.reactions), "n_metabolites": len(model.metabolites),
                        "n_genes": len(model.genes)}
        rep["final"]["n_gene_associated_non_gapfilled"] = sum(
            1 for r in model.reactions if r.gene_reaction_rule and r.id not in gapfilled and r.id not in rep["hand_added"])
        for rid, stage in gapfilled.items():
            if rid in model.reactions:
                model.reactions.get_by_id(rid).notes["gapfilled"] = stage
        notes = json.loads(model.notes.get("salterns_2026", "{}"))
        notes.update({"gapfilled": gapfilled, "final_set_repairs": {
            "nh3_oxidation_block": list(blk), "growth_restoration_added": [x for x, s in gapfilled.items() if s.startswith("repair_")]}})
        model.notes["salterns_2026"] = json.dumps(notes)
        model.objective = "bio1"
        cobra.io.save_json_model(model, str(out / f"{a.tag}.json"))
        rep["status"] = "ok"
    except Exception as e:
        rep["status"] = "error"
        rep["error"] = f"{type(e).__name__}: {e}"
        rep["traceback"] = traceback.format_exc()
    rep["times_s"] = T.t
    with open(out / f"{a.tag}.report.json", "w") as fh:
        json.dump(rep, fh, indent=1, default=str)
    print(a.tag, rep["status"], rep.get("error", ""), flush=True)
    sys.exit(0 if rep["status"] == "ok" else 1)


if __name__ == "__main__":
    main()
