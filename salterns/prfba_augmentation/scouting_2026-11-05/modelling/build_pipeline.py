"""One build/ATP-correction/gapfill pipeline, applied identically to every arm (original RAST annotation,
or an augmented function list).  Steps follow correlational/betaine_gapfilling.ipynb cell 0 (ATP correction
+ MSGapfill + missing exchanges) on a draft built directly from the genome with fliu's builder settings
(cliff_mags/mac/salterns/build_models.ipynb cell 9: build(id, annotate_with_rast=False)).
The MePn module (hand-added transporter, irreversible C-P lyase) is the same for every arm."""
import json
import re
import time

import cobra
from cobra import Metabolite, Reaction

import common as C

PHN_TRANSPORT_ROLES = {   # gapfill.ipynb cell 63 'phnCDE', grouped by subunit
    "PhnC": ["Phosphonate ABC transporter ATP-binding protein PhnC (TC 3.A.1.9.1)"],
    "PhnD": ["Phosphonate ABC transporter substrate-binding protein PhnD (TC 3.A.1.9.1)"],
    "PhnE": ["Phosphonate ABC transporter permease protein PhnE (TC 3.A.1.9.1)",
             "Phosphonate ABC transporter permease protein PhnE1 (TC 3.A.1.9.1)",
             "Phosphonate ABC transporter permease protein PhnE2 (TC 3.A.1.9.1)"],
}
MEPN_RXNS = {"rxn26456_c0": "PhnGHIL synthase", "rxn26457_c0": "PhnM diphosphatase",
             "rxn26459_c0": "PhnJ C-P lyase (template; releases CH4 in c0)", "rxn26447_c0": "PhnP",
             "rxn04651_c0": "PhnN", "rxn10471_c0": "methane transport (spontaneous)",
             "EX_cpd01024_e0": "methane exchange", "rxs00001_c0": "MePn ABC transporter (hand-added)",
             "EX_cpd25960_e0": "MePn exchange", "rxs43120_c0": "fliu C-P lyase (NADH, CH4 to e0)"}


def load_templates(version):
    from modelseedpy.core.mstemplate import MSTemplateBuilder
    gs, core = {"v7": ("v7.0/GramNegModelTemplateV7.json", "v7.0/Core-V6.json"),
                "v6": ("v6.0/GramNegModelTemplateV6.json", "v6.0/Core-V5.2.json")}[version]
    t0 = time.perf_counter()
    # built FRESH per model: modelseedpy ties gapfill state to the template object
    # (codiffusion_bioreactor/build_prfba_community_models.py _load_template docstring)
    gs_t = MSTemplateBuilder.from_dict(json.load(open(C.TPL / gs))).build()
    core_t = MSTemplateBuilder.from_dict(json.load(open(C.TPL / core))).build()
    return gs_t, core_t, {"gs_template": str(C.TPL / gs), "core_template": str(C.TPL / core),
                          "t_load_s": round(time.perf_counter() - t0, 2)}


def role_genes(genome, roles):
    from modelseedpy.core.msgenome import normalize_role
    want = {normalize_role(r) for r in roles}
    return sorted(f.id for f in genome.features
                  if any(normalize_role(t) in want for t in f.ontology_terms.get("RAST", [])))


def add_mepn_module(model, genome):
    """Identical in every arm: (1) the hand-added MePn ABC transporter (fliu's rxs00001, gapfill.ipynb
    cell 47) with a GPR from the annotated PhnC/PhnD/PhnE genes; (2) the template C-P lyase rxn26459
    (direction '=' in the template) made irreversible, as fliu's replacement rxs43120 was."""
    log = {}
    if "cpd25960_c0" not in model.metabolites:
        log["transporter"] = "skipped: no cpd25960_c0 in model (no MePn pathway reaction in draft)"
    else:
        mc = model.metabolites.cpd25960_c0
        if "cpd25960_e0" not in model.metabolites:
            model.add_metabolites([Metabolite("cpd25960_e0", mc.formula, "MePn [e0]", mc.charge, "e0")])
        m = model.metabolites
        r = Reaction("rxs00001_c0", "MePn ABC Transporter", "MePn Degradation", 0, 1000)
        r.add_metabolites({m.cpd00001_c0: -1, m.cpd00002_c0: -1, m.cpd00008_c0: 1, m.cpd00009_c0: 1,
                           m.cpd00067_c0: 1, m.cpd25960_e0: -1, mc: 1})
        groups = {sub: role_genes(genome, roles) for sub, roles in PHN_TRANSPORT_ROLES.items()}
        parts = ["(" + " or ".join(g) + ")" if len(g) > 1 else g[0] for g in groups.values() if g]
        r.gene_reaction_rule = " and ".join(parts)
        model.add_reactions([r])
        from modelseedpy.core.msbuilder import MSBuilder
        MSBuilder.add_exchanges_to_model(model)
        log["transporter"] = {"GPR_subunits": {k: len(v) for k, v in groups.items()},
                              "mass_balance": r.check_mass_balance()}
    if "rxn26459_c0" in model.reactions:
        model.reactions.rxn26459_c0.lower_bound = 0
        log["rxn26459_c0_bounds"] = model.reactions.rxn26459_c0.bounds
    return log


def add_missing_media_exchanges(model, media_dict):
    """betaine_gapfilling.ipynb cell 0: an exchange for every medium compound the model lacks."""
    added = []
    for ex_id in media_dict:
        if ex_id in model.reactions:
            continue
        cpd = ex_id[3:]
        if cpd not in model.metabolites:
            model.add_metabolites([Metabolite(cpd, compartment="e0")])
        ex = Reaction(ex_id, lower_bound=-1000, upper_bound=1000)
        ex.add_metabolites({model.metabolites.get_by_id(cpd): -1})
        model.add_reactions([ex])
        added.append(ex_id)
    return added


def normalize_exchanges(model):
    """Gapfilling leaves the media package's bounds on the model's exchanges (here ub = 100 on 57 EX_).
    Reset every EX_ to (-1000, 1000), as in the original models; the simulation's model.medium then closes
    every uptake not in the medium.  Applied identically to every arm."""
    changed = [r.id for r in model.reactions if r.id.startswith("EX_") and r.bounds != (-1000, 1000)]
    for r in model.reactions:
        if r.id.startswith("EX_"):
            r.bounds = (-1000, 1000)
    return len(changed)


def to_msmedia(media_dict):
    from modelseedpy import MSMedia
    return MSMedia.from_dict({k[3:-3]: v for k, v in media_dict.items()})


def build_and_gapfill(genome, model_id, version="v7", gapfill_media=None, reaction_scores=None,
                      second_gapfill_media=None, second_min_obj=0.02):
    """Returns (draft, gapfilled model, report).  gapfill_media: dict of EX_ ids (default: cell-5 medium)."""
    from modelseedpy.core.msbuilder import MSBuilder
    from modelseedpy.core.msatpcorrection import MSATPCorrection, load_default_medias
    from modelseedpy import MSGapfill
    gapfill_media = gapfill_media or C.MEDIA
    rep, T = {"model_id": model_id, "template_version": version}, C.Timer()
    # modelseedpy/config.cfg has biochem path "../ModelSEEDDatabase" relative to the CWD; MSATPCorrection's
    # expansion test needs it (msmodelutl.assign_reliability_scores_to_reactions -> ModelSEEDBiochem.get())
    from modelseedpy.biochem import modelseed_biochem as MB
    if not MB.ModelSEEDBiochem.default_biochemistry:
        with T("load_biochem"):
            MB.ModelSEEDBiochem.default_biochemistry = MB.from_local("/home/freiburger/Documents/ModelSEEDDatabase")
    gs_t, core_t, tinfo = load_templates(version)
    rep.update(tinfo)
    with T("build_base_model"):
        model = MSBuilder(genome, gs_t, name=model_id).build_base_model(model_id, "0", False, False)
    rep["draft_n_reactions_before_atpm"] = len(model.reactions)
    MSBuilder.add_atpm(model)
    rep["mepn_module"] = add_mepn_module(model, genome)
    draft = model.copy()
    rep["draft_n_reactions"] = len(model.reactions)
    rep["draft_n_genes"] = len(model.genes)
    with T("atp_correction"):
        atp = MSATPCorrection(model, core_t, load_default_medias(), compartment="c0",
                              atp_hydrolysis_id="ATPM_c0", load_default_medias=False)
        atp.evaluate_growth_media()
        atp.determine_growth_media()
        atp.apply_growth_media_gapfilling()
        atp.expand_model_to_genome_scale()
        tests = atp.build_tests()
    rep["atp_tests"] = [{"media": t["media"].id, "threshold": t["threshold"]} for t in tests]
    rep["atp_core_gapfilling_n"] = len(getattr(atp, "cumulative_core_gapfilling", []) or [])
    rep["n_reactions_after_atp_correction"] = len(model.reactions)
    sols = []
    for label, md, mo in (("gapfill_1", gapfill_media, None), ("gapfill_2", second_gapfill_media, second_min_obj)):
        if md is None:
            continue
        with T(label):
            gf = MSGapfill(model, default_gapfill_templates=[gs_t], test_conditions=tests,
                           default_target="bio1", reaction_scores=reaction_scores or {})
            sol = gf.run_gapfilling(media=to_msmedia(md), target="bio1", minimum_obj=mo,
                                    binary_check=False, prefilter=True)
            if sol:
                gf.integrate_gapfill_solution(sol, cumulative_solution=[])
        rep[label] = {"medium": "cell-5 simulation medium" if md is C.MEDIA else
                      ("fliu Pi-free M_GENOME_SCALE" if md is C.M_GENOME_SCALE else "custom"),
                      "solution": None if not sol else {"new": dict(sol.get("new", {})),
                                                        "reversed": dict(sol.get("reversed", {}))}}
        sols.append(sol)
    rep["added_media_exchanges"] = add_missing_media_exchanges(model, C.MEDIA)
    rep["exchange_bounds_reset"] = normalize_exchanges(model)
    rep["times_s"] = T.t
    rep["final_n_reactions"] = len(model.reactions)
    rep["final_n_genes"] = len(model.genes)
    return draft, model, rep


def growth_report(model):
    out = {}
    for lab, md in (("cell5_medium", C.MEDIA), ("fliu_PiFree_lactate1", C.M_GENOME_SCALE)):
        m = model.copy()
        m.medium = {k: v for k, v in md.items() if k in m.reactions}
        m.objective = "bio1"
        out[f"max_bio1_{lab}"] = m.slim_optimize(error_value=float("nan"))
        try:
            s = cobra.flux_analysis.pfba(m)
            out[f"pfba_bio1_{lab}"] = float(s.fluxes["bio1"])
            out[f"pfba_CH4_{lab}"] = float(s.fluxes.get("EX_cpd01024_e0", 0.0))
            out[f"pfba_MePn_{lab}"] = float(s.fluxes.get("EX_cpd25960_e0", 0.0))
        except Exception as e:  # infeasible
            out[f"pfba_{lab}"] = f"error: {str(e)[:80]}"
    return out


def mepn_presence(model):
    return {rid: (rid in model.reactions and {
        "bounds": model.reactions.get_by_id(rid).bounds,
        "GPR": re.sub(r"Salt_Pond_\S+?RAST\.", "", model.reactions.get_by_id(rid).gene_reaction_rule)[:120]})
        for rid in MEPN_RXNS}
