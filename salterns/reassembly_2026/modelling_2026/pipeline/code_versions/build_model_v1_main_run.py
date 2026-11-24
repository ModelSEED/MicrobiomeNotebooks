"""Build, ATP-correct, gap-fill and QC one genome-scale model: one genome x one guild protocol x one annotation arm.

Same pipeline for every arm of the pre-registered comparison (../ANNOTATION_CHOICE_CRITERIA.md); only the gene ->
reaction evidence differs. The winning arm's outputs are the final models.

Protocol sources (all read, none imported):
  - draft: ModelSEEDpy MSBuilder.build_base_model (311f720) replicated step by step, so that the extra (non-RAST) gene
    evidence of arms B-D can be added as gene associations of TEMPLATE reactions only, before spontaneous/universal
    reactions are chosen. The template itself is never changed by an arm.
  - ATPM: gapfill.ipynb cell 3 (= MSBuilder.add_atpm).
  - ATP correction: MSATPCorrection with Core-V6 and ModelSEEDpy's 54 default ATP media (evaluate, determine, apply,
    expand, build tests), as in correlational/betaine_gapfilling.ipynb cell 0 and the scouting build_pipeline.py; the
    tests are passed to every MSGapfill.
  - MePn consumers (cliff_mags/mac/salterns/gapfill.ipynb): cell 5 disables on the draft; cell 72 template additions
    (cpd25960_e, rxs43120 hand-built C-P lyase, rxs00001 MePn ABC transporter, rxn26447 H+ fix, disable lists) on the
    in-memory gap-filling template; cell 74 rxn05602; cell 76 MePn module; cells 92-93 gap-fill on lactate 1 + MePn
    1000, no Pi, minimum biomass 0.02 -- preceded, per the scouting README ("gap-fill first on a medium without MePn"),
    by a gap-fill on the cell-5 simulation medium without MePn (lactate 100 + Pi 100, aerobic).
  - Betaine reducers (correlational/betaine_gapfilling.ipynb): gap-fill on aerobic lactate (cell-0 medium, without MePn
    per the scouting caution), then add_betaine_reductase_pathway() (cell 2, verbatim logic), then the guild check on
    lactate + betaine 0.1 (the simulation medium); only if the model cannot export TMA (net betaine reduction), a
    targeted gap-fill (target EX_cpd00441_e0 >= 0.01).
  - Methanogen: Archaea template, gap-filled on a methanol + methylamines medium (anaerobic; see METHYL_MEDIUM), then,
    only if it cannot export methane, a targeted gap-fill (target EX_cpd01024_e0 >= 0.01) on the same medium.
QC: guild-medium growth, growth_triplet (clauded_insights/mepn_obligacy_check/analysis/vcommon.py, verbatim logic),
free ATP on three carbon-free media, ATP-test pass rate, blocked / dead-end reactions, gap-filled reaction list,
guild-pathway gene support, weak-evidence reactions.

usage: build_model.py --root R --mag MAG --guild G --arm A --out DIR [--msdb PATH] [--solver gurobi]
  R holds templates/, evidence/<mag>.json.gz, modelling_species.tsv
"""
import argparse
import gzip
import hashlib
import json
import logging
import os
import re
import sys
import time
import traceback
from collections import defaultdict
from pathlib import Path

ARMS = {"A": [], "B": ["kofam"], "C": ["kofam", "eggnog"], "D": ["kofam", "eggnog", "bakta"]}

# ---------------------------------------------------------------- media
# correlational/analysis(1).ipynb cell 5 (simulation medium), verbatim
MEDIA = {'EX_cpd00067_e0': 1000.0, 'EX_cpd00058_e0': 1000.0, 'EX_cpd00013_e0': 1000.0,
         'EX_cpd00244_e0': 1000.0, 'EX_cpd00205_e0': 1000.0, 'EX_cpd00034_e0': 1000.0,
         'EX_cpd11574_e0': 1000.0, 'EX_cpd00971_e0': 1000.0, 'EX_cpd00048_e0': 1000.0,
         'EX_cpd00030_e0': 1000.0, 'EX_cpd00305_e0': 100.0, 'EX_cpd00001_e0': 1000.0,
         'EX_cpd10516_e0': 1000.0, 'EX_cpd00007_e0': 1000.0, 'EX_cpd00159_e0': 100.0,
         'EX_cpd25960_e0': 1000.0, "EX_cpd00009_e0": 100, 'EX_cpd00063_e0': 1000.0,
         'EX_cpd00149_e0': 1000.0, 'EX_cpd00254_e0': 1000.0, 'EX_cpd00099_e0': 1000.0}
# the same without MePn: aerobic lactate 100 + Pi 100 (first gap-fill of MePn consumers and betaine reducers)
MEDIA_NO_MEPN = {k: v for k, v in MEDIA.items() if k != "EX_cpd25960_e0"}
# gapfill.ipynb cell 92 (fliu's MePn gap-fill medium: no Pi, lactate 1)
M_GENOME_SCALE = {'EX_cpd00159_e0': 1, 'EX_cpd00007_e0': 1000, 'EX_cpd00067_e0': 1000,
                  'EX_cpd00001_e0': 1000, 'EX_cpd25960_e0': 1000, 'EX_cpd00149_e0': 1000,
                  'EX_cpd00099_e0': 1000, 'EX_cpd00063_e0': 1000, 'EX_cpd00058_e0': 1000,
                  'EX_cpd00048_e0': 1000, 'EX_cpd00034_e0': 1000, 'EX_cpd00030_e0': 1000,
                  'EX_cpd00013_e0': 1000, 'EX_cpd00205_e0': 1000, 'EX_cpd00254_e0': 1000,
                  'EX_cpd00971_e0': 1000, 'EX_cpd10515_e0': 1000, 'EX_cpd10516_e0': 1000,
                  'EX_cpd11574_e0': 1000, 'EX_cpd00244_e0': 1000, 'EX_cpd00305_e0': 100.0}
# betaine simulation medium (correlational/betaine_grdH_r232/README.md section 1: the lactate medium with MePn
# removed, betaine uptake 0.1)
BETAINE_MEDIUM = {**MEDIA_NO_MEPN, "EX_cpd00540_e0": 0.1}
# methanogen (no original code exists; the manuscript states "dimethylamine + H2"): methanol + mono-/di-/trimethylamine
# as carbon and energy sources, CO2, H2S as sulfur source (no sulfate, no O2, no nitrate: no external electron
# acceptor), NH3, Pi and the cell-5 trace elements and thiamin.
_MINERALS = {k: v for k, v in MEDIA.items() if k not in (
    "EX_cpd00007_e0", "EX_cpd00159_e0", "EX_cpd25960_e0", "EX_cpd00048_e0")}
METHYL_MEDIUM = {**_MINERALS, "EX_cpd00116_e0": 100.0, "EX_cpd00187_e0": 100.0, "EX_cpd00425_e0": 100.0,
                 "EX_cpd00441_e0": 100.0, "EX_cpd00011_e0": 100.0, "EX_cpd00239_e0": 1000.0}
# the manuscript's stated methanogen medium ("dimethylamine + H2"), QC only: DMA + H2 + CO2, H2S, NH3, Pi, trace elements
MANUSCRIPT_METHANOGEN_MEDIUM = {**_MINERALS, "EX_cpd00425_e0": 100.0, "EX_cpd11640_e0": 1000.0, "EX_cpd00011_e0": 100.0,
                                "EX_cpd00239_e0": 1000.0}
# free-ATP media (templates/PROVENANCE.md "Validation"): no carbon, no electron donor
_FA_SALTS = ["cpd00013", "cpd00009", "cpd00048", "cpd00205", "cpd00971", "cpd00254", "cpd00063", "cpd00099"]
FREE_ATP_MEDIA = {
    "water_H": ["cpd00001", "cpd00067"],
    "salts_O2": ["cpd00001", "cpd00067", "cpd00007"] + _FA_SALTS,
    "salts_NO3": ["cpd00001", "cpd00067", "cpd00209"] + _FA_SALTS,
}

# ---------------------------------------------------------------- MePn protocol constants (gapfill.ipynb)
DISABLE = ['rxn23850_c0', 'rxn09193_c0', 'rxn05315_c0', 'rxn05206_c0']                      # cell 5 / 72
DISABLE_LB = ['rxn09188_c0', 'rxn00929_c0', 'rxn00145_c0', 'rxn00499_c0', 'rxn00146_c0',
              'rxn08934_c0', 'rxn08941_c0', 'rxn08935_c0', 'rxn08942_c0', 'rxn11942_c0', 'rxn11941_c0',
              'rxn43329_c0', 'rxn11940_c0', 'rxn40559_c0', 'rxn04159_c0', 'rxn04705_c0', 'rxn14415_c0',
              'sul00003_c0', 'rxn48579_c0', 'rxn39860_c0']                                       # cell 5 / 72
BALANCE_CHECK = ['rxn15395_c0', 'rxn28338_c0']                                                   # cell 5
MEPN_MODULE = ['rxs00001_c', 'rxn26456_c', 'rxn26457_c', 'rxs43120_c', 'rxn26447_c', 'rxn04651_c']  # cell 76
PHN_TRANSPORT_ROLES = {   # gapfill.ipynb cell 63 'phnCDE', grouped by subunit (scouting build_pipeline.py)
    "PhnC": ["Phosphonate ABC transporter ATP-binding protein PhnC (TC 3.A.1.9.1)"],
    "PhnD": ["Phosphonate ABC transporter substrate-binding protein PhnD (TC 3.A.1.9.1)"],
    "PhnE": ["Phosphonate ABC transporter permease protein PhnE (TC 3.A.1.9.1)",
             "Phosphonate ABC transporter permease protein PhnE1 (TC 3.A.1.9.1)",
             "Phosphonate ABC transporter permease protein PhnE2 (TC 3.A.1.9.1)"],
}
PHN_TRANSPORT_KOS = {"PhnC": ["K02041"], "PhnD": ["K02044"], "PhnE": ["K02042"]}

# guild pathways for metric 3 (fraction of reactions gene-associated in the draft, before hand-adding/gap-filling)
GUILD_PATHWAY = {
    # C-P lyase phnGHIJKLM: PhnGHIL complex, PhnM, PhnJ lyase (rxn26459; replaced by the hand-built rxs43120)
    "mepn_consumer": {"rxn26456_c0": "PhnGHIL", "rxn26457_c0": "PhnM", "rxn26459_c0": "PhnJ C-P lyase"},
    # betaine reductase grdH/grdI = EC 1.21.4.4; the template's reaction is rxn11678 (hand-added form rxn17220)
    "betaine_reducer": {"rxn11678_c0": "betaine reductase (EC 1.21.4.4)"},
    # methylotrophic methanogenesis core
    "methanogen": {"rxn03127_c0": "methyl-CoM reductase (mcr)", "rxn03126_c0": "heterodisulfide reductase (hdr)",
                   "rxn03020_c0": "methyl-H4MPT:CoM methyltransferase (mtr)",
                   "rxn24606_c0": "methanol:corrinoid methyltransferase (mtaB)",
                   "rxn24607_c0": "methyl-corrinoid:CoM methyltransferase, methanol (mtaA)",
                   "rxn24608_c0": "monomethylamine methyltransferase (mtmB)",
                   "rxn24609_c0": "methyl-corrinoid:CoM methyltransferase, MMA (mtbA)",
                   "rxn24610_c0": "dimethylamine methyltransferase (mtbB)",
                   "rxn24611_c0": "methyl-corrinoid:CoM methyltransferase, DMA (mtbA)",
                   "rxn24612_c0": "trimethylamine methyltransferase (mttB)",
                   "rxn24613_c0": "methyl-corrinoid:CoM methyltransferase, TMA (mtbA)"},
}
GRAMPOS_PHYLA = re.compile(r"^(Bacillota(_[A-Z]+)?|Actinomycetota|Chloroflexota|Patescibacteriota)$")
GUILD_DIR = {"mepn_consumer": "mepn_consumer", "betaine_reducer": "betaine_reducer", "methanogen": "methanogen"}


def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


def template_for(domain, phylum):
    if domain == "Archaea":
        return "ArchaeaTemplateV6.json", "Archaea -> ArchaeaTemplateV6"
    if GRAMPOS_PHYLA.match(str(phylum)):
        return "GramPosModelTemplateV7.json", f"Bacteria, Gram-positive phylum ({phylum}) -> GramPosModelTemplateV7"
    return "GramNegModelTemplateV7.json", f"Bacteria, other phylum ({phylum}) -> GramNegModelTemplateV7"


class Timer:
    def __init__(self):
        self.t = {}

    def __call__(self, label):
        timer = self

        class _C:
            def __enter__(s):
                s.t0 = time.perf_counter()

            def __exit__(s, *a):
                timer.t[label] = round(time.perf_counter() - s.t0, 2)
        return _C()


def load_template(path):
    from modelseedpy.core.mstemplate import MSTemplateBuilder
    return MSTemplateBuilder.from_dict(json.load(open(path))).build()


def to_msmedia(media_dict, mid):
    from modelseedpy import MSMedia
    m = MSMedia.from_dict({k[3:-3]: v for k, v in media_dict.items()})
    m.id = mid
    return m


def set_medium(model, medium):
    model.medium = {k: v for k, v in medium.items() if k in model.reactions}


def max_obj(model, medium, objective, closed=()):
    with model:
        set_medium(model, medium)
        for rid in closed:
            if rid in model.reactions:
                model.reactions.get_by_id(rid).lower_bound = 0
        model.objective = objective
        v = model.slim_optimize(error_value=float("nan"))
        st = model.solver.status
    return (0.0 if (v == v and abs(v) < 1e-9) else v), st


# ---------------------------------------------------------------- vcommon.py growth_triplet (verbatim logic)
def growth_triplet(model):
    def mg(closed=()):
        with model:
            missing = [k for k in MEDIA if k not in model.reactions]
            if missing:
                raise KeyError(f"{model.id}: medium exchanges absent from model: {missing}")
            model.medium = MEDIA
            for rid in closed:
                model.reactions.get_by_id(rid).lower_bound = 0
            model.objective = "bio1"
            v = model.slim_optimize(error_value=float("nan"))
            st = model.solver.status
        return (0.0 if (v == v and abs(v) < 1e-12) else v), st
    a, sa = mg()
    b, sb = mg(("EX_cpd25960_e0",))
    c, sc = mg(("EX_cpd00009_e0",))
    return {"a_as_is": a, "b_MePn_blocked": b, "c_Pi_blocked_MePn_only_P": c,
            "status": [sa, sb, sc], "obligate": bool(a > 0 and not (b >= 1e-6))}


# ---------------------------------------------------------------- draft with injected evidence
def role_index(template):
    """normalized role name -> template reaction ids where the role is a triggering role of a complex."""
    from modelseedpy.core.msgenome import normalize_role
    idx = defaultdict(set)
    for trxn in template.reactions:
        for cpx in trxn.get_complexes():
            for role, (triggering, optional) in cpx.roles.items():
                if triggering:
                    idx[normalize_role(role.name)].add(trxn.id)
    return idx


def arm_evidence(ev, arm, template):
    """template reaction id -> gene -> set('source:type')"""
    from modelseedpy.core.msgenome import normalize_role
    tmpl = {r.id.rsplit("_", 1)[0]: r.id for r in template.reactions}
    extra = defaultdict(lambda: defaultdict(set))
    for s in ARMS[arm]:
        for g, items in ev["sources"][s].items():
            for typ, term, rxns in items:
                for r in rxns:
                    if r in tmpl:
                        extra[tmpl[r]][g].add(f"{s}:{typ}")
        if s == "bakta":
            ridx = role_index(template)
            for g, prods in ev["bakta_products"].items():
                for p in prods:
                    for trid in ridx.get(normalize_role(p), ()):
                        extra[trid][g].add("bakta:product")
    return extra


def build_draft(genome, template, model_id, extra):
    """MSBuilder.build_base_model(model_id, '0', allow_all_non_grp_reactions=False, annotate_with_rast=False),
    step by step (311f720), with `extra` merged into the reaction -> complex sets as a pseudo-complex 'xev'."""
    from cobra.core import Model
    from modelseedpy.core.msbuilder import MSBuilder
    b = MSBuilder(genome, template, name=model_id)
    b.index = "0"
    cobra_model = Model(model_id)
    b.base_model = cobra_model
    b.generate_reaction_complex_sets()
    rast_rxns = set(b.reaction_to_complex_sets)
    for trid, genes in extra.items():
        have = {g for cpx in b.reaction_to_complex_sets.get(trid, {}).values() for gs in cpx.values() for g in gs}
        new = set(genes) - have
        if new:
            b.reaction_to_complex_sets.setdefault(trid, {})["xev"] = {"xev": new}
    complex_groups = b.build_complex_groups(b.reaction_to_complex_sets.values())
    cobra_model.add_reactions(b.build_metabolic_reactions())
    cobra_model.add_reactions(b.build_non_metabolite_reactions(cobra_model, False))
    cobra_model.add_groups(list(complex_groups.values()))
    b.add_exchanges_to_model(cobra_model)
    biomass_reactions = []
    for rxn_biomass in template.biomasses:
        reaction = rxn_biomass.build_biomass(cobra_model, "0", False, 0.5)
        for m in reaction.metabolites:
            if "modelseed_template_id" in m.notes:
                b.template_species_to_model_species[m.notes["modelseed_template_id"]] = m
        biomass_reactions.append(reaction)
    for rxn in biomass_reactions:
        if rxn.id not in cobra_model.reactions:
            cobra_model.add_reactions([rxn])
    cobra_model.objective = biomass_reactions[0].id
    cobra_model.add_reactions(b.build_drains())
    comp = {}
    for cmp_id, data in b.compartments.items():
        comp[f"{cmp_id}0"] = data.name
        cobra_model.notes[f"kbase_compartment_data_{cmp_id}0"] = {"pH": data.ph, "potential": 0, "compartmentIndex": "0"}
    cobra_model.compartments = comp
    return cobra_model, {r + "0" for r in rast_rxns}


# ---------------------------------------------------------------- MePn protocol pieces
def mepn_template_additions(t, log):
    """gapfill.ipynb cell 72, applied to the in-memory template copy of this build only."""
    from modelseedpy.core.mstemplate import MSTemplateReaction, MSTemplateSpecies
    cc = t.compcompounds
    if "cpd25960_e" not in cc:
        t.add_comp_compounds([MSTemplateSpecies("cpd25960_e", cc.cpd25960_c.charge, "e", "cpd25960")])
    r1 = MSTemplateReaction("rxs43120_c", "rxs43120", "rxs43120", "MePn Degradation", 0, 1000, gapfill_direction=">")
    r1.add_metabolites({cc.cpd00017_c: -1, cc.cpd26519_c: -1, cc.cpd01024_e: 1, cc.cpd03091_c: 1, cc.cpd00060_c: 1,
                        cc.cpd26508_c: 1, cc.cpd00003_c: +1, cc.cpd00004_c: -1})
    r2 = MSTemplateReaction("rxs00001_c", "rxs00001", "MePn ABC Transporter", "MePn Degradation", 0, 1000,
                            gapfill_direction=">")
    r2.add_metabolites({cc.cpd00001_c: -1, cc.cpd00002_c: -1, cc.cpd00008_c: 1, cc.cpd00009_c: 1, cc.cpd00067_c: 1,
                        cc.cpd25960_e: -1, cc.cpd25960_c: 1})
    log["rxs43120_c_balance"] = r1.check_mass_balance()
    log["rxs00001_c_balance"] = r2.check_mass_balance()
    t.add_reactions([r1, r2])
    # cell 72 adds -2 H+ to rxn26447_c (unbalanced in fliu's V6 test template). Applied only if the curated template's
    # rxn26447_c is still H/charge-unbalanced in a way that H+ fixes.
    r = t.reactions.rxn26447_c
    imb = r.check_mass_balance()
    log["rxn26447_c_balance_before"] = dict(imb)
    if imb and set(imb) <= {"H", "charge"} and imb.get("H", 0) == imb.get("charge", 0) != 0:
        r.add_metabolites({cc.cpd00067_c: -imb["H"]})
        log["rxn26447_c_H_fix"] = -imb["H"]
        log["rxn26447_c_balance_after"] = dict(r.check_mass_balance())
    else:
        log["rxn26447_c_H_fix"] = None
    removed = [x[:-1] for x in DISABLE if x[:-1] in t.reactions]
    t.reactions -= removed
    log["template_removed"] = removed
    lb0 = []
    for x in DISABLE_LB:
        if x[:-1] in t.reactions:
            t.reactions.get_by_id(x[:-1]).lower_bound = 0
            lb0.append(x[:-1])
    log["template_lb0"] = lb0


def apply_cell5(model, log):
    off, lb0, hfix = [], [], []
    for rid in DISABLE:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).bounds = (0, 0)
            off.append(rid)
    for rid in DISABLE_LB:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).lower_bound = 0
            lb0.append(rid)
    for rid in BALANCE_CHECK:
        if rid in model.reactions:
            rxn = model.reactions.get_by_id(rid)
            if rxn.check_mass_balance() == {"charge": -1, "H": -1}:
                rxn.add_metabolites({model.metabolites.cpd00067_c0: +1})
                hfix.append(rid)
    log.update(cell5_disabled=off, cell5_lb0=lb0, cell5_H_fixed=hfix)


def role_genes(rast, roles):
    from modelseedpy.core.msgenome import normalize_role
    want = {normalize_role(r) for r in roles}
    return sorted(g for g, rs in rast.items() if any(normalize_role(x) in want for x in rs))


def ko_genes(ev, arm, kos):
    out = set()
    for s in ARMS[arm]:
        for g, items in ev["sources"][s].items():
            if any(t == "KO" and term in kos for t, term, _ in items):
                out.add(g)
    return sorted(out)


def add_mepn_module(model, template, ev, arm, log):
    """cells 74 + 76: rxn05602 and the MePn module from the modified template, as hand-added reactions."""
    from modelseedpy.core.msbuilder import MSBuilder
    added = []
    for trid in ["rxn05602_c"] + MEPN_MODULE:
        if trid + "0" not in model.reactions:
            model.add_reactions([template.reactions.get_by_id(trid).to_reaction(model)])
            added.append(trid + "0")
    MSBuilder.add_exchanges_to_model(model)
    # gene associations for the hand-added reactions (identical rule in every arm)
    groups = {sub: sorted(set(role_genes(ev["rast"], roles)) | set(ko_genes(ev, arm, PHN_TRANSPORT_KOS[sub])))
              for sub, roles in PHN_TRANSPORT_ROLES.items()}
    if "rxs00001_c0" in added:
        parts = ["(" + " or ".join(g) + ")" if len(g) > 1 else g[0] for g in groups.values() if g]
        model.reactions.rxs00001_c0.gene_reaction_rule = " and ".join(parts)
    if "rxs43120_c0" in added and "rxn26459_c0" in model.reactions:
        model.reactions.rxs43120_c0.gene_reaction_rule = model.reactions.rxn26459_c0.gene_reaction_rule
    if "rxn26459_c0" in model.reactions:          # scouting: template lyase made irreversible, as rxs43120 is
        model.reactions.rxn26459_c0.lower_bound = 0
    log.update(mepn_hand_added=added, phn_transport_genes={k: len(v) for k, v in groups.items()},
               rxs00001_gpr=model.reactions.rxs00001_c0.gene_reaction_rule if "rxs00001_c0" in model.reactions else None,
               rxs43120_gpr=model.reactions.rxs43120_c0.gene_reaction_rule if "rxs43120_c0" in model.reactions else None)
    return added


# ---------------------------------------------------------------- betaine protocol (betaine_gapfilling.ipynb cell 2)
def add_betaine_reductase_pathway(model):
    """Verbatim logic of cell 2 (print statements dropped)."""
    from cobra import Metabolite, Reaction
    added = []
    required_mets = {
        'cpd00540_c0': ('Glycine betaine', 'c0', 'C5H11NO2'), 'cpd00441_c0': ('Trimethylamine', 'c0', 'C3H9N'),
        'cpd00196_c0': ('Acetyl phosphate', 'c0', 'C2H3O5P'), 'cpd28060_c0': ('Reduced thioredoxin', 'c0', None),
        'cpd27735_c0': ('Oxidized thioredoxin', 'c0', None), 'cpd00009_c0': ('Phosphate', 'c0', 'HO4P'),
        'cpd00067_c0': ('H+', 'c0', 'H'), 'cpd00001_c0': ('H2O', 'c0', 'H2O'),
        'cpd00033_c0': ('Glycine', 'c0', 'C2H5NO2'), 'cpd00013_c0': ('NH3', 'c0', 'H4N'),
        'cpd00441_e0': ('Trimethylamine', 'e0', 'C3H9N')}
    for met_id, (name, comp, formula) in required_mets.items():
        if met_id in model.metabolites:
            continue
        met = Metabolite(met_id, name=name, compartment=comp)
        if formula:
            met.formula = formula
        model.add_metabolites([met])
    M = model.metabolites
    if 'rxn17220_c0' not in model.reactions:
        rxn = Reaction('rxn17220_c0')
        rxn.name, rxn.lower_bound, rxn.upper_bound = 'Betaine reductase', 0, 1000
        rxn.add_metabolites({M.cpd00540_c0: -1, M.cpd11421_c0: -1, M.cpd00009_c0: -1, M.cpd00067_c0: 0,
                             M.cpd00441_c0: 1, M.cpd00196_c0: 1, M.cpd11420_c0: 1, M.cpd00001_c0: 1})
        imbalance = rxn.check_mass_balance()
        if imbalance:
            h_needed = -imbalance.get('H', 0)
            if h_needed != 0:
                rxn.add_metabolites({M.cpd00067_c0: h_needed})
        model.add_reactions([rxn])
        added.append('rxn17220_c0')
    if 'rxn07207_c0' not in model.reactions:
        rxn = Reaction('rxn07207_c0')
        rxn.name, rxn.lower_bound, rxn.upper_bound = 'Glycine reductase', 0, 1000
        rxn.add_metabolites({M.cpd00033_c0: -1, M.cpd11421_c0: -1, M.cpd00009_c0: -1, M.cpd00067_c0: -1,
                             M.cpd00196_c0: 1, M.cpd00013_c0: 1, M.cpd11420_c0: 1, M.cpd00001_c0: 1})
        model.add_reactions([rxn])
        added.append('rxn07207_c0')
    if 'EX_cpd00441_e0' not in model.reactions:
        ex = Reaction('EX_cpd00441_e0')
        ex.name, ex.lower_bound, ex.upper_bound = 'TMA exchange', -1000, 1000
        ex.add_metabolites({M.cpd00441_e0: -1})
        model.add_reactions([ex])
        added.append('EX_cpd00441_e0')
    if not any(r for r in model.reactions if 'cpd00441_c0' in [m.id for m in r.metabolites]
               and 'cpd00441_e0' in [m.id for m in r.metabolites]):
        tp = Reaction('rxn_TMA_transport_c0')
        tp.name, tp.lower_bound, tp.upper_bound = 'TMA transport', -1000, 1000
        tp.add_metabolites({M.cpd00441_c0: -1, M.cpd00441_e0: 1})
        model.add_reactions([tp])
        added.append('rxn_TMA_transport_c0')
    return added


def ensure_thioredoxin(model, template):
    """cell 2 requires cpd11420_c0 / cpd11421_c0 (it would raise otherwise); add them from the template if absent."""
    out = []
    for cid in ("cpd11420_c", "cpd11421_c"):
        if cid + "0" not in model.metabolites:
            model.add_metabolites([template.compcompounds.get_by_id(cid).to_metabolite("0")])
            out.append(cid + "0")
    return out


# ---------------------------------------------------------------- generic pieces
def add_missing_media_exchanges(model, media_dict):
    """betaine_gapfilling.ipynb cell 0: an exchange for every medium compound the model lacks."""
    from cobra import Metabolite, Reaction
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
    """scouting build_pipeline.py: gap-filling leaves the media package's bounds on the exchanges; reset every EX_ to
    (-1000, 1000) as in the original models (model.medium then closes every uptake not in a medium)."""
    n = 0
    for r in model.reactions:
        if r.id.startswith("EX_") and r.bounds != (-1000, 1000):
            r.bounds = (-1000, 1000)
            n += 1
    return n


def is_boundary(r):
    return r.id.startswith(("EX_", "SK_", "DM_")) or r.id in ("bio1", "ATPM_c0") or len(r.metabolites) == 1


def gapfill_stage(model, template, tests, medium, mid, target, min_obj, label, rep, T):
    from modelseedpy import MSGapfill
    before = {r.id: r.bounds for r in model.reactions}
    with T(label):
        gf = MSGapfill(model, default_gapfill_templates=[template], test_conditions=tests, default_target="bio1")
        sol = gf.run_gapfilling(media=to_msmedia(medium, mid), target=target, minimum_obj=min_obj,
                                binary_check=False, prefilter=True)
        integ = None
        fail_status = None if sol else str(gf.gfmodel.solver.status)
        if sol:
            integ = gf.integrate_gapfill_solution(sol, cumulative_solution=[], remove_unneeded_reactions=True)
    after = {r.id: r.bounds for r in model.reactions}
    added = sorted(r for r in after if r not in before and not r.startswith("EX_"))
    widened = sorted(r for r in after if r in before and not r.startswith("EX_") and
                     (after[r][0] < before[r][0] or after[r][1] > before[r][1]))
    rep[label] = {"medium": mid, "target": target, "minimum_obj": min_obj,
                  "solution_found": bool(sol), "gfmodel_solver_status_if_failed": fail_status,
                  "solution": None if not sol else {"new": dict(sol.get("new", {})), "reversed": dict(sol.get("reversed", {}))},
                  "integrated": None if not integ else {"new": dict(integ.get("new", {})),
                                                        "reversed": dict(integ.get("reversed", {})),
                                                        "growth": integ.get("growth")},
                  "added_reactions": {r: [after[r][0], after[r][1]] for r in added},
                  "reversed_reactions": {r: {"before": list(before[r]), "after": list(after[r])} for r in widened}}
    return added, widened


def atp_test_pass(model, tests):
    from modelseedpy.core.msmodelutl import MSModelUtil
    m = model.copy()
    u = MSModelUtil.get(m)
    res = []
    for t in tests:
        with m:
            ok = u.test_single_condition(t, apply_condition=True, model=m)
            res.append({"media": t["media"].id, "threshold": t["threshold"], "pass": bool(ok)})
    return res


def free_atp(model):
    out = {}
    for lab, cpds in FREE_ATP_MEDIA.items():
        med = {f"EX_{c}_e0": 1000 for c in cpds}
        v, st = max_obj(model, med, "ATPM_c0")
        out[lab] = v
    return out


def dead_ends(model):
    internal = [m for m in model.metabolites if not any(is_boundary(r) and r.id.startswith(("EX_", "SK_", "DM_"))
                                                         for r in m.reactions)]
    dead = []
    for m in internal:
        prod = cons = False
        for r in m.reactions:
            c = r.metabolites[m]
            if (c > 0 and r.upper_bound > 0) or (c < 0 and r.lower_bound < 0):
                prod = True
            if (c < 0 and r.upper_bound > 0) or (c > 0 and r.lower_bound < 0):
                cons = True
        if not (prod and cons):
            dead.append(m.id)
    dead_set = set(dead)
    dead_rxns = sorted({r.id for mid in dead for r in model.metabolites.get_by_id(mid).reactions if not is_boundary(r)})
    return {"n_internal_metabolites": len(internal), "dead_end_metabolites": len(dead),
            "dead_end_metabolite_frac": len(dead) / max(1, len(internal)), "dead_end_reactions": len(dead_rxns),
            "dead_end_metabolite_ids": sorted(dead_set)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--mag", required=True)
    ap.add_argument("--guild", required=True, choices=list(GUILD_PATHWAY))
    ap.add_argument("--arm", required=True, choices=list(ARMS))
    ap.add_argument("--out", required=True)
    ap.add_argument("--msdb", default=os.environ.get("MODELSEED_DB_PATH", str(Path.home() / "Documents/ModelSEEDDatabase")))
    ap.add_argument("--solver", default="gurobi")
    a = ap.parse_args()
    root, out = Path(a.root), Path(a.out)
    tag = f"{a.mag}__{a.guild}__{a.arm}"
    out.mkdir(parents=True, exist_ok=True)
    rep_path = out / f"{tag}.report.json"
    logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(name)s %(levelname)s %(message)s")

    import cobra
    import pandas as pd
    import modelseedpy
    cobra.Configuration().solver = a.solver
    T = Timer()
    rep = {"tag": tag, "mag_id": a.mag, "guild": a.guild, "arm": a.arm, "arm_sources": ARMS[a.arm],
           "status": "started", "solver": a.solver, "versions": {"cobra": cobra.__version__,
           "modelseedpy": modelseedpy.__version__, "modelseedpy_file": modelseedpy.__file__, "python": sys.version}}
    try:
        sp = pd.read_csv(root / "modelling_species.tsv", sep="\t")
        sp["mag_id"] = sp.rep_genome.str.split("|").str[1]
        row = sp.set_index("mag_id").loc[a.mag]
        tfile, rule = template_for(row.domain, row.phylum)
        if a.guild == "methanogen":
            assert tfile.startswith("Archaea")
        rep.update(domain=row.domain, phylum=row.phylum, family=row.family, genus=row.genus,
                   template_file=tfile, template_rule=rule, template_md5=md5(root / "templates" / tfile),
                   core_template_md5=md5(root / "templates/Core-V6.json"))
        with gzip.open(root / "evidence" / f"{a.mag}.json.gz", "rt") as fh:
            ev = json.load(fh)
        from modelseedpy.biochem import modelseed_biochem as MB
        with T("load_biochem"):
            MB.ModelSEEDBiochem.default_biochemistry = MB.from_local(a.msdb)
        rep["msdb"] = a.msdb
        with T("load_templates"):
            gs_t = load_template(root / "templates" / tfile)
            core_t = load_template(root / "templates/Core-V6.json")
        # genome: Prodigal genes with RAST roles
        from modelseedpy.core.msgenome import MSGenome, MSFeature
        feats = []
        for g, roles in ev["rast"].items():
            f = MSFeature(g, "")
            for r in roles:
                f.add_ontology_term("RAST", r)
            feats.append(f)
        genome = MSGenome()
        genome.id = a.mag
        genome.add_features(feats)
        extra = arm_evidence(ev, a.arm, gs_t)
        with T("draft"):
            model, rast_rxns = build_draft(genome, gs_t, a.mag, extra)
        model.id = a.mag
        draft_ids = {r.id for r in model.reactions}
        ext_rxns = {k + "0" for k in extra}
        evidence_by_rxn = {}
        for r in model.reactions:
            if r.id in rast_rxns or r.id in ext_rxns:
                ex = extra.get(r.id[:-1], {})
                evidence_by_rxn[r.id] = {"rast": r.id in rast_rxns,
                                         "extra": {g: sorted(v) for g, v in ex.items()},
                                         "template_type": str(gs_t.reactions.get_by_id(r.id[:-1]).type)}
        gene_assoc_draft = {r.id for r in model.reactions if r.gene_reaction_rule}
        rep["draft"] = {"n_reactions": len(model.reactions), "n_metabolites": len(model.metabolites),
                        "n_genes": len(model.genes), "n_gene_associated": len(gene_assoc_draft),
                        "n_rast_supported": len(rast_rxns & draft_ids),
                        "n_extra_supported": len(ext_rxns & draft_ids),
                        "n_extra_only": len((ext_rxns - rast_rxns) & draft_ids),
                        "n_template_rxns_with_extra_evidence": len(extra),
                        "n_bakta_extra_genes_used": len({g for v in extra.values() for g in v if g.startswith("bakta_")})}
        # guild pathway support (metric 3), from the draft
        pw = GUILD_PATHWAY[a.guild]
        pws = {rid: bool(rid in model.reactions and model.reactions.get_by_id(rid).gene_reaction_rule) for rid in pw}
        rep["guild_pathway"] = {"reactions": pw, "gene_supported": pws,
                                "fraction": sum(pws.values()) / len(pws),
                                "evidence": {rid: evidence_by_rxn.get(rid) for rid in pw}}
        # ATPM (cell 3)
        from modelseedpy.core.msbuilder import MSBuilder
        MSBuilder.add_atpm(model)
        prot = {}
        hand_added = []
        if a.guild == "mepn_consumer":
            apply_cell5(model, prot)
            mepn_template_additions(gs_t, prot)
            hand_added += add_mepn_module(model, gs_t, ev, a.arm, prot)
        rep["protocol"] = prot
        rep["draft_free_atp_before_atp_correction"] = free_atp(model)
        # ATP correction
        from modelseedpy.core.msatpcorrection import MSATPCorrection, load_default_medias
        n_before_atp = {r.id: r.bounds for r in model.reactions}
        with T("atp_correction"):
            atp = MSATPCorrection(model, core_t, load_default_medias(), compartment="c0",
                                  atp_hydrolysis_id="ATPM_c0", load_default_medias=False)
            atp.evaluate_growth_media()
            atp.determine_growth_media()
            atp.apply_growth_media_gapfilling()
            atp.expand_model_to_genome_scale()
            tests = atp.build_tests()
        after_atp = {r.id: r.bounds for r in model.reactions}
        rep["atp_correction"] = {
            "selected_media": [m.id for m in (atp.selected_media or [])],
            "n_tests": len(tests), "tests": [{"media": t["media"].id, "threshold": t["threshold"]} for t in tests],
            "core_gapfilling": [list(x) for x in (atp.cumulative_core_gapfilling or [])],
            "expansion_filtered": [list(x)[:2] for x in (atp.filtered_noncore or [])],
            "reactions_added": sorted(set(after_atp) - set(n_before_atp)),
            "reactions_removed": sorted(set(n_before_atp) - set(after_atp)),
            "free_atp_after": free_atp(model)}
        gapfilled = {}
        if a.guild == "mepn_consumer":
            for lab, med, mid, mo in (("gapfill_1", MEDIA_NO_MEPN, "lactate_Pi_noMePn", None),
                                      ("gapfill_2", M_GENOME_SCALE, "fliu_PiFree_lactate1_MePn", 0.02)):
                ad, rv = gapfill_stage(model, gs_t, tests, med, mid, "bio1", mo, lab, rep, T)
                gapfilled.update({r: lab for r in ad + rv if r not in gapfilled})
            sim_media = [MEDIA, M_GENOME_SCALE]
        elif a.guild == "betaine_reducer":
            ad, rv = gapfill_stage(model, gs_t, tests, MEDIA_NO_MEPN, "lactate_Pi_aerobic", "bio1", None,
                                   "gapfill_1", rep, T)
            gapfilled.update({r: "gapfill_1" for r in ad + rv})
            prot["thioredoxin_added"] = ensure_thioredoxin(model, gs_t)
            bet_added = add_betaine_reductase_pathway(model)
            hand_added += [r for r in bet_added if not r.startswith("EX_")]
            MSBuilder.add_exchanges_to_model(model)
            # gene associations for the hand-added reductases: the arm's evidence for the template forms
            # (betaine reductase rxn11678, glycine reductase rxn11676), which are written in the reverse direction
            for hid, tid in (("rxn17220_c0", "rxn11678_c0"), ("rxn07207_c0", "rxn11676_c0")):
                if hid in bet_added and tid in model.reactions:
                    model.reactions.get_by_id(hid).gene_reaction_rule = model.reactions.get_by_id(tid).gene_reaction_rule
            prot["betaine_hand_added"] = bet_added
            prot["rxn17220_gpr"] = model.reactions.rxn17220_c0.gene_reaction_rule
            add_missing_media_exchanges(model, BETAINE_MEDIUM)
            # betaine-reduction capacity = maximum TMA export (net betaine -> TMA). rxn17220_c0 flux alone is not used:
            # with KO evidence the template's rxn11678 (the same reaction written TMA -> betaine, irreversible) enters
            # the draft and forms a null cycle with rxn17220_c0 (test build R2_sediment__metabat2.bin.435, arm B).
            with model:
                set_medium(model, BETAINE_MEDIUM)
                model.reactions.bio1.lower_bound = 0
                model.objective = "EX_cpd00441_e0"
                cap0 = model.slim_optimize(error_value=float("nan"))
            prot["TMA_export_capacity_before_gf2"] = cap0
            if not (cap0 == cap0 and cap0 >= 0.01):
                ad, rv = gapfill_stage(model, gs_t, tests, BETAINE_MEDIUM, "lactate_betaine", "EX_cpd00441_e0", 0.01,
                                       "gapfill_2", rep, T)
                gapfilled.update({r: "gapfill_2" for r in ad + rv if r not in gapfilled})
            sim_media = [BETAINE_MEDIUM, MEDIA_NO_MEPN]
        else:
            ad, rv = gapfill_stage(model, gs_t, tests, METHYL_MEDIUM, "methanol_methylamines", "bio1", None,
                                   "gapfill_1", rep, T)
            gapfilled.update({r: "gapfill_1" for r in ad + rv})
            # the bio1 gap-fill does not require methanogenesis (test build: the model grew as a methyl-to-acetate
            # acetogen, CH4 = 0); stage 2 requires methane export >= 0.01 on the same medium, biomass unconstrained
            add_missing_media_exchanges(model, {"EX_cpd01024_e0": 0})
            with model:
                set_medium(model, METHYL_MEDIUM)
                model.objective = "EX_cpd01024_e0"
                ch4_0 = model.slim_optimize(error_value=float("nan"))
            prot["methane_capacity_before_gf2"] = ch4_0
            if not (ch4_0 == ch4_0 and ch4_0 >= 0.01):
                ad, rv = gapfill_stage(model, gs_t, tests, METHYL_MEDIUM, "methanol_methylamines", "EX_cpd01024_e0",
                                       0.01, "gapfill_2", rep, T)
                gapfilled.update({r: "gapfill_2" for r in ad + rv if r not in gapfilled})
            sim_media = [METHYL_MEDIUM, MANUSCRIPT_METHANOGEN_MEDIUM]
        rep["hand_added"] = hand_added
        for med in sim_media:
            rep.setdefault("added_media_exchanges", []).extend(add_missing_media_exchanges(model, med))
        rep["exchange_bounds_reset"] = normalize_exchanges(model)
        rep["gapfilled"] = gapfilled
        # ------------------------------------------------ QC
        with T("qc"):
            qc = {}
            if a.guild == "mepn_consumer":
                qc["growth_sim_medium"] = max_obj(model, MEDIA, "bio1")
                qc["growth_PiFree_MePn"] = max_obj(model, M_GENOME_SCALE, "bio1")
                qc["growth_no_MePn"] = max_obj(model, MEDIA_NO_MEPN, "bio1")
                qc["growth_triplet"] = growth_triplet(model)
                guild_growth = qc["growth_PiFree_MePn"][0]
                with model:   # methane at max growth on the Pi-free MePn medium
                    set_medium(model, M_GENOME_SCALE)
                    model.objective = "bio1"
                    s = model.optimize()
                    qc["CH4_at_max_growth_PiFree"] = float(s.fluxes.get("EX_cpd01024_e0", 0.0)) if s.status == "optimal" else None
                    qc["MePn_at_max_growth_PiFree"] = float(s.fluxes.get("EX_cpd25960_e0", 0.0)) if s.status == "optimal" else None
            elif a.guild == "betaine_reducer":
                qc["growth_betaine_medium"] = max_obj(model, BETAINE_MEDIUM, "bio1")
                qc["growth_lactate_aerobic"] = max_obj(model, MEDIA_NO_MEPN, "bio1")
                guild_growth = qc["growth_betaine_medium"][0]
                with model:
                    set_medium(model, BETAINE_MEDIUM)
                    model.reactions.bio1.lower_bound = 0.01 if guild_growth >= 0.01 else 0
                    model.objective = "rxn17220_c0"
                    qc["betaine_reductase_max_at_bio_0.01"] = model.slim_optimize(error_value=float("nan"))
                    model.objective = "EX_cpd00441_e0"
                    qc["TMA_export_max_at_bio_0.01"] = model.slim_optimize(error_value=float("nan"))
                qc["growth_no_O2"] = max_obj(model, {k: v for k, v in BETAINE_MEDIUM.items() if k != "EX_cpd00007_e0"}, "bio1")
            else:
                qc["growth_methyl_medium"] = max_obj(model, METHYL_MEDIUM, "bio1")
                guild_growth = qc["growth_methyl_medium"][0]
                from cobra.flux_analysis import flux_variability_analysis
                with model:
                    set_medium(model, METHYL_MEDIUM)
                    model.objective = "bio1"
                    s = model.optimize()
                    qc["CH4_at_max_growth"] = float(s.fluxes.get("EX_cpd01024_e0", 0.0)) if s.status == "optimal" else None
                    if s.status == "optimal" and "EX_cpd01024_e0" in model.reactions:
                        f = flux_variability_analysis(model, ["EX_cpd01024_e0"], fraction_of_optimum=0.99, processes=1)
                        qc["CH4_FVA_at_99pct_growth"] = [float(f.loc["EX_cpd01024_e0", "minimum"]),
                                                         float(f.loc["EX_cpd01024_e0", "maximum"])]
                    model.objective = "EX_cpd01024_e0"
                    qc["CH4_max_biomass_free"] = model.slim_optimize(error_value=float("nan"))
                    model.objective = "bio1"
                    if "EX_cpd00029_e0" in model.reactions:
                        with model:
                            model.reactions.EX_cpd00029_e0.upper_bound = 0
                            qc["growth_acetate_export_blocked"] = model.slim_optimize(error_value=float("nan"))
                    for sub in ("EX_cpd00116_e0", "EX_cpd00187_e0", "EX_cpd00425_e0", "EX_cpd00441_e0"):
                        qc[f"growth_only_{sub}"] = max_obj(model, {k: v for k, v in METHYL_MEDIUM.items()
                                                                    if k not in {"EX_cpd00116_e0", "EX_cpd00187_e0",
                                                                                 "EX_cpd00425_e0", "EX_cpd00441_e0"} - {sub}},
                                                            "bio1")[0]
                qc["growth_manuscript_DMA_H2_medium"] = max_obj(model, MANUSCRIPT_METHANOGEN_MEDIUM, "bio1")
                with model:
                    set_medium(model, MANUSCRIPT_METHANOGEN_MEDIUM)
                    model.objective = "EX_cpd01024_e0"
                    qc["CH4_max_manuscript_medium"] = model.slim_optimize(error_value=float("nan"))
            qc["guild_growth"] = guild_growth
            qc["grows"] = bool(guild_growth == guild_growth and guild_growth >= 1e-6)
            qc["free_atp"] = free_atp(model)
            qc["free_atp_pass"] = all(v == v and v <= 1e-6 for v in qc["free_atp"].values())
            at = atp_test_pass(model, tests)
            qc["atp_tests"] = at
            qc["atp_test_pass_rate"] = sum(x["pass"] for x in at) / max(1, len(at))
            from cobra.flux_analysis import find_blocked_reactions
            internal = [r.id for r in model.reactions if not is_boundary(r)]
            blocked = find_blocked_reactions(model, open_exchanges=True, processes=1)
            bi = sorted(set(blocked) & set(internal))
            qc["n_internal_reactions"] = len(internal)
            qc["blocked_internal"] = len(bi)
            qc["blocked_frac"] = len(bi) / max(1, len(internal))
            qc["blocked_ids"] = bi
            qc["dead_ends"] = dead_ends(model)
        rep["qc"] = qc
        # ------------------------------------------------ metrics
        final_ids = {r.id for r in model.reactions}
        m1 = sorted(r for r in (gene_assoc_draft & final_ids)
                    if model.reactions.get_by_id(r).gene_reaction_rule and r not in gapfilled)
        gf_count = len([r for r in gapfilled if r not in hand_added])
        weak = []
        for rid, evd in evidence_by_rxn.items():
            if rid not in final_ids or evd["rast"]:
                continue
            tags = {t for v in evd["extra"].values() for t in v}
            tools = {t.split(":")[0] for t in tags}
            if len(tools) == 1:
                tool = next(iter(tools))
                types = {t.split(":")[1] for t in tags}
                if types <= {"EC"} or (tool == "eggnog" and "KO" not in types):
                    weak.append(rid)
        rep["metrics"] = {
            "m1_gene_associated_non_gapfilled": len(m1),
            "m2_gapfilled": gf_count,
            "m2_atp_core_gapfill": len(rep["atp_correction"]["core_gapfilling"]),
            "m2_gapfilled_by_stage": {s: sum(1 for v in gapfilled.values() if v == s) for s in set(gapfilled.values())},
            "m3_guild_pathway_fraction": rep["guild_pathway"]["fraction"],
            "m4_free_atp_pass": qc["free_atp_pass"], "m4_atp_test_pass_rate": qc["atp_test_pass_rate"],
            "m5_blocked_frac": qc["blocked_frac"], "m5_dead_end_metabolite_frac": qc["dead_ends"]["dead_end_metabolite_frac"],
            "m6_weak_only_reactions": len(weak), "m6_weak_ids": sorted(weak),
            "extra_only_reactions_in_final": len([r for r, e in evidence_by_rxn.items() if r in final_ids and not e["rast"]]),
            "extra_only_by_template_type": {t: sum(1 for r, e in evidence_by_rxn.items() if r in final_ids and not e["rast"]
                                                   and e["template_type"] == t)
                                            for t in sorted({e["template_type"] for e in evidence_by_rxn.values()})},
            "grows": qc["grows"]}
        rep["final"] = {"n_reactions": len(model.reactions), "n_metabolites": len(model.metabolites),
                        "n_genes": len(model.genes), "n_internal_reactions": qc["n_internal_reactions"]}
        rep["evidence_by_reaction"] = {r: e for r, e in evidence_by_rxn.items() if r in final_ids}
        model.notes["salterns_2026"] = json.dumps({"guild": a.guild, "arm": a.arm, "template": tfile,
                                                   "gapfilled": gapfilled, "hand_added": hand_added})
        for rid, stage in gapfilled.items():
            model.reactions.get_by_id(rid).notes["gapfilled"] = stage
        for rid in hand_added:
            if rid in model.reactions:
                model.reactions.get_by_id(rid).notes["hand_added"] = a.guild + " protocol"
        model.objective = "bio1"
        cobra.io.save_json_model(model, str(out / f"{tag}.json"))
        rep["status"] = "ok"
    except Exception as e:  # record and exit non-zero
        rep["status"] = "error"
        rep["error"] = f"{type(e).__name__}: {e}"
        rep["traceback"] = traceback.format_exc()
    rep["times_s"] = T.t
    json.dump(rep, open(rep_path, "w"), indent=1, default=str)
    print(tag, rep["status"], rep.get("error", ""), json.dumps(rep.get("metrics", {}), default=str)[:400], flush=True)
    sys.exit(0 if rep["status"] == "ok" else 1)


if __name__ == "__main__":
    main()
