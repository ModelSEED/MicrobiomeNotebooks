"""Shared inputs and the MePn simulation protocol (verbatim logic of correlational/analysis(1).ipynb
cell 27 at 22f4629, cell 5 medium), for the smoke test.  Reads only; writes nothing outside the
scratch directory."""
import json
import os
import re
import time
from pathlib import Path

import pandas as pd

SCRATCH = Path(__file__).resolve().parent
S = Path("/home/freiburger/Documents/MicrobiomeNotebooks/salterns")
CORR = S / "correlational"
RAST_DIR = S / "jplfaria" / "Cliff MAGs" / "RASTified_MAG_json"
BASE_MODELS = CORR / "models"                      # symlink -> cliff_mags/mac/salterns/models
MEPN_MODELS = S / "models" / "MePn_consumers"      # what analysis(1).ipynb cell 7 loaded
AB_AS_RUN = CORR / "betaine_grdH_r232" / "inputs" / "averaged_normalized_abundances_unrestored_22f4629.csv"
STORED = CORR / "saltern_fluxes.json"
TPL = Path("/home/freiburger/Documents/ModelSEEDTemplates/templates")

MAG = "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.55"   # Roseovarius.9 (iterativeIDs.json, guilds_r232.csv)
CORES = ["R1_A", "R1_B", "R1_C", "R2_A", "R2_B", "R2_C"]

# analysis(1).ipynb cell 5, verbatim
MEDIA = {'EX_cpd00067_e0': 1000.0, 'EX_cpd00058_e0': 1000.0, 'EX_cpd00013_e0': 1000.0,
         'EX_cpd00244_e0': 1000.0, 'EX_cpd00205_e0': 1000.0, 'EX_cpd00034_e0': 1000.0,
         'EX_cpd11574_e0': 1000.0, 'EX_cpd00971_e0': 1000.0, 'EX_cpd00048_e0': 1000.0,
         'EX_cpd00030_e0': 1000.0, 'EX_cpd00305_e0': 100.0, 'EX_cpd00001_e0': 1000.0,
         'EX_cpd10516_e0': 1000.0, 'EX_cpd00007_e0': 1000.0, 'EX_cpd00159_e0': 100.0,
         'EX_cpd25960_e0': 1000.0, "EX_cpd00009_e0": 100, 'EX_cpd00063_e0': 1000.0,
         'EX_cpd00149_e0': 1000.0, 'EX_cpd00254_e0': 1000.0, 'EX_cpd00099_e0': 1000.0}

# cliff_mags/mac/salterns/gapfill.ipynb cell 92 (fliu's MePn gapfill medium: no Pi, lactate 1)
M_GENOME_SCALE = {'EX_cpd00159_e0': 1, 'EX_cpd00007_e0': 1000, 'EX_cpd00067_e0': 1000,
                  'EX_cpd00001_e0': 1000, 'EX_cpd25960_e0': 1000, 'EX_cpd00149_e0': 1000,
                  'EX_cpd00099_e0': 1000, 'EX_cpd00063_e0': 1000, 'EX_cpd00058_e0': 1000,
                  'EX_cpd00048_e0': 1000, 'EX_cpd00034_e0': 1000, 'EX_cpd00030_e0': 1000,
                  'EX_cpd00013_e0': 1000, 'EX_cpd00205_e0': 1000, 'EX_cpd00254_e0': 1000,
                  'EX_cpd00971_e0': 1000, 'EX_cpd10515_e0': 1000, 'EX_cpd10516_e0': 1000,
                  'EX_cpd11574_e0': 1000, 'EX_cpd00244_e0': 1000, 'EX_cpd00305_e0': 100.0}


def abundances(mag=MAG):
    ab = pd.read_csv(AB_AS_RUN, index_col=0)
    return ab.loc[mag, CORES]


def add_ratio_constraint(model, multiple):
    """analysis(1).ipynb cell 27: MePn uptake >= multiple * Pi uptake (forward/reverse variables)."""
    from optlang.symbolics import Zero
    c = model.problem.Constraint(Zero, lb=0, ub=None, name="Phosphate_ratio")
    model.add_cons_vars(c)
    model.solver.update()
    c.set_linear_coefficients(coefficients={
        model.reactions.EX_cpd00009_e0.forward_variable: multiple,
        model.reactions.EX_cpd00009_e0.reverse_variable: -multiple,
        model.reactions.EX_cpd25960_e0.forward_variable: -1,
        model.reactions.EX_cpd25960_e0.reverse_variable: 1})
    model.solver.update()
    return c


def prepare_mepn(model, multiple):
    """copy -> ratio constraint -> medium (cell 27 call order); medium keys absent from the model are dropped."""
    m = model.copy()
    add_ratio_constraint(m, multiple)
    m.medium = {k: v for k, v in MEDIA.items() if k in m.reactions}
    m.solver.configuration.timeout = 15
    return m


def set_core(m, a):
    m.reactions.bio1.bounds = (a, a)
    ex = m.reactions.EX_cpd00159_e0
    ex.bounds = (-1000, 0)
    m.objective = ex


def simulate_mepn(model, ab, multiple):
    """Per core: fix bio1 to the abundance, maximise EX_cpd00159_e0 (= minimise lactate uptake)."""
    m = prepare_mepn(model, multiple)
    out = {}
    for core, a in ab.items():
        set_core(m, a)
        sol = m.optimize()
        out[core] = {"status": sol.status, "abundance": float(a),
                     "methane": float(sol.fluxes.EX_cpd01024_e0) if sol.status == "optimal" else None,
                     "lactate": float(sol.fluxes.EX_cpd00159_e0) if sol.status == "optimal" else None,
                     "phosphate": float(sol.fluxes.EX_cpd00009_e0) if sol.status == "optimal" else None,
                     "MePn": float(sol.fluxes.EX_cpd25960_e0) if sol.status == "optimal" else None,
                     "objective": float(sol.objective_value) if sol.status == "optimal" else None}
        m.reactions.bio1.bounds = (0, 0)
    return out


def fva_mepn(model, ab, multiple, rxns=("EX_cpd01024_e0", "EX_cpd25960_e0", "EX_cpd00009_e0")):
    from cobra.flux_analysis import flux_variability_analysis
    m = prepare_mepn(model, multiple)
    out = {}
    for core, a in ab.items():
        set_core(m, a)
        r = flux_variability_analysis(m, list(rxns), fraction_of_optimum=1.0, processes=1)
        out[core] = {k: [float(r.loc[k, "minimum"]), float(r.loc[k, "maximum"])] for k in rxns}
        m.reactions.bio1.bounds = (0, 0)
    return out


def split_functions(fn):
    """fliu's split (cliff_mags/mac/salterns/build_models.ipynb cell 25; same set as
    cobrakbase KBaseGenomeFeature.split_annotation)."""
    return [x for x in re.split("; | / | @", fn) if x]


def msgenome_from_rast(path, extra_features=()):
    """MSGenome whose features carry RAST ontology terms from the KBase Genome JSON 'functions'.
    extra_features: iterable of (feature_id, [functions]) appended as protein-less features."""
    from modelseedpy.core.msgenome import MSGenome, MSFeature
    g = json.load(open(path))
    feats, n_skip = [], 0
    for f in g["features"]:
        fns = [p for fn in f.get("functions", []) for p in split_functions(fn)]
        if not fns:
            n_skip += 1
            continue
        mf = MSFeature(f["id"], f.get("protein_translation", ""))
        for p in fns:
            mf.add_ontology_term("RAST", p)
        feats.append(mf)
    for fid, fns in extra_features:
        mf = MSFeature(fid, "")
        for fn in fns:
            for p in split_functions(fn):
                mf.add_ontology_term("RAST", p)
        feats.append(mf)
    genome = MSGenome()
    genome.id = g["id"]
    genome.add_features(feats)
    return genome, {"n_features_json": len(g["features"]), "n_with_function": len(feats) - len(extra_features),
                    "n_no_function_skipped": n_skip, "n_extra": len(extra_features)}


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
