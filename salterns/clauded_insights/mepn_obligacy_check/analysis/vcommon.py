"""Shared inputs + the analysis(1).ipynb cell-27 protocol (commit 22f4629), re-implemented here so that nothing
is imported from (or written into) the salterns repository.  READ-ONLY on every repo path.

Cell 5 medium, cell 7 model paths, cell 24 growth check and cell 27 loop are reproduced verbatim (see
../analysis_1_22f4629.ipynb, extracted with `git show 22f4629:"salterns/correlational/analysis(1).ipynb"`).
"""
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True

import cobra  # noqa: E402
import pandas as pd  # noqa: E402

cobra.Configuration().solver = "glpk"

OUT = Path(__file__).resolve().parent
S = Path("/home/freiburger/Documents/MicrobiomeNotebooks/salterns")
CORR = S / "correlational"
MEPN_MODELS = S / "models" / "MePn_consumers"                    # analysis(1).ipynb cell 7
DRAFTS = S / "cliff_mags" / "mac" / "salterns" / "models"       # fliu's pre-gapfill drafts  <id>__.RAST.json
RAST_DIR = S / "jplfaria" / "Cliff MAGs" / "RASTified_MAG_json"
AB_AS_RUN = CORR / "betaine_grdH_r232" / "inputs" / "averaged_normalized_abundances_unrestored_22f4629.csv"
STORED = CORR / "saltern_fluxes.json"
IDS = json.load(open(CORR / "iterativeIDs_r232.json"))
TPL_DIR = Path("/home/freiburger/Documents/ModelSEEDTemplates/templates")
MSDB = Path("/home/freiburger/Documents/ModelSEEDDatabase")

# analysis(1).ipynb cell 5, verbatim (dict order preserved)
MEDIA = {'EX_cpd00067_e0': 1000.0, 'EX_cpd00058_e0': 1000.0, 'EX_cpd00013_e0': 1000.0,
         'EX_cpd00244_e0': 1000.0, 'EX_cpd00205_e0': 1000.0, 'EX_cpd00034_e0': 1000.0,
         'EX_cpd11574_e0': 1000.0, 'EX_cpd00971_e0': 1000.0, 'EX_cpd00048_e0': 1000.0,
         'EX_cpd00030_e0': 1000.0, 'EX_cpd00305_e0': 100.0, 'EX_cpd00001_e0': 1000.0,
         'EX_cpd10516_e0': 1000.0, 'EX_cpd00007_e0': 1000.0, 'EX_cpd00159_e0': 100.0,
         'EX_cpd25960_e0': 1000.0, "EX_cpd00009_e0": 100, 'EX_cpd00063_e0': 1000.0,
         'EX_cpd00149_e0': 1000.0, 'EX_cpd00254_e0': 1000.0, 'EX_cpd00099_e0': 1000.0}
MEPN_EX, PI_EX, CH4_EX, LAC_EX = "EX_cpd25960_e0", "EX_cpd00009_e0", "EX_cpd01024_e0", "EX_cpd00159_e0"
RATIOS = (0.07, 0.33)


def model_ids():
    return sorted(p.name[:-5] for p in MEPN_MODELS.glob("*.json"))


def mag_of(model_id):
    return model_id.replace(".contigs", "")


def iterative_id(model_id):
    return IDS[mag_of(model_id)]


def load(model_id):
    return cobra.io.load_json_model(str(MEPN_MODELS / f"{model_id}.json"))


def load_draft(model_id):
    return cobra.io.load_json_model(str(DRAFTS / f"{model_id}__.RAST.json"))


def rast_path(model_id):
    return RAST_DIR / f"{model_id}__.RAST.json"


def abundances_as_run():
    """rows = MAG (no .contigs), columns in the CSV's (= the notebook's row.items()) order."""
    return pd.read_csv(AB_AS_RUN, index_col=0)


def stored():
    return json.load(open(STORED))


def set_medium_strict(model, medium=MEDIA):
    missing = [k for k in medium if k not in model.reactions]
    if missing:
        raise KeyError(f"{model.id}: medium exchanges absent from model: {missing}")
    model.medium = medium


def add_ratio_constraint(model, multiple):
    """analysis(1).ipynb cell 27, verbatim: MePn uptake >= multiple * Pi uptake."""
    from optlang.symbolics import Zero
    constraint = model.problem.Constraint(Zero, lb=0, ub=None, name="Phosphate_ratio")
    model.add_cons_vars(constraint)
    model.solver.update()
    constraint.set_linear_coefficients(coefficients={
        model.reactions.EX_cpd00009_e0.forward_variable: multiple,
        model.reactions.EX_cpd00009_e0.reverse_variable: -multiple,
        model.reactions.EX_cpd25960_e0.forward_variable: -1,
        model.reactions.EX_cpd25960_e0.reverse_variable: 1})
    model.solver.update()
    return constraint


def cell27(model_in, row, multiple, extra=None):
    """One model, one ratio: the cell-27 loop body, verbatim order (copy -> constraint -> medium -> per core:
    bio fixed, EX lactate (-1000, 0) objective, optimize, record (even if infeasible), bio (0,0))."""
    model = model_in.copy()
    add_ratio_constraint(model, multiple)
    set_medium_strict(model)
    bio = model.reactions.get_by_id("bio1")
    model.solver.configuration.timeout = 15
    out = {}
    for sample, abundance in row.items():
        bio.bounds = (abundance, abundance)
        EX = model.reactions.get_by_id("EX_cpd00159_e0")
        EX.bounds = (-1000, 0)
        model.objective = EX
        sol = model.optimize()
        rec = {"status": sol.status, "abundance": float(abundance),
               "methane": float(sol.fluxes.EX_cpd01024_e0), "lactate": float(sol.fluxes.EX_cpd00159_e0),
               "phosphate": float(sol.fluxes.EX_cpd00009_e0), "MePn": float(sol.fluxes.EX_cpd25960_e0),
               "bio1": float(sol.fluxes.bio1)}
        if extra is not None:
            rec.update(extra(model, sol))
        out[sample] = rec
        bio.bounds = (0, 0)
    return out


def max_growth(model, medium=MEDIA, closed=(), objective="bio1"):
    with model:
        set_medium_strict(model, medium)
        for rid in closed:
            model.reactions.get_by_id(rid).lower_bound = 0
        model.objective = objective
        v = model.slim_optimize(error_value=float("nan"))
        st = model.solver.status
    return (0.0 if (v == v and abs(v) < 1e-12) else v), st


def growth_triplet(model):
    a, sa = max_growth(model)
    b, sb = max_growth(model, closed=(MEPN_EX,))
    c, sc = max_growth(model, closed=(PI_EX,))
    return {"a_as_is": a, "b_MePn_blocked": b, "c_Pi_blocked_MePn_only_P": c,
            "status": [sa, sb, sc], "obligate": bool(a > 0 and not (b >= 1e-6))}
