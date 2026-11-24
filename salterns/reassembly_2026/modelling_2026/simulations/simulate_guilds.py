#!/usr/bin/env python
"""Community methane simulations of the 2026 guild models (final MAG set), with the original protocols:

MePn consumers (38 models): analysis(1).ipynb cell 27 at commit 22f4629, as re-implemented in
  clauded_insights/mepn_obligacy_check/analysis/vcommon.py (cell27, add_ratio_constraint, MEDIA):
  model copy -> ratio constraint (MePn uptake >= r x Pi uptake) -> cell-5 medium -> per core: bio1 fixed to the
  core's abundance, EX_cpd00159_e0 (lactate) bounds (-1000, 0) as objective (maximised = least lactate uptake),
  optimise, record CH4 export (EX_cpd01024_e0), lactate, Pi, MePn; r = 0.07 and 0.33.
Betaine reducers (18 models): analysis(1).ipynb cells 37/40 (correlational/betaine_grdH_r232/README.md section 1):
  the cell-5 medium without MePn plus betaine EX_cpd00540_e0 = 0.1 (uptake bound), bio1 fixed, the same lactate
  objective; record TMA export (EX_cpd00441_e0) and its FVA range at the lactate optimum (TMA is often not unique).
Abundance (as in the original input): a(core) = mean over D1, D2 of relative abundance / DNA concentration
  (../../abundance_final/final585_relabund_24sediment.tsv, data/SaltPondsDNA.csv), for the six unrestored cores.
Sensitivity 'mmo_blocked': MePn models carrying the RAST methane monooxygenase reaction rxn00844_c0 rerun with it
  closed (models/README.md, repair 1b).
The methanogen model is not simulated (not methanogenic as built; models/README.md caveat 1).
Reference run ('--reference <models_dir> <out.json>'): every MePn model at a fixed biomass of 1e-3 (both ratios, and
with rxn00844 closed where present), giving each model's CH4 per unit biomass k directly; in this protocol k is a
per-model constant, and at the smallest core abundances (about 1e-4 and below) the fluxes fall under GLPK's tolerance.

Run (TF, ~/tf_runs/salterns_sim_final): ~/venvs/salterns_modelling/bin/python simulate_guilds.py <models_dir> <relabund.tsv> <dna.csv> <out.json>
"""
import json
import sys
import time
from pathlib import Path

import cobra
import pandas as pd
from optlang.symbolics import Zero

cobra.Configuration().solver = "glpk"
CORES = ["R1_A", "R1_B", "R1_C", "R2_A", "R2_B", "R2_C"]
# analysis(1).ipynb cell 5, verbatim (dict order preserved)
MEDIA = {'EX_cpd00067_e0': 1000.0, 'EX_cpd00058_e0': 1000.0, 'EX_cpd00013_e0': 1000.0,
         'EX_cpd00244_e0': 1000.0, 'EX_cpd00205_e0': 1000.0, 'EX_cpd00034_e0': 1000.0,
         'EX_cpd11574_e0': 1000.0, 'EX_cpd00971_e0': 1000.0, 'EX_cpd00048_e0': 1000.0,
         'EX_cpd00030_e0': 1000.0, 'EX_cpd00305_e0': 100.0, 'EX_cpd00001_e0': 1000.0,
         'EX_cpd10516_e0': 1000.0, 'EX_cpd00007_e0': 1000.0, 'EX_cpd00159_e0': 100.0,
         'EX_cpd25960_e0': 1000.0, "EX_cpd00009_e0": 100, 'EX_cpd00063_e0': 1000.0,
         'EX_cpd00149_e0': 1000.0, 'EX_cpd00254_e0': 1000.0, 'EX_cpd00099_e0': 1000.0}
MEPN_EX, PI_EX, CH4_EX, LAC_EX, BET_EX, TMA_EX = ("EX_cpd25960_e0", "EX_cpd00009_e0", "EX_cpd01024_e0",
                                                   "EX_cpd00159_e0", "EX_cpd00540_e0", "EX_cpd00441_e0")
RATIOS = (0.07, 0.33)
MMO = "rxn00844_c0"


def abundance(rel, dna, mag):
    return {c: float((rel.loc[mag, f"{c}_D1"] / dna[f"{c}_D1"] + rel.loc[mag, f"{c}_D2"] / dna[f"{c}_D2"]) / 2) for c in CORES}


def set_medium_strict(model, medium):
    missing = [k for k in medium if k not in model.reactions]
    if missing:
        raise KeyError(f"{model.id}: medium exchanges absent from model: {missing}")
    model.medium = medium


def add_ratio_constraint(model, multiple):
    """analysis(1).ipynb cell 27, verbatim: MePn uptake >= multiple * Pi uptake."""
    constraint = model.problem.Constraint(Zero, lb=0, ub=None, name="Phosphate_ratio")
    model.add_cons_vars(constraint)
    model.solver.update()
    constraint.set_linear_coefficients(coefficients={
        model.reactions.get_by_id(PI_EX).forward_variable: multiple,
        model.reactions.get_by_id(PI_EX).reverse_variable: -multiple,
        model.reactions.get_by_id(MEPN_EX).forward_variable: -1,
        model.reactions.get_by_id(MEPN_EX).reverse_variable: 1})
    model.solver.update()


def flux(sol, rid):
    return float(sol.fluxes[rid]) if sol.status == "optimal" and rid in sol.fluxes.index else float("nan")


def cell27(model_in, ab, multiple, closed=()):
    model = model_in.copy()
    for rid in closed:
        model.reactions.get_by_id(rid).bounds = (0, 0)
    add_ratio_constraint(model, multiple)
    set_medium_strict(model, MEDIA)
    bio = model.reactions.get_by_id("bio1")
    model.solver.configuration.timeout = 60
    out = {}
    for core, a in ab.items():
        bio.bounds = (a, a)
        ex = model.reactions.get_by_id(LAC_EX)
        ex.bounds = (-1000, 0)
        model.objective = ex
        sol = model.optimize()
        out[core] = {"status": sol.status, "abundance": a, "methane": flux(sol, CH4_EX), "lactate": flux(sol, LAC_EX),
                     "phosphate": flux(sol, PI_EX), "MePn": flux(sol, MEPN_EX), "bio1": flux(sol, "bio1"),
                     "mmo": flux(sol, MMO) if MMO in model.reactions else None}
        bio.bounds = (0, 0)
    return out


def betaine(model_in, ab):
    model = model_in.copy()
    med = {k: v for k, v in MEDIA.items() if k != MEPN_EX}
    med[BET_EX] = 0.1
    set_medium_strict(model, med)
    bio = model.reactions.get_by_id("bio1")
    model.solver.configuration.timeout = 60
    out = {}
    for core, a in ab.items():
        bio.bounds = (a, a)
        ex = model.reactions.get_by_id(LAC_EX)
        ex.bounds = (-1000, 0)
        model.objective = ex
        sol = model.optimize()
        rec = {"status": sol.status, "abundance": a, "TMA": flux(sol, TMA_EX), "lactate": flux(sol, LAC_EX),
               "betaine": flux(sol, BET_EX), "bio1": flux(sol, "bio1")}
        if sol.status == "optimal":
            fva = cobra.flux_analysis.flux_variability_analysis(model, [TMA_EX], fraction_of_optimum=1.0, processes=1)
            rec["TMA_fva_min"], rec["TMA_fva_max"] = float(fva.loc[TMA_EX, "minimum"]), float(fva.loc[TMA_EX, "maximum"])
        out[core] = rec
        bio.bounds = (0, 0)
    return out


def main(models_dir, rel_path, dna_path, out_path):
    models_dir = Path(models_dir)
    rel = pd.read_csv(rel_path, sep="\t", index_col=0)
    dna = pd.read_csv(dna_path).set_index("Sample")["DNA conc (µg/g soil)"]
    res = {"protocol": __doc__.split("Run (TF")[0].strip(), "cobra": cobra.__version__, "solver": "glpk",
           "mepn": {}, "betaine": {}, "abundance": {}}
    t0 = time.time()
    for f in sorted((models_dir / "mepn_consumer").glob("*.json")):
        mag = f.stem
        m = cobra.io.load_json_model(str(f))
        ab = abundance(rel, dna, mag)
        res["abundance"][mag] = ab
        res["mepn"][mag] = {str(r): cell27(m, ab, r) for r in RATIOS}
        if MMO in m.reactions:
            res["mepn"][mag]["mmo_blocked"] = {str(r): cell27(m, ab, r, closed=(MMO,)) for r in RATIOS}
        print(f"mepn {mag} {time.time() - t0:.0f}s", flush=True)
    for f in sorted((models_dir / "betaine_reducer").glob("*.json")):
        mag = f.stem
        m = cobra.io.load_json_model(str(f))
        ab = abundance(rel, dna, mag)
        res["abundance"][mag] = ab
        res["betaine"][mag] = betaine(m, ab)
        print(f"betaine {mag} {time.time() - t0:.0f}s", flush=True)
    json.dump(res, open(out_path, "w"), indent=1)
    bad = [(g, m, c) for g in ("mepn", "betaine") for m, d in res[g].items()
           for runs in ([d[r] for r in ("0.07", "0.33")] if g == "mepn" else [d]) for c, x in runs.items() if x["status"] != "optimal"]
    print(f"done in {time.time() - t0:.0f}s; non-optimal cells: {len(bad)} {bad[:10]}")


def reference(models_dir, out_path, a=1e-3):
    out = {}
    for f in sorted((Path(models_dir) / "mepn_consumer").glob("*.json")):
        m = cobra.io.load_json_model(str(f))
        rec = {str(r): cell27(m, {"ref": a}, r)["ref"] for r in RATIOS}
        if MMO in m.reactions:
            rec["mmo_blocked"] = {str(r): cell27(m, {"ref": a}, r, closed=(MMO,))["ref"] for r in RATIOS}
        out[f.stem] = rec
    json.dump({"reference_biomass": a, "solver": "glpk", "mepn": out}, open(out_path, "w"), indent=1)
    bad = [m for m, d in out.items() for r in ("0.07", "0.33") if d[r]["status"] != "optimal"]
    print(f"reference run: {len(out)} models; non-optimal: {bad}")


if __name__ == "__main__":
    if sys.argv[1] == "--reference":
        reference(*sys.argv[2:4])
    else:
        main(*sys.argv[1:5])
