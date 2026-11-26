#!/usr/bin/env python
"""
Re-run of the betaine -> TMA simulations behind ``correlational/saltern_fluxes_betaine.json``,
with a diagnosis of every model and a re-simulation restricted to the 9 *grdH* carriers
(GTDB R232 gene-based guilds, ``correlational/r232/guilds_r232.csv``).

Provenance of the original numbers (read from the repo, not assumed)
--------------------------------------------------------------------
* Producer: ``correlational/analysis(1).ipynb`` at commit 22f4629 (2026-04-09), cells 37/40;
  the same code is in ``analysis(2).ipynb`` cells 33/36 and ``analysis_2_marimo.py``
  (the marimo copy has a renamed-variable bug in ``_add_atpm`` and would not run as is).
* Models: ``correlational/models/<MAG>.contigs__.RAST_gapfilled_betaine.json``, built in
  ``betaine_gapfilling.ipynb``: RAST draft -> MSATPCorrection + MSGapfill on the aerobic
  lactate medium (``*_gapfilled.json``) -> ``add_betaine_reductase_pathway()`` (cell 2), which
  hand-adds rxn17220_c0 (betaine reductase, empty GPR), rxn07207_c0, EX_cpd00441_e0 and
  rxn_TMA_transport_c0 to all 15 models.
* Settings: ATPM_c0 added if absent (lb 0); medium = the fixed lactate/MePn medium without
  MePn, plus EX_cpd00540_e0 = 0.1 (a betaine UPTAKE bound; TMA export is uncapped);
  bio1 fixed to the core's abundance; objective = EX_cpd00159_e0 with bounds (-1000, 0),
  maximised (cobra default) = minimise lactate uptake at fixed biomass.
* Abundances: the ``averaged_normalized_abundances_unrestored.csv`` that existed at 22f4629
  (depth MEAN of DNA-normalised relative abundance), saved here as
  ``inputs/averaged_normalized_abundances_unrestored_22f4629.csv``.  The current
  ``data/averaged_normalized_abundances_unrestored.csv`` (DNA-weighted depth merge, written
  after 2026-04-22) does NOT reproduce the logged abundances in ``sols_log_betaine.txt``.

What this script does
---------------------
1. Re-simulates all 15 models with GLPK and (if licensed) Gurobi on the as-run abundances
   and compares with saltern_fluxes_betaine.json (objective, TMA, betaine).
2. Flux variability of TMA export at the lactate optimum (fraction_of_optimum = 1) for all
   15 x 6 cores: is the stored TMA unique, saturated at the 0.1 uptake bound, or one
   arbitrary point of an alternative-optimum range?
3. Diagnoses each model (provenance and GPR of the betaine reactions, betaine transporters,
   annotated grd genes in the RAST genome and whether any reaches a model GPR, marginal
   lactate cost of betaine reduction, O2 dependence).
4. Re-simulates the 9 grdH carriers with the same settings (GLPK, as-run abundances), and
   as a labelled sensitivity with the current abundance file.

No file under correlational/models/ or any existing flux file is modified.
Run:  ~/Documents/py_venv/bin/python rerun_betaine_grdH_r232.py
Outputs (this folder): fluxes_grdH_r232.json, reproduction_check_r232.csv,
diagnostics_grdH_r232.csv, fva_tma_by_core_r232.csv
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import cobra
from cobra.core import Reaction
from cobra.flux_analysis import flux_variability_analysis

HERE = Path(__file__).resolve().parent
CORR = HERE.parent
REPO = CORR.parent
MODELS = CORR / "models"
RAST = REPO / "jplfaria" / "Cliff MAGs" / "RASTified_MAG_json"
AB_AS_RUN = HERE / "inputs" / "averaged_normalized_abundances_unrestored_22f4629.csv"
AB_CURRENT = REPO / "data" / "averaged_normalized_abundances_unrestored.csv"
EXISTING = CORR / "saltern_fluxes_betaine.json"
GUILDS = CORR / "r232" / "guilds_r232.csv"
PHYLO_2020 = REPO / "data" / "Saltern_phylogeny.json"

# verbatim from analysis(1).ipynb / analysis_2_marimo.py (MePn exchange popped before the betaine run)
MEDIA = {'EX_cpd00067_e0': 1000.0, 'EX_cpd00058_e0': 1000.0, 'EX_cpd00013_e0': 1000.0,
         'EX_cpd00244_e0': 1000.0, 'EX_cpd00205_e0': 1000.0, 'EX_cpd00034_e0': 1000.0,
         'EX_cpd11574_e0': 1000.0, 'EX_cpd00971_e0': 1000.0, 'EX_cpd00048_e0': 1000.0,
         'EX_cpd00030_e0': 1000.0, 'EX_cpd00305_e0': 100.0, 'EX_cpd00001_e0': 1000.0,
         'EX_cpd10516_e0': 1000.0, 'EX_cpd00007_e0': 1000.0, 'EX_cpd00159_e0': 100.0,
         'EX_cpd00009_e0': 100, 'EX_cpd00063_e0': 1000.0, 'EX_cpd00149_e0': 1000.0,
         'EX_cpd00254_e0': 1000.0, 'EX_cpd00099_e0': 1000.0}
BETAINE_UPTAKE = 0.1
TOL = 1e-9          # reproduction tolerance
SAT_TOL = 1e-6      # "at the 0.1 bound" / "unique" tolerance


def add_atpm(model):
    """verbatim _add_atpm from analysis(1).ipynb cell 37"""
    if 'ATPM_c0' not in model.reactions:
        atpm = Reaction('ATPM_c0', 'ATPM', 'ATPM', 0, 1000)
        m = model.metabolites
        atpm.add_metabolites({m.cpd00001_c0: -1, m.cpd00002_c0: -1, m.cpd00008_c0: 1,
                              m.cpd00009_c0: 1, m.cpd00067_c0: 1})
        model.add_reactions([atpm])
    return model


def load_prepared(mag_contigs: str, solver: str, return_growth: bool = False):
    """Same call order as analysis(1).ipynb cell 40: copy -> slim_optimize() on the stored
    medium/objective (bio1) -> set medium.  The slim_optimize() matters only for which
    alternative optimum the simplex returns (warm-start basis), not for any optimum value."""
    cobra.Configuration().solver = solver
    model = add_atpm(cobra.io.load_json_model(MODELS / f"{mag_contigs}__.RAST_gapfilled_betaine.json"))
    model = model.copy()
    growth = model.slim_optimize()
    model.medium = {k: v for k, v in {**MEDIA, 'EX_cpd00540_e0': BETAINE_UPTAKE}.items() if k in model.reactions}
    model.solver.configuration.timeout = 15
    return (model, growth) if return_growth else model


def set_core(model, abundance):
    model.reactions.bio1.bounds = (abundance, abundance)
    ex = model.reactions.EX_cpd00159_e0
    ex.bounds = (-1000, 0)
    model.objective = ex


def simulate(mags, abund: pd.DataFrame, solver: str) -> dict:
    out = {k: {} for k in ("TMA", "betaine", "lactate", "phosphate", "status")}
    out["growth_on_stored_medium"] = {}
    for mag in mags:
        model, growth = load_prepared(mag, solver, return_growth=True)
        out["growth_on_stored_medium"][mag] = float(growth)
        for k in ("TMA", "betaine", "lactate", "phosphate", "status"):
            out[k][mag] = {}
        for core, a in abund.loc[mag.replace(".contigs", "")].items():
            set_core(model, a)
            sol = model.optimize()
            out["status"][mag][core] = sol.status
            out["TMA"][mag][core] = float(sol.fluxes.EX_cpd00441_e0)
            out["betaine"][mag][core] = float(sol.fluxes.EX_cpd00540_e0)
            out["lactate"][mag][core] = float(sol.fluxes.EX_cpd00159_e0)
            out["phosphate"][mag][core] = float(sol.fluxes.EX_cpd00009_e0)
            model.reactions.bio1.bounds = (0, 0)
        print(f"  [{solver}] {mag.split('bins')[1]} done", flush=True)
    return out


def fva_tma(mags, abund: pd.DataFrame, solver: str = "glpk") -> dict:
    lo, hi = {}, {}
    for mag in mags:
        model = load_prepared(mag, solver)
        lo[mag], hi[mag] = {}, {}
        for core, a in abund.loc[mag.replace(".contigs", "")].items():
            set_core(model, a)
            r = flux_variability_analysis(model, ["EX_cpd00441_e0"], fraction_of_optimum=1.0, processes=1)
            lo[mag][core] = float(r.loc["EX_cpd00441_e0", "minimum"])
            hi[mag][core] = float(r.loc["EX_cpd00441_e0", "maximum"])
            model.reactions.bio1.bounds = (0, 0)
        print(f"  [FVA {solver}] {mag.split('bins')[1]} done", flush=True)
    return {"TMA_min": lo, "TMA_max": hi}


def classify(lo, hi, v):
    if hi < SAT_TOL:
        return "zero_forced"
    if hi - lo < SAT_TOL:
        return "saturated_at_0.1_uptake_bound" if abs(hi - BETAINE_UPTAKE) < SAT_TOL else "unique_below_bound"
    return "degenerate_alternative_optima"


def marginal_cost(mag, abund, solver="glpk"):
    """lactate uptake (mmol/gDW/h, positive) at TMA fixed 0 and 0.05 for each core;
    cost per unit TMA = (uptake@0.05 - uptake@0)/0.05; negative = betaine reduction saves lactate."""
    model = load_prepared(mag, solver)
    res = {}
    for core, a in abund.loc[mag.replace(".contigs", "")].items():
        set_core(model, a)
        vals = []
        for t in (0.0, 0.05):
            with model:
                model.reactions.EX_cpd00441_e0.bounds = (t, t)
                vals.append(model.slim_optimize(error_value=np.nan))
        res[core] = (-vals[1] - (-vals[0])) / 0.05 if np.all(np.isfinite(vals)) else np.nan
        model.reactions.bio1.bounds = (0, 0)
    return res


def block_rev_glycine_reductase(mag, abund, solver="glpk", core="R1_A"):
    """sensitivity: close rxn11676_c0 (draft glycine reductase, glycine-synthesis direction) and
    report the FVA max of TMA at the new lactate optimum in one core (NOT applied as a fix)."""
    model = load_prepared(mag, solver)
    if "rxn11676_c0" not in model.reactions:
        return "rxn11676 absent"
    model.reactions.rxn11676_c0.bounds = (0, 0)
    set_core(model, abund.loc[mag.replace(".contigs", ""), core])
    if not np.isfinite(model.slim_optimize(error_value=np.nan)):
        return "infeasible"
    r = flux_variability_analysis(model, ["EX_cpd00441_e0"], fraction_of_optimum=1.0, processes=1)
    return f"TMA[{r.loc['EX_cpd00441_e0', 'minimum']:.4f},{r.loc['EX_cpd00441_e0', 'maximum']:.4f}]"


def anaerobic_feasible(mag, abund, solver="glpk"):
    model = load_prepared(mag, solver)
    model.reactions.EX_cpd00007_e0.lower_bound = 0
    ok = {}
    for core, a in abund.loc[mag.replace(".contigs", "")].items():
        set_core(model, a)
        ok[core] = bool(np.isfinite(model.slim_optimize(error_value=np.nan)))
        model.reactions.bio1.bounds = (0, 0)
    return ok


GRD_PATTERNS = {
    "betaine_reductase_compB_EC1.21.4.4": re.compile(r"Betaine reductase component B", re.I),
    "glycine_reductase_compB_EC1.21.4.2": re.compile(r"Glycine reductase component B", re.I),
    "grd_protein_A": re.compile(r"Glycine/sarcosine/betaine reductase protein A", re.I),
    "grd_component_C": re.compile(r"Glycine/sarcosine/betaine reductase component C", re.I),
}


def rxn_table(path):
    return {r["id"]: r for r in json.load(open(path))["reactions"]}


def diagnose_structure(mag):
    base = rxn_table(MODELS / f"{mag}__.RAST.json")
    gf = rxn_table(MODELS / f"{mag}__.RAST_gapfilled.json")
    gb = rxn_table(MODELS / f"{mag}__.RAST_gapfilled_betaine.json")
    gb_only = rxn_table(MODELS / f"{mag}__.RAST_betaine.json")

    def prov(rid):
        if rid not in gb:
            return "absent"
        if rid in base:
            return "draft (RAST)"
        if rid in gf:
            return "MSGapfill"
        return "hand-added: add_betaine_reductase_pathway()"

    def gpr(rid):
        return gb[rid].get("gene_reaction_rule", "") if rid in gb else ""

    gene2rxn = {}
    for r in gb.values():
        for gid in re.findall(r"[\w.]+CDS\.\d+", r.get("gene_reaction_rule", "")):
            gene2rxn.setdefault(gid, set()).add(r["id"])
    rast = json.load(open(RAST / f"{mag}__.RAST.json"))
    genes = {k: [] for k in GRD_PATTERNS}
    for f in rast["features"]:
        fn = " / ".join(f.get("functions", []))
        for k, p in GRD_PATTERNS.items():
            if p.search(fn):
                genes[k].append(f["id"])
    short = lambda ids: ";".join("CDS." + i.split("CDS.")[-1] for i in ids)
    bet_genes = genes["betaine_reductase_compB_EC1.21.4.4"]
    bet_in_model = sorted({r for g in bet_genes for r in gene2rxn.get(g, ())})
    bet_transport = [rid for rid, r in gb.items()
                     if "cpd00540_e0" in r["metabolites"] and "cpd00540_c0" in r["metabolites"]]
    h_coef = gb["rxn17220_c0"]["metabolites"].get("cpd00067_c0", 0.0) if "rxn17220_c0" in gb else np.nan
    return {
        "betaine_json_same_rxns_as_gapfilled_betaine": set(gb_only) == set(gb),
        "rxn17220_provenance": prov("rxn17220_c0"),
        "rxn17220_GPR": gpr("rxn17220_c0") or "(empty)",
        "rxn17220_H+_coef": h_coef,
        "rxn07207_provenance": prov("rxn07207_c0"),
        "TMA_transport_provenance": prov("rxn_TMA_transport_c0"),
        "EX_TMA_provenance": prov("EX_cpd00441_e0"),
        "rxn11676_rev_glycine_reductase": prov("rxn11676_c0") + (" (GPR set)" if gpr("rxn11676_c0") else ""),
        "betaine_transporters": ";".join(f"{rid}[{prov(rid)}{', GPR' if gpr(rid) else ', no GPR'}]" for rid in bet_transport),
        "n_betaine_reductase_compB_genes": len(bet_genes),
        "betaine_reductase_compB_genes": short(bet_genes),
        "betaine_reductase_compB_genes_in_any_model_GPR": ";".join(bet_in_model) or "none",
        "n_glycine_reductase_compB_genes": len(genes["glycine_reductase_compB_EC1.21.4.2"]),
        "n_grd_protein_A_genes": len(genes["grd_protein_A"]),
        "n_grd_component_C_genes": len(genes["grd_component_C"]),
    }


def main():
    existing = json.load(open(EXISTING))
    mags15 = list(existing["TMA"])
    guilds = pd.read_csv(GUILDS).set_index("mag")
    grdh9 = [m for m in mags15 if bool(guilds.loc[m.replace(".contigs", ""), "betaine_reducer"])]
    assert len(grdh9) == 9, grdh9
    ab_run = pd.read_csv(AB_AS_RUN, index_col=0)
    ab_cur = pd.read_csv(AB_CURRENT, index_col=0)
    cores = list(next(iter(existing["TMA"].values())).keys())
    ab_run, ab_cur = ab_run[cores], ab_cur[cores]

    # ------------------------------------------------------------- 1. reproduction
    print("1. reproduction of saltern_fluxes_betaine.json (15 models, as-run abundances)")
    resim = {"glpk": simulate(mags15, ab_run, "glpk")}
    try:
        import gurobipy
        gurobipy.Model()          # raises without a licence
        resim["gurobi"] = simulate(mags15, ab_run, "gurobi")
        fva_solver = "gurobi"     # tighter tolerances than GLPK for the fraction_of_optimum=1 FVA
    except Exception as e:  # licence or package missing
        print("  gurobi unavailable:", e)
        fva_solver = "glpk"

    # ------------------------------------------------------------- 2. FVA envelope
    print(f"2. FVA of TMA export at the lactate optimum (15 models x 6 cores, {fva_solver})")
    fva_run = fva_tma(mags15, ab_run, fva_solver)

    rows = []
    for solver, sim in resim.items():
        for m in mags15:
            d = lambda k: max(abs(sim[k][m][c] - existing[k][m][c]) for c in cores)
            inside = all(fva_run["TMA_min"][m][c] - SAT_TOL <= existing["TMA"][m][c] <= fva_run["TMA_max"][m][c] + SAT_TOL
                         for c in cores)
            rows.append(dict(solver=solver, mag=m.replace(".contigs", ""),
                             iterativeID_r232=guilds.loc[m.replace(".contigs", ""), "iterativeID_r232"],
                             grdH=m in grdh9,
                             max_abs_diff_lactate_objective=d("lactate"), max_abs_diff_TMA=d("TMA"),
                             max_abs_diff_betaine=d("betaine"),
                             TMA_reproduced=d("TMA") < TOL, objective_reproduced=d("lactate") < 1e-8,
                             existing_TMA_within_FVA_range=inside,
                             all_optimal=all(sim["status"][m][c] == "optimal" for c in cores)))
    repro = pd.DataFrame(rows)
    repro.to_csv(HERE / "reproduction_check_r232.csv", index=False)
    print(repro.groupby("solver")[["TMA_reproduced", "objective_reproduced", "existing_TMA_within_FVA_range"]].sum())

    # ------------------------------------------------------------- 3. diagnostics
    print("3. per-model diagnostics")
    phylo20 = json.load(open(PHYLO_2020))
    fva_rows, diag_rows = [], []
    for m in mags15:
        k = m.replace(".contigs", "")
        cls = {c: classify(fva_run["TMA_min"][m][c], fva_run["TMA_max"][m][c], existing["TMA"][m][c]) for c in cores}
        cost = marginal_cost(m, ab_run)
        anaer = anaerobic_feasible(m, ab_run)
        for c in cores:
            fva_rows.append(dict(mag=k, iterativeID_r232=guilds.loc[k, "iterativeID_r232"], grdH=m in grdh9, core=c,
                                 abundance_as_run=ab_run.loc[k, c], TMA_existing=existing["TMA"][m][c],
                                 TMA_resim_glpk=resim["glpk"]["TMA"][m][c],
                                 TMA_resim_gurobi=resim.get("gurobi", {}).get("TMA", {}).get(m, {}).get(c, np.nan),
                                 TMA_fva_min=fva_run["TMA_min"][m][c], TMA_fva_max=fva_run["TMA_max"][m][c],
                                 regime=cls[c], lactate_cost_per_TMA=cost[c], feasible_without_O2=anaer[c]))
        g = guilds.loc[k]
        diag_rows.append(dict(
            mag=k, iterativeID_r232=g["iterativeID_r232"], iterativeID_2020=g["iterativeID_2020"],
            family_r232=g["family"], family_2020=phylo20.get(m + "__", {}).get("Family", ""),
            genus_r232=g["genus"], completeness=g["completeness"], contamination=g["contamination"],
            grdH_carrier=bool(g["betaine_reducer"]), ev_grdH=g["ev_grdH"] if isinstance(g["ev_grdH"], str) else "",
            existing_TMA_min_over_cores=min(existing["TMA"][m].values()),
            existing_TMA_max_over_cores=max(existing["TMA"][m].values()),
            n_cores_existing_at_0p1=sum(abs(v - BETAINE_UPTAKE) < SAT_TOL for v in existing["TMA"][m].values()),
            regimes=";".join(f"{c}:{cls[c]}" for c in cores),
            fva_max_over_cores=max(fva_run["TMA_max"][m].values()),
            lactate_cost_per_TMA_R1_A=cost["R1_A"],
            n_cores_feasible_without_O2=int(sum(anaer.values())),
            R1_A_block_rev_glycine_reductase=block_rev_glycine_reductase(m, ab_run, fva_solver),
            **diagnose_structure(m)))
    pd.DataFrame(fva_rows).to_csv(HERE / "fva_tma_by_core_r232.csv", index=False)
    diag = pd.DataFrame(diag_rows)
    diag.to_csv(HERE / "diagnostics_grdH_r232.csv", index=False)

    # ------------------------------------------------------------- 4. grdH-9 re-simulation
    print("4. 9 grdH carriers: same settings (GLPK, as-run abundances) + current-abundance sensitivity")
    grd_run = {k: {m: v[m] for m in grdh9} for k, v in resim["glpk"].items()}
    grd_cur = simulate(grdh9, ab_cur, "glpk")
    fva_cur = fva_tma(grdh9, ab_cur, fva_solver)

    def sub(d, mags):
        return {k: {m: d[k][m] for m in mags} for k in d}

    out = {
        "meta": {
            "description": "Betaine->TMA FBA fluxes (mmol/gDW/h per MAG and unrestored core). Modelled CH4 = 0.06 x TMA (Methods).",
            "models": "correlational/models/<MAG>.contigs__.RAST_gapfilled_betaine.json + ATPM_c0 (lb 0) if absent",
            "medium": {**MEDIA, "EX_cpd00540_e0": BETAINE_UPTAKE},
            "constraints": "bio1 fixed to core abundance; EX_cpd00159_e0 bounds (-1000,0) as objective (maximised = minimise lactate uptake)",
            "abundance_as_run": str(AB_AS_RUN.relative_to(REPO)) + "  (git show 22f4629:salterns/correlational/averaged_normalized_abundances_unrestored.csv)",
            "abundance_current": str(AB_CURRENT.relative_to(REPO)),
            "grdH9": [m.replace(".contigs", "") for m in grdh9],
            "cobra_version": cobra.__version__,
            "solvers": {"resim_primary": "glpk", "fva": fva_solver},
            "note": "No model was structurally changed: no annotation-supported fix exists for the zero-TMA grdH models (see README).",
        },
        "existing_grdH9": sub({"TMA": existing["TMA"], "betaine": existing["betaine"], "lactate": existing["lactate"]}, grdh9),
        "resim_glpk_asrun_15": resim["glpk"],
        "resim_gurobi_asrun_15": resim.get("gurobi"),
        "fva_asrun_15": fva_run,
        "grdH9_resim_glpk_asrun": grd_run,
        "grdH9_resim_glpk_current_abundance": grd_cur,
        "grdH9_fva_current_abundance": fva_cur,
    }
    json.dump(out, open(HERE / "fluxes_grdH_r232.json", "w"), indent=2)
    print("wrote fluxes_grdH_r232.json, reproduction_check_r232.csv, diagnostics_grdH_r232.csv, fva_tma_by_core_r232.csv")


if __name__ == "__main__":
    main()
