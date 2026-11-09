"""B6. Biology-anchored perturbations of the ORIGINAL models (in memory only; nothing written back), GLPK.
Growth = max bio1 on the cell-5 medium (verified from `git show 22f4629:salterns/correlational/analysis(1).ipynb`,
cell index 5). a = as is, b = MePn exchange closed, c = Pi exchange closed.
Perturbations, each motivated by biochemistry rather than by the template:
  P0  none (reproduce obligacy)
  P1  close the 5'-dA sink SK_cpd03091_c0 (in vivo 5'-dA is salvaged, not excreted at ~0.47 mmol/gDW per unit growth)
  P2  Fhs rxn00690 reversible (formyl-THF synthetase is physiologically reversible; gene present: R9 CDS.706, Alb CDS.2814)
  P3  glycine cleavage rxn06377 + rxn06600 reversible (GCS is reversible in vivo; GcvPTH genes present)
  P4  Amn rxn00131 forward only (AMP + H2O -> adenine + R5P, a hydrolysis; ModelSEED DB direction '>')
  P5  add homoserine kinase (pathway hole: ThrC/Hom present, HSK unannotated)   [ATP + homoserine -> ADP + H+ + O-P-homoserine]
  P6  add PRPP synthetase rxn00770 (R11 only; prs absent from the 81%-complete MAG)
For a perturbation that rescues growth without MePn, the cell-27 single-point methane per unit biomass (bio1 = 1,
ratio constraint MePn >= r * Pi, maximise EX_cpd00159_e0) is reported at r = 0.07 and 0.33."""
import sys, json
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_bio import *
import cobra
from cobra import Reaction
from optlang.symbolics import Zero
cobra.Configuration().solver = "glpk"

MEDIA = {'EX_cpd00067_e0': 1000.0, 'EX_cpd00058_e0': 1000.0, 'EX_cpd00013_e0': 1000.0, 'EX_cpd00244_e0': 1000.0,
         'EX_cpd00205_e0': 1000.0, 'EX_cpd00034_e0': 1000.0, 'EX_cpd11574_e0': 1000.0, 'EX_cpd00971_e0': 1000.0,
         'EX_cpd00048_e0': 1000.0, 'EX_cpd00030_e0': 1000.0, 'EX_cpd00305_e0': 100.0, 'EX_cpd00001_e0': 1000.0,
         'EX_cpd10516_e0': 1000.0, 'EX_cpd00007_e0': 1000.0, 'EX_cpd00159_e0': 100.0, 'EX_cpd25960_e0': 1000.0,
         "EX_cpd00009_e0": 100, 'EX_cpd00063_e0': 1000.0, 'EX_cpd00149_e0': 1000.0, 'EX_cpd00254_e0': 1000.0,
         'EX_cpd00099_e0': 1000.0}
MEPN, PI, CH4, LAC = "EX_cpd25960_e0", "EX_cpd00009_e0", "EX_cpd01024_e0", "EX_cpd00159_e0"

def grow(m, closed=()):
    with m:
        m.medium = MEDIA
        for r in closed:
            m.reactions.get_by_id(r).lower_bound = 0
        m.objective = "bio1"
        v = m.slim_optimize(error_value=float("nan"))
    return 0.0 if (v != v or abs(v) < 1e-9) else round(v, 4)

def triplet(m):
    return grow(m), grow(m, (MEPN,)), grow(m, (PI,))

def cell27_point(m, r):
    with m:
        c = m.problem.Constraint(Zero, lb=0, ub=None, name="Phosphate_ratio")
        m.add_cons_vars(c); m.solver.update()
        c.set_linear_coefficients({m.reactions.get_by_id(PI).forward_variable: r, m.reactions.get_by_id(PI).reverse_variable: -r,
                                   m.reactions.get_by_id(MEPN).forward_variable: -1, m.reactions.get_by_id(MEPN).reverse_variable: 1})
        m.medium = MEDIA
        m.reactions.bio1.bounds = (1, 1)
        m.reactions.get_by_id(LAC).bounds = (-1000, 0)
        m.objective = LAC
        s = m.optimize()
        if s.status != "optimal":
            return {"status": s.status}
        return {"status": s.status, "CH4": round(s.fluxes[CH4], 6), "MePn": round(s.fluxes[MEPN], 6), "Pi": round(s.fluxes[PI], 6),
                "MePn_per_Pi": round(s.fluxes[MEPN] / s.fluxes[PI], 3) if s.fluxes[PI] else None,
                "SK_5dA": round(s.fluxes.get("SK_cpd03091_c0", float("nan")), 6)}

def add_hsk(m):
    hs, php = m.metabolites.get_by_id("cpd00227_c0"), m.metabolites.get_by_id("cpd00809_c0")
    r = Reaction("HSK_hole_c0", "homoserine kinase (pathway hole)", lower_bound=0, upper_bound=1000)
    r.add_metabolites({m.metabolites.cpd00002_c0: -1, hs: -1, m.metabolites.cpd00008_c0: 1, m.metabolites.cpd00067_c0: 1, php: 1})
    m.add_reactions([r])

def add_prs(m):
    r = Reaction("PRS_c0", "PRPP synthetase", lower_bound=0, upper_bound=1000)
    r.add_metabolites({m.metabolites.cpd00002_c0: -1, m.metabolites.cpd00101_c0: -1, m.metabolites.cpd00018_c0: 1,
                       m.metabolites.cpd00067_c0: 1, m.metabolites.cpd00103_c0: 1})
    m.add_reactions([r])

def P(kind):
    def f(m):
        if kind == "P1": m.reactions.SK_cpd03091_c0.bounds = (0, 0)
        if kind == "P2" and "rxn00690_c0" in m.reactions: m.reactions.rxn00690_c0.lower_bound = -1000
        if kind == "P3":
            for x in ("rxn06377_c0", "rxn06600_c0"):
                if x in m.reactions: m.reactions.get_by_id(x).lower_bound = -1000
        if kind == "P2+P3":
            P("P2")(m); P("P3")(m)
        if kind == "P4" and "rxn00131_c0" in m.reactions: m.reactions.rxn00131_c0.lower_bound = 0
        if kind == "P5": add_hsk(m)
        if kind == "P5+P1": add_hsk(m); P("P1")(m)
        if kind == "P2+P1": P("P2")(m); P("P1")(m)
        if kind == "P6": add_prs(m)
        if kind == "P6+P4": add_prs(m); P("P4")(m)
        if kind == "P1+P4": P("P1")(m); P("P4")(m)
    return f

inv = {iid(m): m for m in model_ids()}
plan = {"Roseovarius.9": ["P0", "P1", "P2", "P3", "P2+P3", "P4", "P5", "P5+P1", "P2+P1"],
        "Albimonas.1": ["P0", "P1", "P2", "P3", "P2+P3", "P5", "P5+P1", "P2+P1"],
        "Roseovarius.11": ["P0", "P1", "P4", "P1+P4", "P2", "P3", "P6", "P6+P4"],
        "Roseovarius.8": ["P0", "P1", "P4"]}
out = {}
for t, kinds in plan.items():
    base = cobra.io.load_json_model(str(MEPN_MODELS / f"{inv[t]}.json"))
    out[t] = {}
    for k in kinds:
        m = base.copy()
        if k != "P0":
            P(k)(m)
        a, b, c = triplet(m)
        rec = {"a_as_is": a, "b_no_MePn": b, "c_no_Pi": c}
        if a > 0:
            rec["cell27_r0.07"] = cell27_point(m, 0.07)
            rec["cell27_r0.33"] = cell27_point(m, 0.33)
        out[t][k] = rec
        print(f"{t:15s} {k:6s} growth a/b/c {a:7.4f} {b:7.4f} {c:7.4f}  " +
              (f"r=0.07 {rec['cell27_r0.07']}  r=0.33 {rec['cell27_r0.33']}" if a > 0 else ""))
json.dump(out, open(VB / "b6_bio_perturbations.json", "w"), indent=1)
