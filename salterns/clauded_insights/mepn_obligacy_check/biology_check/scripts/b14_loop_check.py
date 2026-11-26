"""B14. Energy-loop check for the biology-motivated repairs of B6: cell-5 medium with NO carbon source (lactate,
MePn closed; MePn is a C1 donor), bio1 fixed at 0, maximise ATPM_c0. Any value > 0 = free-energy-generating cycle."""
import sys
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib.util
spec = importlib.util.spec_from_file_location("b6", str(Path(__file__).resolve().parent / "b6_bio_perturbations.py"))
# reuse only definitions from b6 (its main loop is at module level, so copy the needed bits instead)
from common_bio import *
import cobra
from cobra import Reaction
cobra.Configuration().solver = "glpk"
MEDIA = {'EX_cpd00067_e0': 1000.0, 'EX_cpd00058_e0': 1000.0, 'EX_cpd00013_e0': 1000.0, 'EX_cpd00244_e0': 1000.0,
         'EX_cpd00205_e0': 1000.0, 'EX_cpd00034_e0': 1000.0, 'EX_cpd11574_e0': 1000.0, 'EX_cpd00971_e0': 1000.0,
         'EX_cpd00048_e0': 1000.0, 'EX_cpd00030_e0': 1000.0, 'EX_cpd00305_e0': 100.0, 'EX_cpd00001_e0': 1000.0,
         'EX_cpd10516_e0': 1000.0, 'EX_cpd00007_e0': 1000.0, "EX_cpd00009_e0": 100, 'EX_cpd00063_e0': 1000.0,
         'EX_cpd00149_e0': 1000.0, 'EX_cpd00254_e0': 1000.0, 'EX_cpd00099_e0': 1000.0}   # cell-5 minus lactate and MePn
def add_hsk(m):
    r = Reaction("HSK_hole_c0", lower_bound=0, upper_bound=1000)
    r.add_metabolites({m.metabolites.cpd00002_c0: -1, m.metabolites.cpd00227_c0: -1, m.metabolites.cpd00008_c0: 1,
                       m.metabolites.cpd00067_c0: 1, m.metabolites.cpd00809_c0: 1}); m.add_reactions([r])
def add_prs(m):
    r = Reaction("PRS_c0", lower_bound=0, upper_bound=1000)
    r.add_metabolites({m.metabolites.cpd00002_c0: -1, m.metabolites.cpd00101_c0: -1, m.metabolites.cpd00018_c0: 1,
                       m.metabolites.cpd00067_c0: 1, m.metabolites.cpd00103_c0: 1}); m.add_reactions([r])
def fhs_rev(m): m.reactions.rxn00690_c0.lower_bound = -1000
def gcs_rev(m):
    for x in ("rxn06377_c0", "rxn06600_c0"): m.reactions.get_by_id(x).lower_bound = -1000
inv = {iid(x): x for x in model_ids()}
plan = {"Roseovarius.9": {"P0": None, "P2 Fhs rev": fhs_rev, "P3 GCS rev": gcs_rev, "P5 HSK": add_hsk},
        "Albimonas.1": {"P0": None, "P2 Fhs rev": fhs_rev, "P3 GCS rev": gcs_rev, "P5 HSK": add_hsk},
        "Roseovarius.11": {"P0": None, "P6 Prs": add_prs}}
for t, d in plan.items():
    base = cobra.io.load_json_model(str(MEPN_MODELS / f"{inv[t]}.json"))
    for lab, f in d.items():
        m = base.copy()
        if f: f(m)
        with m:
            m.medium = MEDIA
            m.reactions.bio1.bounds = (0, 0)
            m.objective = "ATPM_c0"
            v = m.slim_optimize(error_value=float("nan"))
        print(f"{t:15s} {lab:12s} max ATPM with no carbon source = {0.0 if abs(v) < 1e-9 else v:.6g}")
