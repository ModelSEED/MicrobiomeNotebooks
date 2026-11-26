"""Step 1b: is methane unique at the cell-27 lactate optimum?  For each model x ratio x core: solve cell 27,
then fix EX_cpd00159_e0 at the optimum (tolerance 1e-9 relative) and min/max EX_cpd01024_e0, EX_cpd25960_e0,
EX_cpd00009_e0 (GLPK, 60 s per LP).  Models: argv (iterative IDs) or the 3 non-binding + 3 ratio-bound controls.
Output: step1b_fva.json"""
import json
import sys
from multiprocessing import Pool

sys.dont_write_bytecode = True
sys.path.insert(0, "/tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns/a24fbb2e-9d33-5e05-aeea-a07d777ccde7/scratchpad/mepn_verify/analysis")
import vcommon as V  # noqa: E402

TARGETS = ["Roseovarius.9", "Roseovarius.11", "Albimonas.1", "Roseovarius.10", "Roseibacterium.1", "Thiolinea.1"]


def fva_model(mid, model=None, tag="original"):
    model = model if model is not None else V.load(mid)
    row = V.abundances_as_run().loc[V.mag_of(mid)]
    out = {}
    for r in V.RATIOS:
        m = model.copy()
        V.add_ratio_constraint(m, r)
        V.set_medium_strict(m)
        m.solver.configuration.timeout = 60
        out[str(r)] = {}
        for core, a in row.items():
            m.reactions.bio1.bounds = (a, a)
            ex = m.reactions.EX_cpd00159_e0
            ex.bounds = (-1000, 0)
            m.objective = ex
            opt = m.slim_optimize(error_value=float("nan"))
            rec = {"lactate_opt": opt, "status": m.solver.status}
            if opt == opt:
                with m:
                    tol = max(abs(opt) * 1e-9, 1e-12)
                    ex.bounds = (opt - tol, min(0.0, opt + tol))
                    for rid in (V.CH4_EX, V.MEPN_EX, V.PI_EX):
                        rr = m.reactions.get_by_id(rid)
                        m.objective = rr
                        m.objective_direction = "min"
                        lo = m.slim_optimize(error_value=float("nan"))
                        m.objective_direction = "max"
                        hi = m.slim_optimize(error_value=float("nan"))
                        rec[rid] = [lo, hi]
                    m.objective_direction = "max"
                rec["CH4_per_abundance_range"] = [rec[V.CH4_EX][0] / a, rec[V.CH4_EX][1] / a]
            out[str(r)][core] = rec
            m.reactions.bio1.bounds = (0, 0)
    return {"model": mid, "iterative_id": V.iterative_id(mid), "tag": tag, "fva": out}


if __name__ == "__main__":
    want = sys.argv[1:] or TARGETS
    inv = {V.iterative_id(m): m for m in V.model_ids()}
    with Pool(min(6, len(want))) as p:
        res = p.map(fva_model, [inv[w] for w in want])
    out = {r["iterative_id"]: r for r in res}
    json.dump(out, open(V.OUT / "step1b_fva.json", "w"), indent=1)
    for iid, r in out.items():
        for rat, cores in r["fva"].items():
            rng = [c["CH4_per_abundance_range"] for c in cores.values() if "CH4_per_abundance_range" in c]
            lo = min(x[0] for x in rng); hi = max(x[1] for x in rng)
            spread = max(x[1] - x[0] for x in rng)
            print(f"{iid:18s} r={rat}: CH4/abundance FVA over cores [{lo:.6f}, {hi:.6f}], max within-core spread {spread:.2e}; "
                  f"statuses {set(c['status'] for c in cores.values())}")
