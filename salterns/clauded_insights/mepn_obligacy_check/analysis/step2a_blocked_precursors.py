"""Step 2a: which biomass precursors cannot be produced without MePn (condition (b): cell-5 medium,
EX_cpd25960_e0 lower bound 0) in the MePn-obligate models.
  demand test   : bio1 = 0, add DM_X (X -> ), maximise; blocked if max < 1e-9 (under (b)); also under (a) as control
  expansion     : for blocked pseudo-metabolites (Protein, RNA, DNA, lipids, cell wall ...) the reactants of the
                  reaction(s) that make them are demand-tested in the same way (recursively, depth <= 2)
  rescue test   : under (b), add a free source of X (-> X, <= 1000) and maximise bio1 (> 1e-6 = rescued)
  removal test  : under (b), drop X from bio1 (or from its pseudo-reaction) and maximise bio1
Output: step2a_blocked_precursors.json"""
import json
import sys
from multiprocessing import Pool

sys.dont_write_bytecode = True
sys.path.insert(0, "/tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns/a24fbb2e-9d33-5e05-aeea-a07d777ccde7/scratchpad/mepn_verify/analysis")
from cobra import Reaction  # noqa: E402

import vcommon as V  # noqa: E402

OBLIGATE = ["Roseovarius.9", "Roseovarius.11", "Albimonas.1"]
EPS = 1e-9


def condition(model, mepn_open):
    V.set_medium_strict(model)
    if not mepn_open:
        model.reactions.EX_cpd25960_e0.lower_bound = 0


def demand(model, met, mepn_open):
    with model:
        condition(model, mepn_open)
        model.reactions.bio1.bounds = (0, 0)
        dm = Reaction("DMtest", lower_bound=0, upper_bound=1000)
        dm.add_metabolites({met: -1})
        model.add_reactions([dm])
        model.objective = "DMtest"
        v = model.slim_optimize(error_value=0.0)
    return v


def rescue(model, met):
    with model:
        condition(model, False)
        src = Reaction("SRCtest", lower_bound=0, upper_bound=1000)
        src.add_metabolites({met: 1})
        model.add_reactions([src])
        model.objective = "bio1"
        v = model.slim_optimize(error_value=0.0)
    return v


def removal(model, rxn, met):
    with model:
        condition(model, False)
        r = model.reactions.get_by_id(rxn.id)
        c = r.metabolites[met]
        r.add_metabolites({met: -c})
        model.objective = "bio1"
        v = model.slim_optimize(error_value=0.0)
    return v


def analyse(iid):
    inv = {V.iterative_id(m): m for m in V.model_ids()}
    mid = inv[iid]
    model = V.load(mid)
    bio = model.reactions.bio1
    res = {"model": mid, "iterative_id": iid, "bio1_reactants": {}, "expanded": {}}
    blocked_top = []
    for met, coef in bio.metabolites.items():
        if coef >= 0:
            continue
        b = demand(model, met, False)
        a = demand(model, met, True)
        res["bio1_reactants"][met.id] = {"name": met.name, "coef": coef, "demand_max_without_MePn": b,
                                         "demand_max_with_MePn": a, "blocked_without_MePn": b < EPS}
        if b < EPS:
            blocked_top.append(met)
            res["bio1_reactants"][met.id]["rescue_bio1_if_supplied"] = rescue(model, met)
            res["bio1_reactants"][met.id]["bio1_if_removed_from_bio1"] = removal(model, bio, met)
    # expand blocked pseudo-metabolites through their producing reactions (depth 2)
    frontier = [(m, 0) for m in blocked_top]
    seen = set(m.id for m in blocked_top)
    while frontier:
        met, depth = frontier.pop(0)
        if depth >= 2:
            continue
        producers = [r for r in met.reactions if r.metabolites[met] > 0 and r.id != "bio1" and not r.boundary]
        if len(producers) != 1:      # only expand single-producer (pseudo-)metabolites
            res["expanded"][met.id] = {"name": met.name, "n_producers": len(producers),
                                       "producers": [r.id for r in producers][:20]}
            continue
        pr = producers[0]
        ent = {"name": met.name, "producer": pr.id, "producer_name": pr.name, "reactants": {}}
        for m2, c2 in pr.metabolites.items():
            if c2 >= 0:
                continue
            b = demand(model, m2, False)
            ent["reactants"][m2.id] = {"name": m2.name, "coef": c2, "demand_max_without_MePn": b,
                                       "demand_max_with_MePn": demand(model, m2, True), "blocked_without_MePn": b < EPS}
            if b < EPS:
                ent["reactants"][m2.id]["rescue_bio1_if_supplied"] = rescue(model, m2)
                ent["reactants"][m2.id]["bio1_if_removed_from_producer"] = removal(model, pr, m2)
                if m2.id not in seen:
                    seen.add(m2.id)
                    frontier.append((m2, depth + 1))
        res["expanded"][met.id] = ent
    return res


if __name__ == "__main__":
    with Pool(3) as p:
        out = {r["iterative_id"]: r for r in p.map(analyse, OBLIGATE)}
    json.dump(out, open(V.OUT / "step2a_blocked_precursors.json", "w"), indent=1)
    for iid, r in out.items():
        print(f"== {iid} ({r['model']})")
        for mid_, d in r["bio1_reactants"].items():
            if d["blocked_without_MePn"]:
                print(f"  bio1 reactant BLOCKED w/o MePn: {mid_} {d['name']} coef {d['coef']:.5f} | demand with MePn "
                      f"{d['demand_max_with_MePn']:.4g} | rescue-if-supplied bio1 {d['rescue_bio1_if_supplied']:.4f} | "
                      f"bio1 if removed {d['bio1_if_removed_from_bio1']:.4f}")
        for mid_, e in r["expanded"].items():
            if "reactants" not in e:
                print(f"  {mid_} {e['name']}: {e['n_producers']} producers {e['producers']}")
                continue
            bl = {k: v for k, v in e["reactants"].items() if v["blocked_without_MePn"]}
            print(f"  {mid_} {e['name']} <- {e['producer']} ({e['producer_name']}): {len(e['reactants'])} reactants, "
                  f"{len(bl)} blocked w/o MePn")
            for k, v in bl.items():
                print(f"     BLOCKED {k} {v['name']} coef {v['coef']:.5f} | demand with MePn {v['demand_max_with_MePn']:.4g} | "
                      f"rescue {v['rescue_bio1_if_supplied']:.4f} | removed-from-producer bio1 {v['bio1_if_removed_from_producer']:.4f}")
