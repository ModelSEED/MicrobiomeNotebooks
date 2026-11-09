"""Step 2c: (i) sink scan under (b): which single metabolite drain (X ->) restores growth without MePn;
(ii) presence / bounds / GPR / draft / template / RAST support of the canonical reactions around each blocked node
(glycine, threonine, one-carbon disposal, PRPP); (iii) a parsimonious flux trace with MePn (pFBA at bio1 = 1,
lactate uptake at its cell-27 minimum) through the MePn module, the THF one-carbon pool, purine synthesis and the
5'-deoxyadenosine drain.
Output: step2c_sinks_and_c1.json"""
import json
import sys
from multiprocessing import Pool

sys.dont_write_bytecode = True
sys.path.insert(0, "/tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns/a24fbb2e-9d33-5e05-aeea-a07d777ccde7/scratchpad/mepn_verify/analysis")
from cobra import Reaction  # noqa: E402
from cobra.flux_analysis import pfba  # noqa: E402

import annot as A  # noqa: E402
import vcommon as V  # noqa: E402

OBLIGATE = ["Roseovarius.9", "Roseovarius.11", "Albimonas.1"]
MODULE = ["rxs00001_c0", "rxn26456_c0", "rxn26457_c0", "rxs43120_c0", "rxn26447_c0", "rxn04651_c0"]
CANON = {
    # glycine formation
    "rxn00692": "GlyA serine hydroxymethyltransferase (EC 2.1.2.1)",
    "rxn00541": "threonine aldolase (EC 4.1.2.5)", "rxn04142": "L-allo-threonine aldolase (EC 4.1.2.48/49)",
    "rxn00908": "glycine cleavage system, overall (EC 1.4.4.2/2.1.2.10/1.8.1.4)",
    "rxn06377": "GcvP glycine decarboxylase (lipoylprotein)", "rxn06600": "GcvT aminomethyltransferase",
    "rxn06493": "Lpd dihydrolipoyl dehydrogenase (lipoylprotein)",
    "rxn00272": "alanine-glyoxylate aminotransferase (EC 2.6.1.44)",
    "rxn00275": "glycine:2-oxoglutarate aminotransferase (EC 2.6.1.4)",
    "rxn00424": "serine-glyoxylate aminotransferase (EC 2.6.1.45)",
    "rxn00276": "glycine:oxaloacetate aminotransferase (EC 2.6.1.35)",
    "rxn01068": "threonine 3-dehydrogenase (EC 1.1.1.103)", "rxn00274": "2-amino-3-ketobutyrate CoA ligase Kbl (EC 2.3.1.29)",
    "rxn00433": "sarcosine oxidase (EC 1.5.3.1)", "rxn05870": "sarcosine dehydrogenase (EC 1.5.8.3)",
    # threonine from aspartate
    "rxn00337": "aspartate kinase (EC 2.7.2.4)", "rxn01643": "aspartate-semialdehyde dehydrogenase (EC 1.2.1.11)",
    "rxn01301": "homoserine dehydrogenase NAD (EC 1.1.1.3)", "rxn01302": "homoserine dehydrogenase NADP (EC 1.1.1.3)",
    "rxn01300": "homoserine kinase ThrB (EC 2.7.1.39)", "rxn01069": "threonine synthase ThrC (EC 4.2.3.1)",
    # serine
    "rxn01101": "phosphoglycerate dehydrogenase (EC 1.1.1.95)", "rxn02914": "phosphoserine aminotransferase (EC 2.6.1.52)",
    "rxn00420": "phosphoserine phosphatase (EC 3.1.3.3)",
    # one-carbon (THF) pool and its outlets
    "rxn00907": "MTHFD NADP (EC 1.5.1.5)", "rxn00906": "MTHFD NAD (EC 1.5.1.15)", "rxn01211": "MTHFC cyclohydrolase (EC 3.5.4.9)",
    "rxn00691": "PurU 10-formyl-THF deformylase (EC 3.5.1.10)", "rxn00690": "formate-THF ligase (EC 6.3.4.3)",
    "rxn00371": "formate dehydrogenase NAD (EC 1.17.1.9/1.2.1.2)", "rxn00103": "formate dehydrogenase NADP (EC 1.2.1.43)",
    "rxn04954": "MTHFR NAD (EC 1.5.1.20)", "rxn00910": "MTHFR NADP (EC 1.5.1.20)", "rxn00909": "MTHFR FADH2",
    "rxn00693": "MetH methionine synthase (EC 2.1.1.13)", "rxn03052": "MetE (EC 2.1.1.14)",
    "rxn37797": "thymidylate synthase", "rxn00686": "dihydrofolate reductase",
    # PRPP and nucleoside salvage
    "rxn00770": "ribose-phosphate diphosphokinase PrsA (EC 2.7.6.1)", "rxn04651": "PhnN ribose-1,5-bisphosphokinase (EC 2.7.4.23)",
    "rxn00131": "AMP nucleosidase (EC 3.2.2.4)", "rxn00139": "adenine phosphoribosyltransferase (EC 2.4.2.7)",
    "rxn01138": "purine nucleoside phosphorylase, adenosine (EC 2.4.2.1)", "rxn00134": "adenosine kinase (EC 2.7.1.20)",
    "rxn01021": "MTA nucleosidase (EC 3.2.2.16)", "rxn00143": "SAH nucleosidase (EC 3.2.2.9)",
    "rxn01022": "MTA phosphorylase (EC 2.4.2.28)", "rxn00777": "ribose-5-phosphate isomerase",
    "rxn00790": "amidophosphoribosyltransferase PurF", "rxn02895": "GAR synthetase PurD",
}
C1_METS = ["cpd00125_c0", "cpd00347_c0", "cpd00201_c0", "cpd00345_c0", "cpd00087_c0", "cpd00047_c0", "cpd00033_c0",
           "cpd00054_c0", "cpd00161_c0", "cpd03091_c0", "cpd00128_c0", "cpd00017_c0", "cpd00103_c0", "cpd00847_c0"]


def sink_scan(model):
    with model:
        V.set_medium_strict(model)
        model.reactions.EX_cpd25960_e0.lower_bound = 0
        model.objective = "bio1"
        model.solver.configuration.timeout = 20
        sk = Reaction("SINKtest", lower_bound=0, upper_bound=0)
        model.add_reactions([sk])
        hits = {}
        for met in model.metabolites:
            if met.compartment != "c0":
                continue
            sk.add_metabolites({met: -1})
            sk.upper_bound = 1000
            v = model.slim_optimize(error_value=0.0)
            if v > 1e-6:
                hits[met.id] = {"name": met.name, "bio1": v}
            sk.upper_bound = 0
            sk.subtract_metabolites({met: -1})
    return hits


def trace(model):
    with model:
        V.set_medium_strict(model)
        model.reactions.bio1.bounds = (1, 1)
        ex = model.reactions.EX_cpd00159_e0
        ex.bounds = (-1000, 0)
        model.objective = ex
        opt = model.slim_optimize()
        ex.bounds = (opt - 1e-9, opt + 1e-9)
        s = pfba(model).fluxes
        out = {"lactate_per_biomass": opt, "MePn": s["EX_cpd25960_e0"], "Pi": s["EX_cpd00009_e0"],
               "CH4": s["EX_cpd01024_e0"], "module": {r: s[r] for r in MODULE if r in model.reactions},
               "SK_5dA": s.get("SK_cpd03091_c0", None), "by_metabolite": {}}
        for mid in C1_METS:
            if mid not in model.metabolites:
                continue
            met = model.metabolites.get_by_id(mid)
            out["by_metabolite"][mid] = {"name": met.name, "reactions": {
                r.id: {"net_production": s[r.id] * r.metabolites[met], "flux": s[r.id], "bounds": r.bounds,
                       "GPR_empty": not r.gene_reaction_rule, "equation": r.build_reaction_string(True)[:160]}
                for r in met.reactions if abs(s[r.id]) > 1e-9}}
    return out


def analyse(iid):
    inv = {V.iterative_id(m): m for m in V.model_ids()}
    mid = inv[iid]
    model = V.load(mid)
    draft = json.load(open(V.DRAFTS / f"{mid}__.RAST.json"))
    draft_rx = {r["id"]: r for r in draft["reactions"]}
    prov = json.load(open(V.OUT / "step0_provenance.json"))[mid]
    canon = {}
    for base, desc in CANON.items():
        rid = base + "_c0"
        rec = {"desc": desc, "in_model": rid in model.reactions, "in_draft": rid in draft_rx,
               "gapfilled": rid in prov["gapfilled"]}
        if rid in model.reactions:
            r = model.reactions.get_by_id(rid)
            rec.update({"bounds": r.bounds, "GPR": r.gene_reaction_rule.replace(mid + "__.RAST.", ""),
                        "equation": r.build_reaction_string(True)})
        rec["template_and_RAST"] = A.support(mid, base + "_c")
        canon[base] = rec
    return {"model": mid, "iterative_id": iid, "sink_rescuers": sink_scan(model), "canonical": canon,
            "trace_with_MePn_bio1_eq_1": trace(model)}


if __name__ == "__main__":
    with Pool(3) as p:
        out = {r["iterative_id"]: r for r in p.map(analyse, OBLIGATE)}
    json.dump(out, open(V.OUT / "step2c_sinks_and_c1.json", "w"), indent=1, default=str)
    for iid, r in out.items():
        print(f"== {iid}: {len(r['sink_rescuers'])} single-sink rescuers")
        for k, v in sorted(r["sink_rescuers"].items(), key=lambda x: -x[1]["bio1"])[:40]:
            print(f"   sink-rescuer {k:14s} {v['name'][:50]:50s} bio1 {v['bio1']:.3f}")
        print("   canonical reactions:")
        for base, c in r["canonical"].items():
            ta = c["template_and_RAST"]
            print(f"     {base} {c['desc'][:52]:52s} model={str(c['in_model'])[0]} draft={str(c['in_draft'])[0]} GF={str(c['gapfilled'])[0]} "
                  f"bounds={c.get('bounds')} GPR={'yes' if c.get('GPR') else ('-' if not c['in_model'] else 'EMPTY')} "
                  f"tpl={ta.get('template_type')} RAST={ta.get('rast_supported')} {list(ta.get('rast_supported_roles', {}).items())[:1]}")
        t = r["trace_with_MePn_bio1_eq_1"]
        print(f"   trace (bio1=1): lactate {t['lactate_per_biomass']:.4f} MePn {t['MePn']:.4f} Pi {t['Pi']:.4f} CH4 {t['CH4']:.4f} "
              f"SK_5dA {t['SK_5dA']} module {({k: round(v, 4) for k, v in t['module'].items()})}")
