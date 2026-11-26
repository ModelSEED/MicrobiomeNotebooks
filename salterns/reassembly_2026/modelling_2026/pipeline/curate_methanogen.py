"""Curate the 2026 Methanolobus (R1_A_D1__concoct.bin.007) arm-A model into a gene-supported methylotrophic methanogen,
then validate it. Run locally (a few LPs; ~1 min):

    ~/Documents/py_venv/bin/python -I pipeline/curate_methanogen.py

Inputs : models/methanogen/R1_A_D1__concoct.bin.007.json (arm A, protocol of record; md5 checked)
         models/reports/methanogen/R1_A_D1__concoct.bin.007.report.json (ATP-test thresholds)
         annotations/prodigal/R1_A_D1__concoct.bin.007.genes.tsv (gene ids)
Outputs: models/methanogen/R1_A_D1__concoct.bin.007__curated.json
         models/reports/methanogen/R1_A_D1__concoct.bin.007__curated.validation.json
Every change, with its evidence, is listed in models/methanogen/CURATION.md.
"""
import hashlib
import importlib.util
import json
import logging
import sys
from pathlib import Path

import cobra
import pandas as pd
from cobra import Metabolite, Reaction
from cobra.flux_analysis import flux_variability_analysis

ROOT = Path(__file__).resolve().parents[1]
MAG = "R1_A_D1__concoct.bin.007"
SRC = ROOT / "models/methanogen" / f"{MAG}.json"
SRC_MD5 = "e8fd1042225c7b0e906a5cd11c8c2b05"          # models/MANIFEST.md5 entry of the kept arm-A model
OUT = ROOT / "models/methanogen" / f"{MAG}__curated.json"
VAL = ROOT / "models/reports/methanogen" / f"{MAG}__curated.validation.json"

# Rnf ion coupling (Na+ and H+ translocated per 2 e- from Fd to methanophenazine); see CURATION.md "Ion stoichiometry".
RNF_NA, RNF_H = 3, 0

spec = importlib.util.spec_from_file_location("bm", ROOT / "pipeline/build_model.py")
bm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bm)

SUBS = {"methanol": "EX_cpd00116_e0", "MMA": "EX_cpd00187_e0", "DMA": "EX_cpd00425_e0", "TMA": "EX_cpd00441_e0"}
METHYLS = {"EX_cpd00116_e0": 1, "EX_cpd00187_e0": 1, "EX_cpd00425_e0": 2, "EX_cpd00441_e0": 3}
# textbook stoichiometry per mol substrate (Thauer 1998; Zinder 1993): CH4, CO2, NH3
STOICH = {"methanol": (0.75, 0.25, 0.0), "MMA": (0.75, 0.25, 1.0), "DMA": (1.5, 0.5, 1.0), "TMA": (2.25, 0.75, 1.0)}
H2CO2 = {**bm._MINERALS, "EX_cpd11640_e0": 1000.0, "EX_cpd00011_e0": 100.0, "EX_cpd00239_e0": 1000.0}

_genes = pd.read_csv(ROOT / f"annotations/prodigal/{MAG}.genes.tsv", sep="\t").gene.str.split("|").str[1]
FULL = {g.split("_length_")[0] + "_" + g.rsplit("_", 1)[1]: g for g in _genes}   # NODE_296_14 -> model gene id


def G(*short):
    return [FULL[s] for s in short]


def only(sub):
    return {k: v for k, v in bm.METHYL_MEDIUM.items() if k not in set(SUBS.values()) - {SUBS[sub]}}


def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


# ------------------------------------------------------------------------------------------------ curation
def curate(m, rnf_na=RNF_NA, rnf_h=RNF_H):
    met = m.metabolites.get_by_id
    log = []

    def eq(r):
        return r.build_reaction_string()

    def add(rid, name, stoich, bounds, gpr, ref):
        r = Reaction(rid, name=name, lower_bound=bounds[0], upper_bound=bounds[1])
        r.add_metabolites({met(k): v for k, v in stoich.items()})
        m.add_reactions([r])
        r.gene_reaction_rule = gpr
        r.notes["curation"] = ref
        bal = r.check_mass_balance()
        assert not bal, (rid, bal)
        log.append({"action": "added", "id": rid, "equation": eq(r), "bounds": list(r.bounds), "gpr": gpr, "note": ref})

    def remove(rid, why):
        r = m.reactions.get_by_id(rid)
        log.append({"action": "removed", "id": rid, "equation": eq(r), "bounds": list(r.bounds),
                    "gpr": r.gene_reaction_rule, "gapfilled": r.notes.get("gapfilled"), "note": why})
        m.remove_reactions([r])

    def rebound(rid, b, why):
        r = m.reactions.get_by_id(rid)
        log.append({"action": "bounds", "id": rid, "equation": eq(r), "before": list(r.bounds), "after": list(b), "note": why})
        r.bounds = b
        r.notes["curation"] = why

    def regpr(rid, gpr, why):
        r = m.reactions.get_by_id(rid)
        log.append({"action": "gpr", "id": rid, "equation": eq(r), "before": r.gene_reaction_rule, "after": gpr, "note": why})
        r.gene_reaction_rule = gpr
        r.notes["curation"] = why

    # -- 1. Fpo: template rxn24614 runs on the generic F420 pool (cpd28067/cpd27742) that no H4MPT-branch enzyme uses,
    #       is '<' only, and translocates no H+. Replaced by the same enzyme on coenzyme F420 (cpd00792/cpd00649), 2 H+/2e-.
    fpo_gpr = m.reactions.rxn24614_c0.gene_reaction_rule
    remove("rxn24614_c0", "superseded by MCUR_FPO_c0 (wrong F420 pool, '<' only, no H+ translocation)")
    add("MCUR_FPO_c0", "F420H2:methanophenazine oxidoreductase Fpo (EC 1.5.98.3), 2 H+ translocated per 2 e-",
        {"cpd00792_c0": -1, "cpd08701_c0": -1, "cpd00067_c0": -2, "cpd00649_c0": 1, "cpd08702_c0": 1, "cpd00067_e0": 2},
        (-1000, 1000), fpo_gpr, "ModelSEED rxn24614 written on coenzyme F420 cpd00792/cpd00649 + 2 H+ c->e (Baumer et al. 2000)")
    # -- 2. Fwd/Fmd: rxn11938 is NAD-linked (no such enzyme); the ferredoxin form (DB rxn46184) is on generic pools.
    fwd_gpr = (f"{G('NODE_4070_15')[0]} and {G('NODE_4070_16')[0]} and ({' or '.join(G('NODE_404_55', 'NODE_192_9'))}) "
               f"and ({' or '.join(G('NODE_404_56', 'NODE_192_10'))})")
    remove("rxn11938_c0", "NAD-dependent formylmethanofuran dehydrogenase does not exist; replaced by ferredoxin form MCUR_FWD_c0")
    add("MCUR_FWD_c0", "formylmethanofuran dehydrogenase (ferredoxin) Fmd/Fwd (EC 1.2.7.12)",
        {"cpd00735_c0": -1, "cpd00001_c0": -1, "cpd11621_c0": -2, "cpd00011_c0": 1, "cpd00643_c0": 1, "cpd11620_c0": 2,
         "cpd00067_c0": 1}, (-1000, 1000), fwd_gpr,
        "ModelSEED rxn46184 written on the model's cpd00735/cpd00643 and ferredoxin cpd11620/cpd11621")
    # -- 3. HdrDE: rxn03126 carried the soluble HdrA2B2C2 genes; the membrane HdrDE genes were mapped to nothing.
    hdrde = f"({' or '.join(G('NODE_664_33', 'NODE_1064_33'))}) and {G('NODE_664_32')[0]}"
    remove("rxn03126_c0", "superseded by MCUR_HDRDE_c0 (GPR held HdrA2B2C2 genes; no H+ translocation; reversible)")
    add("MCUR_HDRDE_c0", "CoB-CoM heterodisulfide reductase HdrDE (EC 1.8.98.1), methanophenazine, 2 H+ translocated",
        {"cpd02935_c0": -1, "cpd08702_c0": -1, "cpd00067_c0": -3, "cpd02246_c0": 1, "cpd02817_c0": 1, "cpd08701_c0": 1,
         "cpd00067_e0": 2}, (0, 1000), hdrde, "ModelSEED rxn03126 + 2 H+ c->e (Ide et al. 1999), forward only")
    # -- 4. Mcr: methane-forming direction only (dG0' = -30 kJ/mol; Thauer 1998)
    rebound("rxn03127_c0", (0, 1000), "methyl-CoM reductase, methane-forming direction only")
    # -- 5. Rnf: re-enable Fd -> methanophenazine (filtered by ATP correction), set its ion coupling, block reverse
    r = m.reactions.rxa45615_c0
    before = {"equation": eq(r), "bounds": list(r.bounds), "gpr": r.gene_reaction_rule}
    r.subtract_metabolites(dict(r.metabolites))
    r.add_metabolites({met("cpd11620_c0"): -2, met("cpd08701_c0"): -1, met("cpd00067_c0"): -(2 + rnf_h),
                       met("cpd00971_c0"): -rnf_na, met("cpd11621_c0"): 2, met("cpd08702_c0"): 1,
                       **({met("cpd00067_e0"): rnf_h} if rnf_h else {}), met("cpd00971_e0"): rnf_na})
    assert not r.check_mass_balance()
    r.bounds = (0, 1000)
    rnf_gpr = " and ".join(G("NODE_835_25", "NODE_835_24", "NODE_835_29", "NODE_835_28", "NODE_835_26", "NODE_835_27",
                             "NODE_835_30"))
    r.gene_reaction_rule = rnf_gpr
    r.notes["curation"] = (f"Rnf (rnfABCDEG + mmcA) Fd -> MP, {rnf_na} Na+/{rnf_h} H+ per 2 e-, forward only "
                           "(template: 3 Na+ + 2 H+, reverse only after ATP correction)")
    log.append({"action": "restoichiometried+bounds+gpr", "id": "rxa45615_c0", "before": before,
                "after": {"equation": eq(r), "bounds": list(r.bounds), "gpr": rnf_gpr}, "note": r.notes["curation"]})
    # -- 6. F420H2 -> Fd without energy input is uphill; keep only Fd -> F420
    rebound("rxn40355_c0", (-1000, 0), "F420:ferredoxin oxidoreductase, Fd(red) -> F420 direction only (uphill reverse)")
    # -- 7. gene-less hydrogenase / formate hydrogenlyase artefacts from ATP correction (core) and the CH4 gap-fill
    remove("rxn05759_c0", "no Ech or Fd-reducing hydrogenase gene; added by ATP-correction core gap-fill (H2.CO2 test)")
    remove("rxn06299_c0", "no FrhA/FrhG (only FrhB/FpoF/FdhB-like proteins); added by ATP-correction core gap-fill")
    remove("rxn08518_c0", "no FdhAB/Hyc genes; added by the CH4 gap-fill (gapfill_2) only to make CH4 possible apart from growth")
    # -- 7b. complex I on ubiquinone, GPR = FpoK and FpoL (RAST named them "NADH-ubiquinone oxidoreductase chain K/L"):
    #        Fpo subunits mapped to a quinone pump in an organism without quinones; with succinate it is an energy source
    remove("rxn10122_c0", "GPR is Fpo subunits FpoK/FpoL (NODE_863_35/36), which belong to MCUR_FPO; methanogens have no "
                          "ubiquinone; with the gene-associated quinone SDH rxn09272 it lets succinate alone yield ATP")
    # -- 8. Fno: present (RAST + eggNOG K06988) but the template role name differs, so never mapped
    add("rxn16596_c0", "reduced coenzyme F420:NADP+ oxidoreductase Fno (EC 1.5.1.40)",
        {"cpd00006_c0": -1, "cpd00792_c0": -1, "cpd00005_c0": 1, "cpd00067_c0": 1, "cpd00649_c0": 1}, (0, 1000),
        G("NODE_209_20")[0], "template rxn16596 (identical equation); F420H2 -> NADPH direction only")
    # -- 9. GPR corrections (no flux change)
    mtmb = f"(({' and '.join(G('NODE_296_14', 'NODE_296_15'))}) or ({' and '.join(G('NODE_296_29', 'NODE_296_30'))}))"
    mtmc = f"({' or '.join(G('NODE_296_13', 'NODE_296_28'))})"
    regpr("rxn24608_c0", f"{mtmb} and {mtmc}", "MtmB is one Pyl protein split at the in-frame UAG (codon 202): pieces ANDed")
    regpr("rxn24610_c0", f"({' and '.join(G('NODE_296_26', 'NODE_296_27'))}) and {G('NODE_296_25')[0]}",
          "MtbB split at the in-frame UAG (codon 376): pieces ANDed; MtbC")
    regpr("rxn24612_c0", f"({' and '.join(G('NODE_296_19', 'NODE_296_20'))}) and {G('NODE_296_21')[0]}",
          "MttB split at the in-frame UAG (codon 334): pieces ANDed; MttC")
    atp = " and ".join(G("NODE_2103_1", "NODE_863_1", "NODE_2103_3", "NODE_863_2", "NODE_2103_4", "NODE_2103_2",
                         "NODE_2103_6", "NODE_2103_5"))
    regpr("rxn08173_c0", atp, "A1AO ATP synthase subunits A,B,C,D,E,F,I,K are one complex (was OR)")
    na = m.reactions.rxn05209_c0.gene_reaction_rule
    mrp = " and ".join(G("NODE_192_92", "NODE_8305_8", "NODE_8305_7", "NODE_8305_6", "NODE_8305_5", "NODE_8305_4"))
    regpr("rxn05209_c0", f"{na} or ({mrp})", "add Mrp Na+/H+ antiporter mrpACDEFG (Jasso-Chavez et al. 2013)")
    m.notes["salterns_2026_curation"] = json.dumps({"source": SRC.name, "rnf_na": rnf_na, "rnf_h": rnf_h,
                                                    "n_changes": len(log), "doc": "models/methanogen/CURATION.md"})
    return m, log


# ------------------------------------------------------------------------------------------------ validation
def opt(m, med, objective="bio1", fix=None, closed_ex_up=()):
    with m:
        bm.set_medium(m, med)
        for rid, b in (fix or {}).items():
            m.reactions.get_by_id(rid).bounds = b
        for rid in closed_ex_up:
            if rid in m.reactions:
                m.reactions.get_by_id(rid).upper_bound = 0
        m.objective = objective
        s = m.optimize()
        if s.status != "optimal":
            return {"status": s.status}
        ex = {r.id: float(s.fluxes[r.id]) for r in m.exchanges if abs(s.fluxes[r.id]) > 1e-9}
        return {"status": "optimal", "objective": float(s.objective_value), "exchanges": ex}


def growth_row(m, med, closed=()):
    r = opt(m, med, closed_ex_up=closed)
    if r["status"] != "optimal":
        return {"growth": 0.0, "status": r["status"]}
    ex = r["exchanges"]
    out = {"growth": r["objective"], "CH4": ex.get("EX_cpd01024_e0", 0.0), "acetate": ex.get("EX_cpd00029_e0", 0.0),
           "CO2": ex.get("EX_cpd00011_e0", 0.0), "H2": ex.get("EX_cpd11640_e0", 0.0),
           "formate": ex.get("EX_cpd00047_e0", 0.0),
           "uptake": {k: -ex.get(k, 0.0) for k in SUBS.values() if med.get(k, 0) > 0}}
    up = sum(METHYLS[k] * v for k, v in out["uptake"].items())
    out["CH4_per_methyl"] = out["CH4"] / up if up > 1e-9 else None
    return out


def catabolic(m, sub, block_ch4=False):
    """substrate uptake fixed at 1, no biomass, maximise ATP: products per mol substrate"""
    r = opt(m, only(sub), "ATPM_c0", fix={SUBS[sub]: (-1, -1), "bio1": (0, 0),
                                         **({"EX_cpd01024_e0": (0, 0)} if block_ch4 else {})})
    if r["status"] != "optimal":
        return {"status": r["status"]}
    ex = r["exchanges"]
    return {"ATP": r["objective"], "CH4": ex.get("EX_cpd01024_e0", 0.0), "CO2": ex.get("EX_cpd00011_e0", 0.0),
            "NH3": ex.get("EX_cpd00013_e0", 0.0), "acetate": ex.get("EX_cpd00029_e0", 0.0),
            "H2": ex.get("EX_cpd11640_e0", 0.0), "formate": ex.get("EX_cpd00047_e0", 0.0)}


def kvalue(m, med, mu=1e-3, weights=None):
    """bio1 fixed at mu, substrate uptake minimised (weighted), k = CH4 / bio1 (mmol CH4 per g biomass); FVA at optimum"""
    with m:
        bm.set_medium(m, med)
        m.reactions.bio1.bounds = (mu, mu)
        subs = [k for k in METHYLS if med.get(k, 0) > 0 and k in m.reactions]
        w = weights or {k: METHYLS[k] for k in subs}
        m.objective = {m.reactions.get_by_id(k): w[k] for k in subs}
        m.objective_direction = "max"                      # exchanges are negative for uptake
        s = m.optimize()
        if s.status != "optimal":
            return {"status": s.status}
        o = s.objective_value
        res = {"status": "optimal", "objective": "minimise " + " + ".join(f"{w[k]}*uptake({k})" for k in subs),
               "uptake_per_g_biomass": {k: float(-s.fluxes[k] / mu) for k in subs},
               "methyl_uptake_per_g_biomass": float(sum(-METHYLS[k] * s.fluxes[k] for k in subs) / mu),
               "k_mmol_CH4_per_g": float(s.fluxes.EX_cpd01024_e0 / mu)}
        m.add_cons_vars(m.problem.Constraint(m.objective.expression, lb=o - 1e-7 * abs(o)))
        m.objective = "bio1"
        ids = [i for i in ("EX_cpd01024_e0", "EX_cpd00029_e0", "EX_cpd00047_e0", "EX_cpd11640_e0", "EX_cpd00011_e0")
               if i in m.reactions]
        f = flux_variability_analysis(m, ids, fraction_of_optimum=0, processes=1)
        res["FVA_per_g_biomass"] = {i: [float(f.loc[i, "minimum"] / mu), float(f.loc[i, "maximum"] / mu)] for i in ids}
        res["yield_g_per_mol_CH4"] = 1000.0 / res["k_mmol_CH4_per_g"] if res["k_mmol_CH4_per_g"] > 0 else None
        return res


def atp_tests(m):
    from modelseedpy.core.msatpcorrection import load_default_medias
    from modelseedpy.core.msmodelutl import MSModelUtil
    rep = json.load(open(ROOT / f"models/reports/methanogen/{MAG}.report.json"))["build"]
    meds = {x.id: x for x, _ in load_default_medias()}
    mc = m.copy()
    u = MSModelUtil.get(mc)
    out = {}
    for t in rep["atp_correction"]["tests"]:
        cond = {"media": meds[t["media"]], "is_max_threshold": True, "threshold": t["threshold"], "objective": "ATPM_c0"}
        with mc:
            u.apply_test_condition(cond, model=mc)
            v = mc.slim_optimize(error_value=float("nan"))
        out[t["media"]] = {"ATP": v, "threshold": t["threshold"], "pass": bool(v <= t["threshold"] + 1e-9)}
    return out


def validate(m, full=True):
    v = {"free_atp": bm.free_atp(m)}
    v["free_atp_pass"] = all(x == x and x <= 1e-6 for x in v["free_atp"].values())
    v["growth_max"] = {"mixed": growth_row(m, bm.METHYL_MEDIUM), **{s: growth_row(m, only(s)) for s in SUBS}}
    v["growth_acetate_export_blocked"] = {"mixed": growth_row(m, bm.METHYL_MEDIUM, closed=("EX_cpd00029_e0",)),
                                          **{s: growth_row(m, only(s), closed=("EX_cpd00029_e0",)) for s in SUBS}}
    v["manuscript_DMA_H2"] = growth_row(m, bm.MANUSCRIPT_METHANOGEN_MEDIUM)
    v["H2_CO2_only"] = growth_row(m, H2CO2)
    v["k"] = {**{s: kvalue(m, only(s)) for s in SUBS},
              "mixed_methyl_weighted": kvalue(m, bm.METHYL_MEDIUM),
              "mixed_unweighted": kvalue(m, bm.METHYL_MEDIUM, weights={k: 1 for k in METHYLS})}
    if full:
        v["catabolic_per_mol_substrate"] = {s: {"methanogenic": catabolic(m, s), "CH4_export_blocked": catabolic(m, s, True),
                                                "textbook_CH4_CO2_NH3": STOICH[s]} for s in SUBS}
        v["atp_tests"] = atp_tests(m)
        v["alt_substrates_ATP_per_mmol"] = alt_substrates(m)
        v["aerobic_community_media"] = aerobic_check(m)
    return v


def set_coupling(m, fpo_h, hdr_h, vht_h, mtr_na, rnf_na):
    """sensitivity only: rewrite the ion terms of Fpo, HdrDE, Vht, Mtr and Rnf (scalar H+ kept balanced)"""
    met = m.metabolites.get_by_id
    hc, he, nac, nae = (met("cpd00067_c0"), met("cpd00067_e0"), met("cpd00971_c0"), met("cpd00971_e0"))

    def setion(rid, ion_c, ion_e, n, scalar_c):
        r = m.reactions.get_by_id(rid)
        r.subtract_metabolites({k: v for k, v in r.metabolites.items() if k in (hc, he, nac, nae)})
        add = {ion_c: -n, ion_e: n}
        if scalar_c:
            add[hc] = add.get(hc, 0) - scalar_c
        r.add_metabolites(add)
        assert all(abs(x) < 1e-9 for x in r.check_mass_balance().values()), (rid, r.check_mass_balance())
    setion("MCUR_FPO_c0", hc, he, fpo_h, 0)
    setion("MCUR_HDRDE_c0", hc, he, hdr_h, 1)
    setion("rxn15961_c0", hc, he, vht_h, 0)
    setion("rxn03020_c0", nac, nae, mtr_na, 0)
    setion("rxa45615_c0", nac, nae, rnf_na, 2)


def alt_substrates(m):
    """ATP (and CH4) from 1 mmol of each single substrate on the anaerobic minerals, biomass 0"""
    out = {}
    cands = {"acetate": "EX_cpd00029_e0", "pyruvate": "EX_cpd00020_e0", "succinate": "EX_cpd00036_e0",
             "L-lactate": "EX_cpd00159_e0", "formate": "EX_cpd00047_e0", "ethanol": "EX_cpd00363_e0",
             "glycine betaine": "EX_cpd00540_e0", "H2 (+CO2)": "EX_cpd11640_e0", "methanol": "EX_cpd00116_e0"}
    for lab, ex in cands.items():
        if ex not in m.reactions:
            out[lab] = "no exchange"
            continue
        med = {**bm._MINERALS, ex: 1.0, "EX_cpd00011_e0": 100.0}
        r = opt(m, med, "ATPM_c0", fix={ex: (-1, -1), "bio1": (0, 0)})
        out[lab] = ({"ATP": r["objective"], "CH4": r["exchanges"].get("EX_cpd01024_e0", 0.0),
                     "acetate": r["exchanges"].get("EX_cpd00029_e0", 0.0)} if r["status"] == "optimal" else r["status"])
    return out


def aerobic_check(m):
    """community MEDIA (aerobic lactate + MePn + sulfate, no H2S) with MMA added; with and without H2S"""
    out = {}
    for lab, med in (("MEDIA + MMA", {**bm.MEDIA, "EX_cpd00187_e0": 100.0}),
                     ("MEDIA + MMA + H2S", {**bm.MEDIA, "EX_cpd00187_e0": 100.0, "EX_cpd00239_e0": 1000.0}),
                     ("MEDIA + H2S (no methylamine)", {**bm.MEDIA, "EX_cpd00239_e0": 1000.0})):
        g = growth_row(m, med)
        r = opt(m, med)
        g["O2_uptake"] = -r.get("exchanges", {}).get("EX_cpd00007_e0", 0.0) if r["status"] == "optimal" else None
        out[lab] = g
    return out


def methanol_what_if(m):
    """NOT part of the curated model: MtaB/MtaC are absent; adds a lumped methanol:CoM methyltransferase to a copy"""
    mc = m.copy()
    mc.add_metabolites([Metabolite("cpd00116_c0", formula="CH4O", charge=0, name="Methanol", compartment="c0")])
    d = Reaction("whatif_methanol_diffusion", lower_bound=-1000, upper_bound=1000)
    d.add_metabolites({mc.metabolites.cpd00116_e0: -1, mc.metabolites.cpd00116_c0: 1})
    t = Reaction("whatif_mtaCBA", lower_bound=0, upper_bound=1000)
    t.add_metabolites({mc.metabolites.cpd00116_c0: -1, mc.metabolites.cpd02246_c0: -1, mc.metabolites.cpd02438_c0: 1,
                       mc.metabolites.cpd00001_c0: 1})
    mc.add_reactions([d, t])
    return {"growth": growth_row(mc, only("methanol")), "k": kvalue(mc, only("methanol")),
            "catabolic": catabolic(mc, "methanol")}


def main():
    logging.disable(logging.WARNING)
    cobra.Configuration().solver = "gurobi"
    src_md5 = md5(SRC)
    assert src_md5 == SRC_MD5, f"arm-A model changed: {src_md5}"
    arm_a = cobra.io.load_json_model(str(SRC))
    m = cobra.io.load_json_model(str(SRC))
    m, log = curate(m)
    m.id = f"{MAG}__curated"
    cobra.io.save_json_model(m, str(OUT))
    res = {"source_model": str(SRC.relative_to(ROOT)), "source_md5": src_md5, "code_md5": md5(__file__),
           "curated_model": str(OUT.relative_to(ROOT)), "curated_md5": md5(OUT), "rnf": {"Na": RNF_NA, "H": RNF_H},
           "n_reactions": [len(arm_a.reactions), len(m.reactions)], "n_genes": [len(arm_a.genes), len(m.genes)],
           "changes": log, "validation_curated": validate(m), "validation_arm_A": validate(arm_a, full=True)}
    sens = {}
    for na, h in ((3, 2), (3, 0), (2, 0), (1, 0)):
        mm = cobra.io.load_json_model(str(SRC))
        mm, _ = curate(mm, na, h)
        sens[f"Rnf_{na}Na_{h}H"] = {"free_atp": bm.free_atp(mm), "k_MMA": kvalue(mm, only("MMA")),
                                    "growth_MMA": growth_row(mm, only("MMA")),
                                    "catabolic_MMA": catabolic(mm, "MMA"),
                                    "catabolic_MMA_CH4_blocked": catabolic(mm, "MMA", True)}
        rev = mm.copy()
        rev.reactions.rxa45615_c0.lower_bound = -1000
        sens[f"Rnf_{na}Na_{h}H"]["reversible_H2_CO2"] = growth_row(rev, H2CO2)
        sens[f"Rnf_{na}Na_{h}H"]["reversible_DMA_H2"] = growth_row(rev, bm.MANUSCRIPT_METHANOGEN_MEDIUM)
    res["rnf_sensitivity"] = sens
    # coupling sensitivity (Fpo H+, HdrDE H+, Vht H+, Mtr Na+, Rnf Na+ per 2 e- or per methyl)
    coup = {"used (2, 2, 2, 2, 3)": (2, 2, 2, 2, 3),
            "iAF692 (1.8, 1.8, 1.8, 1.7) + iMB745 Rnf 3": (1.8, 1.8, 1.8, 1.7, 3),
            "iMAC868-like Mtr 2, Rnf 4": (2, 2, 2, 2, 4),
            "half coupling (vesicle H+/2e- 0.9; all x0.5)": (1, 1, 1, 1, 1.5)}
    cs = {}
    for lab, c in coup.items():
        mm = cobra.io.load_json_model(str(OUT))
        set_coupling(mm, *c)
        cs[lab] = {"couplings": c, "free_atp": bm.free_atp(mm), "k_MMA": kvalue(mm, only("MMA")),
                   "catabolic_MMA": catabolic(mm, "MMA"), "catabolic_MMA_CH4_blocked": catabolic(mm, "MMA", True),
                   "growth_MMA": growth_row(mm, only("MMA"))}
    res["coupling_sensitivity"] = cs
    res["methanol_what_if"] = methanol_what_if(m)
    VAL.write_text(json.dumps(res, indent=1, default=float))
    print("wrote", OUT, md5(OUT))
    print("wrote", VAL)


if __name__ == "__main__":
    sys.exit(main())
