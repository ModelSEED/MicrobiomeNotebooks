"""Assemble the final model set from the repaired winning-arm models.

Re-runnable. Reads ../models/_repair_runs/<arm>/<tag>.{json,report.json} (repair_model.py output, copied back from
TF) and the comparison build reports in ../annotation_comparison/runs/<arm>/, writes
  ../models/<guild>/<mag_id>.json                   kept models (cobra JSON)
  ../models/reports/<guild>/<mag_id>.report.json    build report + repair report of each kept model
  ../models/dropped/<guild>/<mag_id>.{json,report.json}   dropped models, kept for inspection
  ../models/models_summary.tsv                       one row per (species, guild) build: QC, gap-fill lists, kept/why
  ../models/MANIFEST.md5
Drop rule: no growth on the guild medium after gap-filling and repair (as the original dropped its 4 non-growing
models), or ATP from a carbon-free medium. Guild-function flags (MePn obligacy, no TMA export, non-methanogenic
growth) are reported, not dropped on.
usage: python finalize.py [--arm A]
"""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
M26 = HERE.parent
AC = M26 / "annotation_comparison"
MOD = M26 / "models"


def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="A")
    arm = ap.parse_args().arm
    sp = pd.read_csv(M26.parent / "comparison_rollout/modelling_species.tsv", sep="\t")
    sp["mag_id"] = sp.rep_genome.str.split("|").str[1]
    sp = sp.set_index("mag_id")
    for d in ("mepn_consumer", "betaine_reducer", "methanogen", "reports", "dropped"):
        if (MOD / d).exists():
            shutil.rmtree(MOD / d)
    rows, man = [], []
    for p in sorted((MOD / "_repair_runs" / arm).glob("*.report.json")):
        r = json.load(open(p))
        b = json.load(open(AC / "runs" / arm / p.name))           # comparison build report
        g, mag = r["guild"], r["mag_id"]
        row = {"guild": g, "mag_id": mag, "species_cluster": sp.loc[mag, "species_cluster"],
               "domain": sp.loc[mag, "domain"], "phylum": sp.loc[mag, "phylum"], "family": sp.loc[mag, "family"],
               "genus": sp.loc[mag, "genus"], "arm": arm, "repair_status": r["status"]}
        if r["status"] != "ok":
            row.update(kept=False, drop_reason=f"repair error: {r.get('error')}")
            rows.append(row)
            continue
        q, rp, gf = r["qc"], r["repair"], r["gapfilled"]
        stages = {}
        for v in gf.values():
            stages[v] = stages.get(v, 0) + 1
        row.update({
            "template": r["template_file"], "template_rule": r["template_rule"],
            "n_reactions": r["final"]["n_reactions"], "n_metabolites": r["final"]["n_metabolites"],
            "n_genes": r["final"]["n_genes"],
            "n_gene_associated_non_gapfilled": r["final"]["n_gene_associated_non_gapfilled"],
            "n_gapfilled": len(gf), "n_gapfilled_by_stage": json.dumps(stages, sort_keys=True),
            "atp_core_gapfill": len(b["atp_correction"]["core_gapfilling"]),
            "hand_added": ";".join(r.get("hand_added", [])),
            "repair_nh3_oxidation_blocked": ";".join(f"{k}({v['gapfilled'] or 'draft'})" for k, v in rp["nh3_oxidation_block"].items()),
            "repair_regapfill_added": ";".join(x for x, s in gf.items() if s.startswith("repair_")),
            "methane_consuming_reactions": ";".join(f"{k}:{v['action']}" for k, v in rp["methane_consumption"]["model"].items()),
            "free_atp_before_repair": json.dumps(rp["free_atp_before"]),
            "guild_growth": q["guild_growth"], "grows": q["grows"],
            "free_atp_pass": q["free_atp_pass"], "free_atp": json.dumps(q["free_atp"]),
            "atp_tests_passed": f'{sum(t["pass"] for t in q["atp_tests"])}/{len(q["atp_tests"])}',
            "atp_tests_failed": ";".join(t["media"] for t in q["atp_tests"] if not t["pass"]),
            "blocked_frac": round(q["blocked_frac"], 4),
            "dead_end_metabolite_frac": round(q["dead_ends"]["dead_end_metabolite_frac"], 4),
            "guild_pathway_fraction_template_mapped": round(r["guild_pathway"]["fraction"], 3),
            "guild_pathway_gene_supported": json.dumps(r["guild_pathway"]["gene_supported"], sort_keys=True),
            "gapfilled_reactions": ";".join(f"{k}:{v}" for k, v in sorted(gf.items())),
        })
        if g == "mepn_consumer":
            t = q["growth_triplet"]
            row.update(growth_sim_medium=q["growth_sim_medium"][0], growth_no_MePn=q["growth_no_MePn"][0],
                       triplet_a_as_is=t["a_as_is"], triplet_b_MePn_blocked=t["b_MePn_blocked"],
                       triplet_c_Pi_blocked=t["c_Pi_blocked_MePn_only_P"], mepn_obligate=t["obligate"],
                       CH4_at_max_growth_PiFree=q.get("CH4_at_max_growth_PiFree"),
                       rxs00001_gpr=b["protocol"].get("rxs00001_gpr"), rxs43120_gpr=b["protocol"].get("rxs43120_gpr"),
                       cp_lyase_gene_supported=all(r["guild_pathway"]["gene_supported"].values()))
            fl = []
            if t["obligate"]:
                fl.append("MePn-obligate")
            if not (q.get("CH4_at_max_growth_PiFree") or 0) > 1e-6:
                fl.append("no CH4 export at max growth on the Pi-free MePn medium")
            flag = "; ".join(fl)
        elif g == "betaine_reducer":
            row.update(growth_lactate_aerobic=q["growth_lactate_aerobic"][0], growth_no_O2=q["growth_no_O2"][0],
                       TMA_export_max_at_bio_0_01=q["TMA_export_max_at_bio_0.01"],
                       rxn17220_gpr=rp.get("rxn17220_gpr_after"),
                       betaine_reductase_gene_supported=("annotation-supported (RAST roles), not template-mapped"
                                                         if rp.get("rxn17220_gpr_after") else "no"))
            flag = "" if (q["TMA_export_max_at_bio_0.01"] or 0) > 1e-6 else "no TMA export"
        else:
            fpo = rp.get("methanogen_fpo_test", {})
            row.update(CH4_at_max_growth=q.get("CH4_at_max_growth"), CH4_FVA_99=json.dumps(q.get("CH4_FVA_at_99pct_growth")),
                       CH4_max_biomass_free=q.get("CH4_max_biomass_free"),
                       growth_acetate_export_blocked=q.get("growth_acetate_export_blocked"),
                       growth_manuscript_DMA_H2=q["growth_manuscript_DMA_H2_medium"][0],
                       CH4_max_manuscript_medium=q.get("CH4_max_manuscript_medium"),
                       growth_single_substrate=json.dumps({k[12:-3]: v for k, v in q.items() if k.startswith("growth_only_")}),
                       fpo_freed_test=json.dumps(fpo))
            meth = (q.get("growth_acetate_export_blocked") or 0) > 1e-6 and (q.get("CH4_at_max_growth") or 0) > 1e-6
            flag = "" if meth else "growth not methanogenic (methyl -> acetate); CH4 possible only apart from growth"
        row["guild_flag"] = flag
        reasons = []
        if not q["grows"]:
            reasons.append("no growth on the guild medium after gap-filling and repair")
        if not q["free_atp_pass"]:
            reasons.append(f"ATP from a carbon-free medium {q['free_atp']}")
        row["kept"] = not reasons
        row["drop_reason"] = "; ".join(reasons)
        src = MOD / "_repair_runs" / arm / f"{r['tag']}.json"
        sub = g if row["kept"] else f"dropped/{g}"
        (MOD / sub).mkdir(parents=True, exist_ok=True)
        dst = MOD / sub / f"{mag}.json"
        shutil.copyfile(src, dst)
        rd = MOD / ("reports/" + g if row["kept"] else sub)
        rd.mkdir(parents=True, exist_ok=True)
        dstr = rd / f"{mag}.report.json"
        json.dump({"build": b, "repair": r}, open(dstr, "w"), indent=1, default=str)
        row["model_path"] = str(dst.relative_to(M26))
        row["model_md5"] = md5(dst)
        man += [f"{md5(dst)}  {dst.relative_to(MOD)}", f"{md5(dstr)}  {dstr.relative_to(MOD)}"]
        rows.append(row)
    df = pd.DataFrame(rows).sort_values(["guild", "kept", "genus", "mag_id"], ascending=[True, False, True, True])
    df.to_csv(MOD / "models_summary.tsv", sep="\t", index=False)
    (MOD / "MANIFEST.md5").write_text("\n".join(sorted(man, key=lambda x: x.split("  ")[1])) + "\n")
    print(df.groupby("guild").agg(n=("mag_id", "size"), kept=("kept", "sum")).to_string())
    print(df[~df.kept][["guild", "mag_id", "genus", "drop_reason"]].to_string())
    print(df[df.guild_flag.fillna("") != ""][["guild", "mag_id", "genus", "guild_flag"]].to_string())


if __name__ == "__main__":
    main()
