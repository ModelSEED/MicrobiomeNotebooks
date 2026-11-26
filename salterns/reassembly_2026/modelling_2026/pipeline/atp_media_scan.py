#!/usr/bin/env python
"""Figure S10 data for the 2026 guild models: the core (ATP) gap-filling each arm-A draft needs on every ATP test
medium, i.e. MSATPCorrection.evaluate_growth_media() on the draft exactly as build_model.py makes it before ATP
correction (RAST draft, ATPM, and for MePn consumers the cell-5 disables and the MePn module). Per medium: 0 = the
draft makes ATP without gap-filling, n = reactions added or reversed by the core gap-fill, null = no solution.
This is what the original Figure S10 ('Core Gap Filling - MePn consumers / BET consumers') showed for the 2020 models.

Run (TF, ~/tf_runs/salterns_modelling, same inputs and code as the build):
  ~/venvs/salterns_modelling/bin/python code/atp_media_scan.py --root inputs --mag <mag> --guild <guild> --out <dir>
"""
import argparse
import gzip
import json
import logging
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_model as B  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--mag", required=True)
    ap.add_argument("--guild", required=True, choices=["mepn_consumer", "betaine_reducer"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--msdb", default=os.environ.get("MODELSEED_DB_PATH"))
    ap.add_argument("--solver", default="gurobi")
    a = ap.parse_args()
    logging.basicConfig(level=logging.WARNING)
    import cobra
    import pandas as pd
    cobra.Configuration().solver = a.solver
    root, out = Path(a.root), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    sp = pd.read_csv(root / "modelling_species.tsv", sep="\t")
    sp["mag_id"] = sp.rep_genome.str.split("|").str[1]
    row = sp.set_index("mag_id").loc[a.mag]
    tfile, _ = B.template_for(row.domain, row.phylum)
    with gzip.open(root / "evidence" / f"{a.mag}.json.gz", "rt") as fh:
        ev = json.load(fh)
    from modelseedpy.biochem import modelseed_biochem as MB
    MB.ModelSEEDBiochem.default_biochemistry = MB.from_local(a.msdb)
    gs_t = B.load_template(root / "templates" / tfile)
    core_t = B.load_template(root / "templates/Core-V6.json")
    from modelseedpy.core.msgenome import MSGenome, MSFeature
    feats = []
    for g, roles in ev["rast"].items():
        f = MSFeature(g, "")
        for r in roles:
            f.add_ontology_term("RAST", r)
        feats.append(f)
    genome = MSGenome()
    genome.id = a.mag
    genome.add_features(feats)
    extra = B.arm_evidence(ev, "A", gs_t)
    model, _ = B.build_draft(genome, gs_t, a.mag, extra)
    model.id = a.mag
    from modelseedpy.core.msbuilder import MSBuilder
    MSBuilder.add_atpm(model)
    prot = {}
    if a.guild == "mepn_consumer":
        B.apply_cell5(model, prot)
        B.mepn_template_additions(gs_t, prot)
        B.add_mepn_module(model, gs_t, ev, "A", prot)
    from modelseedpy.core.msatpcorrection import MSATPCorrection, load_default_medias
    atp = MSATPCorrection(model, core_t, load_default_medias(), compartment="c0", atp_hydrolysis_id="ATPM_c0",
                          load_default_medias=False)
    atp.evaluate_growth_media()
    res = {}
    for media, sol in atp.media_gapfill_stats.items():
        res[media.id] = None if sol is None else len(sol.get("new", {})) + len(sol.get("reversed", {}))
    rec = {"mag_id": a.mag, "guild": a.guild, "template": tfile, "core_gapfill_by_medium": res,
           "seconds": round(time.time() - t0, 1), "build_model_md5": B.md5(B.__file__)}
    json.dump(rec, open(out / f"{a.mag}__{a.guild}.json", "w"), indent=1)
    print(f"{a.mag} {a.guild} {len(res)} media, {sum(v is None for v in res.values())} without solution, "
          f"{rec['seconds']} s")


if __name__ == "__main__":
    main()
