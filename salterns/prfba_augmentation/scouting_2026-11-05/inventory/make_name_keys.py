#!/usr/bin/env python
"""One row per dereplicated MAG with every name form used across the salterns files, joined by key() (read-only)."""
import json, os, re, sys
from pathlib import Path
import pandas as pd
S = Path("/home/freiburger/Documents/MicrobiomeNotebooks/salterns")
PAT = re.compile(r"Salt_Pond_MetaG_?(.+?)_MG_DASTool_bins[._](concoct_out|metabat|maxbin)\.(\d+)")
key = lambda n: (lambda m: f"{m[1]}|{m[2]}|{int(m[3])}" if m else None)(PAT.search(str(n)))
g = pd.read_csv(S / "correlational/r232/guilds_r232.csv")
ids20 = json.load(open(S / "correlational/iterativeIDs.json"))
cm = pd.read_csv(S / "MAG_QC/checkm/summary/CheckM_summary_table.tsv", sep="\t")
cmk = {key(b): b for b in cm["Bin Name"]}
gt = pd.concat([pd.read_csv(S / f"MAG_QC/gtdbtk_local_r232/gtdbtk.{d}.summary.tsv", sep="\t") for d in ("ar53", "bac120")])
gtk = {key(u): u for u in gt.user_genome}
fl = json.load(open(S / "correlational/saltern_fluxes.json"))["0.07"]["methane"]; flk = {key(k): k for k in fl}
fb = json.load(open(S / "correlational/saltern_fluxes_betaine.json"))["TMA"]; fbk = {key(k): k for k in fb}
mep = {key(p.name): p.name for p in (S / "models/MePn_consumers").glob("*.json")}
rows = []
for r in g.itertuples():
    k = key(r.mag)
    base = S / "correlational/models" / f"{r.mag}.contigs__.RAST.json"
    rows.append(dict(
        mag=r.mag, key=k, iterativeID_r232=r.iterativeID_r232, iterativeID_2020=ids20.get(r.mag, ""),
        name_variant=("bins." if "bins." in r.mag else "bins_") + ("|MetaG_" if "MetaG_" in r.mag else "|MetaG") +
                     ("|zero-padded" if re.search(r"\.0\d+$", r.mag) else ""),
        rast_json=f"jplfaria/Cliff MAGs/RASTified_MAG_json/{r.mag}.contigs__.RAST.json",
        rast_genome_id=f"{r.mag}.contigs__.RAST", rast_feature_id_form=f"{r.mag}.contigs__.RAST.CDS.<n>",
        contigs_fna=f"MAG_contigs/{r.mag}.contigs__.fna", contigs_exists=(S / f"MAG_contigs/{r.mag}.contigs__.fna").exists(),
        base_model=str(base.relative_to(S)) if base.exists() else "",
        base_model_id=json.load(open(base))["id"] if base.exists() else "",
        mepn_model=f"models/MePn_consumers/{mep[k]}" if k in mep else "",
        betaine_model=(f"correlational/models/{r.mag}.contigs__.RAST_gapfilled_betaine.json"
                       if (S / f"correlational/models/{r.mag}.contigs__.RAST_gapfilled_betaine.json").exists() else ""),
        mepn_flux_key=flk.get(k, ""), betaine_flux_key=fbk.get(k, ""),
        checkm_bin_name=cmk.get(k, ""), gtdbtk_user_genome=gtk.get(k, ""), phylogeny_key=f"{r.mag}.contigs__",
    ))
df = pd.DataFrame(rows)
df.to_csv(sys.argv[1], index=False)
print(df.name_variant.value_counts().to_string()); print("all exist:", df.contigs_exists.all(), "base models:", (df.base_model != "").sum(),
      "MePn models:", (df.mepn_model != "").sum(), "betaine models:", (df.betaine_model != "").sum(),
      "MePn flux keys:", (df.mepn_flux_key != "").sum(), "betaine flux keys:", (df.betaine_flux_key != "").sum(),
      "base_model_id == rast_genome_id:", (df.base_model_id == df.rast_genome_id).all())
