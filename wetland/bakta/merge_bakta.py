#!/usr/bin/env python3
"""Merge the full Bakta 1.11.4 annotations of the 1,184 DRAM2 representatives (run_bakta_reps.sh) and compare them
with DRAM2 on the same MAGs (per Andrew, 2026-11-18).

Input: <runs dir>/<mag>/<mag>.tsv (Bakta's feature table: Sequence Id, Type, Start, Stop, Strand, Locus Tag, Gene,
Product, DbXrefs), one per MAG with a DONE marker.
Outputs (results/):
  bakta_genes_all.tsv.gz      every feature of every MAG, coding and non-coding (cds, sorf, tRNA, tmRNA, rRNA, ncRNA,
                              ncRNA-region, crispr, oriC/oriV/oriT, gap) + parsed KO / EC / GO / COG columns (local
                              only if > 95 MB)
  bakta_per_mag.tsv           per MAG feature counts: CDS (hypothetical / with KO / EC / GO / COG / Pfam / UniRef90 hit /
                              UniRef50 hit only), pseudogenes, sORF,
                              tRNA, tmRNA, rRNA 5S/16S/23S, ncRNA, ncRNA-region, CRISPR, oriC
  bakta_vs_dram2_per_mag.tsv  per MAG: Bakta vs DRAM2 CDS, tRNA, rRNA by kind, KO-annotated CDS (DRAM2: KOfam hit
                              passing the KO threshold)
  summary.json                totals and medians
Note: Bakta attaches KEGG/EC/GO/COG xrefs at the UniRef90 (PSC) level only; a CDS whose best hit is a UniRef50
cluster (PSCC) gets a product name (and Pfam domains) but no KO. Divergent environmental MAGs therefore show low KO
coverage here; DRAM2's KOfam HMMs are the better KO source for them.
Usage: ~/Documents/py_venv/bin/python -I merge_bakta.py <runs dir>
"""
import glob, json, os, sys
import pandas as pd

RUNS = sys.argv[1]
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)
reps = [l.strip() for l in open(os.path.join(HERE, "reps.txt")) if l.strip()]
COLS = ["contig", "type", "start", "stop", "strand", "locus_tag", "gene", "product", "dbxrefs"]

parts, missing = [], []
for m in reps:
    f = os.path.join(RUNS, m, f"{m}.tsv")
    if not (os.path.exists(os.path.join(RUNS, m, "DONE")) and os.path.exists(f)):
        missing.append(m)
        continue
    d = pd.read_csv(f, sep="\t", comment="#", header=None, names=COLS, dtype=str, keep_default_na=False)
    d.insert(0, "mag", m)
    parts.append(d)
g = pd.concat(parts, ignore_index=True)
x = g.dbxrefs.str.split(r",\s*")


def pick(prefix):
    return x.map(lambda v: ",".join(sorted({t[len(prefix):] for t in v if t.startswith(prefix)})))


g["KO"], g["EC"], g["GO"] = pick("KEGG:"), pick("EC:"), pick("GO:")
g["COG"] = x.map(lambda v: ",".join(sorted({t[4:] for t in v if t.startswith("COG:COG")})))
g["PFAM"] = x.map(lambda v: ",".join(sorted({t[5:] for t in v if t.startswith("PFAM:")})))
g["UniRef"] = x.map(lambda v: ",".join(sorted({t[7:] for t in v if t.startswith("UniRef:")})))
g["pseudo"] = g["product"].str.contains(r"\(pseudo\)", regex=True) | g["product"].str.startswith("(pseudo)")
out_all = os.path.join(OUT, "bakta_genes_all.tsv.gz")
g.to_csv(out_all, sep="\t", index=False, compression="gzip")

cds = g[g.type == "cds"]
by = lambda df: df.groupby("mag").size()
per = pd.DataFrame({
    "CDS": by(cds), "CDS_hypothetical": by(cds[cds["product"].str.contains("hypothetical protein", case=False)]),
    "CDS_with_KO": by(cds[cds.KO != ""]), "CDS_with_EC": by(cds[cds.EC != ""]), "CDS_with_GO": by(cds[cds.GO != ""]),
    "CDS_with_COG": by(cds[cds.COG != ""]), "CDS_with_PFAM": by(cds[cds.PFAM != ""]),
    "CDS_with_UniRef90": by(cds[cds.UniRef.str.contains("UniRef90_")]),
    "CDS_with_UniRef50_only": by(cds[cds.UniRef.str.contains("UniRef50_") & ~cds.UniRef.str.contains("UniRef90_")]),
    "pseudogenes": by(cds[cds.pseudo]),
    "sORF": by(g[g.type == "sorf"]), "tRNA": by(g[g.type == "tRNA"]), "tmRNA": by(g[g.type == "tmRNA"]),
    "ncRNA": by(g[g.type == "ncRNA"]), "ncRNA_region": by(g[g.type == "ncRNA-region"]),
    "CRISPR": by(g[g.type == "crispr"]), "oriC": by(g[g.type == "oriC"]),
})
r = g[g.type == "rRNA"]
for kind in ("5S", "16S", "23S"):
    per[f"rRNA_{kind}"] = by(r[r["product"].str.startswith(kind)])
per = per.reindex([m for m in reps if m not in missing]).fillna(0).astype(int)
per.index.name = "mag"
per.to_csv(os.path.join(OUT, "bakta_per_mag.tsv"), sep="\t")

dram = pd.read_csv(os.path.join(HERE, "..", "dram2", "results", "dram2_per_mag_genes.tsv"), sep="\t").set_index("mag")
cmp_ = pd.DataFrame({
    "bakta_CDS": per.CDS, "dram2_CDS": dram.CDS, "bakta_tRNA": per.tRNA, "dram2_tRNA": dram.tRNA,
    "bakta_rRNA_5S": per.rRNA_5S, "dram2_rRNA_5S": dram.rRNA_5S, "bakta_rRNA_16S": per.rRNA_16S,
    "dram2_rRNA_16S": dram.rRNA_16S, "bakta_rRNA_23S": per.rRNA_23S, "dram2_rRNA_23S": dram.rRNA_23S,
    "bakta_CDS_with_KO": per.CDS_with_KO, "dram2_CDS_KO_passing_threshold": dram.CDS_KO_passing_threshold,
}).loc[per.index]
cmp_.to_csv(os.path.join(OUT, "bakta_vs_dram2_per_mag.tsv"), sep="\t")
tot = per.sum()
summ = {
    "mags": len(per), "missing_or_failed": missing, "features": int(len(g)),
    "feature_types": {k: int(v) for k, v in g.type.value_counts().items()},
    "totals": {k: int(v) for k, v in tot.items()},
    "mags_with_at_least_one": {k: int((per[k] > 0).sum()) for k in per.columns},
    "share_CDS_hypothetical": round(tot.CDS_hypothetical / tot.CDS, 4),
    "bakta_vs_dram2_totals": {k: int(v) for k, v in cmp_.sum().items()},
    "median_CDS_ratio_bakta_over_dram2": round(float((cmp_.bakta_CDS / cmp_.dram2_CDS).median()), 4),
    "bakta_all_genes_bytes": os.path.getsize(out_all),
}
json.dump(summ, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in summ.items() if k != "missing_or_failed"}, indent=1),
      f"\nmissing/failed: {len(missing)}")
