#!/usr/bin/env python
"""Reference tables for the min-contig DAS Tool experiment.

results/originals.tsv : the 310 original (2020, dereplicated) MAGs with their 2026 CheckM, R232 genus, species cluster
                        (comparison/mags_old_vs_new.tsv), FASTA path (MAG_contigs/index.tsv), the re-assembly unit
                        name of their sample, and two 'lost' flags:
                          lost58 = species cluster has no re-assembly ('new|') member (compare_vs_original.py)
                          lost30 = lost58, MIMAG medium/high, and no re-assembly MAG above 80 % ANI in skani_edges.tsv
results/current_mags.tsv : the re-assembly's current MAGs (CheckM >= 50 / < 10) of the experiment units
runs/<unit>/skani_originals_query.txt : same-sample originals + all lost58 originals (for skani vs the assembly)
Usage: ~/Documents/py_venv/bin/python -I prep_refs.py
"""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)
RA = os.path.dirname(os.path.dirname(EXP))
S = os.path.dirname(RA)
sys.path.insert(0, HERE)
from units import UNITS, unit_samples, orig_to_unit_sample  # noqa: E402


def main():
    d = pd.read_csv(os.path.join(RA, "comparison/mags_old_vs_new.tsv"), sep="\t", keep_default_na=False)
    idx = pd.read_csv(os.path.join(S, "MAG_contigs/index.tsv"), sep="\t").set_index("assembly_name")
    old, new = d[d.genome.str.startswith("old|")].copy(), d[d.genome.str.startswith("new|")].copy()
    has_new = d.groupby("species_cluster").genome.agg(lambda s: s.str.startswith("new|").any())
    old["lost58"] = ~old.species_cluster.map(has_new)
    an = old.genome.str[4:]
    old["fasta"] = [os.path.join(S, "MAG_contigs", idx.loc[n, "file"]) for n in an]
    assert all(os.path.exists(f) for f in old.fasta)
    e = pd.read_csv(os.path.join(RA, "comparison/skani_edges.tsv"), sep="\t")
    fa2g = dict(zip(old.fasta, old.genome))
    for f in pd.read_csv(os.path.join(RA, "comparison/skani_inputs.txt"), header=None)[0]:
        if "/work/checkm/" in f:
            fa2g[f] = "new|" + os.path.basename(f)[:-4]
    assert set(e.Ref_file) | set(e.Query_file) <= set(fa2g), "unmapped skani edge files"
    near = {}
    for a, b, ani in zip(e.Ref_file.map(fa2g), e.Query_file.map(fa2g), e.ANI):
        for x, y in ((a, b), (b, a)):
            if str(x).startswith("old|") and str(y).startswith("new|"):
                near[x] = max(near.get(x, 0), ani)
    old["best_new_ani_any_af30"] = old.genome.map(near).fillna(0)
    old["lost30"] = old.lost58 & (old.mimag != "below MQ") & (old.best_new_ani_any_af30 < 80)
    old["unit_sample"] = old["sample"].map(orig_to_unit_sample)
    cols = ["genome", "sample", "unit_sample", "completeness", "contamination", "mimag", "phylum", "family", "genus",
            "species", "species_cluster", "n_contigs", "size_bp", "n50", "lost58", "lost30", "best_new_ani_any_af30", "fasta"]
    old[cols].to_csv(os.path.join(EXP, "results/originals.tsv"), sep="\t", index=False)
    print("originals", len(old), "lost58", int(old.lost58.sum()), "lost30", int(old.lost30.sum()))

    cur = new[new["sample"].isin(UNITS)].copy()
    cur["fasta"] = [os.path.join(RA, f"work/checkm/{u}/fasta/{g[4:]}.fna") for u, g in zip(cur["sample"], cur.genome)]
    assert all(os.path.exists(f) for f in cur.fasta)
    cur[["genome", "sample", "completeness", "contamination", "mimag", "genus", "species", "species_cluster",
         "n_contigs", "size_bp", "n50", "fasta"]].to_csv(os.path.join(EXP, "results/current_mags.tsv"), sep="\t", index=False)
    print("current MAGs in experiment units", len(cur))

    for u in UNITS:
        os.makedirs(os.path.join(EXP, "runs", u), exist_ok=True)
        same = old[old.unit_sample.isin(unit_samples(u))]
        q = pd.concat([same.fasta, old[old.lost58].fasta]).drop_duplicates()
        q.to_csv(os.path.join(EXP, "runs", u, "skani_originals_query.txt"), index=False, header=False)
        print(u, "same-sample originals", len(same), "MQ+", int((same.mimag != "below MQ").sum()),
              "lost58", int(same.lost58.sum()), "lost30", int(same.lost30.sum()), "query", len(q))


if __name__ == "__main__":
    main()
