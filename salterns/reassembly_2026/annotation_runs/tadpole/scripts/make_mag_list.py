#!/usr/bin/env python
"""List the original MAGs for the tadpole-pilot genome check (run locally, ~/Documents/py_venv/bin/python -I).

Writes genome_check_mags.tsv next to this script, one row per MAG:
  mag, population, source_sample, file (in salterns/MAG_contigs/), map_to (comma-separated pilot samples),
  plus genus / completeness / contamination / n_contigs / size_bp / n50 from comparison/mags_old_vs_new.tsv where known.

Populations:
  derep310     original dereplicated MAGs binned from a pilot sample's own (iii) assembly -> mapped to that sample
  nonderep140  the non-dereplicated original MAGs from the same samples (also (iii) bins) -> mapped to that sample
  lost30       the 30 medium/high-quality original MAGs whose species the 2026 re-assembly did not recover
               (comparison/README.md: no 2026 MAG above 80 % ANI); none is from a pilot sample, so each is mapped to
               both pilot samples to see whether its genome is present there at all.
The lost30 selection repeats compare_vs_original.py: species cluster with original MAGs only, MIMAG not 'below MQ',
best skani ANI to any 2026 MAG < 80 (0 when there is no edge).
"""
import csv
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
R26 = os.path.abspath(os.path.join(HERE, "..", "..", ".."))           # salterns/reassembly_2026
S = os.path.dirname(R26)                                              # salterns
PILOT = ["R1_A_D1", "SF2_C_D1"]


def sample_of(name):
    # Salt_Pond_MetaG_R1_A_D1_MG_... / Salt_Pond_MetaGSF2_C_D1_MG_...
    core = name.split("Salt_Pond_MetaG", 1)[1].lstrip("_").split("_MG_")[0]
    return core


def main():
    idx = pd.read_csv(os.path.join(S, "MAG_contigs", "index.tsv"), sep="\t", keep_default_na=False)
    d = pd.read_csv(os.path.join(R26, "comparison", "mags_old_vs_new.tsv"), sep="\t", keep_default_na=False)
    meta = {g.split("|", 1)[1]: r for g, r in zip(d.genome, d.to_dict("records")) if g.startswith("old|")}

    # lost30, as in compare_vs_original.py
    old = d[d.genome.str.startswith("old|")]
    kind = {c: "".join(sorted({x.split("|")[0] for x in g.genome}, reverse=True)) for c, g in d.groupby("species_cluster")}
    lost = old[old.species_cluster.map(kind) == "old"]
    e = pd.read_csv(os.path.join(R26, "comparison", "skani_edges.tsv"), sep="\t")

    def gname(p):
        b = os.path.basename(p)[:-4]
        return ("old|" if "/MAG_contigs/" in p else "new|") + b

    L, near = set(lost.genome), {}
    for a, b, ani in zip(e.Ref_file.map(gname), e.Query_file.map(gname), e.ANI):
        for x, y in ((a, b), (b, a)):
            if x in L and y.startswith("new|"):
                near[x] = max(near.get(x, 0), ani)
    good = lost[lost.mimag != "below MQ"]
    lost30 = [g.split("|", 1)[1] for g in good.genome if near.get(g, 0) < 80]
    assert len(lost30) == 30, len(lost30)

    rows = []
    for r in idx.itertuples():
        smp = sample_of(r.mag)
        if smp in PILOT:
            rows.append(dict(mag=r.mag, population=r.population, source_sample=smp, file=r.file, map_to=smp))
    for m in lost30:
        r = idx[idx.mag == m].iloc[0]
        smp = sample_of(m)
        assert smp not in PILOT
        rows.append(dict(mag=m, population="lost30", source_sample=smp, file=r.file, map_to=",".join(PILOT)))
    for row in rows:
        mr = meta.get(row["mag"], {})
        for k in ("genus", "completeness", "contamination", "n_contigs", "size_bp", "n50"):
            row[k] = mr.get(k, "")
        assert os.path.exists(os.path.join(S, "MAG_contigs", row["file"])), row["file"]
    out = os.path.join(HERE, "genome_check_mags.tsv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("wrote %s: %d MAGs (%s)" % (out, len(rows), pd.Series([r["population"] + ":" + r["source_sample"]
                                                                 if r["population"] != "lost30" else "lost30"
                                                                 for r in rows]).value_counts().to_dict()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
