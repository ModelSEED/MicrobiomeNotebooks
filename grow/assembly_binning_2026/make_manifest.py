#!/usr/bin/env python
"""Write manifest.tsv (one row per SRA run, 171 runs / 133 samples) from ../raw_metagenomes_sra/runs.tsv."""
import csv, os, re
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
rows = list(csv.DictReader(open(os.path.join(HERE, "..", "raw_metagenomes_sra", "runs.tsv")), delimiter="\t"))
by = defaultdict(list)
for r in rows:
    by[r["sample"]].append(r)
ok = re.compile(r"^[A-Za-z0-9|._-]+$")
assert all(ok.match(s) for s in by), "sample names must be Workspace object-name safe"
cols = ["sample", "run", "source", "bioproject", "biosample", "spots", "bases", "sample_n_runs", "sample_bases",
        "n_mags_from_sample"]
with open(os.path.join(HERE, "manifest.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(cols)
    for s in sorted(by):
        tot = sum(int(r["bases"]) for r in by[s])
        for r in sorted(by[s], key=lambda r: r["Run"]):
            w.writerow([s, r["Run"], r["source"], r["BioProject"], r["BioSample"], r["spots"], r["bases"], len(by[s]), tot,
                        r["n_mags_from_sample"]])
print(len(rows), "runs,", len(by), "samples,", sum(len(v) > 1 for v in by.values()), "multi-run")
