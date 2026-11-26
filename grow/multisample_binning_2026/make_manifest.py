#!/usr/bin/env python
"""Write manifest.tsv: the 67 grouped samples of site_groups.tsv (25 site groups, 8 co-assembly pilot groups), with
their SRA runs from ../assembly_binning_2026/manifest.tsv (the per-sample run's manifest, read only)."""
import csv, os, re
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sg = [r for r in csv.DictReader(open(os.path.join(HERE, "site_groups.tsv")), delimiter="\t") if r["site_group"]]
ps = defaultdict(list)
for r in csv.DictReader(open(os.path.join(HERE, "..", "assembly_binning_2026", "manifest.tsv")), delimiter="\t"):
    ps[r["sample"]].append(r["run"])
groups = defaultdict(list)
for r in sg:
    groups[r["site_group"]].append(r)
ok = re.compile(r"^[A-Za-z0-9|._-]+$")
assert len(sg) == 67 and len(groups) == 25, (len(sg), len(groups))
assert sum(1 for g in groups.values() if g[0]["coassembly_pilot"] == "True") == 8
assert all(ok.match(g) for g in groups) and all(r["sample"] in ps for r in sg)
assert all(int(r["group_size"]) == len(groups[r["site_group"]]) for r in sg)
assert all(len({x["coassembly_pilot"] for x in g}) == 1 for g in groups.values())
cols = ["site_group", "sample", "group_size", "coassembly_pilot", "program", "river", "station", "date", "Gbp", "runs"]
with open(os.path.join(HERE, "manifest.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(cols)
    for g in sorted(groups):
        for r in sorted(groups[g], key=lambda r: (r["date"], r["sample"])):
            w.writerow([g, r["sample"], r["group_size"], r["coassembly_pilot"], r["program"], r["river"], r["station"],
                        r["date"], r["Gbp"], ",".join(sorted(ps[r["sample"]]))])
print(len(sg), "samples,", len(groups), "groups,", sum(1 for g in groups.values() if g[0]["coassembly_pilot"] == "True"),
      "co-assembly pilots")
