#!/usr/bin/env python3
"""Assign MAG iterativeIDs to the wetland MAGs (per Andrew, 2026-11-15), following the MAG rule used in EmilyKin
(scripts/recompute_mag_ids.py, the GTDB branch; the wetland study has no ASVs, so there is no ASV inheritance):

  root = the GTDB R232 genus if it is named and contains no digits; otherwise the deepest R232 rank (domain..genus)
         whose label contains no digits (placeholder names such as JABFXC01 are skipped); 'Bacteria'/'Archaea' if none.
  id   = f"{root}.{N}_m", N = 1, 2, ... within each root.

Ordering of N (deterministic): within a root, MAGs by descending quality score (CheckM2 completeness - 5 x contamination),
then MAG name, so '.1' is the best MAG of that root. IDs are FROZEN once written: on a re-run, MAGs already in
mag_iterativeIDs.tsv keep their id and new MAGs (e.g. from re-assembled libraries) take the next free N of their root.
A MAG's id is not changed by later taxonomy updates unless --reassign is given (then the old->new map is written).
Inputs: census outputs (mag_census.tsv). Output: mag_iterativeIDs.tsv (mag, iterativeID, root, root_rank, kbase_ref,
species_cluster, r232_classification), plus the iterativeID column added to mag_census.tsv, rep_worklist.tsv and
dram2_hq_multisample.tsv in the same directory.
Usage: python -I iterative_ids.py <census outputs dir> [--reassign]
"""
import os, re, sys
import pandas as pd

D = sys.argv[1]
REASSIGN = "--reassign" in sys.argv
RANKS = [("d", "Domain"), ("p", "Phylum"), ("c", "Class"), ("o", "Order"), ("f", "Family"), ("g", "Genus")]
HASDIGIT = re.compile(r"\d")
OUT = os.path.join(D, "mag_iterativeIDs.tsv")

m = pd.read_csv(os.path.join(D, "mag_census.tsv"), sep="\t")


def root_of(lineage):
    d = dict(x.split("__", 1) for x in str(lineage).split(";") if "__" in x)
    g = d.get("g", "")
    if g and not HASDIGIT.search(g):
        return g, "Genus"
    best = ("", "")
    for k, name in RANKS:
        v = d.get(k, "")
        if v and not HASDIGIT.search(v):
            best = (v, name)
    return best if best[0] else (d.get("d", "Bacteria") or "Bacteria", "Domain")


m["root"], m["root_rank"] = zip(*m.r232_classification.map(root_of))
old = pd.read_csv(OUT, sep="\t") if os.path.exists(OUT) and not REASSIGN else None
ids = dict(zip(old.mag, old.iterativeID)) if old is not None else {}
used = {}
for i in ids.values():
    r, n = i[:-2].rsplit(".", 1)
    used.setdefault(r, set()).add(int(n))
for root, g in m.sort_values(["root", "score", "mag"], ascending=[True, False, True]).groupby("root", sort=False):
    n = 0
    taken = used.get(root, set())
    for mag in g.mag:
        if mag in ids:
            continue
        n += 1
        while n in taken:
            n += 1
        ids[mag] = f"{root}.{n}_m"
        taken.add(n)
    used[root] = taken
m["iterativeID"] = m.mag.map(ids)
assert m.iterativeID.notna().all() and m.iterativeID.is_unique, "iterativeIDs must be unique and complete"
m[["mag", "iterativeID", "root", "root_rank", "kbase_ref", "species_cluster", "r232_classification"]].to_csv(
    OUT, sep="\t", index=False)
for f in ("mag_census.tsv", "rep_worklist.tsv", "dram2_hq_multisample.tsv"):
    p = os.path.join(D, f)
    if not os.path.exists(p):
        continue
    t = pd.read_csv(p, sep="\t")
    t = t.drop(columns=[c for c in ("iterativeID",) if c in t.columns])
    t.insert(1, "iterativeID", t.mag.map(ids))
    t.to_csv(p, sep="\t", index=False)
print(f"{len(m)} MAGs, {m.root.nunique()} roots; root rank: {m.root_rank.value_counts().to_dict()}")
print("examples:", m.sort_values("score", ascending=False).head(5)[["mag", "iterativeID"]].values.tolist())
print("largest roots:", m.root.value_counts().head(8).to_dict())
