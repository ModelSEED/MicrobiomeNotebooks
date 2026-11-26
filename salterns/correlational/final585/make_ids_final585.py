#!/usr/bin/env python
"""IterativeIDs and phylogeny for the 585 species representatives of the final re-assembled MAG set (2026-11-19):
935 second-pass (Arm C) MAGs + 520 kept first-pass site co-assembly MAGs, joint skani species clusters
(../../reassembly_2026/comparison_rollout/mags_three_way.tsv, final_species_rep == True).

Convention as in ../make_iterativeIDs_r232.py (analysis_2_marimo.py "IterativeIDs" cell): walk Domain -> Genus, stop at
the first empty rank (species is never used), name the representative "<deepest named taxon>.<n>". Numbering:
  * a species cluster that contains one of the 310 original MAGs keeps that MAG's R232 ID when the representative's
    deepest taxon is the same taxon (several originals in a cluster: the lowest number), so names cited in the
    manuscript keep pointing at the same organism;
  * every other cluster gets the next free number of its taxon, in order of decreasing mean relative abundance over the
    24 sediment samples (../../reassembly_2026/abundance_final/final585_relabund_24sediment.tsv).
Numbers of lost originals are therefore not reused.
Keys are the abundance-table row names (the representative MAG name without the 'rollout|' / 'coassembly|' prefix).
Writes, in this folder:
  phylogeny_final585.json            MAG -> {Kingdom..Species} (the representative's GTDB R232 lineage)
  iterativeIDs_final585.json         MAG -> iterativeID
  iterativeID_levels_final585.json   iterativeID -> rank used
  iterativeIDs_final585.csv          per species: representative, set, ID, inherited ID, originals in the cluster
Run: ~/Documents/py_venv/bin/python make_ids_final585.py
"""
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
COR = HERE.parent
S = COR.parent
RA = S / "reassembly_2026"
LEVELS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
COLS = ["domain", "phylum", "cls", "order", "family", "genus", "species"]

d = pd.read_csv(RA / "comparison_rollout/mags_three_way.tsv", sep="\t")
ab = pd.read_csv(RA / "abundance_final/final585_relabund_24sediment.tsv", sep="\t", index_col=0)
reps = d[d.final_species_rep == True].copy().rename(columns={"class": "cls"})
reps["name"] = reps.genome.str.split("|", n=1).str[1]
assert len(reps) == 585 and reps.name.is_unique and set(reps.name) <= set(ab.index)
reps["mean_rel"] = ab.loc[reps.name].mean(axis=1).values


def lineage(r):
    return {lv: ("" if pd.isna(getattr(r, c)) else str(getattr(r, c))) for lv, c in zip(LEVELS, COLS)}


def deepest(t):
    best = "Kingdom"
    for lv in LEVELS[:-1]:
        if not t[lv]:
            break
        best = lv
    return best


phylo, level_of, taxon_of = {}, {}, {}
for r in reps.itertuples():
    t = lineage(r)
    phylo[r.name] = t
    level_of[r.name] = deepest(t)
    taxon_of[r.name] = t[level_of[r.name]]

# originals per cluster, with their R232 IDs (key() join as in final_set_guilds.py)
ids_r232 = json.load(open(COR / "iterativeIDs_r232.json"))
o = d[d.set == "original"].copy()
o["mag"] = o.genome.str.split("|", n=1).str[1].str.replace(".contigs__", "", regex=False)
assert o.mag.isin(ids_r232).all()
o["id_r232"] = o.mag.map(ids_r232)
orig_ids = o.groupby("species_cluster").id_r232.apply(lambda s: sorted(s, key=lambda i: (i.rsplit(".", 1)[0], int(i.rsplit(".", 1)[1]))))

ids, inherited = {}, {}
used = {}                                    # taxon -> set of numbers taken (inherited + lost originals' numbers)
for i in ids_r232.values():
    t, n = i.rsplit(".", 1)
    used.setdefault(t, set()).add(int(n))    # never reuse any 2020/R232 number, recovered or not
for r in reps.itertuples():
    cands = [i for i in orig_ids.get(r.species_cluster, []) if i.rsplit(".", 1)[0] == taxon_of[r.name]]
    if cands:
        ids[r.name] = inherited[r.name] = cands[0]
for r in reps.sort_values("mean_rel", ascending=False).itertuples():
    if r.name in ids:
        continue
    t = taxon_of[r.name]
    n = max(used.get(t, {0})) + 1
    used.setdefault(t, set()).add(n)
    ids[r.name] = f"{t}.{n}"
assert len(set(ids.values())) == 585

json.dump(phylo, open(HERE / "phylogeny_final585.json", "w"), indent=1)
json.dump(ids, open(HERE / "iterativeIDs_final585.json", "w"))
json.dump({ids[m]: level_of[m] for m in ids}, open(HERE / "iterativeID_levels_final585.json", "w"))
rows = []
for r in reps.itertuples():
    oi = orig_ids.get(r.species_cluster, [])
    rows.append({"mag": r.name, "set": r.set, "species_cluster": r.species_cluster, "iterativeID": ids[r.name],
                 "level": level_of[r.name], "inherited_from_r232": inherited.get(r.name, ""),
                 "originals_in_cluster_r232": " ".join(oi), "n_originals": len(oi),
                 "mean_rel_abund_24": r.mean_rel, "gtdb_r232": r.classification})
t = pd.DataFrame(rows).sort_values("iterativeID")
t.to_csv(HERE / "iterativeIDs_final585.csv", index=False)
print(f"585 iterativeIDs; levels {pd.Series(level_of).value_counts().to_dict()}")
print(f"clusters with an original: {int((t.n_originals > 0).sum())}; ID inherited: {int((t.inherited_from_r232 != '').sum())}; "
      f"with an original but a different deepest taxon (new ID): {int(((t.n_originals > 0) & (t.inherited_from_r232 == '')).sum())}")
print(f"R232 IDs whose cluster is absent from the final set: {len(set(ids_r232.values()) - set(o[o.species_cluster.isin(reps.species_cluster)].id_r232))}")
