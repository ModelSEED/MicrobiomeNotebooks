#!/usr/bin/env python
"""Compare the 2026 re-QC of the 310 dereplicated salterns MAGs against the 2020-era tables.

CheckM: MAG_QC/checkm/summary/CheckM_summary_table.tsv (KBase CheckM 1.0.18, reduced tree)
        vs cliff_mags/checkm_summary_table.tsv
GTDB:   MAG_QC/gtdbtk_local_r232/gtdbtk.{bac120,ar53}.summary.tsv (GTDB-Tk 2.8.0, r232, local)
        vs cliff_mags/gtdb_table_{bac,arc}.csv (2020-era GTDB-Tk 1.x output)

Each GTDB rank difference is labelled:
    renamed   the 2020 name no longer exists at that rank anywhere in r232 (renamed or split)
    moved     the 2020 name still exists in r232, but the MAG is now placed elsewhere
    resolved  2020 had no name at that rank, r232 has one
    lost      2020 had a name at that rank, r232 has none
Writes MAG_QC/comparison_vs_2020.csv.  Run: ~/Documents/py_venv/bin/python MAG_QC/compare_vs_2020.py
"""
import collections
import os
import re
from pathlib import Path

import pandas as pd

S = Path(__file__).resolve().parents[1]
QC = S / "MAG_QC"
R232_TAX = Path.home() / "Documents/gtdbtk_data/release232/taxonomy/gtdb_taxonomy.tsv"
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
PAT = re.compile(r"Salt_Pond_MetaG_?(.+?)_MG_DASTool_bins[._](concoct_out|metabat|maxbin)\.(\d+)")


def key(name):
    m = PAT.search(str(name))
    return f"{m[1]}|{m[2]}|{int(m[3])}" if m else None


def split_lineage(lineage):
    parts = dict(p.split("__", 1) for p in str(lineage).split(";") if "__" in p)
    return [parts.get(r[0], "") for r in RANKS]


# ---- CheckM ------------------------------------------------------------------------------
old_cm = pd.read_csv(S / "cliff_mags/checkm_summary_table.tsv", sep="\t")
new_cm = pd.read_csv(QC / "checkm/summary/CheckM_summary_table.tsv", sep="\t")
for d in (old_cm, new_cm):
    d["mag_key"] = d["Bin Name"].map(key)
cm = old_cm.merge(new_cm, on="mag_key", suffixes=("_2020", "_2026"))
for c in ("Completeness", "Contamination"):
    cm[f"d_{c}"] = (cm[f"{c}_2026"] - cm[f"{c}_2020"]).round(2)
for y in ("2020", "2026"):
    cm[f"MIMAG_MQplus_{y}"] = (cm[f"Completeness_{y}"] >= 50) & (cm[f"Contamination_{y}"] < 10)
cm = cm[["mag_key", "Bin Name_2026", "Marker Lineage_2026", "Completeness_2020", "Completeness_2026",
         "d_Completeness", "Contamination_2020", "Contamination_2026", "d_Contamination",
         "MIMAG_MQplus_2020", "MIMAG_MQplus_2026"]].rename(columns={"Bin Name_2026": "assembly_name"})

# ---- GTDB --------------------------------------------------------------------------------
old = pd.concat([pd.read_csv(S / f"cliff_mags/gtdb_table_{x}.csv") for x in ("bac", "arc")])
old["mag_key"] = old["User Genome"].map(key)
new = pd.concat([pd.read_csv(QC / f"gtdbtk_local_r232/gtdbtk.{x}.summary.tsv", sep="\t")
                 for x in ("bac120", "ar53")])
new["mag_key"] = new["user_genome"].map(key)
r232_names = [set() for _ in RANKS]
with open(R232_TAX) as fh:
    for line in fh:
        for i, n in enumerate(split_lineage(line.rstrip("\n").split("\t")[1])):
            r232_names[i].add(n)

g = old[["mag_key", "Classification"]].merge(
    new[["mag_key", "classification", "classification_method", "closest_genome_reference",
         "closest_genome_ani", "red_value", "warnings"]], on="mag_key", how="outer")
rows, renames, moves = [], collections.Counter(), []
for _, r in g.iterrows():
    a, b = split_lineage(r["Classification"]), split_lineage(r["classification"])
    kinds = {}
    for i, (x, y) in enumerate(zip(a, b)):
        if x == y:
            continue
        kind = ("resolved" if not x else "lost" if not y else
                "renamed" if x not in r232_names[i] else "moved")
        kinds[RANKS[i]] = kind
        if kind == "renamed":
            renames[(RANKS[i], x, y)] += 1
        elif kind == "moved":
            moves.append((r["mag_key"], RANKS[i], x, y))
    first = next((rk for rk in RANKS if rk in kinds), "")
    rows.append({"mag_key": r["mag_key"], "gtdb_2020": r["Classification"], "gtdb_r232": r["classification"],
                 "first_changed_rank": first,
                 "change_kinds": ";".join(f"{rk}:{k}" for rk, k in kinds.items()) or "unchanged",
                 "any_moved": any(k == "moved" for k in kinds.values()),
                 "r232_method": r["classification_method"], "r232_closest_ref": r["closest_genome_reference"],
                 "r232_closest_ani": r["closest_genome_ani"], "r232_red": r["red_value"],
                 "r232_warnings": r["warnings"]})
gt = pd.DataFrame(rows)

out = cm.merge(gt, on="mag_key", how="outer").sort_values("mag_key")
out.to_csv(QC / "comparison_vs_2020.csv", index=False)

# ---- report ------------------------------------------------------------------------------
derep = {key(l) for l in open(S / "jplfaria/Cliff MAGs/310_MAGS.txt") if l.strip()}
dom = lambda col: collections.Counter(split_lineage(x)[0] for x in g[col])
print(f"MAGs: 2020 GTDB {old.mag_key.nunique()}, r232 {new.mag_key.nunique()}, "
      f"all 310 derep in r232: {set(new.mag_key) == derep}")
print("domains 2020:", dict(dom("Classification")), "| r232:", dict(dom("classification")))
arc_old = {k for k, c in zip(g.mag_key, g.Classification) if split_lineage(c)[0] == "Archaea"}
arc_new = {k for k, c in zip(g.mag_key, g.classification) if split_lineage(c)[0] == "Archaea"}
print("same 30 archaea:", arc_old == arc_new)
print("first changed rank:", dict(collections.Counter(gt.first_changed_rank.replace("", "unchanged"))))
kc = collections.Counter(k.split(":")[1] for ks in gt.change_kinds if ks != "unchanged" for k in ks.split(";"))
print("rank-level differences by kind:", dict(kc))
print("MAGs with any real move (old name still exists in r232):", int(gt.any_moved.sum()))
for i in (5, 6):
    named = lambda col: sum(bool(split_lineage(x)[i]) for x in g[col])
    print(f"{RANKS[i]} named: 2020 {named('Classification')} -> r232 {named('classification')}")
print("top renames:")
for (rk, x, y), n in renames.most_common(12):
    print(f"   {rk:8s} {x} -> {y} ({n})")
print(f"moves ({len(moves)}):")
for k, rk, x, y in moves[:25]:
    print(f"   {k:28s} {rk:8s} {x} -> {y}")
print(f"\nwrote {QC / 'comparison_vs_2020.csv'} {out.shape}")
