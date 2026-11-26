#!/usr/bin/env python
"""IterativeIDs for the 310 dereplicated MAGs from the GTDB r232 taxonomy (2026-10 re-QC).

Same convention as the "IterativeIDs" cell of analysis_2_marimo.py: walk Kingdom -> Genus, stop
at the first empty rank (Species is never used), and name the MAG "<deepest named taxon>.<n>",
with n counting up in MAG order. Differences from the 2020 run:
  * taxonomy = GTDB-Tk 2.8.0 / r232 (MAG_QC/gtdbtk_local_r232/), not the 2020 GTDB-Tk 1.x table;
  * all 310 MAGs are named, including the 30 archaea the 2020 phylogeny file left out;
  * order = the 2020 order (data/Saltern_phylogeny.json, i.e. the 2020 GTDB table rows) for the
    280 bacteria, then the archaea in GTDB-Tk's ar53 summary order -- so a taxon whose member
    MAGs did not change keeps its old numbering.
Writes (the 2020 files are left untouched):
  ../data/Saltern_phylogeny_r232.json        MAG -> {Kingdom..Species}, the convention's input
  iterativeIDs_r232.json                     MAG (".contigs__" stripped) -> iterativeID
  iterativeID_levels_r232.json               iterativeID -> rank used
  iterativeID_phylums_r232.json              iterativeID -> phylum
  iterativeIDs_r232_vs_2020.csv              old vs new ID per MAG
Run from this folder: ~/Documents/py_venv/bin/python make_iterativeIDs_r232.py
"""
from json import dump, load
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
S = HERE.parent
GTDB = S / "MAG_QC/gtdbtk_local_r232"
LEVELS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

# ---- r232 taxonomy in the convention's {Kingdom..Species} shape ----------------------------
summary = pd.concat([pd.read_csv(GTDB / f"gtdbtk.{d}.summary.tsv", sep="\t") for d in ("bac120", "ar53")])
taxa = {}
for genome, lineage in zip(summary.user_genome, summary.classification):
    parts = dict(p.split("__", 1) for p in lineage.split(";"))
    taxa[genome] = {lvl: parts.get(lvl[0].lower() if lvl != "Kingdom" else "d", "") for lvl in LEVELS}

order_2020 = [m for m in load(open(S / "data/Saltern_phylogeny.json")) if m in taxa]
order = order_2020 + [m for m in summary.user_genome if m not in set(order_2020)]
assert len(order) == len(taxa) == 310, (len(order), len(taxa))
phylogeny = {m: taxa[m] for m in order}
dump(phylogeny, open(S / "data/Saltern_phylogeny_r232.json", "w"), indent=1)

# ---- the IterativeID convention (as in analysis_2_marimo.py) -----------------------------
ids, levels, phylums = {}, {}, {}
for MAG, _taxa in phylogeny.items():
    best_level = "Kingdom"
    for level, taxon in _taxa.items():
        if not isinstance(taxon, str) or level == "Species" or taxon == "":
            break
        best_level = level
    taxon = _taxa[best_level] + ".1"
    while taxon in ids:
        taxon, count = taxon.rsplit(".", 1)
        taxon = f"{taxon}.{int(count) + 1}"
    ids[taxon] = MAG.replace(".contigs__", "")
    levels[taxon] = best_level
    phylums[taxon] = _taxa["Phylum"]
mag_to_id = {v: k for k, v in ids.items()}
dump(mag_to_id, open(HERE / "iterativeIDs_r232.json", "w"))
dump(levels, open(HERE / "iterativeID_levels_r232.json", "w"))
dump(phylums, open(HERE / "iterativeID_phylums_r232.json", "w"))

# ---- old vs new ----------------------------------------------------------------------------
old = load(open(HERE / "iterativeIDs.json"))
rows = []
for mag, new_id in mag_to_id.items():
    old_id = old.get(mag, "")
    rows.append({"mag": mag, "iterativeID_2020": old_id, "iterativeID_r232": new_id,
                 "level_r232": levels[new_id], "phylum_r232": phylums[new_id],
                 "taxon_changed": bool(old_id) and old_id.rsplit(".", 1)[0] != new_id.rsplit(".", 1)[0],
                 "id_changed": old_id != new_id,
                 "gtdb_r232": ";".join(f"{lvl[0].lower() if lvl != 'Kingdom' else 'd'}__{t}"
                                       for lvl, t in phylogeny[mag + ".contigs__"].items())})
cmp_ = pd.DataFrame(rows)
cmp_.to_csv(HERE / "iterativeIDs_r232_vs_2020.csv", index=False)

print(f"{len(mag_to_id)} iterativeIDs ({cmp_.iterativeID_2020.eq('').sum()} MAGs new to the naming)")
print("levels used:", pd.Series(levels).value_counts().to_dict())
had = cmp_[cmp_.iterativeID_2020 != ""]
print(f"of the 280 named in 2020: identical ID {int((~had.id_changed).sum())}, "
      f"same taxon/new number {int((had.id_changed & ~had.taxon_changed).sum())}, "
      f"new taxon {int(had.taxon_changed.sum())}")
