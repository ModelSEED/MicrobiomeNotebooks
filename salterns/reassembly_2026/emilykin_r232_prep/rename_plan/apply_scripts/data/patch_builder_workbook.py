"""B6 (builder workbook only): bin-keyed cell patch of ~/Documents/EmilyKin/MAG_target_gene_filled.xlsx,
sheets S1, S1b, S2a, S3, S4, S5: MAG iterativeID -> R232 map (mag_iterativeID_old_to_new.json; also
overwrites the stale 6a04588-era IDs, D4), Phylum..Genus / genus -> R232 rank, Taxonomy lineage -> R232
classification (gtdbtk.bac120.summary.tsv).  A blank R232 rank is written as an empty cell (None), as the
sheet stores blanks.  S6 is rebuilt by scripts/build_table_s6_ani.py and S7 by
scripts/build_table_s7_classification.py.  build_xlsx_291.py is NOT used (gate: tier C, S1 layout drift).
Writes a change list TSV (sheet, cell, bin, column, old, new, kind) to apply_scratch/data/."""
import csv, json, os, sys
import openpyxl
csv.field_size_limit(sys.maxsize)
R = os.path.expanduser("~/Documents/EmilyKin")
OUT = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan/apply_scratch/data/builder_workbook_changes.tsv"
RANKS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
idmap = json.load(open(f"{R}/mag_iterativeID_old_to_new.json"))
gt = {r["user_genome"]: r["classification"] for r in csv.DictReader(open(f"{R}/gtdbtk.bac120.summary.tsv", newline=""), delimiter="\t")}
lin = {b: dict(zip(RANKS, [p.split("__", 1)[1] for p in c.split(";")])) for b, c in gt.items()}
assert len(idmap) == len(gt) == 291
SHEETS = [  # sheet, header row, bin column, {column: kind}
    ("Table S1_MAG quality", 3, "MAG bin name", {"MAG iteration ID": "id", "Phylum": "Phylum", "Class": "Class",
        "Order": "Order", "Family": "Family", "Genus": "Genus", "Taxonomy ": "lineage"}),
    ("Table S1b_relabun_totalreads", 3, "MAG bin name", {"MAG iteration ID": "id", "Phylum": "Phylum", "Class": "Class",
        "Order": "Order", "Family": "Family", "Genus": "Genus", "Taxonomy ": "lineage"}),
    ("Table S2a_RPKM per MAG", 3, "MAGs", {"iterative ID": "id"}),
    ("Table S3_gene copies", 3, "MAGs", {"iterative ID": "id", "Taxonomy": "lineage"}),
    ("Table S4_nosZ_genecopy_type", 3, "MAG", {"MAG iteration ID": "id", "genus": "Genus"}),
    ("Table S5 _ppk_ACC", 3, "MAG", {"MAG iteration ID": "id"}),
]
wb = openpyxl.load_workbook(f"{R}/MAG_target_gene_filled.xlsx")
changes = []
summary = []
for sh, hrow, bcol, cols in SHEETS:
    ws = wb[sh]
    hdr = [ws.cell(hrow, c).value for c in range(1, ws.max_column + 1)]
    assert hdr.count(bcol) == 1, (sh, bcol)
    bi = hdr.index(bcol) + 1
    ci = {}
    for c in cols:
        assert hdr.count(c) == 1, (sh, c)
        ci[c] = hdr.index(c) + 1
    nrows = 0
    seen = set()
    for r in range(hrow + 1, ws.max_row + 1):
        b = ws.cell(r, bi).value
        if b is None or b not in idmap:
            continue
        assert b not in seen, (sh, b); seen.add(b)
        nrows += 1
        for c, kind in cols.items():
            cell = ws.cell(r, ci[c])
            if kind == "id":
                new = idmap[b]
            elif kind == "lineage":
                new = gt[b]
            else:
                new = lin[b][kind] or None
            if cell.value != new:
                changes.append([sh, cell.coordinate, b, c, "" if cell.value is None else cell.value,
                                "" if new is None else new, kind])
                cell.value = new
    summary.append((sh, nrows))
wb.save(f"{R}/MAG_target_gene_filled.xlsx")
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["sheet", "cell", "bin", "column", "old", "new", "kind"])
    w.writerows(changes)
from collections import Counter
print("rows per sheet:", summary)
for (sh, c), n in sorted(Counter((x[0], x[3]) for x in changes).items()):
    print(f"  {sh} | {c}: {n}")
print("total cells changed:", len(changes), "->", OUT)
