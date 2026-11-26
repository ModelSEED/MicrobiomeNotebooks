"""Docs stage, step A: bin-keyed cell patch of the three docs-owned workbooks in ~/Documents/EmilyKin
(the builder workbook MAG_target_gene_filled.xlsx was fully patched by the data stage and is NOT touched).

  MAG_target_gene_filled_UPDATED (1).xlsx  (submission workbook)
      S1 / S2a / S3 / S4 / S5 / S7: MAG iteration ID -> R232 map (mag_iterativeID_old_to_new.json; also
      overwrites the stale 6a04588-era IDs, decision D4); GTDB rank columns -> R232 rank; lineage columns
      -> R232 classification (gtdbtk.bac120.summary.tsv).  Per-cell, s16 and 5-class columns untouched.
      S6_ANI: columns B..T copied by bin from the rebuilt builder sheet (MAG_target_gene_filled.xlsx,
      'Table S6_ANI', data stage B8a); columns U..AA (skani block) copied by bin from the rebuilt
      Table_S6_skani_novelty.tsv (data stage B8c).  At HEAD this copy rule reproduces the submission sheet
      exactly (gate mode --verify-head).
  cladeexcel.xlsx  S4: MAG iteration ID + genus.
  abundance and coverage (1).xlsx  Sheet1: MAG iteration ID.

Every write is keyed by the bin in the row's bin column; nothing is looked up by an ID string.
openpyxl (data_only=False) keeps formulas and styles; only planned cells get new values.
Writes a change list TSV (file, sheet, cell, bin, column, old, new, kind)."""
import csv, json, os, re, sys, subprocess, io
import openpyxl
csv.field_size_limit(sys.maxsize)
R = os.path.expanduser("~/Documents/EmilyKin")
P = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan"
OUT = f"{P}/apply_scratch/docs/workbook_changes.tsv"
RANKS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
VERIFY_HEAD = "--verify-head" in sys.argv
DRY = "--dry-run" in sys.argv or VERIFY_HEAD


def load_lineage(text):
    gt = {r["user_genome"]: r["classification"] for r in csv.DictReader(io.StringIO(text), delimiter="\t")}
    return gt, {b: dict(zip(RANKS, [p.split("__", 1)[1] for p in c.split(";")])) for b, c in gt.items()}


if VERIFY_HEAD:
    def git(path):
        return subprocess.check_output(["git", "-C", R, "show", f"HEAD:{path}"]).decode()
    idmap = json.loads(git("mag_iterativeID_old_to_new.json"))
    gt, lin = load_lineage(git("gtdbtk.bac120.summary.tsv").replace("\r\n", "\n"))
    S6TSV = list(csv.DictReader(io.StringIO(git("Table_S6_skani_novelty.tsv").replace("\r\n", "\n")), delimiter="\t"))
    BUILDER = f"{P}/apply_scratch/data/backups/MAG_target_gene_filled.HEAD.xlsx"
else:
    idmap = json.load(open(f"{R}/mag_iterativeID_old_to_new.json"))
    gt, lin = load_lineage(open(f"{R}/gtdbtk.bac120.summary.tsv", newline="").read().replace("\r\n", "\n"))
    S6TSV = list(csv.DictReader(open(f"{R}/Table_S6_skani_novelty.tsv", newline=""), delimiter="\t"))
    BUILDER = f"{R}/MAG_target_gene_filled.xlsx"
assert len(idmap) == len(gt) == 291 and len(S6TSV) == 291
s6tsv = {r["MAG"]: r for r in S6TSV}
SKANI_COLS = ["skani_best_ref", "skani_best_taxonomy", "skani_ANI", "skani_AF_query", "skani_AF_ref", "ANI_class",
              "Novelty_interpretation"]
NUMERIC = {"skani_ANI", "skani_AF_query", "skani_AF_ref"}


def typed(col, s):
    if s == "":
        return None
    if col in NUMERIC:
        return int(s) if re.fullmatch(r"-?\d+", s) else float(s)
    return s


# builder S6 rows by bin (columns B..T = index 1..19)
bws = openpyxl.load_workbook(BUILDER, data_only=False)["Table S6_ANI"]
bhdr = [bws.cell(1, c).value for c in range(1, 21)]
assert bhdr[0] == "MAG bin name" and bws.max_column == 20
brow = {}
for r in range(2, bws.max_row + 1):
    b = bws.cell(r, 1).value
    if b in idmap:
        assert b not in brow
        brow[b] = [bws.cell(r, c).value for c in range(1, 21)]
assert len(brow) == 291

# (file, sheet, header row, bin column header, {column header: kind})
PLAN = [
    ("MAG_target_gene_filled_UPDATED (1).xlsx", [
        ("Table S1_MAG quality", 3, "MAG", {"MAG iteration ID": "id", "Phylum": "Phylum", "Class": "Class",
                                            "Order": "Order", "Family": "Family", "Genus": "Genus",
                                            "Taxonomy ": "lineage"}),
        ("Table S2a_RPKM per MAG", 3, "MAG", {"MAG iteration ID": "id"}),
        ("Table S3_gene copies", 3, "MAG", {"MAG iteration ID": "id", "Taxonomy": "lineage"}),
        ("Table S4_nosZ_genecopy_type", 3, "MAG", {"MAG iteration ID": "id", "genus": "Genus"}),
        ("Table S5 _ppk_ACC", 3, "MAG", {"MAG iteration ID": "id"}),
        ("Table S6_ANI", 1, "MAG", "S6"),
        ("Table S7_functional", 1, "MAG", {"MAG iteration ID": "id", "Genus": "Genus", "Family": "Family",
                                          "Phylum": "Phylum"}),
    ]),
    ("cladeexcel.xlsx", [
        ("Table S4_nosZ_genecopy_type", 3, "MAG", {"MAG iteration ID": "id", "genus": "Genus"}),
    ]),
    ("abundance and coverage (1).xlsx", [
        ("Sheet1", 1, None, {"MAG iteration ID": "id"}),   # bin column = column A (no header)
    ]),
]

changes = []
for fname, sheets in PLAN:
    wb = openpyxl.load_workbook(f"{R}/{fname}", data_only=False)
    for sh, hrow, bcol, cols in sheets:
        ws = wb[sh]
        hdr = [ws.cell(hrow, c).value for c in range(1, ws.max_column + 1)]
        bi = 1 if bcol is None else hdr.index(bcol) + 1
        assert bcol is None or hdr.count(bcol) == 1
        if cols == "S6":
            assert hdr[1:20] == bhdr[1:20] and hdr[20:27] == SKANI_COLS, "S6 header layout"
            ci = None
        else:
            ci = {}
            for c in cols:
                assert hdr.count(c) == 1, (fname, sh, c)
                ci[c] = hdr.index(c) + 1
        seen = set()
        for r in range(hrow + 1, ws.max_row + 1):
            b = ws.cell(r, bi).value
            if b is None or b not in idmap:
                continue
            assert b not in seen, (fname, sh, b)
            seen.add(b)
            if cols == "S6":
                targets = [(c + 1, hdr[c], brow[b][c], "s6_builder") for c in range(1, 20)]
                targets += [(21 + j, col, typed(col, s6tsv[b][col]), "s6_skani") for j, col in enumerate(SKANI_COLS)]
            else:
                targets = []
                for c, kind in cols.items():
                    if kind == "id":
                        new = idmap[b]
                    elif kind == "lineage":
                        new = gt[b]
                    else:
                        new = lin[b][kind] or None
                    targets.append((ci[c], c, new, kind))
            for cidx, cname, new, kind in targets:
                cell = ws.cell(r, cidx)
                old = cell.value
                same = old == new or (old in (None, "") and new in (None, ""))
                if not same and isinstance(old, (int, float)) and isinstance(new, (int, float)):
                    same = float(old) == float(new) and type(old) is type(new)
                if not same:
                    changes.append([fname, sh, cell.coordinate, b, cname, "" if old is None else old,
                                    "" if new is None else new, kind])
                    if not DRY:
                        cell.value = new
        print(f"{fname} | {sh}: {len(seen)} bin rows")
    if not DRY:
        wb.save(f"{R}/{fname}")

from collections import Counter
for (f, sh, c, k), n in sorted(Counter((x[0], x[1], x[4], x[7]) for x in changes).items()):
    print(f"  {f} | {sh} | {c} [{k}]: {n}")
print("total cells changed:", len(changes), "(dry run)" if DRY else "")
if VERIFY_HEAD:
    s6 = [x for x in changes if x[1] == "Table S6_ANI"]
    print("VERIFY-HEAD: S6 cells that the copy rule would change at HEAD:", len(s6))
    print("  of which ID-only (the D4 drift the builder rebuild also fixes):",
          sum(1 for x in s6 if x[4] == "MAG iteration ID"))
elif not DRY:
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["file", "sheet", "cell", "bin", "column", "old", "new", "kind"])
        w.writerows(changes)
    print("->", OUT)
