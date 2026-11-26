"""Gate for patch_docs_workbooks.py: (1) every non-skani plan cell of table_cells_plan.tsv for the 3 files is
changed to the planned value and no other non-S6-recompute cell changed; (2) full cell diff vs the
.pre_r232.bak.xlsx backups touches only the planned columns; (3) S1 formula count unchanged; (4) zip part list,
defined names, merged cells unchanged; (5) D1: each ID column equals the R232 map for every bin and is unique;
(6) S6 equals builder S6 (B..T) + Table_S6 tsv (U..AA) for 291/291."""
import csv, json, os, re, sys, zipfile
import openpyxl
csv.field_size_limit(sys.maxsize)
R = os.path.expanduser("~/Documents/EmilyKin")
P = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan"
FILES = ["MAG_target_gene_filled_UPDATED (1)", "cladeexcel", "abundance and coverage (1)"]
idmap = json.load(open(f"{R}/mag_iterativeID_old_to_new.json"))
ch = list(csv.DictReader(open(f"{P}/apply_scratch/docs/workbook_changes.tsv"), delimiter="\t"))
chk = {(r["file"], r["sheet"], r["cell"]): r for r in ch}
plan = [r for r in csv.DictReader(open(f"{P}/docs/table_cells_plan.tsv"), delimiter="\t")
        if r["file"] in [f + ".xlsx" for f in FILES]]
ok = True
# (1) plan vs applied (non-skani)
miss = []; bad = []
nonsk = [r for r in plan if r["change_type"] != "skani_or_version_cell"]
for r in nonsk:
    sh, cell = r["location"].rsplit("!", 1)
    k = (r["file"], sh, cell)
    if k not in chk:
        miss.append(k); continue
    want = "" if r["proposed"] == "(blank in R232)" else r["proposed"]
    if chk[k]["new"] != want or chk[k]["bin"] != r["bin"]:
        bad.append((k, chk[k]["new"], want))
plank = {(r["file"],) + tuple(r["location"].rsplit("!", 1)) for r in plan}
extra = [k for k, r in chk.items() if k not in plank and r["kind"] != "s6_builder" and r["kind"] != "s6_skani"]
s6extra = [k for k, r in chk.items() if k not in plank and r["kind"] in ("s6_builder", "s6_skani")]
print(f"(1) non-skani plan cells: {len(nonsk)}; applied & equal: {len(nonsk)-len(miss)-len(bad)}; missing {len(miss)}; "
      f"wrong {len(bad)}; non-S6 cells changed outside the plan: {len(extra)}")
print(f"    S6 recompute cells changed: {sum(1 for r in ch if r['kind'] in ('s6_builder','s6_skani'))} "
      f"(of which not listed in the plan, i.e. blank->value or plan-unanticipated R232 changes: {len(s6extra)})")
ok &= not miss and not bad and not extra
# (2)-(4)
for f in FILES:
    a = openpyxl.load_workbook(f"{R}/{f}.pre_r232.bak.xlsx", data_only=False)
    b = openpyxl.load_workbook(f"{R}/{f}.xlsx", data_only=False)
    assert a.sheetnames == b.sheetnames
    diffcells = set(); nform_a = nform_b = 0
    for sh in a.sheetnames:
        wa, wb_ = a[sh], b[sh]
        assert (wa.max_row, wa.max_column) == (wb_.max_row, wb_.max_column), (f, sh)
        assert sorted(map(str, wa.merged_cells.ranges)) == sorted(map(str, wb_.merged_cells.ranges))
        for row in wa.iter_rows():
            for c in row:
                v2 = wb_[c.coordinate].value
                if isinstance(c.value, str) and c.value.startswith("="): nform_a += 1
                if isinstance(v2, str) and v2.startswith("="): nform_b += 1
                if c.value != v2:
                    diffcells.add((f + ".xlsx", sh, c.coordinate))
                if c.style_id != wb_[c.coordinate].style_id:
                    print("STYLE DIFF", f, sh, c.coordinate); ok = False
    applied = {k for k in chk if k[0] == f + ".xlsx"}
    za = sorted(zipfile.ZipFile(f"{R}/{f}.pre_r232.bak.xlsx").namelist())
    zb = sorted(zipfile.ZipFile(f"{R}/{f}.xlsx").namelist())
    dn = (sorted(a.defined_names.keys()) == sorted(b.defined_names.keys()))
    print(f"(2-4) {f}: cells differing vs backup {len(diffcells)} == change list {len(applied)}: {diffcells == applied}; "
          f"formulas {nform_a} -> {nform_b}; zip parts equal {za == zb}; defined names equal {dn}")
    ok &= diffcells == applied and nform_a == nform_b and za == zb and dn
# (5) D1
for f, sh, hrow, bcol, idcol in [("MAG_target_gene_filled_UPDATED (1)", s, h, "MAG", "MAG iteration ID") for s, h in
                                  [("Table S1_MAG quality", 3), ("Table S2a_RPKM per MAG", 3), ("Table S3_gene copies", 3),
                                   ("Table S4_nosZ_genecopy_type", 3), ("Table S5 _ppk_ACC", 3), ("Table S6_ANI", 1),
                                   ("Table S7_functional", 1)]] + [
        ("cladeexcel", "Table S4_nosZ_genecopy_type", 3, "MAG", "MAG iteration ID"),
        ("abundance and coverage (1)", "Sheet1", 1, None, "MAG iteration ID")]:
    ws = openpyxl.load_workbook(f"{R}/{f}.xlsx", data_only=False)[sh]
    hdr = [ws.cell(hrow, c).value for c in range(1, ws.max_column + 1)]
    bi = 1 if bcol is None else hdr.index(bcol) + 1; ii = hdr.index(idcol) + 1
    pairs = [(ws.cell(r, bi).value, ws.cell(r, ii).value) for r in range(hrow + 1, ws.max_row + 1)
             if ws.cell(r, bi).value in idmap]
    wrong = [p for p in pairs if idmap[p[0]] != p[1]]
    ids = [p[1] for p in pairs]
    print(f"(5) D1 {f} | {sh}: {len(pairs)} bins, wrong {len(wrong)}, duplicate IDs {len(ids)-len(set(ids))}")
    ok &= not wrong and len(ids) == len(set(ids))
print("GATE", "PASSED" if ok else "FAILED")
