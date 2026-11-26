#!/usr/bin/env python3
"""Docs step 6 gate: what an openpyxl load/save round trip (the mechanism of a bin-keyed surgical cell edit) does to
each workbook the docs part will edit, read from the HEAD worktree and saved to gate_scratch (never into the repo).
Reports: cell values/formulas (all sheets), formula cells with/without cached values before and after (data_only
readers), zip parts that disappear (drawings, images, charts, comments, custom XML), merged cells, data validations,
conditional formatting, defined names, column widths."""
import io, json, sys, zipfile, warnings
from pathlib import Path
import openpyxl
warnings.filterwarnings("ignore")
RP = Path("/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan")
WT = RP / "head_worktree"
OUT = RP / "gate_scratch" / "D6_xlsx_roundtrip"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(RP / "gate"))
from harness import cmp_xlsx

BOOKS = ["MAG_target_gene_filled_UPDATED (1).xlsx", "MAG_target_gene_filled.xlsx", "cladeexcel.xlsx",
         "abundance and coverage (1).xlsx"]


def formula_stats(raw):
    wf = openpyxl.load_workbook(io.BytesIO(raw), data_only=False)
    wv = openpyxl.load_workbook(io.BytesIO(raw), data_only=True)
    out = {}
    for sn in wf.sheetnames:
        n = cached = 0
        for row in wf[sn].iter_rows():
            for c in row:
                if isinstance(c.value, str) and c.value.startswith("="):
                    n += 1
                    if wv[sn][c.coordinate].value is not None:
                        cached += 1
        if n:
            out[sn] = (n, cached)
    return out


def features(raw):
    wb = openpyxl.load_workbook(io.BytesIO(raw), data_only=False)
    f = {}
    for ws in wb.worksheets:
        f[ws.title] = dict(merged=len(ws.merged_cells.ranges), dv=len(ws.data_validations.dataValidation),
                           cf=len(ws.conditional_formatting), widths=len(ws.column_dimensions),
                           freeze=str(ws.freeze_panes), filt=str(ws.auto_filter.ref))
    f["_defined_names"] = len(list(wb.defined_names))
    return f


res = []
for b in BOOKS:
    src = WT / b
    raw = src.read_bytes()
    wb = openpyxl.load_workbook(io.BytesIO(raw))
    buf = io.BytesIO(); wb.save(buf); new = buf.getvalue()
    (OUT / (b + ".roundtrip.xlsx")).write_bytes(new)
    tier, how = cmp_xlsx(raw, new)
    pa = set(zipfile.ZipFile(io.BytesIO(raw)).namelist()); pb = set(zipfile.ZipFile(io.BytesIO(new)).namelist())
    lost = sorted(pa - pb); gained = sorted(pb - pa)
    fs_before, fs_after = formula_stats(raw), formula_stats(new)
    fa, fb = features(raw), features(new)
    feat_diff = {k: (fa.get(k), fb.get(k)) for k in set(fa) | set(fb) if fa.get(k) != fb.get(k)}
    r = dict(workbook=b, cells=f"{tier}: {how}", parts_lost=lost, parts_gained=gained,
             formula_cells_before={k: f"{v[0]} formulas, {v[1]} cached" for k, v in fs_before.items()},
             formula_cells_after={k: f"{v[0]} formulas, {v[1]} cached" for k, v in fs_after.items()},
             feature_diff=feat_diff)
    res.append(r)
    print(json.dumps(r, indent=1))
(OUT / "result.json").write_text(json.dumps(res, indent=1))
