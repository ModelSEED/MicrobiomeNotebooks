#!/usr/bin/env python
"""Compare the kofam_rollout validation run (8 first-pass control bins, this track's re-pointed call_genes.py /
make_chunks.py / postprocess.py, 68 marker profiles, -Z = the 935-MAG total) with the kofam track's marker test on
the same bins (kofam/test/markers: chunk 01 of the 1,324 run, -Z 3,911,146).

Expected: every call identical; only E-values differ, by the factor Z_new / Z_old (up to the 2-digit rounding of
the tblout).

Inputs (copied from poplar):  validation/new/{chunks,per_mag}   <- kofam_rollout/validation/results/
                              validation/orig/{chunks,per_mag}  <- kofam/test/markers/results/ (8 control MAGs only)
usage: ~/Documents/py_venv/bin/python -I scripts/compare_validation.py   (from annotation_runs/kofam_rollout/)
writes validation/validation_comparison.tsv (cell-level summary) and prints the verdict.
"""
from pathlib import Path

import numpy as np
import pandas as pd

K = Path(__file__).resolve().parents[1]
V = K / "validation"
ctrl = [l.strip() for l in open(K / "meta/controls_validation.txt") if l.strip()]
z = {}
for line in open(V / "Z.txt"):
    k, v = line.split()
    z[k] = int(v)
ratio = z["Z_new"] / z["Z_old"]
rows = []


def report(table, n_rows, cols, cells, differing, note=""):
    rows.append(dict(table=table, rows_compared=n_rows, columns_compared=cols, cells_compared=cells,
                     cells_differing=differing, note=note))
    print(f"{table:28s} rows {n_rows:6d}  cells {cells:8d}  differing {differing}  {note}")


# 1. per-MAG marker counts + guild columns
a = pd.read_csv(V / "orig/chunks/chunk_01.mag_markers_guilds.tsv", sep="\t", dtype=str).set_index("mag")
b = pd.read_csv(V / "new/chunks/chunk_01.mag_markers_guilds.tsv", sep="\t", dtype=str).set_index("mag")
assert set(b.index) == set(ctrl) and set(ctrl) <= set(a.index)
a = a.loc[ctrl]
b = b.loc[ctrl]
assert list(a.columns) == list(b.columns), set(a.columns) ^ set(b.columns)
cols = [c for c in a.columns if c != "chunk"]
d = (a[cols].fillna("") != b[cols].fillna(""))
report("mag_markers_guilds (all)", len(ctrl), len(cols), d.size, int(d.values.sum()), "every column except 'chunk'")
mk = [c for c in cols if c.startswith("n_") and c not in ("n_genes", "n_genes_with_ko",
                                                          "n_marker_genes_on_short_maxbin2_contigs")]
report("  of which n_<marker>", len(ctrl), len(mk), len(ctrl) * len(mk), int(d[mk].values.sum()),
       "68 marker best-KO counts (the 544 cells of kofam/test/control_comparison.tsv)")

# 1b. against kofam/test/control_comparison.tsv, column poplar_2026_10
cc = pd.read_csv(K.parent / "kofam/test/control_comparison.tsv", sep="\t")
nb = b[mk].astype(int).rename(columns=lambda c: c[2:]).stack().rename("new").reset_index()
nb.columns = ["mag", "marker", "new"]
m = cc.merge(nb, on=["mag", "marker"], how="outer", indicator=True)
assert (m._merge == "both").all(), m._merge.value_counts()
report("control_comparison.tsv", len(ctrl), 68, len(m), int((m.poplar_2026_10 != m.new).sum()),
       "kofam test (KOfam 2026-10) vs this run")


def gene_level(name, f, key, evalue_col):
    a = pd.read_csv(V / "orig/chunks" / f, sep="\t", dtype=str)
    b = pd.read_csv(V / "new/chunks" / f, sep="\t", dtype=str)
    a = a[a.mag.isin(ctrl)].set_index(key).sort_index()
    b = b.set_index(key).sort_index()
    assert b.mag.isin(ctrl).all() and a.index.is_unique and b.index.is_unique
    same_keys = a.index.equals(b.index)
    cols = [c for c in a.columns if c != evalue_col]
    if not same_keys:
        report(name, len(b), len(cols), 0, -1, f"KEY SETS DIFFER: {len(a.index.difference(b.index))} only orig, "
               f"{len(b.index.difference(a.index))} only new")
        return
    dd = (a[cols].fillna("") != b[cols].fillna(""))
    report(name, len(b), len(cols), dd.size, int(dd.values.sum()), f"all columns except {evalue_col}")
    ea, eb = a[evalue_col].astype(float), b[evalue_col].astype(float)
    nz = (ea > 0) & (eb > 0)
    r = (eb[nz] / ea[nz]) / ratio
    report(f"  {evalue_col}", len(b), 1, len(b), int((ea != eb).sum()),
           f"differ as expected; new/old / (Z_new/Z_old): median {np.median(r):.4f}, range {r.min():.3f}-{r.max():.3f}; "
           f"zero in both {int(((ea == 0) & (eb == 0)).sum())}, zero in one only {int(((ea == 0) ^ (eb == 0)).sum())}")


gene_level("marker_hits (gene x KO)", "chunk_01.marker_hits.tsv.gz", ["gene", "ko"], "full_evalue")
gene_level("genes_ko (gene)", "chunk_01.genes_ko.tsv.gz", ["gene"], "best_full_evalue")

# per-MAG gene tables (every gene, best KO, all passing KOs)
n = cells = diff = ediff = 0
for m in ctrl:
    a = pd.read_csv(V / f"orig/per_mag/{m}.kofam.tsv.gz", sep="\t", dtype=str)
    b = pd.read_csv(V / f"new/per_mag/{m}.kofam.tsv.gz", sep="\t", dtype=str)
    assert list(a.gene) == list(b.gene)
    cols = [c for c in a.columns if c != "kofam_full_evalue"]
    dd = (a[cols].fillna("") != b[cols].fillna(""))
    n += len(a); cells += dd.size; diff += int(dd.values.sum())
    ediff += int((a.kofam_full_evalue.fillna("") != b.kofam_full_evalue.fillna("")).sum())
report("per_mag tables (gene)", n, len(cols), cells, diff, f"all columns except kofam_full_evalue "
       f"({ediff} E-values differ)")

t = pd.DataFrame(rows)
t.to_csv(V / "validation_comparison.tsv", sep="\t", index=False)
calls = t[~t.table.str.contains("evalue")]
bad = calls[(calls.cells_differing != 0)]
print("\nVERDICT:", "PASS - all calls identical; only E-values differ" if bad.empty else f"FAIL\n{bad.to_string()}")
