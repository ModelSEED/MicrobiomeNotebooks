#!/usr/bin/env python
"""Compare the kofam_coasm38 validation control chunk (2 rollout MAGs, meta/controls_validation.txt; this track's
call_genes.py / make_chunks.py / postprocess.py; all 27,981 profiles; -Z = the 38-MAG total) with kofam_rollout's
production output for the same 2 MAGs (-Z 2,732,576).

Expected: identical gene calls (byte-identical files, checked on poplar by run_validation.sh genes) and identical
calls in every table; only E-values differ, by the factor Z_new / Z_old (up to the 2-significant-digit rounding of
the tblout). A best-KO / ordering difference would be acceptable only as an E-value rounding tie, and is listed.

Inputs:  validation/new/{chunks,per_mag,all_*}  <- poplar kofam_coasm38/validation/results/
         validation/logs/{genes_md5_compare.txt,Z.txt}  <- poplar kofam_coasm38/validation/logs/
         ../kofam_rollout/results/{per_mag,all_*}  (kofam_rollout production, copied back with md5)
usage: ~/Documents/py_venv/bin/python -I scripts/compare_validation.py   (from annotation_runs/kofam_coasm38/)
writes validation/validation_comparison.tsv (cell-level summary), validation/differing_cells.tsv, prints the verdict.
kofam_coasm38 rewrite of kofam_rollout/scripts/compare_validation.py (whose reference was a marker-only test run).
"""
from pathlib import Path

import numpy as np
import pandas as pd

K = Path(__file__).resolve().parents[1]
V = K / "validation"
RO = K.parent / "kofam_rollout/results"
ctrl = [l.strip() for l in open(K / "meta/controls_validation.txt") if l.strip()]
z = dict(line.split() for line in open(V / "logs/Z.txt"))
ratio = int(z["Z_new"]) / int(z["Z_old"])
rows, diffs = [], []


def report(table, n_rows, cols, cells, differing, note=""):
    rows.append(dict(table=table, rows_compared=n_rows, columns_compared=cols, cells_compared=cells,
                     cells_differing=differing, note=note))
    print(f"{table:34s} rows {n_rows:6d}  cells {cells:8d}  differing {differing}  {note}")


def cmp_frames(name, a, b, cols):
    """a = kofam_rollout, b = this run; same index required."""
    assert a.index.is_unique and b.index.is_unique
    if not a.index.equals(b.index):
        report(name, len(b), len(cols), 0, -1, f"KEY SETS DIFFER: {len(a.index.difference(b.index))} only rollout, "
               f"{len(b.index.difference(a.index))} only new")
        return
    dd = a[cols].fillna("") != b[cols].fillna("")
    report(name, len(b), len(cols), dd.size, int(dd.values.sum()), f"columns: all except {sorted(set(a.columns) - set(cols))}")
    for key, col in zip(*np.nonzero(dd.values)):
        diffs.append(dict(table=name, key=str(a.index[key]), column=cols[col], kofam_rollout=a[cols[col]].iloc[key],
                          kofam_coasm38=b[cols[col]].iloc[key]))


def evalue_ratio(name, ea, eb):
    ea, eb = ea.astype(float), eb.astype(float)
    nz = (ea > 0) & (eb > 0)
    r = (eb[nz] / ea[nz]) / ratio
    report(f"  {name}", len(ea), 1, len(ea), int((ea != eb).sum()),
           f"E-values differ as expected; new/old / (Z_new/Z_old = {ratio:.5f}): median {np.median(r):.4f}, "
           f"range {r.min():.3f}-{r.max():.3f}; zero in both {int(((ea == 0) & (eb == 0)).sum())}, "
           f"zero in one only {int(((ea == 0) ^ (eb == 0)).sum())}")


# 0. gene calls: byte-identical files (md5 on poplar)
g = pd.read_csv(V / "logs/genes_md5_compare.txt", sep="\t", header=None, names=["mag", "file", "new", "rollout", "verdict"])
assert set(g.mag) == set(ctrl) and len(g) == 4 * len(ctrl)
report("gene calls (.faa/.gff/.genes.tsv/.done)", len(g), 1, len(g), int((g.new != g.rollout).sum()),
       "md5 of this run's files vs kofam_rollout/genes/")

# 1. per-MAG gene tables (every gene; best KO = kofam_id; all passing KOs)
pa = pd.concat([pd.read_csv(RO / f"per_mag/{m}.kofam.tsv.gz", sep="\t", dtype=str) for m in ctrl]).set_index("gene")
pb = pd.concat([pd.read_csv(V / f"new/per_mag/{m}.kofam.tsv.gz", sep="\t", dtype=str) for m in ctrl]).set_index("gene")
cols = [c for c in pa.columns if c != "kofam_full_evalue"]
cmp_frames("per_mag tables (gene)", pa, pb, cols)
dd = pa.kofam_id.fillna("") != pb.kofam_id.reindex(pa.index).fillna("")
report("  of which kofam_id (best KO)", len(pa), 1, len(pa), int(dd.sum()),
       f"{int(pa.kofam_id.notna().sum())} genes with a best KO")
evalue_ratio("kofam_full_evalue", pa.kofam_full_evalue.dropna(), pb.kofam_full_evalue.reindex(pa.kofam_full_evalue.dropna().index))


def load(name, key, ro_file, new_file):
    a = pd.read_csv(RO / ro_file, sep="\t", dtype=str)
    a = a[a.mag.isin(ctrl)].set_index(key).sort_index()
    b = pd.read_csv(V / "new" / new_file, sep="\t", dtype=str)
    assert b.mag.isin(ctrl).all()
    return a, b.set_index(key).sort_index()


# 2. genes with >= 1 passing KO
a, b = load("genes_ko", ["gene"], "all_genes_ko.tsv.gz", "all_genes_ko.tsv.gz")
cmp_frames("genes_ko (gene)", a, b, [c for c in a.columns if c != "best_full_evalue"])
evalue_ratio("best_full_evalue", a.best_full_evalue, b.best_full_evalue.reindex(a.index))
# 3. marker hits (gene x KO, any reported score)
a, b = load("marker_hits", ["gene", "ko"], "all_marker_hits.tsv.gz", "all_marker_hits.tsv.gz")
cmp_frames("marker_hits (gene x KO)", a, b, [c for c in a.columns if c != "full_evalue"])
evalue_ratio("full_evalue", a.full_evalue, b.full_evalue.reindex(a.index))
# 4. MAG x KO counts (best KO, any passing KO)
a, b = load("mag_ko_counts", ["mag", "ko"], "all_mag_ko_counts.tsv", "all_mag_ko_counts.tsv")
cmp_frames("mag_ko_counts (MAG x KO)", a, b, list(a.columns))
# 5. per-MAG marker counts + guild calls
a = pd.read_csv(RO / "all_mag_markers_guilds.tsv", sep="\t", dtype=str).set_index("mag").loc[ctrl]
b = pd.read_csv(V / "new/all_mag_markers_guilds.tsv", sep="\t", dtype=str).set_index("mag").loc[ctrl]
assert list(a.columns) == list(b.columns)
cmp_frames("mag_markers_guilds (MAG)", a, b, [c for c in a.columns if c != "chunk"])
# 6. no-threshold hits (informational: postprocess.py keeps rows with full E <= 1e-5, which depends on -Z)
a, b = load("no_threshold_hits", ["gene", "ko"], "all_no_threshold_hits.tsv.gz", "all_no_threshold_hits.tsv.gz")
common = a.index.intersection(b.index)
c6 = [c for c in a.columns if c != "full_evalue"]
d6 = (a.loc[common, c6].fillna("") != b.loc[common, c6].fillna(""))
rows.append(dict(table="no_threshold_hits (informational)", rows_compared=len(common), columns_compared=len(c6),
                 cells_compared=d6.size, cells_differing=int(d6.values.sum()),
                 note=f"Z-dependent row set (E <= 1e-5): rollout {len(a)} rows, this run {len(b)}; rollout rows missing "
                      f"here {len(a.index.difference(b.index))}; extra here {len(b.index.difference(a.index))}"))
print(f"{'no_threshold_hits (informational)':34s} {rows[-1]['note']}; common rows differing cells {rows[-1]['cells_differing']}")

t = pd.DataFrame(rows)
t.to_csv(V / "validation_comparison.tsv", sep="\t", index=False)
pd.DataFrame(diffs, columns=["table", "key", "column", "kofam_rollout", "kofam_coasm38"]).to_csv(
    V / "differing_cells.tsv", sep="\t", index=False)
calls = t[~t.table.str.contains("evalue|informational")]
n_cells, n_diff = int(calls[calls.table.str.match(r"^\S")].cells_compared.sum()), int(calls[calls.table.str.match(r"^\S")].cells_differing.sum())
print(f"\ncall cells compared (top-level tables) {n_cells}, differing {n_diff}")
bad = calls[calls.cells_differing != 0]
print("VERDICT:", "PASS - all calls identical; only E-values differ" if bad.empty else
      f"DIFFERENCES (check for E-rounding ties, validation/differing_cells.tsv)\n{bad.to_string()}")
