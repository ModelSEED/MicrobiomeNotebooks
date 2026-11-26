#!/usr/bin/env python
"""Compare the poplar marker-profile test on chunk 01 (KOfam 2026-10-26, Prodigal bioconda build, contigs >= 1,000 bp)
with the local control test (annotation_choice/kofam_control_test/marker_calls_by_mag.tsv, KOfam 2026-05-28,
gtdbtk_venv Prodigal, min_contig 1000) for the 8 re-assembly control bins. Both runs searched only the marker
profiles, so 'best KO' is best among markers in both.

usage: ~/Documents/py_venv/bin/python compare_controls.py <test results dir with chunks/chunk_01.mag_markers_guilds.tsv>
writes test/control_comparison.tsv next to the results dir and prints the differing cells.
"""
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
CT = RUN.parents[1] / "annotation_choice" / "kofam_control_test"
res = Path(sys.argv[1])
LABEL = {
    "new_Roseovarius_maxbin2_R1_A_D2": "maxbin2.Bin.048_sub.fasta__R1_A_D2",
    "new_Roseovarius_maxbin2_R1_B_D2": "maxbin2.Bin.032_sub.fasta__R1_B_D2",
    "new_Methanolobus_concoct_R1_sed": "concoct.bin.101.fasta__R1_sediment",
}
for b in ["concoct.bin.061.fasta__R1_C_D1", "concoct.bin.157_sub.fasta__R1_sediment", "concoct.bin.018.fasta__R2_B_D2",
          "metabat2.bin.062.fasta__R2_restored_C_black", "concoct.bin.058_sub.fasta__R2_sediment"]:
    LABEL[f"new_Albimonas_{b.split('.fasta__')[0].replace('.', '')}_{b.split('__')[1]}"] = b

old = pd.read_csv(CT / "marker_calls_by_mag.tsv", sep="\t")
old = old[(old["min_contig"] == 1000) & old["mag"].isin(LABEL)].copy()
old["mag"] = old["mag"].map(LABEL)
old = old.drop(columns="min_contig").set_index("mag")
new = pd.read_csv(res / "chunks" / "chunk_01.mag_markers_guilds.tsv", sep="\t").set_index("mag")
new = new.loc[list(LABEL.values()), [c for c in new.columns if c.startswith("n_") and c not in
                                     ("n_genes", "n_genes_with_ko", "n_marker_genes_on_short_maxbin2_contigs")]]
new.columns = [c[2:] for c in new.columns]
genes = sorted(set(old.columns) | set(new.columns))
rows = []
for m in LABEL.values():
    for g in genes:
        o = int(old.at[m, g]) if g in old.columns else 0
        n = int(new.at[m, g]) if g in new.columns else 0
        rows.append({"mag": m, "marker": g, "local_2026_05": o, "poplar_2026_10": n, "diff": n - o})
t = pd.DataFrame(rows)
out = RUN / "test" / "control_comparison.tsv"
out.parent.mkdir(exist_ok=True)
t.to_csv(out, sep="\t", index=False)
d = t[t["diff"] != 0]
print(f"{len(t)} MAG x marker cells, {len(d)} differ")
print(d.to_string(index=False))
