"""Old manuscript Table 3 (2020 GTDB families) next to its corrected values.

The old table came from correlational/abundances.ipynb, which had three faults (fixed 2026-10-21):
  1. `'R1' in col or 'R2' in col` also matched the six R2A (reference) cores, so "unrestored" pooled R1 + R2 + R2A;
  2. cell 0 read only the bacterial GTDB-Tk summary, so the 30 archaea had no taxonomy, `fillna("unbinned_reads")`
     put them in with the unbinned reads, and they were dropped with that row ("68 % unbinned" = unbinned + archaea
     over those 18 cores);
  3. the notebook emitted unrestored / reference / restored, but the table header reads Unrestored / Restored / Reference.

The corrected values are the fixed notebook's (sum over a site type's cores, drop the unbinned row, normalise to the
binned total, group by family), recomputed from its output data/mag_abundances_taxonomy.csv. MAGs without a family,
bacterial or archaeal, are pooled as "undetermined". The R232 version of the table is in the R232 review docx.

Usage: ~/Documents/py_venv/bin/python table3_2020_corrected.py
"""
from pathlib import Path

import docx
import pandas as pd

S = Path(__file__).resolve().parents[1]
OUT = S / "Table 3 - corrected values (2020 GTDB families).csv"
SITES = {"Unrestored (R1, R2)": ("R1_", "R2_"), "Restored (SF2)": ("SF2_",), "Reference (R2A)": ("R2A_",)}


def main():
    mt = pd.read_csv(S / "data" / "mag_abundances_taxonomy.csv")
    unbinned = mt["MAG"].eq("unbinned_reads")
    assert mt.loc[~unbinned, "domain"].value_counts().to_dict() == {"Bacteria": 280, "Archaea": 30}
    fam = mt["family"].fillna("undetermined")
    named_archaeal = mt["domain"].eq("Archaea") & mt["family"].notna()
    fam[named_archaeal] = fam[named_archaeal] + " (archaeal)"

    cores = [c for c in mt.columns if c.startswith(("R1_", "R2_", "SF2_", "R2A_"))]
    new, unb = {}, {}
    for site, pre in SITES.items():
        cs = [c for c in cores if c.startswith(pre)]
        v = mt.loc[~unbinned, cs].sum(axis=1)
        new[site] = (v / v.sum()).groupby(fam[~unbinned]).sum()
        unb[site] = mt.loc[unbinned, cs].values.mean()
    new = pd.DataFrame(new).sort_values("Unrestored (R1, R2)", ascending=False)

    t = docx.Document(S / "Metabolic Manuscript.docx").tables[2]
    old = pd.DataFrame([[c.text.strip() for c in r.cells] for r in t.rows[1:]],
                       columns=[c.text.strip() for c in t.rows[0].cells]).set_index("GTDB-Tk Family").astype(float)
    old.columns = ["Old 'Unrestored total' (R1+R2+R2A)", "Old 'Restored total' (actually R2A)",
                   "Old 'Reference total' (actually SF2)"]

    rows = list(new.index[:10]) + [f for f in old.index if f not in new.index[:10]]
    tab = pd.concat([old.reindex(rows), new.loc[rows]], axis=1).round(6)
    tab.insert(0, "In corrected top 10", ["yes" if f in new.index[:10] else "no" for f in rows])
    tab.index.name = "Family (2020 GTDB)"
    tab.loc["Unbinned share of all reads (caption)"] = ["", None, None, None, *[round(unb[s], 4) for s in SITES]]
    tab.to_csv(OUT)
    print(tab.to_string())
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
