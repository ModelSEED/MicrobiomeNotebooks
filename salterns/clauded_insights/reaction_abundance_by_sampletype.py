#!/usr/bin/env python
"""
Subselect the per-reaction / per-pathway read-abundance %s by SAMPLE TYPE
(Unrestored = R1+R2, Restored = SF2, Reference = R2A) and compare to the manual
attempt in Sheet2 of cliff_mags.xlsx, which computed the same pathway split by a
different method.

Outputs:
  read_abundance_by_reaction_totMAGs_byType.csv   -- my %, aggregated by the 3 types
  reaction_abundance_vs_sheet2.csv                -- side-by-side MePn/Bet comparison

Method note (the crux of the arbitration):
  MINE   : for each reaction, % of the mapped community's READ ABUNDANCE carried by
           MAGs that possess it (denominator = whole mapped community; MAG-abundance
           weighted; presence from the curated MePn table + all MAG models).
  Sheet2 : GENE COUNTS per function, normalised to the MePn+Bet gene pool only
           (MePn% + Bet% = 100%); gene-copy weighted, not abundance weighted.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd

REPO = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent


def sample_type(col: str) -> str:
    depth = "D1" if col.endswith("_D1") else "D2" if col.endswith("_D2") else None
    if depth is None:
        return "other"                       # (H2O etc. — not counted; only D1/D2)
    if col.startswith("R2A"):            t = "Reference_R2A"
    elif col.startswith("SF2"):          t = "Restored_SF2"
    elif col.startswith(("R1_", "R2_")): t = "Unrestored_R1R2"
    else:                                return "other"
    return f"{t}_{depth}"


# ---- 1) my %s, aggregated by sample type -----------------------------------------
mine = pd.read_csv(HERE / "read_abundance_by_reaction_totMAGs.csv", index_col=0)
sample_cols = [c for c in mine.columns if c not in ("n_MAGs", "MEAN_across_samples", "pathway")]
groups = {}
for c in sample_cols:
    groups.setdefault(sample_type(c), []).append(c)
bytype = pd.DataFrame({t: mine[cols].mean(axis=1) for t, cols in groups.items() if t != "other"})
bytype.insert(0, "n_MAGs", mine["n_MAGs"])
bytype.insert(0, "pathway", mine["pathway"])
order = ["pathway", "n_MAGs",
         "Unrestored_R1R2_D1", "Unrestored_R1R2_D2",
         "Restored_SF2_D1", "Restored_SF2_D2",
         "Reference_R2A_D1", "Reference_R2A_D2"]
bytype = bytype[order]
bytype.to_csv(HERE / "read_abundance_by_reaction_totMAGs_byType.csv")
print("== MY % (total % abundance of possessing MAGs), by sample type ==")
print(bytype.round(3).to_string())
print(f"\n(samples per type: " + ", ".join(f"{t}={len(v)}" for t, v in groups.items() if t!='other') + ")")

# ---- 2) Sheet2's values, same sediment samples, by type --------------------------
s2 = pd.ExcelFile(REPO / "cliff_mags.xlsx").parse("Sheet2", header=None)
scols = {c: s2.iloc[3, c] for c in range(1, 35) if isinstance(s2.iloc[3, c], str) and s2.iloc[3, c].endswith("_MG")}
lab = {i: s2.iloc[i, 0] for i in range(s2.shape[0]) if isinstance(s2.iloc[i, 0], str)}
S = pd.DataFrame({lab[i]: {scols[c]: pd.to_numeric(s2.iloc[i, c], errors="coerce") for c in scols}
                  for i in lab}).T
# keep only D1/D2 sediment samples (drop H2O) to match my 24-sample set
sed = [s for s in S.columns if ("_D1_" in s or "_D2_" in s)]
def s2type(s):
    return "Reference_R2A" if s.startswith("R2A") else "Restored_SF2" if s.startswith("SF2") else \
           "Unrestored_R1R2" if s.startswith(("R1_","R2_")) else "other"
s2groups = {}
for s in sed:
    t = s2type(s)
    if t != "other": s2groups.setdefault(t, []).append(s)
def s2mean(rows):
    return {t: S.loc[rows, cols].mean(axis=1).sum() if isinstance(rows,list) else S.loc[rows, cols].mean()
            for t, cols in s2groups.items()}

print("\n== Sheet2 (gene-count method), by sample type ==")
s2summary = pd.DataFrame({t: S.loc[["methylphosphonate","betaine","total","MePn %","Bet %",
                                    "MePn degradation %","Bet → TMA %","Methanogenesis %"], cols].mean(axis=1)
                          for t, cols in s2groups.items()})
print(s2summary.round(3).to_string())
print(f"(Sheet2 sediment samples per type: " + ", ".join(f"{t}={len(v)}" for t,v in s2groups.items()) + ")")

# ---- 3) side-by-side pathway comparison ------------------------------------------
# type-level (mean over both depths) for the Sheet2 comparison
typelevel = pd.DataFrame({
    "Unrestored_R1R2": bytype[["Unrestored_R1R2_D1", "Unrestored_R1R2_D2"]].mean(axis=1),
    "Restored_SF2":    bytype[["Restored_SF2_D1", "Restored_SF2_D2"]].mean(axis=1),
    "Reference_R2A":   bytype[["Reference_R2A_D1", "Reference_R2A_D2"]].mean(axis=1),
})
mp = typelevel.loc["MePn_pathway_union"]
bp = typelevel.loc["Bet_pathway_union"]
cmp = pd.DataFrame({
    "MINE_MePn_%_of_community": mp,
    "MINE_Bet_%_of_community":  bp,
    "MINE_MePn_share_of_(MePn+Bet)_%": 100*mp/(mp+bp),
    "SHEET2_MePn_%_of_pool": s2summary.loc["MePn %"]*100,
    "SHEET2_Bet_%_of_pool":  s2summary.loc["Bet %"]*100,
})
cmp.to_csv(HERE / "reaction_abundance_vs_sheet2.csv")
print("\n== COMPARISON: MePn vs Bet split, mine vs Sheet2 ==")
print(cmp.round(1).to_string())
print("\nSaved read_abundance_by_reaction_totMAGs_byType.csv and reaction_abundance_vs_sheet2.csv")
