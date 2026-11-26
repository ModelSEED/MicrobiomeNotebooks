"""Subset the merged DRAM2 outputs to one MAG set, for a distill/product run over exactly that set.
Reads merged/raw-annotations.all.tsv (every chunk) and writes <out>/raw-annotations.tsv with only the MAGs listed (column
1 of a TSV with a 'mag' header), plus <out>/rrnas and <out>/trnas linking those MAGs' per-MAG tables.
Usage: conda_cache/dram2_beta17_distill/bin/python -I bin/subset_merged.py <set.tsv> <out dir>
"""
import os, sys
from pathlib import Path
import pandas as pd
B = Path("/scratch1/afreiburger/wetland_dram2")
want = set(pd.read_csv(sys.argv[1], sep="\t", usecols=["mag"]).mag)
out = Path(sys.argv[2]); (out / "rrnas").mkdir(parents=True, exist_ok=True); (out / "trnas").mkdir(exist_ok=True)
nm = pd.read_csv(B / "inputs" / "name_map.tsv", sep="\t", dtype=str)
missing = want - set(nm.mag)
assert not missing, f"{len(missing)} MAGs of the set were never run: {sorted(missing)[:3]}"
seen, first = set(), True
with open(out / "raw-annotations.tsv", "w") as fh:
    for ch in pd.read_csv(B / "merged" / "raw-annotations.all.tsv", sep="\t", chunksize=500_000, low_memory=False):
        ch = ch[ch.mag.isin(want)]
        seen |= set(ch.mag)
        ch.to_csv(fh, sep="\t", index=False, header=first); first = False
assert seen == want, f"annotations missing for {len(want - seen)} MAGs"
names = set(nm[nm.mag.isin(want)].dram2_name)
n = {"rrnas": 0, "trnas": 0}
for kind in n:
    for f in (B / "merged" / kind).glob("*.tsv"):
        stem = f.name.split("_processed_")[0]
        if stem in names:
            dst = out / kind / f.name
            if not dst.exists():
                os.symlink(os.path.realpath(f), dst)
            n[kind] += 1
print(f"{len(seen)} MAGs; rRNA tables {n['rrnas']}, tRNA tables {n['trnas']} -> {out}")
