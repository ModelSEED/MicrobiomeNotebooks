#!/usr/bin/env python
"""Step 10d (2026-11-25): genus of the methanol:cobalamin methyltransferase (mtaB) and its corrinoid-protein partner
(mtaA) reads, which 10b did not cover. Same rule as 10b: genus of the best-hit representative (refs/rep_info.tsv.gz),
identity bin of the best hit. All 39 samples. Output results/methanol_taxonomy.tsv (sample, label, domain, phylum,
family, genus, genus_source, ident_bin, reads). Run on poplar, 1 CPU: python scripts/10d_methanol_taxonomy.py
"""
import os
import pandas as pd

W = "/scratch1/afreiburger/salterns_2026/read_markers"
R = "/scratch1/afreiburger/salterns_2026/rollout"
os.chdir(W)
samples = sorted(os.listdir(f"{R}/reads"))
ri = pd.read_csv("refs/rep_info.tsv.gz", sep="\t", dtype=str).fillna("").set_index("rep")
LABELS = {"mtaB", "mtaA"}


def ibin(p):
    return ">=95" if p >= 95 else ("85-95" if p >= 85 else ("70-85" if p >= 70 else "<70"))


rows = []
for s in samples:
    rd = pd.read_csv(f"results/per_sample/{s}.reads.tsv.gz", sep="\t", dtype={"read": str, "mate": str})
    rd = rd[rd.label.isin(LABELS)].join(ri[["genus", "family", "phylum", "domain", "genus_source"]], on="rep")
    rd["ident_bin"] = rd.pident.map(ibin)
    for key, gg in rd.groupby(["label", "domain", "phylum", "family", "genus", "genus_source", "ident_bin"], dropna=False):
        rows.append(dict(zip(["label", "domain", "phylum", "family", "genus", "genus_source", "ident_bin"], key),
                         sample=s, reads=len(gg)))
    print(s, flush=True)
pd.DataFrame(rows).to_csv("results/methanol_taxonomy.tsv", sep="\t", index=False)
