#!/usr/bin/env python
"""Relative abundance of the 585 final-set species representatives in the 39 salterns metagenomes (2026-11-20).

Inputs: coverm/<sample>.tsv (CoverM 0.7.0 genome mode on poplar, run_coverm.sh v2: minimap2-sr against one index of the
585 genomes, >= 95 % identity, >= 75 % aligned, min-covered-fraction 0; read counts) and qc_read_totals.tsv
(read_counts.json['qtrim']['result_reads'] = the reads CoverM mapped, i.e. the total it reports).
Relative abundance = genome read count / total QC'd reads, with an 'unbinned_reads' row = 1 - sum, i.e. the basis
of data/Cliff_310MAG_relabund.tsv (whose counts + unbinned_reads sum exactly to its denominator).
Writes final585_read_count.tsv (39 samples), final585_relabund.tsv (39 samples) and final585_relabund_24sediment.tsv
(the 24 D1/D2 columns, in the column order of Cliff_310MAG_relabund.tsv).
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
S = os.path.dirname(os.path.dirname(HERE))
tot = pd.read_csv(os.path.join(HERE, "qc_read_totals.tsv"), sep="\t", header=None, index_col=0)[1]
counts = {}
for s in tot.index:
    t = pd.read_csv(os.path.join(HERE, "coverm", f"{s}.tsv"), sep="\t", index_col=0)
    c = [x for x in t.columns if x.endswith("Read Count")]
    assert len(c) == 1 and len(t) == 585, (s, t.shape)
    counts[s] = t[c[0]]
rc = pd.DataFrame(counts)
assert (rc.sum() <= tot).all()
rc.loc["unbinned_reads"] = tot - rc.sum()
ra = rc / tot
assert ((ra.sum() - 1).abs() < 1e-12).all()
rc.astype(int).to_csv(os.path.join(HERE, "final585_read_count.tsv"), sep="\t")
ra.to_csv(os.path.join(HERE, "final585_relabund.tsv"), sep="\t")
cols = list(pd.read_csv(os.path.join(S, "data", "Cliff_310MAG_relabund.tsv"), sep="\t", index_col=0, nrows=1).columns)
ra[cols].to_csv(os.path.join(HERE, "final585_relabund_24sediment.tsv"), sep="\t")
print("binned fraction per sample (24 sediment):")
print((1 - ra.loc["unbinned_reads", cols]).round(3).to_string())
