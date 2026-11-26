#!/usr/bin/env python3
"""Threshold test: share of each sample's READS in each contig class (map_reads_tf.sh + compare_reassembly.py).

Reads are mapped once, to the 1 kb re-assembly (its contigs >= 2 kb are the KBase assembly). Classes: binned, unb_2.5,
unb_2.0, unb_1.0 (as in threshold_sensitivity.py), plus "unmapped" (reads on no contig >= 1 kb: shorter contigs,
unassembled sequence, low-quality reads). Shares are of all reads in the fastq files. Also the mean read depth of
each class relative to binned (reads per bp, normalised), the unbiased counterpart of the k-mer coverage ratio.
Inputs: ~/gtdbtk_kbase/reasm/<run>.classes.tsv.gz, <run>.contig_reads.tsv and <run>.reads_total (copied from TF).
Output: unbinned/threshold/reassembly_reads.json
Usage: ~/Documents/py_venv/bin/python -I reads_by_class.py RUN [RUN ...]
"""
import json, os, subprocess, sys
import pandas as pd

RE = os.path.expanduser("~/gtdbtk_kbase/reasm")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "threshold")
p = f"{OUT}/reassembly_reads.json"; summ = json.load(open(p)) if os.path.exists(p) else {}
for run in sys.argv[1:]:
    for f in (f"{run}.contig_reads.tsv", f"{run}.reads_total"):
        if not os.path.exists(f"{RE}/{f}"):
            subprocess.run(["rsync", "-a", f"TF:tf_runs/wetland_reasm/{run}/{f}", f"{RE}/"], check=True)
    c = pd.read_csv(f"{RE}/{run}.classes.tsv.gz", sep="\t", keep_default_na=False)
    r = pd.read_csv(f"{RE}/{run}.contig_reads.tsv", sep="\t", header=None, names=["contig", "reads"])
    tot = int(open(f"{RE}/{run}.reads_total").read())
    c = c.merge(r, on="contig", how="left").fillna({"reads": 0})
    g = c.groupby("class").agg(reads=("reads", "sum"), bp=("length", "sum"))
    o = {"reads_total": tot, "reads_mapped_pct": round(g.reads.sum() / tot * 100, 2)}
    dens_b = g.loc["binned", "reads"] / g.loc["binned", "bp"] if "binned" in g.index else float("nan")
    for k in ("binned", "unb_2.5", "unb_2.0", "unb_1.0"):
        if k in g.index:
            o[f"{k}_reads_pct"] = round(g.loc[k, "reads"] / tot * 100, 2)
            o[f"{k}_depth_rel_binned"] = round((g.loc[k, "reads"] / g.loc[k, "bp"]) / dens_b, 3)
    o["unmapped_pct"] = round(100 - g.reads.sum() / tot * 100, 2)
    o["kbase_2kb_floor_reads_pct"] = round(g.loc[[k for k in ("binned", "unb_2.5", "unb_2.0") if k in g.index], "reads"].sum() / tot * 100, 2)
    o["gain_from_1kb_floor_reads_pct"] = o.get("unb_1.0_reads_pct", 0.0)
    summ[run] = o; print(run, json.dumps(o), flush=True)
json.dump(summ, open(p, "w"), indent=1)
