#!/usr/bin/env python3
"""Smoke test before the launch: the first 1,000 proteins of kofam_rollout chunk 01, searched with this track's
hmmsearch command (--cpu 0, tblout only, -o /dev/null) at the kofam_rollout -Z (2,732,576), must reproduce
(a) kofam_rollout's tblout rows for those proteins (all fields except the description) and
(b) kofam_rollout's best-KO calls (results/all_genes_ko.tsv.gz) after 02_postprocess_chunk.py.
usage: 00_smoke_compare.py <smoke_W> <kofam_rollout_W>
"""
import gzip
import sys

SW, KW = sys.argv[1], sys.argv[2]
names = {l[1:].split()[0] for l in open(f"{SW}/chunks/chunk_9999.faa") if l.startswith(">")}


def rows(fh):
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split(None, 18)
        if f[0] in names:
            yield tuple(f[:18])


new = sorted(rows(open(f"{SW}/hmm/chunk_9999.tblout")))
old = sorted(rows(gzip.open(f"{KW}/hmm/chunk_01.tblout.gz", "rt")))
print("tblout rows new", len(new), "old", len(old), "identical", new == old)
if new != old:
    sn, so = set(new), set(old)
    print("only new", len(sn - so), list(sn - so)[:3])
    print("only old", len(so - sn), list(so - sn)[:3])
# calls
oc = {}
with gzip.open(f"{KW}/results/all_genes_ko.tsv.gz", "rt") as fh:
    h = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        r = dict(zip(h, line.rstrip("\n").split("\t")))
        if r["gene"] in names:
            oc[r["gene"]] = (r["best_ko"], r["best_score"], r["all_passing"])
nc = {}
with gzip.open(f"{SW}/results/chunks/chunk_9999.calls.tsv.gz", "rt") as fh:
    h = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        r = dict(zip(h, line.rstrip("\n").split("\t")))
        if r["best_ko"]:
            nc[r["rep"]] = (r["best_ko"], r["best_score"], r["all_passing"])
print("genes with passing KO new", len(nc), "old", len(oc), "identical (best KO, score, all passing)", nc == oc)
for g in sorted(set(nc) | set(oc)):
    if nc.get(g) != oc.get(g):
        print("DIFF", g, nc.get(g), oc.get(g))
print(open(f"{SW}/results/chunks/chunk_9999.qc.txt").read())
