#!/usr/bin/env python
"""Binning-free check of the coverage dependence: is the 2026 contig-length ceiling a property of the assembly, or
only of what the binners kept?
Contigs >= 2,500 bp are put into classes of their own SPAdes k-mer coverage (k = 127, from the NODE name):
  2026 assembly, ALL contigs (lengths_2026.tsv.gz, 39 per-sample metaSPAdes assemblies)
  2026 MAG contigs only (work/checkm/<unit>/fasta, selected MAGs)
  original (iii) MAG contigs (salterns/MAG_contigs, 450 MAGs); (iii) has ~2.9x the k-mer coverage for the same
     genome (kmer_cov_pairs.tsv), so its class limits are the 2026 limits x 2.9.
For each class: contigs, Mbp, N50, Mbp in contigs >= 20 kb, pooled over the 39 samples.
Writes contiguity_by_coverage.tsv
Usage: ~/Documents/py_venv/bin/python -I contiguity_by_coverage.py
"""
import csv, glob, gzip, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SALTERNS = os.path.normpath(os.path.join(ROOT, ".."))
NODE = re.compile(r">?NODE_(\d+)_length_(\d+)_cov_([\d.]+)")
EDGES = (0, 1.0, 1.5, 2.5, 5.0, float("inf"))      # 2026 k-mer coverage classes
SCALE_III = 2.9                                     # median (iii)/2026 k-mer coverage, same genome & sample


def cls(cov, scale=1.0):
    c = cov / scale
    for lo, hi in zip(EDGES, EDGES[1:]):
        if lo <= c < hi:
            return f"{lo}-{hi}" if hi != float("inf") else f">={lo}"


def summarise(lens):
    lens = sorted(lens, reverse=True)
    tot, acc, n50 = sum(lens), 0, 0
    for l in lens:
        acc += l
        if acc >= tot / 2:
            n50 = l
            break
    return {"contigs": len(lens), "Mbp": round(tot / 1e6, 1), "N50": n50,
            "Mbp_ge20kb": round(sum(l for l in lens if l >= 20000) / 1e6, 1),
            "pct_in_ge20kb": round(100 * sum(l for l in lens if l >= 20000) / tot, 1) if tot else 0}


def headers(files):
    for f in files:
        for line in open(f):
            if line.startswith(">"):
                yield line


def main():
    groups = {}
    with gzip.open(os.path.join(HERE, "lengths_2026.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            u, cid, l = line.rstrip("\n").split("\t")
            m = NODE.match(cid)
            if m and int(l) >= 2500:
                groups.setdefault(("2026 assembly, all contigs", cls(float(m.group(3)))), []).append(int(l))
    for line in headers(glob.glob(os.path.join(ROOT, "work", "checkm", "*", "fasta", "*.fna"))):
        m = NODE.match(line)
        if m and int(m.group(2)) >= 2500:          # co-assembly MAGs (MEGAHIT ids) do not match NODE and drop out
            groups.setdefault(("2026 MAG contigs", cls(float(m.group(3)))), []).append(int(m.group(2)))
    for line in headers(glob.glob(os.path.join(SALTERNS, "MAG_contigs", "*.fna"))):
        m = NODE.match(line)
        groups.setdefault(("original (iii) MAG contigs, cov/2.9", cls(float(m.group(3)), SCALE_III)), []).append(
            int(m.group(2)))
    order = [f"{lo}-{hi}" if hi != float("inf") else f">={lo}" for lo, hi in zip(EDGES, EDGES[1:])]
    rows = []
    for src in ("2026 assembly, all contigs", "2026 MAG contigs", "original (iii) MAG contigs, cov/2.9"):
        for c in order:
            if (src, c) in groups:
                rows.append({"source": src, "kmer_cov_class_2026_scale": c, **summarise(groups[(src, c)])})
                print(rows[-1])
    with open(os.path.join(HERE, "contiguity_by_coverage.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
