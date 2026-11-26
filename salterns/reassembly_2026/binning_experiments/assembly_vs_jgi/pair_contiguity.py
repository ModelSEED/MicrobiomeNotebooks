#!/usr/bin/env python
"""Same genome, same sample, two assemblies: does the 2026 contig-length deficit depend on the genome's coverage?
Input: kmer_cov_pairs.tsv (from kmer_coverage_check.py: 190 original/2026 MAG pairs of one species cluster from the
same sample's own assembly). For each MAG, N50 and bp over its contigs >= 2,500 bp only (so MaxBin2's 1-2.5 kb
contigs do not bias the 2026 side). Pairs are split into quartiles of the 2026 k-mer coverage (k = 127).
Writes pair_contiguity.tsv (per pair) and prints the quartile table.
Usage: ~/Documents/py_venv/bin/python -I pair_contiguity.py
"""
import csv, os, re, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SALTERNS = os.path.normpath(os.path.join(ROOT, ".."))
LEN = re.compile(r"_length_(\d+)_")


def n50_ge2500(path):
    lens = sorted((int(m.group(1)) for line in open(path) if line.startswith(">")
                   for m in [LEN.search(line)] if m and int(m.group(1)) >= 2500), reverse=True)
    tot, acc = sum(lens), 0
    for l in lens:
        acc += l
        if acc >= tot / 2:
            return l, tot, len(lens)
    return 0, 0, 0


def main():
    pairs = list(csv.DictReader(open(os.path.join(HERE, "kmer_cov_pairs.tsv")), delimiter="\t"))
    for p in pairs:
        p["orig_N50_ge2500"], p["orig_bp_ge2500"], _ = n50_ge2500(os.path.join(SALTERNS, "MAG_contigs", p["orig_mag"] + ".fna"))
        p["new_N50_ge2500"], p["new_bp_ge2500"], _ = n50_ge2500(
            os.path.join(ROOT, "work", "checkm", p["sample"], "fasta", p["new_mag"] + ".fna"))
        p["new_over_orig_N50"] = round(p["new_N50_ge2500"] / p["orig_N50_ge2500"], 3)
        p["new_binner"] = p["new_mag"].split(".")[0]
    with open(os.path.join(HERE, "pair_contiguity.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pairs[0]), delimiter="\t")
        w.writeheader()
        w.writerows(pairs)
    pairs.sort(key=lambda p: float(p["new_median_kmer_cov"]))
    n = len(pairs)
    print("quartile of 2026 k-mer cov | n | median 2026 cov | median orig cov | cov ratio | orig N50 | 2026 N50 | "
          "N50 ratio")
    for q in range(4):
        g = pairs[q * n // 4:(q + 1) * n // 4]
        med = lambda k: st.median(float(p[k]) for p in g)  # noqa: E731
        print(f"Q{q + 1} | {len(g)} | {med('new_median_kmer_cov'):.2f} | {med('orig_median_kmer_cov'):.2f} | "
              f"{med('ratio_orig_over_new'):.2f} | {med('orig_N50_ge2500'):.0f} | {med('new_N50_ge2500'):.0f} | "
              f"{med('new_over_orig_N50'):.2f}")
    print("all | median N50 ratio", st.median(p["new_over_orig_N50"] for p in pairs))


if __name__ == "__main__":
    main()
