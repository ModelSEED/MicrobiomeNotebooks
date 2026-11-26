#!/usr/bin/env python
"""Did the MAG-source assembly (iii) see the same reads as the 2026 assembly? Two checks that need no read mapping.

1. Final k of (iii). In SPAdes output, contigs that are adjacent in the assembly graph overlap by exactly the final k.
   For each candidate k (15..151) count contig-end k-mers (first/last k bp, canonical) shared by two contig ends of
   the sample's MAG contigs; the true k stands out far above the background. Run on the original MAG contigs and,
   as a positive control, on the 2026 MAG contigs of the same sample (known final k = 127).
2. k-mer coverage of the same genome. SPAdes names carry cov = k-mer coverage at the final k, which for read depth d
   and read length L is d * (L - k + 1) / L. For original MAGs and 2026 MAGs of the same species cluster
   (comparison/mags_old_vs_new.tsv, skani ANI >= 95 %) from the SAME sample's own assembly, both at k = 127, the
   ratio of their length-weighted median contig cov is ~ (d_iii / d_2026) * (f_iii / f_2026), f = (L - 126) / L.
   Same reads and read lengths -> ratio ~ 1 (1.0-1.15 allowing for the 2026 trimming: mean read length 148.3 vs
   <= 151). A ratio near 2 would mean (iii) was assembled from about twice the reads.
Writes final_k_overlaps.tsv, kmer_cov_pairs.tsv
Usage: ~/Documents/py_venv/bin/python -I kmer_coverage_check.py
"""
import collections, csv, glob, os, re, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SALTERNS = os.path.normpath(os.path.join(ROOT, ".."))
NODE = re.compile(r">NODE_(\d+)_length_(\d+)_cov_([\d.]+)")
RC = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(RC)[::-1]


def read_fa(files):
    out = []
    for f in files:
        name, buf = None, []
        for line in open(f):
            if line.startswith(">"):
                if name:
                    out.append((name, "".join(buf)))
                name, buf = line.strip(), []
            else:
                buf.append(line.strip().upper())
        if name:
            out.append((name, "".join(buf)))
    return out


def overlap_profile(seqs):
    prof = {}
    for k in range(15, 152):
        c = collections.Counter()
        for _, s in seqs:
            if len(s) >= 2 * k:
                for km in (s[:k], s[-k:]):
                    c[min(km, rc(km))] += 1
        prof[k] = sum(v for v in c.values() if v >= 2)
    return prof


def wmedian_cov(files):
    pts = []
    for f in files:
        for line in open(f):
            if line.startswith(">"):
                m = NODE.match(line)
                if m and int(m.group(2)) >= 2500:
                    pts.append((float(m.group(3)), int(m.group(2))))
    pts.sort()
    tot, acc = sum(l for _, l in pts), 0
    for c, l in pts:
        acc += l
        if acc >= tot / 2:
            return c, len(pts)
    return None, 0


def main():
    orig_files = lambda s: glob.glob(os.path.join(SALTERNS, "MAG_contigs",  # noqa: E731
                                                  f"Salt_Pond_MetaG*{s.replace('R2_restored_D_shore', 'R2_restored_DShore')}_MG_DASTool*.fna"))
    # 1. final k, a spread of samples (sediment/water, rich/poor)
    rows = []
    for s in ("R1_A_D1", "R1_B_H2O", "R2_B_D2", "R2A_B_D1", "SF2_C_D1", "SF2_C_H2O", "R2_restored_D_shore"):
        for which, files in (("original (iii)", orig_files(s)),
                             ("2026", glob.glob(os.path.join(ROOT, "work", "checkm", s, "fasta", "*.fna")))):
            seqs = read_fa(files)
            prof = overlap_profile(seqs)
            best = max(prof, key=prof.get)
            bg = st.median(prof.values())
            rows.append({"sample": s, "assembly": which, "n_contigs": len(seqs), "best_k": best,
                         "shared_ends_at_best_k": prof[best], "median_background": bg,
                         "shared_at_k55": prof[55], "shared_at_k77": prof[77], "shared_at_k99": prof[99],
                         "shared_at_k127": prof[127]})
            print(rows[-1], flush=True)
    with open(os.path.join(HERE, "final_k_overlaps.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # 2. k-mer coverage of the same species in the same sample
    t = list(csv.DictReader(open(os.path.join(ROOT, "comparison", "mags_old_vs_new.tsv")), delimiter="\t"))
    by = collections.defaultdict(lambda: {"old": [], "new": []})
    for r in t:
        side, name = r["genome"].split("|", 1)
        by[(r["sample"], r["species_cluster"])][side].append(name)
    pairs = []
    for (s, cl), d in sorted(by.items()):
        if not d["old"] or not d["new"] or s.endswith("_sediment"):
            continue
        for o in d["old"]:
            co, no = wmedian_cov([os.path.join(SALTERNS, "MAG_contigs", o + ".fna")])
            for n in d["new"]:
                cn, nn = wmedian_cov([os.path.join(ROOT, "work", "checkm", s, "fasta", n + ".fna")])
                if co and cn:
                    pairs.append({"sample": s, "species_cluster": cl, "orig_mag": o, "new_mag": n,
                                  "orig_median_kmer_cov": co, "new_median_kmer_cov": cn,
                                  "orig_contigs_ge2500": no, "new_contigs_ge2500": nn,
                                  "ratio_orig_over_new": round(co / cn, 3)})
    with open(os.path.join(HERE, "kmer_cov_pairs.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pairs[0]), delimiter="\t")
        w.writeheader()
        w.writerows(pairs)
    rat = sorted(p["ratio_orig_over_new"] for p in pairs)
    q = lambda f: rat[int(f * (len(rat) - 1))]  # noqa: E731
    print(f"{len(pairs)} same-sample pairs; ratio orig/new k-mer cov: q10 {q(.1)} q25 {q(.25)} median {q(.5)} "
          f"q75 {q(.75)} q90 {q(.9)}")


if __name__ == "__main__":
    main()
