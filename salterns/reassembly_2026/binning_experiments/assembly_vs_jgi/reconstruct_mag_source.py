#!/usr/bin/env python
"""Reconstruct the contig-length distribution (contigs >= 2,500 bp) of the non-public assembly the original MAGs
came from (assembly (iii) in provenance/PROVENANCE.md: Tringe-lab Tadpole + SPAdes, ~2018), from the MAG contigs.

Why this works: SPAdes numbers its output NODE_1, NODE_2, ... in order of decreasing length, and the original MAG
contigs keep those names (NODE_<n>_length_<L>_cov_<c>). So each MAG contig is an exact point (rank n, length L) of the
assembly's rank-length curve: exactly n-1 contigs are at least as long. (Checked: in every sample the MAG contigs
have no rank-order violation, and the same holds for all 39 2026 KBase assemblies; see assemblies_2026.tsv.)
Between two sampled ranks r_i < r_j every unsampled contig has a length in [L_j, L_i], which gives hard lower/upper
bounds on Mbp in contigs >= 2,500 bp; the point estimate interpolates linearly.
  - head (ranks before the first MAG contig): counted at L_first (lower bound; no upper bound exists);
  - tail (ranks after the last MAG contig >= 2,500 bp): extrapolated from the local rank density of the MAG contigs
    with L in [2500, 2600]; the lower bound adds nothing there.
Validation: the same reconstruction is run on the 2026 assemblies from their own selected MAG contigs (>= 2,500 bp;
work/checkm/<unit>/fasta), once with all MAG contigs and once subsampled to the same number of points the original
MAGs give for that sample, and compared with the exact values from lengths_2026.tsv.gz.

Inputs: salterns/MAG_contigs/*.fna (all 450 original MAGs; 310 dereplicated + 140 unretained, same per-sample
assemblies), lengths_2026.tsv.gz, reassembly_2026/work/checkm/<unit>/fasta/*.fna.
Writes: mag_source_reconstruction.tsv, reconstruction_validation_2026.tsv
Usage: ~/Documents/py_venv/bin/python -I reconstruct_mag_source.py
"""
import collections, csv, glob, gzip, os, random, re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SALTERNS = os.path.normpath(os.path.join(ROOT, ".."))
NODE = re.compile(r">?NODE_(\d+)_length_(\d+)_cov_")
T = 2500
THRESH = (2500, 5000, 10000, 20000)


def norm_sample(s):
    return "R2_restored_D_shore" if s == "R2_restored_DShore" else s


def read_points(files):
    pts = set()
    for f in files:
        for line in open(f):
            if line.startswith(">"):
                m = NODE.match(line)
                assert m, (f, line[:80])
                pts.add((int(m.group(1)), int(m.group(2))))
    return sorted(pts)


def check_monotone(pts):
    return sum(1 for a, b in zip(pts, pts[1:]) if b[1] > a[1])


def reconstruct(pts, t=T):
    """pts: sorted unique (rank, length), all from one SPAdes output. Returns estimates for contigs >= t."""
    p = [x for x in pts if x[1] >= t]
    if not p:
        return None
    # rank density near the threshold -> tail extrapolation
    near = [x for x in p if x[1] < t + 100]
    if len(near) >= 3 and near[0][1] != near[-1][1]:
        dens = (near[-1][0] - near[0][0]) / (near[0][1] - near[-1][1])   # ranks per bp of length
    else:
        dens = 0.0
    r_last, l_last = p[-1]
    tail_n = int(round(dens * (l_last - t + 1))) if l_last >= t else 0
    tail_n = max(tail_n, 0)
    # integrate: lower/upper/point over ranks 1..r_last, then the tail
    r1, l1 = p[0]
    lo = up = est = (r1 - 1) * l1           # head: >= l1 each (lower bound only)
    curve = [(r, float(l1)) for r in range(1, r1)]   # point-estimate curve for N50 (head at l1)
    for (ra, la), (rb, lb) in zip(p, p[1:]):
        lo += la
        up += la
        est += la
        curve.append((ra, float(la)))
        gap = rb - ra - 1
        if gap > 0:
            lo += gap * lb
            up += gap * la
            # linear interpolation over the gap
            est += sum(la + (lb - la) * (k / (gap + 1)) for k in range(1, gap + 1))
            curve.extend((ra + k, la + (lb - la) * (k / (gap + 1))) for k in range(1, gap + 1))
    lo += l_last
    up += l_last
    est += l_last
    curve.append((r_last, float(l_last)))
    est += tail_n * (t + l_last) / 2
    up += tail_n * l_last
    curve.extend((r_last + k, (t + l_last) / 2) for k in range(1, tail_n + 1))
    lens = [l for _, l in curve]
    n50 = nx(lens, 0.5)
    out = {"n_points": len(p), "first_rank": r1, "first_len": l1, "last_rank": r_last, "last_len": l_last,
           "n_ge_lower": r_last, "n_ge_est": r_last + tail_n, "bp_lower": lo, "bp_est": est, "bp_upper_excl_head": up,
           "N50_est": n50, "L50_est": lx(lens, 0.5)}
    for th in THRESH[1:]:
        q = [x for x in p if x[1] >= th]
        if q:
            out[f"bp_ge{th}_est"] = sum(l for l in lens if l >= th)
            out[f"n_ge{th}_lower"] = q[-1][0]
    return out


def nx(lens, f):
    lens = sorted(lens, reverse=True)
    tot, acc = sum(lens), 0
    for l in lens:
        acc += l
        if acc >= f * tot:
            return int(round(l))
    return 0


def lx(lens, f):
    lens = sorted(lens, reverse=True)
    tot, acc = sum(lens), 0
    for i, l in enumerate(lens, 1):
        acc += l
        if acc >= f * tot:
            return i
    return 0


def exact(lens, t=T):
    l = sorted((x for x in lens if x >= t), reverse=True)
    return {"n": len(l), "bp": sum(l), "N50": nx(l, 0.5)}


def main():
    # ---- original MAG-source assemblies (iii)
    by_sample = collections.defaultdict(list)
    for f in glob.glob(os.path.join(SALTERNS, "MAG_contigs", "*.fna")):
        m = re.match(r"Salt_Pond_MetaG_?(.+?)_MG_DASTool", os.path.basename(f))
        assert m, f
        by_sample[norm_sample(m.group(1))].append(f)
    assert len(by_sample) == 39, len(by_sample)
    rows = []
    npts = {}
    for s in sorted(by_sample):
        pts = read_points(by_sample[s])
        assert check_monotone(pts) == 0, s
        r = reconstruct(pts)
        npts[s] = r["n_points"]
        rows.append({"sample": s, "n_mag_files": len(by_sample[s]), **r})
    with open(os.path.join(HERE, "mag_source_reconstruction.tsv"), "w", newline="") as fh:
        keys = list(dict.fromkeys(k for r in rows for k in r))
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"reconstructed {len(rows)} MAG-source assemblies")

    # ---- validation on the 2026 assemblies, where the truth is known
    lens = collections.defaultdict(list)
    with gzip.open(os.path.join(HERE, "lengths_2026.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            u, _, l = line.rstrip("\n").split("\t")
            lens[u].append(int(l))
    rnd = random.Random(20261111)
    val = []
    for s in sorted(by_sample):
        files = glob.glob(os.path.join(ROOT, "work", "checkm", s, "fasta", "*.fna"))
        if not files:
            continue
        pts = read_points(files)
        assert check_monotone(pts) == 0, s
        truth = exact(lens[s])
        p25 = [x for x in pts if x[1] >= T]
        for mode in ("all_2026_MAG_contigs", "subsampled_to_original_n"):
            q = p25 if mode.startswith("all") else sorted(rnd.sample(p25, min(len(p25), npts[s])))
            r = reconstruct(q)
            val.append({"unit": s, "mode": mode, "n_points": r["n_points"], "true_n_ge2500": truth["n"],
                        "est_n_ge2500": r["n_ge_est"], "true_bp_ge2500": truth["bp"], "lower_bp": r["bp_lower"],
                        "est_bp": round(r["bp_est"]), "upper_bp_excl_head": r["bp_upper_excl_head"],
                        "est_over_true_bp": round(r["bp_est"] / truth["bp"], 4),
                        "lower_over_true_bp": round(r["bp_lower"] / truth["bp"], 4),
                        "true_N50_ge2500": truth["N50"], "est_N50_ge2500": r["N50_est"],
                        "est_over_true_N50": round(r["N50_est"] / truth["N50"], 4)})
    with open(os.path.join(HERE, "reconstruction_validation_2026.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(val[0]), delimiter="\t")
        w.writeheader()
        w.writerows(val)
    for mode in ("all_2026_MAG_contigs", "subsampled_to_original_n"):
        v = [x for x in val if x["mode"] == mode]
        eb = sorted(x["est_over_true_bp"] for x in v)
        lb = sorted(x["lower_over_true_bp"] for x in v)
        en = sorted(x["est_over_true_N50"] for x in v)
        print(mode, "units", len(v), "est/true bp: min %.3f median %.3f max %.3f" % (eb[0], eb[len(eb) // 2], eb[-1]),
              "| lower/true bp: min %.3f median %.3f max %.3f" % (lb[0], lb[len(lb) // 2], lb[-1]),
              "| est/true N50: min %.3f median %.3f max %.3f" % (en[0], en[len(en) // 2], en[-1]))


if __name__ == "__main__":
    main()
