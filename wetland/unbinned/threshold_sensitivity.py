#!/usr/bin/env python3
"""Contig-length threshold sensitivity: 2.5 kb (the binners' minimum) vs lower floors, per Andrew 2026-11-25 ("is
2.5kb a reasonable threshold, or would lowering it to 1kb better capture the system?").
The KBase assemblies (narrative 261700) were made by the KBase MEGAHIT app with min_contig_length = 2000 (provenance),
and every assembly's shortest contig is 2,000-2,285 bp: the 1-2 kb band is NOT in them. It is measured on re-assemblies
with a 1 kb floor (reassembly/), with the same classes. Binning floor checked: all 4,086,693 contigs of the 13,164 MAG FASTAs are >= 2,500 bp (shortest
exactly 2,500), so no 1-2.5 kb contig can be in a MAG.
Denominator: assembled sequence in the file (KBase: contigs >= 2 kb), not the sample's reads.

Part 1 (assemblies only, this script): for every downloaded assembly, the share of contigs, bases and approximate
abundance (contig length x MEGAHIT k-mer coverage `multi=`) that falls in
  binned   = contigs of the assembly's admitted MAGs (all >= 2.5 kb by construction)
  unb_2.5  = unbinned contigs >= 2,500 bp          (what the 2.5 kb paradigm adds over MAGs)
  unb_2.0  = contigs 2,000-2,499 bp               (what a 2 kb floor adds; never binnable)
  unb_1.0  = contigs 1,000-1,999 bp               (what a 1 kb floor adds; empty in the KBase assemblies)
plus the median length of the 2.0-2.5 kb band and its coverage relative to the >= 2.5 kb contigs.
Summary: assemblies with and without admitted MAGs are reported apart; shares are pooled over all bases (sum over
assemblies, so a 5 kb assembly does not weigh like a 2.6 Gbp one), with the per-assembly median and IQR beside them.
Abundance shares are not pooled across samples (k-mer coverage units differ with depth): their pooled figure is the
per-assembly share weighted by assembly bp.
Part 2: threshold_assign.py (species assignability per band) and threshold_genes.py (genes, partial genes, functions
and non-coding RNAs per band, from Bakta).
Outputs: unbinned/threshold/per_assembly.tsv and summary.json.
Usage: ~/Documents/py_venv/bin/python -I threshold_sensitivity.py
"""
import glob, json, os, re
from collections import defaultdict
import pandas as pd

ASM = os.path.expanduser("~/gtdbtk_kbase/assemblies"); MAGS = os.path.expanduser("~/gtdbtk_kbase/genomes")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "threshold"); os.makedirs(OUT, exist_ok=True)
binned = defaultdict(set)
for f in glob.glob(f"{MAGS}/*.fna"):
    run = os.path.basename(f)[:-4].split("__")[-1]
    for l in open(f):
        if l.startswith(">"):
            binned[run].add(l[1:].split()[0])
rows = []
for f in sorted(glob.glob(f"{ASM}/*.fna")):
    run = os.path.basename(f).split("_")[0]
    agg = defaultdict(lambda: [0, 0, 0.0]); lens_short = []; mins = []
    has_mags = run in binned; mine = binned.get(run, set())   # read before the loop: binned is a defaultdict
    name = L = cov = None; min_len = 10**12
    def flush():
        if name is None: return
        cls = "binned" if name in mine else ("unb_2.5" if L >= 2500 else "unb_2.0" if L >= 2000 else ("unb_1.0" if L >= 1000 else "short"))
        a = agg[cls]; a[0] += 1; a[1] += L; a[2] += L * cov; mins.append(L)
        if cls == "unb_2.0": lens_short.append(L)
    for l in open(f):
        if l.startswith(">"):
            flush()
            name = l[1:].split()[0]; m = re.search(r"multi=([\d.]+)", l); cov = float(m.group(1)) if m else 0.0; L = 0
        else:
            L += len(l.strip())
    flush(); min_len = min(mins) if mins else 0
    tot = [sum(v[i] for v in agg.values()) for i in range(3)]
    r = {"run": run, "contigs": tot[0], "bp": tot[1], "has_mags": has_mags,
         "len_median_2to2.5": int(pd.Series(lens_short).median()) if lens_short else 0}
    for c in ("binned", "unb_2.5", "unb_2.0", "unb_1.0", "short"):
        for i, k in enumerate(("contigs", "bp", "abund")):
            r[f"{c}_{k}_pct"] = round(agg[c][i] / tot[i] * 100, 2) if tot[i] else 0.0
    long_cov = (agg["binned"][2] + agg["unb_2.5"][2]) / max(1, agg["binned"][1] + agg["unb_2.5"][1])
    r["cov_ratio_2to2.5_vs_long"] = round((agg["unb_2.0"][2] / max(1, agg["unb_2.0"][1])) / max(1e-9, long_cov), 3)
    r["min_contig_bp"] = min_len
    rows.append(r)
d = pd.DataFrame(rows); d.to_csv(f"{OUT}/per_assembly.tsv", sep="\t", index=False)
CL = ("binned", "unb_2.5", "unb_2.0", "unb_1.0")
def block(x):
    o = {"assemblies": len(x), "bp": int(x.bp.sum()), "contigs": int(x.contigs.sum())}
    for c in CL:
        for k in ("contigs", "bp"):
            o[f"{c}_{k}_pct_pooled"] = round(float((x[f"{c}_{k}_pct"] * x[k if k == "bp" else "contigs"]).sum() / x[k if k == "bp" else "contigs"].sum()), 2)
        o[f"{c}_abund_pct_bpweighted"] = round(float((x[f"{c}_abund_pct"] * x.bp).sum() / x.bp.sum()), 2)
        for k in ("contigs", "bp", "abund"):
            q = x[f"{c}_{k}_pct"].quantile([.25, .5, .75]).round(2).tolist(); o[f"{c}_{k}_pct_median_iqr"] = q
    q = x["cov_ratio_2to2.5_vs_long"].quantile([.25, .5, .75]).round(3).tolist(); o["cov_ratio_2to2.5_vs_long_median_iqr"] = q
    o["cov_ratio_2to2.5_vs_long_bpweighted"] = round(float((x["cov_ratio_2to2.5_vs_long"] * x.bp).sum() / x.bp.sum()), 3)
    o["len_median_2to2.5_median"] = float(x["len_median_2to2.5"].median())
    o["min_contig_bp_range"] = [int(x.min_contig_bp.min()), int(x.min_contig_bp.max())]
    return o
summ = {"assemblies_of_781": len(d), "with_mags": block(d[d.has_mags]),
        "without_mags": {"assemblies": int((~d.has_mags).sum()), "bp": int(d[~d.has_mags].bp.sum()),
                         "unb_2.0_bp_pct_pooled": round(float((d[~d.has_mags]["unb_2.0_bp_pct"] * d[~d.has_mags].bp).sum() / max(1, d[~d.has_mags].bp.sum())), 2)}}
json.dump(summ, open(f"{OUT}/summary.json", "w"), indent=1); print(json.dumps(summ, indent=1))
