#!/usr/bin/env python3
"""Threshold test: check the 1 kb re-assemblies against the KBase assemblies and class their contigs.

For each re-assembled sample (run_reassembly_tf.sh: MEGAHIT 1.2.9, KBase app parameters except --min-contig-len 1000):
 1. Reproduction. The KBase assembly (min_contig_length 2000) should be the re-assembly's contigs >= 2 kb. Checked by
    exact sequence identity (either strand), by count and bp in both directions. A re-assembly contig that is not
    identical is aligned to the KBase assembly (minimap2 2.28 asm5, base-level) and counted as reproduced if it is
    >= 99% identical over >= 95% of its length.
 2. Are the new 1-2 kb contigs new sequence? Same alignment: the share of 1-2 kb contigs (and bp) contained (>= 99% id,
    >= 95% of length) in a KBase contig, i.e. already represented in the 2 kb assembly.
 3. Class per re-assembly contig, as in threshold_sensitivity.py: binned (identical to, or contained in, a contig of
    one of the sample's admitted MAGs), unb_2.5 (>= 2,500 bp), unb_2.0 (2,000-2,499), unb_1.0 (1,000-1,999).
    Shares of contigs, bp and abundance (contig length x MEGAHIT k-mer coverage `multi=` of the re-assembly).
Inputs: ~/gtdbtk_kbase/reasm/<run>.min1000.fa (copied back from TF), ~/gtdbtk_kbase/assemblies/<run>_*.fna,
~/gtdbtk_kbase/genomes/*__<run>.fna.
Outputs: unbinned/threshold/reassembly_summary.json; ~/gtdbtk_kbase/reasm/<run>.classes.tsv.gz (contig, length, multi,
class, kbase_match) and <run>.1to2kb.fa (the 1-2 kb contigs, for Bakta and the assignment pass).
Usage: ~/Documents/py_venv/bin/python -I compare_reassembly.py RUN [RUN ...]
"""
import glob, gzip, hashlib, json, os, re, subprocess, sys
from collections import Counter
import pandas as pd

MM2 = os.path.expanduser("~/Documents/EmilyKin/meta/mag_gene_ab/tools/minimap2-2.28_x64-linux/minimap2")
RE = os.path.expanduser("~/gtdbtk_kbase/reasm"); ASM = os.path.expanduser("~/gtdbtk_kbase/assemblies")
MAGS = os.path.expanduser("~/gtdbtk_kbase/genomes")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "threshold")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
band = lambda L: "unb_2.5" if L >= 2500 else "unb_2.0" if L >= 2000 else "unb_1.0" if L >= 1000 else "short"


def fasta(path):
    name = None; seq = []
    for l in open(path):
        if l[0] == ">":
            if name is not None: yield name, "".join(seq)
            name = l[1:].rstrip("\n"); seq = []
        else: seq.append(l.strip())
    if name is not None: yield name, "".join(seq)


def key(s):
    s = s.upper(); r = s.translate(COMP)[::-1]
    return hashlib.md5(min(s, r).encode()).digest()


summ = json.load(open(f"{OUT}/reassembly_summary.json")) if os.path.exists(f"{OUT}/reassembly_summary.json") else {}
for run in sys.argv[1:]:
    ra = f"{RE}/{run}.min1000.fa"; ka = glob.glob(f"{ASM}/{run}_*.fna")[0]
    binned = set()
    for f in glob.glob(f"{MAGS}/*__{run}.fna"):
        binned |= {l[1:].split()[0] for l in open(f) if l[0] == ">"}
    kb = {}                                   # digest -> (KBase contig, length)
    for n, s in fasta(ka): kb[key(s)] = (n.split()[0], len(s))
    kb_names = {v[0]: v[1] for v in kb.values()}
    rows, hit_kb = [], set()
    for n, s in fasta(ra):
        cid = n.split()[0]; m = re.search(r"multi=([\d.]+)", n); cov = float(m.group(1)) if m else 0.0
        k = kb.get(key(s)); rows.append([cid, len(s), cov, k[0] if k else ""])
        if k: hit_kb.add(k[0])
    d = pd.DataFrame(rows, columns=["contig", "length", "multi", "kbase_match"])
    d["match"] = d.kbase_match.ne("").map({True: "identical", False: ""})
    # align the non-identical contigs to the KBase assembly
    todo = d[d.kbase_match == ""]
    q = f"{RE}/{run}.nonident.fa"
    want = set(todo.contig)
    with open(q, "w") as fo:
        for n, s in fasta(ra):
            if n.split()[0] in want: fo.write(f">{n.split()[0]}\n{s}\n")
    paf = subprocess.run(["bash", "-c", f"{MM2} -x asm5 -c -t 8 --secondary=no {ka} {q} 2>/dev/null"],
                         capture_output=True, text=True, check=True).stdout
    best = {}
    for l in paf.splitlines():
        p = l.split("\t"); ql, qs, qe, mt, bl = int(p[1]), int(p[2]), int(p[3]), int(p[9]), int(p[10])
        ident, cov = mt / bl, (qe - qs) / ql
        if ident >= 0.99 and cov >= 0.95 and cov > best.get(p[0], (0, 0, ""))[0]: best[p[0]] = (cov, ident, p[5])
    for i in todo.index:
        b = best.get(d.at[i, "contig"])
        if b: d.at[i, "kbase_match"] = b[2]; d.at[i, "match"] = "contained"
    d["class"] = [("binned" if km in binned else band(L)) if km else band(L) for km, L in zip(d.kbase_match, d.length)]
    os.remove(q)
    d.to_csv(f"{RE}/{run}.classes.tsv.gz", sep="\t", index=False)
    one = set(d[d.length < 2000].contig)
    with open(f"{RE}/{run}.1to2kb.fa", "w") as fo:
        for n, s in fasta(ra):
            if n.split()[0] in one: fo.write(f">{n}\n{s}\n")
    ge2 = d[d.length >= 2000]; lt2 = d[d.length < 2000]
    kb_bp = sum(kb_names.values()); kb_hit_bp = sum(kb_names[c] for c in hit_kb)
    o = {"kbase_contigs": len(kb_names), "kbase_bp": kb_bp,
         "reasm_contigs": len(d), "reasm_bp": int(d.length.sum()),
         "reasm_ge2kb_contigs": len(ge2), "reasm_ge2kb_bp": int(ge2.length.sum()),
         "reasm_ge2kb_identical_pct_bp": round(ge2[ge2.match == "identical"].length.sum() / ge2.length.sum() * 100, 2),
         "reasm_ge2kb_identical_or_contained_pct_bp": round(ge2[ge2.match != ""].length.sum() / ge2.length.sum() * 100, 2),
         "kbase_bp_with_identical_reasm_contig_pct": round(kb_hit_bp / kb_bp * 100, 2),
         "reasm_1to2kb_contigs": len(lt2), "reasm_1to2kb_bp": int(lt2.length.sum()),
         "reasm_1to2kb_contained_in_kbase_pct_bp": round(lt2[lt2.match != ""].length.sum() / max(1, lt2.length.sum()) * 100, 2)}
    tot = {"contigs": len(d), "bp": d.length.sum(), "abund": (d.length * d.multi).sum()}
    for c in ("binned", "unb_2.5", "unb_2.0", "unb_1.0"):
        x = d[d["class"] == c]
        o[f"{c}_contigs_pct"] = round(len(x) / tot["contigs"] * 100, 2)
        o[f"{c}_bp_pct"] = round(x.length.sum() / tot["bp"] * 100, 2)
        o[f"{c}_abund_pct"] = round((x.length * x.multi).sum() / tot["abund"] * 100, 2)
    lng = d[d.length >= 2500]; o["cov_ratio_1to2_vs_ge2.5"] = round(
        ((lt2.length * lt2.multi).sum() / max(1, lt2.length.sum())) / ((lng.length * lng.multi).sum() / lng.length.sum()), 3)
    o["len_median_1to2"] = int(lt2.length.median()) if len(lt2) else 0
    summ[run] = o; print(run, json.dumps(o, indent=1), flush=True)
json.dump(summ, open(f"{OUT}/reassembly_summary.json", "w"), indent=1)
