#!/usr/bin/env python
"""Rank the pre-scanned bins by raw grdH (K21579) / grdI (K21578) / mtgB / mttB scores (not pass/fail: a Sec-truncated GrdH
scores below its threshold), and check whether each best grdH gene stops at TGA with the reading frame continuing.
Usage: parse_grdscan.py <grd.domtblout> <ko_list> <faa_dir_with_prodigal_headers> <mags_old_vs_new.tsv> <checkm_root> <out.tsv>"""
import gzip, re, sys
from pathlib import Path
import pandas as pd

dom, kol, faadir, table, checkm, out = sys.argv[1:7]
ko = pd.read_csv(kol, sep="\t", index_col=0)
rows = []
for line in open(dom):
    if line.startswith("#"):
        continue
    f = line.split()
    rows.append(dict(gene=f[0], glen=int(f[2]), ko=f[3], hmm_len=int(f[5]), full_e=float(f[6]), full=float(f[7]),
                     dom=float(f[13]), hf=int(f[15]), ht=int(f[16])))
h = pd.DataFrame(rows)
h["bin"] = h.gene.str.split("|").str[0]
h["thr"] = pd.to_numeric(h.ko.map(ko.threshold), errors="coerce")
h["stype"] = h.ko.map(ko.score_type)
h["score"] = h.full.where(h.stype == "full", h.dom)
g = h.groupby(["gene", "ko"]).agg(bin=("bin", "first"), glen=("glen", "first"), score=("score", "max"), thr=("thr", "first"),
                                  full_e=("full_e", "min"), hf=("hf", "min"), ht=("ht", "max"), hmm_len=("hmm_len", "first")).reset_index()
g["hmm_cov"] = (g.ht - g.hf + 1) / g.hmm_len
# a gene counts for the KO it scores best on (by E-value) among the 8 profiles
best = g.sort_values("full_e").drop_duplicates("gene")
t = pd.read_csv(table, sep="\t")
t = t[t.genome.str.startswith("new|")].assign(bin=lambda d: d.genome.str[4:]).set_index("bin")
summ = []
for b, x in best.groupby("bin"):
    r = dict(bin=b)
    for k, name in [("K21579", "grdH"), ("K21578", "grdI"), ("K10672", "grdB"), ("K10670", "grdA"), ("K25220", "mtgB"), ("K14083", "mttB")]:
        y = x[x.ko == k].sort_values("score", ascending=False)
        r[f"{name}_n_sig"] = int((y.score > y.thr).sum())
        r[f"{name}_best"] = round(y.score.iloc[0], 1) if len(y) else None
        r[f"{name}_best_len"] = int(y.glen.iloc[0]) if len(y) else None
        r[f"{name}_best_gene"] = y.gene.iloc[0] if len(y) else ""
    summ.append(r)
s = pd.DataFrame(summ).set_index("bin")
s = s.join(t[["sample", "completeness", "contamination", "phylum", "family", "genus", "n50", "quality"]])


def read_fasta(p):
    seqs, name = {}, None
    for line in open(p):
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        else:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


comp = str.maketrans("ACGTacgt", "TGCAtgca")
# Sec check on the best grdH-scoring gene per bin: stop codon used, and length of the in-frame read-through to the next stop
sec = {}
for b, r in s.iterrows():
    gname = r.grdH_best_gene
    if not gname:
        continue
    hdr = None
    for line in open(Path(faadir) / f"{b}.faa"):
        if line.startswith(">") and line[1:].split()[0] == gname:
            hdr = line; break
    p = [x.strip() for x in hdr[1:].split("#")]
    contig = gname.split("|")[1].rsplit("_", 1)[0]
    beg, end, strand, partial = int(p[1]), int(p[2]), int(p[3]), re.search(r"partial=(\d\d)", hdr)[1]
    fa = list(Path(checkm).glob(f"*/fasta/{b}.fna"))[0]
    seq = read_fasta(fa)[contig]
    if strand == 1:
        stop = seq[end - 3:end]; down = seq[end:]
    else:
        stop = seq[beg - 1:beg + 2].translate(comp)[::-1]; down = seq[:beg - 1].translate(comp)[::-1]
    ext = 0
    for i in range(0, len(down) - 2, 3):
        if down[i:i + 3].upper() in ("TAA", "TAG", "TGA"):
            break
        ext += 1
    sec[b] = dict(grdH_stop=stop.upper(), grdH_partial=partial, grdH_readthrough_codons=ext)
s = s.join(pd.DataFrame(sec).T)
s.sort_values(["grdH_best"], ascending=False).to_csv(out, sep="\t")
pd.set_option("display.width", 250)
print(s.sort_values("grdH_best", ascending=False)[["phylum", "family", "genus", "completeness", "contamination", "grdH_best", "grdH_best_len",
      "grdH_n_sig", "grdH_stop", "grdH_readthrough_codons", "grdI_best", "grdI_n_sig", "grdB_n_sig", "mtgB_n_sig", "mttB_n_sig"]].head(40).to_string())
