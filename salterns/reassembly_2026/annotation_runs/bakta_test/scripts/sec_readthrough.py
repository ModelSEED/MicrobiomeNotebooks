#!/usr/bin/env python
"""Sec-aware rescue test on the grd pre-scan: for each bin whose best grdH (K21579) gene ends at TGA with the frame
continuing (grdscan_ranking.tsv), translate from the Prodigal start through the TGA (as U) to the next in-frame stop,
and write the joined proteins plus a 2 kb window around the TGA (for the SECIS cmscan test).
Usage: sec_readthrough.py <grdscan_ranking.tsv> <prodigal faa dir> <checkm root> <out joined.faa> <out windows.fna>"""
import re, sys
from pathlib import Path
import pandas as pd

rank, faadir, checkm, out_faa, out_win = sys.argv[1:6]
CODE = {}
bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i = 0
for a in bases:
    for b in bases:
        for c in bases:
            CODE[a + b + c] = aas[i]; i += 1
comp = str.maketrans("ACGTacgtN", "TGCAtgcaN")


def rf(p):
    s, n = {}, None
    for line in open(p):
        if line.startswith(">"):
            n = line[1:].split()[0]; s[n] = []
        else:
            s[n].append(line.strip())
    return {k: "".join(v).upper() for k, v in s.items()}


s = pd.read_csv(rank, sep="\t", index_col=0)
s = s[(s.grdH_stop == "TGA") & (s.grdH_readthrough_codons >= 40) & (s.grdH_best >= 500)]
with open(out_faa, "w") as fa, open(out_win, "w") as fw:
    for b, r in s.iterrows():
        g = r.grdH_best_gene
        hdr = next(l for l in open(Path(faadir) / f"{b}.faa") if l.startswith(">") and l[1:].split()[0] == g)
        p = [x.strip() for x in hdr[1:].split("#")]
        beg, end, strand = int(p[1]), int(p[2]), int(p[3])
        contig = g.split("|")[1].rsplit("_", 1)[0]
        seq = rf(next(Path(checkm).glob(f"*/fasta/{b}.fna")))[contig]
        if strand == 1:
            nt = seq[beg - 1:]
            tga = end - 3                      # 0-based start of the TGA
            win = seq[max(0, tga - 1000):tga + 1000]
        else:
            nt = seq[:end].translate(comp)[::-1]
            tga = beg - 1
            win = seq[max(0, tga - 1000):tga + 1000].translate(comp)[::-1]
        prot, n_stop = [], 0
        for k in range(0, len(nt) - 2, 3):
            aa = CODE.get(nt[k:k + 3], "X")
            if aa == "*":
                n_stop += 1
                if n_stop == 1 and nt[k:k + 3] == "TGA":
                    prot.append("U"); continue
                break
            prot.append(aa)
        fa.write(f">{b}|joined|{contig}|{beg}-{end}|{strand}|U{len(''.join(prot).split('U')[0]) + 1}\n{''.join(prot)}\n")
        fw.write(f">{b}|{contig}|tga_at_{tga + 1}|strand_{strand}\n{win}\n")
print(len(s), "joined proteins written")
