#!/usr/bin/env python
"""KOfam marker calls on the 310 original (2020, dereplicated) MAGs, with the same gene caller, release and
reporting settings as the 1,324-MAG kofam track, so original-vs-re-assembly guild flips can be attributed.

  * genes: Prodigal 2.6.3 (local gtdbtk_venv build) `-p meta -g 11 -q`; every original contig is >= 2,500 bp
    (MAG_contigs/README.md), so the >= 1,000 bp filter of the kofam track removes nothing;
  * search: the 68 marker profiles of KOfam 2026-10-26 (`../db/markers.hmm`, copied from poplar
    kofam/db/markers.hmm), `hmmsearch -T 14.7 --domT 14.7 -Z 3911146` (= the kofam track's settings);
  * rule: DRAM 1.5 `sig_scores` (score by ko_list score_type > threshold; '-' never called); best KO per gene =
    lowest full E-value, ties by higher full score then KO id -- among the 68 MARKERS only (the kofam track
    took the best over all 27,981 profiles; see RESULTS.md for how often that differs on the new MAGs).

Writes results/originals_310/: genes.tsv.gz (same columns as genes_all.tsv.gz), proteins.faa.gz,
marker_hits.tsv.gz (same columns as all_marker_hits.tsv.gz, gene_best_ko = best among markers),
markers.{tblout,domtblout}.gz.
Run: ~/Documents/py_venv/bin/python -I originals_markers.py   (nice 10, <= 16 threads)
"""
import gzip
import os
import re
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
K = HERE.parent                                        # annotation_runs/kofam
S = K.parents[2]                                       # salterns/
OUT = K / "results/originals_310"
TMP = OUT / "tmp"
BIN = Path.home() / "Documents/gtdbtk_venv/bin"
PAT = re.compile(r"Salt_Pond_MetaG_?(.+?)_MG_DASTool_bins[._](concoct_out|metabat|maxbin)\.(\d+)")
SHORT = 2500
THREADS = 16


def key(name):
    m = PAT.search(str(name))
    return f"{m[1]}|{m[2]}|{int(m[3])}" if m else None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(exist_ok=True)
    t = pd.read_csv(S / "reassembly_2026/comparison/mags_old_vs_new.tsv", sep="\t")
    old = t[t.genome.str.startswith("old|")].genome.str[4:].tolist()
    assert len(old) == 310
    files = {key(p.name): p for p in (S / "MAG_contigs").glob("*.contigs__.fna")}
    fa = {m: files[key(m)] for m in old}
    assert len(fa) == 310 and len(set(fa.values())) == 310

    def prodigal(m):
        faa, gff = TMP / f"{m}.faa", TMP / f"{m}.gff"
        if not faa.exists():
            subprocess.run(["nice", "-n", "10", BIN / "prodigal", "-p", "meta", "-g", "11", "-q", "-i", fa[m],
                            "-a", str(faa) + ".part", "-f", "gff", "-o", gff], check=True)
            os.replace(str(faa) + ".part", faa)
        return m

    with ThreadPoolExecutor(THREADS) as ex:
        list(ex.map(prodigal, old))

    # gene table + pooled proteins (IDs <MAG>|<contig>_<n>, as in the kofam track)
    rows = []
    clen, alen = {}, {}
    with gzip.open(OUT / "proteins.faa.gz", "wt") as fo, open(TMP / "all.faa", "w") as fa_all:
        for m in old:
            for line in open(fa[m]):
                if line.startswith(">"):
                    c = line[1:].split()[0]
                    clen[(m, c)] = 0
                else:
                    clen[(m, c)] += len(line.strip())
            binner = {"concoct_out": "concoct", "metabat": "metabat2", "maxbin": "maxbin2"}[key(m).split("|")[1]]
            name, seq = None, []
            for line in open(TMP / f"{m}.faa"):
                if line.startswith(">"):
                    if name:
                        s = "".join(seq).rstrip("*")
                        alen[name] = len(s)
                        fo.write(f">{name}\n{s}\n"); fa_all.write(f">{name}\n{s}\n")
                    p = [x.strip() for x in line[1:].split("#")]
                    gid = p[0]
                    contig = gid.rsplit("_", 1)[0]
                    name = f"{m}|{gid}"
                    part = re.search(r"partial=(\d\d)", line)[1]
                    rbs = re.search(r"rbs_motif=([^;]+)", line)[1]
                    seq = []
                    rows.append([name, m, binner, contig, clen[(m, contig)], int(p[1]), int(p[2]),
                                 "+" if p[3] == "1" else "-", part, None, rbs])
                else:
                    seq.append(line.strip())
            s = "".join(seq).rstrip("*")
            alen[name] = len(s)
            fo.write(f">{name}\n{s}\n"); fa_all.write(f">{name}\n{s}\n")
            name = None
    g = pd.DataFrame(rows, columns=["gene", "mag", "binner", "contig", "contig_len", "start", "end", "strand",
                                    "partial", "aa_len", "rbs_motif"])
    g["aa_len"] = g.gene.map(alen)
    assert g.contig_len.min() >= 1000, g.contig_len.min()
    g.to_csv(OUT / "genes.tsv.gz", sep="\t", index=False)
    print("originals:", len(old), "MAGs", len(g), "genes; min contig", g.contig_len.min())

    tbl, dom = TMP / "markers.tblout", TMP / "markers.domtblout"
    if not dom.exists():
        subprocess.run(["nice", "-n", "10", BIN / "hmmsearch", "--cpu", str(THREADS), "--noali", "-T", "14.7",
                        "--domT", "14.7", "-Z", "3911146", "--tblout", tbl, "--domtblout", str(dom) + ".part",
                        "-o", "/dev/null", K / "db/markers.hmm", TMP / "all.faa"], check=True)
        os.replace(str(dom) + ".part", dom)

    info = pd.read_csv(K / "results/hmm_info.tsv", sep="\t", dtype={"threshold": str}).set_index("ko")
    mk = pd.read_csv(K / "meta/marker_ko_table.tsv", sep="\t").set_index("ko").gene.to_dict()
    gi = g.set_index("gene")
    # coverage from domtblout
    cov = defaultdict(lambda: [10 ** 9, 0, 0])
    for line in open(dom):
        if line.startswith("#"):
            continue
        f = line.split(None, 17)
        a = cov[(f[0], f[3])]
        a[0], a[1], a[2] = min(a[0], int(f[15])), max(a[1], int(f[16])), a[2] + 1
    hits = []
    for line in open(tbl):
        if line.startswith("#"):
            continue
        f = line.split(None, 18)
        gene, ko = f[0], f[2]
        full_e, full_s, dom_s = float(f[4]), float(f[5]), float(f[8])
        st, thr = info.at[ko, "score_type"], info.at[ko, "threshold"]
        if st in ("full", "domain"):
            score = full_s if st == "full" else dom_s
            sig = score > float(thr)
        else:
            score, sig = None, False
        a = cov.get((gene, ko))
        hc = (a[1] - a[0] + 1) / int(info.at[ko, "leng"]) if a else 0.0
        hits.append([gene, ko, mk[ko], st, thr, score, full_s, dom_s, full_e, round(hc, 2), a[2] if a else 0, int(sig)])
    h = pd.DataFrame(hits, columns=["gene", "ko", "marker", "score_type", "threshold", "score", "full_score",
                                    "best_domain_score", "full_evalue", "hmm_cov", "n_domains", "passes"])
    sig = h[h.passes == 1].sort_values(["gene", "full_evalue", "full_score", "ko"], ascending=[True, True, False, True])
    best = sig.drop_duplicates("gene").set_index("gene").ko
    h["gene_best_ko"] = h.gene.map(best).fillna("")
    h["gene_best_marker"] = h.gene_best_ko.map(mk).fillna("")
    h["is_best_ko"] = (h.gene_best_ko == h.ko).astype(int)
    for c in ("mag", "binner", "contig", "contig_len", "aa_len", "partial"):
        h[c] = h.gene.map(gi[c])
    h["short_maxbin2_contig"] = ((h.binner == "maxbin2") & (h.contig_len < SHORT)).astype(int)
    h = h[["gene", "mag", "binner", "contig", "contig_len", "short_maxbin2_contig", "aa_len", "partial", "ko", "marker",
           "score_type", "threshold", "score", "full_score", "best_domain_score", "full_evalue", "hmm_cov", "n_domains",
           "passes", "is_best_ko", "gene_best_ko", "gene_best_marker"]]
    h.to_csv(OUT / "marker_hits.tsv.gz", sep="\t", index=False)
    for p in (tbl, dom):
        subprocess.run(["gzip", "-kf", p], check=True)
        os.replace(str(p) + ".gz", OUT / (p.name + ".gz"))
    print("marker hits", len(h), "passing", int(h.passes.sum()), "genes with a marker KO", best.size)


if __name__ == "__main__":
    main()
