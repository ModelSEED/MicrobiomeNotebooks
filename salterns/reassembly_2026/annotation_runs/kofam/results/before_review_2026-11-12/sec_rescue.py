#!/usr/bin/env python
"""Selenocysteine (Sec) read-through rescue for the glycine/sarcosine/betaine-reductase selenoproteins.

Prodigal stops at the Sec UGA codon, so a Sec protein is called as an N-terminal piece ending at TGA (+ often a
separate C-terminal piece in the same frame). For GrdH (betaine reductase B beta, K21579) the N-terminal piece
scores just under the KOfam threshold and is often labelled grdB (bakta_test/README.md). This script joins the
pieces through the TGA (translated as U) and re-scores the joined protein.

Sec KOs (marker table; the four Sec subunits of the Grd complexes):
  K21579 grdH (betaine reductase B beta), K10672 grdB (glycine reductase B gamma),
  K21584 grdF (sarcosine reductase B beta), K10670 grdA (Gly/Sar/betaine reductase A; Sec near the N-terminus).

build:  for every Prodigal gene with a hit (score >= 14.7) to a Sec KO (the trigger), look on the DNA for
        (a) a TGA stop codon of the gene with the frame open downstream (junction after the gene), and
        (b) an in-frame TGA as the first stop upstream of the gene's start (junction before the gene; the case of
            grdA, whose ~40-aa N-terminal piece Prodigal often does not call).
        The joined ORF runs from the upstream piece's Prodigal start (or, if that piece was not called, the first
        ATG/GTG/TTG after the previous in-frame stop) through the TGA (U) to the next in-frame stop. The called
        partner gene, if any, is recorded (pair = same strand, same frame, no stop between). One record per
        junction (mag, contig, strand, TGA). Also writes a 603-nt window around each TGA for cmscan (SECIS).
score:  reads hmmsearch of the joined proteins against the 68 marker profiles (local) and against all 27,981
        KOfam profiles (poplar), and cmscan of the windows; applies the DRAM rule and the acceptance gate:
          1. best KO of the joined protein over ALL profiles (significant, lowest full E-value, ties: higher full
             score, KO id) is a Sec KO, and it passes that KO's threshold/score_type;
          2. an aligned domain for that KO spans the U (ali_from <= U-10 and ali_to >= U+10);
          3. the joined score for that KO beats every piece's score for it.
        SECIS (Rfam SECIS_* hit starting -10..+100 nt from the TGA's first base, same strand) is recorded as evidence, not required.
        Status: 'sec_readthrough' = accepted and no piece already had that KO as its best KO (a new call);
        'confirmed' = accepted, a piece already had it; 'rejected' otherwise.

usage: sec_rescue.py build <set_name> <genes.tsv.gz> <marker_hits.tsv.gz> <fasta_map.tsv> <genes_ko.tsv.gz|-> <outdir>
       sec_rescue.py score <outdir> <hmm_info.tsv> <markers.domtblout> <full.tblout.gz> <full.domtblout.gz> <cmscan.tbl>
"""
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

SEC_KOS = {"K21579": "grdH", "K10672": "grdB", "K21584": "grdF", "K10670": "grdA"}
MIN_EXT = 10          # codons that the join must add on the other side of the TGA
FLANK = 10            # aligned residues required on each side of the U
WIN = 300             # nt either side of the TGA for the SECIS window
SECIS_MAX_DIST = 100  # nt from the TGA end to the SECIS start

CODE = {}
_b, _a = "TCAG", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_i = 0
for x in _b:
    for y in _b:
        for z in _b:
            CODE[x + y + z] = _a[_i]; _i += 1
STOPS = {"TAA", "TAG", "TGA"}
STARTS = {"ATG", "GTG", "TTG"}
COMP = str.maketrans("ACGTN", "TGCAN")


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def read_fasta(p, want=None):
    seqs, name, buf = {}, None, []
    with opener(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None and (want is None or name in want):
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None and (want is None or name in want):
        seqs[name] = "".join(buf).upper()
    return seqs


def tr(nt):
    return "".join(CODE.get(nt[k:k + 3], "X") for k in range(0, len(nt) - 2, 3))


def oriented(start, end, strand, L):
    """Prodigal 1-based inclusive coords -> 0-based half-open coords on the gene-oriented strand."""
    return (start - 1, end) if strand == "+" else (L - end, L - start + 1)


def genomic(j, strand, L):
    """oriented 0-based codon start j -> genomic 1-based (lo, hi) of that codon."""
    return (j + 1, j + 3) if strand == "+" else (L - j - 2, L - j)


def build(set_name, genes_p, hits_p, fmap_p, genes_ko_p, outdir):
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    h = pd.read_csv(hits_p, sep="\t", usecols=["gene", "mag", "ko", "score", "full_score", "passes", "gene_best_ko"])
    trig = h[h.ko.isin(SEC_KOS)]
    trig_genes = set(trig.gene)
    trig_mags = set(trig.mag)
    # piece scores for every Sec KO (gene, ko) -> score by score_type (full for all four)
    piece = {(r.gene, r.ko): r.score for r in trig.itertuples()}
    gbest = h.drop_duplicates("gene").set_index("gene").gene_best_ko.fillna("").to_dict()
    if genes_ko_p != "-":                                    # full-DB best KO (kofam track) for every gene
        gk = pd.read_csv(genes_ko_p, sep="\t", usecols=["gene", "mag", "best_ko", "best_score"])
        gk = gk[gk.mag.isin(trig_mags)]
        gbest.update(gk.set_index("gene").best_ko.to_dict())
    fmap = pd.read_csv(fmap_p, sep="\t", header=None, names=["mag", "path"]).set_index("mag").path.to_dict()

    cols = ["gene", "mag", "binner", "contig", "contig_len", "start", "end", "strand", "partial"]
    parts = []
    for ch in pd.read_csv(genes_p, sep="\t", usecols=cols, chunksize=500_000, dtype={"partial": str}):
        parts.append(ch[ch.mag.isin(trig_mags)])
    g = pd.concat(parts)
    tcont = set(zip(g[g.gene.isin(trig_genes)].mag, g[g.gene.isin(trig_genes)].contig))
    g = g[[(m, c) in tcont for m, c in zip(g.mag, g.contig)]]
    print(f"[{set_name}] trigger genes {len(trig_genes)} in {len(trig_mags)} MAGs, {len(tcont)} contigs, "
          f"{len(g)} genes on them")

    recs, seen = [], {}
    n = dict(trigger=len(trig_genes), tga_stop=0, upstream_tga=0, partial_skip=0)
    for mag, gm in g.groupby("mag"):
        seqs = read_fasta(fmap[mag], want=set(gm.contig))
        for contig, gc in gm.groupby("contig"):
            seq = seqs[contig]
            L = len(seq)
            assert L == int(gc.contig_len.iloc[0]), (mag, contig, L)
            O = {"+": seq, "-": seq.translate(COMP)[::-1]}
            ori = {r.gene: (*oriented(int(r.start), int(r.end), r.strand, L), r.strand, r.partial)
                   for r in gc.itertuples()}
            by_end = {(v[2], v[1]): k for k, v in ori.items()}          # (strand, oriented end) -> gene
            for gene in gc.gene[gc.gene.isin(trig_genes)]:
                s0, e0, strand, partial = ori[gene]
                o = O[strand]
                p5 = partial[0] == "1" if strand == "+" else partial[1] == "1"   # 5' incomplete (no start)
                p3 = partial[1] == "1" if strand == "+" else partial[0] == "1"   # 3' incomplete (no stop)
                junctions = []
                # (a) the gene itself ends at TGA
                if not p3 and o[e0 - 3:e0] == "TGA":
                    n["tga_stop"] += 1
                    junctions.append(("after", e0 - 3, gene, s0))
                elif p3:
                    n["partial_skip"] += 1
                # (b) the first in-frame stop upstream of the gene's start is TGA
                if not p5:
                    k = s0 - 3
                    while k >= 0 and o[k:k + 3] not in STOPS:
                        k -= 3
                    if k >= 0 and o[k:k + 3] == "TGA":
                        n["upstream_tga"] += 1
                        up = by_end.get((strand, k + 3))
                        if up is not None:
                            junctions.append(("before", k, up, ori[up][0]))
                        else:                                   # piece not called: longest ORF start
                            k2 = k - 3
                            while k2 >= 0 and o[k2:k2 + 3] not in STOPS:
                                k2 -= 3
                            st = next((q for q in range(k2 + 3, k - 2, 3) if o[q:q + 3] in STARTS), None)
                            if st is not None:
                                junctions.append(("before", k, None, st))
                for side, j, up_gene, start in junctions:
                    jkey = (mag, contig, strand, j)
                    if jkey in seen:
                        seen[jkey]["trigger_genes"].add(gene)
                        continue
                    # read through to the next in-frame stop
                    q = j + 3
                    while q + 3 <= L and o[q:q + 3] not in STOPS:
                        q += 3
                    three_partial = q + 3 > L
                    stop_end = q + 3 if not three_partial else q
                    down_gene = by_end.get((strand, stop_end)) if not three_partial else None
                    if down_gene is not None and ori[down_gene][0] <= j:
                        down_gene = None
                    if side == "before" and down_gene is None:
                        down_gene = gene
                    ext_up, ext_down = (j - start) // 3, (q - j - 3) // 3
                    if ext_up < MIN_EXT or ext_down < MIN_EXT:
                        continue
                    prot = tr(o[start:j]) + "U" + tr(o[j + 3:q])
                    assert "*" not in prot, (mag, contig, j)
                    upos = (j - start) // 3 + 1
                    lo, hi = genomic(j, strand, L)
                    gap = (ori[down_gene][0] - (j + 3)) if down_gene else None
                    pieces = [x for x in (up_gene, down_gene) if x]
                    r = dict(set=set_name, mag=mag, contig=contig, contig_len=L, binner=gc.binner.iloc[0],
                             strand=strand, tga_lo=lo, tga_hi=hi, up_gene=up_gene or "", down_gene=down_gene or "",
                             up_called=up_gene is not None, down_called=down_gene is not None,
                             gap_tga_to_down_start_nt=gap, joined_len=len(prot), u_pos=upos,
                             ext_upstream_aa=ext_up, ext_downstream_aa=ext_down, three_prime_partial=three_partial,
                             trigger_genes={gene},
                             pieces_best_ko=";".join(f"{x.split('|', 1)[1]}:{gbest.get(x, '')}" for x in pieces),
                             _pieces=pieces, _prot=prot,
                             _win=O[strand][max(0, j - WIN):j + 3 + WIN], _woff=min(j, WIN))
                    for ko, sym in SEC_KOS.items():
                        sc = [piece.get((x, ko)) for x in pieces]
                        sc = [v for v in sc if v is not None and v == v]
                        r[f"piece_max_{sym}"] = max(sc) if sc else None
                    seen[jkey] = r
                    recs.append(r)
    for i, r in enumerate(recs, 1):
        r["junction_id"] = f"{set_name}_sec{i:04d}"
        r["trigger_genes"] = ";".join(sorted(x.split("|", 1)[1] for x in r["trigger_genes"]))
    with open(out / f"{set_name}.joined.faa", "w") as fa, open(out / f"{set_name}.windows.fna", "w") as fw:
        for r in recs:
            fa.write(f">{r['junction_id']}\n{r['_prot']}\n")
            fw.write(f">{r['junction_id']} tga_offset={r['_woff']}\n{r['_win']}\n")
    d = pd.DataFrame(recs)
    d["pieces"] = d._pieces.map(lambda x: ";".join(x))
    d["tga_offset_in_window"] = d._woff
    d.drop(columns=[c for c in d.columns if c.startswith("_")]).to_csv(out / f"{set_name}.junctions.tsv", sep="\t",
                                                                        index=False)
    n["junctions"] = len(recs)
    n["junctions_with_called_partner_pair"] = int((d.up_called & d.down_called).sum()) if len(d) else 0
    with open(out / f"{set_name}.build_counts.txt", "w") as fo:
        for k, v in n.items():
            fo.write(f"{k}\t{v}\n")
    print(f"[{set_name}]", n)


def read_dom(p):
    rows = []
    with opener(p) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split(None, 22)
            rows.append(dict(q=f[0], ko=f[3], full_e=float(f[6]), full=float(f[7]), dom=float(f[13]),
                             hf=int(f[15]), ht=int(f[16]), ali_from=int(f[17]), ali_to=int(f[18])))
    return pd.DataFrame(rows)


def read_tbl(p):
    rows = []
    with opener(p) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split(None, 18)
            rows.append(dict(q=f[0], ko=f[2], full_e=float(f[4]), full=float(f[5]), dom_best=float(f[8])))
    return pd.DataFrame(rows)


def score(outdir, info_p, mk_dom_p, full_tbl_p, full_dom_p, cm_p):
    out = Path(outdir)
    info = pd.read_csv(info_p, sep="\t", dtype={"threshold": str}).set_index("ko")
    j = pd.concat([pd.read_csv(p, sep="\t") for p in sorted(out.glob("*.junctions.tsv"))], ignore_index=True)
    j = j.set_index("junction_id")

    def significant(t):
        t = t.join(info[["threshold", "score_type", "leng"]], on="ko")
        num = t.score_type.isin(["full", "domain"])
        t["score"] = t.full.where(t.score_type == "full", t.dom_best)
        t["thr"] = pd.to_numeric(t.threshold.where(num), errors="coerce")
        t["sig"] = num & (t.score > t.thr)
        return t

    ft = significant(read_tbl(full_tbl_p))
    fd = read_dom(full_dom_p)
    md = read_dom(mk_dom_p)
    # best KO over all profiles (DRAM rule)
    s = ft[ft.sig].sort_values(["q", "full_e", "full", "ko"], ascending=[True, True, False, True])
    best = s.drop_duplicates("q").set_index("q")
    j["best_ko_all"] = best.ko
    j["best_ko_all_score"] = best.score.round(1)
    j["best_ko_all_threshold"] = best.thr
    # best Sec KO hit (by score) and second-best KO
    for ko, sym in SEC_KOS.items():
        x = ft[ft.ko == ko].set_index("q")
        j[f"joined_{sym}"] = x.score.round(1)
    j["all_sig_kos"] = s.groupby("q").ko.apply(lambda v: ";".join(v))       # in DRAM best-KO order
    # best among the 68 markers (local marker search) for cross-check
    mdt = md.groupby(["q", "ko"]).agg(full_e=("full_e", "min"), full=("full", "max"), dom_best=("dom", "max")).reset_index()
    mdt = significant(mdt)
    mb = mdt[mdt.sig].sort_values(["q", "full_e", "full", "ko"], ascending=[True, True, False, True]).drop_duplicates("q")
    j["best_ko_markers"] = mb.set_index("q").ko
    # U spanned by an aligned domain of the called KO
    span = {}
    for (q, ko), x in fd.groupby(["q", "ko"]):
        u = j.at[q, "u_pos"] if q in j.index else None
        if u is None:
            continue
        span[(q, ko)] = bool(((x.ali_from <= u - FLANK) & (x.ali_to >= u + FLANK)).any())
    j["u_spanned"] = [span.get((q, ko), False) for q, ko in zip(j.index, j.best_ko_all.fillna(""))]
    # SECIS evidence from cmscan --tblout (window strand '+' = gene strand)
    secis = {}
    with open(cm_p) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            model, q, strand, sfrom, sto, sc, ev = f[0], f[2], f[9], int(f[7]), int(f[8]), float(f[14]), float(f[15])
            if "SECIS" not in model.upper() or strand != "+" or q not in j.index:
                continue
            dist = sfrom - (int(j.at[q, "tga_offset_in_window"]) + 1)      # nt from the TGA's first base
            if -10 <= dist <= SECIS_MAX_DIST and (q not in secis or sc > secis[q][1]):
                secis[q] = (model, sc, ev, dist)
    j["secis_model"] = [secis[q][0] if q in secis else "" for q in j.index]
    j["secis_score"] = [secis[q][1] if q in secis else None for q in j.index]
    j["secis_evalue"] = [secis[q][2] if q in secis else None for q in j.index]
    j["secis_dist_nt"] = [secis[q][3] if q in secis else None for q in j.index]
    # gate
    is_sec = j.best_ko_all.isin(SEC_KOS)
    called_sym = j.best_ko_all.map(SEC_KOS)
    joined_sc = [j.at[q, f"joined_{s}"] if isinstance(s, str) else None for q, s in zip(j.index, called_sym)]
    piece_sc = [j.at[q, f"piece_max_{s}"] if isinstance(s, str) else None for q, s in zip(j.index, called_sym)]
    j["called_marker"] = called_sym.fillna("")
    j["joined_score_called"] = joined_sc
    j["piece_max_score_called"] = piece_sc
    better = [(a is not None and a == a) and (b is None or b != b or a > b) for a, b in zip(joined_sc, piece_sc)]
    j["beats_pieces"] = better
    acc = is_sec & j.u_spanned & pd.Series(better, index=j.index)
    already = [any(p.endswith(":" + ko) for p in str(pb).split(";")) if isinstance(ko, str) else False
               for pb, ko in zip(j.pieces_best_ko, j.best_ko_all)]
    j["piece_already_called"] = already
    j["status"] = "rejected"
    j.loc[acc & ~pd.Series(already, index=j.index), "status"] = "sec_readthrough"
    j.loc[acc & pd.Series(already, index=j.index), "status"] = "confirmed"
    j.reset_index().to_csv(out / "sec_readthrough_junctions_scored.tsv", sep="\t", index=False)
    print(j.groupby(["set", "status", "called_marker"]).size().to_string())


if __name__ == "__main__":
    if sys.argv[1] == "build":
        build(*sys.argv[2:8])
    elif sys.argv[1] == "score":
        score(*sys.argv[2:8])
    else:
        sys.exit(__doc__)
