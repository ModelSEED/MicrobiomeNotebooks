#!/usr/bin/env python3
"""Genome-level check of the tadpole-pilot assemblies: are the ORIGINAL (iii) MAG genomes there, and in how many pieces?

Stdlib only (Python >= 3.8; on poplar: /scratch1/afreiburger/salterns_2026/envs/spades3153/bin/python -I).

Query  = original MAG contigs (every one is a (iii) contig, NODE_n_length_L_cov_C), listed in genome_check_mags.tsv
         (make_mag_list.py): the pilot samples' own MAGs (derep310 + nonderep140) and the 30 'lost' MAGs.
Target = each assembly of the same sample, sequences >= 1 kb (the 2026 KBase assembly was saved at >= 1 kb, so every
         target is cut to the same floor): <root>/<sample>/arm{A,B,C}/spades/scaffolds.fasta and the 2026 KBase
         scaffolds (--kbase-dir/<sample>_kbase2026.fna, fetched with local/fetch_assembly.py, md5-checked).
Mapping: minimap2 -x asm10 -c --secondary=no (PAF, base-level identity = matches / block length). asm10 so that
         95-99 % identity hits (strain differences) are still reported; the 99 % cut is applied here.

Per MAG x target (genome_check_mags_by_target.tsv):
  frac_cov99 / frac_cov95   fraction of MAG bp covered by alignments >= 99 / 95 % identity (block >= 500 bp), merged
  frac_cov99_tge2500        the same, counting only target contigs >= 2.5 kb (what MetaBAT2/CONCOCT bin)
  n_tcontigs99              target contigs carrying >= 1 kb of the MAG at >= 99 % identity (the genome's pieces)
  nga50_99 / lga90_99       contiguity of the same genome in the target: target contigs sorted by MAG bp they carry;
                            NGA50 = that bp at 50 % of MAG size (0 if < 50 % covered); LGA90 = pieces to 90 % (blank if < 90 %)
  n_tcontigs95 / nga50_95 / lga90_95   the same at >= 95 % identity: a (iii) contig that is a strain consensus
                            aligns to the arm at 95-99 %, and the 99 % versions count such stretches as gaps
  mag_n_contigs / mag_n50   the (iii) reference: the MAG's own contigs (all >= 2.5 kb)
  joins / splits            target contigs carrying >= 2 MAG contigs (each >= 1 kb) / MAG contigs >= 5 kb whose
                            >= 99 % sequence lies on >= 2 target contigs (each >= 1 kb)
  n_identical               MAG contigs reproduced end to end: one alignment covering the whole MAG contig and the whole
                            target contig (<= 5 bp overhang each side), >= 99.9 % identity -> the same (iii) contig
  cov_ratio_wmed            length-weighted median of SPAdes k-mer coverage, target contig / MAG contig, over >= 99 %
                            alignments >= 5 kb (2026 KBase vs (iii) was 2.9x lower for the same genome; ~1 means the
                            read set fed to SPAdes behaves like (iii)'s)
Summaries per sample x target x population: genome_check_summary.tsv. Re-runnable: PAFs and >= 1 kb targets are cached.

Usage: genome_check.py --root /scratch1/afreiburger/salterns_2026/tadpole_pilot --mags scripts/genome_check_mags.tsv
         --mag-dir eval/genome_check/mags --kbase-dir eval/genome_check/kbase2026 --out eval/genome_check
         --minimap2 /scratch1/afreiburger/metag-hybrid/env/miniforge/envs/assembly/bin/minimap2 --threads 8
"""
import argparse
import csv
import gzip
import os
import re
import statistics
import subprocess
import sys
from collections import defaultdict

MIN_TARGET = 1000
MIN_BLOCK = 500
PIECE_MIN = 1000


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_fasta(path):
    name, seq = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def cov_of(name):
    m = re.search(r"_cov_([0-9.]+)", name)
    return float(m.group(1)) if m else None


def n50(lens):
    s, acc = sorted(lens, reverse=True), 0
    half = sum(s) / 2.0
    for x in s:
        acc += x
        if acc >= half:
            return x
    return 0


def wmedian(pairs):
    pairs = sorted(p for p in pairs if p[0] is not None)
    if not pairs:
        return None
    half, acc = sum(w for _, w in pairs) / 2.0, 0
    for v, w in pairs:
        acc += w
        if acc >= half:
            return v
    return None


def union_len(iv):
    tot, cs, ce = 0, None, None
    for s, e in sorted(iv):
        if cs is None or s > ce:
            if cs is not None:
                tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    if cs is not None:
        tot += ce - cs
    return tot


def filter_target(src, dst):
    if os.path.exists(dst) and os.path.getmtime(dst) >= os.path.getmtime(src):
        return
    tmp = dst + ".tmp"
    with open(tmp, "w") as out:
        for n, s in read_fasta(src):
            if len(s) >= MIN_TARGET:
                out.write(">%s\n%s\n" % (n, s))
    os.replace(tmp, dst)


def run_minimap2(mm2, target, query, paf, threads):
    if os.path.exists(paf) and os.path.getmtime(paf) >= max(os.path.getmtime(target), os.path.getmtime(query)):
        return
    tmp = paf + ".tmp"
    with open(tmp, "w") as out:
        subprocess.run(["nice", "-n", "10", mm2, "-x", "asm10", "-c", "--secondary=no", "-t", str(threads),
                        target, query], stdout=out, stderr=subprocess.DEVNULL, check=True)
    os.replace(tmp, paf)


def read_paf(paf):
    hits = []
    with open(paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            tags = dict(t.split(":", 1) for t in f[12:] if t.count(":") >= 2)
            if tags.get("tp", "A:P")[2:] not in ("P",):
                continue
            q, ql, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
            t, tl, ts, te = f[5], int(f[6]), int(f[7]), int(f[8])
            nm, blk = int(f[9]), int(f[10])
            hits.append(dict(q=q, ql=ql, qs=qs, qe=qe, t=t, tl=tl, ts=ts, te=te, ident=nm / blk if blk else 0,
                             blk=blk, alen=qe - qs))
    return hits


def summarise_mag(mag_contigs, hits):
    """mag_contigs: {qname: (len, cov)}; hits: PAF records whose q is in mag_contigs."""
    size = sum(l for l, _ in mag_contigs.values())
    out = {"mag_bp": size, "mag_n_contigs": len(mag_contigs), "mag_n50": n50([l for l, _ in mag_contigs.values()]),
           "mag_kmer_cov_wmed": wmedian([(c, l) for l, c in mag_contigs.values()])}
    iv99, iv95, iv99g = defaultdict(list), defaultdict(list), defaultdict(list)
    per_t = defaultdict(int)                      # MAG bp aligned (>= 99 %) on each target contig
    per_t95 = defaultdict(int)                    # the same at >= 95 % (strain-consensus differences tolerated)
    per_t_q = defaultdict(lambda: defaultdict(int))
    per_q_t = defaultdict(lambda: defaultdict(int))
    tlen, ratios, identical, ident_bp = {}, [], set(), 0
    for h in hits:
        if h["blk"] < MIN_BLOCK:
            continue
        if h["ident"] >= 0.95:
            iv95[h["q"]].append((h["qs"], h["qe"]))
            per_t95[h["t"]] += h["alen"]
        if h["ident"] >= 0.99:
            iv99[h["q"]].append((h["qs"], h["qe"]))
            if h["tl"] >= 2500:
                iv99g[h["q"]].append((h["qs"], h["qe"]))
            per_t[h["t"]] += h["alen"]
            per_t_q[h["t"]][h["q"]] += h["alen"]
            per_q_t[h["q"]][h["t"]] += h["alen"]
            tlen[h["t"]] = h["tl"]
            if h["alen"] >= 5000:
                cq, ct = mag_contigs[h["q"]][1], cov_of(h["t"])
                if cq and ct:
                    ratios.append((ct / cq, h["alen"]))
        if (h["ident"] >= 0.999 and h["qs"] <= 5 and h["ql"] - h["qe"] <= 5 and h["ts"] <= 5
                and h["tl"] - h["te"] <= 5):
            identical.add(h["q"])
    for q in identical:
        ident_bp += mag_contigs[q][0]
    c99 = sum(union_len(v) for v in iv99.values())
    c95 = sum(union_len(v) for v in iv95.values())
    c99g = sum(union_len(v) for v in iv99g.values())
    def contiguity(pt):
        pcs = sorted((bp for bp in pt.values() if bp >= PIECE_MIN), reverse=True)
        acc, g50, l90 = 0, 0, ""
        for i, bp in enumerate(pcs, 1):
            acc += bp
            if not g50 and acc >= 0.5 * size:
                g50 = bp
            if l90 == "" and acc >= 0.9 * size:
                l90 = i
        return pcs, g50, l90
    pieces, nga50, lga90 = contiguity(per_t)
    pieces95, nga50_95, lga90_95 = contiguity(per_t95)
    big_t = [t for t, bp in per_t.items() if bp >= PIECE_MIN]
    joins = sum(1 for t in big_t if sum(1 for bp in per_t_q[t].values() if bp >= PIECE_MIN) >= 2)
    splits = sum(1 for q, d in per_q_t.items()
                 if mag_contigs[q][0] >= 5000 and sum(1 for bp in d.values() if bp >= PIECE_MIN) >= 2)
    cr = wmedian(ratios)
    out.update(frac_cov99=round(c99 / size, 4), frac_cov95=round(c95 / size, 4),
               frac_cov99_tge2500=round(c99g / size, 4), n_tcontigs99=len(pieces),
               tcontig_n50_99=n50([tlen[t] for t in big_t]) if big_t else 0,
               nga50_99=nga50, lga90_99=lga90, n_tcontigs95=len(pieces95), nga50_95=nga50_95,
               lga90_95=lga90_95, joins=joins, splits=splits,
               n_identical=len(identical), frac_bp_identical=round(ident_bp / size, 4),
               cov_ratio_wmed=round(cr, 3) if cr is not None else "")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True)
    ap.add_argument("--mags", required=True)
    ap.add_argument("--mag-dir", required=True)
    ap.add_argument("--kbase-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--minimap2", required=True)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()

    mags = list(csv.DictReader(open(a.mags), delimiter="\t"))
    os.makedirs(os.path.join(a.out, "asm"), exist_ok=True)
    os.makedirs(os.path.join(a.out, "paf"), exist_ok=True)
    samples = sorted({s for m in mags for s in m["map_to"].split(",")})
    rows = []
    for smp in samples:
        mine = [m for m in mags if smp in m["map_to"].split(",")]
        # query: MAG contigs renamed <i>|<contig> (lost30 MAGs from other samples share NODE names)
        query = os.path.join(a.out, "asm", "%s.query.fna" % smp)
        contigs = {}
        tmp = query + ".tmp"
        with open(tmp, "w") as out:
            for i, m in enumerate(mine):
                d = {}
                for n, s in read_fasta(os.path.join(a.mag_dir, m["file"])):
                    qn = "%d|%s" % (i, n)
                    d[qn] = (len(s), cov_of(n))
                    out.write(">%s\n%s\n" % (qn, s))
                contigs[i] = d
        if os.path.exists(query) and open(query).read() == open(tmp).read():
            os.remove(tmp)                         # keep the old mtime so cached PAFs stay valid
        else:
            os.replace(tmp, query)
        targets = []
        for arm in ("armA", "armB", "armC"):
            src = os.path.join(a.root, smp, arm, "spades", "scaffolds.fasta")
            if os.path.exists(src) and os.path.exists(os.path.join(a.root, smp, ".done", "C1_spades_default"
                                                                   if arm == "armC" else
                                                                   {"armA": "A4_metaspades",
                                                                    "armB": "B3_metaspades"}[arm])):
                targets.append((arm, src))
        kb = os.path.join(a.kbase_dir, "%s_kbase2026.fna" % smp)
        if os.path.exists(kb):
            targets.append(("kbase2026", kb))
        for name, src in targets:
            tgt = os.path.join(a.out, "asm", "%s.%s.ge1000.fna" % (smp, name))
            filter_target(src, tgt)
            paf = os.path.join(a.out, "paf", "%s.%s.paf" % (smp, name))
            run_minimap2(a.minimap2, tgt, query, paf, a.threads)
            byq = defaultdict(list)
            for h in read_paf(paf):
                byq[int(h["q"].split("|", 1)[0])].append(h)
            for i, m in enumerate(mine):
                r = dict(sample=smp, target=name, mag=m["mag"], population=m["population"],
                         source_sample=m["source_sample"], genus=m.get("genus", ""),
                         completeness=m.get("completeness", ""), contamination=m.get("contamination", ""))
                r.update(summarise_mag(contigs[i], byq.get(i, [])))
                rows.append(r)
            print("%s %s: %d MAGs" % (smp, name, len(mine)), flush=True)

    cols = list(rows[0])
    with open(os.path.join(a.out, "genome_check_mags_by_target.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # summary per sample x target x population (own-sample MAGs: derep310, nonderep140, both; lost30)
    groups = defaultdict(list)
    for r in rows:
        groups[(r["sample"], r["target"], r["population"])].append(r)
        if r["population"] in ("derep310", "nonderep140"):
            groups[(r["sample"], r["target"], "own_all")].append(r)
    summ = []
    for (smp, tgt, pop), rs in sorted(groups.items()):
        bp = sum(r["mag_bp"] for r in rs)
        summ.append(dict(
            sample=smp, target=tgt, population=pop, n_mags=len(rs), mag_Mbp=round(bp / 1e6, 2),
            pct_cov99=round(100 * sum(r["frac_cov99"] * r["mag_bp"] for r in rs) / bp, 1),
            pct_cov95=round(100 * sum(r["frac_cov95"] * r["mag_bp"] for r in rs) / bp, 1),
            pct_cov99_tge2500=round(100 * sum(r["frac_cov99_tge2500"] * r["mag_bp"] for r in rs) / bp, 1),
            mag_contigs=sum(r["mag_n_contigs"] for r in rs), target_pieces99=sum(r["n_tcontigs99"] for r in rs),
            median_mag_n50=statistics.median(r["mag_n50"] for r in rs),
            median_nga50_99=statistics.median(r["nga50_99"] for r in rs),
            n_mags_nga50_ge_mag_n50=sum(1 for r in rs if r["nga50_99"] >= r["mag_n50"]),
            n_mags_pieces_le_mag_contigs=sum(1 for r in rs if 0 < r["n_tcontigs99"] <= r["mag_n_contigs"]),
            target_pieces95=sum(r["n_tcontigs95"] for r in rs),
            median_nga50_95=statistics.median(r["nga50_95"] for r in rs),
            n_mags_nga50_95_ge_mag_n50=sum(1 for r in rs if r["nga50_95"] >= r["mag_n50"]),
            joins=sum(r["joins"] for r in rs), splits=sum(r["splits"] for r in rs),
            pct_contigs_identical=round(100 * sum(r["n_identical"] for r in rs) / sum(r["mag_n_contigs"] for r in rs), 1),
            pct_bp_identical=round(100 * sum(r["frac_bp_identical"] * r["mag_bp"] for r in rs) / bp, 1),
            median_cov_ratio=(statistics.median(r["cov_ratio_wmed"] for r in rs if r["cov_ratio_wmed"] != "")
                              if any(r["cov_ratio_wmed"] != "" for r in rs) else "")))
    with open(os.path.join(a.out, "genome_check_summary.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summ[0]), delimiter="\t")
        w.writeheader()
        w.writerows(summ)
    for s in summ:
        print("\t".join(str(v) for v in s.values()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
