#!/usr/bin/env python3
"""KOfam calls (DRAM's KO layer) + marker / guild calls for one finished hmmsearch chunk, or combine all chunks.

Rule (DRAM 1.5 `sig_scores` + `generic_hmmscan_formater`, sources in annotation_choice/sources/dram1_v1.5.0/):
  * a (gene, KO) hit is significant when score > KO threshold, with score = full-sequence score if ko_list
    score_type == 'full', or the domain score if 'domain' (DRAM tests every domain row; here: the best domain
    score of the hit, which is equivalent);
  * KOs whose ko_list threshold is '-' are searched but never called (DRAM 1.5 skips them);
  * best KO per gene = significant hit with the lowest full-sequence E-value (DRAM); ties (e.g. E = 0) broken by
    the higher full score, then KO id;
  * 'all passing KOs' = every significant KO of the gene (cf. DRAM2's kofam_ids; DRAM2 uses >=, this uses >).
hmmsearch was run with -Z <total proteins> so E-values do not depend on chunking, and -T/--domT at the lowest KO
threshold so that every hit that could pass is reported. Each chunk searched one representative per distinct
protein sequence; hits are copied to identical proteins (chunks/chunk_XX.members.tsv).

Guild rules (revised, annotation_choice/README.md section 6), counted on best-KO calls (primary) and on
'any passing KO' (column suffix _anypass):
  methanogen        = mcrA (K00399)
  mepn_consumer     = phnJ (K06163); + counts of phnG/H/I/K/L/M, n distinct of phnGHIJKLM, n contigs carrying them
  betaine_reducer   = grdH (K21579) or grdI (K21578); grdA/C/D reported as support
  methylamine_user  = mttB (K14083) / mtbB (K16178) / mtmB (K16176) only in an mcrA genome
  tma_producer      = cutC+cutD (K20038+K20037) or cntA+cntB (K22443+K22444); one-sided calls -> tma_candidate
Flags marker genes on MaxBin2 contigs < 2,500 bp (DRAM's default minimum) and whether a guild call depends on them.

stdlib only.
usage: postprocess.py chunk <workdir> <chunk_id>      (writes results/chunks/chunk_XX.* and results/per_mag/<MAG>.kofam.tsv.gz)
       postprocess.py combine <workdir>                (concatenates finished chunks into results/all_*.tsv[.gz])
"""
import gzip
import os
import sys
from collections import defaultdict

SHORT = 2500
GUILD_MARKERS = {   # symbol -> KO, for the guild rules
    "mcrA": "K00399", "phnG": "K06166", "phnH": "K06165", "phnI": "K06164", "phnJ": "K06163", "phnK": "K05781",
    "phnL": "K05780", "phnM": "K06162", "grdH": "K21579", "grdI": "K21578", "grdA": "K10670", "grdC": "K21577",
    "grdD": "K21576", "mttB": "K14083", "mtbB": "K16178", "mtmB": "K16176", "mtgB": "K25220", "cutC": "K20038",
    "cutD": "K20037", "cntA": "K22443", "cntB": "K22444",
}
PHN7 = ["phnG", "phnH", "phnI", "phnJ", "phnK", "phnL", "phnM"]


def opener(p):
    if os.path.exists(p):
        return open(p)
    return gzip.open(p + ".gz", "rt")


def read_tsv(p):
    with opener(p) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(hdr, line.rstrip("\n").split("\t")))


def load_meta(W):
    info = {r["ko"]: r for r in read_tsv(os.path.join(W, "db/hmm_info.tsv"))}
    markers = {r["ko"]: r["gene"] for r in read_tsv(os.path.join(W, "meta/marker_ko_table.tsv"))}
    for sym, ko in GUILD_MARKERS.items():
        assert ko in markers, (sym, ko)
    return info, markers


def guilds(cnt, contigs_of):
    """cnt: symbol -> n genes; contigs_of: symbol -> set(contigs)."""
    g = {}
    g["methanogen"] = cnt["mcrA"] > 0
    g["mepn_consumer"] = cnt["phnJ"] > 0
    g["phnGHIJKLM_n_distinct"] = sum(cnt[s] > 0 for s in PHN7)
    g["phnGHIJKLM_n_contigs"] = len(set().union(*(contigs_of[s] for s in PHN7)))
    g["betaine_reducer"] = cnt["grdH"] > 0 or cnt["grdI"] > 0
    g["grdACD_n_distinct"] = sum(cnt[s] > 0 for s in ("grdA", "grdC", "grdD"))
    mt = cnt["mttB"] + cnt["mtbB"] + cnt["mtmB"]
    g["methylamine_user"] = mt > 0 and cnt["mcrA"] > 0
    g["mttB_family_without_mcrA"] = (mt + cnt["mtgB"]) > 0 and cnt["mcrA"] == 0
    cut, cnt_ = cnt["cutC"] > 0 and cnt["cutD"] > 0, cnt["cntA"] > 0 and cnt["cntB"] > 0
    g["tma_producer"] = cut or cnt_
    g["tma_candidate_one_sided"] = (not (cut or cnt_)) and (cnt["cutC"] + cnt["cntA"]) > 0
    return g


GUILD_BOOL = ["methanogen", "mepn_consumer", "betaine_reducer", "methylamine_user", "tma_producer"]


def do_chunk(W, cid):
    info, markers = load_meta(W)
    R = os.path.join(W, "results")
    for d in ("chunks", "per_mag"):
        os.makedirs(os.path.join(R, d), exist_ok=True)
    mags = [r["mag"] for r in read_tsv(os.path.join(W, "chunks/chunk_assignment.tsv")) if r["chunk"] == cid]
    genes = {}
    for m in mags:
        for r in read_tsv(os.path.join(W, "genes", m + ".genes.tsv")):
            genes[r["gene"]] = r
    pref = os.path.join(W, "hmm", f"chunk_{cid}")

    # ---- pass 1, tblout: one line per reported (gene, KO) with the full-sequence score and the best single
    #      domain score (that domain's score even if it was below --domT) -> significance per DRAM's rule
    passing = defaultdict(list)          # gene -> [ko, score, thr, stype, full_e, full_s, dom_s, cov]
    marker_rows, nothr_rows = [], []
    need_cov = {}                        # (gene, ko) -> [hmm_from, hmm_to, n_domains], filled from domtblout
    n_tbl = 0
    with opener(pref + ".tblout") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split(None, 18)
            n_tbl += 1
            gene, ko = f[0], f[2]
            full_e, full_s, dom_s = float(f[4]), float(f[5]), float(f[8])
            inf = info[ko]
            thr, stype = inf["threshold"], inf["score_type"]
            if stype in ("full", "domain"):
                t = float(thr)
                score = full_s if stype == "full" else dom_s
                sig = score > t
            else:
                t, score, sig = None, None, False
            if sig:
                passing[gene].append([ko, score, t, stype, full_e, full_s, dom_s, None])
                need_cov[(gene, ko)] = None
            if ko in markers:
                marker_rows.append([gene, ko, markers[ko], stype, thr, score, full_s, dom_s, full_e, None, None, sig])
                need_cov[(gene, ko)] = None
            elif stype not in ("full", "domain") and full_e <= 1e-5:
                nothr_rows.append([gene, ko, full_s, dom_s, full_e, None, None])
                need_cov[(gene, ko)] = None
    # ---- pass 2, domtblout: HMM coverage (min hmm_from .. max hmm_to over reported domains) for the kept pairs
    n_rows, dom_pairs = 0, set()
    with opener(pref + ".domtblout") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split(None, 17)
            n_rows += 1
            key = (f[0], f[3])
            dom_pairs.add(hash(key))
            if key in need_cov:
                hf, ht = int(f[15]), int(f[16])
                a = need_cov[key]
                if a is None:
                    need_cov[key] = [hf, ht, 1]
                else:
                    a[0], a[1], a[2] = min(a[0], hf), max(a[1], ht), a[2] + 1
    n_pairs = len(dom_pairs)
    del dom_pairs

    def cov_of(gene, ko):
        a = need_cov.get((gene, ko))
        return ((a[1] - a[0] + 1) / int(info[ko]["leng"]), a[2]) if a else (0.0, 0)

    pass_nocov = 0                       # significant pairs whose domains were all below --domT (coverage 0)
    for gene, hits in passing.items():
        for h in hits:
            h[7] = cov_of(gene, h[0])[0]
            pass_nocov += need_cov.get((gene, h[0])) is None
    for r in marker_rows:
        r[9], r[10] = cov_of(r[0], r[1])
    for r in nothr_rows:
        r[5], r[6] = cov_of(r[0], r[1])
    # ---- copy each representative's hits to the genes with an identical protein (chunk-internal de-duplication)
    members = defaultdict(list)
    n_members = 0
    for r in read_tsv(os.path.join(W, "chunks", f"chunk_{cid}.members.tsv")):
        members[r["representative_gene"]].append(r["member_gene"])
        n_members += 1
    n_searched = len(genes) - n_members
    for rep in [g for g in passing if g in members]:
        for mem in members[rep]:
            passing[mem] = [list(h) for h in passing[rep]]
    marker_rows += [[mem] + r[1:] for r in marker_rows if r[0] in members for mem in members[r[0]]]
    nothr_rows += [[mem] + r[1:] for r in nothr_rows if r[0] in members for mem in members[r[0]]]
    unknown = [g for g in passing if g not in genes]
    assert not unknown, unknown[:5]

    best = {}
    for gene, hits in passing.items():
        best[gene] = sorted(hits, key=lambda h: (h[4], -h[5], h[0]))[0]

    # ---- per-gene table for genes with >= 1 passing KO
    with gzip.open(os.path.join(R, "chunks", f"chunk_{cid}.genes_ko.tsv.gz"), "wt") as fo:
        fo.write("gene\tmag\tbinner\tcontig\tcontig_len\tbest_ko\tbest_score\tbest_threshold\tbest_score_type\t"
                 "best_full_evalue\tbest_hmm_cov\tn_passing\tall_passing\n")
        for gene in sorted(passing):
            g, b = genes[gene], best[gene]
            allp = ";".join(f"{h[0]}:{h[1]:.1f}" for h in sorted(passing[gene], key=lambda h: (h[4], -h[5], h[0])))
            fo.write(f"{gene}\t{g['mag']}\t{g['binner']}\t{g['contig']}\t{g['contig_len']}\t{b[0]}\t{b[1]:.1f}\t"
                     f"{b[2]}\t{b[3]}\t{b[4]:.3g}\t{b[7]:.2f}\t{len(passing[gene])}\t{allp}\n")

    # ---- per-MAG gene-level annotation (all genes, DRAM-like)
    by_mag = defaultdict(list)
    for gene, g in genes.items():
        by_mag[g["mag"]].append(gene)
    for m in mags:
        with gzip.open(os.path.join(R, "per_mag", f"{m}.kofam.tsv.gz"), "wt") as fo:
            fo.write("gene\tcontig\tcontig_len\tstart\tend\tstrand\tpartial\taa_len\tkofam_id\tkofam_score\t"
                     "kofam_threshold\tkofam_score_type\tkofam_full_evalue\tkofam_hmm_cov\tkofam_definition\t"
                     "kofam_all_passing\n")
            for gene in by_mag[m]:
                g = genes[gene]
                row = [gene, g["contig"], g["contig_len"], g["start"], g["end"], g["strand"], g["partial"], g["aa_len"]]
                if gene in best:
                    b = best[gene]
                    row += [b[0], f"{b[1]:.1f}", b[2], b[3], f"{b[4]:.3g}", f"{b[7]:.2f}", info[b[0]]["definition"],
                            ";".join(h[0] for h in sorted(passing[gene], key=lambda h: (h[4], -h[5], h[0])))]
                else:
                    row += [""] * 8
                fo.write("\t".join(map(str, row)) + "\n")

    # ---- MAG x KO counts (best KO, any passing KO)
    with open(os.path.join(R, "chunks", f"chunk_{cid}.mag_ko_counts.tsv"), "w") as fo:
        fo.write("mag\tko\tn_best\tn_anypass\n")
        for m in mags:
            nb, na = defaultdict(int), defaultdict(int)
            for gene in by_mag[m]:
                if gene in best:
                    nb[best[gene][0]] += 1
                    for h in passing[gene]:
                        na[h[0]] += 1
            for ko in sorted(na):
                fo.write(f"{m}\t{ko}\t{nb.get(ko, 0)}\t{na[ko]}\n")

    # ---- marker hits (any reported score) with gene context
    sym_of = {v: k for k, v in GUILD_MARKERS.items()}
    with gzip.open(os.path.join(R, "chunks", f"chunk_{cid}.marker_hits.tsv.gz"), "wt") as fo:
        fo.write("gene\tmag\tbinner\tcontig\tcontig_len\tshort_maxbin2_contig\taa_len\tpartial\tko\tmarker\t"
                 "score_type\tthreshold\tscore\tfull_score\tbest_domain_score\tfull_evalue\thmm_cov\tn_domains\t"
                 "passes\tis_best_ko\tgene_best_ko\tgene_best_marker\n")
        for r in sorted(marker_rows, key=lambda r: (genes[r[0]]["mag"], r[0], r[8])):
            g = genes[r[0]]
            b = best.get(r[0])
            short = g["binner"] == "maxbin2" and int(g["contig_len"]) < SHORT
            fo.write("\t".join(map(str, [
                r[0], g["mag"], g["binner"], g["contig"], g["contig_len"], int(short), g["aa_len"], g["partial"],
                r[1], r[2], r[3], r[4], "" if r[5] is None else f"{r[5]:.1f}", f"{r[6]:.1f}", f"{r[7]:.1f}",
                f"{r[8]:.3g}", f"{r[9]:.2f}", r[10], int(r[11]), int(b is not None and b[0] == r[1]),
                b[0] if b else "", markers.get(b[0], "") if b else ""])) + "\n")
    with gzip.open(os.path.join(R, "chunks", f"chunk_{cid}.no_threshold_hits.tsv.gz"), "wt") as fo:
        fo.write("gene\tmag\tko\tfull_score\tbest_domain_score\tfull_evalue\thmm_cov\tn_domains\tgene_best_ko\n")
        for r in nothr_rows:
            fo.write("\t".join(map(str, [r[0], genes[r[0]]["mag"], r[1], f"{r[2]:.1f}", f"{r[3]:.1f}", f"{r[4]:.3g}",
                                         f"{r[5]:.2f}", r[6], best[r[0]][0] if r[0] in best else ""])) + "\n")

    # ---- per-MAG marker counts + guild calls
    mk_syms = list(dict.fromkeys(markers.values()))
    cols_counts = [f"n_{s}" for s in mk_syms]
    rows = []
    for m in mags:
        binner = m.split(".")[0]
        cb = defaultdict(int)                                       # marker-table symbol -> n genes (best KO)
        cg, ca, cg_long = defaultdict(int), defaultdict(int), defaultdict(int)   # guild symbol -> n genes
        con_b, con_a, con_bl = defaultdict(set), defaultdict(set), defaultdict(set)
        short_list = []
        for gene in by_mag[m]:
            if gene not in best:
                continue
            g = genes[gene]
            short = binner == "maxbin2" and int(g["contig_len"]) < SHORT
            ko = best[gene][0]
            if ko in markers:
                cb[markers[ko]] += 1
                if ko in sym_of:
                    s = sym_of[ko]
                    cg[s] += 1
                    con_b[s].add(g["contig"])
                    if not short:
                        cg_long[s] += 1
                        con_bl[s].add(g["contig"])
                if short:
                    short_list.append(f"{markers[ko]}@{gene.split('|', 1)[1]}({g['contig_len']})")
            for h in passing[gene]:
                if h[0] in sym_of:
                    ca[sym_of[h[0]]] += 1
                    con_a[sym_of[h[0]]].add(g["contig"])
        gsym = cg
        gb = guilds(cg, con_b)
        ga = guilds(ca, con_a)
        gl = guilds(cg_long, con_bl)
        row = {"mag": m, "binner": binner, "unit": m.split("__", 1)[1], "chunk": cid,
               "n_genes": len(by_mag[m]), "n_genes_with_ko": sum(1 for x in by_mag[m] if x in best)}
        row.update({k: int(v) if isinstance(v, bool) else v for k, v in gb.items()})
        row.update({f"{k}_anypass": int(v) if isinstance(v, bool) else v for k, v in ga.items() if k in GUILD_BOOL})
        row["n_marker_genes_on_short_maxbin2_contigs"] = len(short_list)
        row["short_maxbin2_marker_genes"] = ";".join(short_list)
        row["guilds_lost_without_short_maxbin2_contigs"] = ";".join(k for k in GUILD_BOOL if gb[k] and not gl[k])
        row["phn_accessory_counts"] = ";".join(f"{s}:{gsym[s]}" for s in PHN7 if s != "phnJ")
        for s, c in zip(mk_syms, cols_counts):
            row[c] = cb[s]
        rows.append(row)
    with open(os.path.join(R, "chunks", f"chunk_{cid}.mag_markers_guilds.tsv"), "w") as fo:
        keys = list(rows[0].keys())
        fo.write("\t".join(keys) + "\n")
        for r in rows:
            fo.write("\t".join(str(r[k]) for k in keys) + "\n")

    with open(os.path.join(R, "chunks", f"chunk_{cid}.qc.txt"), "w") as fo:
        fo.write(f"chunk\t{cid}\nmags\t{len(mags)}\ngenes\t{len(genes)}\tsearched_representatives\t{n_searched}\t"
                 f"copied_to_identical_members\t{n_members}\ndomtblout_rows\t{n_rows}\n"
                 f"gene_ko_pairs_domtblout\t{n_pairs}\ngene_ko_pairs_tblout\t{n_tbl}\n"
                 f"tblout_pairs_without_domtblout_row\t{n_tbl - n_pairs}\t"
                 f"passing_pairs_without_domtblout_row\t{pass_nocov}\n"
                 f"genes_with_ko\t{len(best)}\t"
                 f"fraction\t{len(best) / max(1, len(genes)):.4f}\n"
                 f"marker_hit_rows\t{len(marker_rows)}\nno_threshold_hits_e1e-5\t{len(nothr_rows)}\n")
    with open(os.path.join(R, "chunks", f"chunk_{cid}.postprocessed"), "w") as fo:
        fo.write("ok\n")
    print(open(os.path.join(R, "chunks", f"chunk_{cid}.qc.txt")).read())


def do_combine(W):
    R = os.path.join(W, "results")
    done = sorted(f.split("_")[1].split(".")[0] for f in os.listdir(os.path.join(R, "chunks"))
                  if f.endswith(".postprocessed"))
    for name, gz in (("mag_markers_guilds.tsv", False), ("mag_ko_counts.tsv", False), ("marker_hits.tsv.gz", True),
                     ("genes_ko.tsv.gz", True), ("no_threshold_hits.tsv.gz", True)):
        out = os.path.join(R, "all_" + name)
        op = gzip.open if gz else open
        with op(out, "wt") as fo:
            for i, c in enumerate(done):
                with op(os.path.join(R, "chunks", f"chunk_{c}.{name}"), "rt") as fh:
                    hdr = fh.readline()
                    if i == 0:
                        fo.write(hdr)
                    for line in fh:
                        fo.write(line)
    n_mags = sum(1 for r in read_tsv(os.path.join(W, "chunks/chunk_assignment.tsv")) if r["chunk"] in done)
    with open(os.path.join(R, "combine_status.txt"), "w") as fo:
        fo.write(f"chunks_combined\t{len(done)}\t{','.join(done)}\nmags_covered\t{n_mags}\n")
    print(open(os.path.join(R, "combine_status.txt")).read())


if __name__ == "__main__":
    if sys.argv[1] == "chunk":
        do_chunk(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "combine":
        do_combine(sys.argv[2])
    else:
        sys.exit(__doc__)
