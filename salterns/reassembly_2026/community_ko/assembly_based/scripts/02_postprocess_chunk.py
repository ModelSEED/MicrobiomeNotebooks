#!/usr/bin/env python3
"""Step 2b: KOfam calls for one finished hmmsearch chunk (representative proteins only; identical proteins get the
same call in step 3).

Call rule = kofam_rollout/scripts/postprocess.py (DRAM 1.5 `sig_scores`):
  * a (gene, KO) hit passes when score > KO threshold (score = full-sequence score if ko_list score_type == 'full',
    best single-domain score (tblout column 9) if 'domain');
  * KOs with threshold '-' are never called;
  * best KO = passing hit with the lowest full-sequence E-value; ties -> higher full score -> KO id.
Marker rescue (for the genome-equivalent denominator only; read_markers/scripts/04_build_db.py 'marker_relaxed_0.5thr'):
  a gene with NO passing KO whose top-scoring profile (highest full score) is one of the 10 marker KOs and whose
  score (by that KO's score_type) is >= 0.5 x its threshold counts for that marker. Two readings are recorded:
    rescue_all   : top-scoring profile over ALL reported KOfam hits (used here; full DB searched),
    rescue_panel : top-scoring profile over the 62 KOs of the v1 read_markers panel (v1 searched only those).
Also recorded for every gene with any hit: the top-scoring KO and its full score / E-value (used for the
no-threshold KOs K18049 mpnS / K12905 hepD, which the rule never calls).
stdlib only.  usage: 02_postprocess_chunk.py <W> <chunk_id>
"""
import gzip
import os
import sys

W, cid = sys.argv[1], sys.argv[2]
info = {}
with open(os.path.join(W, "db/hmm_info.tsv")) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        r = dict(zip(hdr, line.rstrip("\n").split("\t")))
        st = r["score_type"]
        info[r["ko"]] = (float(r["threshold"]), st) if st in ("full", "domain") else (None, st)
MARKERS = {"K02886", "K02867", "K02982", "K02996", "K02948", "K02950", "K02926", "K02930", "K02878", "K02866"}
V1_PANEL = {l.strip().split("\t")[0] for l in open(os.path.join(W, "meta/v1_panel_kos.txt")) if l.strip()}
assert MARKERS <= V1_PANEL and len(V1_PANEL) == 62, len(V1_PANEL)

pref = os.path.join(W, "hmm", f"chunk_{cid}.tblout")
src = pref if os.path.exists(pref) else pref + ".gz"
op = open if src == pref else (lambda p: gzip.open(p, "rt"))

passing = {}       # gene -> list of (E, -full, ko, score, thr, stype)
top = {}           # gene -> (full, ko, dom, E)   highest full score over all hits (tie -> lower KO id)
top_panel = {}     # gene -> same over the v1 panel KOs
n_lines, last = 0, ""
with op(src) as fh:
    for line in fh:
        if line.startswith("#"):
            last = line.rstrip("\n")
            continue
        f = line.split(None, 18)
        n_lines += 1
        gene, ko = f[0], f[2]
        full_e, full_s, dom_s = float(f[4]), float(f[5]), float(f[8])
        t, st = info[ko]
        if t is not None:
            score = full_s if st == "full" else dom_s
            if score > t:
                passing.setdefault(gene, []).append((full_e, -full_s, ko, score, t, st))
        cand = (full_s, ko, dom_s, full_e)
        cur = top.get(gene)
        if cur is None or full_s > cur[0] or (full_s == cur[0] and ko < cur[1]):
            top[gene] = cand
        if ko in V1_PANEL:
            cur = top_panel.get(gene)
            if cur is None or full_s > cur[0] or (full_s == cur[0] and ko < cur[1]):
                top_panel[gene] = cand
assert last == "# [ok]", f"tblout not complete: {src} last={last!r}"


def rescue(t):
    if t is None or t[1] not in MARKERS:
        return ""
    thr, st = info[t[1]]
    sc = t[0] if st == "full" else t[2]
    return t[1] if sc >= 0.5 * thr else ""


n_reps = sum(1 for l in open(os.path.join(W, "chunks", f"chunk_{cid}.faa")) if l.startswith(">"))
out = os.path.join(W, "results/chunks", f"chunk_{cid}.calls.tsv.gz")
n_pass = n_res_all = n_res_panel = 0
with gzip.open(out + ".part", "wt", compresslevel=4) as fo:
    fo.write("rep\tbest_ko\tbest_score\tbest_threshold\tbest_score_type\tbest_full_evalue\tn_passing\tall_passing\t"
             "top_ko\ttop_full_score\ttop_full_evalue\trescue_all\trescue_panel\n")
    for gene in sorted(top):
        hits = sorted(passing.get(gene, []))
        tp = top[gene]
        ra = rp = ""
        if hits:
            b = hits[0]
            n_pass += 1
            bcols = [b[2], f"{b[3]:.1f}", f"{b[4]:.2f}", b[5], f"{b[0]:.3g}", str(len(hits)),
                     ";".join(f"{h[2]}:{h[3]:.1f}" for h in hits)]
        else:
            ra, rp = rescue(tp), rescue(top_panel.get(gene))
            n_res_all += bool(ra)
            n_res_panel += bool(rp)
            bcols = ["", "", "", "", "", "0", ""]
        fo.write("\t".join([gene] + bcols + [tp[1], f"{tp[0]:.1f}", f"{tp[3]:.3g}", ra, rp]) + "\n")
os.replace(out + ".part", out)
with open(os.path.join(W, "results/chunks", f"chunk_{cid}.qc.txt"), "w") as fo:
    fo.write(f"chunk\t{cid}\nreps\t{n_reps}\ntblout_rows\t{n_lines}\nreps_with_any_hit\t{len(top)}\n"
             f"reps_with_passing_ko\t{n_pass}\nmarker_rescue_all\t{n_res_all}\nmarker_rescue_panel\t{n_res_panel}\n")
open(os.path.join(W, "results/chunks", f"chunk_{cid}.postprocessed"), "w").write("ok\n")
