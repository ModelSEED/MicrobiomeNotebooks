#!/usr/bin/env python3
"""Equivalence check of the two-pass design against the full search, on the 4 chunks that finished all 27,981
profiles before the full run was stopped (full_run_stopped/: chunks 0027, 0046, 0049, 0063).

(a) pass-1 tblout rows (148 profiles) vs the full run's tblout rows for the same (gene, KO): all 18 fields;
(b) per representative of the chunk: panel/marker label (best KO if it is one of the 148, else none), marker rescue
    (both readings), and for every candidate the complete call row (best KO, score, threshold, E-value, n passing,
    all passing KOs with scores, top-scoring KO, its score and E-value) from pass 2 vs the full run;
(c) completeness of the candidate set: every representative that passes one of the 148 KOs or is rescued in the full
    run is a candidate.
Writes results/two_pass_equivalence_check.txt and results/two_pass_equivalence_diffs.tsv.
usage: 09_two_pass_check.py <W>
"""
import gzip
import os
import sys

W = sys.argv[1]
os.chdir(W)
P148 = {l.strip() for l in open("meta/pass1_kos.txt") if l.strip()}
FULL = sorted(f.split("_")[1].split(".")[0] for f in os.listdir("full_run_stopped/results_chunks")
              if f.endswith(".postprocessed"))


def read_calls(path):
    out = {}
    with gzip.open(path, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = line.rstrip("\n").split("\t")
            out[r[0]] = dict(zip(h, r))
    return out


cand = set()
with gzip.open("pass2/candidates.tsv.gz", "rt") as fh:
    fh.readline()
    for line in fh:
        cand.add(line.split("\t", 1)[0])
p2 = {}
for f in os.listdir("pass2/results/chunks"):
    if f.endswith(".calls.tsv.gz"):
        p2.update(read_calls(f"pass2/results/chunks/{f}"))
COLS = ["best_ko", "best_score", "best_threshold", "best_score_type", "best_full_evalue", "n_passing", "all_passing",
        "top_ko", "top_full_score", "top_full_evalue", "rescue_all", "rescue_panel"]
lines, diffs = [], []
tot = dict(reps=0, cand=0, label_same=0, rescue_same=0, cand_rows_same=0, missed=0, rows_full148=0, rows_p1=0,
           rows_same=0, panel_pass_nonpanel_best=0, panel_labelled=0)
for c in FULL:
    reps = [l[1:].split()[0] for l in open(f"chunks/chunk_{c}.faa") if l[0] == ">"]
    full = read_calls(f"full_run_stopped/results_chunks/chunk_{c}.calls.tsv.gz")
    # (a) tblout rows
    fr, pr = set(), set()
    with gzip.open(f"full_run_stopped/hmm/chunk_{c}.tblout.gz", "rt") as fh:
        for line in fh:
            if line[0] != "#":
                x = tuple(line.split(None, 18)[:18])
                if x[2] in P148:
                    fr.add(x)
    for line in open(f"pass1/hmm/chunk_{c}.tblout"):
        if line[0] != "#":
            pr.add(tuple(line.split(None, 18)[:18]))
    tot["rows_full148"] += len(fr)
    tot["rows_p1"] += len(pr)
    tot["rows_same"] += len(fr & pr)
    # (b), (c)
    for r in reps:
        tot["reps"] += 1
        f = full.get(r, {})
        q = p2.get(r, {}) if r in cand else {}
        isc = r in cand
        tot["cand"] += isc
        fl = f.get("best_ko", "") if f.get("best_ko", "") in P148 else ""
        ql = q.get("best_ko", "") if q.get("best_ko", "") in P148 else ""
        tot["panel_labelled"] += fl != ""
        same_l = fl == ql
        same_r = (f.get("rescue_all", ""), f.get("rescue_panel", "")) == (q.get("rescue_all", ""), q.get("rescue_panel", ""))
        tot["label_same"] += same_l
        tot["rescue_same"] += same_r
        passes148 = any(x.split(":")[0] in P148 for x in f.get("all_passing", "").split(";") if x)
        if (passes148 or f.get("rescue_all") or f.get("rescue_panel")) and not isc:
            tot["missed"] += 1
        if passes148 and f.get("best_ko", "") not in P148:
            tot["panel_pass_nonpanel_best"] += 1
        if isc:
            same_row = all(f.get(k, "") == q.get(k, "") for k in COLS)
            tot["cand_rows_same"] += same_row
        else:
            same_row = True
        if not (same_l and same_r and same_row):
            diffs.append((c, r, isc, fl, ql, f.get("best_ko", ""), q.get("best_ko", ""), f.get("best_score", ""),
                          q.get("best_score", ""), f.get("all_passing", ""), q.get("all_passing", "")))
    lines.append(f"chunk_{c}\treps\t{len(reps)}\tcandidates\t{sum(1 for r in reps if r in cand)}\t"
                 f"tblout_rows_148_full\t{len(fr)}\tpass1\t{len(pr)}\tidentical\t{len(fr & pr)}\n")
with open("results/two_pass_equivalence_diffs.tsv", "w") as fo:
    fo.write("chunk\trep\tcandidate\tlabel_full\tlabel_two_pass\tbest_full\tbest_two_pass\tscore_full\tscore_two_pass\t"
             "all_passing_full\tall_passing_two_pass\n")
    for d in diffs:
        fo.write("\t".join(map(str, d)) + "\n")
with open("results/two_pass_equivalence_check.txt", "w") as fo:
    fo.write(f"full-run chunks checked\t{','.join(FULL)}\n")
    fo.writelines(lines)
    fo.write(f"tblout rows of the 148 profiles: full run {tot['rows_full148']}, pass 1 {tot['rows_p1']}, "
             f"identical (all 18 fields) {tot['rows_same']}\n"
             f"representatives\t{tot['reps']}\tcandidates\t{tot['cand']}\n"
             f"panel/marker label identical\t{tot['label_same']}\tof\t{tot['reps']}\t"
             f"(reps with a panel/marker best KO in the full run: {tot['panel_labelled']})\n"
             f"marker rescue (both readings) identical\t{tot['rescue_same']}\tof\t{tot['reps']}\n"
             f"candidate call rows identical (all {len(COLS)} columns)\t{tot['cand_rows_same']}\tof\t{tot['cand']}\n"
             f"reps passing a 148 KO or rescued in the full run but not candidates\t{tot['missed']}\n"
             f"reps passing a 148 KO whose best KO is another KO (counted for that KO, as in the full run)\t"
             f"{tot['panel_pass_nonpanel_best']}\n"
             f"differing reps\t{len(diffs)}\n")
print(open("results/two_pass_equivalence_check.txt").read())
