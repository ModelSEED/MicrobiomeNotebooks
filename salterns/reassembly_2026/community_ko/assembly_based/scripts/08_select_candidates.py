#!/usr/bin/env python3
"""Two-pass design, step between the passes: select the representatives that need the full-database search (pass 2)
from the pass-1 screen (148 profiles: 138 panel KOs + 10 marker KOs), and write the pass-2 chunks.

A representative is a candidate if any of its pass-1 hits has
  (a) score >= 0.5 x the KO's threshold (score = full or best-domain score by the KO's ko_list score_type).
      This includes every passing hit (score > threshold: a panel / marker KO can only be the best KO of a gene that
      passes it), the marker rescue (top-scoring profile is a marker KO at >= 0.5 x threshold), and, as a safety
      margin requested by the coordinator, panel hits at 0.5-1 x threshold;
  (b) or is one of the five no-threshold panel KOs with full E <= 1e-30 (supplement rule, 07_nothreshold_panel.py).
Why exact for every panel and marker KO: a gene's hmmsearch scores and E-values against a profile do not depend on
which other profiles are searched (same -Z, same flags), so pass 2 (all 27,981 profiles on the candidates) gives the
same best KO / passing KOs / top-scoring profile as the full run; a non-candidate passes no panel or marker KO and
cannot be rescued, so its panel/marker assignment is 'none' in both designs.
Writes pass1/hits.tsv.gz (all pass-1 hits), pass2/candidates.tsv.gz, pass2/chunks/chunk_NNNN.faa, pass2/chunk_stats.txt.
usage: 08_select_candidates.py <W> [n_chunks=180]
"""
import gzip
import os
import sys

W = sys.argv[1]
NCH = int(sys.argv[2]) if len(sys.argv) > 2 else 180
os.chdir(W)
info = {}
with open("db/hmm_info.tsv") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        r = dict(zip(hdr, line.rstrip("\n").split("\t")))
        info[r["ko"]] = (float(r["threshold"]), r["score_type"]) if r["score_type"] in ("full", "domain") else (None, "-")
kos = [l.strip() for l in open("meta/pass1_kos.txt") if l.strip()]
assert len(kos) == 148
NOTHR = {k for k in kos if info[k][0] is None}
assert NOTHR == {"K16967", "K23996", "K05915", "K18049", "K12905"}, NOTHR
MARKERS = {"K02886", "K02867", "K02982", "K02996", "K02948", "K02950", "K02926", "K02930", "K02878", "K02866"}

reason = {}                      # rep -> set of reasons
n_rows = 0
files = sorted(f for f in os.listdir("pass1/hmm") if f.endswith(".tblout"))
nch1 = sum(1 for f in os.listdir("chunks") if f.startswith("chunk_") and f.endswith(".faa"))
assert len(files) == nch1, (len(files), nch1)
with gzip.open("pass1/hits.tsv.gz.part", "wt", compresslevel=4) as fo:
    fo.write("rep\tko\tfull_evalue\tfull_score\tbest_domain_score\tscore\tthreshold\tratio\n")
    for f in files:
        last = ""
        for line in open(f"pass1/hmm/{f}"):
            if line[0] == "#":
                last = line.rstrip("\n")
                continue
            x = line.split(None, 9)
            gene, ko, e, fs, ds = x[0], x[2], float(x[4]), float(x[5]), float(x[8])
            n_rows += 1
            t, st = info[ko]
            if t is not None:
                sc = fs if st == "full" else ds
                ratio = sc / t
                if sc > t:
                    reason.setdefault(gene, set()).add("marker_pass" if ko in MARKERS else "panel_pass")
                elif ratio >= 0.5:
                    reason.setdefault(gene, set()).add("marker_half" if ko in MARKERS else "panel_half")
                fo.write(f"{gene}\t{ko}\t{e:.3g}\t{fs:.1f}\t{ds:.1f}\t{sc:.1f}\t{t:.2f}\t{ratio:.3f}\n")
            else:
                if e <= 1e-30:
                    reason.setdefault(gene, set()).add("nothr_e30")
                fo.write(f"{gene}\t{ko}\t{e:.3g}\t{fs:.1f}\t{ds:.1f}\t\t\t\n")
        assert last == "# [ok]", f
os.replace("pass1/hits.tsv.gz.part", "pass1/hits.tsv.gz")
os.makedirs("pass2/chunks", exist_ok=True)
with gzip.open("pass2/candidates.tsv.gz", "wt") as fo:
    fo.write("rep\treasons\n")
    for g in sorted(reason):
        fo.write(f"{g}\t{','.join(sorted(reason[g]))}\n")
# sequences of the candidates, in original chunk order
seqs, tot = [], 0
for f in sorted(os.listdir("chunks")):
    if not (f.startswith("chunk_") and f.endswith(".faa")):
        continue
    name, buf = None, []
    for line in open(f"chunks/{f}"):
        if line[0] == ">":
            if name in reason:
                s = "".join(buf)
                seqs.append((name, s))
                tot += len(s)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name in reason:
        s = "".join(buf)
        seqs.append((name, s))
        tot += len(s)
assert len(seqs) == len(reason), (len(seqs), len(reason))
target = tot / NCH
cid, cres, lines, out = 1, 0, [], None
cum = 0
out = open(f"pass2/chunks/chunk_{cid:04d}.faa", "w")
n_in = 0
for name, s in seqs:
    if cum >= target * cid and n_in > 0:
        out.close()
        lines.append(f"chunk_{cid:04d}\tproteins\t{n_in}\tresidues\t{cres}\n")
        cid, cres, n_in = cid + 1, 0, 0
        out = open(f"pass2/chunks/chunk_{cid:04d}.faa", "w")
    out.write(f">{name}\n")
    for k in range(0, len(s), 80):
        out.write(s[k:k + 80] + "\n")
    cres += len(s)
    cum += len(s)
    n_in += 1
out.close()
lines.append(f"chunk_{cid:04d}\tproteins\t{n_in}\tresidues\t{cres}\n")
from collections import Counter
rc = Counter(r for v in reason.values() for r in v)
strict = sum(1 for v in reason.values() if v & {"panel_pass", "marker_pass"})
with open("pass2/chunk_stats.txt", "w") as fo:
    fo.write(f"pass1_hit_rows\t{n_rows}\ncandidates\t{len(reason)}\tresidues\t{tot}\n"
             f"candidates_passing_a_panel_or_marker_KO\t{strict}\n"
             + "".join(f"reason\t{k}\t{v}\n" for k, v in sorted(rc.items()))
             + f"chunks\t{cid}\n")
    fo.writelines(lines)
print(open("pass2/chunk_stats.txt").read()[:800])
