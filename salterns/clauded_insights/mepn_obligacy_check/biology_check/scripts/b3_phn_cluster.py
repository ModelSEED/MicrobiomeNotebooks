"""B3. Gene structure of the C-P lyase (phn) cluster in each of the 3 obligate MAGs (plus 2 non-obligate references):
CDS order, contig coordinates, protein length, RAST function, best C-P-lyase KOfam hit with HMM coverage.
A split (frameshifted / truncated) gene shows as two adjacent CDSs hitting different halves of one HMM."""
import sys, re
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_bio import *
import pyhmmer

PHN_KOS = ["K06162", "K06163", "K06164", "K06165", "K06166", "K06167", "K05780", "K05781", "K02041", "K02042",
           "K02044", "K05774"]
LAB = {"K06162": "PhnM", "K06163": "PhnJ", "K06164": "PhnI", "K06165": "PhnH", "K06166": "PhnG", "K06167": "PhnP",
       "K05780": "PhnL", "K05781": "PhnK", "K02041": "PhnC", "K02042": "PhnE", "K02044": "PhnD", "K05774": "PhnN"}
import csv
thr = {}
with open("/home/freiburger/Documents/EmilyKin/kofam_new/ko_list") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["knum"] in PHN_KOS:
            thr[r["knum"]] = (float(r["threshold"]), r["score_type"])
alpha = pyhmmer.easel.Alphabet.amino()
hmms = {}
for k in PHN_KOS:
    with pyhmmer.plan7.HMMFile(str(VB / "kofam_subset" / f"{k}.hmm")) as hf:
        hmms[k] = hf.read()
rx = re.compile(r"phosphonate|\bPhn[A-Z]\b|C-P lyase|1,2-cyclic|ribose 1,5-bisphosph", re.I)
targets = ["Roseovarius.9", "Roseovarius.11", "Albimonas.1", "Roseovarius.8", "JAVJPJ01.1"]
inv = {iid(m): m for m in model_ids()}
for t in targets:
    mid = inv[t]
    g = rast(mid)
    feats = g["features"]
    idx = {short(f["id"]): i for i, f in enumerate(feats)}
    sel = [f for f in feats if any(rx.search(x) for x in (f.get("functions") or []))]
    # widen to +-2 neighbours to catch unannotated fragments
    keep = set()
    for f in sel:
        i = idx[short(f["id"])]
        for j in range(max(0, i - 2), min(len(feats), i + 3)):
            keep.add(j)
    seqs = []
    for j in sorted(keep):
        f = feats[j]
        if f.get("protein_translation"):
            seqs.append(pyhmmer.easel.TextSequence(name=short(f["id"]).encode(), sequence=f["protein_translation"].replace("*", "")).digitize(alpha))
    best = {}
    for k, h in hmms.items():
        for top in pyhmmer.hmmsearch([h], seqs, E=10, Z=len(seqs) or 1):
            for hit in top:
                n = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                d = max(hit.domains, key=lambda d: d.score)
                sc = hit.score if thr[k][1] == "full" else d.score
                rec = (sc, k, round(sc / thr[k][0], 2), d.alignment.hmm_from, d.alignment.hmm_to, h.M)
                if n not in best or rec[0] > best[n][0]:
                    best[n] = rec
    print(f"===== {t}  ({mid})")
    for j in sorted(keep):
        f = feats[j]
        cid, start, strand, ln = f["location"][0]
        plen = len(f.get("protein_translation") or "")
        fn = "; ".join(f.get("functions") or [])
        b = best.get(short(f["id"]))
        bs = f"{LAB[b[1]]}({b[1]}) {b[2]}x thr, HMM {b[3]}-{b[4]}/{b[5]}" if b else "-"
        print(f"  {short(f['id']):9s} {cid[:34]:34s} {start:>8} {strand} {plen:4d} aa | {bs:38s} | {fn[:95]}")
