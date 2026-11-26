"""B12. Roseovarius.9 SHMT (glyA, K00600): the 4 RAST 'Serine hydroxymethyltransferase' CDSs are 60-68 aa; where do
they sit on the HMM and on their contigs, and is a full-length glyA present anywhere in the six-frame translation?"""
import sys, csv
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_bio import *
import pyhmmer
from Bio.Seq import Seq
alpha = pyhmmer.easel.Alphabet.amino()
with pyhmmer.plan7.HMMFile(str(VB / "kofam_subset" / "K00600.hmm")) as hf:
    hmm = hf.read()
inv = {iid(m): m for m in model_ids()}
for t in ["Roseovarius.9", "Albimonas.1", "Roseovarius.8"]:
    mid = inv[t]
    g = rast(mid)
    seqs = {}
    name = None
    for line in open(S / "MAG_contigs" / f"{mid}__.fna"):
        line = line.rstrip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        else:
            seqs[name].append(line)
    seqs = {k: "".join(v).upper() for k, v in seqs.items()}
    print(f"===== {t}  (K00600 HMM length {hmm.M})")
    prots = []
    for f in g["features"]:
        if any("hydroxymethyltransferase" in x.lower() for x in (f.get("functions") or [])):
            cid, st, sd, ln = f["location"][0]
            prots.append(pyhmmer.easel.TextSequence(name=short(f["id"]).encode(), sequence=f["protein_translation"]).digitize(alpha))
            print(f"  RAST {short(f['id'])}: {len(f['protein_translation'])} aa, contig {cid} (len {len(seqs[cid])}), start {st} {sd}, gene length {ln} nt")
    for top in pyhmmer.hmmsearch([hmm], prots, E=10):
        for h in top:
            for d in h.domains:
                print(f"     {h.name.decode() if isinstance(h.name, bytes) else h.name}: HMM {d.alignment.hmm_from}-{d.alignment.hmm_to}, {d.score:.1f} bits")
    # six-frame, whole MAG
    orfs = []
    for cid, s in seqs.items():
        L = len(s)
        for strand, sq in (("+", s), ("-", str(Seq(s).reverse_complement()))):
            for fr in range(3):
                sub = sq[fr:]; sub = sub[: len(sub) - len(sub) % 3]
                pos = 0
                for piece in str(Seq(sub).translate(table=11)).split("*"):
                    if len(piece) >= 40:
                        orfs.append(pyhmmer.easel.TextSequence(name=f"{cid}|{strand}{fr}|{fr+3*pos}".encode(), sequence=piece.replace("X", "")).digitize(alpha))
                    pos += len(piece) + 1
    best = []
    for top in pyhmmer.hmmsearch([hmm], orfs, E=1e-5, Z=len(orfs), cpus=8):
        for h in top:
            d = max(h.domains, key=lambda d: d.score)
            best.append((h.score, h.name.decode() if isinstance(h.name, bytes) else h.name, d.alignment.hmm_from, d.alignment.hmm_to))
    for sc, n, a, b in sorted(best, reverse=True)[:6]:
        print(f"  six-frame: {sc:7.1f} bits (KOfam thr 84.6 full)  {n}  HMM {a}-{b}")
