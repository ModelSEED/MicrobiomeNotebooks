"""B4. Gene-calling-independent check on the contigs (MAG_contigs/<mag>.contigs__.fna, read-only): six-frame
translation of every contig, split at stop codons into ORF fragments >= 25 aa, hmmsearch with the decisive KOfam
HMMs. Answers (i) is a prs/HSK/PurU homolog present but uncalled (e.g. at a contig edge)? (ii) is Albimonas.1's
phnGHIJKL locus frameshifted (one HMM matched piecewise across frames)?"""
import sys, csv, json
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_bio import *
import pyhmmer
from Bio.Seq import Seq

KOS = ["K00948", "K00872", "K02204", "K02203", "K01433", "K06163", "K06164", "K06165", "K06166", "K05780",
       "K05781", "K06162", "K06167", "K05774", "K02041", "K02042", "K02044", "K01938", "K00772", "K00759"]
thr = {}
with open("/home/freiburger/Documents/EmilyKin/kofam_new/ko_list") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["knum"] in KOS:
            thr[r["knum"]] = (float(r["threshold"]), r["score_type"], r["definition"][:60])
alpha = pyhmmer.easel.Alphabet.amino()
hmms = []
for k in KOS:
    with pyhmmer.plan7.HMMFile(str(VB / "kofam_subset" / f"{k}.hmm")) as hf:
        hmms.append(hf.read())

def read_fasta(p):
    name, buf = None, []
    for line in open(p):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield name, "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name:
        yield name, "".join(buf)

def orfs(cid, s, minlen=25):
    L = len(s)
    out = []
    for strand, seq in (("+", s), ("-", str(Seq(s).reverse_complement()))):
        for fr in range(3):
            sub = seq[fr:]
            sub = sub[: len(sub) - len(sub) % 3]
            prot = str(Seq(sub).translate(table=11))
            pos = 0
            for piece in prot.split("*"):
                if len(piece) >= minlen:
                    a_start = fr + 3 * pos
                    a_end = a_start + 3 * len(piece)
                    if strand == "+":
                        nt = (a_start + 1, a_end)
                    else:
                        nt = (L - a_end + 1, L - a_start)
                    out.append((f"{cid}|{strand}{fr}|{nt[0]}-{nt[1]}", piece.replace("X", "")))
                pos += len(piece) + 1
    return out

targets = ["Roseovarius.9", "Roseovarius.11", "Albimonas.1"]
inv = {iid(m): m for m in model_ids()}
summary = {}
for t in targets:
    mid = inv[t]
    fa = S / "MAG_contigs" / f"{mid}__.fna"
    contigs = dict(read_fasta(fa))
    seqs = []
    for cid, s in contigs.items():
        for name, p in orfs(cid, s.upper()):
            seqs.append(pyhmmer.easel.TextSequence(name=name.encode(), sequence=p).digitize(alpha))
    print(f"===== {t}: {len(contigs)} contigs, {sum(map(len, contigs.values()))} bp, {len(seqs)} six-frame ORF fragments >=25 aa")
    res = {}
    for top in pyhmmer.hmmsearch(hmms, seqs, cpus=8, E=1e-3, Z=len(seqs)):
        ko = top.query.name.decode() if isinstance(top.query.name, bytes) else top.query.name
        for h in top:
            n = h.name.decode() if isinstance(h.name, bytes) else h.name
            d = max(h.domains, key=lambda d: d.score)
            sc = h.score if thr[ko][1] == "full" else d.score
            res.setdefault(ko, []).append((round(sc, 1), round(sc / thr[ko][0], 2), n, d.alignment.hmm_from, d.alignment.hmm_to, top.query.M if hasattr(top.query, 'M') else None))
    summary[t] = {}
    for ko in KOS:
        hs = sorted(res.get(ko, []), key=lambda x: -x[0])[:4]
        summary[t][ko] = hs
        if not hs:
            print(f"  {ko} {thr[ko][2]:60s}: no six-frame fragment with E<=1e-3")
        for h in hs:
            print(f"  {ko} {thr[ko][2]:60s}: score {h[0]:7.1f} ({h[1]}x thr)  {h[2]}  HMM {h[3]}-{h[4]}")
json.dump(summary, open(VB / "b4_sixframe.json", "w"), indent=1)
