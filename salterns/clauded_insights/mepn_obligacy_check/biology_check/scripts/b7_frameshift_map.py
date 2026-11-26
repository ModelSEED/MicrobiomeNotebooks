"""B7. Frameshift map of the phnGHIJKL locus: for each C-P lyase HMM, every domain hit on the six-frame ORF
fragments of the locus contig, with exact nucleotide coordinates of the aligned segment, frame, HMM range and score.
Colinear segments that jump frame while the HMM coordinate keeps increasing = a gene broken by frameshifts.
Albimonas.1 NODE_10693 vs the intact loci of Roseovarius.9 (NODE_52) and Roseovarius.11 (NODE_6588) as controls."""
import sys, csv
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_bio import *
import pyhmmer
from Bio.Seq import Seq

KOS = {"K06166": "PhnG", "K06165": "PhnH", "K06164": "PhnI", "K06163": "PhnJ", "K05781": "PhnK", "K05780": "PhnL", "K06162": "PhnM"}
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
            if name: yield name, "".join(buf)
            name, buf = line[1:].split()[0], []
        else: buf.append(line)
    if name: yield name, "".join(buf)

def frame_orfs(s, minlen=15):
    L = len(s); out = []
    for strand, seq in (("+", s), ("-", str(Seq(s).reverse_complement()))):
        for fr in range(3):
            sub = seq[fr:]; sub = sub[: len(sub) - len(sub) % 3]
            prot = str(Seq(sub).translate(table=11)); pos = 0
            for piece in prot.split("*"):
                if len(piece) >= minlen:
                    out.append((strand, fr, fr + 3 * pos, piece))
                pos += len(piece) + 1
    return out, L

inv = {iid(m): m for m in model_ids()}
for t, contig in [("Albimonas.1", "NODE_10693_length_5378_cov_3.50314"), ("Roseovarius.9", "NODE_52_length_131214_cov_6.03563"),
                  ("Roseovarius.11", "NODE_6588_length_6200_cov_1.47456")]:
    seqs = dict(read_fasta(S / "MAG_contigs" / f"{inv[t]}__.fna"))
    s = seqs[contig].upper()
    if t == "Roseovarius.9":
        s = s[:12000]   # the phn cluster sits at the start of this 131 kb contig (CDS.5-13, 3.3-10.2 kb)
    orfs, L = frame_orfs(s)
    dseqs, meta = [], {}
    for i, (strand, fr, a0, p) in enumerate(orfs):
        n = f"o{i}"
        meta[n] = (strand, fr, a0, len(p))
        dseqs.append(pyhmmer.easel.TextSequence(name=n.encode(), sequence=p.replace("X", "")).digitize(alpha))
    rows = []
    for top in pyhmmer.hmmsearch(hmms, dseqs, E=1e-2, Z=len(dseqs)):
        ko = top.query.name.decode() if isinstance(top.query.name, bytes) else top.query.name
        for h in top:
            n = h.name.decode() if isinstance(h.name, bytes) else h.name
            strand, fr, a0, plen = meta[n]
            for d in h.domains:
                if d.score < 15:
                    continue
                ali = d.alignment
                nt_from = a0 + 3 * (ali.target_from - 1) + 1
                nt_to = a0 + 3 * ali.target_to
                if strand == "-":
                    nt_from, nt_to = L - nt_to + 1, L - nt_from + 1
                rows.append((KOS[ko], nt_from, nt_to, f"{strand}{fr}", ali.hmm_from, ali.hmm_to, top.query.M, round(d.score, 1)))
    print(f"===== {t} {contig} (analysed length {L} bp)")
    for g in KOS.values():
        rs = sorted([r for r in rows if r[0] == g], key=lambda r: r[1])
        if not rs:
            print(f"  {g}: no domain >= 15 bits"); continue
        frames = sorted({r[3] for r in rs})
        cov = sum(r[5] - r[4] + 1 for r in rs)
        print(f"  {g}: {len(rs)} segment(s) in frame(s) {frames}; HMM positions covered (sum) {cov}/{rs[0][6]}")
        for r in rs:
            print(f"      nt {r[1]:>6}-{r[2]:<6} frame {r[3]}  HMM {r[4]:>3}-{r[5]:<3}/{r[6]}  {r[7]:6.1f} bits")
