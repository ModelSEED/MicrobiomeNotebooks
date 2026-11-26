"""B8. Context for the Albimonas.1 phn locus: (i) contig coverage (from the SPAdes header) relative to the MAG's
length-weighted distribution; (ii) homopolymer runs (>=5 nt) at the frameshift junctions found in B7."""
import sys, re, statistics
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_bio import *

def read_fasta(p):
    name, buf = None, []
    for line in open(p):
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield name, "".join(buf)
            name, buf = line[1:].split()[0], []
        else: buf.append(line)
    if name: yield name, "".join(buf)

inv = {iid(m): m for m in model_ids()}
for t in ["Albimonas.1", "Roseovarius.9", "Roseovarius.11"]:
    seqs = dict(read_fasta(S / "MAG_contigs" / f"{inv[t]}__.fna"))
    covs = []
    for n, s in seqs.items():
        m = re.search(r"cov_([0-9.]+)", n)
        covs.append((float(m.group(1)), len(s), n))
    tot = sum(l for _, l, _ in covs)
    wmean = sum(c * l for c, l, _ in covs) / tot
    srt = sorted(covs)
    acc, wmed = 0, None
    for c, l, _ in srt:
        acc += l
        if acc >= tot / 2:
            wmed = c; break
    print(f"{t}: {len(covs)} contigs, length-weighted mean cov {wmean:.2f}, length-weighted median {wmed:.2f}")
    for key in ["NODE_10693_", "NODE_2963_", "NODE_190_", "NODE_52_", "NODE_6588_", "NODE_8872_", "NODE_7647_"]:
        for c, l, n in covs:
            if n.startswith(key):
                below = sum(ll for cc, ll, _ in covs if cc < c) / tot
                print(f"   {n}: cov {c}, {l} bp, {100*below:.1f}% of MAG bp on lower-coverage contigs")
seq = dict(read_fasta(S / "MAG_contigs" / f"{inv['Albimonas.1']}__.fna"))["NODE_10693_length_5378_cov_3.50314"].upper()
print("\nAlbimonas.1 NODE_10693 frameshift junctions (from b7): +-15 nt, homopolymers >=5 in caps context")
for lab, pos in [("PhnH +1->+2", 1192), ("PhnI +1->+0", 1684), ("PhnI +0->+1", 2178), ("PhnI +1->+2", 2241),
                 ("PhnJ stop in +1", 2550), ("PhnJ +1->+0", 2628), ("PhnJ +0->+2", 2946), ("PhnJ +2->+1", 3116)]:
    w = seq[pos - 16: pos + 15]
    hp = [m.group(0) for m in re.finditer(r"A{5,}|C{5,}|G{5,}|T{5,}", w)]
    print(f"   {lab:16s} ~{pos}: {w}  homopolymers>=5: {hp or 'none'}")
gc = (seq.count("G") + seq.count("C")) / len(seq)
print(f"   NODE_10693 GC {gc:.3f}; homopolymer runs >=5 nt in whole contig: {len(re.findall(r'A{5,}|C{5,}|G{5,}|T{5,}', seq))}")
