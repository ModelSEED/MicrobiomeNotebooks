#!/usr/bin/env python3
"""Step 1: global de-duplication of the 12,055,118 Prodigal proteins of the 39 per-sample rollout assemblies
(read_markers/genes/<s>/<s>.faa) and hmmsearch chunks of the unique sequences.

* One representative per distinct protein sequence (md5 of the sequence exactly as Prodigal wrote it, trailing '*'
  included, as in the kofam_rollout chunks), the first occurrence in sorted sample order. hmmsearch scores do not
  depend on the target set and E-values depend on it only through the fixed -Z (= all 12,055,118 genes), so
  searching one copy and copying its hits to identical proteins is exact (as kofam_rollout/make_chunks.py).
* Chunks: unique sequences in file order, cut every ~RES_PER_CHUNK residues; run through a queue (xargs -P), so
  chunk balance does not matter much.
* Writes meta/genes/<s>.tsv.gz (gene, contig, partial, aa_len, rep) for every gene, chunks/chunk_NNNN.faa,
  chunks/chunk_stats.txt (proteins = -Z). stdlib only.
usage: 01_make_chunks.py <W> <genes_dir> [res_per_chunk=4500000]
"""
import gzip
import hashlib
import os
import re
import sys

W, GD = sys.argv[1], sys.argv[2]
RES = int(sys.argv[3]) if len(sys.argv) > 3 else 4_500_000
samples = sorted(d for d in os.listdir(GD) if os.path.isfile(os.path.join(GD, d, "DONE")) and d != "coasm520")
assert len(samples) == 39, len(samples)
PART = re.compile(r"partial=(\d\d)")


def read_faa(p):
    name, hdr, seq = None, None, []
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, hdr, "".join(seq)
                hdr = line[1:].rstrip("\n")
                name, seq = hdr.split()[0], []
            else:
                seq.append(line.strip())
    if name:
        yield name, hdr, "".join(seq)


# pass 1: representatives, gene tables
rep_of = {}           # md5 digest -> representative gene name
n_all = res_all = n_uni = res_uni = 0
per_sample = []
os.makedirs(os.path.join(W, "meta/genes"), exist_ok=True)
for s in samples:
    ns = nu = 0
    tmp = os.path.join(W, "meta/genes", f"{s}.tsv.gz.part")
    with gzip.open(tmp, "wt", compresslevel=4) as fo:
        fo.write("gene\tcontig\tpartial\taa_len\trep\n")
        for name, hdr, seq in read_faa(os.path.join(GD, s, f"{s}.faa")):
            smp, rest = name.split("|", 1)
            assert smp == s, name
            contig = rest.rsplit("_", 1)[0]
            partial = PART.search(hdr).group(1)
            L = len(seq.rstrip("*"))
            h = hashlib.md5(seq.encode()).digest()
            r = rep_of.get(h)
            if r is None:
                rep_of[h] = r = name
                n_uni += 1
                res_uni += len(seq)
                nu += 1
            ns += 1
            n_all += 1
            res_all += len(seq)
            fo.write(f"{name}\t{contig}\t{partial}\t{L}\t{r}\n")
    os.replace(tmp, tmp[:-5])
    per_sample.append((s, ns, nu))
    print(s, ns, nu, flush=True)

# pass 2: chunks of the representatives
os.makedirs(os.path.join(W, "chunks"), exist_ok=True)
n_chunks_target = max(1, round(res_uni / RES))
target = res_uni / n_chunks_target
cid, cres, cn, cum = 1, 0, 0, 0
lines = []
fo = open(os.path.join(W, "chunks", f"chunk_{cid:04d}.faa.part"), "w")


def close(cid, cn, cres):
    fo.close()
    p = os.path.join(W, "chunks", f"chunk_{cid:04d}.faa")
    os.replace(p + ".part", p)
    lines.append(f"chunk_{cid:04d}\tproteins\t{cn}\tresidues\t{cres}\n")


for s in samples:
    for name, hdr, seq in read_faa(os.path.join(GD, s, f"{s}.faa")):
        if rep_of[hashlib.md5(seq.encode()).digest()] != name:
            continue
        if cum >= target * cid and cn > 0:
            close(cid, cn, cres)
            cid, cres, cn = cid + 1, 0, 0
            fo = open(os.path.join(W, "chunks", f"chunk_{cid:04d}.faa.part"), "w")
        fo.write(f">{name}\n")
        for k in range(0, len(seq), 80):
            fo.write(seq[k:k + 80] + "\n")
        cres += len(seq)
        cn += 1
        cum += len(seq)
close(cid, cn, cres)
with open(os.path.join(W, "chunks", "chunk_stats.txt"), "w") as f:
    f.write(f"samples\t{len(samples)}\nproteins\t{n_all}\nresidues\t{res_all}\n"
            f"unique_sequences\t{n_uni}\tfraction_duplicate\t{1 - n_uni / n_all:.4f}\n"
            f"searched_residues\t{res_uni}\tfraction_saved\t{1 - res_uni / res_all:.4f}\n"
            f"chunks\t{cid}\tres_per_chunk_target\t{target:.0f}\n")
    for s, ns, nu in per_sample:
        f.write(f"sample\t{s}\tgenes\t{ns}\tnew_unique\t{nu}\n")
    f.writelines(lines)
print(open(os.path.join(W, "chunks", "chunk_stats.txt")).read()[:600])
