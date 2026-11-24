#!/usr/bin/env python3
"""Prodigal gene calls for one MAG, DRAM-style (KOfam track, salterns 2026 re-assembly MAGs).

  prodigal -p meta -g 11 on contigs >= MIN_CONTIG (1,000 bp; every contig of the 1,324 MAGs is >= 1,000 bp)
  protein IDs  <MAG>|<contig>_<n>   (n = prodigal's per-contig gene number, as in run_kofam_controls.py)
  gene table   gene, mag, binner, contig, contig_len, start, end, strand, partial, aa_len, rbs_motif

stdlib only (runs with any python3 on poplar).
usage: call_genes.py <mag.fna> <out_dir> [min_contig=1000]
writes <out_dir>/<MAG>.faa, <MAG>.gff, <MAG>.genes.tsv, and <MAG>.done last.
"""
import os
import re
import subprocess
import sys
import tempfile

fna, out_dir = sys.argv[1], sys.argv[2]
min_contig = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
mag = os.path.basename(fna)
for suf in (".fna", ".fa", ".fasta"):
    if mag.endswith(suf):
        mag = mag[: -len(suf)]
        break
binner = mag.split(".")[0]          # maxbin2 / metabat2 / concoct
PRODIGAL = os.environ.get("PRODIGAL", "prodigal")

seqs, order, name = {}, [], None
with open(fna) as fh:
    for line in fh:
        if line.startswith(">"):
            name = line[1:].split()[0]
            assert name not in seqs and "|" not in name, (mag, name)
            seqs[name] = []
            order.append(name)
        else:
            seqs[name].append(line.strip())
seqs = {k: "".join(v) for k, v in seqs.items()}
keep = [k for k in order if len(seqs[k]) >= min_contig]

os.makedirs(out_dir, exist_ok=True)
base = os.path.join(out_dir, mag)
# local scratch, not NFS: TemporaryDirectory cleanup on NFS failed with ENOTEMPTY (.nfs silly-renames) for 8 MAGs
with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR", "/tmp"), prefix="afr_kofam_",
                                 ignore_cleanup_errors=True) as tmp:
    flt = os.path.join(tmp, "in.fna")
    with open(flt, "w") as fo:
        for k in keep:
            fo.write(f">{k}\n{seqs[k]}\n")
    raw_faa, gff = os.path.join(tmp, "raw.faa"), os.path.join(tmp, "out.gff")
    subprocess.run([PRODIGAL, "-p", "meta", "-g", "11", "-q", "-i", flt, "-a", raw_faa, "-f", "gff", "-o", gff],
                   check=True)
    # gene table from the GFF (ID=<seqnum>_<n>)
    rows = []
    with open(gff) as fh, open(base + ".gff.part", "w") as go:
        for line in fh:
            go.write(line)
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            n = re.search(r"ID=\d+_(\d+);", f[8]).group(1)
            partial = re.search(r"partial=(\d\d)", f[8]).group(1)
            rbs = re.search(r"rbs_motif=([^;]*)", f[8])
            start, end = int(f[3]), int(f[4])
            rows.append((f"{mag}|{f[0]}_{n}", mag, binner, f[0], len(seqs[f[0]]), start, end, f[6], partial,
                         (end - start + 1) // 3 - 1, rbs.group(1) if rbs else ""))
    nprot = 0
    with open(raw_faa) as fi, open(base + ".faa.part", "w") as fo:
        for line in fi:
            if line.startswith(">"):
                nprot += 1
                fo.write(f">{mag}|{line[1:].split()[0]}\n")
            else:
                fo.write(line)
    assert nprot == len(rows), (mag, nprot, len(rows))
    with open(base + ".genes.tsv.part", "w") as fo:
        fo.write("gene\tmag\tbinner\tcontig\tcontig_len\tstart\tend\tstrand\tpartial\taa_len\trbs_motif\n")
        for r in rows:
            fo.write("\t".join(map(str, r)) + "\n")
for suf in (".faa", ".gff", ".genes.tsv"):
    os.replace(base + suf + ".part", base + suf)
with open(base + ".done", "w") as fo:
    fo.write(f"{mag}\tcontigs_in={len(order)}\tcontigs_used={len(keep)}\tgenes={len(rows)}\n")
