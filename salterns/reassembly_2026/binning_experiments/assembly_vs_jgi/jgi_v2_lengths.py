#!/usr/bin/env python
"""Contig lengths of JGI assembly (ii), the IMG-team 'SPAdes version 2' reassembly (2018, IMG 3300026094-251,
Ga0209xxx), read from the @SQ headers of the Henry-lab CoverM BAMs (salterns/reads/coverm/<sample>_MG/bam_metagenome/
*.bam). Those BAMs were mapped against '<submission>.assembled.fna.filtered' (contigs >= 500 bp), so the header holds
every contig >= 500 bp of the assembly. 33 of 39 samples have a BAM.
Only the BAM header is read (pysam, no records).
Writes lengths_jgi_v2.tsv.gz: sample, contig_id, length
Usage: ~/Documents/py_venv/bin/python -I jgi_v2_lengths.py
"""
import glob, gzip, os

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
COVERM = os.path.normpath(os.path.join(HERE, "..", "..", "..", "reads", "coverm"))


def main():
    out = gzip.open(os.path.join(HERE, "lengths_jgi_v2.tsv.gz"), "wt")
    out.write("sample\tcontig_id\tlength\n")
    for d in sorted(glob.glob(os.path.join(COVERM, "*_MG"))):
        bams = glob.glob(os.path.join(d, "bam_metagenome", "*.bam"))
        if not bams:
            continue
        assert len(bams) == 1, bams
        sample = os.path.basename(d)[:-3]                 # strip '_MG'
        with pysam.AlignmentFile(bams[0], "rb") as bam:   # header only
            refs = [(n, l) for n, l in zip(bam.references, bam.lengths) if "~Ga0" in n]
        assert refs, sample
        for n, l in refs:
            out.write(f"{sample}\t{n}\t{l}\n")
        print(sample, len(refs), sum(l for _, l in refs), flush=True)
    out.close()


if __name__ == "__main__":
    main()
