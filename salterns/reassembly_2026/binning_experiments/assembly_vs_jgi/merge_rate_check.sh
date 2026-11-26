#!/usr/bin/env bash
# How many R1_A_D1 read pairs overlap and merge? (Zhou et al. 2022: "reads overlapped and merged" before the
# Tadpole + SPAdes re-assembly; the 2026 KBase run fed unmerged Trimmomatic pairs to metaSPAdes.)
# 300,000 pairs (reads 100,001-400,000; the first 100k are skipped because tile-edge reads carry N runs) from the
# first 40 MB of the ENA FASTQ files of SRR5209943, merged with BBMerge 39.13 (BBTools tarball, not installed
# system-wide; Java 21). Two settings: BBMerge defaults, and bbmerge-auto with Tadpole extension and error
# correction (rem k=62 extend2=50 ecct), the BBTools assembly-pipeline style. On a 300k-pair subsample the Tadpole
# extension has little k-mer coverage, so the second number is a lower bound for the full read set.
# Results (2026-11-11) are in merge_rate_R1_A_D1.txt.
# Usage: bash merge_rate_check.sh <scratch dir>     (downloads ~95 MB into <scratch dir>)
set -euo pipefail
S=${1:?scratch dir}
mkdir -p "$S/reads" "$S/bbtools"
cd "$S/bbtools"
[ -x bbmap/bbmerge.sh ] || { curl -sL -o BBMap.tar.gz "https://sourceforge.net/projects/bbmap/files/BBMap_39.13.tar.gz/download"; tar -xzf BBMap.tar.gz; }
cd "$S/reads"
for i in 1 2; do
  [ -s part_$i.gz ] || curl -s -r 0-40000000 -o part_$i.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR520/003/SRR5209943/SRR5209943_$i.fastq.gz
  (gzip -dc part_$i.gz 2>/dev/null || true) | head -n 1600000 | tail -n 1200000 > sub_$i.fq
done
BB="$S/bbtools/bbmap"
nice -n 10 "$BB/bbmerge.sh" in1=sub_1.fq in2=sub_2.fq ihist=ihist_default.txt t=8 -Xmx4g
nice -n 10 "$BB/bbmerge-auto.sh" in1=sub_1.fq in2=sub_2.fq ihist=ihist_rem.txt rem k=62 extend2=50 ecct t=8 -Xmx8g
