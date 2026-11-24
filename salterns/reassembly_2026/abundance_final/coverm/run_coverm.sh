#!/usr/bin/env bash
# Abundance of the 585 final-set species representatives in all 39 salterns metagenomes (2026-11-20, v2).
# v1 failed: CoverM concatenated + indexed the genomes per sample in /tmp (TMPDIR unset) and filled poplar's shared /tmp.
# v2: one concatenated reference (contigs renamed '<genome>~<contig>') and one minimap2 -x sr index, both on /scratch1;
# CoverM genome mode with --separator '~' and the pre-built index; TMPDIR on /scratch1 per sample; /tmp is never used.
# Settings: minimap2-sr, >= 95 % identity, >= 75 % aligned, min-covered-fraction 0, methods count/covered_fraction/mean.
# Relative abundance (computed afterwards) = count / read_counts.json['qtrim']['result_reads'] (QC'd map_* reads).
set -uo pipefail
A=/scratch1/afreiburger/salterns_2026/abundance_final; cd $A
B=/scratch1/afreiburger/emilykin/processed/.snakemake_envs
CM=$(ls -d $B/8c3854a8*_)/bin; MM=$(ls -d $B/d501a828*_)/bin
export PATH=$CM:$MM:$PATH
mkdir -p out logs tmp ref
coverm --version > coverm.version 2>&1; minimap2 --version >> coverm.version 2>&1
if [ ! -s ref/reps585.mmi ]; then
  : > ref/reps585.fa.part
  for g in genomes/*.fa; do n=$(basename $g .fa); sed "s/^>/>${n}~/" $g >> ref/reps585.fa.part; done
  mv ref/reps585.fa.part ref/reps585.fa
  [ "$(grep -c '^>' ref/reps585.fa)" = "$(grep -c '^>' genomes/*.fa | awk -F: '{s+=$2} END{print s}')" ] || { echo "contig count mismatch"; exit 1; }
  [ "$(grep '^>' ref/reps585.fa | sort | uniq -d | wc -l)" = 0 ] || { echo "duplicate contig names"; exit 1; }
  TMPDIR=$A/tmp nice -n 10 minimap2 -x sr -t 8 -d ref/reps585.mmi.part ref/reps585.fa > logs/minimap2_index.log 2>&1 && mv ref/reps585.mmi.part ref/reps585.mmi
fi
R=/scratch1/afreiburger/salterns_2026/rollout/reads
run1() { s=$1; [ -s out/$s.tsv ] && return 0
  mkdir -p tmp/$s
  TMPDIR=$A/tmp/$s nice -n 10 coverm genome --coupled $R/$s/map_1.fq.gz $R/$s/map_2.fq.gz \
     --reference ref/reps585.mmi --minimap2-reference-is-index --separator '~' \
     --mapper minimap2-sr --methods count covered_fraction mean --min-read-percent-identity 95 --min-read-aligned-percent 75 \
     --min-covered-fraction 0 --threads 16 --output-file out/$s.tsv.part > logs/$s.log 2>&1 \
     && mv out/$s.tsv.part out/$s.tsv || echo "$s FAILED" >> FAILED.txt
  rm -rf tmp/$s; }
export -f run1; export R A
cat ${1:-libs39.txt} | xargs -P ${2:-4} -I{} bash -c 'run1 {}'
[ "$(ls out/*.tsv 2>/dev/null | wc -l)" = 39 ] && touch DONE
