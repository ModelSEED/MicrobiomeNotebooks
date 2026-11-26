#!/bin/bash
# Prodigal (-p meta -g 11, contigs >= 1,000 bp) on all 1,324 MAGs, 40 parallel single-thread jobs, nice 10.
set -uo pipefail
W=/scratch1/afreiburger/salterns_2026
ANNOT=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_
export PRODIGAL=$ANNOT/bin/prodigal
PY=$ANNOT/bin/python3
cd $W/kofam
date -u +"start %FT%TZ"
ls $W/mags/*.fna | while read f; do b=$(basename $f .fna); [ -e genes/$b.done ] || echo $f; done \
  | xargs -P 40 -I{} nice -n 10 $PY scripts/call_genes.py {} genes 1000
n=$(ls genes/*.done | wc -l); echo "done files: $n"
date -u +"end %FT%TZ"
[ "$n" -eq 1324 ] && touch genes/ALL_DONE
