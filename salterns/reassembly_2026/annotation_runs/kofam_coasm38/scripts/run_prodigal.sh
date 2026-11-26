#!/bin/bash
# Prodigal (-p meta -g 11, contigs >= 1,000 bp) on the 38 kept co-assembly MAGs without KOfam calls
# (coasm_final/kofam_missing.tsv -> meta/mags_38.txt), 16 parallel single-thread jobs, nice 10.
# kofam_coasm38 copy of kofam_rollout/scripts/run_prodigal.sh: input = the 38 listed MAGs read in place from
# coasm_qc/mags_kept/<mag>.fa (NOT a glob of that dir: it holds all 520 kept co-assembly MAGs), work dir
# kofam_coasm38, ALL_DONE at 38, xargs -P 16.
set -uo pipefail
W=/scratch1/afreiburger/salterns_2026
ANNOT=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_
export PRODIGAL=$ANNOT/bin/prodigal
PY=$ANNOT/bin/python3
NMAG=38
cd $W/kofam_coasm38
mkdir -p genes
date -u +"start %FT%TZ"
[ "$(sort -u meta/mags_38.txt | wc -l)" -eq $NMAG ] || { echo "meta/mags_38.txt does not list $NMAG MAGs"; exit 1; }
sed "s#^#$W/coasm_qc/mags_kept/#; s#\$#.fa#" meta/mags_38.txt | while read f; do b=$(basename $f .fa); [ -e genes/$b.done ] || echo $f; done \
  | xargs -P 16 -I{} nice -n 10 $PY scripts/call_genes.py {} genes 1000
n=$(ls genes/*.done | wc -l); echo "done files: $n"
date -u +"end %FT%TZ"
[ "$n" -eq $NMAG ] && touch genes/ALL_DONE
