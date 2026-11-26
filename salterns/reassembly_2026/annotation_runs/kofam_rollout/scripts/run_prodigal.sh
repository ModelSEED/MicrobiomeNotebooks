#!/bin/bash
# Prodigal (-p meta -g 11, contigs >= 1,000 bp) on the 935 kept MAGs of the Arm C rollout, 40 parallel single-thread jobs, nice 10.
# kofam_rollout copy of kofam/scripts/run_prodigal.sh: input rollout/mags_kept/*.fa (read in place), work dir kofam_rollout, ALL_DONE at 935.
set -uo pipefail
W=/scratch1/afreiburger/salterns_2026
ANNOT=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_
export PRODIGAL=$ANNOT/bin/prodigal
PY=$ANNOT/bin/python3
NMAG=935
cd $W/kofam_rollout
mkdir -p genes
date -u +"start %FT%TZ"
ls $W/rollout/mags_kept/*.fa | while read f; do b=$(basename $f .fa); [ -e genes/$b.done ] || echo $f; done \
  | xargs -P 40 -I{} nice -n 10 $PY scripts/call_genes.py {} genes 1000
n=$(ls genes/*.done | wc -l); echo "done files: $n"
date -u +"end %FT%TZ"
[ "$n" -eq $NMAG ] && touch genes/ALL_DONE
