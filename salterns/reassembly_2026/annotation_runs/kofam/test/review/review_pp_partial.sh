#!/bin/bash
# Review test: run postprocess.py on a truncated copy of a running chunk's real full-DB output (complete KOs only).
set -euo pipefail
W=/scratch1/afreiburger/salterns_2026/kofam
T=$W/test/review_partial_pp
P=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin/python3
c=${1:-01}
rm -rf $T; mkdir -p $T/hmm
for d in db meta chunks genes; do ln -s $W/$d $T/$d; done
# last complete KO = the second-to-last distinct query in each table
lastko_t=$(grep -v '^#' $W/hmm/chunk_$c.tblout | awk '{print $3}' | uniq | tail -2 | head -1)
lastko_d=$(grep -v '^#' $W/hmm/chunk_$c.domtblout | awk '{print $4}' | uniq | tail -2 | head -1)
cut=$(printf "%s\n%s\n" $lastko_t $lastko_d | sort | head -1)
echo "keep KOs <= $cut"
grep -v '^#' $W/hmm/chunk_$c.tblout | awk -v k=$cut '$3<=k' > $T/hmm/chunk_$c.tblout || true
grep -v '^#' $W/hmm/chunk_$c.domtblout | awk -v k=$cut '$4<=k' > $T/hmm/chunk_$c.domtblout || true
wc -l $T/hmm/chunk_$c.*
/usr/bin/time -v nice -n 10 $P $W/scripts/postprocess.py chunk $T $c 2> $T/time.txt
grep -E "Elapsed|Maximum resident" $T/time.txt
ls -la $T/results/chunks; ls $T/results/per_mag | wc -l
