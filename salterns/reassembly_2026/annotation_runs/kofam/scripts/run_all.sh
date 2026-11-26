#!/bin/bash
# Launch the 40 single-CPU hmmsearch chunk jobs in parallel (detached driver), then combine.
set -uo pipefail
W=/scratch1/afreiburger/salterns_2026/kofam
ANNOT=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_
cd $W
mkdir -p results/chunks
echo "$(date -u +%FT%TZ) run_all start pid $$" >> logs/run_all.log
: > logs/chunk_pids.tsv
for c in $(seq -w 1 40); do
    bash scripts/run_chunk.sh $c &
    echo -e "$c\t$!" >> logs/chunk_pids.tsv
done
wait
nd=$(ls hmm/chunk_*.done 2>/dev/null | wc -l); nf=$(ls hmm/chunk_*.failed 2>/dev/null | wc -l)
echo "$(date -u +%FT%TZ) chunks done $nd failed $nf" >> logs/run_all.log
$ANNOT/bin/python3 scripts/postprocess.py combine $W >> logs/run_all.log 2>&1
[ "$nd" -eq 40 ] && touch ALL_DONE
echo "$(date -u +%FT%TZ) run_all end" >> logs/run_all.log
