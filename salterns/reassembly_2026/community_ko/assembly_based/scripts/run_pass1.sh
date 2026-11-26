#!/bin/bash
# Pass 1 driver: all 521 chunks x 148 profiles through xargs -P NP; then 08_select_candidates.py (pass-2 chunks).
# Launch: nohup setsid bash scripts/run_pass1.sh > logs/run_pass1.nohup 2>&1 < /dev/null &
set -uo pipefail
W=/scratch1/afreiburger/salterns_2026/community_ko
A=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin
NP=${NP:-90}
cd $W
export TMPDIR=$W/tmp
mkdir -p pass1/hmm logs/pass1
echo "$(date -u +%FT%TZ) run_pass1 start pid $$ NP $NP load $(cut -d' ' -f1-3 /proc/loadavg)" >> logs/run_pass1.log
ls chunks/ | sed -n 's/^chunk_\([0-9]*\)\.faa$/\1/p' | xargs -P $NP -I{} bash scripts/p1_run_chunk.sh {}
nch=$(ls chunks/ | grep -c '^chunk_[0-9]*\.faa$'); nd=$(ls pass1/hmm | grep -c '\.done$'); nf=$(ls pass1/hmm | grep -c 'failed$')
echo "$(date -u +%FT%TZ) pass1 search finished: chunks $nch done $nd failed $nf" >> logs/run_pass1.log
if [ "$nd" -eq "$nch" ] && [ "$nf" -eq 0 ]; then
    touch pass1/PASS1_DONE
    nice -n 10 $A/python3 -I scripts/08_select_candidates.py $W >> logs/08_select_candidates.log 2>&1 && touch pass1/SELECT_DONE
    echo "$(date -u +%FT%TZ) select rc=$?" >> logs/run_pass1.log
fi
echo "$(date -u +%FT%TZ) run_pass1 end" >> logs/run_pass1.log
