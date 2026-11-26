#!/bin/bash
# Pass 2 driver: every candidate chunk x all 27,981 profiles through xargs -P NP; then (all postprocessed, none
# failed) 03_aggregate.py -> ALL_DONE. Restart-safe. Launch:
#   nohup setsid bash scripts/run_pass2.sh > logs/run_pass2.nohup 2>&1 < /dev/null &
set -uo pipefail
W=/scratch1/afreiburger/salterns_2026/community_ko
A=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin
NP=${NP:-90}
cd $W
export TMPDIR=$W/tmp
mkdir -p pass2/hmm pass2/results/chunks logs/pass2
[ -e pass2/db ] || ln -s ../db pass2/db
[ -e pass2/meta ] || ln -s ../meta pass2/meta
echo "$(date -u +%FT%TZ) run_pass2 start pid $$ NP $NP load $(cut -d' ' -f1-3 /proc/loadavg)" >> logs/run_pass2.log
ls pass2/chunks/ | sed -n 's/^chunk_\([0-9]*\)\.faa$/\1/p' | xargs -P $NP -I{} bash scripts/p2_run_chunk.sh {}
nch=$(ls pass2/chunks/ | grep -c '^chunk_[0-9]*\.faa$')
nd=$(ls pass2/results/chunks/ | grep -c '\.postprocessed$')
nf=$(ls pass2/hmm/ pass2/results/chunks/ | grep -c 'failed$')
echo "$(date -u +%FT%TZ) pass2 finished: chunks $nch postprocessed $nd failed $nf" >> logs/run_pass2.log
if [ "$nd" -eq "$nch" ] && [ "$nf" -eq 0 ]; then
    touch pass2/PASS2_DONE
    echo "$(date -u +%FT%TZ) aggregate start" >> logs/run_pass2.log
    nice -n 10 $A/python3 -I scripts/03_aggregate.py $W >> logs/03_aggregate.log 2>&1 && touch ALL_DONE
    echo "$(date -u +%FT%TZ) aggregate end rc=$?" >> logs/run_pass2.log
fi
echo "$(date -u +%FT%TZ) run_pass2 end" >> logs/run_pass2.log
