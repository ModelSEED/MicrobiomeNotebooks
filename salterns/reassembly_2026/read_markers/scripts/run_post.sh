#!/usr/bin/env bash
# Per-sample post-processing as soon as both mates of a sample are through DIAMOND:
#   08_classify.py sample (1 CPU) -> results/per_sample/<s>.{counts.tsv,reads.tsv.gz}
#   09_map_assigned.sh <s> (bowtie2/samtools, 4 threads)  -> map/<s>.hits.tsv.gz
# Thread budget: a post job (5 threads) starts only while 8 x (DIAMOND slots still needed) + 5 x (running post
# jobs) + 5 <= 64, so the track never exceeds 64 threads.
set -uo pipefail
W=/scratch1/afreiburger/salterns_2026/read_markers
A=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin
cd $W; mkdir -p results/per_sample logs/post
SAMPLES="$@"
post() {
  s=$1
  nice -n 10 $A/python3 -I scripts/08_classify.py sample cutoffs.tsv results/per_sample/$s dmnd/${s}_1.tsv.gz dmnd/${s}_2.tsv.gz \
    > logs/post/$s.classify.log 2>&1 || { echo "$s classify FAILED" > logs/post/$s.FAILED; return 1; }
  T=4 bash scripts/09_map_assigned.sh $s > logs/post/$s.map.log 2>&1 || { echo "$s map FAILED" > logs/post/$s.FAILED; return 1; }
  touch logs/post/$s.done
}
# DIAMOND slots still needed = samples (of all 39) with a mate not yet through DIAMOND, at most 8 (6 + 1 + 1 drivers);
# counted from the done markers, not from running processes, so the short gaps between mates do not open the budget.
ALL39=$(ls /scratch1/afreiburger/salterns_2026/rollout/reads)
ndm() { n=0; for x in $ALL39; do { [ -f dmnd/${x}_1.done ] && [ -f dmnd/${x}_2.done ]; } || n=$((n+1)); done; [ $n -gt 8 ] && n=8; echo $n; }
declare -A started
while true; do
  left=0
  for s in $SAMPLES; do
    [ -f logs/post/$s.done ] && continue
    [ -f logs/post/$s.FAILED ] && continue
    left=$((left+1))
    [ -n "${started[$s]:-}" ] && continue
    [ -f dmnd/${s}_1.done ] && [ -f dmnd/${s}_2.done ] || continue
    npost=$(jobs -rp | wc -l)
    if [ $((8 * $(ndm) + 5 * npost + 5)) -le 64 ]; then
      post $s &
      started[$s]=1
      echo "$(date -u +%FT%TZ) started $s" >> logs/post/driver.log
    fi
  done
  [ $left -eq 0 ] && break
  sleep 60
done
wait
echo "$(date -u +%FT%TZ) all post-processing finished" >> logs/post/driver.log
