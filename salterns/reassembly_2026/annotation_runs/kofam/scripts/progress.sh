#!/bin/bash
# Progress of the 40 chunk jobs: profiles done (Query: lines) -> fraction of total profile length (LENG),
# elapsed since the chunk started, and a linear ETA. usage: bash scripts/progress.sh
W=/scratch1/afreiburger/salterns_2026/kofam
cd $W
awk -F'\t' 'NR>1{c+=$2; print $5"\t"c}' db/hmm_info.tsv > /tmp/afr_kofam_cumleng.tsv
now=$(date +%s)
for c in $(seq -w 1 40); do
  if [ -e hmm/chunk_$c.done ]; then echo -e "$c\tdone"; continue; fi
  [ -e hmm/chunk_$c.failed ] && { echo -e "$c\tFAILED"; continue; }
  n=$(grep -c "^Query:" hmm/chunk_$c.out 2>/dev/null || echo 0)
  st=$(head -1 logs/chunk_$c.log | cut -d' ' -f1); s=$(date -d "$st" +%s)
  awk -v n=$n -v c=$c -v el=$((now-s)) -F'\t' '$1==(n>1?n-1:1){f=$2/16466034; printf "%s\tprofiles %d/27981\tLENG %.2f%%\telapsed %.2f h\tETA %.1f h\n", c, n, 100*f, el/3600, (f>0? el/f/3600 : -1)}' /tmp/afr_kofam_cumleng.tsv
done
rm -f /tmp/afr_kofam_cumleng.tsv
