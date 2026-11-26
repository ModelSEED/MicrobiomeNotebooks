#!/usr/bin/env bash
# Step 6 (v2: v1 script, v2 paths; run as 6 jobs x 16 threads): DIAMOND 2.2.8 blastx of every QC read (map_1 and map_2 of each sample, each mate a separate query) against
# db/survey.dmnd. Permissive output (E <= 1e-3, top 5 targets per read); cut-offs are applied in 07_classify.py.
# Usage: bash 06_diamond_reads.sh [samples...]   env: J (concurrent samples, default 4), T (threads each, default 16),
#        MODE (DIAMOND sensitivity flag, default --mid-sensitive), DB (default db/survey), OUT (default dmnd)
set -euo pipefail
W=/scratch1/afreiburger/salterns_2026/read_markers_v2
R=/scratch1/afreiburger/salterns_2026/rollout
DM=/scratch1/afreiburger/salterns_2026/envs/bt_dastool/bin/diamond
export TMPDIR=$W/tmp
J=${J:-4}; T=${T:-16}; MODE=${MODE:---mid-sensitive}; DB=${DB:-db/survey}; OUT=${OUT:-dmnd}
cd $W; mkdir -p $OUT logs
export W R DM T MODE DB OUT
one() {
  s=$1
  for m in 1 2; do
    o=$OUT/${s}_$m
    [ -f $o.done ] && continue
    rm -f $o.tsv.gz
    /usr/bin/time -v -o $o.time nice -n 10 $DM blastx -d $DB -q $R/reads/$s/map_$m.fq.gz -o $o.tsv.gz --compress 1 \
      $MODE -k 5 -e 1e-3 --threads $T -b 4 -c 1 -t $TMPDIR --quiet \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
      && touch $o.done || { echo "$(date -u +%FT%TZ) $s mate $m FAILED" >> logs/06_diamond.log; return 1; }
    echo "$(date -u +%FT%TZ) $s mate $m done" >> logs/06_diamond.log
  done
}
export -f one
for s in "$@"; do echo $s; done | xargs -P $J -I{} bash -c 'one {}'
echo "$(date -u +%FT%TZ) all requested samples finished" >> logs/06_diamond.log
