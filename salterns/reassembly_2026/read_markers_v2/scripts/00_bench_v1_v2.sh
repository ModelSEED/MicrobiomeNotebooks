#!/usr/bin/env bash
# DIAMOND benchmark (v2): the same 8 M reads (first 8 M of R1_A_D1 map_1) against the v1 and the v2 survey DB, run at
# the same time under the same load, production flags (--mid-sensitive -k K -e 1e-3 -b 4 -c 1), 16 threads each.
# Usage: bash 00_bench_v1_v2.sh <tag> <db> <k>
set -euo pipefail
W=/scratch1/afreiburger/salterns_2026/read_markers_v2
DM=/scratch1/afreiburger/salterns_2026/envs/bt_dastool/bin/diamond
export TMPDIR=$W/tmp
tag=$1; db=$2; k=$3
cd $W/bench
/usr/bin/time -v -o bench_$tag.time nice -n 10 $DM blastx -d $db -q q8M.fq.gz -o bench_$tag.tsv.gz --compress 1 --mid-sensitive \
  -k $k -e 1e-3 --threads 16 -b 4 -c 1 -t $TMPDIR --quiet \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
echo done > bench_$tag.done
