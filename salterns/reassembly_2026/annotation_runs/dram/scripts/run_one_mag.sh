#!/usr/bin/env bash
# DRAM.py annotate on one MAG. Usage: run_one_mag.sh MAG_NAME   (skips a MAG that already has annotations.tsv)
set -uo pipefail
source /scratch1/afreiburger/salterns_2026/dram/scripts/dram_env.sh
n=$1; out=$DRAM_OUT/$n; lg=$DRAM_LOGS/mags/$n.log
if [ -s "$out/annotations.tsv" ]; then echo "$(date '+%F %T') skip $n"; exit 0; fi
rm -rf "$out"
t0=$(date +%s)
nice -n 10 /usr/bin/time -v DRAM.py annotate -i "$MAGS/$n.fna" -o "$out" --min_contig_size $MIN_CONTIG \
  --threads $THREADS --gtdb_taxonomy "$TAX" --checkm_quality "$QUAL" > "$lg" 2>&1
rc=$?
if [ $rc -eq 0 ] && [ -s "$out/annotations.tsv" ]; then
  echo "$(date '+%F %T') done $n genes=$(($(wc -l < "$out/annotations.tsv")-1)) sec=$(( $(date +%s)-t0 ))"
else
  echo "$(date '+%F %T') FAILED $n rc=$rc log=$lg"; exit 1
fi
