#!/usr/bin/env bash
# 1. Copy the Bakta 6.0-full DB to poplar's local disk in parallel and verify it (copy_db_parallel.sh: sizes for every
#    file, byte-for-byte cmp for files > 2 GB).
# 2. Smoke test: Bakta on the smallest representative alone.
# 3. Launch run_bakta_reps.sh on all 1,184 representatives (16 jobs x 6 threads).
# Usage (poplar): cd /scratch1/afreiburger/wetland_bakta && setsid nohup bash launch_after_db_copy.sh > logs/launcher.log 2>&1 < /dev/null &
set -uo pipefail
SRC=/scratch1/fliu/hub_scratch/data/bakta/db
DST=/scratch/afreiburger/bakta_db_6.0_full
B=/scratch1/afreiburger/wetland_bakta
ev() { echo "[$(date '+%F %T')] $*" >> "$B/bakta.events"; }
cd "$B"
ev "LAUNCHER copying the DB to $DST (parallel)"
t0=$(date +%s)
if bash copy_db_parallel.sh "$SRC" "$DST" >> logs/db_copy.log 2>&1; then
  ev "LAUNCHER DB copy verified in $(( $(date +%s) - t0 )) s ($(du -sh "$DST" | cut -f1); $(tail -1 logs/db_copy.log))"
else
  ev "LAUNCHER FAIL: DB copy not verified (logs/db_copy.log); not launching"
  exit 1
fi
m=$(while read -r x; do printf "%s %s\n" "$(stat -c %s "/scratch1/afreiburger/wetland_coverm_pilot/all_mags/$x.fna")" "$x"; done < reps.txt | sort -n | head -1 | cut -d" " -f2)
echo "$m" > smoke.txt
REPS=smoke.txt bash run_bakta_reps.sh 1 16 >> logs/driver.log 2>&1
if [ -e "runs/$m/DONE" ]; then
  ev "LAUNCHER smoke test passed on $m ($(cat "runs/$m/DONE") s); launching all"
  exec bash run_bakta_reps.sh 16 6 >> logs/driver.log 2>&1
else
  ev "LAUNCHER FAIL: smoke test on $m failed (runs/$m.FAILED); not launching"
  exit 1
fi
