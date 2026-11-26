#!/usr/bin/env bash
# Driver: run every inputs/chunk_* through run_dram2.sh, ONE AT A TIME (each launch is capped at
# 64 CPUs by launch/wetland_dram2.config, so the total never exceeds 64).
# Restartable: chunks with runs/<chunk>/DONE are skipped; anything else is (re)launched with
# -resume, so finished tasks are reused.  A failed chunk is retried once (resume), then the
# driver moves on and reports it at the end.
#
#   setsid nohup bin/run_all_chunks.sh > logs/driver.log 2>&1 < /dev/null &
#
# CLEAN_WORK=1 (default) deletes a chunk's local work dir once its DONE marker is written
# (published outputs are copies; -resume is not needed after DONE).
set -uo pipefail
source /scratch1/afreiburger/wetland_dram2/bin/dram2_settings.sh
CLEAN_WORK=${CLEAN_WORK:-1}
cd "$B"

exec 9> "$B/logs/driver.lock"
flock -n 9 || { echo "[$(date -Is)] another driver holds logs/driver.lock; exiting"; exit 1; }
echo "[$(date -Is)] driver pid=$$ start; dbs='$DRAM2_DBS'"

failed=()
for dir in "$B"/inputs/chunk_*/; do
  c=$(basename "$dir")
  if [ -e "runs/$c/DONE" ]; then echo "[$(date -Is)] $c already DONE, skip"; continue; fi
  if [ -e "runs/$c/LOWDISK" ]; then echo "[$(date -Is)] runs/$c/LOWDISK present: stopping driver"; exit 2; fi
  for attempt in 1 2; do
    echo "[$(date -Is)] $c attempt $attempt start ($(ls "$dir" | wc -l) MAGs)"
    if bin/run_dram2.sh "$c" "$dir"; then
      echo "[$(date -Is)] $c DONE"
      if [ "$CLEAN_WORK" = 1 ] && [ -e "runs/$c/DONE" ]; then
        rm -rf "${WORK_ROOT:?}/$c" && echo "[$(date -Is)] $c work dir removed"
      fi
      break
    fi
    echo "[$(date -Is)] $c attempt $attempt FAILED: $(cat "runs/$c/FAILED" 2>/dev/null)"
    if [ -e "runs/$c/LOWDISK" ]; then echo "[$(date -Is)] low disk on $WORK_ROOT: stopping driver (rerun after freeing space and deleting runs/$c/LOWDISK)"; exit 2; fi
    [ $attempt = 2 ] && failed+=("$c")
  done
done

if [ ${#failed[@]} -eq 0 ] && [ "$(ls -d runs/chunk_*/DONE 2>/dev/null | wc -l)" = "$(ls -d inputs/chunk_*/ | wc -l)" ]; then
  date -Is > runs/ALL_CHUNKS_DONE
  echo "[$(date -Is)] ALL CHUNKS DONE -> next: bin/merge_and_distill.sh"
else
  echo "[$(date -Is)] finished with failed chunks: ${failed[*]:-none}; rerun this driver to resume them"
  exit 1
fi
