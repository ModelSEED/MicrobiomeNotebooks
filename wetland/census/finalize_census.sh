#!/usr/bin/env bash
# Wait for the poplar census compute (CheckM2), copy outputs back, run census.py + select_reps.py with CheckM2.
# Usage: setsid nohup bash finalize_census.sh > ~/gtdbtk_kbase/census_final.log 2>&1 < /dev/null &
set -uo pipefail
W=/scratch1/afreiburger/wetland_census
L=$HOME/gtdbtk_kbase/census_compute
OUT=$HOME/gtdbtk_kbase/census_final
HERE=$(cd "$(dirname "$0")" && pwd)
PY="$HOME/Documents/py_venv/bin/python -I"
log() { echo "[$(date '+%F %T')] $*"; }
deadline=$(( $(date +%s) + 24*3600 ))
until timeout 60 ssh -o BatchMode=yes poplar "grep -q DONE $W/census.events" < /dev/null; do
  [ "$(date +%s)" -gt "$deadline" ] && { log "FAIL: compute not done within 24 h"; exit 1; }
  sleep 300
done
log "poplar: $(timeout 60 ssh -o BatchMode=yes poplar "tail -1 $W/census.events" < /dev/null)"
mkdir -p "$L/checkm2"
scp -q poplar:$W/checkm2/quality_report.tsv "$L/checkm2/" && scp -q "poplar:$W/time_checkm2.txt" "poplar:$W/census.events" "$L/" || { log "FAIL: copy"; exit 1; }
n=$(( $(wc -l < "$L/checkm2/quality_report.tsv") - 1 ))
log "CheckM2 rows: $n"
[ "$n" -eq 12859 ] || { log "FAIL: expected 12859 CheckM2 rows"; exit 1; }
$PY "$HERE/census.py" "$L" "$OUT" > "$OUT.census.log" 2>&1 || { log "FAIL: census.py ($OUT.census.log)"; exit 1; }
$PY "$HERE/select_reps.py" "$OUT" > "$OUT.reps.log" 2>&1 || { log "FAIL: select_reps.py ($OUT.reps.log)"; exit 1; }
log "DONE census_final written: $OUT"
