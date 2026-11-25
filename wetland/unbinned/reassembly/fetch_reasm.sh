#!/usr/bin/env bash
# Copy finished TF re-assemblies back (~/gtdbtk_kbase/reasm/<run>.min1000.fa), then run compare_reassembly.py and push the
# run's 1-2 kb contigs to poplar for priority Bakta. Polls every 5 min until all given runs are handled.
set -uo pipefail
L=$HOME/gtdbtk_kbase/reasm; mkdir -p "$L"; W=$HOME/Documents/MicrobiomeNotebooks/wetland/unbinned/reassembly
for r in "$@"; do
  until timeout 30 ssh -o BatchMode=yes TF "test -e tf_runs/wetland_reasm/$r/DONE" < /dev/null; do
    timeout 30 ssh -o BatchMode=yes TF "grep -q 'FAIL $r' tf_runs/wetland_reasm/events.log" < /dev/null && { echo "[$(date '+%F %T')] $r FAILED on TF"; continue 2; }
    sleep 300
  done
  [ -s "$L/$r.classes.tsv.gz" ] && continue
  rsync -a "TF:tf_runs/wetland_reasm/$r/$r.min1000.fa" "$L/" && echo "[$(date '+%F %T')] $r copied"
  ~/Documents/py_venv/bin/python -I "$W/compare_reassembly.py" "$r" && echo "[$(date '+%F %T')] $r compared" || { echo "$r compare FAILED"; continue; }
  rsync -a "$L/$r.1to2kb.fa" poplar:/scratch1/afreiburger/wetland_bakta_asm/reasm_inputs/ &&
    timeout 60 ssh -o BatchMode=yes poplar "cd /scratch1/afreiburger/wetland_bakta_asm && (setsid nohup bash bakta_priority.sh 1 reasm:/scratch1/afreiburger/wetland_bakta_asm/reasm_inputs/$r.1to2kb.fa >> logs/priority.log 2>&1 < /dev/null &)" < /dev/null &&
    echo "[$(date '+%F %T')] $r 1-2 kb Bakta launched"
done
echo "[$(date '+%F %T')] fetch_reasm finished"
