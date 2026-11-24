#!/bin/bash
# Re-runnable TF driver for the final-set repairs of the winning arm (default A): one repair_model.py per model,
# N in parallel, nice 10, Gurobi 1 thread (gurobi.env). Skips models whose repair report says ok.
set -uo pipefail
cd ~/tf_runs/salterns_modelling
ARM=${1:-A}; N=${2:-12}; LIST=${3:-}
export MODELSEED_DB_PATH=$HOME/tf_runs/salterns_modelling/src/ModelSEEDDatabase
PY=$HOME/venvs/salterns_modelling/bin/python
mkdir -p final/$ARM logs_final
one() {
  tag=$1; rep=final/$ARM/$tag.report.json
  if [ -f "$rep" ] && grep -q '"status": "ok"' "$rep"; then echo "skip $tag"; return 0; fi
  /usr/bin/time -v -o logs_final/$tag.time nice -n 10 $PY code/repair_model.py --root inputs --runs out/$ARM --tag $tag \
     --out final/$ARM --msdb $MODELSEED_DB_PATH > logs_final/$tag.log 2>&1
  rc=$?; echo "$(date -u +%FT%TZ) done $tag exit $rc"
}
export -f one; export PY MODELSEED_DB_PATH ARM
if [ -n "$LIST" ]; then cat $LIST; else ls out/$ARM/*.report.json | xargs -n1 basename | sed 's/.report.json$//'; fi \
  | xargs -P $N -I{} bash -c 'one {}' >> logs_final/driver.log 2>&1
echo "$(date -u +%FT%TZ) ALL_DONE" >> logs_final/driver.log
