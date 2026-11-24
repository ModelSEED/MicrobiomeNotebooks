#!/bin/bash
# Re-runnable TF driver: builds every (mag, guild, arm) job in jobs.tsv not yet reported 'ok', N in parallel,
# nice 10, one Gurobi thread each (gurobi.env). Run detached:  nohup setsid bash code/run_tf.sh jobs.tsv 18 &
set -uo pipefail
cd ~/tf_runs/salterns_modelling
JOBS=${1:-jobs.tsv}; N=${2:-18}
export MODELSEED_DB_PATH=$HOME/tf_runs/salterns_modelling/src/ModelSEEDDatabase
PY=$HOME/venvs/salterns_modelling/bin/python
run_one() {
  mag=$1; guild=$2; arm=$3; tag=${mag}__${guild}__${arm}
  rep=out/$arm/$tag.report.json
  if [ -f "$rep" ] && grep -q '"status": "ok"' "$rep"; then echo "skip $tag"; return 0; fi
  /usr/bin/time -v -o logs/$tag.time nice -n 10 $PY code/build_model.py --root inputs --mag $mag --guild $guild \
     --arm $arm --out out/$arm --msdb $MODELSEED_DB_PATH > logs/$tag.log 2>&1
  echo "$(date -u +%FT%TZ) done $tag exit $?"
}
export -f run_one; export PY MODELSEED_DB_PATH
echo "$(date -u +%FT%TZ) start $(wc -l < $JOBS) jobs, $N parallel" >> logs/driver.log
tr '\t' ' ' < $JOBS | xargs -P $N -L 1 bash -c 'run_one $0 $1 $2' >> logs/driver.log 2>&1
echo "$(date -u +%FT%TZ) ALL_DONE" >> logs/driver.log
