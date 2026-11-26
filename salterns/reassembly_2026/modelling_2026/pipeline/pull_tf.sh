#!/bin/bash
# Re-runnable: copy the TF build outputs back (md5-checked) into ../annotation_comparison/{runs,tf_logs,tf_env}.
set -euo pipefail
M26=/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/modelling_2026
RUN=tf_runs/salterns_modelling
AC=$M26/annotation_comparison
mkdir -p $AC/runs $AC/tf_logs $AC/tf_env
ssh TF "cd ~/$RUN && find out logs -type f -print0 | sort -z | xargs -0 md5sum > pulled.md5"
rsync -a TF:$RUN/out/ $AC/runs/
rsync -a TF:$RUN/logs/ $AC/tf_logs/
rsync -a TF:$RUN/src/COMMITS.tsv TF:$RUN/src/pip_freeze.txt TF:$RUN/gurobi.env TF:$RUN/inputs.md5 TF:$RUN/pulled.md5 TF:$RUN/jobs.tsv $AC/tf_env/
( cd $AC && sed -e 's#  out/#  runs/#' -e 's#  logs/#  tf_logs/#' tf_env/pulled.md5 | md5sum -c --quiet ) && echo "PULL_MD5_OK $(wc -l < $AC/tf_env/pulled.md5) files"
