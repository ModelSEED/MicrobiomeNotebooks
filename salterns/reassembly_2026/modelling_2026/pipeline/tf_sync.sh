#!/bin/bash
# Re-runnable: push code + inputs to TF (~/tf_runs/salterns_modelling), checked by md5.
set -euo pipefail
M26=/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/modelling_2026
RUN=tf_runs/salterns_modelling
ssh TF "mkdir -p ~/$RUN/inputs/templates ~/$RUN/inputs/evidence ~/$RUN/code ~/$RUN/out ~/$RUN/logs"
rsync -a $M26/templates/{GramNegModelTemplateV7,GramPosModelTemplateV7,ArchaeaTemplateV6,Core-V6}.json TF:$RUN/inputs/templates/
rsync -a $M26/annotation_comparison/evidence/ TF:$RUN/inputs/evidence/
rsync -a $M26/../comparison_rollout/modelling_species.tsv TF:$RUN/inputs/
rsync -a $M26/pipeline/build_model.py $M26/pipeline/run_tf.sh $M26/pipeline/repair_model.py $M26/pipeline/run_repair_tf.sh TF:$RUN/code/
printf "Threads 1\nTimeLimit 1200\n" | ssh TF "cat > ~/$RUN/gurobi.env"
( cd $M26/templates && md5sum {GramNegModelTemplateV7,GramPosModelTemplateV7,ArchaeaTemplateV6,Core-V6}.json | sed 's#  #  inputs/templates/#'
  cd $M26/annotation_comparison/evidence && md5sum *.json.gz | sed 's#  #  inputs/evidence/#'
  cd $M26/../comparison_rollout && md5sum modelling_species.tsv | sed 's#  #  inputs/#'
  cd $M26/pipeline && md5sum build_model.py run_tf.sh repair_model.py run_repair_tf.sh | sed 's#  #  code/#' ) > $M26/pipeline/tf_inputs.md5
scp -q $M26/pipeline/tf_inputs.md5 TF:$RUN/inputs.md5
ssh TF "cd ~/$RUN && md5sum -c --quiet inputs.md5 && echo TF_INPUTS_MD5_OK \$(wc -l < inputs.md5) files"
