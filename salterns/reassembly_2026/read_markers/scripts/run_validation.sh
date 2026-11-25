#!/usr/bin/env bash
# Validation driver: hold-out simulation (07), DIAMOND blastx of the simulated reads against the validation DB in three
# sensitivity modes, blastp of the held-out proteins against the validation DB (divergence strata), the cut-off grid (08).
set -euo pipefail
W=/scratch1/afreiburger/salterns_2026/read_markers
DM=/scratch1/afreiburger/salterns_2026/envs/bt_dastool/bin/diamond
A=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin
export TMPDIR=$W/tmp
T=${T:-16}
cd $W
[ -f val/sim_reads.fa ] || nice -n 10 $A/python3 -I scripts/07_simulate.py > logs/07_simulate.log 2>&1
nice -n 10 $DM makedb --in val/val_survey.faa -d val/val_survey --threads $T --quiet
for mode in default mid-sensitive sensitive; do
  flag=""; [ $mode != default ] && flag="--$mode"
  [ -f val/sim_$mode.tsv.gz ] || nice -n 10 $DM blastx -d val/val_survey -q val/sim_reads.fa -o val/sim_$mode.tsv.gz --compress 1 \
    $flag -k 5 -e 1e-3 --threads $T -b 4 -c 1 -t $TMPDIR --quiet \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
done
cut -f1 val/holdout_genes.tsv | tail -n +2 > val/holdout.ids
$A/seqkit grep -f val/holdout.ids refs/labeled.faa > val/holdout.faa
nice -n 10 $DM blastp -d val/val_survey -q val/holdout.faa -o val/holdout_vs_valdb.tsv -k 200 --sensitive --threads $T --quiet \
  --outfmt 6 qseqid sseqid pident length evalue bitscore
for mode in default mid-sensitive sensitive; do
  nice -n 10 $A/python3 -I scripts/08_classify.py validate val/sim_$mode.tsv.gz val/sim_truth.tsv val/grid_$mode.tsv &
done
wait
echo "$(date -u +%FT%TZ) validation done" > val/VAL_DONE
