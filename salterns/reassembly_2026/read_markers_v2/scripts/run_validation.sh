#!/usr/bin/env bash
# Validation driver (v2): v1 design (07_simulate.py hold-out of 20 % of each label's 90 % clusters, validation DB
# without them, simulated 148-nt reads with 0.3 % errors, background genes, KOfam-twilight homologs).
# DIAMOND 2.2.8 blastx of the simulated reads against the validation DB with the production flags
# (--mid-sensitive -k 5 -e 1e-3 -b 4 -c 1) and, as a test of top-k saturation in the larger v2 DB, with -k 25;
# blastp of the held-out proteins against the validation DB (divergence strata); the cut-off grid (08 validate,
# reduced to EMAX 1e-5, LMIN 30, MARGIN 0/2, PMIN 0/50/60/70: v1 showed EMAX and LMIN insensitive).
set -euo pipefail
W=/scratch1/afreiburger/salterns_2026/read_markers_v2
DM=/scratch1/afreiburger/salterns_2026/envs/bt_dastool/bin/diamond
A=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin
export TMPDIR=$W/tmp
T=${T:-32}
cd $W; mkdir -p val
[ -f val/sim_reads.fa ] || nice -n 10 $A/python3 -I scripts/07_simulate.py > logs/07_simulate.log 2>&1
nice -n 10 $DM makedb --in val/val_survey.faa -d val/val_survey --threads $T --quiet
for k in 5 25; do
  [ -f val/sim_mid-sensitive_k$k.tsv.gz ] || nice -n 10 $DM blastx -d val/val_survey -q val/sim_reads.fa -o val/sim_mid-sensitive_k$k.tsv.gz \
    --compress 1 --mid-sensitive -k $k -e 1e-3 --threads $T -b 4 -c 1 -t $TMPDIR --quiet \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
done
ln -sf sim_mid-sensitive_k5.tsv.gz val/sim_mid-sensitive.tsv.gz
cut -f1 val/holdout_genes.tsv | tail -n +2 > val/holdout.ids
$A/seqkit grep -f val/holdout.ids refs/labeled.faa > val/holdout.faa
nice -n 10 $DM blastp -d val/val_survey -q val/holdout.faa -o val/holdout_vs_valdb.tsv -k 200 --sensitive --threads $T --quiet \
  --outfmt 6 qseqid sseqid pident length evalue bitscore
for k in 5 25; do
  nice -n 10 $A/python3 -I scripts/08_classify.py validate val/sim_mid-sensitive_k$k.tsv.gz val/sim_truth.tsv val/grid_mid-sensitive_k$k.tsv &
done
wait
cp val/grid_mid-sensitive_k5.tsv val/grid_mid-sensitive.tsv
echo "$(date) validation done" > val/VAL_DONE
