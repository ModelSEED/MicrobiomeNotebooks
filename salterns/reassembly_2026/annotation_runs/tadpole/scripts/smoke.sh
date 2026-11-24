#!/usr/bin/env bash
# End-to-end smoke test of run_sample.sh + eval_assemblies.py on 200,000 R1_A_D1 pairs (reads 100,001-300,000),
# with small threads/memory and k = 33,55, including Arm C. Run on poplar from tadpole_pilot/smoke.
set -euo pipefail
cd /scratch1/afreiburger/salterns_2026/tadpole_pilot/smoke
mkdir -p reads
for i in 1 2; do head -n 1200000 ../reads/SRR5209943_$i.fastq | tail -n 800000 > reads/SMOKE_$i.fastq; done
touch reads/SMOKE.ok
READS_DIR=$PWD/reads OUT_DIR=$PWD TPRE=8 TASM=4 XMX_BIG=8g XMX_SMALL=4g SPADES_MEM=32 KLIST=33,55 NOWAIT=1 \
  bash ../scripts/run_sample.sh R1_A_D1 SMOKE armC
/scratch1/afreiburger/salterns_2026/envs/spades3153/bin/python -I ../scripts/eval_assemblies.py \
  --root $PWD --targets ../scripts/targets.tsv --out-prefix $PWD/eval/smoke
