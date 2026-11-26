#!/usr/bin/env bash
# Summaries after all 39 samples are classified and re-mapped: 10 (per_sample_KO, GE, guilds), 10b (placement, read
# taxonomy, MAG attribution, assembly complement), 10c (key tables, mcrA / C-P lyase taxonomy, split per gene,
# methanogenesis-pathway archaeal share), 10d (methanol taxonomy), 12 (v1 -> v2 regression) and the lost-read
# attribution for three samples. 1 CPU each, sequential, nice 10.
set -euo pipefail
W=/scratch1/afreiburger/salterns_2026/read_markers_v2
PY=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin/python3
cd $W
n=$(ls logs/post/*.done | wc -l); [ "$n" = 39 ] || { echo "only $n of 39 samples post-processed" >&2; exit 1; }
for s in 10_summarize.py 10b_split_and_taxonomy.py 10c_key_tables.py 10d_methanol_taxonomy.py 12_regression_v1_v2.py; do
  echo "$(date) $s start"; nice -n 10 $PY -I scripts/$s > logs/$s.log 2>&1; echo "$(date) $s done"
done
nice -n 10 $PY -I scripts/12_regression_v1_v2.py attribute R1_A_D1 R2_A_D1 SF2_C_H2O > logs/12_attribute.log 2>&1
echo "$(date) summaries done" > results/SUMMARIES_DONE
