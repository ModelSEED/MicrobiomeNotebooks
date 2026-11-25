#!/usr/bin/env bash
# DRAM2 for the run-4 representative set (1,206): merge all 6 chunks, subset to the set, distill + product on it.
set -uo pipefail
B=/scratch1/afreiburger/wetland_dram2; cd $B
PY=$B/conda_cache/dram2_beta17_distill/bin/python
ev() { echo "[$(date '+%F %T')] $*" >> $B/logs/finalize1206.events; }
ev START
$PY -I bin/merge_dram2.py >> logs/finalize1206.log 2>&1 || { ev "FAIL merge"; exit 1; }
ev "merged: $(tail -n +2 merged/raw-annotations.all.tsv | wc -l) genes"
$PY -I bin/subset_merged.py dram2_hq_multisample_run4.tsv merged_set1206 >> logs/finalize1206.log 2>&1 || { ev "FAIL subset"; exit 1; }
ev "subset: $(tail -1 logs/finalize1206.log)"
RUN=all_mags_v2 SKIP_MERGE=1 ANN=$B/merged_set1206/raw-annotations.tsv RRNAS=$B/merged_set1206/rrnas TRNAS=$B/merged_set1206/trnas \
  bash bin/merge_and_distill.sh >> logs/finalize1206.log 2>&1 && ev "DONE distill+product (runs/all_mags_v2)" || ev "FAIL distill"
