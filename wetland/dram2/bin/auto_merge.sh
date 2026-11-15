#!/usr/bin/env bash
# Wait for all 5 DRAM2 chunks, then run the tested merge (merged annotations + one distill/product over all 1,184 MAGs).
cd /scratch1/afreiburger/wetland_dram2
until [ -e runs/ALL_CHUNKS_DONE ] || ls runs/chunk_0*/FAILED >/dev/null 2>&1; do sleep 900; done
if ls runs/chunk_0*/FAILED >/dev/null 2>&1; then echo "[$(date -Is)] a chunk FAILED; merge not started" >> logs/auto_merge.log; exit 1; fi
echo "[$(date -Is)] all chunks done; starting merge" >> logs/auto_merge.log
bin/merge_and_distill.sh > logs/merge_and_distill.log 2>&1
echo "[$(date -Is)] merge exit $?" >> logs/auto_merge.log
