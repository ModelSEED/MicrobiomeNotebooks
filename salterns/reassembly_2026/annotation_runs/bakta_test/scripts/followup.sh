#!/bin/bash
# After PIPELINE_DONE: re-score the 35 Sec read-through grdH proteins (results/interim/sec_readthrough_grdH_joined.faa)
# against ALL 27,983 KOfam 2026-10 profiles (re-uses kofam/chunks), 16 x hmmsearch --cpu 1, so "K21579 is the best KO"
# is checked beyond the 68 marker profiles. Output: kofam/sec_joined.domtblout.gz, FOLLOWUP_DONE.
set -uo pipefail
W=/scratch1/afreiburger/salterns_2026/bakta_test
ANNOT=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_
export PATH=$ANNOT/bin:$PATH
cd $W
while [ ! -e PIPELINE_DONE ]; do sleep 120; done
mkdir -p kofam/sec_out
ls kofam/chunks/chunk_*.hmm | xargs -P 16 -I{} bash -c 'c=$(basename {} .hmm); nice -n 10 hmmsearch --cpu 1 --domtblout kofam/sec_out/$c.domtblout -o /dev/null {} grdscan/sec_joined.faa'
cat kofam/sec_out/chunk_*.domtblout | grep -v '^#' | gzip > kofam/sec_joined.domtblout.gz
touch FOLLOWUP_DONE
echo "FOLLOWUP_DONE $(date +%FT%T)"
