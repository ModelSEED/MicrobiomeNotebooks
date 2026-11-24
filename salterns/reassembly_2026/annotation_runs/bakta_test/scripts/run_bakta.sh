#!/bin/bash
# Bakta 1.12.0 (full DB v6.0, 2025-02-24), default options, --threads 8 --skip-plot, one MAG per job.
# Usage (inside /scratch1/afreiburger/salterns_2026/bakta_test): bash scripts/run_bakta.sh <list of labels> <parallel jobs>
# Each job is timed with /usr/bin/time -v -> logs/bakta_<label>.time
set -uo pipefail
W=/scratch1/afreiburger/salterns_2026/bakta_test
ANNOT=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_
DB=/scratch1/afreiburger/emilykin/processed/refdata/bakta/db
export PATH=$ANNOT/bin:$PATH W DB
cd $W; mkdir -p bakta logs tmp
one(){
  l=$1
  if [ -s bakta/$l/$l.json ]; then echo "SKIP $l"; return; fi
  rm -rf bakta/$l; mkdir -p tmp/bakta_$l
  echo "START $l $(date +%FT%T)"
  /usr/bin/time -v -o logs/bakta_$l.time nice -n 10 bakta --db $DB --threads 8 --skip-plot \
      --tmp-dir tmp/bakta_$l --output bakta/$l --prefix $l fasta/$l.fna > logs/bakta_$l.log 2>&1 \
    && echo "OK $l $(date +%FT%T)" || echo "FAIL $l $(date +%FT%T)"
  rm -rf tmp/bakta_$l
}
export -f one
xargs -P ${2:-4} -I{} bash -c 'one {}' < $1
echo "BAKTA_DONE $(date +%FT%T)"
