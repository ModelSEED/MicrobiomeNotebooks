#!/bin/bash
# D3: Table S6 skani search, R232, with the Table S6 settings of reprocess_291/run_skani.sh
# (skani search --min-af 0 -n 3), run locally against the GTDB-Tk R232 skani database.
set -uo pipefail
cd "$(dirname "$0")"
DB=/home/freiburger/Documents/gtdbtk_data/release232/skani/database
SKANI=/home/freiburger/Documents/gtdbtk_venv/bin/skani
echo "[skani] start $(date)  db=$DB  $($SKANI --version)"
/usr/bin/time -v systemd-run --user --scope --unit=emilykin-skani-r232-s6 -p MemoryMax=60G -p MemorySwapMax=0 \
  nice -n 10 $SKANI search -d "$DB" -q /home/freiburger/Documents/EmilyKin/dereplicated_genomes/*.fa \
  -o mag_vs_gtdb_r232_skani.tsv -t 24 --min-af 0 -n 3
rc=$?
echo "[skani] rows: $(wc -l < mag_vs_gtdb_r232_skani.tsv 2>/dev/null)"
echo "SKANI_DONE rc=$rc $(date)"
