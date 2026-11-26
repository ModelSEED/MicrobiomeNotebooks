#!/bin/bash
# usage: run.sh <logname> <cmd...>   (cwd = repo root)
SC=/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan/apply_scratch/figures
name=$1; shift
cd /home/freiburger/Documents/EmilyKin
export OMP_NUM_THREADS=24 OPENBLAS_NUM_THREADS=24 MKL_NUM_THREADS=24 NUMEXPR_NUM_THREADS=24 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
{ echo "# $(date '+%F %T') cmd: $*"; } > $SC/logs/$name.log
/usr/bin/time -v nice -n 10 "$@" >> $SC/logs/$name.log 2> $SC/logs/$name.stderr; rc=$?
echo "# rc=$rc" >> $SC/logs/$name.log
grep -E "Elapsed|Maximum resident" $SC/logs/$name.stderr | sed 's/^\s*/# /' >> $SC/logs/$name.log
exit $rc
