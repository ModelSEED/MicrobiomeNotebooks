#!/bin/bash
# Runs validate_build.py for the three test MAGs on both template arms (6 processes, 2 Gurobi threads each).
cd "$(dirname "$0")"
G="/home/freiburger/Documents/MicrobiomeNotebooks/salterns/jplfaria/Cliff MAGs/RASTified_MAG_json"
PY=/home/freiburger/Documents/py_venv/bin/python
declare -A GEN=( [Roseovarius_R2BD2_mb38]="Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.38.contigs__.RAST.json gramneg"
                 [BM714_R1AD1_mb51]="Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs__.RAST.json grampos"
                 [Methanolobus_R1AD2_cc21]="Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_concoct_out.21.contigs__.RAST.json archaea" )
for lab in "${!GEN[@]}"; do
  set -- ${GEN[$lab]}
  for arm in combined upstream; do
    $PY validate_build.py --genome "$G/$1" --template-kind $2 --arm $arm --label ${lab}__$arm \
        --out results/${lab}__$arm.json > logs/${lab}__$arm.log 2>&1 &
  done
done
wait
