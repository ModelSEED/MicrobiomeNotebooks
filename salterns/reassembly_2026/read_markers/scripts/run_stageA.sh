#!/usr/bin/env bash
# Waits for 01_prodigal.sh and 01b_prodigal_coasm.sh, then runs the panel hmmsearch on the 39 assemblies + 520 co-assembly MAGs.
set -euo pipefail
W=/scratch1/afreiburger/salterns_2026/read_markers
cd $W
until [ -f genes/N_DONE ] && [ -f genes/coasm520/DONE ]; do sleep 60; done
n=$(for s in $(ls /scratch1/afreiburger/salterns_2026/rollout/reads); do [ -f genes/$s/DONE ] && echo $s; done | wc -l); [ "$n" = 39 ] || { echo "only $n of 39 samples have genes" >&2; exit 1; }
inputs=$(for s in $(ls /scratch1/afreiburger/salterns_2026/rollout/reads); do echo genes/$s/$s.faa; done; echo genes/coasm520/coasm520.faa)
NP=${NP:-24} bash scripts/03_panel_hmmsearch.sh $inputs
n=$(ls hmm/*.done | wc -l)
echo "$(date -u +%FT%TZ) stage A done, $n hmm inputs" > hmm/STAGE_A_DONE
