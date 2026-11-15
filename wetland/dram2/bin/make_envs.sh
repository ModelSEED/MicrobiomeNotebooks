#!/usr/bin/env bash
# Build the DRAM2 v2.0.0-beta17 module conda envs that cannot be reused from chenry's cache:
#   call / collect_rna / rename : chenry's copies (nxf_conda/env-8393b9af..., env-21ab46b9...,
#       env-7a8564ac...) were built on popjup under /storage/chenry/DRAM2, a path that does not
#       exist on poplar (quast.py shebang and tRNAscan-SE.conf bin_dir point there).
#   distill / product           : never built in chenry's tree.
# Each env is built from the exact environment.yml of the pinned beta17 module (copied to
# envs_yml/), with the user's own conda, into conda_cache/ under this run tree.
# Usage: make_envs.sh <name> [<name> ...]   (names: call collect_rna rename distill product)
set -uo pipefail
B=/scratch1/afreiburger/wetland_dram2
CONDA=/scratch/afreiburger/metag-hybrid/env/miniforge/condabin/conda
for y in "$@"; do
  p=$B/conda_cache/dram2_beta17_$y
  if [ -e "$p/conda-meta/history" ]; then echo "[$(date)] $y exists, skip"; continue; fi
  echo "[$(date)] creating $y -> $p"
  timeout 3h "$CONDA" env create -q -p "$p" -f "$B/envs_yml/$y.environment.yml"
  rc=$?
  echo "[$(date)] $y rc=$rc"
  if [ $rc -eq 0 ]; then
    "$CONDA" list -p "$p" --explicit --md5 > "$B/envs_yml/$y.explicit.txt"
    "$CONDA" list -p "$p" > "$B/envs_yml/$y.list.txt"
  fi
done
echo "[$(date)] ALL DONE"
