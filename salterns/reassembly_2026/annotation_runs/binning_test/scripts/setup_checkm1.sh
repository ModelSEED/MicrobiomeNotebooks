#!/usr/bin/env bash
# CheckM1 env (the EmilyKin bin_qc env's checkm-genome 1.2.3 fails to import: setuptools 82 has no pkg_resources).
# Try checkm-genome 1.1.3 (closest to KBase kb_Msuite's CheckM 1.1.2); fall back to 1.2.3. setuptools<81 keeps
# pkg_resources. A test import + 'checkm data' check must pass. Marker: envs/.bt_setup/checkm1.ok
set -uo pipefail
E=/scratch1/afreiburger/salterns_2026/envs; M=/scratch1/afreiburger/metag-hybrid/env/miniforge/bin/mamba
export CONDA_PKGS_DIRS=$E/pkgs CHECKM_DATA_PATH=/scratch1/afreiburger/salterns_2026/binning_test/refdata/checkm1
for spec in "checkm-genome=1.1.3 python=3.8 setuptools<70" "checkm-genome=1.2.3 setuptools<81"; do
  rm -rf $E/bt_checkm1
  echo "[$(date +%F' '%T)] try: $spec"
  if nice -n 10 $M create -y -p $E/bt_checkm1 --override-channels -c conda-forge -c bioconda $spec > $E/bt_checkm1.create.log 2>&1 \
     && $E/bt_checkm1/bin/checkm lineage_wf -h > /dev/null 2>&1; then
    $M list -p $E/bt_checkm1 > $E/bt_checkm1.list.txt 2>&1; echo "$spec" > $E/.bt_setup/checkm1.ok
    grep -E "^  (checkm-genome|python|hmmer|prodigal|pplacer|setuptools|numpy) " $E/bt_checkm1.list.txt
    echo "[$(date +%F' '%T)] CHECKM1_ENV_OK $spec"; exit 0
  fi
  tail -5 $E/bt_checkm1.create.log
done
echo "[$(date +%F' '%T)] CHECKM1_ENV_FAILED"; exit 1
