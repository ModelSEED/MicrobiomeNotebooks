#!/usr/bin/env bash
# Build the binning-test conda envs under /scratch1/afreiburger/salterns_2026/envs/ (new envs only) and fetch the
# CheckM1 data and the GUNC database. Run once on poplar, detached; each env has its own log; markers in envs/.bt_setup/.
# Usage: nohup setsid bash setup_envs.sh > /scratch1/afreiburger/salterns_2026/binning_test/logs/setup_envs.log 2>&1 < /dev/null &
set -uo pipefail
E=/scratch1/afreiburger/salterns_2026/envs
REF=/scratch1/afreiburger/salterns_2026/binning_test/refdata
M=/scratch1/afreiburger/metag-hybrid/env/miniforge/bin/mamba
export CONDA_PKGS_DIRS=$E/pkgs
mkdir -p $E/.bt_setup $REF
ts() { date +%F' '%T; }
mk() {  # mk <name> <specs...>
  local name=$1; shift
  if [ -e $E/.bt_setup/$name.ok ]; then echo "[$(ts)] $name exists"; return 0; fi
  echo "[$(ts)] create $name: $*"
  if nice -n 10 $M create -y -p $E/bt_$name --override-channels -c conda-forge -c bioconda "$@" > $E/bt_$name.create.log 2>&1; then
    $M list -p $E/bt_$name > $E/bt_$name.list.txt 2>&1; touch $E/.bt_setup/$name.ok; echo "[$(ts)] $name ok"
  else
    echo "[$(ts)] $name FAILED (see $E/bt_$name.create.log)"; tail -20 $E/bt_$name.create.log; touch $E/.bt_setup/$name.FAIL
  fi
}
mk map python=3.12 bowtie2 samtools pigz seqkit &
mk metabat2 metabat2=2.12.1 &
mk maxbin2 maxbin2=2.2.7 &
mk concoct concoct=1.1.0 samtools &
mk dastool das_tool=1.1.7 diamond=2.2.8 prodigal=2.6.3 &
mk gunc gunc=1.1.1 &
wait
# CheckM1 reference data (checkm_data_2015_01_16, the data set CheckM 1.x uses; KBase kb_Msuite uses the same)
if [ ! -e $REF/checkm1/.ok ]; then
  mkdir -p $REF/checkm1 && cd $REF/checkm1
  for u in https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz \
           https://zenodo.org/records/7401545/files/checkm_data_2015_01_16.tar.gz; do
    echo "[$(ts)] checkm1 data from $u"
    nice -n 10 curl -sSfL --retry 3 -C - -o checkm_data_2015_01_16.tar.gz "$u" && break
  done
  md5sum checkm_data_2015_01_16.tar.gz
  nice -n 10 tar xzf checkm_data_2015_01_16.tar.gz && rm -f checkm_data_2015_01_16.tar.gz && touch .ok && echo "[$(ts)] checkm1 data ok"
fi
# GUNC database (default progenomes 2.1 diamond db) via gunc download_db
if [ -e $E/.bt_setup/gunc.ok ] && [ ! -e $REF/gunc/.ok ]; then
  mkdir -p $REF/gunc
  echo "[$(ts)] gunc download_db"
  nice -n 10 $E/bt_gunc/bin/gunc download_db $REF/gunc > $REF/gunc/download.log 2>&1 && touch $REF/gunc/.ok && echo "[$(ts)] gunc db ok"
  ls -la $REF/gunc
fi
echo "[$(ts)] setup finished"
