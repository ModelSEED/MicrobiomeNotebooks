#!/usr/bin/env bash
# Salterns rollout: GTDB-Tk 2.8.0 env + GTDB R232 package (and the 26 MB R232 mockup DB used by the smoke test).
# Run once on poplar, detached:
#   nohup setsid nice -n 10 bash /scratch1/afreiburger/salterns_2026/rollout/scripts/setup_gtdbtk.sh \
#       >> /scratch1/afreiburger/salterns_2026/rollout/logs/setup_gtdbtk.log 2>&1 < /dev/null &
# Markers: envs/.ro_setup/gtdbtk280.ok, refdata/gtdbtk_r232_mockup/.ok, refdata/gtdbtk_r232/.ok (holds the data path).
# Integrity: every split part and the whole tarball are checked against GTDB's MD5SUM.txt
# (full package 25a59e0352b1fd150c589f56559767d4 for release232/232.0), then `gtdbtk check_install`.
set -uo pipefail
S=/scratch1/afreiburger/salterns_2026
E=$S/envs; ENV=$E/gtdbtk280; R=$S/refdata/gtdbtk_r232; MK=$S/refdata/gtdbtk_r232_mockup
M=/scratch1/afreiburger/metag-hybrid/env/miniforge/bin/mamba
PIGZ=$E/bt_map/bin/pigz
BASE=https://data.gtdb.ecogenomic.org/releases/release232/232.0
export CONDA_PKGS_DIRS=$E/pkgs
mkdir -p $E/.ro_setup $R/parts $MK
ts() { date +%F' '%T; }
echo "[$(ts)] setup_gtdbtk start (pid $$)"

env_build() {
  if [ -e $E/.ro_setup/gtdbtk280.ok ]; then echo "[$(ts)] env exists"; return 0; fi
  echo "[$(ts)] create env gtdbtk=2.8.0 at $ENV"
  rm -rf $ENV
  if nice -n 10 $M create -y -p $ENV --override-channels -c conda-forge -c bioconda gtdbtk=2.8.0 > $E/gtdbtk280.create.log 2>&1 \
     && $ENV/bin/gtdbtk --version > $E/gtdbtk280.version.txt 2>&1 && grep -q "2.8.0" $E/gtdbtk280.version.txt; then
    $M list -p $ENV > $E/gtdbtk280.list.txt 2>&1
    touch $E/.ro_setup/gtdbtk280.ok; echo "[$(ts)] env ok: $(head -3 $E/gtdbtk280.version.txt | tr '\n' ' ')"
  else
    echo "[$(ts)] env FAILED (see $E/gtdbtk280.create.log)"; tail -20 $E/gtdbtk280.create.log; touch $E/.ro_setup/gtdbtk280.FAIL
  fi
}

fetch_md5sum() {
  for k in 1 2 3 4 5; do
    curl -sSf --retry 3 -o $R/MD5SUM.txt.tmp $BASE/MD5SUM.txt && mv $R/MD5SUM.txt.tmp $R/MD5SUM.txt && return 0
    sleep 30
  done
  return 1
}

getpart() {  # getpart <url> <out> <md5>: resume until the md5 matches (at most 8 tries)
  local u=$1 o=$2 want=$3 k got
  for k in 1 2 3 4 5 6 7 8; do
    if [ -s "$o" ]; then got=$(md5sum "$o" | cut -d' ' -f1); [ "$got" = "$want" ] && { echo "[$(ts)] ok $(basename $o)"; return 0; }; fi
    [ "$k" -gt 1 ] && { echo "[$(ts)] $(basename $o) md5 ${got:-none} != $want, try $k"; [ "$k" -gt 3 ] && rm -f "$o"; }
    nice -n 10 curl -sS --retry 5 --retry-delay 20 -C - -o "$o" "$u" || sleep 20
  done
  echo "[$(ts)] FAILED $(basename $o)"; return 1
}

mockup() {
  [ -e $MK/.ok ] && { echo "[$(ts)] mockup exists"; return 0; }
  local m1 m2
  m1=$(awk '$2=="./auxillary_files/gtdbtk_package/mockup_db/r232_mockup.tar.gz"{print $1}' $R/MD5SUM.txt)
  m2=$(awk '$2=="./auxillary_files/gtdbtk_package/mockup_db/test_genome_G021783475.fna.gz"{print $1}' $R/MD5SUM.txt)
  cd $MK
  getpart $BASE/auxillary_files/gtdbtk_package/mockup_db/r232_mockup.tar.gz r232_mockup.tar.gz $m1 || return 1
  getpart $BASE/auxillary_files/gtdbtk_package/mockup_db/test_genome_G021783475.fna.gz test_genome_G021783475.fna.gz $m2 || return 1
  tar xzf r232_mockup.tar.gz && d=$(dirname "$(find $MK -maxdepth 3 -type d -name metadata | head -1)") && echo "$d" > DATA_PATH && touch .ok \
    && echo "[$(ts)] mockup ok: $(cat DATA_PATH)"
}

full() {
  [ -e $R/.ok ] && { echo "[$(ts)] R232 exists: $(cat $R/.ok)"; return 0; }
  local full_md5 parts p
  full_md5=$(awk '$2=="./auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz"{print $1}' $R/MD5SUM.txt)
  [ "$full_md5" = 25a59e0352b1fd150c589f56559767d4 ] || { echo "unexpected full md5 '$full_md5' in MD5SUM.txt"; return 1; }
  grep "split_package/gtdbtk_r232_data.tar.gz.part_" $R/MD5SUM.txt | awk '{print $1, $2}' > $R/parts.md5
  echo "[$(ts)] $(wc -l < $R/parts.md5) split parts; full md5 $full_md5"
  if [ ! -e $R/.extracted ]; then
    # parts in parallel, 6 at a time
    local pids=() n=0
    while read -r md5 path; do
      p=$(basename $path)
      getpart "$BASE/${path#./}" "$R/parts/$p" "$md5" > $R/parts/$p.log 2>&1 & pids+=($!); n=$((n+1))
      if [ $((n % 6)) -eq 0 ]; then for x in "${pids[@]}"; do wait $x; done; pids=(); fi
    done < $R/parts.md5
    for x in "${pids[@]}"; do wait $x; done
    while read -r md5 path; do p=$(basename $path); tail -1 $R/parts/$p.log; grep -q "ok $p" $R/parts/$p.log || { echo "part $p failed"; return 1; }; done < $R/parts.md5
    echo "[$(ts)] all parts ok; checking the whole-tarball md5"
    got=$(cd $R/parts && cat $(awk '{n=split($2,a,"/"); print a[n]}' $R/parts.md5) | md5sum | cut -d' ' -f1)
    echo "[$(ts)] whole tarball md5 $got (want $full_md5)"
    [ "$got" = "$full_md5" ] || { echo "WHOLE MD5 MISMATCH"; return 1; }
    echo "$got  gtdbtk_r232_data.tar.gz (cat of $(wc -l < $R/parts.md5) split parts)" > $R/full_md5.verified
    echo "[$(ts)] extracting"
    rm -rf $R/extract.tmp; mkdir -p $R/extract.tmp
    (cd $R/parts && cat $(awk '{n=split($2,a,"/"); print a[n]}' $R/parts.md5)) | nice -n 10 $PIGZ -dc -p 4 | (cd $R/extract.tmp && nice -n 10 tar xf -) || { echo "extract FAILED"; return 1; }
    touch $R/.extracted; echo "[$(ts)] extracted"
  fi
  d=$(dirname "$(find $R/extract.tmp -maxdepth 3 -type d -name metadata | head -1)")
  [ -d "$d/markers" ] || { echo "no GTDB-Tk data dir found under $R/extract.tmp"; return 1; }
  echo "[$(ts)] data dir $d; waiting for the env"
  while [ ! -e $E/.ro_setup/gtdbtk280.ok ]; do [ -e $E/.ro_setup/gtdbtk280.FAIL ] && { echo "env failed"; return 1; }; sleep 60; done
  if GTDBTK_DATA_PATH=$d nice -n 10 $ENV/bin/gtdbtk check_install > $R/check_install.log 2>&1; then
    echo "$d" > $R/DATA_PATH; echo "$d" > $R/.ok; rm -rf $R/parts; echo "[$(ts)] R232 ok (check_install passed); parts deleted"
  else
    echo "[$(ts)] check_install FAILED"; tail -30 $R/check_install.log; return 1
  fi
}

env_build &
EPID=$!
fetch_md5sum || { echo "[$(ts)] cannot fetch MD5SUM.txt"; wait $EPID; exit 1; }
mockup || echo "[$(ts)] mockup FAILED"
full || { echo "[$(ts)] R232 FAILED"; touch $R/.FAIL; }
wait $EPID
echo "[$(ts)] setup_gtdbtk finished"
