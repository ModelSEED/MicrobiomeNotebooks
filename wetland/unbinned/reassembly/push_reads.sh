#!/usr/bin/env bash
# Download ENA reads for the re-assembly test on the workstation, push them to TF and drop READS_OK there for
# run_reassembly_tf.sh (DL=0). ENA serves ~1 MB/s per stream but honours HTTP ranges, so each file is fetched as
# NPART parallel ranges and re-joined; md5 checked against ENA here and again on TF. The local copy is deleted after.
# Usage: setsid nohup bash push_reads.sh SRR... > ~/gtdbtk_kbase/reasm_push.log 2>&1 < /dev/null &
set -uo pipefail
L=$HOME/gtdbtk_kbase/reasm_reads; R=tf_runs/wetland_reasm; NPART=${NPART:-16}
part() {   # url partfile start end: fetch bytes start..end into partfile, resuming
  local url=$1 p=$2 s=$3 e=$4 want have t
  want=$(( e - s + 1 ))
  for t in 1 2 3 4 5 6; do
    have=$(stat -c %s "$p" 2>/dev/null || echo 0)
    [ "$have" -eq "$want" ] && return 0
    [ "$have" -gt "$want" ] && { rm -f "$p"; have=0; }
    curl -s -L --max-time 7200 -r $(( s + have ))-"$e" "$url" >> "$p"; sleep $(( t * 10 ))
  done
  [ "$(stat -c %s "$p" 2>/dev/null || echo 0)" -eq "$want" ]
}
export -f part
pget() {   # url size md5 out
  local url=$1 size=$2 md5=$3 out=$4 i
  mkdir -p "$out.parts"
  for ((i = 0; i < NPART; i++)); do
    echo "$url $out.parts/$i $(( size * i / NPART )) $(( size * (i + 1) / NPART - 1 ))"
  done | xargs -P "$NPART" -n 4 bash -c 'part "$@"' _
  for ((i = 0; i < NPART; i++)); do cat "$out.parts/$i"; done > "$out"
  if [ "$(stat -c %s "$out")" -eq "$size" ] && [ "$(md5sum "$out" | cut -d' ' -f1)" = "$md5" ]; then rm -rf "$out.parts"; return 0; fi
  rm -f "$out"; return 1
}
for r in "$@"; do
  mkdir -p "$L/$r"; cd "$L/$r" || exit 1
  ssh -o BatchMode=yes TF "test -e $R/$r/READS_OK" < /dev/null && { echo "$r already on TF"; continue; }
  rep=$(curl -s --retry 5 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$r&result=read_run&fields=fastq_ftp,fastq_bytes,fastq_md5&format=tsv" | tail -1)
  IFS=';' read -r -a urls <<< "$(echo "$rep" | cut -f2)"; IFS=';' read -r -a sizes <<< "$(echo "$rep" | cut -f3)"
  IFS=';' read -r -a md5s <<< "$(echo "$rep" | cut -f4)"
  [ "${#urls[@]}" -eq 2 ] || { echo "FAIL $r: $rep"; continue; }
  : > md5.txt; ok=1; t0=$(date +%s)
  for i in 0 1; do
    f=$(basename "${urls[$i]}"); echo "${md5s[$i]}  $f" >> md5.txt
    [ -s "$f" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "${md5s[$i]}" ] && continue
    rm -f "$f"
    pget "https://${urls[$i]}" "${sizes[$i]}" "${md5s[$i]}" "$f" || pget "https://${urls[$i]}" "${sizes[$i]}" "${md5s[$i]}" "$f" || ok=0
  done
  [ $ok -eq 1 ] || { echo "[$(date '+%F %T')] FAIL $r download"; continue; }
  echo "[$(date '+%F %T')] $r downloaded in $(( $(date +%s) - t0 ))s"
  rsync -a --partial ./*.fastq.gz md5.txt "TF:$R/$r/" &&
    ssh -o BatchMode=yes TF "cd $R/$r && md5sum -c --quiet md5.txt && touch READS_OK" < /dev/null &&
    { echo "[$(date '+%F %T')] $r pushed"; rm -f "$L/$r"/*.fastq.gz; } || echo "FAIL $r push"
done
echo "[$(date '+%F %T')] ALL pushed"
