#!/usr/bin/env bash
# Threshold test (per Andrew 2026-11-25, "would lowering it to 1kb better capture the system?"): the KBase assemblies
# (narrative 261700) were made by the KBase MEGAHIT app with min_contig_length = 2000, so they hold no contig < 2 kb.
# This re-assembles a few samples with MEGAHIT 1.2.9 and the same parameters (MEGAHIT defaults: k-list 21..141,
# min-count 2) except --min-contig-len 1000, so the 1-2 kb contigs can be measured. Contigs >= 2 kb are then checked
# against the KBase assembly (compare_reassembly.py) to show the re-assembly reproduces it.
# Reads: ENA fastq (md5-checked), all samples downloaded in parallel. MEGAHIT one sample at a time (TF has ~23 GB RAM free); waits until MemAvailable > MEM + 4 GB.
# Resumable: <run>/DONE. Events: events.log.
# Usage (TF): cd ~/tf_runs/wetland_reasm && DL=0 setsid nohup bash run_reassembly_tf.sh SRR... > driver.log 2>&1 < /dev/null &
set -uo pipefail
B=$HOME/tf_runs/wetland_reasm; MH=$B/megahit/bin/megahit   # wrapper is "#!/usr/bin/env python"; TF has only python3
THREADS=${THREADS:-16}; MEM=${MEM:-18000000000}
DL=${DL:-1}   # DL=0: reads are pushed by push_reads.sh (workstation; ENA is ~6x faster there), which drops READS_OK
ev() { echo "[$(date '+%F %T')] $*" >> "$B/events.log"; }
exec 9> "$B/driver.lock"; flock -n 9 || { echo "another driver is running"; exit 1; }
dl() {   # download one sample's reads (md5-checked); marker <run>/READS_OK or READS_FAIL
  local r=$1 d=$B/$1 rep f i t ok=1; mkdir -p "$d"; cd "$d" || return
  [ -e READS_OK ] && return
  rep=$(curl -s --retry 5 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$r&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv" | tail -1)
  IFS=';' read -r -a urls <<< "$(echo "$rep" | cut -f2)"; IFS=';' read -r -a md5s <<< "$(echo "$rep" | cut -f3)"
  [ "${#urls[@]}" -eq 2 ] || { ev "FAIL $r: expected 2 fastq files, got '${rep}'"; touch READS_FAIL; return; }
  for i in 0 1; do
    f=$(basename "${urls[$i]}")
    if [ ! -e "$f.ok" ]; then
      for t in 1 2 3 4; do
        curl -s -L -C - --retry 5 -o "$f" "https://${urls[$i]}" && [ "$(md5sum "$f" | cut -d' ' -f1)" = "${md5s[$i]}" ] && { touch "$f.ok"; break; }
        [ $t -ge 2 ] && rm -f "$f"; sleep $((t * 30))
      done
    fi
    [ -e "$f.ok" ] || { ok=0; ev "FAIL $r: download of $f"; }
  done
  if [ $ok -eq 1 ]; then touch READS_OK; ev "READS $r ok"; else touch READS_FAIL; fi
}
[ "$DL" = 1 ] && for r in "$@"; do rm -f "$B/$r/READS_FAIL"; [ -e "$B/$r/DONE" ] || dl "$r" & done   # all downloads in parallel
for r in "$@"; do
  d=$B/$r; mkdir -p "$d"; cd "$d" || exit 1
  [ -e DONE ] && continue
  until [ -e READS_OK ] || [ -e READS_FAIL ]; do sleep 60; done
  [ -e READS_OK ] || continue
  until [ "$(awk '/MemAvailable/{print $2 * 1024}' /proc/meminfo)" -gt $((MEM + 4000000000)) ]; do sleep 300; done
  rm -rf megahit_out tmp; mkdir -p tmp; t0=$(date +%s)
  python3 "$MH" -1 "${r}_1.fastq.gz" -2 "${r}_2.fastq.gz" --min-contig-len 1000 -t "$THREADS" -m "$MEM" \
        -o megahit_out --tmp-dir "$d/tmp" > megahit.stdout 2>&1
  rc=$?; rm -rf tmp
  if [ $rc -eq 0 ] && [ -s megahit_out/final.contigs.fa ]; then
    cp megahit_out/final.contigs.fa "$r.min1000.fa" && rm -rf megahit_out/intermediate_contigs
    echo "$(( $(date +%s) - t0 ))" > DONE; ev "DONE $r $(( $(date +%s) - t0 ))s"
  else
    ev "FAIL $r megahit rc=$rc"; tail -5 megahit.stdout >> "$B/events.log"
  fi
done
wait; ev "ALL finished"
