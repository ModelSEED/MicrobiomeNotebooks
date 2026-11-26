#!/usr/bin/env bash
# Full Bakta 1.11.4 (DB 6.0-full) on all 781 metagenome assemblies (per Andrew, 2026-11-25: "run bakta on all
# metagenomes now, even if it takes weeks"). Coding and non-coding genes on every contig KBase kept (>= 2 kb: the KBase
# MEGAHIT app ran with min_contig_length = 2000), so the same output answers what MAGs miss and the 2.0 vs 2.5 kb part
# of the threshold question (annotations split by contig length afterwards). The 1-2 kb band is not in these
# assemblies; it is tested on re-assemblies (unbinned/reassembly/).
# --meta (pyrodigal metagenome mode: mixed communities, no single-genome training); otherwise the MAG run's flags.
# Inputs arrive while it runs: assemblies are downloaded on the workstation and pushed here by sync_assemblies.sh;
# the driver keeps polling inputs/ until inputs/ALL_INPUTS_SYNCED exists and every input is DONE or FAILED.
# Resumable: runs/<asm>/DONE; a failure leaves runs/<asm>.FAILED and is retried once on the next pass.
# Temp files on poplar's local /scratch (97% full, so at most JOBS assemblies' temp at a time); outputs on /scratch1.
# Usage (poplar): cd /scratch1/afreiburger/wetland_bakta_asm && setsid nohup bash run_bakta_assemblies.sh [JOBS] [THREADS] \
#                   > logs/driver.log 2>&1 < /dev/null &
set -uo pipefail
B=/scratch1/afreiburger/wetland_bakta_asm
E=/scratch1/afreiburger/tools/bakta_1.11.4
DB=/scratch/afreiburger/bakta_db_6.0_full
TMP=/scratch/afreiburger/bakta_tmp_asm
JOBS=${1:-6}; THREADS=${2:-16}
cd "$B"; mkdir -p inputs runs logs "$TMP"
export PATH=$E/bin:$PATH B DB TMP THREADS
ev() { echo "[$(date '+%F %T')] $*" >> "$B/bakta_asm.events"; }
[ "$(bakta --version 2>&1 | tail -1)" = "bakta 1.11.4" ] || { ev "FAIL wrong bakta"; exit 1; }
exec 9> "$B/logs/driver.lock"; flock -n 9 || { echo "another driver is running"; exit 1; }
one() {
  local f=$1 a; a=$(basename "$f" .fna); local out=$B/runs/$a t0 rc
  [ -e "$out/DONE" ] && return 0
  [ -e "$out.FAILED" ] && [ "$(grep -c '^attempt' "$out.FAILED")" -ge 2 ] && return 0
  local att=$(( $( [ -e "$out.FAILED" ] && grep -c '^attempt' "$out.FAILED" || echo 0) + 1 ))
  rm -rf "$out" "$out.tmp" "$TMP/$a"; mkdir -p "$TMP/$a"; t0=$(date +%s)
  bakta --db "$DB" --output "$out.tmp" --prefix "$a" --keep-contig-headers --meta --translation-table 11 \
        --threads "$THREADS" --skip-plot --tmp-dir "$TMP/$a" "$f" > "$B/logs/$a.stdout" 2>&1
  rc=$?; rm -rf "$TMP/$a"
  if [ $rc -eq 0 ] && [ -s "$out.tmp/$a.tsv" ] && [ -s "$out.tmp/$a.json" ]; then
    mv "$out.tmp" "$out" && echo "$(( $(date +%s) - t0 ))" > "$out/DONE"; rm -f "$out.FAILED"
    echo "[$(date '+%F %T')] DONE $a $(( $(date +%s) - t0 ))s" >> "$B/bakta_asm.events"
  else
    { echo "attempt $att rc=$rc"; tail -5 "$B/logs/$a.stdout"; } >> "$out.FAILED"
    echo "[$(date '+%F %T')] FAILED $a rc=$rc attempt $att" >> "$B/bakta_asm.events"
  fi
}
export -f one
ev "START jobs=$JOBS threads=$THREADS"
while :; do
  todo=$(for f in inputs/*.fna; do [ -e "$f" ] || continue; a=$(basename "$f" .fna); [ -e "runs/$a/DONE" ] && continue
         [ -e "runs/$a.FAILED" ] && [ "$(grep -c '^attempt' "runs/$a.FAILED")" -ge 2 ] && continue; echo "$B/$f"; done)
  if [ -n "$todo" ]; then
    echo "$todo" | xargs -P "$JOBS" -I{} bash -c 'one "$@"' _ {}
  elif [ -e inputs/ALL_INPUTS_SYNCED ]; then
    ev "ALL done=$(ls runs/*/DONE 2>/dev/null | wc -l) failed=$(ls runs/*.FAILED 2>/dev/null | wc -l) inputs=$(ls inputs/*.fna | wc -l)"; break
  else
    sleep 300
  fi
done
