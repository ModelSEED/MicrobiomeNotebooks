#!/usr/bin/env bash
# Threshold test: Bakta (same version, DB and flags as run_bakta_assemblies.sh) ahead of the full run, for
#   kbase <run>  -> the KBase assembly of a re-assembled sample, written to prio_runs/<asm> (its own directory: the
#                   full run may reach the same assembly while this runs and would clear runs/<asm>);
#   reasm <fa>   -> the 1-2 kb contigs of a re-assembly (min 1000), written to reasm_runs/<name>.
# Usage (poplar): cd /scratch1/afreiburger/wetland_bakta_asm && setsid nohup bash bakta_priority.sh JOBS kbase:SRR.. reasm:/path/x.fa ... \
#                   > logs/priority.log 2>&1 < /dev/null &
set -uo pipefail
B=/scratch1/afreiburger/wetland_bakta_asm; E=/scratch1/afreiburger/tools/bakta_1.11.4
DB=/scratch/afreiburger/bakta_db_6.0_full; TMP=/scratch/afreiburger/bakta_tmp_asm; THREADS=16
JOBS=$1; shift
export PATH=$E/bin:$PATH B DB TMP THREADS
ev() { echo "[$(date '+%F %T')] $*" >> "$B/priority.events"; }; export -f ev
run1() {
  local kind=${1%%:*} arg=${1#*:} f a out t0 rc
  if [ "$kind" = kbase ]; then f=$(ls "$B"/inputs/"$arg"_*.fna 2>/dev/null | head -1); [ -n "$f" ] || { ev "MISSING input $arg"; return; }
    a=$(basename "$f" .fna); out=$B/prio_runs/$a; mkdir -p "$B/prio_runs"
  else f=$arg; a=$(basename "$f" .fa); out=$B/reasm_runs/$a; mkdir -p "$B/reasm_runs"; fi
  [ -e "$out/DONE" ] && { ev "SKIP $a (done)"; return; }
  rm -rf "$out.tmp" "$TMP/$a.p"; mkdir -p "$TMP/$a.p"; t0=$(date +%s)
  bakta --db "$DB" --output "$out.tmp" --prefix "$a" --keep-contig-headers --meta --translation-table 11 \
        --threads "$THREADS" --skip-plot --tmp-dir "$TMP/$a.p" "$f" > "$B/logs/$a.prio.stdout" 2>&1
  rc=$?; rm -rf "$TMP/$a.p"
  if [ $rc -eq 0 ] && [ -s "$out.tmp/$a.tsv" ] && [ -s "$out.tmp/$a.json" ]; then
    rm -rf "$out"; mv "$out.tmp" "$out" && echo "$(( $(date +%s) - t0 ))" > "$out/DONE"; rm -f "$out.FAILED"; ev "DONE $a $(( $(date +%s) - t0 ))s"
  else ev "FAILED $a rc=$rc"; fi
}
export -f run1
printf '%s\n' "$@" | xargs -P "$JOBS" -I{} bash -c 'run1 "$@"' _ {}
ev "batch finished: $*"
