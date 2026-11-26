#!/usr/bin/env bash
# Full Bakta annotation (coding + non-coding) of the 1,184 MAGs that DRAM2 annotated (per Andrew, 2026-11-18: "run
# bakta on the same representative set that was used for the DRAM annotations"; census/outputs/dram2_hq_multisample.tsv).
#
# Bakta 1.11.4 (conda env /scratch1/afreiburger/tools/bakta_1.11.4, = KBDL's pinned version; AMRFinderPlus pinned to
# 4.0.23 because 4.2.7, the solver's default, requires AMR DB >= 2025-09-22.2 and DB 6.0-full ships 2025-06-03.1) with DB 6.0-full
# (/scratch1/fliu/hub_scratch/data/bakta/db, = KBDL's BaktaDB ref:2fdfb17d...), on the same FASTAs DRAM2 read
# (/scratch1/afreiburger/wetland_coverm_pilot/all_mags/<mag>.fna). Flags match KBDL's structural driver where it
# matters for gene calls (--keep-contig-headers --min-contig-length 1 --translation-table 11, no --meta/--complete/
# --compliant), but this is a FULL run: functional annotation of every CDS (UPS/IPS/PSC/PSCC, Pfam, AMRFinderPlus,
# expert systems, pseudogenes) and every non-coding caller (tRNA, tmRNA, rRNA, ncRNA, ncRNA regions, CRISPR, sORF,
# gaps, oriC/oriV/oriT). Only the circular plot is skipped.
#
# Resumable: runs/<mag>/DONE marks a finished MAG; a failure leaves runs/<mag>.FAILED (rc + log tail) and is retried
# on the next launch. Temp files on poplar's local disk (/scratch/afreiburger/bakta_tmp).
# Usage (poplar): cd /scratch1/afreiburger/wetland_bakta && setsid nohup bash run_bakta_reps.sh [JOBS] [THREADS] \
#                   > logs/driver.log 2>&1 < /dev/null &
set -uo pipefail
B=/scratch1/afreiburger/wetland_bakta
E=/scratch1/afreiburger/tools/bakta_1.11.4
DB=${BAKTA_DB:-/scratch/afreiburger/bakta_db_6.0_full}   # local copy of /scratch1/fliu/hub_scratch/data/bakta/db
#   (launch_after_db_copy.sh verifies it file by file; the NFS share reads at ~12 MB/s and psc.dmnd alone is 44 GB)
MAGS=/scratch1/afreiburger/wetland_coverm_pilot/all_mags
TMP=/scratch/afreiburger/bakta_tmp
JOBS=${1:-16}
THREADS=${2:-6}
cd "$B"
mkdir -p runs logs "$TMP"
export PATH=$E/bin:$PATH B E DB MAGS TMP THREADS
ev() { echo "[$(date '+%F %T')] $*" >> "$B/bakta.events"; }
[ "$(bakta --version 2>&1 | tail -1)" = "bakta 1.11.4" ] || { ev "FAIL wrong bakta: $(bakta --version 2>&1)"; exit 1; }
[ "$(python3 -c "import json; d = json.load(open('$DB/version.json')); print(d['major'], d['minor'], d['type'])")" = "6 0 full" ] \
  || { ev "FAIL DB at $DB is not Bakta DB 6.0 full"; exit 1; }
REPS=${REPS:-reps.txt}             # one MAG name per line (REPS=smoke.txt for the launcher's one-MAG smoke test)
[ -s "$REPS" ] || { echo "$REPS (one MAG name per line) missing"; exit 1; }
n=$(grep -c . "$REPS")
ev "START jobs=$JOBS threads=$THREADS reps=$n done_already=$(ls runs/*/DONE 2>/dev/null | wc -l)"

one() {
  local m=$1 out=$B/runs/$1 t0 rc
  [ -e "$out/DONE" ] && return 0
  rm -rf "$out" "$out.tmp" "$TMP/$m"; rm -f "$out.FAILED"; mkdir -p "$TMP/$m"
  t0=$(date +%s)
  bakta --db "$DB" --output "$out.tmp" --prefix "$m" --keep-contig-headers --min-contig-length 1 \
        --translation-table 11 --threads "$THREADS" --skip-plot --tmp-dir "$TMP/$m" \
        "$MAGS/$m.fna" > "$B/logs/$m.stdout" 2>&1
  rc=$?
  rm -rf "$TMP/$m"
  if [ $rc -eq 0 ] && [ -s "$out.tmp/$m.tsv" ] && [ -s "$out.tmp/$m.gbff" ] && [ -s "$out.tmp/$m.json" ]; then
    mv "$out.tmp" "$out" && echo "$(( $(date +%s) - t0 ))" > "$out/DONE"
    echo "[$(date '+%F %T')] DONE $m $(( $(date +%s) - t0 ))s" >> "$B/bakta.events"
  else
    { echo "rc=$rc"; tail -20 "$B/logs/$m.stdout"; } > "$out.FAILED"
    echo "[$(date '+%F %T')] FAILED $m rc=$rc" >> "$B/bakta.events"
  fi
}
export -f one
grep . "$REPS" | xargs -P "$JOBS" -I{} bash -c 'one "$@"' _ {}
d=$(ls runs/*/DONE 2>/dev/null | wc -l); f=$(ls runs/*.FAILED 2>/dev/null | wc -l)
ev "ALL done=$d/$n failed=$f"
