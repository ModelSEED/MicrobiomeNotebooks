#!/usr/bin/env bash
# Launch the EmilyKin R232 GTDB-Tk run the way reassembly_2026/pipeline.py (gtdbtk_stage) launches the salterns one:
# a transient systemd USER scope with a hard memory cap and no swap, plus ../../local/mem_guard.sh, which freezes ONLY this
# scope if the machine's MemAvailable drops below 12 GB and thaws it above 30 GB (never touches other users' jobs).
# Refuses to start if MemAvailable < 150 GB (the R232 class-level pplacer trees peak at ~133 GB pplacer RSS / 138 GB
# scope anon with 8 pplacer threads: salterns run 2026-11-11), if another GTDB-Tk classify_wf is running, or if the
# scope already exists. WAIT=1 re-checks memory every 10 min for up to WAIT_HOURS (default 12) instead of refusing.
# Must be started from a shell that has the user systemd bus (a normal login / tmux shell on frankensteinsmonster).
# Usage: [WORK=~/gtdbtk_kbase/emilykin_r232] [WAIT=1] bash 20_launch_gtdbtk.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=${WORK:-$HOME/gtdbtk_kbase/emilykin_r232}
SCOPE=gtdbtk-emilykin-r232
MIN_GB=${MIN_GB:-150}
[ -s "$WORK/batchfile.tsv" ] || { echo "no $WORK/batchfile.tsv - run 10_stage_genomes.sh first"; exit 1; }
[ -e "$WORK/gtdbtk_out.exit_code" ] && { echo "$WORK/gtdbtk_out.exit_code exists: a run already finished here"; exit 1; }
st=$(systemctl --user is-active "$SCOPE.scope" 2>/dev/null || true)
case "$st" in active|activating|reloading) echo "scope $SCOPE.scope is already $st"; exit 1;; esac
pgrep -af "gtdbtk classify_w[f]" && { echo "another GTDB-Tk classify_wf is running (above) - not starting a second one"; exit 1; }
avail() { awk '/MemAvailable/ {printf "%d", $2/1048576}' /proc/meminfo; }
deadline=$(( $(date +%s) + ${WAIT_HOURS:-12} * 3600 ))
while [ "$(avail)" -lt "$MIN_GB" ]; do
  [ "${WAIT:-0}" = 1 ] && [ "$(date +%s)" -lt "$deadline" ] || { echo "MemAvailable $(avail) GB < $MIN_GB GB - not starting"; exit 3; }
  echo "[$(date '+%F %T')] MemAvailable $(avail) GB < $MIN_GB GB; re-checking in 10 min"; sleep 600
done
systemctl --user reset-failed "$SCOPE.scope" 2>/dev/null || true
echo "[$(date '+%F %T')] launching in $SCOPE.scope (MemAvailable $(avail) GB)" | tee -a "$WORK/gtdbtk_run.log"
setsid nohup systemd-run --user --scope --unit="$SCOPE" -p MemoryMax=170G -p MemorySwapMax=0 \
    nice -n 10 bash "$HERE/run_gtdbtk_emilykin.sh" "$WORK/batchfile.tsv" "$WORK/gtdbtk_out" \
    >> "$WORK/gtdbtk_run.log" 2>&1 < /dev/null &
echo $! > "$WORK/launcher.pid"
UNIT="$SCOPE.scope" GUARD_DIR="$WORK" setsid nohup bash "$HERE/../../local/mem_guard.sh" > /dev/null 2>&1 < /dev/null &
echo $! > "$WORK/mem_guard.pid"
echo "started; follow: tail -f $WORK/gtdbtk_run.log $WORK/gtdbtk_out/gtdbtk.log $WORK/mem_guard.events"
echo "done when $WORK/gtdbtk_out.exit_code exists; then: bash $HERE/30_collect_and_diff.sh"
