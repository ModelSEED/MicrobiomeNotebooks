#!/usr/bin/env bash
# GTDB-Tk 2.8.0 / r232 classification of the two Halomonas isolates (I1, I4), whose manuscript
# taxonomy came from GTDB-Tk v1.7.0. Same venv, tools and reference data as run_gtdbtk_local.sh
# (run that first; this script only reuses what it installed).
# Input: data/I1_SPAdes.Assembly.fa, data/I4_SPAdes.Assembly.fa
# Output: MAG_QC/gtdbtk_isolates_r232/ (gtdbtk.bac120.summary.tsv)
# Usage: setsid nohup bash run_gtdbtk_isolates.sh > MAG_QC/gtdbtk_isolates_r232.log 2>&1 &
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
VENV="${VENV:-$HOME/Documents/gtdbtk_venv}"
DATA="${DATA:-$HOME/Documents/gtdbtk_data}"
OUT="$HERE/MAG_QC/gtdbtk_isolates_r232"
TMP="${TMP_DIR:-/tmp/gtdbtk_isolates}"   # short: AF_UNIX socket paths are capped at 107 bytes
CPUS="${CPUS:-32}"
PPLACER_CPUS="${PPLACER_CPUS:-8}"

log() { echo "[$(date '+%F %T')] $*"; }
export PATH="$VENV/bin:$PATH"
[ -f "$DATA/.r232_extracted" ] || { log "r232 data missing; run run_gtdbtk_local.sh first"; exit 1; }
GTDBTK_DATA_PATH="$(find "$DATA" -maxdepth 2 -name metadata -type d -printf '%h\n' -quit)"
export GTDBTK_DATA_PATH

mkdir -p "$OUT" "$TMP"
export TMPDIR="$TMP"
batch="$OUT/batchfile.tsv"
printf '%s\t%s\n' "$HERE/data/I1_SPAdes.Assembly.fa" I1 "$HERE/data/I4_SPAdes.Assembly.fa" I4 > "$batch"
log "gtdbtk classify_wf on the isolates (cpus=$CPUS, pplacer_cpus=$PPLACER_CPUS)"
gtdbtk classify_wf --batchfile "$batch" --out_dir "$OUT" --cpus "$CPUS" \
    --pplacer_cpus "$PPLACER_CPUS" --tmpdir "$TMP"
log "done: $(ls "$OUT"/gtdbtk.*.summary.tsv 2>/dev/null | tr '\n' ' ' || true)"
