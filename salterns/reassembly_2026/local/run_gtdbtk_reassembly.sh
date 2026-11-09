#!/usr/bin/env bash
# GTDB-Tk 2.8.0 / GTDB R232 classify_wf over the medium-quality-or-better re-binned MAGs (CheckM >= 50% complete,
# < 10% contaminated) of KBase narrative 275645. Same environment and settings as the wetland runs
# (~/gtdbtk_kbase/run_gtdbtk_run2.sh): ~/gtdbtk_kbase/env.sh, --cpus 48 --pplacer_cpus 8 --min_perc_aa 10.
# Launched by ../pipeline.py (stage gtdbtk) inside the memory-capped systemd user scope gtdbtk-salterns-reassembly
# (MemoryMax=170G, MemorySwapMax=0) once MemAvailable >= 150 GB. Usage: run_gtdbtk_reassembly.sh BATCHFILE OUTDIR
set -euo pipefail
cd "$(dirname "$0")"
source "$HOME/gtdbtk_kbase/env.sh"
CPUS="${CPUS:-48}"
PPLACER_CPUS="${PPLACER_CPUS:-8}"
batch="$1"
OUT="$2"
mkdir -p "$TMPDIR" "$(dirname "$OUT")"
want=$(wc -l < "$batch")
[ "$(cut -f2 "$batch" | sort | uniq -d | wc -l)" -eq 0 ] || { echo "duplicate genome names"; exit 1; }
while IFS=$'\t' read -r f n; do [ -s "$f" ] || { echo "missing $f"; exit 1; }; done < "$batch"
echo "[$(date '+%F %T')] classify_wf on $want genomes (cpus=$CPUS, pplacer_cpus=$PPLACER_CPUS) in $(cat /proc/self/cgroup)"
free -g | head -2
rc=0
gtdbtk classify_wf --batchfile "$batch" --out_dir "$OUT" --cpus "$CPUS" \
    --pplacer_cpus "$PPLACER_CPUS" --min_perc_aa 10 --tmpdir "$TMPDIR" || rc=$?
echo "[$(date '+%F %T')] classify_wf exit=$rc"
echo "$rc" > "$OUT.exit_code"
exit "$rc"
