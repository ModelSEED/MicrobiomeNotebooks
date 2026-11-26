#!/usr/bin/env bash
# GTDB-Tk 2.8.0 / GTDB R232 classify_wf on the 291 EmilyKin MAGs. Same environment and settings as the salterns
# re-assembly and wetland runs (../../local/run_gtdbtk_reassembly.sh: ~/gtdbtk_kbase/env.sh, --pplacer_cpus 8,
# --min_perc_aa 10), but --cpus 24 (shared-machine rule). Launched by 20_launch_gtdbtk.sh inside the memory-capped
# systemd user scope gtdbtk-emilykin-r232 (MemoryMax=170G, MemorySwapMax=0). Usage: run_gtdbtk_emilykin.sh BATCHFILE OUTDIR
set -euo pipefail
source "$HOME/gtdbtk_kbase/env.sh"          # PATH -> ~/Documents/gtdbtk_venv/bin, GTDBTK_DATA_PATH -> release232, TMPDIR
CPUS="${CPUS:-24}"
PPLACER_CPUS="${PPLACER_CPUS:-8}"
batch="$1"; OUT="$2"
mkdir -p "$TMPDIR" "$(dirname "$OUT")"
[ "$(cut -f2 "$batch" | sort | uniq -d | wc -l)" -eq 0 ] || { echo "duplicate genome names"; exit 1; }
while IFS=$'\t' read -r f n; do [ -s "$f" ] || { echo "missing $f"; exit 1; }; done < "$batch"
grep -q '^VERSION_DATA=r232' "$GTDBTK_DATA_PATH/metadata/metadata.txt" || { echo "GTDBTK_DATA_PATH is not r232"; exit 1; }
echo "[$(date '+%F %T')] $(gtdbtk --version 2>&1 | head -1); data $GTDBTK_DATA_PATH"
echo "[$(date '+%F %T')] classify_wf on $(wc -l < "$batch") genomes (cpus=$CPUS, pplacer_cpus=$PPLACER_CPUS) in $(cat /proc/self/cgroup)"
free -g | head -2
rc=0
extra=()                                     # SCRATCH=<dir>: pplacer memory-maps to disk (much lower RAM, slower)
[ -n "${SCRATCH:-}" ] && { mkdir -p "$SCRATCH"; extra=(--scratch_dir "$SCRATCH"); }
gtdbtk classify_wf --batchfile "$batch" --out_dir "$OUT" --cpus "$CPUS" \
    --pplacer_cpus "$PPLACER_CPUS" --min_perc_aa 10 --tmpdir "$TMPDIR" "${extra[@]}" || rc=$?
echo "[$(date '+%F %T')] classify_wf exit=$rc"
echo "$rc" > "$OUT.exit_code"
exit "$rc"
