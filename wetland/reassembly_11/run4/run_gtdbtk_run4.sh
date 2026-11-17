#!/usr/bin/env bash
# GTDB-Tk 2.8.0 / GTDB R232 classify_wf over every MAG in all_manifest_run4.tsv: all MAGs of KBase narrative 272312
# extracted from the 11 re-assembled libraries, plus 60 controls. Launched by pipeline_run4.sh inside the memory-capped
# systemd scope gtdbtk-wetland-run4. Same settings as run 1 (run_gtdbtk_all.sh) except --cpus 48.
set -euo pipefail
cd "$(dirname "$0")"
source ./env.sh
CPUS="${CPUS:-48}"
PPLACER_CPUS="${PPLACER_CPUS:-8}"
OUT=gtdbtk_run4
batch=batchfile_run4.tsv

mkdir -p "$TMPDIR"                 # gtdbtk rejects a --tmpdir that does not exist
awk -F'\t' -v d="$PWD" 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
                        { print d "/" $c["file"] "\t" $c["name"] }' all_manifest_run4.tsv > "$batch"
want=$(( $(wc -l < all_manifest_run4.tsv) - 1 ))
[ "$(wc -l < "$batch")" -eq "$want" ] || { echo "batchfile has $(wc -l < "$batch") rows, expected $want"; exit 1; }
[ "$(cut -f2 "$batch" | sort | uniq -d | wc -l)" -eq 0 ] || { echo "duplicate genome names"; exit 1; }
while IFS=$'\t' read -r f n; do [ -s "$f" ] || { echo "missing $f"; exit 1; }; done < "$batch"
echo "[$(date '+%F %T')] classify_wf on $want genomes (cpus=$CPUS, pplacer_cpus=$PPLACER_CPUS) in $(cat /proc/self/cgroup)"
free -g | head -2
rc=0
gtdbtk classify_wf --batchfile "$batch" --out_dir "$OUT" --cpus "$CPUS" \
    --pplacer_cpus "$PPLACER_CPUS" --min_perc_aa 10 --tmpdir "$TMPDIR" || rc=$?
echo "[$(date '+%F %T')] classify_wf exit=$rc"
exit "$rc"
