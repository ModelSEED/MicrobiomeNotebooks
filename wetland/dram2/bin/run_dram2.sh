#!/usr/bin/env bash
# One DRAM2 launch: rename -> call (prodigal) -> annotate -> distill -> product, for every
# *.fna in <input_dir>, from its own launch dir (so .nextflow history / -resume are per run).
#
#   run_dram2.sh <run_name> <input_dir> [extra nextflow args]
#
# Writes   launch/<run>/                (launch dir: databases symlink, .nextflow/)
#          runs/<run>/janitor.{log,state} (work-dir trimming by bin/janitor.sh, see there)
#          runs/<run>/out/              (published DRAM2 outputs)
#          runs/<run>/DONE | FAILED     (marker, written after verification)
#          logs/<run>/nextflow.log, logs/<run>/stdout.log, logs/<run>/time.txt
#          /scratch/afreiburger/wetland_dram2_work/<run>/   (Nextflow work dir, local disk)
# Re-running the same command resumes (-resume) from the cached tasks.
set -uo pipefail
source /scratch1/afreiburger/wetland_dram2/bin/dram2_settings.sh

RUN=$1; IN=$(readlink -f "$2"); shift 2
L=$B/launch/$RUN; R=$B/runs/$RUN; O=$R/out; W=$WORK_ROOT/$RUN; G=$B/logs/$RUN
mkdir -p "$L" "$O" "$W" "$G" "$TMPDIR" "$MAMBA_ROOT_PREFIX"
ln -sfn "$D/databases" "$L/databases"
rm -f "$R/DONE" "$R/FAILED"
if [ -e "$R/LOWDISK" ]; then
  echo "[$(date -Is)] $R/LOWDISK exists (work disk ran low earlier); free space, delete it, then rerun" | tee -a "$G/stdout.log"
  exit 2
fi
free=$(df -BG --output=avail "$WORK_ROOT" | tail -1 | tr -dc 0-9)
if [ "$free" -lt "${MIN_START_GB:-200}" ]; then
  echo "$(date -Is) only ${free}GB free on $WORK_ROOT (< ${MIN_START_GB:-200}GB); not starting" | tee "$R/FAILED" | tee -a "$G/stdout.log"
  exit 1
fi

n_in=$(ls "$IN"/*.fna | wc -l)
dbargs=""; for db in $DRAM2_DBS; do dbargs="$dbargs --use_$db"; done

cmd=( nextflow-classic -log "$G/nextflow.log" run "$PIPELINE"
      -profile conda -c "$B/launch/wetland_dram2.config"
      -work-dir "$W" -resume -ansi-log false
      --rename --call --annotate
      --input_fasta "$IN" --fasta_fmt '*.fna'
      --outdir "$O"
      --threads 2
      $dbargs
      --distill_topic "$DISTILL_TOPIC" --distill_ecosystem "$DISTILL_ECOSYSTEM"
      --product
      --custom_config_base "$B/launch/nfcore_offline"
      "$@" )

cd "$L"
{
  echo "[$(date -Is)] host=$(hostname) run=$RUN inputs=$n_in"
  echo "[$(date -Is)] cmd: ${cmd[*]}"
} | tee -a "$G/stdout.log"
printf '%q ' "${cmd[@]}" > "$R/command.sh"; echo >> "$R/command.sh"

# work-dir janitor (trims finished tasks only; stops the run if the disk gets low)
"$B/bin/janitor.sh" "$W" "$R" >> "$G/janitor.stdout" 2>&1 < /dev/null &
jpid=$!

/usr/bin/time -v -o "$G/time.txt" "${cmd[@]}" >> "$G/stdout.log" 2>&1
rc=$?
kill "$jpid" 2>/dev/null; wait "$jpid" 2>/dev/null
echo "[$(date -Is)] nextflow exit=$rc" | tee -a "$G/stdout.log"

# ---- verification -------------------------------------------------------------------------
ok=1; why=""
ann=$O/RAW/raw-annotations.tsv
if [ $rc -ne 0 ]; then ok=0; why="nextflow exit $rc"; fi
if [ -e "$R/LOWDISK" ]; then ok=0; why="stopped by janitor: low disk ($(cat "$R/LOWDISK"))"; fi
if [ $ok = 1 ] && [ ! -s "$ann" ]; then ok=0; why="missing $ann"; fi
if [ $ok = 1 ]; then
  # every input MAG must appear in the input_fasta column
  col=$(head -1 "$ann" | tr '\t' '\n' | grep -nx input_fasta | cut -d: -f1)
  n_out=$(tail -n +2 "$ann" | cut -f"$col" | sort -u | grep -c .)
  missing=$(comm -23 <(ls "$IN"/*.fna | xargs -n1 basename | sed 's/\.fna$//' | sort) \
                      <(tail -n +2 "$ann" | cut -f"$col" | sort -u))
  if [ -n "$missing" ]; then ok=0; why="MAGs absent from annotations: $(echo $missing | wc -w)"; fi
  echo "[$(date -Is)] annotations: rows=$(($(wc -l < "$ann")-1)) mags=$n_out/$n_in" | tee -a "$G/stdout.log"
fi
if [ $ok = 1 ]; then
  # per MAG: annotation rows == Prodigal genes, and >= 1 KOfam hit (combine_annotations.py turns a
  # per-database read error into an empty frame, so presence in input_fasta alone is not enough)
  kcol=$(head -1 "$ann" | tr '\t' '\n' | grep -nx kofam_id | cut -d: -f1)
  awk -F'\t' -v c="$col" -v k="${kcol:-0}" 'NR>1{n[$c]++; if(k && $k!="") h[$c]++} END{for(m in n) print m, n[m], h[m]+0}' "$ann" \
    | while read -r m n h; do
        g=$(grep -c '>' "$O/Prodigal_v2.6.3/${m}_called_genes.faa" 2>/dev/null || echo 0)
        if [ "$n" != "$g" ] || [ "$h" -eq 0 ]; then echo "$m rows=$n genes=$g kofam_hits=$h"; fi
      done > "$G/verify_bad.txt"
  if [ -s "$G/verify_bad.txt" ]; then ok=0; why="rows/genes/kofam check failed for $(wc -l < "$G/verify_bad.txt") MAGs ($G/verify_bad.txt)"; fi
  if grep -q "Error loading" "$O/RAW/combine_annotations.log" 2>/dev/null; then ok=0; why="$why; read errors in RAW/combine_annotations.log"; fi
  nk=$(ls "$O/HMM_search_kofam_raw" 2>/dev/null | wc -l)
  [ "$nk" = "$n_in" ] || echo "[$(date -Is)] WARNING: $nk/$n_in raw KOfam domtblouts in $O/HMM_search_kofam_raw" | tee -a "$G/stdout.log"
fi
if [ $ok = 1 ]; then
  dx=$(ls "$O"/*/distillate.xlsx 2>/dev/null | head -1); px=$(ls "$O"/*/product.tsv 2>/dev/null | head -1)
  [ -s "$dx" ] || { ok=0; why="no distillate.xlsx"; }
  [ -s "$px" ] || { ok=0; why="$why no product.tsv"; }
fi
if [ $ok = 1 ]; then
  date -Is > "$R/DONE"
  echo "[$(date -Is)] VERIFIED -> $R/DONE" | tee -a "$G/stdout.log"
  exit 0
else
  echo "$(date -Is) $why" > "$R/FAILED"
  echo "[$(date -Is)] FAILED: $why" | tee -a "$G/stdout.log"
  exit 1
fi
