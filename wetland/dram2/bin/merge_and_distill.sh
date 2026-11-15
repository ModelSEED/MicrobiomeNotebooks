#!/usr/bin/env bash
# After all chunks are DONE:
#  1) merge the per-chunk raw-annotations.tsv into merged/raw-annotations.all.tsv (+ `mag` column)
#  2) one more DRAM2 launch (no --call/--annotate) that runs DRAM2's own distill + product on the
#     merged table for all MAGs, with the per-MAG rRNA/tRNA tables from every chunk.
#
#   setsid nohup bin/merge_and_distill.sh > logs/merge_and_distill.log 2>&1 < /dev/null &
#
# Outputs: merged/raw-annotations.all.tsv
#          runs/all_mags/out/{distill or DISTILL}/distillate.xlsx, summarized_genomes.tsv
#          runs/all_mags/out/PRODUCT/product.tsv, product.html ; runs/all_mags/DONE
set -uo pipefail
source /scratch1/afreiburger/wetland_dram2/bin/dram2_settings.sh
cd "$B"
RUN=${RUN:-all_mags}
PY=$B/conda_cache/dram2_beta17_distill/bin/python

[ -e runs/ALL_CHUNKS_DONE ] || [ "${FORCE:-0}" = 1 ] || { echo "runs/ALL_CHUNKS_DONE missing (FORCE=1 to override)"; exit 1; }

echo "[$(date -Is)] merging"
"$PY" -I bin/merge_dram2.py ${MERGE_ARGS:-} || { echo "merge failed"; exit 1; }
ANN=$B/merged/raw-annotations.all.tsv

L=$B/launch/$RUN; R=$B/runs/$RUN; O=$R/out; W=$WORK_ROOT/$RUN; G=$B/logs/$RUN
mkdir -p "$L" "$O" "$W" "$G"; ln -sfn "$D/databases" "$L/databases"; rm -f "$R/DONE" "$R/FAILED"
cmd=( nextflow-classic -log "$G/nextflow.log" run "$PIPELINE"
      -profile conda -c "$B/launch/wetland_dram2.config"
      -work-dir "$W" -resume -ansi-log false
      --annotations "$ANN"
      --rrnas "$B/merged/rrnas" --trnas "$B/merged/trnas"
      --distill_topic "$DISTILL_TOPIC" --distill_ecosystem "$DISTILL_ECOSYSTEM"
      --product
      --outdir "$O" --threads 2
      --custom_config_base "$B/launch/nfcore_offline"
      --validate_params false )
# ^ beta17's schema types --rrnas/--trnas as files, but subworkflows/local/collect_rna.nf globs
#   "${params.rrnas}/*.tsv" (a directory of per-MAG tables), so schema validation must be off here.
cd "$L"
echo "[$(date -Is)] cmd: ${cmd[*]}" | tee -a "$G/stdout.log"
printf '%q ' "${cmd[@]}" > "$R/command.sh"; echo >> "$R/command.sh"
/usr/bin/time -v -o "$G/time.txt" "${cmd[@]}" >> "$G/stdout.log" 2>&1; rc=$?
dx=$(ls "$O"/*/distillate.xlsx 2>/dev/null | head -1); px=$(ls "$O"/*/product.tsv 2>/dev/null | head -1)
if [ $rc -eq 0 ] && [ -s "$dx" ] && [ -s "$px" ]; then
  date -Is > "$R/DONE"; echo "[$(date -Is)] DONE: $ANN, $dx, $px"
else
  echo "$(date -Is) rc=$rc distill=$dx product=$px" > "$R/FAILED"; echo "[$(date -Is)] FAILED rc=$rc"; exit 1
fi
