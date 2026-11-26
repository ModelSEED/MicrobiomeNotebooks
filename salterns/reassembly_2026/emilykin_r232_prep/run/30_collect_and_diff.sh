#!/usr/bin/env bash
# After the run: copy the R232 summaries next to this plan (../results/), diff them against the R220 classification
# MAG by MAG, score the desk prediction, and simulate (read-only) what the EmilyKin iterativeID rule would do with R232.
# Writes nothing into ~/Documents/EmilyKin.  Usage: [WORK=~/gtdbtk_kbase/emilykin_r232] bash 30_collect_and_diff.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=${WORK:-$HOME/gtdbtk_kbase/emilykin_r232}
PY=$HOME/Documents/py_venv/bin/python
RES="$HERE/../results"; mkdir -p "$RES"
rc=$(cat "$WORK/gtdbtk_out.exit_code" 2>/dev/null || echo missing)
[ "$rc" = 0 ] || { echo "classify_wf exit code: $rc (see $WORK/gtdbtk_run.log)"; exit 1; }
ins=()
for dom in bac120 ar53; do
  f="$WORK/gtdbtk_out/gtdbtk.$dom.summary.tsv"
  [ -s "$f" ] && cp -p "$f" "$RES/" && ins+=("$RES/gtdbtk.$dom.summary.tsv")
done
cp -p "$WORK/gtdbtk_out/gtdbtk.log" "$RES/gtdbtk.log"
[ -s "$WORK/gtdbtk_out/gtdbtk.failed_genomes.tsv" ] && cp -p "$WORK/gtdbtk_out/gtdbtk.failed_genomes.tsv" "$RES/"
cp -p "$WORK/staged_md5.txt" "$RES/"
"$PY" -I "$HERE/../diff_r220_r232.py" --r232 "${ins[@]}" --out "$RES"
if [ -s "$HERE/../poplar_inventory/gtdbtk_291.bac120.summary.tsv" ]; then   # also against the server-native R220 summary
  "$PY" -I "$HERE/../diff_r220_r232.py" --r220 "$HERE/../poplar_inventory/gtdbtk_291.bac120.summary.tsv" \
      --r232 "${ins[@]}" --out "$RES/vs_server_native_r220" > /dev/null
fi
"$PY" -I "$HERE/../simulate_iterativeids_r232.py" --r232-summary "$RES/gtdbtk.r232.summary.tsv" \
    --out "$RES/iterativeid_simulation_r232.tsv" | tee "$RES/iterativeid_simulation_r232.log"
