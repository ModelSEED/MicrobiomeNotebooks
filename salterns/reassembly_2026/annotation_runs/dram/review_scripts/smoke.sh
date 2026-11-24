#!/usr/bin/env bash
# Reviewer smoke test (2026-11-11): DRAM annotate + distill on the test MAG with the phase-A databases
# (Pfam not yet configured, so DRAM skips it). Writes only to $W/review/. Does not touch dram_out/, state/, logs/.
set -uo pipefail
source /scratch1/afreiburger/salterns_2026/dram/scripts/dram_env.sh
R=$W/review; mkdir -p "$R/tmp"; export TMPDIR=$R/tmp
n=$TEST_MAG; out=$R/annot_$n
echo "[$(date '+%F %T')] smoke start pid $$"
DRAM-setup.py print_config > "$R/print_config.txt" 2>&1; echo "print_config rc=$?"
rm -rf "$out"
nice -n 10 /usr/bin/time -v DRAM.py annotate -i "$MAGS/$n.fna" -o "$out" --min_contig_size $MIN_CONTIG \
  --threads $THREADS --gtdb_taxonomy "$TAX" --checkm_quality "$QUAL" > "$R/annotate.log" 2>&1
echo "[$(date '+%F %T')] annotate rc=$?"
rm -rf "$R/distill"
nice -n 10 DRAM.py distill -i "$out/annotations.tsv" -o "$R/distill" $( [ -s "$out/trnas.tsv" ] && echo --trna_path "$out/trnas.tsv" ) \
  $( [ -s "$out/rrnas.tsv" ] && echo --rrna_path "$out/rrnas.tsv" ) --custom_distillate "$SHEET" > "$R/distill.log" 2>&1
echo "[$(date '+%F %T')] distill rc=$?"
python scripts/check_test.py "$out" "$R/distill" > "$R/check_test.log" 2>&1
echo "[$(date '+%F %T')] check_test rc=$? (expected: only 'missing column pfam_hits')"
echo SMOKE_FINISHED
