#!/usr/bin/env bash
# QC of the 2026 first-pass site co-assembly bins with the rollout's own QC pipeline (Andrew, 2026-11-19):
# all 718 DAS Tool bins of the 4 MEGAHIT sediment co-assemblies (R1/R2/R2A/SF2_sediment; uploaded md5-checked to
# $C/inputs/<unit>/ from reassembly_2026/work/dastool/<unit>/attempt1-*/out/<unit>_DASTool_bins/), not only the 548
# that KBase CheckM had selected. Same commands, tool envs, databases and thread counts as the rollout's task bodies
# (driver.py: checkm2 / gunc / checkm1 / keep / gtdbtk; values from rollout/scripts/config.json):
#   CheckM2 1.1.0, GUNC 1.1.1 (proGenomes 2.1), CheckM1 1.1.3 --reduced_tree (informational), then rolib.py keep
#   (CheckM2 >= 50 / < 10 and GUNC pass) -> $C/mags_kept/<unit>__<bin>.fa, then GTDB-Tk 2.8.0 / R232 classify_wf on
#   the kept MAGs. Per unit, CheckM2 / GUNC / CheckM1 run in parallel (3 x 16 threads). nice 10, detached:
#   cd $C && nohup setsid nice -n 10 bash coasm_qc.sh > logs/run.log 2>&1 < /dev/null &
# Markers: $C/.done/<step>:<unit>, $C/DONE or $C/FAILED.
set -uo pipefail
C=/scratch1/afreiburger/salterns_2026/coasm_qc
QCENV=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/8c703f8abdc13cdc1ef374a2bcf37f14_/bin
CHECKM2_DB=/scratch1/afreiburger/emilykin/processed/refdata/checkm2/CheckM2_database/uniref100.KO.1.dmnd
CHECKM1=/scratch1/afreiburger/salterns_2026/envs/bt_checkm1/bin
CHECKM1_DATA=/scratch1/afreiburger/salterns_2026/binning_test/refdata/checkm1
GUNC=/scratch1/afreiburger/salterns_2026/envs/bt_gunc/bin
GUNC_DB=/scratch1/afreiburger/salterns_2026/binning_test/refdata/gunc/gunc_db_progenomes2.1.dmnd
PY=/scratch1/afreiburger/salterns_2026/envs/bt_map/bin/python
RO=/scratch1/afreiburger/salterns_2026/rollout/scripts/rolib.py
GT_BIN=/scratch1/afreiburger/salterns_2026/envs/gtdbtk280/bin
GT_DATA=$(cat /scratch1/afreiburger/salterns_2026/refdata/gtdbtk_r232/.ok)
TH=16
cd $C; mkdir -p .done logs qc mags_kept
ts() { date -u +%FT%TZ; }
fail() { echo "$(ts) $1" > FAILED; echo "$(ts) FAILED: $1"; exit 1; }
tmpd() { local d=/tmp/cq_$1_$$; rm -rf $d; mkdir -p $d; echo $d; }
echo "$(ts) start (pid $$)"
for u in R1_sediment R2_sediment R2A_sediment SF2_sediment; do
  IN=$C/inputs/$u; test -e $IN/.ok || fail "no md5-verified inputs for $u"
  n=$(find -L $IN -maxdepth 1 -name '*.fa' | wc -l); echo "$(ts) $u: $n bins"
  pids=()
  if [ ! -e .done/checkm2:$u ]; then (
    O=$C/qc/$u/checkm2; rm -rf $O; mkdir -p $O; T=$(tmpd c2_$u); trap "rm -rf $T" EXIT
    export PATH=$QCENV:$PATH OMP_NUM_THREADS=$TH
    checkm2 predict --input $IN -x fa --output-directory $O/out --threads $TH --database_path $CHECKM2_DB --tmpdir $T > $O/checkm2.log 2>&1 \
      && cp $O/out/quality_report.tsv $O/quality_report.tsv && touch $C/.done/checkm2:$u ) & pids+=($!); fi
  if [ ! -e .done/gunc:$u ]; then (
    O=$C/qc/$u/gunc; rm -rf $O; mkdir -p $O; T=$(mktemp -d $C/logs/gunc_tmp_XXXX); trap "rm -rf $T" EXIT
    export PATH=$GUNC:$PATH OMP_NUM_THREADS=$TH
    gunc run --input_dir $IN --file_suffix .fa --db_file $GUNC_DB --threads $TH --out_dir $O --temp_dir $T > $O/gunc.log 2>&1 \
      && ls $O/GUNC.*.maxCSS_level.tsv > /dev/null && rm -rf $O/diamond_output $O/gene_calls && touch $C/.done/gunc:$u ) & pids+=($!); fi
  if [ ! -e .done/checkm1:$u ]; then (
    O=$C/qc/$u/checkm1; rm -rf $O; mkdir -p $O; T=$(tmpd c1_$u); trap "rm -rf $T" EXIT
    export PATH=$CHECKM1:$PATH CHECKM_DATA_PATH=$CHECKM1_DATA OMP_NUM_THREADS=$TH
    checkm lineage_wf --reduced_tree -x fa -t $TH --pplacer_threads 4 --tab_table -f $O/qa.tsv --tmpdir $T $IN $O/out > $O/checkm1.log 2>&1 \
      && rm -rf $O/out/bins $O/out/storage/aai_qa && touch $C/.done/checkm1:$u ) & pids+=($!); fi
  for p in "${pids[@]}"; do wait $p; done
  for k in checkm2 gunc; do test -e .done/$k:$u || fail "$k:$u (qc/$u/$k/$k.log)"; done
  test -e .done/checkm1:$u || echo "$(ts) WARNING checkm1:$u failed (informational; keep runs without it)"
  r=$(($(wc -l < qc/$u/checkm2/quality_report.tsv) - 1)); [ "$r" -eq "$n" ] || fail "CheckM2 has $r rows for $n bins ($u)"
  G=$(ls qc/$u/gunc/GUNC.*.maxCSS_level.tsv | head -1); g=$(($(wc -l < $G) - 1))
  [ "$g" -eq "$n" ] || echo "$(ts) WARNING: GUNC has $g rows for $n bins in $u (missing rows count as not passing)"
  C1=qc/$u/checkm1/qa.tsv; [ -s "$C1" ] || C1=none
  $PY -I $RO keep --sample $u --bins $IN --checkm2 qc/$u/checkm2/quality_report.tsv --gunc "$G" --checkm1 "$C1" \
     --outdir $C/$u/keep --kept-dir $C/mags_kept > logs/keep_$u.log 2>&1 || fail "keep:$u (logs/keep_$u.log)"
  touch .done/keep:$u; echo "$(ts) $u: kept $(wc -l < $C/$u/keep/kept.txt) of $n"
done
# GTDB-Tk 2.8.0 / R232 on all kept co-assembly MAGs, flags as the rollout's gtdbtk task
if [ ! -e .done/gtdbtk ]; then
  O=$C/gtdbtk; rm -rf $O; mkdir -p $O; T=$(tmpd gt all)
  [ "$(df --output=avail -BG /tmp | tail -1 | tr -dc 0-9)" -ge 10 ] || fail "/tmp has < 10 GB free"
  export PATH=$GT_BIN:$PATH GTDBTK_DATA_PATH=$GT_DATA
  gtdbtk --version 2>&1 | head -1
  find $C/mags_kept -maxdepth 1 -name '*.fa' | sort | awk -F/ '{n=$NF; sub(/\.fa$/,"",n); print $0"\t"n}' > $O/batchfile.tsv
  [ "$(cut -f2 $O/batchfile.tsv | sort | uniq -d | wc -l)" -eq 0 ] || fail "duplicate genome names"
  echo "$(ts) GTDB-Tk on $(wc -l < $O/batchfile.tsv) kept MAGs"
  gtdbtk classify_wf --batchfile $O/batchfile.tsv --out_dir $O/out --cpus 48 --pplacer_cpus 8 --min_perc_aa 10 --tmpdir "$T" > $O/gtdbtk.log 2>&1 \
    || { rm -rf $T; fail "gtdbtk (gtdbtk/gtdbtk.log)"; }
  rm -rf $T; cp $O/out/gtdbtk.*.summary.tsv $O/ && touch .done/gtdbtk
  echo "$(ts) GTDB-Tk classified $(cat $O/gtdbtk.*.summary.tsv | awk 'NR>1 && $1!="user_genome"' | wc -l)"
fi
echo "$(ts) units=4 bins=$(cat inputs/*/.md5 | wc -l) kept=$(ls mags_kept/*.fa | wc -l)" > DONE; cat DONE
