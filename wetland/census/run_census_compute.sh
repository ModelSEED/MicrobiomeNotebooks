#!/usr/bin/env bash
# Weisenhorn wetland MAG census (handoff Phases 3-4), compute on poplar. Read-only on shared data; writes only under $W.
# 1. CheckM2 1.1.0 (DB uniref100.KO.1, v3) on all 12,859 MAGs                   -> checkm2/quality_report.tsv
# 2. skani 0.3.1 triangle over all MAGs (sparse, screen 90%)                     -> skani_triangle.tsv  (MAG-vs-MAG clustering)
# 3. skani 0.3.1 search vs KBDL novel_all_representatives (ref:1be8415b...) and novel_hq_v2 (ref:9f4ac5e9...) -> skani_vs_<db>.tsv
# 4. skani 0.3.1 search vs KBDL gtdb_r220_consolidated (ref:cee67d7f...) (compare the handoff's r220 route with local R232)
# skani 0.3.1 = the version the KBDL reference DBs were built with. AF is kept (columns Align_fraction_ref/query, percent).
# Steps are resumable (.<step>.done). Usage: cd $W && setsid nohup bash run_census_compute.sh > run.log 2>&1 < /dev/null &
set -uo pipefail
W=/scratch1/afreiburger/wetland_census
cd "$W"
MAGS=/scratch1/afreiburger/wetland_coverm_pilot/all_mags
QC=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/8c703f8abdc13cdc1ef374a2bcf37f14_
CK2DB=/scratch1/afreiburger/emilykin/processed/refdata/checkm2/CheckM2_database/uniref100.KO.1.dmnd
SKANI=/scratch1/chenry/opt/skani/0.3.1/skani
declare -A DB=( [novel_all_representatives]=/scratch1/chenry/kbdl/skani_dbs/ad1c32979e1b4cc897c5f32631e2cc9e
                [novel_hq_v2]=/scratch1/chenry/kbdl/skani_dbs/c63e6b662bbf402e901974f84cd29d31
                [gtdb_r220_consolidated]=/scratch1/fliu/hub_scratch/fliu/data/ani/skani/gtdb_r220 )
ev() { echo "[$(date '+%F %T')] $*" >> census.events; }
tm() { local tag=$1; shift; /usr/bin/time -v -o "time_${tag}.txt" "$@"; local rc=$?; \
       ev "TIME $tag rc=$rc wall=$(grep 'Elapsed (wall' time_${tag}.txt | awk '{print $NF}') maxrss_GB=$(awk '/Maximum resident/ {printf "%.1f", $NF/1048576}' time_${tag}.txt)"; return $rc; }
ev "START load $(cut -d' ' -f1-3 /proc/loadavg)"
ls "$MAGS"/*.fna > mags.list
n=$(wc -l < mags.list); [ "$n" -eq 12859 ] || { ev "FAIL expected 12859 MAGs, found $n"; exit 1; }
md5sum "$SKANI" > skani_binary.md5; "$SKANI" --version > skani_version.txt

# 1. CheckM2 in the background (its own tmp on local, non-NFS /scratch)
if [ ! -f .checkm2.done ]; then
  ( mkdir -p /scratch/afreiburger/tmp_checkm2_wetland
    tm checkm2 env PATH=$QC/bin:$PATH TMPDIR=/scratch/afreiburger/tmp_checkm2_wetland checkm2 predict -i "$MAGS" -x fna \
        -t 32 --database_path "$CK2DB" --tmpdir /scratch/afreiburger/tmp_checkm2_wetland -o checkm2 --force \
      && [ -s checkm2/quality_report.tsv ] && touch .checkm2.done ) &
  CK_PID=$!
fi

# 2-4. skani (16 threads each, one after another)
if [ ! -f .triangle.done ]; then
  tm triangle "$SKANI" triangle -l mags.list -t 16 -E -s 90 -o skani_triangle.tsv && touch .triangle.done
fi
for name in novel_all_representatives novel_hq_v2 gtdb_r220_consolidated; do
  [ -f ".search_${name}.done" ] && continue
  tm "search_${name}" "$SKANI" search --ql mags.list -d "${DB[$name]}" -t 16 -n 5 -o "skani_vs_${name}.tsv" \
     && touch ".search_${name}.done"
done
[ -n "${CK_PID:-}" ] && wait "$CK_PID"
ev "DONE checkm2=$([ -f .checkm2.done ] && echo ok || echo FAIL) triangle=$([ -f .triangle.done ] && echo ok || echo FAIL) searches=$(ls .search_*.done 2>/dev/null | wc -l)/3"
