#!/usr/bin/env bash
# Census compute for the 305 MAGs that run 4 added (the 11 re-assembled libraries; 13,164 MAGs in narrative 272312).
# Same binaries, databases and parameters as run_census_compute.sh / run_census_rescan.sh:
#  1. CheckM2 1.1.0 (DB uniref100.KO.1) --allmodels on the 305 new MAGs            -> checkm2_new305/quality_report.tsv
#  2. skani 0.3.1 triangle over ALL 13,164 MAGs (-E -s 90; 78 s for 12,859)        -> skani_triangle.tsv
#     (deterministic: the 12,859 x 12,859 part must reproduce the census triangle; checked by the local step)
#  3. skani 0.3.1 search of the 305 vs novel_all_representatives, novel_hq_v2, gtdb_r220_consolidated with
#     --no-marker-index and no -n cap (= the brute-force census re-scan)          -> skani_vs_<db>.new305.tsv
# Writes only under $W4. Resumable (.<step>.done).
# Usage (poplar): cd /scratch1/afreiburger/wetland_census && setsid nohup bash run_census_compute_run4.sh > run4.log 2>&1 < /dev/null &
set -uo pipefail
W=/scratch1/afreiburger/wetland_census
W4=$W/run4
mkdir -p "$W4/new305"
cd "$W4"
MAGS=/scratch1/afreiburger/wetland_coverm_pilot/all_mags
QC=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/8c703f8abdc13cdc1ef374a2bcf37f14_
CK2DB=/scratch1/afreiburger/emilykin/processed/refdata/checkm2/CheckM2_database/uniref100.KO.1.dmnd
SKANI=/scratch1/chenry/opt/skani/0.3.1/skani
declare -A DB=( [novel_all_representatives]=/scratch1/chenry/kbdl/skani_dbs/ad1c32979e1b4cc897c5f32631e2cc9e
                [novel_hq_v2]=/scratch1/chenry/kbdl/skani_dbs/c63e6b662bbf402e901974f84cd29d31
                [gtdb_r220_consolidated]=/scratch1/fliu/hub_scratch/fliu/data/ani/skani/gtdb_r220 )
ev() { echo "[$(date '+%F %T')] $*" >> census_run4.events; }
tm() { local tag=$1; shift; /usr/bin/time -v -o "time_${tag}.txt" "$@"; local rc=$?; \
       ev "TIME $tag rc=$rc wall=$(grep 'Elapsed (wall' time_${tag}.txt | awk '{print $NF}')"; return $rc; }
ev "START load $(cut -d' ' -f1-3 /proc/loadavg)"
[ -s new305.txt ] || { ev "FAIL new305.txt missing"; exit 1; }
[ "$(grep -c . new305.txt)" -eq 305 ] || { ev "FAIL new305.txt is not 305 MAGs"; exit 1; }
ls "$MAGS"/*.fna > mags.list
n=$(wc -l < mags.list); [ "$n" -eq 13164 ] || { ev "FAIL expected 13164 MAGs, found $n"; exit 1; }
sed "s#^#$MAGS/#; s#\$#.fna#" new305.txt > new305.list
while read -r m; do ln -sf "$MAGS/$m.fna" "new305/$m.fna"; done < new305.txt
md5sum "$SKANI" > skani_binary.md5
[ "$(cut -d' ' -f1 skani_binary.md5)" = "$(cut -d' ' -f1 "$W/skani_binary.md5")" ] \
  || { ev "FAIL the skani binary differs from the one the census used"; exit 1; }

if [ ! -f .checkm2.done ]; then
  ( mkdir -p /scratch/afreiburger/tmp_checkm2_wetland_run4
    tm checkm2 env PATH=$QC/bin:$PATH TMPDIR=/scratch/afreiburger/tmp_checkm2_wetland_run4 checkm2 predict -i new305 -x fna \
        -t 32 --allmodels --database_path "$CK2DB" --tmpdir /scratch/afreiburger/tmp_checkm2_wetland_run4 \
        -o checkm2_new305 --force \
      && [ "$(tail -n +2 checkm2_new305/quality_report.tsv | wc -l)" -eq 305 ] && touch .checkm2.done ) &
  CK_PID=$!
fi
if [ ! -f .triangle.done ]; then
  tm triangle "$SKANI" triangle -l mags.list -t 16 -E -s 90 -o skani_triangle.tsv && touch .triangle.done
fi
for name in novel_all_representatives novel_hq_v2 gtdb_r220_consolidated; do
  [ -f ".search_${name}.done" ] && continue
  tm "search_${name}" "$SKANI" search --ql new305.list -d "${DB[$name]}" -t 16 --no-marker-index \
      -o "skani_vs_${name}.new305.tsv" && touch ".search_${name}.done"
done
[ -n "${CK_PID:-}" ] && wait "$CK_PID"
ev "DONE checkm2=$([ -f .checkm2.done ] && echo ok || echo FAIL) triangle=$([ -f .triangle.done ] && echo ok || echo FAIL) searches=$(ls .search_*.done 2>/dev/null | wc -l)/3"
