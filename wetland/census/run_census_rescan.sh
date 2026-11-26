#!/usr/bin/env bash
# Re-run the census's three reference searches (run_census_compute.sh steps 3-4) WITHOUT skani's marker index.
#
# Why: skani 0.3.1 `search` builds an approximate marker-hash index when given > 100 query files (default; see
# `--no-marker-index`), and that index loses hits, more the more queries a call has. Measured 2026-11-17 on
# novel_all_representatives:
#   - the census's single 12,859-query call has 57 rows (-n 5) for 50 test MAGs, while 50- and 101-query calls give 68;
#   - bin.001.fasta__SRR10968274 gets no row at all in the census, but 90.75% / AF 41.7 to spire_mag_00720575 alone;
#   - a 2,000-query call gives 2,476 rows vs 2,372 for the same MAGs inside the census call.
# A re-run of the exact census command reproduced the census output (14,039 rows, same 6,511 queries), so the loss is
# deterministic, not noise. The skani triangle (clustering) is NOT affected: a 300-MAG subset gave the same 684 edges
# with identical ANI.
#
# This run: same skani 0.3.1 binary, DBs, MAG list and screen as the census, plus --no-marker-index, and no -n cap
# (so a species-level hit with AF >= 50 cannot be pushed out of a top-5 by higher-ANI low-AF hits).
# Output: $W/rescan_noindex/skani_vs_<db>.tsv (+ time_*.txt, .<db>.done); the census originals are left untouched.
# Usage (poplar): cd /scratch1/afreiburger/wetland_census && setsid nohup bash run_census_rescan.sh > rescan.log 2>&1 < /dev/null &
set -uo pipefail
W=/scratch1/afreiburger/wetland_census
O=$W/rescan_noindex
mkdir -p "$O"
cd "$W"
SKANI=/scratch1/chenry/opt/skani/0.3.1/skani
declare -A DB=( [novel_all_representatives]=/scratch1/chenry/kbdl/skani_dbs/ad1c32979e1b4cc897c5f32631e2cc9e
                [novel_hq_v2]=/scratch1/chenry/kbdl/skani_dbs/c63e6b662bbf402e901974f84cd29d31
                [gtdb_r220_consolidated]=/scratch1/fliu/hub_scratch/fliu/data/ani/skani/gtdb_r220 )
ev() { echo "[$(date '+%F %T')] $*" >> "$O/rescan.events"; }
run() {
  local name=$1
  [ -f "$O/.${name}.done" ] && return 0
  ev "START $name"
  /usr/bin/time -v -o "$O/time_${name}.txt" "$SKANI" search --ql mags.list -d "${DB[$name]}" -t 12 --no-marker-index \
      -o "$O/skani_vs_${name}.tsv" 2> "$O/${name}.err" && touch "$O/.${name}.done" && ev "DONE $name" || ev "FAIL $name"
}
for name in novel_all_representatives novel_hq_v2 gtdb_r220_consolidated; do run "$name" & done
wait
ev "ALL $(ls "$O"/.*.done 2>/dev/null | wc -l)/3 done"
