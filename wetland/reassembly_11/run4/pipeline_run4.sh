#!/usr/bin/env bash
# Run 4 (per Andrew, 2026-11-17): put the MAGs of the 11 re-assembled libraries (binned by
# ~/Documents/MicrobiomeNotebooks/wetland/reassembly_11/binning_driver.py) through
# the same pipeline as the other 12,611 MAGs of KBase narrative 272312:
#   extraction (narrative cell B v36, POOL_CUTOFF 2026-11-05) -> local GTDB-Tk 2.8.0 / GTDB R232 (+ 60 run-2 controls)
#   -> write R232 onto the MAGs + report -> cell B makes MAG_all.part128+ sets and submits CheckM (stage 2b)
#   -> CheckM pushed onto the MAGs + report.
# Stages are resumable (.run4_<stage>.done); milestones -> pipeline_run4.events, details -> pipeline_run4.log.
# Usage: cd ~/gtdbtk_kbase && setsid nohup bash pipeline_run4.sh > /dev/null 2>&1 < /dev/null &
set -uo pipefail
cd "$(dirname "$0")"
EV=pipeline_run4.events
LOG=pipeline_run4.log
PY="$HOME/Documents/py_venv/bin/python -I -B"
UNIT=gtdbtk-wetland-run4
ev()    { local l="[$(date '+%F %T')] $*"; echo "$l" >> "$EV"; echo "$l" >> "$LOG"; }
die()   { ev "FAIL $*"; exit 1; }
done_() { [ -f ".run4_$1.done" ]; }
mark()  { touch ".run4_$1.done"; }
advance() {   # narrative cells A/B, unchanged (v36), via the local harness; output appended to advance_run4.log
  (cd narrative_run && $HOME/Documents/py_venv/bin/python -B run_cells.py --mode advance) > narrative_run/advance_run4_last.log 2>&1
  local rc=$?; cat narrative_run/advance_run4_last.log >> narrative_run/advance_run4.log; return $rc
}
export KB_AUTH_TOKEN="$(cat "$HOME/.kbase/token")"          # never printed
export WB_LINEAGES="$PWD/mag_gtdb_lineages_run4_merged.tsv" WB_MANIFEST="$PWD/all_manifest_run4_merged.tsv"
export WB_OUT="$PWD/merged_run4/gtdbtk_run4" WB_BATCH="$PWD/batchfile_run4.tsv"
ev "START pipeline_run4 (pid $$)"

# ---- 0. wait for binning (MetaBAT2/MaxBin2/CONCOCT -> DAS Tool), then move cell A POOL_CUTOFF past the new DAS Tool runs
BINEV=/home/freiburger/Documents/MicrobiomeNotebooks/wetland/reassembly_11/binning.events
if ! done_ cutoff; then
  until grep -q ALLDONE "$BINEV" 2>/dev/null; do sleep 1800; done
  ev "BINNING $(grep ALLDONE "$BINEV" | tail -1 | cut -c1-200)"
  CUT=$(date -u -d '+1 hour' '+%Y-%m-%dT%H:00:00+0000')
  ver=$(KB_AUTH_TOKEN="$KB_AUTH_TOKEN" $PY -c "import sys; sys.path.insert(0,'/home/freiburger/gtdbtk_kbase'); import kbase_mags as k; print(k.rpc(k.WS,'Workspace.get_object_info3',[{'objects':[{'ref':'272312/1'}]}])[0]['infos'][0][4])")
  (cd narrative_run && $PY save_cells_cutoff.py --expect-version "$ver" --cutoff "$CUT" --apply) >> "$LOG" 2>&1 || die "narrative cutoff edit (see $LOG)"
  $PY - "$CUT" >> "$LOG" 2>&1 <<'PYEOF' || die "pool cache seed"
import json, sys
p = "/home/freiburger/gtdbtk_kbase/narrative_run/mag_annotation_pool.json"
c = json.load(open(p)); c["thresholds"][0] = sys.argv[1]; json.dump(c, open(p, "w")); print("pool cache thresholds ->", c["thresholds"])
PYEOF
  ev "CUTOFF cell A POOL_CUTOFF -> $CUT (narrative v$ver -> v$((ver+1)))"
  mark cutoff
fi

# ---- 1. extraction: re-run cell B every 30 min (it resubmits failed extractions, up to 5 tries) until stage 1 is
#         done and the new MAGs wait for their GTDB lineage ---------------------------------------------------------
if ! done_ extracted; then
  deadline=$(date -d "2026-12-10" +%s)
  while :; do
    advance || ev "WARN cell B exit $? (narrative_run/advance_run4_last.log)"
    s1=$(grep -E '^1 extract' narrative_run/advance_run4_last.log | tail -1)
    w=$(grep -oE '[0-9]+ newly extracted MAGs wait' narrative_run/advance_run4_last.log | tail -1)
    ev "STATUS ${s1:-no stage-1 line} | ${w:-}"
    grep -q 'gave up' narrative_run/advance_run4_last.log && die "an extraction job failed 5 times (advance_run4_last.log)"
    [ -n "$w" ] && break
    [ "$(date +%s)" -gt "$deadline" ] && die "extraction not done by the deadline"
    sleep 1800
  done
  ev "EXTRACTED $w"
  mark extracted
fi

# ---- 2. list + select (new MAGs + controls) + fetch -----------------------------------------------------------------
if ! done_ fetched; then
  $PY make_all_list.py all_mags_run4_full.tsv >> "$LOG" 2>&1 || die "make_all_list.py"
  $PY run4_tools.py select >> "$LOG" 2>&1 || die "run4_tools.py select"
  ev "SELECTED $(grep '^select:' "$LOG" | tail -1)"
  $PY -u kbase_mags_kbaseus.py fetch --list all_mags_run4.tsv --genome-dir genomes --manifest all_manifest_run4.tsv \
      > fetch_run4.log 2>&1 || die "fetch (fetch_run4.log)"
  $PY - >> "$LOG" 2>&1 <<'EOF' || die "manifest check (see pipeline_run4.log)"
import csv, sys
rows = list(csv.DictReader(open("all_manifest_run4.tsv"), delimiter="\t"))
want = {r["name"] for r in csv.DictReader(open("all_mags_run4.tsv"), delimiter="\t")}
names = [r["name"] for r in rows]
if len(set(names)) != len(names) or set(names) != want:
    sys.exit(f"manifest has {len(names)} rows / {len(set(names))} unique; list has {len(want)}")
bad = []
for r in rows:
    n = bp = 0
    for line in open(r["file"]):
        if line.startswith(">"):
            n += 1
        else:
            bp += len(line.strip())
    if str(n) != r["contigs"] or str(bp) != r["size_bp"]:
        bad.append(r["name"])
if bad:
    sys.exit(f"{len(bad)} FASTAs differ from KBase contig count/bp: {bad[:5]}")
print(f"manifest OK: {len(rows)} MAGs, every FASTA matches KBase's contig count and bp")
EOF
  ev "FETCHED $(grep 'manifest OK' "$LOG" | tail -1)"
  mark fetched
fi

# ---- 3. GTDB-Tk (memory-capped scope + guard; waits for memory and for other sessions' GTDB-Tk runs) -------------
if ! done_ classified; then
  waited=0
  while [ "$(awk '/MemAvailable/ {printf "%d", $2/1048576}' /proc/meminfo)" -lt 150 ] || pgrep -f '[g]tdbtk classify_wf' > /dev/null; do
    [ $(( waited % 1800 )) -eq 0 ] && ev "WAIT MemAvailable $(awk '/MemAvailable/ {printf "%d", $2/1048576}' /proc/meminfo) GB (need 150) / other GTDB-Tk runs: $(pgrep -fc '[g]tdbtk classify_wf')"
    sleep 300; waited=$(( waited + 300 ))
  done
  pgrep -f 'bash mem_guard[.]sh' > /dev/null || UNIT=$UNIT.scope setsid nohup bash mem_guard.sh > /dev/null 2>&1 < /dev/null &
  rm -rf tmp gtdbtk_run4 && mkdir -p tmp
  ev "LAUNCH GTDB-Tk classify_wf on $(( $(wc -l < all_manifest_run4.tsv) - 1 )) MAGs in scope $UNIT (MemoryMax=170G, 48 CPUs, 8 pplacer threads)"
  systemd-run --user --scope --unit=$UNIT -p MemoryMax=170G -p MemorySwapMax=0 bash run_gtdbtk_run4.sh >> gtdbtk_run4.log 2>&1; rc=$?
  [ $rc -eq 0 ] && [ -f gtdbtk_run4/gtdbtk.json ] || die "GTDB-Tk exit $rc (gtdbtk_run4.log, gtdbtk_run4/gtdbtk.log, mem_guard.events)"
  ev "CLASSIFIED $(tr '\r' '\n' < gtdbtk_run4/gtdbtk.log | grep -E 'genomes have a warning|classified using the ANI' | sed 's/ \{2,\}/ /g' | tr '\n' ' ')"
  mark classified
fi

# ---- 4. collect + gate + merge ----------------------------------------------------------------------------------
if ! done_ gated; then
  $PY kbase_mags.py collect --manifest all_manifest_run4.tsv --out-dir gtdbtk_run4 --out mag_gtdb_lineages_run4.tsv \
      > collect_run4.log 2>&1 || die "collect (collect_run4.log)"
  $PY run4_tools.py gate > check_run4.log 2>&1; rc=$?
  ev "GATE $(tr '\n' ' ' < check_run4.log | cut -c1-600)"
  [ $rc -eq 0 ] || die "gate failed (check_run4.log); nothing written to KBase"
  $PY run4_tools.py merge >> "$LOG" 2>&1 || die "merge"
  ev "MERGED $(grep '^merge:' "$LOG" | tail -1)"
  mark gated
fi

# ---- 5. R232 onto the new MAGs -----------------------------------------------------------------------------------
if ! done_ applied; then
  nnew=$(grep -oE 'new MAGs [0-9]+, [0-9]+ with a lineage' check_run4.log | grep -oE '[0-9]+ with' | grep -oE '[0-9]+')
  $PY writeback_r232.py apply > wb_run4_plan.log 2>&1 || die "writeback dry run (wb_run4_plan.log)"
  ev "PLAN $(grep -E 'to write|^HELD' wb_run4_plan.log | cut -c1-300 | tr '\n' ' ')"
  planned=$(grep -oE '[0-9]+ to write' wb_run4_plan.log | head -1 | grep -oE '[0-9]+')
  [ "$planned" = "$nnew" ] || die "plan writes ${planned:-?} MAGs, expected the $nnew new ones (wb_run4_plan.log)"
  $PY -u writeback_r232.py apply --apply > wb_run4_apply.log 2>&1 || die "writeback apply (wb_run4_apply.log)"
  ev "APPLIED $(grep -E '^Done:' wb_run4_apply.log)"
  mark applied
fi
if ! done_ pushed; then
  $PY writeback_r232.py verify > verify_run4.log 2>&1 || die "verify (verify_run4.log)"
  grep -qE 'NOT written|DIFFER:|WRONG:' verify_run4.log && die "verify found problems (verify_run4.log)"
  ev "VERIFIED $(grep -E 'carry exactly|sequence fields|std_lineages.gtdb holds' verify_run4.log | cut -c1-160 | tr '\n' ' ')"
  # no 'sets': MAG_pool sets are kb_gtdbtk inputs and keep pointing at their analysed versions (--allow-partial)
  $PY writeback_r232.py report --apply --allow-partial > report_run4.log 2>&1 || die "report (report_run4.log)"
  ev "REPORT $(grep -E '^Saved report|^Read back' report_run4.log | tr '\n' ' ')"
  mark pushed
fi

# ---- 6. cell B: MAG_all.part128+ sets + CheckM (stage 2b); wait for CheckM ----------------------------------------
if ! done_ checkm; then
  deadline=$(date -d "2026-12-20" +%s)
  while :; do
    advance || ev "WARN cell B exit $? (narrative_run/advance_run4_last.log)"
    s2b=$(grep -E '^2b ' narrative_run/advance_run4_last.log | tail -1)
    ev "STATUS $(grep -E '^  saved MAG_all' narrative_run/advance_run4_last.log | tr '\n' ' ')${s2b}"
    echo "$s2b" | grep -qE 'gave up' && die "a CheckM job failed 5 times (advance_run4_last.log)"
    read -r tasks done_n < <($HOME/Documents/py_venv/bin/python -I -c '
import sys
L = open(sys.argv[1]).read().splitlines()
hdr = next((l.split("\t") for l in L if l.startswith("stage\t")), None)
row = next((l.split("\t") for l in L if l.startswith("2b ")), None)
d = dict(zip(hdr, row)) if hdr and row else {}
print(d.get("tasks") or 0, d.get("done") or 0)' narrative_run/advance_run4_last.log)
    [ -n "$tasks" ] && [ "$tasks" -gt 130 ] && [ "$tasks" = "$done_n" ] && break
    [ "$(date +%s)" -gt "$deadline" ] && die "CheckM not done by the deadline"
    sleep 1800
  done
  ev "CHECKM DONE $s2b"
  mark checkm
fi

# ---- 7. CheckM onto the new MAGs ------------------------------------------------------------------------------------
if ! done_ checkm_pushed; then
  $PY narrative_run/collect_checkm.py > narrative_run/collect_checkm_run4.log 2>&1 || die "collect_checkm (collect_checkm_run4.log)"
  ev "CHECKM TABLE $(tail -2 narrative_run/collect_checkm_run4.log | tr '\n' ' ')"
  (cd checkm_push && $PY writeback_checkm.py apply) > checkm_push/plan_run4.log 2>&1 || die "CheckM plan (checkm_push/plan_run4.log)"
  ev "CHECKM PLAN $(grep -E 'to write' checkm_push/plan_run4.log | cut -c1-300)"
  (cd checkm_push && $PY -u writeback_checkm.py apply --apply) > checkm_push/apply_run4.log 2>&1 || die "CheckM apply (checkm_push/apply_run4.log)"
  (cd checkm_push && $PY writeback_checkm.py verify) > checkm_push/verify_run4.log 2>&1 || die "CheckM verify (checkm_push/verify_run4.log)"
  ev "CHECKM VERIFIED $(tail -4 checkm_push/verify_run4.log | tr '\n' ' ' | cut -c1-400)"
  (cd checkm_push && $PY writeback_checkm.py report --apply) > checkm_push/report_run4.log 2>&1 || die "CheckM report (checkm_push/report_run4.log)"
  ev "CHECKM REPORT $(grep -E '^Saved report|^Read back' checkm_push/report_run4.log | tr '\n' ' ')"
  mark checkm_pushed
fi
ev "DONE run 3: the new MAGs carry GTDB R232 + CheckM like the other MAGs of narrative 272312"
