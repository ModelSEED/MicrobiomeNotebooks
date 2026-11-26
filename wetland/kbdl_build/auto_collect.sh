#!/usr/bin/env bash
# Run the post-completion steps automatically (per Andrew, 2026-11-18: "automatically run this after the job completes").
# Local, detached; polls poplar every 10 min; each task runs once, when its job reports completion:
#  1. KBDLSKANI sweep (skani_events.log "DONE all archives" or "STOPPED"): collect_skani_sweep.py -> kbdl_build/skani_sweep/
#     (merged hits, failed MAGs, pair-level check vs the brute-force census); commit + push. A STOPPED sweep is
#     collected too and the commit says PARTIAL.
#  2. DRAM2 residual bundle (dram2_residual/bundle.events "DONE"): rsync the tarballs to dram2/results/residual/
#     (gitignored), verify SHA256SUMS, write dram2/results/RESIDUAL_MANIFEST.tsv (tarball, bytes, sha256, member
#     files); commit the manifest + the non-coding tables + the per-MAG gene summary; push.
# Gives up after 4 days. Log: ~/gtdbtk_kbase/auto_collect.log
# Usage: setsid nohup bash kbdl_build/auto_collect.sh > /dev/null 2>&1 < /dev/null &
set -uo pipefail
REPO=$HOME/Documents/MicrobiomeNotebooks
WET=$REPO/wetland
PY="$HOME/Documents/py_venv/bin/python -I"
LOG=$HOME/gtdbtk_kbase/auto_collect.log
W=/scratch/afreiburger/wetland_kbdl_build
B=/scratch/afreiburger/dram2_residual
BR=clauded-insights-figure7
log() { echo "[$(date '+%F %T')] $*" >> "$LOG"; }
rsh() { timeout 120 ssh -o BatchMode=yes poplar "$@" < /dev/null; }
push() {
  for i in 1 2 3 4; do git -C "$REPO" push -q -u origin "$BR" >> "$LOG" 2>&1 && return 0; sleep $((2 ** i)); done
  log "WARN push failed 4 times"; return 1
}
commit() {  # commit <message> <paths...>
  local msg=$1; shift
  for i in 1 2 3 4 5; do
    [ -e "$REPO/.git/index.lock" ] && { sleep 20; continue; }
    git -C "$REPO" add -- "$@" >> "$LOG" 2>&1 || return 1
    git -C "$REPO" diff --cached --quiet && { log "nothing to commit for: ${msg%%$'\n'*}"; return 0; }
    git -C "$REPO" commit -q -m "$msg" -- "$@" >> "$LOG" 2>&1 && { push; return $?; }
    sleep 20
  done
  log "WARN commit failed: ${msg%%$'\n'*}"; return 1
}
TRAILER=$'\n\nCo-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>'
sweep_done=0; dram_done=0; deadline=$(( $(date +%s) + 4 * 86400 ))
log "START auto_collect pid $$"
while [ $sweep_done = 0 ] || [ $dram_done = 0 ]; do
  [ "$(date +%s)" -gt "$deadline" ] && { log "GIVE UP after 4 days (sweep_done=$sweep_done dram_done=$dram_done)"; exit 1; }
  if [ $sweep_done = 0 ]; then
    st=$(rsh "grep -E 'DONE all archives|STOPPED' $W/skani_events.log | tail -1")
    if [ -n "$st" ]; then
      log "sweep: $st"
      if (cd "$WET/kbdl_build" && $PY collect_skani_sweep.py > "$WET/kbdl_build/skani_sweep/collect.log" 2>&1); then
        s=$($PY -c "import json; s=json.load(open('$WET/kbdl_build/skani_sweep/summary.json')); p=s['novel_pairs_vs_bruteforce_census']; print(f\"{s['mags_swept']}/{s['of_census_mags']} MAGs swept, {s['failed_mags']} failed; novel pairs vs brute-force census: {p['kbdl_pairs']} vs {p['census_pairs']}, only KBDL {p['only_kbdl']}, only census {p['only_census']}\")")
        kind=$(echo "$st" | grep -q STOPPED && echo "PARTIAL (sweep halted)" || echo "complete")
        log "collected: $s"
        commit "Wetland KBDL skani sweep $kind: R232 + novel_all_representatives hits for all archives"$'\n\n'"$s"$'\n'"(collected automatically by kbdl_build/auto_collect.sh when the sweep reported: ${st#*] })$TRAILER" \
          "wetland/kbdl_build/skani_sweep" && sweep_done=1
      else
        log "WARN collect_skani_sweep.py failed; see skani_sweep/collect.log; retry next poll"
      fi
    fi
  fi
  if [ $dram_done = 0 ] && rsh "grep -q '] DONE' $B/bundle.events"; then
    D=$WET/dram2/results/residual
    mkdir -p "$D"
    if rsync -a "poplar:$B/" "$D/" >> "$LOG" 2>&1 && (cd "$D" && sha256sum -c --quiet SHA256SUMS >> "$LOG" 2>&1); then
      M=$WET/dram2/results/RESIDUAL_MANIFEST.tsv
      { printf "tarball\tbytes\tsha256\tmember_files\tsource_dir_on_poplar\n"
        for f in "$D"/*.tar.gz; do
          n=$(tar -tzf "$f" | grep -vc '/$'); src=$(tar -tzf "$f" | head -1 | cut -d/ -f1-4)
          printf "%s\t%s\t%s\t%s\t%s\n" "$(basename "$f")" "$(stat -c %s "$f")" "$(grep " $(basename "$f")\$" "$D/SHA256SUMS" | cut -d' ' -f1)" "$n" "/scratch1/afreiburger/wetland_dram2/$src"
        done; } > "$M"
      log "DRAM2 residual: $(ls "$D"/*.tar.gz | wc -l) tarballs verified, $(du -sh "$D" | cut -f1)"
      commit "Wetland DRAM2: residual outputs pulled (gene calls, per-DB hits, raw KOfam hits, tRNA/rRNA, QUAST, logs) + manifest"$'\n\n'"Tarballs in dram2/results/residual/ (gitignored, checksums verified); RESIDUAL_MANIFEST.tsv lists each with size, sha256 and file count. Pulled automatically by kbdl_build/auto_collect.sh.$TRAILER" \
        "wetland/dram2/results/RESIDUAL_MANIFEST.tsv" && dram_done=1
    else
      log "WARN DRAM2 residual rsync/checksum failed; retry next poll"
    fi
  fi
  [ $sweep_done = 1 ] && [ $dram_done = 1 ] && break
  sleep 600
done
log "DONE auto_collect"
