#!/usr/bin/env bash
# Keeps a running DRAM2 chunk's local work dir small (/scratch is shared and ~90 % full).
# Touches ONLY tasks that have finished with exit code 0:
#   * HMM_SEARCH_* tasks: hmmsearch prints full alignments to stdout (the module has no -o), which
#     Nextflow keeps twice (.command.out + .command.log, up to ~400 MB each for KOfam). Truncated to 0.
#     Not used downstream (the --domtblout file is the declared output).
#   * MMSEQS_INDEX tasks: `mmseqs createindex` writes a ~0.9 GB <mag>.mmsdb.idx per MAG that the
#     searches never read (MMseqs2 only uses the TARGET db's index). Deleted once every
#     MMSEQS_SEARCH_* task of that MAG has finished OK.
# Safety valve: if free space on the work filesystem drops below LOW_DISK_GB, writes
# runs/<run>/LOWDISK and sends SIGTERM to that run's Nextflow (resumable later with -resume).
#
#   janitor.sh <work_dir> <run_dir> [once]
set -uo pipefail
source /scratch1/afreiburger/wetland_dram2/bin/dram2_settings.sh
W=$1; R=$2; ONCE=${3:-}
INTERVAL=${JANITOR_INTERVAL:-240}
LOW_DISK_GB=${LOW_DISK_GB:-120}
WARN_DISK_GB=${WARN_DISK_GB:-250}
STATE=$R/janitor.state; LOG=$R/janitor.log
n_mm=0; for db in $DRAM2_DBS; do case $db in merops|camper|canthyd|methyl|pfam|viral|kegg|uniref) n_mm=$((n_mm+1));; esac; done
touch "$STATE"
declare -A seen
while IFS=$'\t' read -r d k m; do seen[$d]=1; done < "$STATE"

while true; do
  nh=0
  for d in "$W"/??/*; do
    [ -n "${seen[$d]+x}" ] && continue
    [ -s "$d/.exitcode" ] || continue
    if [ "$(<"$d/.exitcode")" != 0 ]; then seen[$d]=1; printf '%s\tfailed\t-\n' "$d" >> "$STATE"; continue; fi
    kind=other; mag=-
    if grep -q '^ *hmmsearch' "$d/.command.sh" 2>/dev/null; then
      kind=hmm; : > "$d/.command.out"; : > "$d/.command.log"; nh=$((nh+1))
    elif l=$(grep -m1 -o 'mmseqs createdb [^ ]* [^ ]*\.mmsdb' "$d/.command.sh" 2>/dev/null); then
      kind=index; mag=$(basename "${l##* }" .mmsdb)
    elif l=$(grep -m1 -o 'mmseqs search query_database/[^ ]*\.mmsdb' "$d/.command.sh" 2>/dev/null); then
      kind=search; mag=$(basename "${l##*/}" .mmsdb)
    fi
    seen[$d]=1; printf '%s\t%s\t%s\n' "$d" "$kind" "$mag" >> "$STATE"
  done
  ni=0
  if [ "$n_mm" -gt 0 ]; then
    while read -r d; do
      if compgen -G "$d/*.mmsdb.idx*" >/dev/null; then rm -f "$d"/*.mmsdb.idx*; ni=$((ni+1)); fi
    done < <(awk -F'\t' -v n="$n_mm" '$2=="search"{c[$3]++} $2=="index"{dir[$3]=$1} END{for(m in dir) if(c[m]>=n) print dir[m]}' "$STATE")
  fi
  free=$(df -BG --output=avail "$W" | tail -1 | tr -dc 0-9)
  used=$(find "$W" -type f -printf '%s\n' 2>/dev/null | awk '{s+=$1} END{printf "%.1f", s/1e9}')
  echo "[$(date -Is)] trimmed hmm=$nh idx_dirs=$ni work=${used}GB free=${free}GB" >> "$LOG"
  if [ "$free" -lt "$LOW_DISK_GB" ]; then
    echo "[$(date -Is)] free ${free}GB < ${LOW_DISK_GB}GB: stopping Nextflow for $W" >> "$LOG"
    echo "$(date -Is) free=${free}GB" > "$R/LOWDISK"
    pkill -TERM -u "$(id -un)" -f -- "-work-dir $W " ; exit 3
  elif [ "$free" -lt "$WARN_DISK_GB" ]; then
    echo "[$(date -Is)] WARNING free ${free}GB" >> "$LOG"
  fi
  [ -n "$ONCE" ] && exit 0
  sleep "$INTERVAL"
done
