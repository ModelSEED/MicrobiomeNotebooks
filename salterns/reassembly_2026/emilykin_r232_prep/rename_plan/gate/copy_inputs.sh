#!/bin/bash
# copy untracked input files (not symlinks) from the main repo into the HEAD worktree; log each one.
# usage: copy_inputs.sh REASON RELPATH...
set -euo pipefail
RP=/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan
WT=$RP/head_worktree; MAIN=$HOME/Documents/EmilyKin; LST=$RP/gate_scratch/copied_inputs.tsv
reason=$1; shift
[ -f "$LST" ] || printf "# relpath\tbytes\tmd5\treason\n" > "$LST"
for rel in "$@"; do
  src="$MAIN/$rel"
  if [ -d "$src" ]; then
    (cd "$MAIN" && find "$rel" -type f) | while read -r f; do "$0" "$reason" "$f"; done; continue
  fi
  [ -f "$src" ] || { echo "ABSENT in main repo: $rel"; continue; }
  if (cd "$WT" && git ls-files --error-unmatch "$rel" >/dev/null 2>&1); then echo "tracked, not copied: $rel"; continue; fi
  grep -qF -- "$rel	" "$LST" && continue
  mkdir -p "$(dirname "$WT/$rel")"; cp -p "$src" "$WT/$rel"
  printf "%s\t%s\t%s\t%s\n" "$rel" "$(stat -c %s "$src")" "$(md5sum < "$src" | cut -c1-32)" "$reason" >> "$LST"
done
