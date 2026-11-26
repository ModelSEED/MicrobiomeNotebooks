#!/usr/bin/env bash
# Stage the 291 EmilyKin MAG FASTAs for the R232 run and write the GTDB-Tk batchfile.
# Default source = the committed best-of-both catalogue in the local EmilyKin repo (git 6541627; the same sequences as the
# KBase 262319 assemblies). Its md5s were frozen in ../local_catalogue_md5.txt when this plan was written: if any file
# changed since, this script stops (re-run the projection, or pass ACCEPT_CHANGED=1 after checking why).
# SRC=poplar copies poplar's /scratch1/afreiburger/emilykin/EmilyKin/dereplicated_genomes instead (the directory the R220
# gtdbtk_291 run used, per reprocess_291/env.sh MAGDIR); run 00_poplar_inventory.sh first and only use it if the
# catalogue comparison says 0 files differ, or if you deliberately want to re-classify poplar's version.
# Usage: [SRC=local|poplar] [WORK=~/gtdbtk_kbase/emilykin_r232] bash 10_stage_genomes.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
SRC=${SRC:-local}
WORK=${WORK:-$HOME/gtdbtk_kbase/emilykin_r232}
G="$WORK/genomes"
mkdir -p "$G"
if [ "$SRC" = local ]; then
  LOCAL=$HOME/Documents/EmilyKin/dereplicated_genomes
  ( cd "$LOCAL" && md5sum -c --quiet "$HERE/../local_catalogue_md5.txt" ) || {
    [ "${ACCEPT_CHANGED:-0}" = 1 ] || { echo "local catalogue differs from local_catalogue_md5.txt - stopping"; exit 1; }; }
  cp -p "$LOCAL"/*.fa "$G"/
  ( cd "$G" && md5sum -c --quiet "$HERE/../local_catalogue_md5.txt" ) || [ "${ACCEPT_CHANGED:-0}" = 1 ] || {
    echo "staged copy does not match local_catalogue_md5.txt - stopping"; exit 1; }
elif [ "$SRC" = poplar ]; then
  ssh -o BatchMode=yes -o ClearAllForwardings=yes poplar true || { echo "poplar not reachable"; exit 2; }
  rsync -a --include='*.fa' --exclude='*' -e "ssh -o BatchMode=yes -o ClearAllForwardings=yes" \
      poplar:/scratch1/afreiburger/emilykin/EmilyKin/dereplicated_genomes/ "$G"/
  ( cd "$G" && md5sum *.fa ) > "$WORK/staged_md5.txt"
  if [ -s "$HERE/../poplar_inventory/poplar_catalogue_md5.txt" ]; then
    if diff -q <(sort -k2 "$WORK/staged_md5.txt") <(sort -k2 "$HERE/../poplar_inventory/poplar_catalogue_md5.txt") > /dev/null
    then echo "md5 matches the poplar inventory"; else echo "WARNING: staged files differ from the poplar inventory md5s"; fi
  fi
  if ! ( cd "$G" && md5sum -c --quiet "$HERE/../local_catalogue_md5.txt" ) > /dev/null 2>&1; then
    echo "NOTE: poplar's catalogue differs from the local best-of-both one (local_catalogue_md5.txt)"
  fi
else
  echo "SRC must be local or poplar"; exit 1
fi
n=$(ls "$G"/*.fa | wc -l)
[ "$n" -eq 291 ] || { echo "expected 291 MAGs, found $n"; exit 1; }
for f in "$G"/*.fa; do printf '%s\t%s\n' "$f" "$(basename "$f" .fa)"; done > "$WORK/batchfile.tsv"
( cd "$G" && md5sum *.fa ) > "$WORK/staged_md5.txt"
echo "staged $n MAGs from $SRC -> $WORK/batchfile.tsv"
