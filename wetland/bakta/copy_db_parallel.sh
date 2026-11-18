#!/usr/bin/env bash
# Copy the Bakta 6.0-full DB from the NFS share to poplar's local disk with many parallel read streams (one rsync
# stream reads ~9-13 MB/s from the share; the DB is 83 GB, psc.dmnd 46 GB and bakta.db 42 GB). Files > 2 GB are split
# into 8 byte ranges copied concurrently (dd seek/skip, conv=notrunc); smaller files are copied 16 at a time. Every file
# keeps the source mtime/mode; files > 2 GB are then compared byte for byte (cmp) against the source.
# The root-only amrfinderplus-db/2024-12-18.1 is unreadable and not needed (Bakta reads amrfinderplus-db/latest ->
# 2025-06-03.1); it is skipped.
# Exit 0 only if every readable file is present with the source size and every large file is byte-identical.
set -uo pipefail
SRC=${1:-/scratch1/fliu/hub_scratch/data/bakta/db}
DST=${2:-/scratch/afreiburger/bakta_db_6.0_full}
mkdir -p "$DST"
cd "$SRC" || exit 1
find . -type d -readable -printf '%P\n' 2>/dev/null | while read -r d; do mkdir -p "$DST/$d"; done
[ -L "$DST/amrfinderplus-db/latest" ] || cp -P amrfinderplus-db/latest "$DST/amrfinderplus-db/"
find . -type f -readable -printf '%s %P\n' 2>/dev/null | sort -rn > "$DST/.filelist"
seg() {  # seg <file> <i> <n> <size>
  local f=$1 i=$2 n=$3 size=$4 bs=$((16 * 1024 * 1024))
  local blocks=$(( (size + bs - 1) / bs )) per=$(( (size + bs - 1) / bs / n + 1 ))
  local start=$(( i * per ))
  [ $start -ge $blocks ] && return 0
  dd if="$SRC/$f" of="$DST/$f.part" bs=$bs skip=$start seek=$start count=$per conv=notrunc status=none
}
small() {  # small <file>
  local f=$1
  [ -f "$DST/$f" ] && [ "$(stat -c %s "$DST/$f")" = "$(stat -c %s "$SRC/$f")" ] && return 0
  # --preserve=mode fails here ("Operation not supported": NFSv4 ACL -> xfs), so copy, then set mtime/mode separately
  cp "$SRC/$f" "$DST/$f.part" && touch -r "$SRC/$f" "$DST/$f.part" && { chmod --reference="$SRC/$f" "$DST/$f.part" || true; } \
    && mv "$DST/$f.part" "$DST/$f"
}
export -f seg small; export SRC DST
big=$(awk '$1 > 2000000000 {print $2}' "$DST/.filelist")
for f in $big; do
  [ -f "$DST/$f" ] && [ "$(stat -c %s "$DST/$f")" = "$(stat -c %s "$SRC/$f")" ] && continue
  rm -f "$DST/$f.part"
done
{
  for f in $big; do
    [ -f "$DST/$f" ] && continue
    s=$(stat -c %s "$SRC/$f")
    for i in 0 1 2 3 4 5 6 7; do echo "seg $f $i 8 $s"; done
  done
  awk '$1 <= 2000000000 {print "small " $2}' "$DST/.filelist"
} | xargs -P 16 -L 1 bash -c '"$@"' _
rc=0
for f in $big; do
  if [ -f "$DST/$f.part" ]; then
    truncate -s "$(stat -c %s "$SRC/$f")" "$DST/$f.part"
    if cmp -s "$SRC/$f" "$DST/$f.part"; then
      touch -r "$SRC/$f" "$DST/$f.part"; chmod --reference="$SRC/$f" "$DST/$f.part" 2>/dev/null; mv "$DST/$f.part" "$DST/$f"
      echo "verified $f (cmp)"
    else
      echo "MISMATCH $f"; rc=1
    fi
  fi
done
missing=0
while read -r size f; do
  [ "$(stat -c %s "$DST/$f" 2>/dev/null)" = "$size" ] || { echo "BAD $f"; missing=$((missing + 1)); }
done < "$DST/.filelist"
echo "files: $(wc -l < "$DST/.filelist"), bad or missing: $missing"
[ $missing -eq 0 ] && [ $rc -eq 0 ]
