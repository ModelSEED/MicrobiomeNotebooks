#!/bin/bash
# KOfam release download for salterns_2026 (bakta_test track).
# genome.jp serves ~0.4 MB/s per connection, so profiles.tar.gz is fetched as 16 byte ranges in parallel and concatenated.
# Release = ko_list.gz Last-Modified 2026-10-26 05:26 GMT; profiles.tar.gz Last-Modified 2026-10-27 01:41 GMT (1,567,491,980 bytes).
set -euo pipefail
cd /scratch1/afreiburger/salterns_2026/kofam_db
URL=https://www.genome.jp/ftp/db/kofam
date
[ -s ko_list.gz ] || curl -sS -R -o ko_list.gz $URL/ko_list.gz
SIZE=$(curl -sI $URL/profiles.tar.gz | awk "tolower(\$1)==\"content-length:\"{print \$2+0}")
echo "size $SIZE"
N=16; CH=$(( (SIZE + N - 1) / N )); mkdir -p parts
for i in $(seq 0 $((N-1))); do
  a=$((i*CH)); b=$(( a+CH-1 )); [ $b -ge $SIZE ] && b=$((SIZE-1))
  ( for try in 1 2 3; do curl -sS --retry 3 -r $a-$b -o parts/p$(printf %02d $i) $URL/profiles.tar.gz && [ $(stat -c %s parts/p$(printf %02d $i)) -eq $((b-a+1)) ] && break; done ) &
done
wait
cat parts/p* > profiles.tar.gz
[ $(stat -c %s profiles.tar.gz) -eq $SIZE ] || { echo "SIZE MISMATCH"; exit 1; }
rm -rf parts
gzip -t profiles.tar.gz && echo "gzip OK"
md5sum ko_list.gz profiles.tar.gz > md5sums.txt
gunzip -kf ko_list.gz
tar -xzf profiles.tar.gz
echo "profiles: $(ls profiles | wc -l)"
date
touch .download_done
