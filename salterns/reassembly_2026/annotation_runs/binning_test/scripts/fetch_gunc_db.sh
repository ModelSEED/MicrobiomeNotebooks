#!/usr/bin/env bash
# GUNC progenomes_2.1 DIAMOND db (the gunc download_db default), fetched in 8 byte ranges from the EMBL mirror
# (gunc 1.1.1's second URL; zenodo, its first, ran at 1-3 MB/s from poplar; EMBL gives ~12 MB/s per connection but
# answers 403 to some concurrent requests, so each part is retried until it has its exact size) and checked against
# gunc 1.1.1's md5s (gunc/gunc_database.py: gz bc93a855e0760aad5c4e5f2d0e26da46, dmnd 447c9330056b02f29f30fe81fe4af4eb).
set -euo pipefail
D=/scratch1/afreiburger/salterns_2026/binning_test/refdata/gunc
U=https://black.embl.de/~fullam/gunc/gunc_db_progenomes2.1.dmnd.gz
cd $D; rm -f gunc_db_progenomes2.1.dmnd.gz
N=$(curl -sI $U | awk 'tolower($1)=="content-length:"{print $2+0}')
echo "size $N"; P=8; S=$(( (N + P - 1) / P ))
part() {  # part i: fetch until the file has exactly its range size
  local i=$1 a=$(( $1*S )) b=$(( ($1+1)*S - 1 )) k
  [ $b -ge $N ] && b=$((N-1))
  for k in $(seq 1 30); do
    [ -e part.$i ] && [ "$(stat -c %s part.$i)" = "$((b-a+1))" ] && { echo "part $i ok"; return 0; }
    rm -f part.$i; curl -sSf --retry 3 -r $a-$b -o part.$i $U || true
    [ "$(stat -c %s part.$i 2>/dev/null || echo 0)" = "$((b-a+1))" ] || sleep $((5 + RANDOM % 20))
  done
  echo "part $i FAILED"; return 1
}
pids=(); for i in $(seq 0 $((P-1))); do part $i & pids+=($!); sleep 2; done
for p in "${pids[@]}"; do wait $p; done
cat $(for i in $(seq 0 $((P-1))); do echo part.$i; done) > gunc_db_progenomes2.1.dmnd.gz
[ "$(stat -c %s gunc_db_progenomes2.1.dmnd.gz)" = "$N" ]
g=$(md5sum gunc_db_progenomes2.1.dmnd.gz | cut -d' ' -f1); echo "gz md5 $g (want bc93a855e0760aad5c4e5f2d0e26da46)"
[ "$g" = bc93a855e0760aad5c4e5f2d0e26da46 ]; rm -f part.*
/scratch1/afreiburger/salterns_2026/envs/bt_map/bin/pigz -d -p 4 gunc_db_progenomes2.1.dmnd.gz
m=$(md5sum gunc_db_progenomes2.1.dmnd | cut -d' ' -f1); echo "dmnd md5 $m (want 447c9330056b02f29f30fe81fe4af4eb)"
[ "$m" = 447c9330056b02f29f30fe81fe4af4eb ]
touch .ok; ls -la; echo GUNC_DB_OK
