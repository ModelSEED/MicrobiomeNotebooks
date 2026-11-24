#!/usr/bin/env bash
# Status of the tadpole pilot. Run on poplar:  bash /scratch1/afreiburger/salterns_2026/tadpole_pilot/scripts/status.sh
ROOT=/scratch1/afreiburger/salterns_2026/tadpole_pilot
echo "== $(date '+%F %T')  $(uptime | sed 's/.*load/load/')"
echo "== download: $(tail -1 $ROOT/reads/download.log)"; ls $ROOT/reads | grep -E '\.ok$|\.FAIL$|download.done' | tr '\n' ' '; echo
for s in R1_A_D1 SF2_C_D1; do
  D=$ROOT/$s; [ -d $D ] || continue
  pid=$(cat $D/driver.pid 2>/dev/null); alive=no; kill -0 "$pid" 2>/dev/null && alive=yes
  echo "== $s  driver PID $pid alive=$alive  STATUS: $(cat $D/STATUS 2>/dev/null)"
  echo "   done: $(ls $D/.done 2>/dev/null | tr '\n' ' ')"
  tail -3 $D/driver.log | sed 's/^/   /'
  for arm in armA armB armC; do
    L=$D/$arm/spades/spades.log; [ -f $L ] || continue
    fin=$(grep -c "SPAdes pipeline finished" $L)
    echo "   $arm SPAdes: finished=$fin  last: $(grep -E 'INFO|ERROR|==' $L | tail -1 | cut -c1-150)"
  done
done
echo "== my SPAdes/Java processes:"; ps -o pid,etime,pcpu,rss,args -u afreiburger | grep -E "tadpole_pilot" | grep -E "spades-core|java" | grep -v grep \
  | awk '{printf "   %s %s cpu=%s%% rss=%.0fG %s\n",$1,$2,$3,$4/1048576,substr($0,index($0,$5),90)}'
