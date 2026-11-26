#!/usr/bin/env bash
# Restart a rollout driver on poplar without disturbing its tasks (used for the R2A migration, 2026-11-16):
#   1. touch STOP and wait until the driver has written a STATE.json saying so (from then on it launches nothing);
#   2. SIGTERM it only while its main thread sleeps between cycles (wchan hrtimer_nanosleep), so no launch, marker or
#      run-record write can be cut in half; tasks are detached (own sessions) and keep running;
#   3. check no driver runs on that root, start the new one with the README launch command (same --config argument),
#      wait until its log shows the adoptions and a written STATE.json; STOP is left in place: remove it after checking
#      (rm $ROOT/STOP), or pass --rm-stop to remove it here once the new driver has adopted every running task.
# usage: bash restart_driver.sh ROOT CONFIG_ARG [--rm-stop]
#   ROOT=/scratch1/afreiburger/salterns_2026/rollout  CONFIG_ARG=scripts/config.json  (exactly as in the launch command)
set -euo pipefail
R=$1; C=$2; RMSTOP=${3:-}
PY=/scratch1/afreiburger/salterns_2026/envs/bt_map/bin/python
cd "$R"
pid=$(cat driver.pid)
grep -q "scripts/driver.py" /proc/$pid/cmdline || { echo "pid $pid is not a driver"; exit 1; }
nrun=$(ls .run | wc -l)
echo "[$(date -u +%FT%TZ)] driver $pid; $nrun running tasks: $(ls .run | tr '\n' ' ')"
touch STOP
for i in $(seq 1 90); do grep -q '"note": "STOP file present' STATE.json 2>/dev/null && break; sleep 2; done
grep -q '"note": "STOP file present' STATE.json || { echo "driver did not acknowledge STOP within 3 min (STOP left in place)"; exit 1; }
echo "[$(date -u +%FT%TZ)] STOP acknowledged in STATE.json"
killed=0
for i in $(seq 1 600); do
  if [ "$(cat /proc/$pid/wchan 2>/dev/null)" = hrtimer_nanosleep ]; then kill -TERM $pid; killed=1; break; fi
  sleep 0.1
done
[ $killed = 1 ] || { echo "driver never seen sleeping in 60 s; not killed (STOP left in place)"; exit 1; }
for i in $(seq 1 50); do kill -0 $pid 2>/dev/null || break; sleep 0.2; done
kill -0 $pid 2>/dev/null && { echo "driver $pid still alive after SIGTERM"; exit 1; }
echo "[$(date -u +%FT%TZ)] driver $pid stopped while sleeping"
n=$(pgrep -fc "^$PY -I scripts/driver.py run --config $C\$" || true)
[ "$n" = 0 ] || { echo "another driver still runs with --config $C"; pgrep -af "^$PY -I scripts/driver.py run --config $C\$"; exit 1; }
nrun2=$(ls .run | wc -l)
lines=$(wc -l < logs/driver.log)
nohup setsid nice -n 10 $PY -I scripts/driver.py run --config "$C" >> logs/driver.log 2>&1 < /dev/null &
sleep 1
for i in $(seq 1 60); do
  [ -s driver.pid ] && [ "$(cat driver.pid)" != "$pid" ] && grep -q '"driver_pid": '"$(cat driver.pid)" STATE.json 2>/dev/null && break
  sleep 2
done
new=$(cat driver.pid)
echo "[$(date -u +%FT%TZ)] new driver $new; log since restart:"
tail -n +$((lines + 1)) logs/driver.log | cut -c1-200
nad=$(tail -n +$((lines + 1)) logs/driver.log | grep -c "adopted " || true)
echo "adopted $nad of $nrun2 running-task records"
if [ "$RMSTOP" = --rm-stop ]; then
  if [ "$nad" = "$nrun2" ] && ! tail -n +$((lines + 1)) logs/driver.log | grep -q "died without an exit code\|cycle error"; then
    rm -f STOP; echo "[$(date -u +%FT%TZ)] STOP removed"
  else
    echo "NOT removing STOP: check the log above"; exit 1
  fi
fi
