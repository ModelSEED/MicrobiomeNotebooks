#!/usr/bin/env bash
# Memory watchdog for the GTDB-Tk run (systemd user scope gtdbtk-wetland.scope) on this shared box.
# Every 30 s: append a status line to mem_guard.log. Alert lines (also echoed to stdout, which the
# Claude Monitor turns into notifications) go to mem_guard.events:
#   WARN/CRIT  MemAvailable below 40 / 20 GB, with the top memory users outside our scope
#   BIG        a process outside our scope newly above 15 GB RSS (another session's job)
#   CAP        our scope above 85% of its MemoryMax
#   FREEZE     MemAvailable < 12 GB: freeze our scope (stops it allocating; never touches other jobs)
#   THAW       MemAvailable back > 30 GB: unfreeze
#   OOM        oom_kill events counted in our scope's memory.events
#   END        our scope is gone (run finished or died) -> guard exits
#   END        our scope never appeared within 15 min (launch failed) -> guard exits
# Usage: UNIT=<unit>.scope GUARD_DIR=<dir> setsid nohup bash mem_guard.sh > /dev/null 2>&1 &
set -u
cd "${GUARD_DIR:-$(dirname "$0")}"        # logs go next to the run (GUARD_DIR), else next to this script
UNIT=${UNIT:-gtdbtk-wetland.scope}
CG=/sys/fs/cgroup/user.slice/user-$(id -u).slice/user@$(id -u).service/app.slice/$UNIT
LOG=mem_guard.log; EV=mem_guard.events
declare -A seen_big
frozen=0; last_oom=0; started=0; t0=$(date +%s)

ev() { local l="[$(date '+%F %T')] $*"; echo "$l" >> "$EV"; echo "$l"; }
# A guard that is stopped while it holds the scope frozen thaws it if memory allows; else it says so (the watcher
# restarts a guard, which adopts the freeze below). SIGKILL cannot be trapped: that case is the adoption.
thaw_exit() {
  local a; a=$(awk '/MemAvailable/ {printf "%d", $2/1048576}' /proc/meminfo)
  if [ "$frozen" -eq 1 ] && [ -e "$CG/cgroup.freeze" ]; then
    if [ "$a" -gt 30 ]; then echo 0 > "$CG/cgroup.freeze"; ev "THAW guard stopping (avail=${a}G): our run resumed"
    else ev "GUARD STOPPING avail=${a}G: scope left frozen; the watcher restarts a guard"; fi
  fi
  exit 0
}
trap thaw_exit TERM INT HUP
top_others() {
  ps -eo pid=,rss=,user=,cgroup:200=,args= --sort=-rss | grep -v "$UNIT" | head -4 |
    awk '{printf "%s %s %.1fG %s | ", $1, $3, $2/1048576, $5" "$6" "$7}'
}

while true; do
  avail=$(awk '/MemAvailable/ {printf "%d", $2/1048576}' /proc/meminfo)
  swap=$(awk '/SwapTotal/ {t=$2} /SwapFree/ {f=$2} END {printf "%.1f", (t-f)/1048576}' /proc/meminfo)
  if [ -d "$CG" ]; then
    [ "$started" -eq 0 ] && [ "$(cat "$CG/cgroup.freeze" 2>/dev/null)" = 1 ] && frozen=1   # adopt a freeze a dead guard left
    started=1
    cur=$(( $(cat "$CG/memory.current") / 1073741824 ))
    max=$(cat "$CG/memory.max"); [ "$max" = max ] && maxg=0 || maxg=$(( max / 1073741824 ))
    anon=$(awk '/^anon / {printf "%d", $2/1073741824}' "$CG/memory.stat")
    oom=$(awk '/^oom_kill / {print $2}' "$CG/memory.events")
    stage=$(pgrep -a -P "$(pgrep -f 'gtdbtk classify_wf --batchfil[e]' | head -1)" 2>/dev/null |
            awk '{print $2}' | xargs -r -n1 basename | sort | uniq -c | tr '\n' ' ')
    pp=$(pgrep -x pplacer | head -1); pprss=$([ -n "$pp" ] && awk '/VmRSS/ {printf "%.1f", $2/1048576}' /proc/$pp/status 2>/dev/null)
    echo "$(date '+%F %T') avail=${avail}G swap=${swap}G scope=${cur}G anon=${anon}G cap=${maxg}G pplacer=${pprss:-0}G frozen=$frozen children: $stage" >> "$LOG"
    [ "$oom" -gt "$last_oom" ] && ev "OOM our scope recorded oom_kill=$oom (cap ${maxg}G)" && last_oom=$oom
    if [ "$maxg" -gt 0 ] && [ "$anon" -gt $(( maxg * 85 / 100 )) ] && [ $(( $(date +%s) - ${last_cap:-0} )) -ge 600 ]; then
      ev "CAP our anon memory ${anon}G > 85% of ${maxg}G"; last_cap=$(date +%s)
    fi
    if [ "$avail" -lt 12 ] && [ "$frozen" -eq 0 ]; then
      echo 1 > "$CG/cgroup.freeze" && frozen=1 && ev "FREEZE avail=${avail}G: our run paused. Others: $(top_others)"
    elif [ "$avail" -gt 30 ] && [ "$frozen" -eq 1 ]; then
      echo 0 > "$CG/cgroup.freeze" && frozen=0 && ev "THAW avail=${avail}G: our run resumed"
    fi
  else
    echo "$(date '+%F %T') avail=${avail}G swap=${swap}G scope=absent" >> "$LOG"
    if [ "$started" -eq 1 ]; then ev "END scope $UNIT is gone (run finished or died)"; exit 0; fi
    if [ "$started" -eq 0 ] && [ $(( $(date +%s) - t0 )) -gt 900 ]; then ev "END scope $UNIT never appeared (15 min)"; exit 1; fi
  fi
  if [ "$avail" -lt 20 ]; then ev "CRIT avail=${avail}G swap=${swap}G | $(top_others)"
  elif [ "$avail" -lt 40 ]; then ev "WARN avail=${avail}G swap=${swap}G | $(top_others)"; fi
  while read -r pid rss args; do
    [ -n "${seen_big[$pid]:-}" ] && continue
    seen_big[$pid]=1
    ev "BIG pid $pid $(awk -v r="$rss" 'BEGIN {printf "%.1fG", r/1048576}') cwd=$(readlink /proc/$pid/cwd 2>/dev/null) :: ${args:0:120}"
  done < <(ps -eo pid=,rss=,cgroup:200=,args= | grep -v "$UNIT" | awk '$2 > 15728640 {pid=$1; rss=$2; $1=$2=$3=""; print pid, rss, $0}')
  sleep 30 & wait $!                     # so a TERM is handled at once
done
