#!/usr/bin/env bash
# Status of the salterns Arm C rollout. On poplar:
#   bash /scratch1/afreiburger/salterns_2026/rollout/scripts/status.sh [root]
# root defaults to the real run; pass /scratch1/afreiburger/salterns_2026/rollout/smoke for the smoke test.
R=${1:-/scratch1/afreiburger/salterns_2026/rollout}
PY=/scratch1/afreiburger/salterns_2026/envs/bt_map/bin/python
S=/scratch1/afreiburger/salterns_2026
echo "== $(date '+%F %T %Z')  $(uptime | sed 's/.*load/load/')   MemAvailable $(awk '/MemAvailable/{printf "%.0f GB", $2/1048576}' /proc/meminfo)"
pid=$(cat $R/driver.pid 2>/dev/null); alive=no; [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null && alive=yes
echo "== driver PID ${pid:-none} alive=$alive   STOP file: $([ -e $R/STOP ] && echo yes || echo no)"
echo "== end marker: $( { cat $R/DONE 2>/dev/null || cat $R/FAILED 2>/dev/null; } || echo 'none yet')"
echo "== GTDB-Tk setup: env $([ -e $S/envs/.ro_setup/gtdbtk280.ok ] && echo ok || echo pending)  R232 $([ -e $S/refdata/gtdbtk_r232/.ok ] && echo "ok ($(cat $S/refdata/gtdbtk_r232/.ok))" || { [ -e $S/refdata/gtdbtk_r232/.FAIL ] && echo FAILED || echo pending; })  last: $(tail -1 $S/rollout/logs/setup_gtdbtk.log 2>/dev/null | cut -c1-120)"
echo "== disk: /scratch1 $(df -h --output=avail /scratch1 | tail -1 | tr -d ' ') free; /tmp $(df -h --output=avail /tmp | tail -1 | tr -d ' ') free"
# samples run on another host (R2A migration): staged by its push task, imported here by import:<s>
if [ -d $R/imports ] && [ -n "$(ls $R/imports 2>/dev/null)" ]; then
  echo "== staged imports (imports/<s>: READY = pushed and md5-verified; import = .done/import:<s>):"
  for d in $R/imports/*/; do s=$(basename $d); [ "$s" = "*" ] && continue
    printf '     %-12s READY=%s import=%s%s\n' "$s" "$([ -e $d/READY ] && echo yes || echo no)" \
      "$([ -e $R/.done/import:$s ] && echo done || { [ -e $R/.failed/import:$s ] && echo FAILED || echo no; })" \
      "$([ -e $d/FAILED ] && echo "  remote FAILED: $(head -c 120 $d/FAILED)")"
  done
fi
[ -f $R/STATE.json ] && $PY -I - $R/STATE.json <<'PYEOF'
import json, sys, collections
st = json.load(open(sys.argv[1]))
sc = st.get("sched") or {}
print(f"== STATE.json {st['utc']}  note={st.get('note')!r}")
if sc:
    print(f"   scheduler {sc.get('utc')}: load1 {sc.get('load1')}  ours {sc.get('ours')} thr  others {sc.get('others')}  "
          f"cap {sc.get('cap')} (literal {sc.get('cap_literal')}, guard {sc.get('cap_guard')})  "
          f"MemAvailable {sc.get('mem_available_gb')} GB, outstanding reserve {sc.get('outstanding_reserve_gb')} GB  "
          f"rollout usage {sc.get('rollout_usage_gb')} GB  ready {sc.get('n_ready')}  asm lanes {sc.get('asm_lanes')}")
    sk = collections.Counter((sc.get("skipped") or {}).values())
    if sk:
        print("   held back (first 40 ready tasks):", dict(sk))
c = collections.Counter(st["tasks"].values())
print("   tasks:", dict(c))
ss = collections.Counter(v["state"] for v in st["samples"].values())
print("   samples:", dict(ss))
for s, v in st["samples"].items():
    i = v.get("info") or {}
    extra = f"last_done={i.get('last_done')} running={','.join(i.get('running') or [])}" if v["state"] == "running" else json.dumps(i)[:150]
    print(f"     {s:22s} {v['state']:9s} {extra}")
xs = [n for n, x in st["tasks"].items() if x == "failed"]
if xs:
    print("   FAILED tasks:", " ".join(xs)[:800])
print(f"   running ({len(st['running'])}):")
for n, r in sorted(st["running"].items(), key=lambda x: -x[1]["elapsed_h"]):
    print(f"     {n:40s} att {r['attempt']} pid {r['pid']:>8} {r['cost']:3d} thr  {r['elapsed_h']:6.2f} h  rss {r['rss_gb']:6.1f} GB"
          f"{'  ' + json.dumps({k: v for k, v in (r.get('params') or {}).items() if k in ('mode', 'mem_gb', 'heap')}) if r.get('params') else ''}")
PYEOF
echo "== last driver log lines:"; tail -8 $R/logs/driver.log 2>/dev/null | sed 's/^/   /'
[ -f $R/results/REPORT.md ] && echo "== report: $R/results/REPORT.md ($(stat -c %y $R/results/REPORT.md | cut -c1-19))"
exit 0
