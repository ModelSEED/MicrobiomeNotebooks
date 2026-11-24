#!/usr/bin/env bash
# Status of the salterns binning test. On poplar:
#   bash /scratch1/afreiburger/salterns_2026/binning_test/scripts/status.sh [root]
# root defaults to the real run; pass /scratch1/afreiburger/salterns_2026/binning_test/smoke for the smoke test.
R=${1:-/scratch1/afreiburger/salterns_2026/binning_test}
PY=/scratch1/afreiburger/salterns_2026/envs/bt_map/bin/python
echo "== $(date '+%F %T %Z')  $(uptime | sed 's/.*load/load/')"
pid=$(cat $R/driver.pid 2>/dev/null); alive=no; [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null && alive=yes
echo "== driver PID ${pid:-none} alive=$alive   DONE: $(cat $R/DONE 2>/dev/null || echo 'not yet')   STOP file: $([ -e $R/STOP ] && echo yes || echo no)"
P=/scratch1/afreiburger/salterns_2026/tadpole_pilot/SF2_C_D1
echo "== SF2_C_D1 Arm C assembly: marker $([ -e $P/.done/C1_spades_default ] && echo PRESENT || echo absent); spades.py 3634657 $(kill -0 3634657 2>/dev/null && echo alive || echo gone); last: $(grep -E 'INFO|ERROR' $P/armC/spades/spades.log 2>/dev/null | tail -1 | cut -c1-110)"
[ -f $R/STATE.json ] && $PY -I - $R/STATE.json <<'EOF'
import json, sys, collections
st = json.load(open(sys.argv[1]))
print(f"== STATE.json {st['utc']}  note={st.get('note')!r}  threads {st['threads_used']}/{st['thread_cap']}  SF2 deadline {st['sf2_deadline_utc']}")
for s, v in st["samples"].items():
    print(f"   {s}: {v['state']}  {json.dumps(v['info'])[:200] if v['info'] else ''}")
c = collections.Counter(st["tasks"].values())
print("   tasks:", dict(c))
for k in ("failed", "blocked"):
    xs = [n for n, v in st["tasks"].items() if v == k]
    if xs:
        print(f"   {k}: {' '.join(xs)[:600]}")
print("   running:", " ".join(st["running"]))
EOF
echo "== last driver log lines:"; tail -6 $R/logs/driver.log 2>/dev/null | sed 's/^/   /'
[ -f $R/DECISION.json ] && { echo "== DECISION.json:"; head -c 600 $R/DECISION.json; echo; }
exit 0
