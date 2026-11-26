#!/usr/bin/env bash
# Dry run of a config change on a running rollout root (R2 migration, 2026-11-17; same check as the R2A migration):
# copies the root's markers (.done .failed .run imports/) twice under OUT/{old,new}, points a copy of each config at its
# copy, runs `driver.py plan` (state from markers only; no task is launched), then compares the two plans:
#   tasks only in old / only in new, state changes of common tasks, sample state changes, and markers (.done/.failed)
#   that the new task set no longer uses but the old one did.
# usage: bash dryrun_migration.sh ROOT PY OLD_CONFIG NEW_CONFIG OUT
set -euo pipefail
ROOT=$1; PY=$2; OLD=$3; NEW=$4; OUT=$5
rm -rf "$OUT"; mkdir -p "$OUT"
for v in old new; do
  mkdir -p "$OUT/$v/scripts"
  for m in .done .failed .run imports; do [ -d "$ROOT/$m" ] && cp -a "$ROOT/$m" "$OUT/$v/" || mkdir -p "$OUT/$v/$m"; done
  cp -p "$ROOT/scripts/driver.py" "$ROOT/scripts/rolib.py" "$ROOT/scripts/btlib.py" "$OUT/$v/scripts/"
done
"$PY" -I - "$OLD" "$NEW" "$OUT" <<'EOF'
import json, sys
old, new, out = sys.argv[1:4]
for v, p in (("old", old), ("new", new)):
    c = json.load(open(p)); c["root"] = f"{out}/{v}"
    json.dump(c, open(f"{out}/{v}/scripts/config.json", "w"), indent=1)
EOF
for v in old new; do
  ( cd "$OUT/$v" && "$PY" -I scripts/driver.py plan --config scripts/config.json > plan.txt 2> plan.err ) || { echo "plan $v failed"; cat "$OUT/$v/plan.err"; exit 1; }
done
"$PY" -I - "$OUT" <<'EOF'
import json, os, sys
out = sys.argv[1]
S = {v: json.load(open(f"{out}/{v}/STATE.json")) for v in ("old", "new")}
T = {v: S[v]["tasks"] for v in S}
only_old = sorted(set(T["old"]) - set(T["new"])); only_new = sorted(set(T["new"]) - set(T["old"]))
print(f"tasks: old {len(T['old'])}  new {len(T['new'])}  only-old {len(only_old)}  only-new {len(only_new)}")
from collections import Counter
print("  only-old by sample:", dict(Counter(n.split(':')[-1] if n.split(':')[0] != 'map' else n.split(':')[1] for n in only_old)))
print("  only-old by state:", dict(Counter(T['old'][n] for n in only_old)))
print("  only-new:", only_new)
ch = {n: (T["old"][n], T["new"][n]) for n in set(T["old"]) & set(T["new"]) if T["old"][n] != T["new"][n]}
print(f"common tasks whose state changed: {len(ch)}", sorted(ch.items())[:20])
sc = {s: (S['old']['samples'].get(s, {}).get('state'), S['new']['samples'].get(s, {}).get('state'))
      for s in set(S['old']['samples']) | set(S['new']['samples'])
      if S['old']['samples'].get(s, {}).get('state') != S['new']['samples'].get(s, {}).get('state')}
print("sample state changes:", sc)
marks = set(os.listdir(f"{out}/new/.done")) | set(os.listdir(f"{out}/new/.failed"))
orph = sorted(m for m in marks if m in T["old"] and m not in T["new"])
print(f"markers used by the old task set but not the new: {len(orph)}", orph)
running = sorted(r[:-5].replace("__", ":") for r in os.listdir(f"{out}/new/.run") if r.endswith(".json"))
print("running (.run) all still tasks in new:", all(r in T["new"] for r in running), running)
EOF
