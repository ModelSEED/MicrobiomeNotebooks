#!/usr/bin/env python3
"""Move narrative 272312 cell A's POOL_CUTOFF (per Andrew 2026-11-17: put the MAGs of the 11 re-assembled libraries through
the same pipeline). Only the POOL_CUTOFF statement and its comment change; cells B and C and every other cell stay byte-
identical. Refuses unless the narrative is at --expect-version and cell A equals cells_v<expect>/cell_A.py.
Usage: KB_AUTH_TOKEN=... python -I -B save_cells_cutoff.py --expect-version 36 --cutoff 2026-11-20T00:00:00+0000 [--apply]
"""
import argparse, json, os, re, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/home/freiburger/gtdbtk_kbase")
import kbase_mags as k

p = argparse.ArgumentParser()
p.add_argument("--expect-version", type=int, required=True)
p.add_argument("--cutoff", required=True)
p.add_argument("--apply", action="store_true")
a = p.parse_args()
assert re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\+0000", a.cutoff), a.cutoff
obj = k.rpc(k.WS, "Workspace.get_objects2", [{"objects": [{"ref": "272312/1"}]}])[0]["data"][0]
info, data = obj["info"], obj["data"]
if info[4] != a.expect_version:
    sys.exit(f"narrative is at version {info[4]}, expected {a.expect_version}")
src = lambda c: "".join(c["source"]) if isinstance(c.get("source"), list) else (c.get("source") or "")
old_a = open(os.path.join(HERE, f"cells_v{info[4]}", "cell_A.py")).read()
m = re.search(r'^POOL_CUTOFF = "[^"]+"[^\n]*\n(?:[ ]{20,}#[^\n]*\n)*', old_a, re.M)
assert m, "POOL_CUTOFF block not found"
new_block = (f'POOL_CUTOFF = "{a.cutoff}"  # pool = bins from DAS Tool runs submitted before this: every DAS\n'
             '                                          # Tool run in this narrative, incl. the 11 libraries re-assembled\n'
             '                                          # 2026-11-15/16 and binned 2026-11-17. Earlier values: 2026-11-05\n'
             '                                          # (12,859 MAGs), 2026-10-30 (12,611), 2026-10-08 (2,449)\n')
new_a = old_a[:m.start()] + new_block + old_a[m.end():]
compile(new_a, "cell_A", "exec")
hit = 0
for c in data["cells"]:
    if src(c).startswith("# MAG annotation A:"):
        if src(c) != old_a:
            sys.exit("cell A differs from the saved copy; someone edited it")
        c["source"] = new_a.splitlines(True) if isinstance(c.get("source"), list) else new_a
        hit += 1
assert hit == 1
print(f"cell A: {m.group(0).splitlines()[0][:60]} -> POOL_CUTOFF = \"{a.cutoff}\"")
if not a.apply:
    sys.exit("dry run; re-run with --apply")
prov = [{"time": time.strftime("%Y-%m-%dT%H:%M:%S+0000", time.gmtime()), "script": "narrative_run/save_cells_cutoff.py",
         "description": f"cell A POOL_CUTOFF -> {a.cutoff}", "input_ws_objects": [k.upa(info)]}]
saved = k.rpc(k.WS, "Workspace.save_objects", [{"id": k.WSID, "objects": [
    {"type": info[2], "objid": info[0], "data": data, "meta": info[10], "provenance": prov}]}], tries=1)[0][0]
back = k.rpc(k.WS, "Workspace.get_objects2", [{"objects": [{"ref": k.upa(saved)}]}])[0]["data"][0]["data"]
d = os.path.join(HERE, f"cells_v{saved[4]}")
os.makedirs(d, exist_ok=True)
same = 0
for c0, c1 in zip(obj["data"]["cells"], back["cells"]):
    s = src(c1)
    lab = next((l for l in "ABC" if s.startswith(f"# MAG annotation {l}:")), None)
    if lab:
        open(os.path.join(d, f"cell_{lab}.py"), "w").write(s)
    if lab != "A":
        same += json.dumps(c0, sort_keys=True) == json.dumps(c1, sort_keys=True)
print(f"saved {k.upa(saved)}; cell A as intended: {open(os.path.join(d, 'cell_A.py')).read() == new_a}; "
      f"other cells unchanged: {same}/{len(back['cells']) - 1}")
