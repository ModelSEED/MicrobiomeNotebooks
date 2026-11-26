"""Field-level diff of ModelSEED template JSON files between git revisions."""
import json, subprocess, sys
from functools import lru_cache

V7 = ["GramNegModelTemplateV7.json", "GramPosModelTemplateV7.json", "ArchaeaTemplateV6.json", "Core-V6.json"]

@lru_cache(maxsize=64)
def load(repo, rev, path):
    raw = subprocess.run(["git", "-C", repo, "show", f"{rev}:{path}"], capture_output=True, check=True).stdout
    return json.loads(raw)

def index(t):
    out = {}
    for k, v in t.items():
        if isinstance(v, list) and v and all(isinstance(x, dict) and "id" in x for x in v):
            out[k] = {x["id"]: x for x in v}
        else:
            out[k] = v
    return out

def diff(a, b):
    A, B = index(a), index(b)
    res = {}
    for k in sorted(set(A) | set(B)):
        va, vb = A.get(k), B.get(k)
        if isinstance(va, dict) and isinstance(vb, dict) and k not in ("info",):
            added = sorted(set(vb) - set(va)); removed = sorted(set(va) - set(vb))
            changed = {}
            for i in sorted(set(va) & set(vb)):
                if va[i] != vb[i]:
                    fa, fb = va[i], vb[i]
                    changed[i] = {f: (fa.get(f, "<absent>"), fb.get(f, "<absent>"))
                                  for f in sorted(set(fa) | set(fb)) if fa.get(f, "<absent>") != fb.get(f, "<absent>")}
            if added or removed or changed:
                res[k] = {"added": {i: vb[i] for i in added}, "removed": {i: va[i] for i in removed}, "changed": changed}
        elif va != vb:
            res[k] = {"value": (va, vb)}
    return res

def summarize(d):
    lines = []
    for sec, v in d.items():
        if "value" in v:
            lines.append(f"  {sec}: scalar changed"); continue
        if v["added"]: lines.append(f"  {sec}: +{len(v['added'])} {list(v['added'])[:8]}")
        if v["removed"]: lines.append(f"  {sec}: -{len(v['removed'])} {list(v['removed'])[:8]}")
        if v["changed"]:
            fields = {}
            for i, ch in v["changed"].items():
                for f in ch: fields.setdefault(f, []).append(i)
            for f, ids in sorted(fields.items()):
                lines.append(f"  {sec}.{f}: {len(ids)} {ids[:12]}{' ...' if len(ids)>12 else ''}")
    return "\n".join(lines) if lines else "  (no change)"

if __name__ == "__main__":
    repo, parent, child = sys.argv[1:4]
    files = sys.argv[4:] or V7
    for f in files:
        p = f if "/" in f else f"templates/v7.0/{f}"
        try:
            a = load(repo, parent, p); b = load(repo, child, p)
        except subprocess.CalledProcessError:
            print(f"{p}: missing in one revision"); continue
        print(f"{p}:"); print(summarize(diff(a, b)))
