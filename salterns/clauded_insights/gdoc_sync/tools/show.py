import json, sys
sys.path.insert(0, sys.argv[3])
from word_diff import wdiff
old = json.load(open(sys.argv[1])); new = json.load(open(sys.argv[2]))
pairs = [(int(a), int(b)) for a, b in (x.split(":") for x in sys.argv[4:])]
for i, j in pairs:
    o = old[i] if i >= 0 else {"base": "", "proposed": ""}
    n = new[j] if j >= 0 else {"base": "", "proposed": ""}
    print(f"### old {i} -> new {j}")
    print("  base (suggestions rejected):"); print("\n".join(wdiff(o["base"], n["base"])) or "   (same)")
    print("  proposed (accepted):"); print("\n".join(wdiff(o["proposed"], n["proposed"])) or "   (same)")
    print("  old-base vs new-proposed:"); print("\n".join(wdiff(o["base"], n["proposed"])) or "   (same)")
