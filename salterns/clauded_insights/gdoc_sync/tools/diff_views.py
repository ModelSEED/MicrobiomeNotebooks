import json, sys, difflib
old = json.load(open(sys.argv[1])); new = json.load(open(sys.argv[2]))
key = lambda x: (x["where"].split(" r")[0] if x["where"] != "body" else "body", x["base"], x["proposed"])
a = [key(x) for x in old]; b = [key(x) for x in new]
sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
for tag, i1, i2, j1, j2 in sm.get_opcodes():
    if tag == "equal": continue
    print(f"== {tag} old[{i1}:{i2}] new[{j1}:{j2}]")
    for i in range(i1, i2): print(f"  OLD {i} {old[i]['where']}\n    base: {old[i]['base']!r}\n    prop: {old[i]['proposed']!r}")
    for j in range(j1, j2): print(f"  NEW {j} {new[j]['where']}\n    base: {new[j]['base']!r}\n    prop: {new[j]['proposed']!r}")
