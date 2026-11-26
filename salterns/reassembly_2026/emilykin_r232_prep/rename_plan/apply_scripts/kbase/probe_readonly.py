"""Read-only probe of ws 262319: AssemblySet 585, one Assembly, one Genome info, type info. Token never printed."""
import json, os, sys, urllib.request
T = open(os.path.expanduser("~/.kbase/token")).read().strip()
def rpc(m, p):
    assert m.split(".")[1] in ("get_objects2", "get_object_info3", "get_type_info", "get_workspace_info", "ver")
    r = urllib.request.Request("https://kbase.us/services/ws", data=json.dumps({"version": "1.1", "id": "1", "method": m, "params": p}).encode(),
                               headers={"Authorization": T, "Content-Type": "application/json"})
    return json.load(urllib.request.urlopen(r, timeout=300))["result"]
s = rpc("Workspace.get_objects2", [{"objects": [{"ref": "262319/585"}]}])[0]["data"][0]
print("SET info", s["info"][:10]); print("SET meta", s["info"][10])
d = s["data"]; print("SET keys", list(d)); print("description:", d.get("description")); print("items[0:2]", d["items"][:2], "n", len(d["items"]))
print("provenance", json.dumps(s.get("provenance"))[:1500])
a = rpc("Workspace.get_objects2", [{"objects": [{"ref": "262319/10"}], "no_data": 0}])[0]["data"][0]
print("ASM info", a["info"]); print("ASM keys", list(a["data"]))
print("ASM fields", {k: (v if not isinstance(v, (dict, list)) else f"<{type(v).__name__} {len(v)}>") for k, v in a["data"].items()})
print("ASM prov", json.dumps(a.get("provenance"))[:800])
print("ASM size", len(json.dumps(a["data"])))
ti = rpc("Workspace.get_type_info", ["KBaseGenomeAnnotations.Assembly-5.1"])[0]
print("type has std_lineages:", "std_lineages" in ti["spec_def"])
