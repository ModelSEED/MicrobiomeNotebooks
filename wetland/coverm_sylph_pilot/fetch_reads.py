#!/usr/bin/env python3
"""Fetch KBase paired-end read libraries (ws 261700) as their stored blob (interleaved FASTQ, usually gzipped).
Token from stdin (never written to disk). Usage: python3 -I fetch_reads.py OUTDIR NAME [NAME ...]  < token"""
import json, os, sys, urllib.request, hashlib
tok = sys.stdin.read().strip()
if os.fork():                      # detach after reading the token from stdin (it is never written to disk)
    sys.exit(0)
os.setsid()
log = open(os.path.join(sys.argv[1], "..", "fetch.log"), "a", buffering=1)
sys.stdout = sys.stderr = log
WS = "https://kbase.us/services/ws"
def rpc(method, params):
    req = urllib.request.Request(WS, data=json.dumps({"method": method, "params": params, "version": "1.1", "id": "1"}).encode(),
                                 headers={"Authorization": tok})
    return json.load(urllib.request.urlopen(req, timeout=300))["result"]
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
for name in sys.argv[2:]:
    d = rpc("Workspace.get_objects2", [{"objects": [{"wsid": 261700, "name": name}]}])[0]["data"][0]
    lib = d["data"]["lib1"]["file"]; info = d["info"]
    print(name, info[2], "interleaved", d["data"].get("interleaved"), "reads", d["data"].get("read_count"),
          "bases", d["data"].get("total_bases"), "file", lib.get("file_name"), flush=True)
    dest = os.path.join(out, name + ("." + lib["file_name"].split(".", 1)[1] if "." in lib.get("file_name", "") else ".fastq"))
    if os.path.exists(dest + ".done"):
        continue
    req = urllib.request.Request(f"https://kbase.us/services/shock-api/node/{lib['id']}?download_raw",
                                 headers={"Authorization": "OAuth " + tok})
    h = hashlib.md5(); n = 0
    with urllib.request.urlopen(req, timeout=3600) as r, open(dest, "wb") as fh:
        while True:
            b = r.read(1 << 22)
            if not b: break
            fh.write(b); h.update(b); n += len(b)
    print(name, dest, n, "bytes md5", h.hexdigest(), flush=True)
    open(dest + ".done", "w").write(h.hexdigest() + "\n")
