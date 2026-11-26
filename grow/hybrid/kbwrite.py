"""KBase WRITE helper for the GROW hybrid local/KBase mode (user decision 2026-11-23). Token from ~/.kbase/token, never
printed or written anywhere. Separate from the read-only local/kbrpc.py on purpose: everything that can write goes
through the explicit allowlist below, nothing else.

  blob_upload(path, file_name)   -> Blobstore node (raw POST) + persisted handle, in DataFileUtil.file_to_shock's shape
  save_object(ws_id, type, name, data, provenance, meta) -> object_info (Workspace.save_objects, one object)
  blob_delete(node_id)           -> only for nodes this helper made and never referenced (failed saves / self-tests)
  download_blob(node_id, path)   -> streams a Blobstore node to a file (read; md5 checked by the caller)
"""
import hashlib
import http.client
import json
import os
import socket
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

import requests
from requests_toolbelt.multipart.encoder import MultipartEncoder

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "assembly_binning_2026", "local")))
from kbrpc import rpc, token  # noqa: E402  (read-only RPC, shared token loader)

BLOB = "https://kbase.us/services/shock-api"
HANDLE_URL = "https://services.kbase.us/services/shock-api"   # the url KBase's own handles carry
WS_URL = "https://kbase.us/services/ws"
HS_URL = "https://kbase.us/services/handle_service"
WRITE_METHODS = {"Workspace.save_objects", "AbstractHandle.persist_handle"}
GROW_WS = {276048, 277632}                                       # the only workspaces this helper may write to


def _post_rpc(url, method, params, tries=6, timeout=1800):
    if method not in WRITE_METHODS:
        raise SystemExit(f"refusing method {method}")
    body = json.dumps({"version": "1.1", "id": "1", "method": method, "params": params}).encode()
    last = None
    for i in range(tries):
        req = urllib.request.Request(url, data=body, headers={"Content-Type": "application/json",
                                                              "Authorization": token()})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read())["result"]
        except urllib.error.HTTPError as e:
            txt = e.read().decode(errors="replace")
            try:
                err = json.loads(txt).get("error")
            except Exception:
                err = None
            if err:                                              # the service answered: do not blindly retry a write
                raise RuntimeError(f"{method}: {str(err.get('message') if isinstance(err, dict) else err)[:800]}")
            last = f"HTTP {e.code}: {txt[:300]}"
        except (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError,
                http.client.RemoteDisconnected) as e:
            last = repr(e)
        time.sleep(min(600, 10 * 2 ** i))
    raise RuntimeError(f"{method}: transport failure after {tries} tries: {last}")


def md5_file(path, bs=16 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(bs), b""):
            h.update(b)
    return h.hexdigest()


def blob_upload(path, file_name=None, tries=5):
    """Raw-body POST to the Blobstore, then AbstractHandle.persist_handle. Returns DataFileUtil.file_to_shock's dict:
    {shock_id, node_file_name, size, handle: {hid, file_name, id, url, type, remote_md5}}. md5 is verified."""
    file_name = file_name or os.path.basename(path)
    size = os.path.getsize(path)
    local_md5 = md5_file(path)
    last = None
    for i in range(tries):                               # DataFileUtil.file_to_shock's own request shape
        try:
            with open(path, "rb") as fh:
                mpe = MultipartEncoder(fields={"upload": (file_name, fh, None, {"Content-Length": str(size)})})
                r = requests.post(BLOB + "/node", data=mpe, stream=True, allow_redirects=True, timeout=(60, 7200),
                                  headers={"Authorization": "OAuth " + token(), "content-type": mpe.content_type})
            if r.status_code != 200:
                last = f"HTTP {r.status_code}: {r.text[:300]}"
                time.sleep(min(600, 15 * 2 ** i))
                continue
            node = r.json()["data"]
            break
        except (OSError, requests.RequestException) as e:
            last = repr(e)
            time.sleep(min(600, 15 * 2 ** i))
    else:
        raise RuntimeError(f"blob upload of {path} failed: {last}")
    remote_md5 = (node.get("file") or {}).get("checksum", {}).get("md5")
    if remote_md5 != local_md5 or int((node.get("file") or {}).get("size", -1)) != size:
        blob_delete(node["id"])
        raise RuntimeError(f"blob upload of {path}: md5/size mismatch (local {local_md5}/{size}, node {remote_md5})")
    hid = _post_rpc(HS_URL, "AbstractHandle.persist_handle", [{
        "id": node["id"], "file_name": file_name, "type": "shock", "url": HANDLE_URL, "remote_md5": remote_md5}])[0]
    return {"shock_id": node["id"], "node_file_name": file_name, "size": size,
            "handle": {"hid": hid, "file_name": file_name, "id": node["id"], "url": HANDLE_URL, "type": "shock",
                       "remote_md5": remote_md5}}


def blob_delete(node_id):
    req = urllib.request.Request(f"{BLOB}/node/{node_id}", method="DELETE",
                                 headers={"Authorization": "OAuth " + token()})
    with urllib.request.urlopen(req, timeout=300) as r:
        return r.status


def download_blob(node_id, path, expect_md5=None, tries=5):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(f"{BLOB}/node/{node_id}?download_raw",
                                         headers={"Authorization": "OAuth " + token()})
            h = hashlib.md5()
            with urllib.request.urlopen(req, timeout=600) as r, open(path + ".part", "wb") as out:
                for b in iter(lambda: r.read(8 << 20), b""):
                    out.write(b)
                    h.update(b)
            if expect_md5 and h.hexdigest() != expect_md5:
                last = f"md5 {h.hexdigest()} != {expect_md5}"
                continue
            os.replace(path + ".part", path)
            return h.hexdigest()
        except (OSError, urllib.error.URLError, http.client.HTTPException) as e:
            last = repr(e)
            time.sleep(min(600, 15 * 2 ** i))
    raise RuntimeError(f"download of node {node_id} failed: {last}")


def save_object(ws_id, type_, name, data, provenance, meta=None):
    """One Workspace.save_objects call; returns object_info. Only ws 276048 / 277632."""
    if int(ws_id) not in GROW_WS:
        raise SystemExit(f"refusing to write to ws {ws_id}")
    obj = {"type": type_, "data": data, "name": name, "provenance": provenance}
    if meta:
        obj["meta"] = {str(k): str(v) for k, v in meta.items()}
    return _post_rpc(WS_URL, "Workspace.save_objects", [{"id": int(ws_id), "objects": [obj]}])[0][0]


def save_objects(ws_id, objs):
    """Several objects in ONE Workspace.save_objects call (all or nothing); objs: [{type, name, data, provenance,
    meta?}]. Returns the object_infos in order."""
    if int(ws_id) not in GROW_WS:
        raise SystemExit(f"refusing to write to ws {ws_id}")
    for o in objs:
        if o.get("meta"):
            o["meta"] = {str(k): str(v) for k, v in o["meta"].items()}
    return _post_rpc(WS_URL, "Workspace.save_objects", [{"id": int(ws_id), "objects": objs}])[0]


def upa(info):
    return f"{info[6]}/{info[0]}/{info[4]}"


def get_info(ref):
    return rpc("ws", "Workspace.get_object_info3", [{"objects": [{"ref": ref}], "includeMetadata": 1}])[0]["infos"][0]
