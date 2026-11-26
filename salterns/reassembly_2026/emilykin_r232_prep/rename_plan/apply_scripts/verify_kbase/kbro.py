"""Read-only KBase access for the independent kbase verifier.  Standard library only.

Only the methods in READ_ONLY can be called; blob nodes are only ever fetched with GET.
Token: ~/.kbase/token (never printed or written anywhere).
"""
import json, os, socket, time, urllib.error, urllib.request

WS = "https://kbase.us/services/ws"
HS = "https://kbase.us/services/handle_service"
READ_ONLY = {
    (WS, "Workspace.get_workspace_info"), (WS, "Workspace.list_objects"), (WS, "Workspace.get_object_info3"),
    (WS, "Workspace.get_objects2"), (WS, "Workspace.list_referencing_objects"), (WS, "Workspace.get_object_history"),
    (HS, "AbstractHandle.hids_to_handles"),
}
_T = None


def token():
    global _T
    if _T is None:
        _T = open(os.path.expanduser("~/.kbase/token")).read().strip()
    return _T


def rpc(url, method, params, tries=4):
    if (url, method) not in READ_ONLY:
        raise SystemExit(f"refusing non-read-only call {url} {method}")
    body = json.dumps({"version": "1.1", "id": "1", "method": method, "params": params}).encode()
    last = None
    for k in range(tries):
        req = urllib.request.Request(url, data=body, headers={"Authorization": token(), "Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=600) as r:
                res = json.load(r)
        except urllib.error.HTTPError as e:
            try:
                res = json.load(e)
            except Exception:
                res = None
        except (urllib.error.URLError, TimeoutError, socket.timeout, ConnectionError) as e:
            res, last = None, e
        if res is not None:
            if res.get("error"):
                raise RuntimeError(f"{method}: {str(res['error'].get('message') or res['error'])[:400]}")
            return res["result"]
        time.sleep(5 * (k + 1))
    raise RuntimeError(f"{method}: no response ({last})")


def ws(method, params):
    return rpc(WS, "Workspace." + method, params)


def blob_get(url):
    """GET a blobstore node record (metadata only, not the file). Returns (http_status, json or None)."""
    if "/node/" not in url or "?" in url:
        raise SystemExit(f"refusing blob url {url}")
    # report links record the canonical internal host; fetch the same node through the public host
    url = url.replace("https://services.kbase.us/services/shock-api", "https://kbase.us/services/shock-api")
    req = urllib.request.Request(url, method="GET", headers={"Authorization": "OAuth " + token()})
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return r.status, json.load(r)
    except urllib.error.HTTPError as e:
        try:
            return e.code, json.load(e)
        except Exception:
            return e.code, None
    except (urllib.error.URLError, TimeoutError, socket.timeout, ConnectionError) as e:
        return f"neterr {e}", None


def upa(i):
    return f"{i[6]}/{i[0]}/{i[4]}"
