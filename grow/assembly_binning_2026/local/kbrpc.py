"""Minimal read-only KBase JSON-RPC helper (stdlib). Token read from ~/.kbase/token, never printed."""
import http.client, json, os, sys, time, urllib.request, urllib.error, socket

URLS = {
    "ws": "https://kbase.us/services/ws",
    "hs": "https://kbase.us/services/handle_service",
    "nms": "https://kbase.us/services/narrative_method_store/rpc",
    "catalog": "https://kbase.us/services/catalog",
    "ee2": "https://kbase.us/services/ee2",
}
_TOKEN = None

def token():
    global _TOKEN
    if _TOKEN is None:
        _TOKEN = open(os.path.expanduser("~/.kbase/token")).read().strip()
    return _TOKEN

READ_ONLY = {  # whitelist: nothing that writes
    "Workspace.get_objects2", "Workspace.get_object_info3", "Workspace.list_objects", "Workspace.get_workspace_info",
    "Workspace.get_type_info", "Workspace.get_module_info", "Workspace.list_workspace_info", "Workspace.ver",
    "AbstractHandle.hids_to_handles", "AbstractHandle.fetch_handles_by",
    "NarrativeMethodStore.get_method_spec", "NarrativeMethodStore.list_methods", "NarrativeMethodStore.get_method_full_info",
    "NarrativeMethodStore.list_methods_spec",
    "Catalog.get_module_info", "Catalog.list_client_group_configs", "Catalog.get_exec_aggr_stats",
    "Catalog.get_module_version", "Catalog.list_released_module_versions", "Catalog.get_exec_aggr_table",
    "execution_engine2.check_workspace_jobs", "execution_engine2.check_jobs", "execution_engine2.check_job",
    "execution_engine2.get_job_logs", "execution_engine2.status", "execution_engine2.get_client_groups",
}

def rpc(svc, method, params, auth=True, tries=8, timeout=300):
    """Transport failures (no answer, a proxy 5xx without a JSON-RPC body) are retried with exponential backoff
    (5 s .. 10 min, ~20 min in all); a JSON-RPC error is the service's answer and is raised on its 3rd repeat."""
    if method not in READ_ONLY:
        raise SystemExit(f"refusing non-read-only method {method}")
    url = URLS.get(svc, svc)
    body = json.dumps({"version": "1.1", "id": "1", "method": method, "params": params}).encode()
    hdr = {"Content-Type": "application/json"}
    if auth:
        hdr["Authorization"] = token()
    app_errors, err = 0, None
    for k in range(tries):
        req = urllib.request.Request(url, data=body, headers=hdr)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                res = json.load(r)
        except urllib.error.HTTPError as e:
            try:
                res = json.load(e)
            except Exception:
                res = None
                if e.code < 500:                       # a 4xx without a JSON-RPC body will not change on retry
                    raise RuntimeError(f"{method}: HTTP {e.code} from {url}")
        except (urllib.error.URLError, TimeoutError, socket.timeout, ConnectionError, http.client.HTTPException,
                ValueError):
            res = None
        if res is not None:
            if not res.get("error"):
                return res["result"]
            err = res["error"]
            app_errors += 1
            if app_errors >= 3:
                raise RuntimeError(f"{method}: {str((err.get('message') if isinstance(err, dict) else None) or err)[:600]}")
        if k < tries - 1:
            time.sleep(min(5 * 2 ** k, 600))
    raise RuntimeError(f"{method}: no response from {url} after {tries} tries"
                       + (f" (last JSON-RPC error: {str(err)[:300]})" if err else ""))
