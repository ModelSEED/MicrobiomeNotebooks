#!/usr/bin/env python
"""Upload the 57 modelling-species MAGs of the 2026 salterns re-assembly to KBase narrative 275645, RAST-annotate them
there, and download the RAST genomes for model building.

Inputs
  ../../comparison_rollout/modelling_species.tsv   57 species representatives (rep_genome = '<set>|<mag_id>')
  ../../comparison_rollout/mags_three_way.tsv      GTDB R232 lineage (domain, classification)
  ../../work/{rollout,coasm}_mags_kept/<mag_id>.fa  FASTAs, md5-checked against ../../{rollout,coasm}_final/mags_kept.md5
Writes (KBase: ws 275645 only, only these names; never overwrites, never deletes)
  m26_<mag_id>.assembly   AssemblyUtil.save_assemblies_from_fastas (release 3.0.0), from Blobstore nodes of the FASTAs
  m26_<mag_id>.RAST       RAST_SDK/annotate_genome_assembly -> RAST_SDK.rast_genome_assembly (release 1.9.5), one job
                          per genome: default annotation workflow, domain A/B and scientific name from GTDB R232,
                          genetic code 11, create_report 0 (the narrative sets 1, which would also save a report object)
Local
  ledger.json             per genome: Blobstore node, upload job, Assembly ref + checks, RAST attempts, genome ref,
                          download checks, timestamps; resumable (re-running skips what is done)
  genomes/<mag_id>.RAST.json   the Workspace object data of each genome (as jplfaria/Cliff MAGs/RASTified_MAG_json)
  logs/                   EE2 logs of failed jobs and the pilot archaeon; ws snapshot taken before the first write

Safety
  * a write-ahead marker is saved before every submit; a marker without a job id is reconciled against EE2
    (check_workspace_jobs) before anything is resubmitted, so no output name is ever written by two jobs;
  * every output name is checked to be absent from ws 275645 right before its job is submitted;
  * a completed RAST job without output_genome_ref, or an output that is not version 1 of the expected name, fails;
  * a non-transient failure sets ledger['halt']: nothing new is submitted until a person looks;
  * the 55 non-pilot genomes are refused until `pilot-check --apply` has passed.

Usage (Python from ~/Documents/py_venv; KBUtilLib from ~/Documents/KBUtilLib/src; token ~/.kbase/token, never printed)
  python rast_driver.py plan
  python rast_driver.py advance --scope pilot [--apply]     # poll, verify, submit what is ready, download
  python rast_driver.py pilot-check [--apply]
  python rast_driver.py watch --scope auto --interval 600 --apply [--rast-cap 10]
        # auto: the pilot end to end, then pilot-check (automatic; halts on failure), then the other 55
  python rast_driver.py status
  python rast_driver.py readme
"""
import argparse, csv, fcntl, hashlib, http.client, json, os, re, socket, sys, time, traceback, urllib.error
import urllib.parse, urllib.request
from collections import Counter
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))                 # reassembly_2026/modelling_2026/rast
R = os.path.dirname(os.path.dirname(HERE))                         # reassembly_2026
sys.path.insert(0, os.path.join(R, "local"))
from kbrpc import rpc, token  # noqa: E402  read-only JSON-RPC (whitelisted methods) and the token reader

KBUTILLIB_SRC = os.path.expanduser("~/Documents/KBUtilLib/src")
WSID = 275645
WS_NAME = "freiburgermsu:narrative_1793847843029"
SPECIES_TSV = os.path.join(R, "comparison_rollout", "modelling_species.tsv")
LINEAGE_TSV = os.path.join(R, "comparison_rollout", "mags_three_way.tsv")
FASTA_DIR = {"rollout": os.path.join(R, "work", "rollout_mags_kept"), "coassembly": os.path.join(R, "work", "coasm_mags_kept")}
MD5_MANIFEST = {"rollout": os.path.join(R, "rollout_final", "mags_kept.md5"),
                "coassembly": os.path.join(R, "coasm_final", "mags_kept.md5")}
LEDGER = os.path.join(HERE, "ledger.json")
GENOMES = os.path.join(HERE, "genomes")
LOGS = os.path.join(HERE, "logs")
SNAPSHOT = os.path.join(LOGS, "ws275645_objects_before_first_write.json")
PILOT = ("R1_sediment__concoct.bin.036", "R1_A_D1__concoct.bin.007")   # a co-assembly bacterium, a rollout archaeon

BLOB = "https://kbase.us/services/shock-api"
UPLOAD = {"module": "AssemblyUtil", "method": "AssemblyUtil.save_assemblies_from_fastas", "version": "3.0.0",
          "tag": "release", "git_commit_hash": "59855bcb24014c3030e146b98c8d3e70abac1fa1"}
RAST = {"module": "RAST_SDK", "app_id": "RAST_SDK/annotate_genome_assembly", "method": "RAST_SDK.rast_genome_assembly",
        "version": "1.9.5", "tag": "release", "git_commit_hash": "7171090d87fccc8b7ecf1a1d02398995dcc2dd45"}
RAST_CAP = 10                 # live RAST jobs (bigmem is shared with the GROW pipeline)
RAST_TRIES = 3                # attempts per genome, transient failures only
UPLOAD_TRIES = 3
LIVE = {"created", "estimating", "queued", "running"}
TRANSIENT = ("URLError", "timed out", "TimeoutError", "no response from", "RemoteDisconnected", "ConnectionError",
             "ConnectionResetError", "ConnectionRefusedError", "ConnectionAbortedError", "BrokenPipeError",
             "Bad Gateway", "Service Unavailable", "Gateway Time-out", "Gateway Timeout", "HTTP Error 5",
             " Server Error: ", "ServerError: Unknown: 0.", "IncompleteRead", "ChunkedEncodingError", "SSLError",
             "Max retries exceeded", "Connection aborted", "Temporary failure", "502", "503", "504")


def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def log(msg):
    print(f"[{now()}] {msg}", flush=True)


def upa(i):
    return f"{i[6]}/{i[0]}/{i[4]}"


def transient(msg):
    return any(m in (msg or "") for m in TRANSIENT)


# ------------------------------------------------------------------------------------------------ inputs
def sci_name(classification, domain):
    """GTDB R232 species where named; else '<genus> sp.'; else '<lowest named taxon> bacterium|archaeon'."""
    r = {k: v for k, v in (x.split("__", 1) for x in classification.split(";")) if v}
    if r.get("s"):
        return r["s"]
    if r.get("g"):
        return f"{r['g']} sp."
    for k in ("f", "o", "c", "p"):
        if r.get(k):
            return f"{r[k]} {'archaeon' if domain == 'Archaea' else 'bacterium'}"
    return "Unknown species"


def parse_fasta(path):
    """AssemblyUtil's view: contig id = first header word, sequence upper-cased; assembly md5 = md5 of the sorted
    contig md5s joined by ','."""
    seqs, cur, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf).upper()
                cur = line[1:].split()[0]
                assert cur not in seqs, ("duplicate contig id", path, cur)
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf).upper()
    cm = {c: hashlib.md5(s.encode()).hexdigest() for c, s in seqs.items()}
    return {"n_contigs": len(seqs), "bp": sum(map(len, seqs.values())),
            "assembly_md5": hashlib.md5(",".join(sorted(cm.values())).encode()).hexdigest(),
            "contigs": {c: [len(seqs[c]), cm[c]] for c in seqs}}


def load_inputs():
    rows = list(csv.DictReader(open(SPECIES_TSV), delimiter="\t"))
    lin = {r["genome"]: r for r in csv.DictReader(open(LINEAGE_TSV), delimiter="\t")}
    man = {}
    for s, p in MD5_MANIFEST.items():
        for line in open(p):
            h, n = line.split()
            man[(s, n[:-3] if n.endswith(".fa") else n)] = h
    out = {}
    for r in rows:
        s, m = r["rep_genome"].split("|", 1)
        g = lin[r["rep_genome"]]
        assert g["domain"] == r["domain"] and g["domain"] in ("Bacteria", "Archaea"), r["rep_genome"]
        out[m] = {"rep_genome": r["rep_genome"], "set": s, "fasta": os.path.join(FASTA_DIR[s], f"{m}.fa"),
                  "manifest_md5": man[(s, m)], "domain": g["domain"], "rast_domain": g["domain"][0],
                  "gtdb_r232_classification": g["classification"],
                  "scientific_name": sci_name(g["classification"], g["domain"]),
                  "guilds": [k for k in ("mepn_consumer", "betaine_reducer", "methanogen", "methylamine_user")
                             if r[k] == "True"],
                  "assembly_name": f"m26_{m}.assembly", "genome_name": f"m26_{m}.RAST"}
    assert len(out) == 57, len(out)
    return out


# ------------------------------------------------------------------------------------------------ ledger
def load_ledger():
    if os.path.exists(LEDGER):
        return json.load(open(LEDGER))
    return {"workspace_id": WSID, "workspace_name": WS_NAME, "created_utc": now(), "apps": {"upload": UPLOAD, "rast": RAST},
            "rast_params_rule": ("rast_genome_assembly service params: object_ref (Assembly UPA), output_workspace "
                                 f"'{WS_NAME}', output_genome_name, genetic_code 11, domain 'A'/'B' (GTDB R232 domain), "
                                 "scientific_name (GTDB R232: species; else '<genus> sp.'; else '<lowest taxon> "
                                 "bacterium|archaeon'), create_report 0; no ncbi_taxon_id; default RAST workflow"),
            "pilot": {"ids": list(PILOT), "passed": False}, "halt": None, "genomes": {}, "upload_jobs": {}, "events": []}


def save_ledger(L):
    L["updated_utc"] = now()
    tmp = LEDGER + ".part"
    with open(tmp, "w") as fh:
        json.dump(L, fh, indent=1, sort_keys=False)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, LEDGER)


def event(L, msg):
    L["events"].append({"utc": now(), "msg": msg})
    log(msg)


def halt(L, why):
    if not L.get("halt"):
        L["halt"] = {"utc": now(), "why": why}
    event(L, f"HALT: {why}")


def init_genomes(L, inputs):
    for m, d in inputs.items():
        g = L["genomes"].setdefault(m, {})
        if "fasta_md5" not in g:
            raw = open(d["fasta"], "rb").read()
            md5 = hashlib.md5(raw).hexdigest()
            if md5 != d["manifest_md5"]:
                sys.exit(f"{m}: FASTA md5 {md5} != manifest {d['manifest_md5']}")
            p = parse_fasta(d["fasta"])
            g.update({k: v for k, v in d.items() if k != "manifest_md5"})
            g.update(fasta_md5=md5, fasta_bytes=len(raw), n_contigs=p["n_contigs"], bp=p["bp"],
                     assembly_md5_local=p["assembly_md5"], state="new", upload={}, rast={"attempts": []}, download={})


# ------------------------------------------------------------------------------------------------ KBase
class KB:
    """Reads go through kbrpc (read-only whitelist); writes (Blobstore upload, EE2 submit) only through this class."""
    _runner = _jobs = None

    def infos(self, names):
        out = {}
        for i in range(0, len(names), 100):
            part = names[i:i + 100]
            res = rpc("ws", "Workspace.get_object_info3", [{"objects": [{"wsid": WSID, "name": n} for n in part],
                                                            "ignoreErrors": 1, "includeMetadata": 1}])[0]["infos"]
            out.update(zip(part, res))
        return out

    def info(self, ref):
        return rpc("ws", "Workspace.get_object_info3", [{"objects": [{"ref": ref}], "includeMetadata": 1}])[0]["infos"][0]

    def states(self, job_ids):
        out = {}
        for i in range(0, len(job_ids), 50):
            res = rpc("ee2", "execution_engine2.check_jobs", [{"job_ids": job_ids[i:i + 50], "exclude_fields": ["job_input"]}])[0]
            for s in res.get("job_states", []):
                out[s["job_id"]] = s
        return out

    def ws_jobs(self):
        return rpc("ee2", "execution_engine2.check_workspace_jobs",
                   [{"workspace_id": WSID, "return_list": 1, "exclude_fields": ["job_output"]}])[0].get("job_states", [])

    def find_job(self, method, needle, since_utc, exclude=()):
        since_ms = int(datetime.fromisoformat(since_utc).timestamp() * 1000) - 120_000
        for j in sorted(self.ws_jobs(), key=lambda j: j.get("created") or 0, reverse=True):
            ji = j.get("job_input") or {}
            if (ji.get("method") == method and (j.get("created") or 0) >= since_ms and j["job_id"] not in exclude
                    and json.dumps(needle) in json.dumps(ji.get("params"))):
                return j
        return None

    def job_log(self, job_id):
        lines, skip = [], 0
        while True:
            r = rpc("ee2", "execution_engine2.get_job_logs", [{"job_id": job_id, "skip_lines": skip}])[0]
            got = r.get("lines") or []
            lines += [x.get("line", "") for x in got]
            skip += len(got)
            if not got or skip >= (r.get("last_line_number") or 0):
                return lines

    def blob_upload(self, path, name, want_md5):
        payload = open(path, "rb").read()
        assert hashlib.md5(payload).hexdigest() == want_md5, ("local md5 changed", path)
        for attempt in range(4):
            try:
                req = urllib.request.Request(f"{BLOB}/node?filename={urllib.parse.quote(name)}", data=payload, method="POST",
                                             headers={"Authorization": "OAuth " + token(),
                                                      "Content-Type": "application/octet-stream"})
                with urllib.request.urlopen(req, timeout=600) as r:
                    res = json.load(r)
                break
            except (urllib.error.URLError, OSError, http.client.HTTPException) as e:
                if (isinstance(e, urllib.error.HTTPError) and e.code < 500) or attempt == 3:
                    raise
                time.sleep(15 * (attempt + 1))
        f = res["data"].get("file") or {}
        got = (f.get("checksum") or {}).get("md5")
        if got != want_md5 or int(f.get("size") or -1) != len(payload):
            raise RuntimeError(f"Blobstore node {res['data'].get('id')} md5/size {got}/{f.get('size')} != local")
        return {"node": res["data"]["id"], "node_md5": got, "node_size": int(f["size"]), "node_file_name": f.get("name"),
                "blob_utc": now()}

    def _kbutil(self):
        if self._jobs is None:
            if KBUTILLIB_SRC not in sys.path:
                sys.path.insert(0, KBUTILLIB_SRC)
            from kbutillib.kb_app_runner.nms import NMSSpecCache
            from kbutillib.kb_app_runner.runner import AppRunner
            from kbutillib.kb_genome_utils import KBGenomeUtils
            from kbutillib.kb_job_utils import KBJobUtils
            u = KBGenomeUtils(kb_version="prod")
            self._jobs = KBJobUtils(u, kb_version="prod")
            self._nms = NMSSpecCache()
            self._runner = AppRunner(ws=u, kb_version="prod", job_store=self._jobs, nms_cache=self._nms)
        return self._jobs

    def submit_upload(self, params):
        jobs = self._kbutil()
        rec = jobs.run_job(UPLOAD["method"], [params], workspace_id=WSID, service_ver=UPLOAD["git_commit_hash"],
                           meta={"task": "m26 modelling upload"})
        return rec.job_id

    def submit_rast(self, params):
        self._kbutil()
        spec = self._nms.get(RAST["app_id"], None)
        if spec.method != RAST["method"] or spec.service_ver != RAST["git_commit_hash"]:
            raise SystemExit(f"NMS release of {RAST['app_id']} changed: {spec.method} {spec.service_ver}; "
                             "update RAST and record the change before continuing")
        h = self._runner.run_app(RAST["app_id"], [params], WSID, pin_version=RAST["git_commit_hash"],
                                 meta={"task": "m26 modelling RAST"})
        return h.job_id


def err_of(s):
    e = s.get("error")
    return str((e.get("message") if isinstance(e, dict) else None) or s.get("errormsg") or "")[:3000]


def save_log(kb, job_id, tag):
    try:
        lines = kb.job_log(job_id)
    except Exception as e:
        return f"(log unavailable: {e})"
    p = os.path.join(LOGS, f"{tag}.{job_id}.log")
    with open(p, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return os.path.relpath(p, HERE)


# ------------------------------------------------------------------------------------------------ verification
def verify_assembly(kb, g, ref):
    inc = ["/md5", "/num_contigs", "/dna_size", "/assembly_id", "/type", "/fasta_handle_ref",
           "/contigs/*/contig_id", "/contigs/*/length", "/contigs/*/md5"]
    o = rpc("ws", "Workspace.get_objects2", [{"objects": [{"ref": ref, "included": inc}]}])[0]["data"][0]
    i, d = o["info"], o["data"]
    loc = parse_fasta(g["fasta"])
    probs = []
    if i[1] != g["assembly_name"]:
        probs.append(f"name {i[1]}")
    if not i[2].startswith("KBaseGenomeAnnotations.Assembly"):
        probs.append(f"type {i[2]}")
    if i[4] != 1:
        probs.append(f"version {i[4]} (expected 1)")
    if d.get("md5") != loc["assembly_md5"]:
        probs.append(f"md5 {d.get('md5')} != local {loc['assembly_md5']}")
    if d.get("num_contigs") != loc["n_contigs"] or d.get("dna_size") != loc["bp"]:
        probs.append(f"contigs/bp {d.get('num_contigs')}/{d.get('dna_size')} != {loc['n_contigs']}/{loc['bp']}")
    oc = d.get("contigs") or {}
    if set(oc) != set(loc["contigs"]):
        probs.append("contig id sets differ")
    else:
        bad = [c for c, (ln, m5) in loc["contigs"].items() if oc[c]["length"] != ln or oc[c]["md5"] != m5]
        if bad:
            probs.append(f"{len(bad)} contigs differ in length/md5")
    return {"ref": upa(i), "type": i[2], "save_date": i[3], "md5": d.get("md5"), "num_contigs": d.get("num_contigs"),
            "dna_size": d.get("dna_size"), "fasta_handle_ref": d.get("fasta_handle_ref"), "problems": probs}


def genome_counts(d):
    F = d.get("features") or []
    prot = [f for f in F if f.get("protein_translation")]
    funcs = lambda f: [x for x in (f.get("functions") or []) if (x or "").strip()]
    with_fn = [f for f in prot if funcs(f)]
    hypo = [f for f in with_fn if all(x.strip().lower() == "hypothetical protein" for x in funcs(f))]
    sso_terms = sum(1 for f in F if (f.get("ontology_terms") or {}).get("SSO"))
    return {"contigs": d.get("num_contigs"), "dna_size": d.get("dna_size"), "features": len(F),
            "cdss": len(d.get("cdss") or []), "mrnas": len(d.get("mrnas") or []),
            "non_coding_features": len(d.get("non_coding_features") or []),
            "non_coding_types": dict(Counter(f.get("type") for f in d.get("non_coding_features") or [])),
            "protein_features": len(prot), "protein_with_function": len(with_fn), "hypothetical_only": len(hypo),
            "protein_with_non_hypothetical_function": len(with_fn) - len(hypo),
            "distinct_functions": len({x for f in with_fn for x in funcs(f)}),
            "features_with_SSO_terms": sso_terms, "SSO_present": len((d.get("ontologies_present") or {}).get("SSO") or {}),
            "ontology_events": [{k: e.get(k) for k in ("method", "method_version", "id", "timestamp")}
                                for e in d.get("ontology_events") or []],
            "domain": d.get("domain"), "scientific_name": d.get("scientific_name"), "genetic_code": d.get("genetic_code"),
            "assembly_ref": d.get("assembly_ref"), "gc_content": d.get("gc_content"), "id": d.get("id")}


def check_genome(g, c):
    p = []
    if c["domain"] != g["rast_domain"]:
        p.append(f"domain {c['domain']!r} != {g['rast_domain']!r}")
    if c["scientific_name"] != g["scientific_name"]:
        p.append(f"scientific_name {c['scientific_name']!r} != {g['scientific_name']!r}")
    if c["genetic_code"] != 11:
        p.append(f"genetic_code {c['genetic_code']}")
    if (c["assembly_ref"] or "").split(";")[-1] != g["upload"]["assembly_ref"]:
        p.append(f"assembly_ref {c['assembly_ref']} != {g['upload']['assembly_ref']}")
    if c["contigs"] != g["n_contigs"] or c["dna_size"] != g["bp"]:
        p.append(f"contigs/bp {c['contigs']}/{c['dna_size']} != {g['n_contigs']}/{g['bp']}")
    if not c["protein_features"]:
        p.append("no protein features")
    if not c["protein_with_non_hypothetical_function"]:
        p.append("no RAST functions")
    if not any(e.get("id") == "SSO" for e in c["ontology_events"]):
        p.append("no SSO ontology event")
    return p


# ------------------------------------------------------------------------------------------------ advance
def scope_ids(L, scope):
    if scope == "pilot":
        return list(PILOT)
    if not L["pilot"].get("passed"):
        sys.exit("scope 'all' refused: the pilot has not passed `pilot-check --apply`")
    return sorted(L["genomes"])


def snapshot_ws():
    objs, start = [], 0
    while True:
        res = rpc("ws", "Workspace.list_objects", [{"ids": [WSID], "minObjectID": start + 1, "limit": 5000, "showHidden": 1}])[0]
        if not res:
            return objs
        objs += [[i[0], i[1], i[2], i[4]] for i in res]
        start = max(i[0] for i in res)


def poll_uploads(kb, L):
    live = {j: u for j, u in L["upload_jobs"].items() if u["status"] in LIVE}
    if not live:
        return
    st = kb.states(list(live))
    for jid, u in live.items():
        s = st.get(jid)
        if not s:
            continue
        u["status"] = s["status"]
        if s["status"] in LIVE:
            continue
        u["finished_utc"] = now()
        if s["status"] == "completed":
            res = ((s.get("job_output") or {}).get("result") or [{}])[0].get("results") or []
            if len(res) != len(u["members"]):
                halt(L, f"upload job {jid}: {len(res)} results for {len(u['members'])} inputs")
                continue
            for m, r in zip(u["members"], res):
                g = L["genomes"][m]
                v = verify_assembly(kb, g, r["upa"])
                g["upload"].update(assembly_ref=v["ref"], verify=v, uploaded_utc=v["save_date"], verified_utc=now())
                if v["problems"]:
                    g["state"] = "failed"
                    halt(L, f"{m}: uploaded Assembly {v['ref']} fails checks: {v['problems']}")
                else:
                    g["state"] = "uploaded"
                    log(f"{m}: Assembly {v['ref']} verified ({v['num_contigs']} contigs, {v['dna_size']} bp, md5 match)")
        else:
            msg = err_of(s)
            u["error"] = msg
            u["log"] = save_log(kb, jid, "upload")
            existing = {n: i for n, i in kb.infos([L["genomes"][m]["assembly_name"] for m in u["members"]]).items() if i}
            if existing:
                halt(L, f"upload job {jid} {s['status']} but {len(existing)} assembly names now exist: {sorted(existing)[:3]}")
                continue
            tries = len(L["genomes"][u["members"][0]]["upload"].get("jobs", []))
            if (transient(msg) or s["status"] == "terminated") and tries < UPLOAD_TRIES:
                for m in u["members"]:
                    L["genomes"][m]["state"] = "blob"
                event(L, f"upload job {jid} {s['status']} (transient: {msg[:200]}); members will be resubmitted")
            else:
                for m in u["members"]:
                    L["genomes"][m]["state"] = "failed"
                halt(L, f"upload job {jid} {s['status']}: {msg[:400]}")


def submit_uploads(kb, L, ids, apply):
    pend = [m for m in ids if L["genomes"][m]["state"] in ("new", "blob")]
    if not pend:
        return
    if L.get("upload_submitting"):
        mk = L["upload_submitting"]
        j = kb.find_job(UPLOAD["method"], L["genomes"][mk["members"][0]]["assembly_name"], mk["utc"],
                        exclude=set(L["upload_jobs"]))
        if j:
            adopt_upload(L, j["job_id"], mk["members"], mk["utc"], j.get("status", "queued"))
            event(L, f"adopted upload job {j['job_id']} of an interrupted submit")
            return
        L.pop("upload_submitting")
        event(L, "an interrupted upload submit left no EE2 job; submitting afresh")
    exist = {n: i for n, i in kb.infos([L["genomes"][m]["assembly_name"] for m in pend]).items() if i}
    if exist:
        halt(L, f"assembly names already exist in ws {WSID}, not overwriting: {sorted(exist)[:5]}")
        return
    if not apply:
        log(f"would upload {len(pend)} FASTAs to the Blobstore and submit one {UPLOAD['method']} job")
        return
    for m in pend:
        g = L["genomes"][m]
        if not g["upload"].get("node"):
            g["upload"].update(kb.blob_upload(g["fasta"], os.path.basename(g["fasta"]), g["fasta_md5"]))
            g["state"] = "blob"
            save_ledger(L)
            log(f"{m}: Blobstore node {g['upload']['node']} (md5 verified)")
    params = {"workspace_id": WSID,
              "inputs": [{"node": L["genomes"][m]["upload"]["node"], "assembly_name": L["genomes"][m]["assembly_name"]}
                         for m in pend]}
    t0 = now()
    L["upload_submitting"] = {"utc": t0, "members": pend}
    save_ledger(L)                                               # write-ahead
    try:
        jid = kb.submit_upload(params)
    except Exception as e:
        event(L, f"upload submit raised {type(e).__name__}: {str(e)[:300]} (marker kept; reconciled next round)")
        save_ledger(L)
        return
    adopt_upload(L, jid, pend, t0, "queued", params)
    event(L, f"submitted {UPLOAD['method']}@{UPLOAD['git_commit_hash'][:8]} for {len(pend)} assemblies -> {jid}")


def adopt_upload(L, jid, members, t0, status, params=None):
    L["upload_jobs"][jid] = {"members": members, "submitted_utc": t0, "status": status, "params": params}
    for m in members:
        g = L["genomes"][m]
        g["upload"]["job_id"] = jid
        g["upload"].setdefault("jobs", []).append(jid)
        g["state"] = "uploading"
    L.pop("upload_submitting", None)
    save_ledger(L)


def rast_params(g):
    return {"object_ref": g["upload"]["assembly_ref"], "output_workspace": WS_NAME, "genetic_code": 11,
            "domain": g["rast_domain"], "scientific_name": g["scientific_name"],
            "output_genome_name": g["genome_name"], "create_report": 0}


def poll_rast(kb, L):
    live = {a["job_id"]: m for m, g in L["genomes"].items() for a in g["rast"]["attempts"] if a["status"] in LIVE}
    if not live:
        return
    st = kb.states(list(live))
    for jid, m in live.items():
        s = st.get(jid)
        if not s:
            continue
        g = L["genomes"][m]
        a = next(x for x in g["rast"]["attempts"] if x["job_id"] == jid)
        if s["status"] != a["status"]:
            a["status"] = s["status"]
            if s["status"] == "running":
                a["running_utc"] = now()
        if s["status"] in LIVE:
            continue
        a["finished_utc"] = now()
        if s["status"] == "completed":
            out = ((s.get("job_output") or {}).get("result") or [{}])[0] or {}
            ref = out.get("output_genome_ref")
            if not ref:
                a["log"] = save_log(kb, jid, "rast")
                g["state"] = "failed"
                halt(L, f"{m}: RAST job {jid} completed without output_genome_ref (log {a['log']})")
                continue
            i = kb.info(ref)
            probs = []
            if i[1] != g["genome_name"]:
                probs.append(f"name {i[1]}")
            if not i[2].startswith("KBaseGenomes.Genome"):
                probs.append(f"type {i[2]}")
            if i[4] != 1:
                probs.append(f"version {i[4]}")
            if i[6] != WSID:
                probs.append(f"ws {i[6]}")
            a["output_ref"] = upa(i)
            g["rast"].update(output_ref=upa(i), type=i[2], save_date=i[3], job_id=jid, service_ver=None)
            if probs:
                g["state"] = "failed"
                halt(L, f"{m}: RAST output {upa(i)} unexpected: {probs}")
            else:
                g["state"] = "annotated"
                log(f"{m}: RAST {jid} -> {upa(i)}")
        else:
            msg = err_of(s)
            a["error"] = msg
            a["log"] = save_log(kb, jid, "rast")
            ex = kb.infos([g["genome_name"]])[g["genome_name"]]
            if ex:
                g["state"] = "failed"
                halt(L, f"{m}: RAST job {jid} {s['status']} but {g['genome_name']} exists ({upa(ex)})")
            elif (transient(msg) or s["status"] == "terminated") and len(g["rast"]["attempts"]) < RAST_TRIES:
                a["transient"] = True
                g["state"] = "uploaded"
                event(L, f"{m}: RAST job {jid} {s['status']} (transient: {msg[:200]}); will resubmit")
            else:
                g["state"] = "failed"
                halt(L, f"{m}: RAST job {jid} {s['status']}: {msg[:400]} (log {a['log']})")


def submit_rast(kb, L, ids, apply):
    live = sum(1 for g in L["genomes"].values() for a in g["rast"]["attempts"] if a["status"] in LIVE)
    for m in ids:
        g = L["genomes"][m]
        if g["state"] != "uploaded":
            continue
        mk = g["rast"].get("submitting")
        if mk:
            j = kb.find_job(RAST["method"], g["genome_name"], mk["utc"], exclude={a["job_id"] for a in g["rast"]["attempts"]})
            if j:
                g["rast"]["attempts"].append({"job_id": j["job_id"], "submitted_utc": mk["utc"], "status": j.get("status", "queued"),
                                              "params": mk.get("params"), "adopted": True})
                g["rast"].pop("submitting")
                g["state"] = "rast_active"
                live += 1
                event(L, f"{m}: adopted RAST job {j['job_id']} of an interrupted submit")
                save_ledger(L)
                continue
            g["rast"].pop("submitting")
            event(L, f"{m}: an interrupted RAST submit left no EE2 job; submitting afresh")
        if live >= RAST_CAP:
            break
        ex = kb.infos([g["genome_name"]])[g["genome_name"]]
        if ex:
            g["state"] = "failed"
            halt(L, f"{m}: {g['genome_name']} already exists ({upa(ex)}); not resubmitting")
            continue
        params = rast_params(g)
        if not apply:
            log(f"would submit RAST for {m}: {json.dumps(params)}")
            live += 1
            continue
        t0 = now()
        g["rast"]["submitting"] = {"utc": t0, "params": params}
        save_ledger(L)                                           # write-ahead
        try:
            jid = kb.submit_rast(params)
        except SystemExit:
            raise
        except Exception as e:
            event(L, f"{m}: RAST submit raised {type(e).__name__}: {str(e)[:300]} (marker kept; reconciled next round)")
            save_ledger(L)
            continue
        g["rast"]["attempts"].append({"job_id": jid, "submitted_utc": t0, "status": "queued", "params": params,
                                      "app_id": RAST["app_id"], "method": RAST["method"], "tag": RAST["tag"],
                                      "service_ver": RAST["git_commit_hash"]})
        g["rast"].pop("submitting")
        g["state"] = "rast_active"
        live += 1
        save_ledger(L)
        log(f"{m}: submitted {RAST['app_id']}@{RAST['git_commit_hash'][:8]} ({g['rast_domain']}, {g['scientific_name']!r}) -> {jid}")


def download(kb, L, ids, apply):
    for m in ids:
        g = L["genomes"][m]
        if g["state"] != "annotated":
            continue
        if not apply:
            log(f"would download {g['rast']['output_ref']}")
            continue
        o = rpc("ws", "Workspace.get_objects2", [{"objects": [{"ref": g["rast"]["output_ref"]}]}], timeout=900)[0]["data"][0]
        d = o["data"]
        prov = (o.get("provenance") or [{}])[0]
        c = genome_counts(d)
        probs = check_genome(g, c)
        path = os.path.join(GENOMES, f"{m}.RAST.json")
        tmp = path + ".part"
        with open(tmp, "w") as fh:
            json.dump(d, fh)
        raw = open(tmp, "rb").read()
        os.replace(tmp, path)
        g["download"] = {"path": os.path.relpath(path, HERE), "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(),
                         "utc": now(), "counts": c, "problems": probs,
                         "ws_info": o["info"][:10], "provenance": {k: prov.get(k) for k in ("service", "service_ver", "method",
                                                                                            "time", "epoch")}}
        g["rast"]["service_ver"] = prov.get("service_ver")
        if probs:
            g["state"] = "check_failed"
            halt(L, f"{m}: downloaded genome fails checks: {probs}")
        else:
            g["state"] = "downloaded"
            log(f"{m}: saved {g['download']['path']} ({c['features']} features, {c['protein_features']} proteins, "
                f"{c['protein_with_non_hypothetical_function']} with a non-hypothetical RAST function)")
        save_ledger(L)


def advance(L, kb, scope, apply):
    inputs = load_inputs()
    init_genomes(L, inputs)
    ids = scope_ids(L, scope)
    if apply and not os.path.exists(SNAPSHOT):
        snap = snapshot_ws()
        with open(SNAPSHOT, "w") as fh:
            json.dump({"utc": now(), "objects": snap}, fh)
        event(L, f"ws {WSID} snapshot before the first write: {len(snap)} objects")
    poll_uploads(kb, L)
    poll_rast(kb, L)
    if apply:
        save_ledger(L)
    download(kb, L, ids, apply)
    if L.get("halt"):
        log(f"halted ({L['halt']['why'][:300]}); polling and downloads only, nothing new submitted")
    else:
        submit_uploads(kb, L, ids, apply)
        submit_rast(kb, L, ids, apply)
    if apply:
        save_ledger(L)
    c = Counter(L["genomes"][m]["state"] for m in ids)
    log(f"scope {scope}: {dict(c)}")
    return c


def done(c):
    return all(k in ("downloaded", "failed", "check_failed") for k in c)


# ------------------------------------------------------------------------------------------------ pilot check
def pilot_check(L, kb, apply):
    out, ok = {}, True
    for m in PILOT:
        g = L["genomes"].get(m) or {}
        if g.get("state") != "downloaded":
            sys.exit(f"pilot {m} is {g.get('state')}, not downloaded")
        c = g["download"]["counts"]
        d = json.load(open(os.path.join(HERE, g["download"]["path"])))
        log_lines = kb.job_log(g["rast"]["job_id"])
        p = os.path.join(LOGS, f"pilot_rast.{m}.{g['rast']['job_id']}.log")
        open(p, "w").write("\n".join(log_lines) + "\n")
        dom_lines = [x.strip() for x in log_lines if "***** Domain" in x or "Scientific_name =" in x or "Genitic_code" in x]
        f_ok = [f for f in d["features"] if f.get("protein_translation") and f.get("functions")]
        ex = [{"id": f["id"], "functions": f["functions"], "location": f["location"][:1]} for f in f_ok[:3]]
        frac = c["protein_with_non_hypothetical_function"] / max(c["protein_features"], 1)
        res = {"state": g["state"], "assembly_ref": g["upload"]["assembly_ref"], "genome_ref": g["rast"]["output_ref"],
               "domain_recorded": c["domain"], "expected_domain": g["rast_domain"], "scientific_name": c["scientific_name"],
               "genetic_code": c["genetic_code"], "contigs": c["contigs"], "features": c["features"], "cdss": c["cdss"],
               "protein_features": c["protein_features"], "protein_with_function": c["protein_with_function"],
               "protein_with_non_hypothetical_function": c["protein_with_non_hypothetical_function"],
               "frac_non_hypothetical": round(frac, 3), "distinct_functions": c["distinct_functions"],
               "non_coding_types": c["non_coding_types"], "ontology_events": c["ontology_events"],
               "log_lines": dom_lines, "examples": ex, "problems": g["download"]["problems"],
               "job_log": os.path.relpath(p, HERE)}
        doms = [mm.group(1) for x in dom_lines for mm in [re.search(r"\*\*\*\*\* Domain\s*=\s*(\S+)", x)] if mm]
        res["log_domain"] = doms
        if res["problems"] or frac < 0.3 or not doms or any(x != g["rast_domain"] for x in doms):
            ok = False
        out[m] = res
    snap = json.load(open(SNAPSHOT))
    before = {o[0]: o for o in snap["objects"]}
    now_objs = snapshot_ws()
    new = [o for o in now_objs if o[0] not in before]
    changed = [o for o in now_objs if o[0] in before and o[3] != before[o[0]][3]]
    expected = {L["genomes"][m][k] for m in PILOT for k in ("assembly_name", "genome_name")}
    diff = {"new": new, "changed": changed, "expected_new": sorted(expected)}
    if sorted(o[1] for o in new) != sorted(expected) or any(o[3] != 1 for o in new) or \
            [o for o in changed if not o[2].startswith("KBaseNarrative.Narrative")]:
        ok = False
    out["ws_diff"] = diff
    print(json.dumps(out, indent=1))
    log(f"pilot check {'PASSED' if ok else 'FAILED'}")
    if apply:
        L["pilot"].update(passed=ok, checked_utc=now(), checks=out)
        if not ok:
            halt(L, "pilot check failed; see ledger['pilot']['checks']")
        save_ledger(L)
    return ok


# ------------------------------------------------------------------------------------------------ status / readme
def status(L):
    c = Counter(g["state"] for g in L["genomes"].values())
    live = [(m, a["job_id"], a["status"]) for m, g in L["genomes"].items() for a in g["rast"]["attempts"] if a["status"] in LIVE]
    print(f"genomes: {dict(c)}; live RAST jobs: {len(live)} {dict(Counter(s for _, _, s in live))}; "
          f"pilot passed: {L['pilot'].get('passed')}; halt: {L.get('halt')}")


def write_readme(L):
    G = L["genomes"]
    rows = sorted(G.items(), key=lambda kv: (kv[1]["domain"], kv[0]))
    done_ = [m for m, g in rows if g["state"] == "downloaded"]
    failed = [m for m, g in rows if g["state"] in ("failed", "check_failed")]
    retries = [(m, a) for m, g in rows for a in g["rast"]["attempts"] if a["status"] != "completed"]
    up_jobs = L["upload_jobs"]
    first = min((a["submitted_utc"] for g in G.values() for a in g["rast"]["attempts"]), default="")
    last = max((a.get("finished_utc", "") for g in G.values() for a in g["rast"]["attempts"]), default="")
    svc = Counter(g["rast"].get("service_ver") for g in G.values() if g["rast"].get("service_ver"))
    p = L["pilot"].get("checks") or {}
    tot = Counter()
    for m in done_:
        c = G[m]["download"]["counts"]
        for k in ("features", "cdss", "protein_features", "protein_with_function", "protein_with_non_hypothetical_function"):
            tot[k] += c[k]
    L_ = []
    w = L_.append
    w("# RAST genomes of the 57 modelling species (salterns 2026 re-assembly)\n")
    w(f"Written by `rast_driver.py readme` from `ledger.json` ({L.get('updated_utc')}). KBase narrative ws {WSID} "
      f"(`{WS_NAME}`, \"Salterns metagenomes: reassembly and binning (2026)\"). Write approved by Andrew on 2026-11-20 for "
      "these `m26_*` objects only.\n")
    w(f"**Status:** {len(done_)} of {len(G)} genomes annotated and downloaded; {len(failed)} failed.\n")
    if len(done_) + len(failed) < len(G):
        w(f"**In progress** (state counts {dict(Counter(g['state'] for g in G.values()))}). A detached watcher "
          "(`watch --scope auto --apply`, output in `watch.log`) polls KBase every 10 min: pilot end to end, automatic "
          "pilot check, then the other 55. Check with `python rast_driver.py status`; re-write this file with "
          "`python rast_driver.py readme`. To resume after a stop, start the same watch command again (the ledger "
          "is resumable; an interrupted submit is reconciled with EE2 before anything is resubmitted). KBase queues "
          "were saturated when this started (15-28 h waits on njs and bigmem for this account's jobs, 2026-11-20).\n")
    w("## Inputs\n")
    w("- Genomes: `../../comparison_rollout/modelling_species.tsv`, the species representatives of every MePn-consumer, "
      "betaine-reducer, methanogen or methylamine-user species of the final salterns MAG set "
      f"({Counter(g['set'] for g in G.values())['rollout']} rollout, {Counter(g['set'] for g in G.values())['coassembly']} "
      f"co-assembly; {Counter(g['domain'] for g in G.values())['Bacteria']} Bacteria, "
      f"{Counter(g['domain'] for g in G.values())['Archaea']} Archaea).")
    w("- FASTAs: `../../work/rollout_mags_kept/<mag_id>.fa` and `../../work/coasm_mags_kept/<mag_id>.fa`; each md5 equals "
      "`../../rollout_final/mags_kept.md5` / `../../coasm_final/mags_kept.md5`.")
    w("- Lineage: GTDB-Tk 2.8.0 / GTDB R232, `../../comparison_rollout/mags_three_way.tsv` (`domain`, `classification`).\n")
    w("## What was run\n")
    w(f"1. **Upload.** Each FASTA was POSTed to the KBase Blobstore (node md5 checked against the local file), then one EE2 job "
      f"per batch ran `{UPLOAD['method']}` (AssemblyUtil {UPLOAD['version']}, tag {UPLOAD['tag']}, git "
      f"`{UPLOAD['git_commit_hash']}`) with `workspace_id` {WSID} and per input `node` + `assembly_name` "
      "`m26_<mag_id>.assembly`; no `min_contig_length`, `type` left at AssemblyUtil's default `Unknown` (as the "
      "project's other MAG Assemblies). Batches: the 2 pilot genomes, then the other 55. Each saved Assembly was checked "
      "against its FASTA: version 1 of the expected name, contig ids, lengths and md5s, `num_contigs`, `dna_size`, and the "
      "AssemblyUtil assembly md5 (md5 of the sorted contig md5s). The Blobstore nodes of step 1 are left in place "
      "(AssemblyUtil stores its own copy); they hold only the FASTAs and are not referenced by any object.")
    w(f"2. **RAST.** One EE2 job per genome, app `{RAST['app_id']}` (\"Annotate Genome/Assembly with RASTtk - v1.073\"), "
      f"method `{RAST['method']}`, RAST_SDK {RAST['version']}, tag {RAST['tag']}, git `{RAST['git_commit_hash']}` "
      f"(pinned; the genomes' provenance records service_ver {dict(svc)}). Parameters, as the narrative app sends them "
      "except `create_report`:")
    w("   - `object_ref`: the verified Assembly UPA; `output_workspace`: the narrative's workspace name; "
      "`output_genome_name`: `m26_<mag_id>.RAST`;")
    w("   - `genetic_code`: 11 (GTDB-Tk reported table 11 for 23 genomes and none for the 34 placed by ANI screen; no "
      "genome uses 4 or 25);")
    w("   - `domain`: `A` for the 7 Archaea, `B` for the 50 Bacteria (the app's dropdown values; GTDB R232 domain);")
    w("   - `scientific_name`: from GTDB R232: the species name where GTDB names one (34 genomes), else `<genus> sp.` "
      "(22), else `<family> bacterium` (1). No `ncbi_taxon_id`, so RAST takes the name as given and no Relation Engine "
      "lookup is made;")
    w("   - `create_report`: 0. The narrative app sends 1, which also saves a `KBaseReport.Report` under a generated name; "
      "0 keeps the writes to the `m26_*` objects. It changes no annotation step;")
    w("   - annotation workflow: RAST_SDK's default for this app (gene calling, rRNA/tRNA/repeat/CRISPR/selenoprotein/"
      "pyrrolysoprotein calls, k-mer v2 and similarity-based functions, overlap resolution); this app exposes no "
      "workflow options.")
    w(f"   - Concurrency: at most {L.get('rast_cap', RAST_CAP)} live (queued + running) RAST jobs (client group bigmem, shared with the GROW pipeline). "
      f"First submitted {first}, last finished {last}.")
    w("3. **Download.** `Workspace.get_objects2` of each genome; `genomes/<mag_id>.RAST.json` is the object's `data` "
      "(features with locations, translations and `functions`, `cdss`, `mrnas`, `non_coding_features`, `ontology_events`, "
      "...), the same layout as `jplfaria/Cliff MAGs/RASTified_MAG_json/*.RAST.json`. Each file was checked: domain, "
      "scientific name, genetic code 11, `assembly_ref` = the uploaded Assembly, contig count and size = the FASTA, "
      "protein features present, RAST functions present, and an SSO ontology event.\n")
    w("## Pilot (1 bacterium, 1 archaeon, end to end before the rest)\n")
    if p:
        for m in PILOT:
            x = p.get(m) or {}
            w(f"- `{m}` ({G[m]['domain']}, {G[m]['scientific_name']}): Assembly `{x.get('assembly_ref')}`, genome "
              f"`{x.get('genome_ref')}`; domain recorded `{x.get('domain_recorded')}` (job log: "
              f"{'; '.join(x.get('log_lines') or [])}); {x.get('contigs')} contigs, {x.get('features')} features, "
              f"{x.get('protein_features')} protein features, {x.get('protein_with_function')} with a function, "
              f"{x.get('protein_with_non_hypothetical_function')} non-hypothetical "
              f"({100 * (x.get('frac_non_hypothetical') or 0):.0f} %), {x.get('distinct_functions')} distinct functions; "
              f"non-coding {x.get('non_coding_types')}; ontology events "
              f"{[(e.get('method'), e.get('id')) for e in x.get('ontology_events') or []]}.")
        d = p.get("ws_diff") or {}
        w(f"- Workspace diff before/after the pilot: new objects {sorted(o[1] for o in d.get('new', []))} (all version 1), "
          f"changed objects {[(o[1], o[3]) for o in d.get('changed', [])]}. Pilot check "
          f"{'passed' if L['pilot'].get('passed') else 'NOT passed'} ({L['pilot'].get('checked_utc')}).")
        w("- Baseline for comparison, two of the 2025 Cliff RAST genomes: 61-73 % of protein features with a function, "
          "10-13 % `hypothetical protein`, one SSO ontology event.\n")
    w("## Counts per genome\n")
    if done_:
        w(f"Totals over the {len(done_)} downloaded genomes: {tot['features']} features, {tot['cdss']} CDS, "
          f"{tot['protein_features']} protein features, {tot['protein_with_function']} with a RAST function, "
          f"{tot['protein_with_non_hypothetical_function']} non-hypothetical.\n")
    w("`fn` = protein features with a RAST function; `non-hyp` = of those, not only `hypothetical protein`; "
      "`ncRNA/repeat` = non-coding features of type rna / repeat.\n")
    w("| MAG | set | domain | scientific name (GTDB R232) | Assembly | RAST genome | contigs | features | CDS | fn | non-hyp | ncRNA/repeat | state |")
    w("|---|---|---|---|---|---|---:|---:|---:|---:|---:|---|---|")
    for m, g in rows:
        c = (g.get("download") or {}).get("counts") or {}
        nc = c.get("non_coding_types") or {}
        w(f"| {m} | {g['set']} | {g['rast_domain']} | {g['scientific_name']} | {g['upload'].get('assembly_ref', '')} | "
          f"{g['rast'].get('output_ref', '')} | {g['n_contigs']} | {c.get('features', '')} | {c.get('cdss', '')} | "
          f"{c.get('protein_with_function', '')} | {c.get('protein_with_non_hypothetical_function', '')} | "
          f"{nc.get('rna', 0) if c else ''}/{nc.get('repeat', 0) if c else ''} | {g['state']} |")
    w("\n## Jobs, failures and retries\n")
    for jid, u in up_jobs.items():
        w(f"- Upload job `{jid}` ({len(u['members'])} assemblies): {u['status']}, submitted {u['submitted_utc']}, "
          f"finished {u.get('finished_utc', '')}" + (f"; error: {u['error'][:300]}" if u.get("error") else ""))
    if retries:
        for m, a in retries:
            w(f"- RAST `{m}` job `{a['job_id']}`: {a['status']}" + (f" ({'transient, resubmitted' if a.get('transient') else 'not retried'})"
                                                                     if a["status"] not in LIVE else "")
              + (f"; error: {a.get('error', '')[:300]}" if a.get("error") else "") + (f"; log `{a['log']}`" if a.get("log") else ""))
    else:
        w("- RAST: no failed or retried jobs.")
    if failed:
        w(f"- Failed genomes: {failed}")
    if L.get("halt"):
        w(f"- HALT set {L['halt']['utc']}: {L['halt']['why']}")
    w("\n## Files\n")
    w("- `rast_driver.py`: the driver (plan / advance / pilot-check / watch / status / readme).")
    w("- `ledger.json`: per genome the FASTA md5, Blobstore node, upload job, Assembly ref and checks, RAST attempts "
      "(job ids, params, timestamps), genome ref, app version, download checks and counts; `events` log.")
    w("- `genomes/<mag_id>.RAST.json`: the 57 genome objects (Workspace `data`).")
    w("- `logs/`: workspace snapshot before the first write, EE2 logs of the pilot RAST jobs and of any failed job.")
    w("- `watch.log`: the detached watcher's output.")
    open(os.path.join(HERE, "README.md"), "w").write("\n".join(L_) + "\n")
    log("README.md written")


# ------------------------------------------------------------------------------------------------ main
class Lock:
    def __enter__(self):
        self.fh = open(LEDGER + ".lock", "w")
        try:
            fcntl.flock(self.fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            sys.exit("another rast_driver.py holds the ledger lock")
        return self

    def __exit__(self, *a):
        fcntl.flock(self.fh, fcntl.LOCK_UN)
        self.fh.close()


def main():
    global RAST_CAP
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan")
    s = sub.add_parser("advance"); s.add_argument("--scope", choices=("pilot", "all"), required=True); s.add_argument("--apply", action="store_true")
    s = sub.add_parser("watch"); s.add_argument("--scope", choices=("pilot", "all", "auto"), required=True)
    s.add_argument("--interval", type=int, default=600); s.add_argument("--apply", action="store_true")
    s.add_argument("--rast-cap", type=int, default=RAST_CAP, help="live (queued + running) RAST jobs")
    s = sub.add_parser("pilot-check"); s.add_argument("--apply", action="store_true")
    sub.add_parser("status")
    sub.add_parser("readme")
    a = ap.parse_args()
    os.makedirs(GENOMES, exist_ok=True)
    os.makedirs(LOGS, exist_ok=True)
    if a.cmd == "plan":
        inputs = load_inputs()
        names = [d[k] for d in inputs.values() for k in ("assembly_name", "genome_name")]
        ex = {n: i for n, i in KB().infos(names).items() if i}
        for m, d in sorted(inputs.items(), key=lambda kv: (kv[1]["domain"], kv[0])):
            print(f"{d['rast_domain']}  {m:50s} {d['set']:10s} {d['scientific_name']:40s} {','.join(d['guilds'])}")
        print(f"{len(inputs)} genomes; {Counter(d['domain'] for d in inputs.values())}; existing target names: {sorted(ex) or 'none'}")
        return
    if a.cmd in ("status", "readme"):                       # read-only: no lock, the ledger is replaced atomically
        L = load_ledger()
        status(L) if a.cmd == "status" else write_readme(L)
        return
    with Lock():
        L = load_ledger()
        kb = KB()
        if a.cmd == "advance":
            advance(L, kb, a.scope, a.apply)
        elif a.cmd == "watch":
            RAST_CAP = a.rast_cap
            L["rast_cap"] = RAST_CAP
            log(f"watcher pid {os.getpid()} on {socket.gethostname()}: scope {a.scope}, interval {a.interval} s, "
                f"RAST cap {RAST_CAP}, apply {a.apply}")
            while True:
                scope = a.scope if a.scope != "auto" else ("all" if L["pilot"].get("passed") else "pilot")
                try:
                    c = advance(L, kb, scope, a.apply)
                except SystemExit:
                    raise
                except Exception as e:
                    log(f"advance raised {type(e).__name__}: {str(e)[:500]}\n{traceback.format_exc()[-1500:]}")
                    c = {"error": 1}
                if a.scope == "auto" and scope == "pilot" and done(c) and "error" not in c:
                    if all(L["genomes"][m]["state"] == "downloaded" for m in PILOT):
                        try:
                            ok = pilot_check(L, kb, a.apply)
                        except Exception as e:
                            log(f"pilot check raised {type(e).__name__}: {str(e)[:500]}")
                            ok = None
                        if ok:
                            event(L, "pilot passed the automatic check; continuing with all 57")
                            save_ledger(L)
                            continue
                        if ok is None:
                            time.sleep(a.interval)
                            continue
                    else:
                        halt(L, "a pilot genome failed; the other 55 are not started")
                        save_ledger(L)
                    write_readme(L)
                    log("pilot did not pass; watcher exits for a person to look")
                    return
                if done(c) and "error" not in c and scope == ("all" if a.scope == "auto" else a.scope):
                    write_readme(L)
                    log("all genomes in scope are terminal; watcher exits")
                    return
                live = [1 for g in L["genomes"].values() for x in g["rast"]["attempts"] if x["status"] in LIVE] + \
                    [1 for u in L["upload_jobs"].values() if u["status"] in LIVE]
                if L.get("halt") and not live and not any(g["state"] == "annotated" for g in L["genomes"].values()):
                    write_readme(L)
                    log("halted and nothing live; watcher exits for a person to look")
                    return
                time.sleep(a.interval)
        elif a.cmd == "pilot-check":
            pilot_check(L, kb, a.apply)


if __name__ == "__main__":
    main()
