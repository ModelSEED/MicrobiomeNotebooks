#!/usr/bin/env python3
"""Write the local GTDB-Tk 2.8.0 / GTDB R232 classification of the 2026 salterns re-assembly MAGs into KBase
narrative 275645, the way KBase's kb_gtdbtk app records a lineage.

Adapted from ~/gtdbtk_kbase/writeback_r232.py (the wetland narrative 272312 writeback), minus its kb_gtdbtk
chunk logic (no kb_gtdbtk job ever ran in 275645).
  * each MAG Assembly that passed CheckM (>= 50 % complete, < 10 % contaminated) gets a new version whose
    data.std_lineages.gtdb = {lineage, source_ver "232.0", taxon_id}; the Workspace extracts GTDB_lineage /
    GTDB_source_ver / GTDB_taxon_id metadata from that field. Sequence data are unchanged, and the MD5 of the
    version written must equal the one GTDB-Tk classified (work/checkm/<unit>/selected_<unit>.tsv).
  * user metadata adds GTDB_classified_by; provenance records the GTDB-Tk command, versions and per-MAG details.
  * sets: each <unit>.DASTool_bins.AssemblySet is re-pointed at the new versions (kb_gtdbtk re-saves sets too).
  * report: a KBaseReport.Report with the summary tables, trees and run files, and an HTML summary.

Subcommands (every write is a dry run unless --apply):
  apply  [--apply] [--limit N]   new MAG versions; the first two saves are single-object canaries
  verify                         read-only: KBase vs local lineages and sequence identity
  sets   [--apply]               re-point the DAS Tool AssemblySets
  report [--apply]               upload files and save the report (only after apply + sets are complete)
Token: $KB_AUTH_TOKEN or ~/.kbase/token (never printed). Standard library only.
"""
import argparse, csv, glob, hashlib, html, json, os, socket, sys, tempfile, time, urllib.error, urllib.parse
import urllib.request, zipfile
from collections import Counter

R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))           # salterns/reassembly_2026
WSID = 275645
WS, HS = "https://kbase.us/services/ws", "https://kbase.us/services/handle_service"
OUT = os.path.expanduser("~/gtdbtk_kbase/salterns_reassembly_2026/gtdbtk_out")
BATCH = os.path.expanduser("~/gtdbtk_kbase/salterns_reassembly_2026/batchfile.tsv")
WB = os.path.join(R, "writeback")
LOG = os.path.join(WB, "writeback_gtdbtk_r232_log.tsv")
LOG_COLS = ["time", "name", "old_ref", "new_ref", "new_lineage", "assembly_md5", "note"]
SRC_VER = "232.0"
CLASSIFIED_BY = "local GTDB-Tk 2.8.0, GTDB R232 package"
BLOB_PUBLIC = "https://kbase.us/services/shock-api"
BLOB_CANON = "https://services.kbase.us/services/shock-api"
GTDB_URL = ("https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/"
            "full_package/gtdbtk_r232_data.tar.gz")
EXTRACTED = {"GC content", "MD5", "Name", "Size", "N Contigs"} | {
    f"{s}_{f}" for s in ("GTDB", "NCBI") for f in ("lineage", "source_ver", "taxon_id", "source_id")}
SEQ = ("MD5", "Size", "N Contigs")
DETAIL_COLS = ("classification_method", "closest_genome_reference", "closest_genome_ani", "closest_genome_af",
               "closest_placement_reference", "closest_placement_ani", "red_value", "msa_percent",
               "translation_table", "note", "warnings")
_TOKEN = None


def token():
    global _TOKEN
    if _TOKEN is None:
        _TOKEN = os.environ.get("KB_AUTH_TOKEN", "").strip() or open(os.path.expanduser("~/.kbase/token")).read().strip()
    return _TOKEN


def rpc(url, method, params, tries=4):
    body = json.dumps({"version": "1.1", "id": "1", "method": method, "params": params}).encode()
    for k in range(tries):
        req = urllib.request.Request(url, data=body, headers={"Authorization": token(), "Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                res = json.load(r)
        except urllib.error.HTTPError as e:
            try:
                res = json.load(e)
            except Exception:
                res = None
        except (urllib.error.URLError, TimeoutError, socket.timeout, ConnectionError):
            res = None
        if res is not None:
            if res.get("error"):
                sys.exit(f"{method} failed: {str(res['error'].get('message') or res['error'])[:400]}")
            return res["result"]
        if k == tries - 1:
            sys.exit(f"{method}: no usable response after {tries} tries")
        time.sleep(10 * (k + 1))


def upa(i):
    return f"{i[6]}/{i[0]}/{i[4]}"


def now():
    return time.strftime("%Y-%m-%dT%H:%M:%S+0000", time.gmtime())


def script_ver():
    return hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()[:12]


def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def taxon_id(c):
    parts = [x for x in c.split(";") if "__" in x and x.split("__", 1)[1]]
    return parts[-1] if parts else ""


def expected(lineage):
    return {"GTDB_lineage": lineage, "GTDB_source_ver": SRC_VER, "GTDB_taxon_id": taxon_id(lineage)}


def carries(info, lineage):
    m = (info[10] if info else None) or {}
    return all(m.get(k) == v for k, v in expected(lineage).items()) and m.get("GTDB_classified_by") == CLASSIFIED_BY


def run_info():
    j = json.load(open(os.path.join(OUT, "gtdbtk.json")))
    if j.get("version") != "2.8.0" or j.get("database_version") != "r232":
        sys.exit(f"gtdbtk.json: {j.get('version')} / {j.get('database_version')}, expected 2.8.0 / r232")
    if any(s.get("status") != "completed" for s in j.get("steps") or []) or BATCH not in j.get("command_line", ""):
        sys.exit("gtdbtk.json: run incomplete or not over the salterns batchfile")
    return j


def summaries():
    rows = {}
    for f in ("classify/gtdbtk.bac120.summary.tsv", "classify/gtdbtk.ar53.summary.tsv"):
        for r in read_tsv(os.path.join(OUT, f)):
            rows[r["user_genome"]] = r
    return rows


def mags():
    """name -> {ref, md5, unit, completeness, contamination, classification, summary} for every CheckM-passing MAG."""
    summ, batch = summaries(), {l.split("\t")[1].strip() for l in open(BATCH) if l.strip()}
    out = {}
    for f in sorted(glob.glob(os.path.join(R, "work/checkm/*/selected_*.tsv"))):
        unit = os.path.basename(os.path.dirname(f))
        for r in read_tsv(f):
            if r.get("pass_50_10") != "1":
                continue
            n = r["assembly_name"]
            out[n] = {"ref": r["assembly_ref"], "md5": r["md5"], "unit": unit, "completeness": r["completeness"],
                      "contamination": r["contamination"], "summary": summ.get(n, {}),
                      "classification": summ.get(n, {}).get("classification", "")}
    if set(out) != batch:
        sys.exit(f"CheckM-passing MAGs ({len(out)}) != GTDB-Tk batchfile ({len(batch)})")
    return out


def infos_by_name(names):
    out = {}
    for i in range(0, len(names), 100):
        part = names[i:i + 100]
        res = rpc(WS, "Workspace.get_object_info3", [{"objects": [{"wsid": WSID, "name": n} for n in part],
                                                      "includeMetadata": 1}])[0]["infos"]
        out.update(zip(part, res))
    return out


def user_meta(old):
    meta = {k: v for k, v in (old or {}).items() if k not in EXTRACTED}
    meta["GTDB_classified_by"] = CLASSIFIED_BY
    return meta


def new_data(data, lineage):
    d = dict(data)
    sl = dict(d.get("std_lineages") or {})
    sl["gtdb"] = {"lineage": lineage, "source_ver": SRC_VER, "taxon_id": taxon_id(lineage)}
    d["std_lineages"] = sl
    return d


def provenance(info, m, host, ver, cmdline, gtdbtk_cmd):
    s = m["summary"]
    custom = {c: str(s.get(c, ""))[:900] for c in DETAIL_COLS if s.get(c) not in (None, "", "N/A")}
    custom.update(classified_assembly_version=m["ref"], checkm_completeness=m["completeness"],
                  checkm_contamination=m["contamination"])
    return [{
        "time": now(), "script": "writeback_gtdbtk_r232.py", "script_ver": ver, "script_command_line": cmdline,
        "method_params": [{"tool": "GTDB-Tk", "version": "2.8.0", "workflow": "classify_wf", "command": gtdbtk_cmd,
                           "reference_data": "GTDB R232 (gtdbtk_r232_data.tar.gz)", "min_perc_aa": 10,
                           "ani_screen": "skani (GTDB-Tk default)", "host": host}],
        "input_ws_objects": [upa(info)],
        "subactions": [{"name": "GTDB-Tk", "ver": "2.8.0", "code_url": "https://github.com/Ecogenomics/GTDBTk"},
                       {"name": "skani", "ver": "0.4.0", "code_url": "https://github.com/bluenote-1577/skani"},
                       {"name": "pplacer", "ver": "1.1.alpha19", "code_url": "https://github.com/matsen/pplacer"},
                       {"name": "Prodigal", "ver": "2.6.3", "code_url": "https://github.com/hyattpd/Prodigal"},
                       {"name": "HMMER", "ver": "3.4", "code_url": "http://hmmer.org/"}],
        "external_data": [{"resource_name": "Genome Taxonomy Database (GTDB)", "resource_url": "https://gtdb.ecogenomic.org/",
                           "resource_version": "R232", "data_url": GTDB_URL,
                           "description": "GTDB-Tk reference data package r232"}],
        "custom": custom,
        "description": ("GTDB lineage from GTDB-Tk 2.8.0 classify_wf with the GTDB R232 reference package, run "
                        f"locally on {host} over the 1,324 CheckM-passing MAGs of the 2026 salterns re-assembly; "
                        "written to data.std_lineages.gtdb (the GTDB_* metadata are extracted from it). "
                        "Sequence data unchanged.")}]


def append_log(rows):
    os.makedirs(WB, exist_ok=True)
    new = not os.path.exists(LOG)
    with open(LOG, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=LOG_COLS, delimiter="\t", extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerows(rows)
        fh.flush(); os.fsync(fh.fileno())


# ------------------------------------------------------------------ apply
def cmd_apply(a):
    j, M = run_info(), mags()
    ok = {n: m for n, m in M.items() if m["classification"] and not m["classification"].startswith("Unclassified")}
    infos = infos_by_name(sorted(ok))
    todo = [n for n in sorted(ok) if not carries(infos[n], ok[n]["classification"])]
    bad = [n for n in todo if (infos[n][10] or {}).get("MD5") != ok[n]["md5"] or upa(infos[n]) != ok[n]["ref"]]
    for n in bad[:10]:
        print(f"  SKIP {n}: latest {upa(infos[n])} MD5 {(infos[n][10] or {}).get('MD5')} vs classified "
              f"{ok[n]['ref']} MD5 {ok[n]['md5']}")
    todo = [n for n in todo if n not in bad]
    print(f"{len(M)} CheckM-passing MAGs; {len(ok)} with an R232 lineage; {len(todo)} to write; "
          f"{len(ok) - len(todo) - len(bad)} already written; {len(bad)} skipped (latest version != classified)")
    if a.limit:
        todo = todo[:a.limit]
    batches = [b for b in [todo[:1], todo[1:2]] + [todo[i:i + 25] for i in range(2, len(todo), 25)] if b]
    host, ver = socket.gethostname(), script_ver()
    cmdline = "python -I -B writeback_gtdbtk_r232.py " + " ".join(sys.argv[1:])
    if not a.apply:
        for n in todo[:2]:
            print(f"\n== {n} ({upa(infos[n])})\nstd_lineages.gtdb after:",
                  json.dumps(new_data({}, ok[n]["classification"])["std_lineages"]["gtdb"]),
                  "\nuser meta:", json.dumps(user_meta(infos[n][10])),
                  "\nprovenance:", json.dumps(provenance(infos[n], ok[n], host, ver, cmdline, j["command_line"]))[:1200])
        print(f"\nDry run: would save {len(todo)} versions in {len(batches)} batches. Re-run with --apply.")
        return
    saved = 0
    for b, part in enumerate(batches):
        fresh = infos_by_name(part)
        jumped = [n for n in part if fresh[n][4] != infos[n][4]]
        if jumped:
            sys.exit(f"STOPPING: {jumped[:5]} changed since planning; re-run apply")
        data = rpc(WS, "Workspace.get_objects2", [{"objects": [{"ref": upa(fresh[n])} for n in part]}])[0]["data"]
        objs = []
        for n, d in zip(part, data):
            i = fresh[n]
            if d["info"][1] != n or d["data"].get("md5") != ok[n]["md5"]:
                sys.exit(f"object mismatch for {n}: {upa(d['info'])} md5 {d['data'].get('md5')}")
            objs.append({"type": i[2], "objid": i[0], "data": new_data(d["data"], ok[n]["classification"]),
                         "meta": user_meta(i[10]),
                         "provenance": provenance(i, ok[n], host, ver, cmdline, j["command_line"])})
        try:
            res = rpc(WS, "Workspace.save_objects", [{"id": WSID, "objects": objs}], tries=1)[0]
        except BaseException as e:
            after = infos_by_name(part)
            done = [n for n in part if after[n][4] > fresh[n][4] and carries(after[n], ok[n]["classification"])]
            append_log([{"time": now(), "name": n, "old_ref": upa(fresh[n]), "new_ref": upa(after[n]),
                         "new_lineage": ok[n]["classification"], "assembly_md5": ok[n]["md5"],
                         "note": "reconciled after failed save"} for n in done])
            sys.exit(f"save_objects failed ({e!r}); {len(done)}/{len(part)} committed and logged. Re-run apply.")
        rows, problems = [], []
        for n, ni in zip(part, res):
            i = fresh[n]
            if ni[1] != n or ni[0] != i[0] or ni[4] != i[4] + 1:
                problems.append(f"{n}: saved as {upa(ni)} (from {upa(i)})")
            if not carries(ni, ok[n]["classification"]):
                problems.append(f"{n}: metadata " + json.dumps({k: (ni[10] or {}).get(k) for k in expected("x")}))
            if any((ni[10] or {}).get(k) != (i[10] or {}).get(k) for k in SEQ):
                problems.append(f"{n}: sequence metadata changed")
            rows.append({"time": ni[3], "name": n, "old_ref": upa(i), "new_ref": upa(ni),
                         "new_lineage": ok[n]["classification"], "assembly_md5": (ni[10] or {}).get("MD5", ""), "note": ""})
        append_log(rows)
        saved += len(part)
        if problems:
            sys.exit("STOPPING, stored versions do not match:\n  " + "\n  ".join(problems[:10]))
        print(f"  {'canary ' if b < 2 else ''}saved {saved}/{len(todo)} (e.g. {part[0]} -> {upa(res[0])})", flush=True)
    print(f"Done: {saved} new versions; log {LOG}. Next: verify, then sets.")


# ------------------------------------------------------------------ verify
def cmd_verify(a):
    M = mags()
    ok = {n: m for n, m in M.items() if m["classification"] and not m["classification"].startswith("Unclassified")}
    infos = infos_by_name(sorted(ok))
    good = [n for n in ok if carries(infos[n], ok[n]["classification"]) and (infos[n][10] or {}).get("MD5") == ok[n]["md5"]]
    print(f"{len(good)}/{len(ok)} MAGs carry their R232 lineage with unchanged sequence (MD5) in KBase")
    for n in sorted(set(ok) - set(good))[:10]:
        print("  not yet:", n, upa(infos[n]), (infos[n][10] or {}).get("GTDB_lineage"))
    return len(good) == len(ok)


# ------------------------------------------------------------------ sets
def unit_sets():
    units = sorted({os.path.basename(os.path.dirname(f)) for f in glob.glob(os.path.join(R, "work/checkm/*/selected_*.tsv"))})
    return [f"{u}.DASTool_bins.AssemblySet" for u in units]


def latest_by_item(items):
    """item index -> latest info of the item's object (items have empty labels and version-pinned refs)."""
    refs = ["/".join(it["ref"].split("/")[:2]) for it in items]
    out = {}
    for i in range(0, len(refs), 100):
        res = rpc(WS, "Workspace.get_object_info3", [{"objects": [{"ref": r} for r in refs[i:i + 100]],
                                                      "includeMetadata": 1}])[0]["infos"]
        out.update({i + k: v for k, v in enumerate(res)})
    return out


def set_moves(items, latest, ok):
    out = []
    for j, it in enumerate(items):
        i = latest.get(j)
        n = i[1] if i else None
        if i and n in ok and carries(i, ok[n]["classification"]) and it["ref"] != upa(i) \
                and it["ref"].split("/")[:2] == upa(i).split("/")[:2]:
            out.append(j)
    return out


def cmd_sets(a):
    M = mags()
    ok = {n: m for n, m in M.items() if m["classification"] and not m["classification"].startswith("Unclassified")}
    host, ver, total, saved = socket.gethostname(), script_ver(), 0, []
    for nm in unit_sets():
        d = rpc(WS, "Workspace.get_objects2", [{"objects": [{"wsid": WSID, "name": nm}]}])[0]["data"][0]
        items = d["data"]["items"]
        latest = latest_by_item(items)
        moves = set_moves(items, latest, ok)
        total += len(moves)
        if not moves:
            continue
        print(f"  {nm} ({upa(d['info'])}): {len(moves)}/{len(items)} items to re-point")
        if a.apply:
            new_items = [dict(it, ref=upa(latest[j])) if j in moves else dict(it) for j, it in enumerate(items)]
            prov = [{"time": now(), "script": "writeback_gtdbtk_r232.py sets", "script_ver": ver,
                     "script_command_line": "python -I -B writeback_gtdbtk_r232.py " + " ".join(sys.argv[1:]),
                     "input_ws_objects": [upa(d["info"])],
                     "description": ("Items re-pointed at the MAG versions carrying the GTDB R232 lineage from a "
                                     f"local GTDB-Tk 2.8.0 run on {host} (as kb_gtdbtk re-saves the set).")}]
            info = rpc(WS, "Workspace.save_objects", [{"id": WSID, "objects": [
                {"type": d["info"][2], "objid": d["info"][0], "data": dict(d["data"], items=new_items),
                 "provenance": prov}]}], tries=1)[0][0]
            if info[4] != d["info"][4] + 1:
                print(f"     WARNING {nm}: version jumped {d['info'][4]}->{info[4]}")
            saved.append(upa(info))
            print(f"     saved {upa(info)}")
    print(f"{total} set items {'re-pointed' if a.apply else 'would be re-pointed (dry run; --apply to save)'}")
    if saved:
        with open(os.path.join(WB, "sets_saved.txt"), "a") as fh:
            fh.write("\n".join(saved) + "\n")


# ------------------------------------------------------------------ report
def blob_upload(path, name):
    payload = open(path, "rb").read()
    md5 = hashlib.md5(payload).hexdigest()
    for attempt in range(3):
        try:
            req = urllib.request.Request(f"{BLOB_PUBLIC}/node?filename={urllib.parse.quote(name)}", data=payload,
                                         method="POST", headers={"Authorization": "OAuth " + token(),
                                                                 "Content-Type": "application/octet-stream"})
            with urllib.request.urlopen(req, timeout=1800) as r:
                res = json.load(r)
            break
        except (urllib.error.URLError, OSError) as e:
            if attempt == 2:
                sys.exit(f"blob upload failed for {name}: {e}")
            time.sleep(15 * (attempt + 1))
    node = res["data"]["id"]
    if ((res["data"].get("file") or {}).get("checksum") or {}).get("md5") != md5:
        sys.exit(f"blob upload md5 mismatch for {name}")
    hid = rpc(HS, "AbstractHandle.persist_handle", [{"id": node, "type": "shock", "url": BLOB_CANON,
                                                     "file_name": name, "remote_md5": md5}], tries=1)[0]
    return hid, node


def ranks(c):
    return dict(x.split("__", 1) for x in c.split(";") if "__" in x)


def cmd_report(a):
    j, M = run_info(), mags()
    ok = {n: m for n, m in M.items() if m["classification"] and not m["classification"].startswith("Unclassified")}
    infos = infos_by_name(sorted(ok))
    current = [n for n in sorted(ok) if carries(infos[n], ok[n]["classification"])]
    sets = [l.strip() for l in open(os.path.join(WB, "sets_saved.txt"))] if os.path.exists(os.path.join(WB, "sets_saved.txt")) else []
    print(f"{len(current)}/{len(ok)} MAGs carry R232 in KBase; {len(sets)} sets re-pointed")
    if a.apply and len(current) != len(ok):
        sys.exit("apply has not finished; run apply and verify first")
    stage = tempfile.mkdtemp(prefix="report_stage_", dir=WB if os.path.isdir(WB) else None)
    tpath = os.path.join(stage, "salterns_reassembly_MAGs_gtdb_r232_checkm.tsv")
    cols = ["name", "kbase_ref", "unit", "checkm_completeness", "checkm_contamination", "classification",
            "classification_method", "closest_genome_reference", "closest_genome_ani", "red_value", "warnings"]
    with open(tpath, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(cols)
        for n in sorted(M):
            m, s = M[n], M[n]["summary"]
            w.writerow([n, upa(infos[n]) if n in infos else m["ref"], m["unit"], m["completeness"], m["contamination"],
                        m["classification"]] + [s.get(c, "") for c in cols[6:]])
    files = [(tpath, os.path.basename(tpath), "per-MAG GTDB R232 lineage, GTDB-Tk details, CheckM completeness/"
                                              "contamination and KBase object version")]
    for p in sorted(glob.glob(f"{OUT}/classify/*.tree") + glob.glob(f"{OUT}/classify/*.summary.tsv")):
        files.append((p, os.path.basename(p), f"{os.path.basename(p)} (GTDB-Tk 2.8.0, GTDB R232)"))
    zpath = os.path.join(stage, "GTDB-Tk_classify_wf_R232_salterns.zip")
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for p in sorted(glob.glob(f"{OUT}/**/*", recursive=True)):
            rel = os.path.relpath(p, OUT)
            if not os.path.isfile(p) or "intermediate_results" in rel.split(os.sep):
                continue
            if p.endswith(".msa.fasta.gz") and "user_msa" not in p:
                continue                                   # GTDB reference alignment: public, large
            z.write(p, os.path.join("gtdbtk_out", rel))
    files.append((zpath, os.path.basename(zpath), "GTDB-Tk classify_wf output (R232); intermediates and the public "
                                                  "GTDB reference alignment omitted"))
    good = [ok[n]["classification"] for n in current]
    dom = Counter(ranks(c).get("d", "") for c in good)
    phy = Counter(ranks(c).get("p", "") for c in good)
    nov = Counter("no genus" if not ranks(c).get("g") else "no species" if not ranks(c).get("s") else "species" for c in good)
    tab = lambda cnt, h: (f"<table><tr><th>{h}</th><th>MAGs</th></tr>" + "".join(
        f"<tr><td>{html.escape(x or '(none)')}</td><td>{v}</td></tr>" for x, v in cnt.most_common()) + "</table>")
    body = f"""<!doctype html><html><head><meta charset="utf-8"><title>GTDB-Tk R232 (local), salterns re-assembly</title>
<style>body{{font-family:sans-serif;margin:16px}}table{{border-collapse:collapse;margin:8px 0 20px}}
td,th{{border:1px solid #ccc;padding:3px 8px;text-align:left;font-size:13px}}</style></head><body>
<h2>GTDB-Tk 2.8.0 classify_wf, GTDB R232, run locally</h2>
<p>{len(current)} of {len(M)} CheckM-passing MAGs (&ge; 50 % complete, &lt; 10 % contaminated; KBase CheckM
lineage_wf, reduced tree) of the 2026 salterns re-assembly carry an R232 lineage in KBase
(<code>data.std_lineages.gtdb</code>, GTDB_source_ver {SRC_VER}); {len(sets)} DAS Tool AssemblySets re-pointed.
The MAGs are not dereplicated: per-sample assemblies and per-site sediment co-assemblies overlap.
Command: <code>{html.escape(j['command_line'])}</code>.</p>
<h3>Domains</h3>{tab(dom, 'Domain')}<h3>Lowest named rank</h3>{tab(nov, 'Novelty')}<h3>Phyla</h3>{tab(phy, 'Phylum')}
<h3>Per-MAG results</h3><table><tr><th>MAG</th><th>KBase version</th><th>CheckM compl./contam.</th><th>R232 lineage</th></tr>
{''.join(f"<tr><td>{html.escape(n)}</td><td>{upa(infos[n])}</td><td>{float(M[n]['completeness']):.1f} / {float(M[n]['contamination']):.1f}</td><td>{html.escape(ok[n]['classification'])}</td></tr>" for n in current)}
</table></body></html>"""
    hdir = os.path.join(stage, "html"); os.makedirs(hdir)
    open(os.path.join(hdir, "index.html"), "w").write(body)
    hzip = os.path.join(stage, "output.zip")
    with zipfile.ZipFile(hzip, "w", zipfile.ZIP_DEFLATED) as z:
        z.write(os.path.join(hdir, "index.html"), "index.html")
    name = "GTDBTk_R232_local_report_salterns_reassembly"
    print(f"report {name}: {len(files)} files + HTML, {sum(os.path.getsize(p) for p, _, _ in files) / 1e6:.1f} MB; staged {stage}")
    if not a.apply:
        print("Dry run. Re-run with --apply to upload and save.")
        return
    exist = rpc(WS, "Workspace.get_object_info3", [{"objects": [{"wsid": WSID, "name": name}], "ignoreErrors": 1}])[0]["infos"][0]
    if exist:
        sys.exit(f"{name} already exists ({upa(exist)})")
    links = []
    for p, n, d in files:
        hid, node = blob_upload(p, n)
        links.append({"handle": hid, "description": d, "name": n, "label": "", "URL": f"{BLOB_CANON}/node/{node}"})
        print(f"   uploaded {n}", flush=True)
    hid, node = blob_upload(hzip, "output.zip")
    rep = {"text_message": (f"GTDB-Tk 2.8.0 classify_wf with GTDB R232, run locally on {socket.gethostname()}: "
                            f"{len(current)} of {len(M)} CheckM-passing MAGs carry an R232 lineage; "
                            f"{len(sets)} DAS Tool AssemblySets re-pointed."),
           "warnings": [], "file_links": links, "direct_html_link_index": 0,
           "html_links": [{"handle": hid, "description": "GTDB-Tk Classify Report (R232, local)", "name": "index.html",
                           "label": "index.html", "URL": f"{BLOB_CANON}/node/{node}"}],
           "objects_created": [{"ref": upa(infos[n]), "description": "GTDB R232 lineage added (local GTDB-Tk 2.8.0)"}
                               for n in current] + [{"ref": s, "description": "Items re-pointed at the R232 MAG versions"} for s in sets]}
    prov = [{"time": now(), "script": "writeback_gtdbtk_r232.py report", "script_ver": script_ver(),
             "script_command_line": "python -I -B writeback_gtdbtk_r232.py " + " ".join(sys.argv[1:]),
             "input_ws_objects": sets,
             "description": "Report for the local GTDB-Tk 2.8.0 / GTDB R232 classification of the salterns re-assembly MAGs"}]
    info = rpc(WS, "Workspace.save_objects", [{"id": WSID, "objects": [
        {"type": "KBaseReport.Report", "name": name, "data": rep, "provenance": prov}]}], tries=1)[0][0]
    back = rpc(WS, "Workspace.get_objects2", [{"objects": [{"ref": upa(info)}]}])[0]["data"][0]["data"]
    hs = rpc(HS, "AbstractHandle.hids_to_handles", [[l["handle"] for l in back["file_links"] + back["html_links"]]])[0]
    print(f"Saved {name} as {upa(info)}: {len(back['objects_created'])} objects, {len(back['file_links'])} files, "
          f"{len(hs)} handles resolve")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("apply"); s.add_argument("--apply", action="store_true"); s.add_argument("--limit", type=int, default=0)
    sub.add_parser("verify")
    s = sub.add_parser("sets"); s.add_argument("--apply", action="store_true")
    s = sub.add_parser("report"); s.add_argument("--apply", action="store_true")
    a = p.parse_args()
    {"apply": cmd_apply, "verify": cmd_verify, "sets": cmd_sets, "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    main()
