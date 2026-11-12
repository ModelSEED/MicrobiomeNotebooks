#!/usr/bin/env python3
"""Move KBase MAGs that KBase's GTDB-Tk could not classify through a local GTDB-Tk run.

STATUS (2026-11-12): historical. Its list/fetch/collect steps were reused, but the wetland MAGs were classified with
GTDB-Tk 2.8.0 / GTDB R232 and written by ~/gtdbtk_kbase/writeback_r232.py (data.std_lineages.gtdb). Do NOT use
'writeback' here: it hard-codes GTDB-Tk 2.3.2 / R214 labels (GTDB_source_ver 214.1) and writes only user metadata,
which KBase overwrites from std_lineages. See gtdbtk_r232_and_kbutillib_skani.md.

Subcommands (run in order):
  list       find MAG Assembly objects in narrative 272312 that have no GTDB lineage yet
             (--controls N also picks N MAGs KBase already classified, as a known-answer test set)
  fetch      download FASTA files (md5-checked) for a list
  collect    read local GTDB-Tk output for a manifest; compares with KBase's lineages for control sets
  writeback  add the lineage to each MAG's KBase object metadata (dry run unless --apply)

Standard library only. The KBase token is read from $KB_AUTH_TOKEN or ~/.kbase_token and is never printed.
"""
import argparse, csv, glob, hashlib, json, os, random, re, socket, sys, time, urllib.error, urllib.request
from collections import Counter, defaultdict

WSID = 272312                                   # narrative 272312 ("wetland bins")
WS = "https://kbase.us/services/ws"
HS = "https://kbase.us/services/handle_service"
EE2 = "https://kbase.us/services/ee2"
GTDBTK_APP = "kb_gtdbtk/run_kb_gtdbtk_classify_wf"
ACTIVE = {"created", "estimating", "queued", "running"}
MAG_NAME = re.compile(r"__SRR\d+$")             # MAG objects are named <bin id>__<run accession>
CHUNK = re.compile(r"^MAG_pool\.part(\d\d)\.AssemblySet$")
LIST_TSV, MANIFEST, LINEAGES = "mags_to_classify.tsv", "genomes_manifest.tsv", "mag_gtdb_lineages.tsv"
CONTROLS_TSV = "controls_to_check.tsv"


def token():
    t = os.environ.get("KB_AUTH_TOKEN", "").strip()
    if not t:
        p = os.path.expanduser("~/.kbase_token")
        if os.path.exists(p):
            t = open(p).read().strip()
    if not t:
        sys.exit("No KBase token: set KB_AUTH_TOKEN or put it in ~/.kbase_token (chmod 600).")
    return t


def rpc(url, method, params, tries=4):
    """KBase JSON-RPC 1.1 call. Service errors stop the script; network/gateway errors are retried."""
    body = json.dumps({"version": "1.1", "id": "1", "method": method, "params": params}).encode()
    for k in range(tries):
        req = urllib.request.Request(url, data=body, headers={"Authorization": token(),
                                                              "Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                res = json.load(r)
        except urllib.error.HTTPError as e:      # KBase reports service errors as JSON with HTTP 500
            try:
                res = json.load(e)
            except Exception:
                res = None                         # gateway/proxy error: retry
        except (urllib.error.URLError, TimeoutError, socket.timeout, ConnectionError):
            res = None                             # network problem: retry
        if res is not None:
            if res.get("error"):
                err = res["error"]
                sys.exit(f"{method} failed: {str(err.get('message') or err)[:400]}")
            return res["result"]
        if k == tries - 1:
            sys.exit(f"{method}: no usable response from {url} after {tries} tries")
        time.sleep(10 * (k + 1))


def upa(i):
    return f"{i[6]}/{i[0]}/{i[4]}"


def list_objects(type_str, meta=False):
    out, minid = [], 1
    while True:
        batch = rpc(WS, "Workspace.list_objects", [{"ids": [WSID], "type": type_str, "minObjectID": minid,
                                                     "includeMetadata": int(meta), "limit": 5000}])[0]
        out += batch
        if len(batch) < 5000:
            return out
        minid = max(i[0] for i in batch) + 1


def write_tsv(path, rows, cols):
    if any("kbase_lineage" in r for r in rows) and "kbase_lineage" not in cols:
        cols = cols + ["kbase_lineage"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ------------------------------------------------------------------ list
def cmd_list(a):
    mags = {i[1]: i for i in list_objects("KBaseGenomeAnnotations.Assembly", meta=True) if MAG_NAME.search(i[1])}
    chunk_of = {}
    for s in list_objects("KBaseSets.AssemblySet"):
        if CHUNK.match(s[1]):
            items = rpc(WS, "Workspace.get_objects2", [{"objects": [{"ref": upa(s)}]}])[0]["data"][0]["data"]["items"]
            for it in items:
                chunk_of[it.get("label") or ""] = s[1]
    jobs = rpc(EE2, "execution_engine2.check_workspace_jobs",
               [{"workspace_id": WSID, "return_list": 1,
                 "exclude_fields": ["job_output", "authstrat", "condor_job_ads", "scheduler_type", "scheduler_id"]}])[0]
    jobs = jobs["job_states"] if isinstance(jobs, dict) else jobs
    by_chunk = defaultdict(Counter)
    for j in jobs:
        ji = j.get("job_input") or {}
        if ji.get("app_id") != GTDBTK_APP:
            continue
        base = (ji.get("params") or [{}])[0].get("output_tree_basename", "")
        st = j.get("status")
        by_chunk[base.replace(".GTDB_tree", ".AssemblySet")]["active" if st in ACTIVE else st] += 1
    rows, skipped = [], Counter()
    for name, i in sorted(mags.items()):
        if (i[10] or {}).get("GTDB_lineage"):
            continue                                    # KBase already classified it
        chunk = chunk_of.get(name, "")
        st = by_chunk.get(chunk, Counter())
        if st["active"] and not a.include_active:
            skipped[chunk] += 1
            continue
        rows.append({"name": name, "ref": upa(i), "chunk": chunk, "kbase_gtdbtk_failures": st["error"],
                     "size_bp": (i[10] or {}).get("Size", ""), "contigs": (i[10] or {}).get("N Contigs", "")})
    write_tsv(LIST_TSV, rows, ["name", "ref", "chunk", "kbase_gtdbtk_failures", "size_bp", "contigs"])
    print(f"{len(mags)} MAG objects in narrative {WSID}; "
          f"{sum(1 for i in mags.values() if (i[10] or {}).get('GTDB_lineage'))} already have a GTDB lineage.")
    print(f"{len(rows)} MAGs to classify locally -> {LIST_TSV}")
    for c, n in sorted(Counter(r["chunk"] for r in rows).items()):
        print(f"   {c or '(not in a MAG_pool chunk)'}: {n} MAGs, KBase GTDB-Tk failures: {by_chunk[c]['error']}")
    for c, n in sorted(skipped.items()):
        print(f"   skipped {c}: {n} MAGs, a KBase GTDB-Tk job for this chunk is still queued/running")
    if a.controls:
        done = sorted(n for n, i in mags.items() if (i[10] or {}).get("GTDB_lineage"))
        pick = sorted(random.Random(214).sample(done, min(a.controls, len(done))))
        write_tsv(CONTROLS_TSV, [{"name": n, "ref": upa(mags[n]), "chunk": chunk_of.get(n, ""),
                                  "kbase_lineage": mags[n][10]["GTDB_lineage"]} for n in pick],
                  ["name", "ref", "chunk", "kbase_lineage"])
        print(f"{len(pick)} already-classified MAGs picked as known-answer controls -> {CONTROLS_TSV}")


# ------------------------------------------------------------------ fetch
def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def cmd_fetch(a):
    rows = read_tsv(a.list)
    os.makedirs(a.genome_dir, exist_ok=True)
    manifest, done = [], 0
    for k in range(0, len(rows), 50):
        part = rows[k:k + 50]
        objs = rpc(WS, "Workspace.get_objects2", [{"objects": [{"ref": r["ref"], "included": ["/fasta_handle_ref"]}
                                                               for r in part]}])[0]["data"]
        hids = [o["data"]["fasta_handle_ref"] for o in objs]
        handles = rpc(HS, "AbstractHandle.hids_to_handles", [hids])[0]
        hmap = {h["hid"]: h for h in handles}
        for r, hid in zip(part, hids):
            h = hmap[hid]
            path = os.path.join(a.genome_dir, r["name"] + ".fna")
            if not (os.path.exists(path) and md5sum(path) == h.get("remote_md5")):
                url = f'{h["url"]}/node/{h["id"]}?download_raw'
                tmp = path + ".part"
                for attempt in range(4):
                    try:
                        req = urllib.request.Request(url, headers={"Authorization": "OAuth " + token()})
                        with urllib.request.urlopen(req, timeout=600) as resp, open(tmp, "wb") as out:
                            while True:
                                block = resp.read(1 << 20)
                                if not block:
                                    break
                                out.write(block)
                        if h.get("remote_md5") and md5sum(tmp) != h["remote_md5"]:
                            raise IOError("md5 mismatch")
                        os.replace(tmp, path)
                        break
                    except Exception as e:
                        if attempt == 3:
                            sys.exit(f"download failed for {r['name']}: {e}")
                        time.sleep(10 * (attempt + 1))
            with open(path) as fh:
                if not fh.read(1) == ">":
                    sys.exit(f"{path} is not FASTA")
            manifest.append(dict(r, file=path, md5=md5sum(path)))
            done += 1
        print(f"  {done}/{len(rows)} MAGs downloaded and md5-checked")
    write_tsv(a.manifest, manifest, ["name", "ref", "chunk", "file", "md5", "size_bp", "contigs"])
    print(f"Wrote {len(manifest)} FASTA files to {a.genome_dir}/ and {a.manifest}")


# ------------------------------------------------------------------ collect
def first_real(r, *keys):
    return next((r[k] for k in keys if r.get(k) not in (None, "", "N/A")), "")


def cmd_collect(a):
    want = {r["name"]: r for r in read_tsv(a.manifest)}
    got, files = {}, sorted(set(os.path.realpath(p) for p in
                                glob.glob(os.path.join(a.out_dir, "**", "gtdbtk.*.summary.tsv"), recursive=True)))
    files = [f for f in files if re.search(r"gtdbtk\.(bac120|ar53|ar122)\.summary\.tsv$", f)]
    if not files:
        sys.exit(f"No gtdbtk.*.summary.tsv under {a.out_dir}; did classify_wf finish?")
    for f in files:
        for r in read_tsv(f):
            g = r.get("user_genome")
            if g in want and r.get("classification"):
                got[g] = r
    rows, missing = [], []
    for name, m in sorted(want.items()):
        r = got.get(name)
        if not r:
            missing.append(name)
            continue
        row = {"name": name, "ref": m["ref"], "chunk": m["chunk"], "classification": r["classification"],
               "classification_method": r.get("classification_method", ""),
               "closest_reference": first_real(r, "fastani_reference", "closest_placement_reference"),
               "ani": first_real(r, "fastani_ani", "closest_placement_ani"),
               "red_value": r.get("red_value", ""), "warnings": r.get("warnings", ""), "note": r.get("note", "")}
        if m.get("kbase_lineage"):
            row["kbase_lineage"] = m["kbase_lineage"]
        rows.append(row)
    write_tsv(a.out, rows, ["name", "ref", "chunk", "classification", "classification_method", "closest_reference",
                            "ani", "red_value", "warnings", "note"])
    ok = [r for r in rows if not r["classification"].startswith("Unclassified")]
    ranks = lambda c: dict(x.split("__", 1) for x in c.split(";") if "__" in x)
    print(f"Read {len(files)} summary file(s). {len(rows)}/{len(want)} MAGs in the output, {len(ok)} with a lineage.")
    print("Domains:", dict(Counter(ranks(r["classification"]).get("d", "") for r in ok)))
    print("Phyla:", dict(Counter(ranks(r["classification"]).get("p", "") for r in ok).most_common()))
    nov = Counter("no genus" if not ranks(r["classification"]).get("g") else
                  "no species" if not ranks(r["classification"]).get("s") else "species" for r in ok)
    print("Novelty:", dict(nov))
    if missing:
        print(f"{len(missing)} MAGs are not in any summary (filtered for too few markers, or failed). "
              f"See {a.out_dir}/align/*filtered.tsv and {a.out_dir}/identify/*failed*: {', '.join(missing[:10])}"
              + (" ..." if len(missing) > 10 else ""))
    ctrl = [r for r in rows if r.get("kbase_lineage")]
    if ctrl:                                  # known-answer controls: compare with KBase's GTDB-Tk lineages
        exact = [r for r in ctrl if r["classification"] == r["kbase_lineage"]]
        genus = [r for r in ctrl if ranks(r["classification"]).get("g") == ranks(r["kbase_lineage"]).get("g")
                 and ranks(r["classification"]).get("p") == ranks(r["kbase_lineage"]).get("p")]
        print(f"CONTROLS: {len(exact)}/{len(ctrl)} identical to KBase's lineage; {len(genus)}/{len(ctrl)} agree "
              "on phylum and genus.")
        for r in ctrl:
            if r["classification"] != r["kbase_lineage"]:
                print(f"   differs: {r['name']}\n      local: {r['classification']}\n      KBase: {r['kbase_lineage']}")
    print(f"Wrote {a.out}")


# ------------------------------------------------------------------ writeback
def taxon_id(classification):
    parts = [x for x in classification.split(";") if "__" in x and x.split("__", 1)[1]]
    return parts[-1] if parts else ""


def cmd_writeback(a):
    rows = [r for r in read_tsv(a.lineages)
            if r["classification"] and not r["classification"].startswith("Unclassified") and not r.get("kbase_lineage")]
    host, now = socket.gethostname(), time.strftime("%Y-%m-%dT%H:%M:%S+0000", time.gmtime())
    plan, skip = [], []
    for k in range(0, len(rows), 25):
        part = rows[k:k + 25]
        names = [{"wsid": WSID, "name": r["name"]} for r in part]    # latest version of each MAG
        infos = rpc(WS, "Workspace.get_object_info3", [{"objects": names, "includeMetadata": 1}])[0]["infos"]
        for r, info in zip(part, infos):
            if (info[10] or {}).get("GTDB_lineage"):
                skip.append(r["name"])                           # KBase classified it in the meantime
            else:
                plan.append((r, info))
    print(f"{len(plan)} MAGs would get a GTDB lineage; {len(skip)} skipped (already have one in KBase).")
    for r, info in plan[:3]:
        print(f"   {r['name']} ({upa(info)}): {r['classification']}")
    if not a.apply:
        print("Dry run only. Re-run with --apply after Andrew approves writing to KBase.")
        return
    saved = 0
    for k in range(0, len(plan), 25):
        part = plan[k:k + 25]
        data = rpc(WS, "Workspace.get_objects2", [{"objects": [{"ref": upa(i)} for _, i in part]}])[0]["data"]
        objs = []
        for (r, info), d in zip(part, data):
            meta = dict(info[10] or {})
            meta.update({"GTDB_lineage": r["classification"], "GTDB_taxon_id": taxon_id(r["classification"]),
                         "GTDB_source_ver": "214.1", "GTDB_classified_by": "local GTDB-Tk 2.3.2, GTDB R214 package"})
            objs.append({"type": info[2], "objid": info[0], "data": d["data"], "meta": meta,
                         "provenance": [{"time": now, "script": "kbase_mags.py writeback", "script_ver": "1",
                                         "input_ws_objects": [upa(info)],
                                         "description": f"GTDB-Tk 2.3.2 classify_wf (GTDB R214) run locally on {host}; "
                                                        "lineage added to object metadata"}]})
        rpc(WS, "Workspace.save_objects", [{"id": WSID, "objects": objs}])
        saved += len(objs)
        print(f"  saved {saved}/{len(plan)}")
    check = rpc(WS, "Workspace.get_object_info3",
                [{"objects": [{"wsid": WSID, "name": r["name"]} for r, _ in plan[:50]], "includeMetadata": 1}])[0]["infos"]
    print(f"Verified: {sum(1 for i in check if (i[10] or {}).get('GTDB_lineage'))}/{len(check)} checked objects "
          "now carry GTDB_lineage.")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("list")
    s.add_argument("--include-active", action="store_true",
                   help="also list MAGs whose chunk still has a KBase GTDB-Tk job queued or running")
    s.add_argument("--controls", type=int, default=0, metavar="N",
                   help=f"also write N already-classified MAGs to {CONTROLS_TSV} as a known-answer test")
    s = sub.add_parser("fetch")
    s.add_argument("--list", default=LIST_TSV)
    s.add_argument("--genome-dir", default="genomes")
    s.add_argument("--manifest", default=MANIFEST)
    s = sub.add_parser("collect")
    s.add_argument("--manifest", default=MANIFEST)
    s.add_argument("--out-dir", default="gtdbtk_out", help="GTDB-Tk classify_wf --out_dir")
    s.add_argument("--out", default=LINEAGES)
    s = sub.add_parser("writeback")
    s.add_argument("--lineages", default=LINEAGES)
    s.add_argument("--apply", action="store_true", help="actually save new object versions in KBase")
    a = p.parse_args()
    {"list": cmd_list, "fetch": cmd_fetch, "collect": cmd_collect, "writeback": cmd_writeback}[a.cmd](a)


if __name__ == "__main__":
    main()
