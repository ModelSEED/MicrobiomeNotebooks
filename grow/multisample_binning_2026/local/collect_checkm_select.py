#!/usr/bin/env python
"""Collect a kb_Msuite run_checkM_lineage_wf report for an AssemblySet and pick MAGs for local GTDB-Tk (read-only).

Report (kb_Msuite 7495e02e CheckMUtil._build_output_packages, identical in release fae9ff52):
  file_links[name == 'CheckM_summary_table.tsv.zip'] -> zip member 'CheckM_summary_table.tsv'
      columns: Bin Name, Marker Lineage, # Genomes, # Markers, # Marker Sets, 0, 1, 2, 3, 4, 5+, Completeness,
      Contamination (last two rounded to 2 dp); one row per bin, Bin Name = the Assembly object name
  file_links[name == 'full_output.zip'] -> 'storage/bin_stats_ext.tsv': '<bin id>\\t<python dict literal>' per bin
      (full-precision Completeness / Contamination, Genome size, GC, N50 ...); used for the thresholds.
URL host services.kbase.us (unreachable from here) is rewritten to kbase.us; '?download' returns the zip.
Selection: Completeness >= 50 and Contamination < 10 (MIMAG medium+). Output: <out>/checkm_<tag>.tsv,
<out>/selected_<tag>.tsv and a GTDB-Tk --batchfile (FASTA path, genome id) whose FASTA is downloaded from each
selected Assembly's blob (fetch_assembly.py; md5-verified) unless --local-fasta-dir has a file whose AssemblyUtil md5
equals the Assembly's MD5.
Usage: python -I collect_checkm_select.py --report REF --set REF --tag NAME --out DIR [--local-fasta-dir DIR]
"""
import argparse, ast, csv, hashlib, io, json, os, sys, urllib.request, zipfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kbrpc import rpc, token
from fetch_assembly import fetch, _retry


def blob_get(url):
    """A report zip; transport errors and 5xx are retried (fetch_assembly._retry), a 4xx raises at once."""
    u = url.replace("://services.kbase.us/", "://kbase.us/") + "?download"

    def get():
        with urllib.request.urlopen(urllib.request.Request(u, headers={"Authorization": "OAuth " + token()}),
                                    timeout=300) as r:         # small zips; this runs inside the watcher's cycle
            return r.read()                                     # read() without a size raises IncompleteRead if short
    return _retry(get, "CheckM report zip")


def asm_md5(path):
    md5s, cur, buf = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if cur is not None:
                md5s.append(hashlib.md5("".join(buf).upper().encode()).hexdigest())
            cur, buf = line, []
        else:
            buf.append(line.strip())
    if cur is not None:
        md5s.append(hashlib.md5("".join(buf).upper().encode()).hexdigest())
    return hashlib.md5(",".join(sorted(md5s)).encode()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", required=True, help="KBaseReport UPA from the CheckM job's result[0].report_ref")
    ap.add_argument("--set", required=True, help="the AssemblySet UPA given to CheckM as input_ref")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--local-fasta-dir", default=None, help="dir of <bid>.fa (DAS Tool --write_bins) to reuse")
    ap.add_argument("--min-completeness", type=float, default=50.0)
    ap.add_argument("--max-contamination", type=float, default=10.0)
    a = ap.parse_args()
    os.makedirs(os.path.join(a.out, "fasta"), exist_ok=True)
    rep = rpc("ws", "Workspace.get_objects2", [{"objects": [{"ref": a.report}]}])[0]["data"][0]["data"]
    fl = {f["name"]: f for f in rep.get("file_links") or []}
    summ = zipfile.ZipFile(io.BytesIO(blob_get(fl["CheckM_summary_table.tsv.zip"]["URL"]))).read("CheckM_summary_table.tsv")
    rows = list(csv.DictReader(io.StringIO(summ.decode()), delimiter="\t"))
    full = zipfile.ZipFile(io.BytesIO(blob_get(fl["full_output.zip"]["URL"])))
    ext = {}
    for line in full.read("storage/bin_stats_ext.tsv").decode().splitlines():
        if line.strip():
            bid, js = line.split("\t", 1)
            ext[bid] = ast.literal_eval(js)
    head = rpc("ws", "Workspace.get_object_info3", [{"objects": [{"ref": a.set}]}])[0]["infos"][0]
    if head[2].startswith("KBaseGenomeAnnotations.Assembly"):   # a one-bin unit: CheckM ran on a bare Assembly
        items = [{"ref": f"{head[6]}/{head[0]}/{head[4]}"}]
    else:
        items = rpc("ws", "Workspace.get_objects2", [{"objects": [{"ref": a.set}]}])[0]["data"][0]["data"]["items"]
    infos = rpc("ws", "Workspace.get_object_info3", [{"objects": [{"ref": it["ref"]} for it in items],
                                                      "includeMetadata": 1}])[0]["infos"]
    by_name = {i[1]: i for i in infos}
    names = {r["Bin Name"] for r in rows}
    if names != set(by_name) or set(ext) != names or len(rows) != len(items):
        sys.exit(f"CheckM rows {len(rows)} / bin_stats_ext {len(ext)} / set items {len(items)} do not match 1:1: "
                 f"{sorted(names ^ set(by_name))[:5]}")
    table, sel = [], []
    for r in sorted(rows, key=lambda r: r["Bin Name"]):
        e, i = ext[r["Bin Name"]], by_name[r["Bin Name"]]
        comp, cont = float(e["Completeness"]), float(e["Contamination"])
        q = comp >= a.min_completeness and cont < a.max_contamination
        row = {"assembly_name": i[1], "assembly_ref": f"{i[6]}/{i[0]}/{i[4]}", "md5": i[10].get("MD5"),
               "size_bp": i[10].get("Size"), "n_contigs": i[10].get("N Contigs"),
               "completeness": comp, "contamination": cont, "marker_lineage": e.get("marker lineage"),
               "genome_size_checkm": e.get("Genome size"), "n50": e.get("N50 (contigs)"),
               "summary_completeness": r["Completeness"], "summary_contamination": r["Contamination"],
               "pass_50_10": int(q)}
        table.append(row)
        if q:
            sel.append(row)
    cols = list(table[0])
    for path, rr in ((f"checkm_{a.tag}.tsv", table), (f"selected_{a.tag}.tsv", sel)):
        with open(os.path.join(a.out, path), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(rr)
    batch = []
    for r in sel:
        dst = os.path.join(a.out, "fasta", r["assembly_name"] + ".fna")
        src = None
        if a.local_fasta_dir:
            bid = r["assembly_name"].rsplit("__", 1)[0]
            cand = os.path.join(a.local_fasta_dir, bid[:-6] + ".fa" if bid.endswith(".fasta") else bid + ".fa")
            if os.path.isfile(cand) and asm_md5(cand) == r["md5"]:
                src = cand
        if src:
            if not os.path.exists(dst):
                os.link(src, dst)
            r["fasta_source"] = "local DAS Tool bin (AssemblyUtil md5 == KBase MD5)"
        else:
            fetch(r["assembly_ref"], dst, contig_check=False)
            r["fasta_source"] = "KBase blob (md5-verified)"
        batch.append(f"{os.path.abspath(dst)}\t{r['assembly_name']}\n")
    with open(os.path.join(a.out, f"gtdbtk_batchfile_{a.tag}.tsv"), "w") as fh:
        fh.writelines(batch)
    print(json.dumps({"bins": len(table), "selected_50_10": len(sel),
                      "sources": {s: sum(1 for r in sel if r["fasta_source"] == s) for s in {r["fasta_source"] for r in sel}},
                      "batchfile": os.path.join(a.out, f"gtdbtk_batchfile_{a.tag}.tsv"),
                      "example": {k: table[0][k] for k in ("assembly_name", "completeness", "contamination",
                                                           "summary_completeness", "pass_50_10")}}, indent=1))


if __name__ == "__main__":
    main()
