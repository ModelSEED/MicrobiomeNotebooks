#!/usr/bin/env python
"""Per-sample read and metaSPAdes run statistics of the 2026 assemblies (read-only on KBase).

For each of the 39 per-sample units:
  - raw import and Trimmomatic output read objects (ledger stages 'import', 'trim'): read_count, total_bases,
    read_length_mean (Workspace.get_objects2 with included paths);
  - the kb_SPAdes job log (execution_engine2.get_job_logs; job id from the ledger), saved gzipped to spades_logs/,
    parsed for SPAdes version, k list, command line, mean read length after correction, BayesHammer correction
    totals, insert size, % of pairs aligned to long edges (last estimate, K127), whether scaffolds.fasta was saved,
    QUAST N50 and mimalloc allocation warnings.
Writes spades_run_stats.tsv.
Usage: ~/Documents/py_venv/bin/python -I spades_run_stats.py
"""
import csv, gzip, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "local"))
from kbrpc import rpc  # noqa: E402

LOGDIR = os.path.join(HERE, "spades_logs")


def reads_stats(ref):
    d = rpc("ws", "Workspace.get_objects2", [{"objects": [{"ref": ref, "included": [
        "/read_count", "/total_bases", "/read_length_mean", "/read_orientation_outward", "/insert_size_mean"]}]}])[0][
        "data"][0]["data"]
    return d


def job_log(job_id):
    path = os.path.join(LOGDIR, f"{job_id}.log.gz")
    if os.path.exists(path):
        return gzip.open(path, "rt").read().split("\n")
    r = rpc("ee2", "execution_engine2.get_job_logs", [{"job_id": job_id, "skip_lines": 0, "limit": 200000}])[0]
    lines = [l["line"] for l in r["lines"]]
    assert len(lines) == r["count"] or r.get("last_line_number", 0) + 1 == len(lines), (job_id, len(lines))
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(lines))
    return lines


def parse(lines):
    txt = "\n".join(lines)
    out = {}
    m = re.search(r"SPAdes version: (\S+)", txt)
    out["spades_version"] = m.group(1) if m else ""
    m = re.search(r"SPADES CMD:(\[.*?\])", txt)
    out["spades_cmd"] = m.group(1) if m else ""
    m = re.search(r"-k (\d[\d,]+)", out["spades_cmd"])
    out["k_list"] = m.group(1) if m else ""
    m = re.findall(r"Average read length ([\d.]+)", txt)
    out["avg_read_len_K33"] = m[0] if m else ""
    m = re.search(r"Total (\d+) reads processed", txt)
    out["reads_into_bayeshammer"] = m.group(1) if m else ""
    m = re.search(r"Correction done\. Changed (\d+) bases in (\d+) reads", txt)
    out["bh_bases_changed"], out["bh_reads_changed"] = (m.group(1), m.group(2)) if m else ("", "")
    ins = re.findall(r"Insert size = ([\d.]+), deviation = ([\d.]+)", txt)
    out["insert_size"], out["insert_sd"] = ins[-1] if ins else ("", "")
    al = re.findall(r"(\d+) paired reads \(([\d.]+)% of all\) aligned to long edges", txt)
    out["pct_pairs_on_long_edges_last"] = al[-1][1] if al else ""
    out["pct_pairs_on_long_edges_all"] = ";".join(a[1] for a in al)
    out["saved_file"] = "scaffolds.fasta" if "scaffolds.fasta.filtered" in txt else (
        "contigs.fasta" if "contigs.fasta.filtered" in txt else "")
    m = re.search(r"N50 = (\d+), L50 = (\d+), Total length = (\d+)", txt)
    out["quast_N50_ge1000"] = m.group(1) if m else ""
    out["mimalloc_warnings"] = txt.count("mimalloc: warning")
    m = re.search(r"Memory limit \(in Gb\): (\d+)", txt)
    out["mem_limit_gb"] = m.group(1) if m else ""
    m = re.findall(r"K(\d+) started", txt)
    out["k_run"] = ",".join(m)
    return out


def main():
    os.makedirs(LOGDIR, exist_ok=True)
    t = json.load(open(os.path.join(ROOT, "ledger.json")))["tasks"]
    rows = []
    for k in sorted(k for k in t if k.startswith("assemble:")):
        v = t[k]
        s = v["sample"]
        imp, trim = t[f"import:{s}"]["output_ref"], t[f"trim:{s}"]["output_ref"]
        ri, rt = reads_stats(imp), reads_stats(trim)
        att = [a for a in v["attempts"] if a.get("status") == "completed"][-1]
        p = parse(job_log(att["job_id"]))
        row = {"sample": s, "import_ref": imp, "trim_ref": trim, "assembly_job": att["job_id"],
               "raw_read_count": ri.get("read_count"), "raw_bases": ri.get("total_bases"),
               "trim_read_count": rt.get("read_count"), "trim_bases": rt.get("total_bases"),
               "trim_read_len_mean": rt.get("read_length_mean"),
               "pct_bases_kept_by_trim": round(100 * rt["total_bases"] / ri["total_bases"], 2),
               "pct_reads_kept_by_trim": round(100 * rt["read_count"] / ri["read_count"], 2), **p}
        rows.append(row)
        print(s, row["pct_bases_kept_by_trim"], p["spades_version"], p["k_list"], p["insert_size"],
              p["pct_pairs_on_long_edges_last"], p["saved_file"], p["quast_N50_ge1000"], p["mimalloc_warnings"],
              flush=True)
    with open(os.path.join(HERE, "spades_run_stats.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
