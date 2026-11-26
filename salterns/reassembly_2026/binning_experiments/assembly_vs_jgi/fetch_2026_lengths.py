#!/usr/bin/env python
"""Contig lengths of the 43 KBase 2026 assemblies (39 metaSPAdes per-sample + 4 MEGAHIT sediment co-assemblies).

Read-only on KBase: Workspace.get_objects2 with included paths (/contigs/*/length, /base_counts, ...), via
../../local/kbrpc.py. Assembly refs come from ../../ledger.json (stage 'assemble' or 'coassemble', output_ref).
Writes:
  lengths_2026.tsv.gz      unit, contig_id, length   (one row per contig; contigs < 1000 bp were never saved:
                           kb_SPAdes/MEGAHIT ran with min_contig_length 1000)
  assemblies_2026.tsv      unit, ref, assembler, num_contigs, dna_size, n_bases_N, node_rank_check
Usage: ~/Documents/py_venv/bin/python -I fetch_2026_lengths.py
"""
import csv, gzip, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "local"))
from kbrpc import rpc  # noqa: E402  (read-only whitelist)

NODE = re.compile(r"NODE_(\d+)_length_(\d+)_cov_")


def main():
    ledger = json.load(open(os.path.join(ROOT, "ledger.json")))["tasks"]
    units = []
    for k, v in ledger.items():
        if v["stage"] in ("assemble", "coassemble") and v.get("output_ref"):
            app = v["attempts"][-1]["app_id"]
            units.append((v["sample"], v["output_ref"], "metaSPAdes 3.15.3" if "SPAdes" in app else "MEGAHIT 1.2.9"))
    units.sort()
    assert len(units) == 43, len(units)
    out_len = gzip.open(os.path.join(HERE, "lengths_2026.tsv.gz"), "wt")
    out_len.write("unit\tcontig_id\tlength\n")
    rows = []
    for unit, ref, asm in units:
        o = rpc("ws", "Workspace.get_objects2", [{"objects": [{"ref": ref, "included": [
            "/contigs/*/length", "/num_contigs", "/dna_size", "/base_counts", "/assembly_id"]}]}])[0]["data"][0]
        d = o["data"]
        contigs = d["contigs"]
        assert len(contigs) == d["num_contigs"], unit
        assert sum(c["length"] for c in contigs.values()) == d["dna_size"], unit
        # SPAdes numbers NODE_n by descending length: check that rank n really is the length rank
        check = "n/a (MEGAHIT ids)"
        if asm.startswith("metaSPAdes"):
            pts = sorted((int(m.group(1)), int(m.group(2))) for m in (NODE.match(c) for c in contigs) if m)
            viol = sum(1 for a, b in zip(pts, pts[1:]) if b[1] > a[1])
            lenok = all(int(NODE.match(c).group(2)) == v["length"] for c, v in contigs.items())
            check = (f"{len(pts)} NODE ids; max n {pts[-1][0]}; rank-order violations {viol}; "
                     f"name length == length: {lenok}")
        for cid, c in sorted(contigs.items(), key=lambda kv: -kv[1]["length"]):
            out_len.write(f"{unit}\t{cid}\t{c['length']}\n")
        rows.append({"unit": unit, "ref": ref, "assembler": asm, "num_contigs": d["num_contigs"],
                     "dna_size": d["dna_size"], "n_bases_N": d.get("base_counts", {}).get("N", 0),
                     "node_rank_check": check})
        print(unit, ref, d["num_contigs"], d["dna_size"], check, flush=True)
    out_len.close()
    with open(os.path.join(HERE, "assemblies_2026.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
