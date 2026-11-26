#!/usr/bin/env python
"""KBaseMetagenomes.BinnedContigs -> DAS Tool contig2bin TSV (read-only on KBase).

contig2bin = '<contig_id>\t<bin label>' per binned contig, no header (DAS_Tool.R fread(header=F), col 1 trimmed after
the first whitespace). Bin labels get the binner prefix because KBase's MetaBAT2 and CONCOCT both name bins
'bin.NNN.fasta' (DAS Tool would otherwise rename them '<binner>_<bid>' itself).
Checks: object's assembly_ref == the expected Assembly UPA; every contig id is a FASTA id of that assembly and its
'len' equals the FASTA length; no contig in two bins of one binner.
Usage: python -I binnedcontigs_to_contig2bin.py <assembly.fna> <assembly_upa> <label>=<binnedcontigs_ref> ... --out DIR
"""
import argparse, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kbrpc import rpc


def fasta_lengths(path):
    lens, cur, n = {}, None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    lens[cur] = n
                cur, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if cur is not None:
        lens[cur] = n
    return lens


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta")
    ap.add_argument("assembly_upa")
    ap.add_argument("binners", nargs="+", help="label=ref, e.g. metabat2=275645/80/1")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    lens = fasta_lengths(a.fasta)
    report = {"fasta_contigs": len(lens), "binners": {}}
    for spec in a.binners:
        label, ref = spec.split("=", 1)
        o = rpc("ws", "Workspace.get_objects2", [{"objects": [{"ref": ref, "included": ["/assembly_ref", "/bins/[*]/bid",
                                                                                        "/bins/[*]/contigs"]}]}])[0]["data"][0]
        d = o["data"]
        if d["assembly_ref"] != a.assembly_upa:
            sys.exit(f"{label} {ref}: assembly_ref {d['assembly_ref']} != {a.assembly_upa}")
        seen, rows, badlen = {}, [], []
        for b in d.get("bins", []):
            lab = f"{label}.{b['bid'][:-6] if b['bid'].endswith('.fasta') else b['bid']}"
            for cid, ci in b["contigs"].items():
                if cid not in lens:
                    sys.exit(f"{label} {ref}: contig {cid} of {b['bid']} is not in {a.fasta}")
                if cid in seen:
                    sys.exit(f"{label} {ref}: contig {cid} in {seen[cid]} and {b['bid']}")
                if ci.get("len") not in (None, lens[cid]):
                    badlen.append(cid)
                seen[cid] = b["bid"]
                rows.append(f"{cid}\t{lab}\n")
        if badlen:
            sys.exit(f"{label} {ref}: {len(badlen)} contigs with len != FASTA length, e.g. {badlen[:3]}")
        path = os.path.join(a.out, f"{label}.contigs2bin.tsv")
        with open(path, "w") as fh:
            fh.writelines(rows)
        report["binners"][label] = {"ref": f"{o['info'][6]}/{o['info'][0]}/{o['info'][4]}", "name": o["info"][1],
                                    "bins": len(d.get("bins", [])), "contigs": len(rows),
                                    "bp": sum(lens[c] for c in seen), "tsv": path}
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
