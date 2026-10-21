#!/usr/bin/env python
"""Compute assembly QC stats for every downloaded KBase assembly FASTA and
cross-check them against the values KBase recorded in manifest.json.

Reads each FASTA directly (authoritative for the downloaded files), computes
length/contig/N50/GC metrics in parallel, and writes assembly_qc.csv.

Metagenome assemblies vary enormously, so this script only *measures*; the
sensible/complete-vs-redo judgement is applied downstream from the CSV.
"""
import csv
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
ASM_DIR = os.path.join(HERE, "kbase_assemblies_261700")
MANIFEST = os.path.join(ASM_DIR, "manifest.json")
OUT_CSV = os.path.join(HERE, "assembly_qc.csv")

# contig-length thresholds we bin on (bp)
THRESHOLDS = (1_000, 5_000, 10_000, 50_000, 100_000)


def scan_fasta(path):
    """Single-pass stats for one FASTA. Handles wrapped or single-line seqs."""
    lengths = []          # per-contig length
    cur = 0               # current contig length
    gc = 0                # G/C/g/c bytes
    at = 0                # A/T/a/t bytes
    started = False
    with open(path, "rb") as fh:
        for line in fh:
            if line[:1] == b">":
                if started:
                    lengths.append(cur)
                cur = 0
                started = True
            else:
                seq = line.rstrip(b"\r\n")
                cur += len(seq)
                gc += seq.count(b"G") + seq.count(b"C") + seq.count(b"g") + seq.count(b"c")
                at += seq.count(b"A") + seq.count(b"T") + seq.count(b"a") + seq.count(b"t")
    if started:
        lengths.append(cur)

    n = len(lengths)
    total = sum(lengths)
    stats = {
        "file": os.path.basename(path),
        "n_contigs": n,
        "total_bp": total,
        "gc_frac": round(gc / (gc + at), 5) if (gc + at) else 0.0,
        "longest": max(lengths) if lengths else 0,
        "shortest": min(lengths) if lengths else 0,
        "mean_len": round(total / n, 1) if n else 0.0,
    }
    # N50 / L50
    n50 = l50 = 0
    if total:
        half = total / 2.0
        acc = 0
        for i, L in enumerate(sorted(lengths, reverse=True), 1):
            acc += L
            if acc >= half:
                n50, l50 = L, i
                break
    stats["n50"] = n50
    stats["l50"] = l50
    # threshold bins
    for t in THRESHOLDS:
        cnt = sum(1 for L in lengths if L >= t)
        bp = sum(L for L in lengths if L >= t)
        stats[f"n_ge_{t}"] = cnt
        stats[f"bp_ge_{t}"] = bp
    stats["frac_bp_ge_1000"] = round(stats["bp_ge_1000"] / total, 4) if total else 0.0
    return stats


def main():
    manifest = json.load(open(MANIFEST))
    kb = {o["file"]: o for o in manifest["objects"] if o.get("file")}
    files = sorted(os.path.join(ASM_DIR, f) for f in kb)
    print(f"[qc] scanning {len(files)} FASTAs with {min(32, os.cpu_count())} workers", flush=True)

    rows = []
    with ProcessPoolExecutor(max_workers=min(32, os.cpu_count())) as ex:
        futs = {ex.submit(scan_fasta, f): f for f in files}
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            rows.append(r)
            if i % 25 == 0 or i == len(files):
                print(f"[qc] {i}/{len(files)} done", flush=True)

    # merge KBase-recorded values + integrity deltas
    for r in rows:
        o = kb[r["file"]]
        r["sample"] = r["file"].split("_reads.MEGAHIT")[0]
        r["kb_dna_size"] = o.get("dna_size")
        r["kb_num_contigs"] = o.get("num_contigs")
        r["kb_gc"] = o.get("gc_content")
        r["objid"] = o.get("objid")
        r["ref"] = o.get("ref")
        r["bytes"] = o.get("bytes")
        # relative disagreement vs KBase (download/integrity check)
        kd, kc = o.get("dna_size"), o.get("num_contigs")
        r["dna_size_reldiff"] = round(abs(r["total_bp"] - kd) / kd, 6) if kd else None
        r["num_contigs_reldiff"] = round(abs(r["n_contigs"] - kc) / kc, 6) if kc else None

    rows.sort(key=lambda x: x["total_bp"])
    cols = ["sample", "file", "objid", "ref", "n_contigs", "total_bp", "n50", "l50",
            "longest", "shortest", "mean_len", "gc_frac",
            "n_ge_1000", "bp_ge_1000", "frac_bp_ge_1000",
            "n_ge_5000", "n_ge_10000", "n_ge_50000", "n_ge_100000",
            "bp_ge_10000", "bp_ge_50000", "bp_ge_100000",
            "kb_dna_size", "kb_num_contigs", "kb_gc",
            "dna_size_reldiff", "num_contigs_reldiff", "bytes"]
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c) for c in cols})
    print(f"[qc] wrote {OUT_CSV}", flush=True)

    # quick integrity summary
    bad_int = [r for r in rows
               if (r["dna_size_reldiff"] or 0) > 0.001 or (r["num_contigs_reldiff"] or 0) > 0.001]
    print(f"[qc] integrity mismatches vs KBase (>0.1%): {len(bad_int)}")
    for r in bad_int[:20]:
        print(f"     {r['sample']}: total_bp={r['total_bp']} kb={r['kb_dna_size']} "
              f"(dna d={r['dna_size_reldiff']}, ctg d={r['num_contigs_reldiff']})")


if __name__ == "__main__":
    main()
