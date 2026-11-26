"""Merge per-chunk DRAM2 raw-annotations.tsv files into one table for all MAGs.

Run with the run tree's own DRAM2 distill env (pandas 1.5.2), e.g.
    conda_cache/dram2_beta17_distill/bin/python -I bin/merge_dram2.py [--runs chunk_01 ...]

* Adds a `mag` column (original MAG name, from inputs/name_map.tsv) right after `input_fasta`
  (the DRAM2-safe name, dots -> '-').
* Columns are the union over chunks, in DRAM2's order (base columns first, then per-database
  blocks); chunks are streamed one at a time so memory stays at one chunk.
* Checks: every chunk has a DONE marker; every MAG of the name map appears in exactly one chunk;
  no duplicated (input_fasta, query_id); rows per MAG == genes in its Prodigal .faa; every MAG has
  >= 1 KOfam hit; no "Error loading" in the chunk's RAW/combine_annotations.log.
* Adds kofam_ko_threshold / kofam_ko_score_type / kofam_passes_ko_threshold (from the install's
  databases/kofam/kofam_ko_list.tsv). DRAM2 beta17 keeps the best KOfam hit per gene at
  E <= 1e-5 and does NOT apply KOfam's per-KO thresholds; kofam_score_rank is that hit's
  full-sequence score, so pass = kofam_score_rank >= threshold. For score_type "domain" KOs this
  is an upper-bound approximation (the domain score is not kept). DRAM2's own columns are unchanged.
* Also links every per-MAG rRNA / tRNA table into merged/rrnas and merged/trnas, which the final
  all-MAG distill launch (--rrnas / --trnas) consumes.
"""
import argparse
import os
import sys
from pathlib import Path

import pandas as pd

B = Path("/scratch1/afreiburger/wetland_dram2")
KO_LIST = Path("/scratch1/fliu/hub_scratch/chenry/DRAM2/databases/kofam/kofam_ko_list.tsv")
KO_COLS = ["kofam_ko_threshold", "kofam_ko_score_type", "kofam_passes_ko_threshold"]
BASE = ["query_id", "input_fasta", "mag", "start_position", "stop_position",
        "strandedness", "rank", "gene_number"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="*", default=None,
                    help="run names under runs/ (default: every runs/chunk_*)")
    ap.add_argument("--out", default=str(B / "merged" / "raw-annotations.all.tsv"))
    ap.add_argument("--allow-partial", action="store_true",
                    help="do not require every name-map MAG to be present")
    a = ap.parse_args()

    runs = a.runs or sorted(p.name for p in (B / "runs").glob("chunk_*") if p.is_dir())
    nm = pd.read_csv(B / "inputs" / "name_map.tsv", sep="\t", dtype=str)
    to_mag = dict(zip(nm.dram2_name, nm.mag))
    ko = pd.read_csv(KO_LIST, sep="\t", usecols=["knum", "threshold", "score_type"], dtype=str)
    ko_thr = dict(zip(ko.knum, ko.threshold))
    ko_type = dict(zip(ko.knum, ko.score_type))

    paths = []
    for r in runs:
        if not (B / "runs" / r / "DONE").exists():
            sys.exit(f"runs/{r} has no DONE marker")
        p = B / "runs" / r / "out" / "RAW" / "raw-annotations.tsv"
        if not p.is_file():
            sys.exit(f"missing {p}")
        paths.append((r, p))

    # union of columns, DRAM2-style order
    cols = []
    for _, p in paths:
        with open(p) as fh:
            for c in fh.readline().rstrip("\n").split("\t"):
                if c not in cols:
                    cols.append(c)
    cols = [c for c in cols if c != "mag"]
    base = [c for c in BASE if c in cols or c == "mag"]
    rest = [c for c in cols if c not in base]
    order = base + rest + (KO_COLS if "kofam_id" in cols else [])

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".part")
    seen_mags: dict = {}
    total = 0
    with open(tmp, "w", newline="") as fo:
        fo.write("\t".join(order) + "\n")
        for r, p in paths:
            df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
            unknown = sorted(set(df.input_fasta) - set(to_mag))
            if unknown:
                sys.exit(f"{r}: input_fasta values not in name map, e.g. {unknown[:3]}")
            dup = df.duplicated(["input_fasta", "query_id"]).sum()
            if dup:
                sys.exit(f"{r}: {dup} duplicated (input_fasta, query_id) rows")
            for m in df.input_fasta.unique():
                if m in seen_mags:
                    sys.exit(f"{m} appears in {seen_mags[m]} and {r}")
                seen_mags[m] = r
            # rows per MAG == Prodigal genes; >= 1 KOfam hit per MAG; no read errors in combine log
            rows = df.groupby("input_fasta").size()
            kof = (df.kofam_id != "").groupby(df.input_fasta).sum() if "kofam_id" in df else rows * 0
            prod = B / "runs" / r / "out" / "Prodigal_v2.6.3"
            bad = []
            for m, n in rows.items():
                with open(prod / f"{m}_called_genes.faa") as fh:
                    g = sum(1 for line in fh if line.startswith(">"))
                if n != g or kof.get(m, 0) == 0:
                    bad.append(f"{m} rows={n} genes={g} kofam_hits={kof.get(m, 0)}")
            if bad:
                sys.exit(f"{r}: {len(bad)} MAGs fail rows/genes/kofam checks, e.g. {bad[:3]}")
            log = B / "runs" / r / "out" / "RAW" / "combine_annotations.log"
            if log.exists() and "Error loading" in log.read_text():
                sys.exit(f"{r}: read errors in {log}")
            df["mag"] = df.input_fasta.map(to_mag)
            if "kofam_id" in df.columns:
                thr = df.kofam_id.map(ko_thr).fillna("")
                df["kofam_ko_threshold"] = thr
                df["kofam_ko_score_type"] = df.kofam_id.map(ko_type).fillna("")
                t = pd.to_numeric(thr, errors="coerce")
                sc = pd.to_numeric(df.get("kofam_score_rank", ""), errors="coerce")
                passed = (sc >= t).map({True: "True", False: "False"})
                df["kofam_passes_ko_threshold"] = passed.where(t.notna() & sc.notna(), "")
            df = df.reindex(columns=order, fill_value="")
            df.to_csv(fo, sep="\t", index=False, header=False)
            total += len(df)
            print(f"{r}: {len(df):,} rows, {df.input_fasta.nunique()} MAGs", flush=True)

    missing = sorted(set(to_mag) - set(seen_mags))
    if missing and not a.allow_partial:
        os.remove(tmp)
        sys.exit(f"{len(missing)} name-map MAGs absent from the merged runs, e.g. {missing[:3]}")
    os.replace(tmp, out)
    print(f"merged {total:,} rows, {len(seen_mags)} MAGs, {len(order)} columns -> {out}")

    # rRNA / tRNA per-MAG tables for the all-MAG distill launch
    for kind in ("rrnas", "trnas"):
        d = B / "merged" / kind
        d.mkdir(parents=True, exist_ok=True)
        n = 0
        for r, _ in paths:
            for f in (B / "runs" / r / "out" / "tRNA_rRNA").glob(f"*_processed_{kind}.tsv"):
                link = d / f.name
                if link.is_symlink() or link.exists():
                    link.unlink()
                link.symlink_to(f)
                n += 1
        print(f"{kind}: {n} per-MAG tables linked into {d}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
