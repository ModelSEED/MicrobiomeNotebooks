#!/usr/bin/env python
"""Fetch the contig FASTA of every salterns MAG into MAG_contigs/.

The MAG universe is the 450 DASTool bins listed in jplfaria/Cliff MAGs/:
    310_MAGS.txt              310 dereplicated   (Assemblies in ws 77580)
    140_non_dereped_MAGS.txt  140 non-dereplicated (Assemblies in ws 163264)
The local RAST genome JSONs (jplfaria/Cliff MAGs/RASTified_MAG_json/) carry
per-gene DNA and proteins but not the contigs; each names its Assembly through
`assembly_ref`, which is what this script downloads.

MAGs whose contigs BERDL_MAGs/contigs/ already holds (the `.dlt` downloads made
by download_berdl_mags.py, which point at the same Assembly objects) are
hard-linked instead of re-downloaded, so the folder is complete at no extra disk.

Every FASTA is checked against its Assembly's num_contigs and dna_size before it
is kept; MAG_contigs/index.tsv records provenance and an md5 per file.

Usage:
    python download_mag_contigs.py --dry-run    # resolve refs, download nothing
    python download_mag_contigs.py              # fetch whatever is missing

Re-running is safe: verified files already on disk are kept.
"""

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import sys
import time

_KBUTILLIB_SRC = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "KBUtilLib", "src")
)
if os.path.isdir(_KBUTILLIB_SRC) and _KBUTILLIB_SRC not in sys.path:
    sys.path.insert(0, _KBUTILLIB_SRC)

from kbutillib.kb_genome_utils import KBGenomeUtils  # noqa: E402

GZIP_MAGIC = b"\x1f\x8b"
ASSEMBLY_REF = re.compile(r'"assembly_ref"\s*:\s*"([^"]+)"')
MAG_KEY = re.compile(r"Salt_Pond_MetaG_?(.+?)_MG_DASTool_bins[._](concoct_out|metabat|maxbin)\.(\d+)")
FIELDS = ["/fasta_handle_ref", "/dna_size", "/num_contigs"]


def mag_key(name: str):
    """sample|binner|bin# -- collapses the .contigs__.RAST / _.dlt / bins. vs bins_ variants."""
    m = MAG_KEY.search(name)
    return f"{m[1]}|{m[2]}|{int(m[3])}" if m else None


def fasta_stats(path: str):
    """(n_records, total_sequence_length) of a plain or gzipped FASTA."""
    with open(path, "rb") as fh:
        gz = fh.read(2) == GZIP_MAGIC
    n = total = 0
    with (gzip.open(path, "rt") if gz else open(path)) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
            else:
                total += len(line.strip())
    return n, total


def retry(fn, *a, tries: int = 5, **kw):
    """Call fn, backing off 2/4/8/16 s on errors -- KBase's edge returns transient 5xx (520/523)."""
    for attempt in range(tries):
        try:
            return fn(*a, **kw)
        except Exception as e:
            if attempt == tries - 1:
                raise
            wait = 2 ** (attempt + 1)
            print(f"    retry in {wait}s after {type(e).__name__}: {str(e)[:120]}", file=sys.stderr)
            time.sleep(wait)


def md5sum(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    jpl = os.path.join(here, "jplfaria", "Cliff MAGs")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(here, "MAG_contigs"), help="output dir (default ./MAG_contigs)")
    ap.add_argument("--berdl", default=os.path.join(here, "BERDL_MAGs"),
                    help="BERDL_MAGs dir whose contigs/ are hard-linked when they match")
    ap.add_argument("--dry-run", action="store_true", help="resolve refs and report; download nothing")
    args = ap.parse_args()

    # ---- the 450 MAGs and the Assembly each RAST genome was built from -------------
    mags = []
    for list_file, population in (("310_MAGS.txt", "derep310"), ("140_non_dereped_MAGS.txt", "nonderep140")):
        with open(os.path.join(jpl, list_file)) as fh:
            mags += [(line.strip(), population) for line in fh if line.strip()]
    rast_dir = os.path.join(jpl, "RASTified_MAG_json")
    aref_by_key = {}
    for fn in os.listdir(rast_dir):
        if fn.endswith(".json"):
            with open(os.path.join(rast_dir, fn)) as fh:
                m = ASSEMBLY_REF.search(fh.read(4000))  # first key of the genome object
            if m:
                aref_by_key[mag_key(fn)] = m[1]
    missing = [name for name, _ in mags if mag_key(name) not in aref_by_key]
    if missing:
        sys.exit(f"[refs] no assembly_ref for {len(missing)} MAGs, e.g. {missing[:3]}")

    # ---- BERDL contigs already on disk, by assembly_ref ------------------------------
    berdl = {}
    man_path = os.path.join(args.berdl, "manifest.json")
    if os.path.exists(man_path):
        with open(man_path) as fh:
            for g in json.load(fh)["genomes"]:
                p = os.path.join(args.berdl, "contigs", g.get("contigs", ""))
                if g.get("contigs") and os.path.exists(p):
                    berdl[g["assembly_ref"]] = p

    try:
        u = KBGenomeUtils(kb_version="prod")
        ws = u.ws_client()
    except Exception as e:
        sys.exit(f"[auth] could not build KBUtilLib client: {e}")

    # ---- Assembly metadata (handle, size, contig count) in batches ----------------------
    arefs = [aref_by_key[mag_key(name)] for name, _ in mags]
    meta = {}
    for i in range(0, len(arefs), 25):
        chunk = arefs[i : i + 25]
        res = retry(ws.get_objects2, {"objects": [{"ref": r, "included": FIELDS} for r in chunk]})["data"]
        for r, d in zip(chunk, res):
            meta[r] = {"name": d["info"][1], **{f[1:]: d["data"].get(f[1:]) for f in FIELDS}}
    n_link = sum(r in berdl for r in arefs)
    print(f"[refs] {len(mags)} MAGs -> {len(set(arefs))} Assemblies; "
          f"{n_link} reusable from {args.berdl}/contigs, {len(mags) - n_link} to download")
    if args.dry_run:
        print("[dry-run] stopping before download.")
        return

    os.makedirs(args.out, exist_ok=True)
    rows, linked, fetched, kept, failed = [], 0, 0, 0, 0
    t0 = time.time()
    for i, ((name, population), aref) in enumerate(zip(mags, arefs), 1):
        a = meta[aref]
        dest = os.path.join(args.out, a["name"] + ".fna")
        tag = f"[{i}/{len(mags)}]"
        row = {"mag": re.sub(r"\.RAST$", "", name), "population": population, "assembly_ref": aref,
               "assembly_name": a["name"], "fasta_handle_ref": a["fasta_handle_ref"],
               "num_contigs": a["num_contigs"], "dna_size": a["dna_size"]}
        try:
            have = [p for p in (dest, dest + ".gz") if os.path.exists(p)]
            if have:
                dest, source = have[0], "existing"
            elif aref in berdl:
                src = berdl[aref]
                dest += ".gz" if src.endswith(".gz") else ""
                try:
                    os.link(src, dest)
                    source = "hardlink:BERDL_MAGs/contigs"
                except OSError:  # different filesystem -> plain copy
                    import shutil
                    shutil.copy2(src, dest)
                    source = "copy:BERDL_MAGs/contigs"
            else:
                tmp = dest + ".part"
                if not retry(u.download_blob_file, a["fasta_handle_ref"], tmp):
                    raise RuntimeError(f"blob download failed for {a['fasta_handle_ref']}")
                with open(tmp, "rb") as fh:
                    if fh.read(2) == GZIP_MAGIC:
                        dest += ".gz"
                os.replace(tmp, dest)
                source = "kbase"
            n, total = fasta_stats(dest)
            if (n, total) != (a["num_contigs"], a["dna_size"]):
                if source == "kbase":
                    os.remove(dest)
                raise RuntimeError(f"FASTA has {n} contigs / {total} bp; Assembly says "
                                   f"{a['num_contigs']} / {a['dna_size']}")
            row.update(file=os.path.basename(dest), bytes=os.path.getsize(dest), md5=md5sum(dest),
                       source=source, status="ok")
            if source == "kbase":
                fetched += 1
                print(f"{tag} got  {a['name'][:70]:<70} {n:>5} contigs {total:>9} bp")
            elif source == "existing":
                kept += 1
            else:
                linked += 1
        except Exception as e:
            print(f"{tag} FAILED {name}: {type(e).__name__}: {e}", file=sys.stderr)
            row.update(status=f"failed: {type(e).__name__}: {e}")
            failed += 1
        rows.append(row)

    cols = ["mag", "population", "assembly_ref", "assembly_name", "fasta_handle_ref", "num_contigs",
            "dna_size", "file", "bytes", "md5", "source", "status"]
    with open(os.path.join(args.out, "index.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\n[done] downloaded={fetched} hardlinked={linked} already_present={kept} failed={failed} "
          f"in {time.time() - t0:.0f}s -> {args.out}  (index.tsv written)")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
