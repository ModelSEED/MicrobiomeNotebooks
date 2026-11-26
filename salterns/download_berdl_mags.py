#!/usr/bin/env python
"""Download recently-updated genomes from a KBase narrative into BERDL_MAGs/.

Targets narrative https://narrative.kbase.us/narrative/155805 (workspace 155805,
cliff_bueno's Salt Pond narrative) and pulls every KBaseGenomes.Genome whose
latest version was saved within the last --days days (default 14).

Uses KBUtilLib (../../KBUtilLib) for all KBase access: KBGenomeUtils gives the
Workspace client, ref-path object retrieval, and Blobstore/Shock downloads.

For each genome three artifacts are written:
    genomes/<name>.json     full Genome object -- lossless, keeps RAST
                            functions, ontology terms, taxonomy, per-feature DNA
    proteins/<name>.faa     protein translations
    contigs/<name>.fna      assembly contigs, fetched from the Assembly object
                            the genome points at (usually a different workspace,
                            reached through a ref-path so the ACL resolves)

Auth: KBUtilLib reads ~/.kbase/token (or $KB_AUTH_TOKEN). Tokens expire; a
"login failed"/permission error usually just means you need a fresh one from
https://narrative.kbase.us -> Account -> Developer Tokens.

Usage:
    python download_berdl_mags.py --dry-run     # inventory only
    python download_berdl_mags.py               # ws 155805, last 14 days
    python download_berdl_mags.py --days 30 --limit 5

Re-running is safe: artifacts already on disk are kept and the genome is
skipped, so an interrupted run resumes where it left off.
"""

import argparse
import json
import os
import re
import sys
import time
from datetime import datetime, timedelta, timezone

# KBUtilLib is used from its checkout rather than site-packages (it is not
# installed into the venv); keep the repo-relative path working from anywhere.
_KBUTILLIB_SRC = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "KBUtilLib", "src")
)
if os.path.isdir(_KBUTILLIB_SRC) and _KBUTILLIB_SRC not in sys.path:
    sys.path.insert(0, _KBUTILLIB_SRC)

from kbutillib.kb_genome_utils import KBGenomeUtils  # noqa: E402

GENOME_TYPE_PREFIX = "KBaseGenomes.Genome"
GZIP_MAGIC = b"\x1f\x8b"
FASTA_WRAP = 60


def sanitize(name: str) -> str:
    """Turn a KBase object name into a safe-ish filename stem."""
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._-")
    return stem or "genome"


def human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f}B" if unit == "B" else f"{n:.1f}{unit}"
        n /= 1024.0


def parse_ws_date(ts: str) -> datetime:
    """Workspace save_date -> aware datetime ('2026-07-21T00:12:10+0000')."""
    for fmt in ("%Y-%m-%dT%H:%M:%S%z", "%Y-%m-%dT%H:%M:%S.%f%z"):
        try:
            return datetime.strptime(ts, fmt)
        except ValueError:
            continue
    raise ValueError(f"unrecognized workspace timestamp: {ts!r}")


def write_fasta(path: str, records, wrap: int = FASTA_WRAP) -> int:
    """Write (header, sequence) pairs; atomic via .part. Returns bytes written."""
    tmp = path + ".part"
    written = 0
    with open(tmp, "w") as fh:
        for header, seq in records:
            written += fh.write(f">{header}\n")
            for i in range(0, len(seq), wrap):
                written += fh.write(seq[i : i + wrap] + "\n")
    os.replace(tmp, path)
    return written


def genome_proteins(data: dict):
    """Yield (header, protein) for every translated feature in a Genome object.

    Genomes carry translations on `features`, on `cdss`, or on both (duplicated).
    Prefer `features`; fall back to `cdss` only when features yielded nothing.
    """
    def harvest(entries):
        out = []
        for ftr in entries or []:
            seq = ftr.get("protein_translation")
            if not seq:
                continue
            funcs = ftr.get("functions") or ([ftr["function"]] if ftr.get("function") else [])
            desc = "; ".join(f for f in funcs if f)
            out.append((f"{ftr.get('id', 'unknown')} {desc}".rstrip(), seq))
        return out

    return harvest(data.get("features")) or harvest(data.get("cdss"))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--ws", type=int, default=155805, help="workspace id (default 155805)")
    ap.add_argument("--days", type=int, default=14,
                    help="only genomes saved within this many days (default 14)")
    ap.add_argument("--out", default=os.path.join(here, "BERDL_MAGs"),
                    help="output dir (default ./BERDL_MAGs)")
    ap.add_argument("--limit", type=int, default=None, help="process only the first N (testing)")
    ap.add_argument("--dry-run", action="store_true", help="inventory only; no downloads")
    ap.add_argument("--skip-contigs", action="store_true", help="do not fetch assembly FASTA")
    ap.add_argument("--skip-proteins", action="store_true", help="do not write protein FASTA")
    args = ap.parse_args()

    cutoff = datetime.now(timezone.utc) - timedelta(days=args.days)

    try:
        u = KBGenomeUtils(kb_version="prod")
    except Exception as e:
        sys.exit(f"[auth] could not build KBUtilLib client: {e}")
    if not u.get_token(namespace="kbase"):
        sys.exit("[auth] no KBase token found -- put one in ~/.kbase/token or set $KB_AUTH_TOKEN.")

    try:
        u.set_ws(args.ws)
        wsinfo = u.ws_client().get_workspace_info({"id": args.ws})
    except Exception as e:
        sys.exit(f"[auth] workspace {args.ws} not readable: {e}\n"
                 f"       Most likely an expired token -- refresh ~/.kbase/token from "
                 f"https://narrative.kbase.us (Account -> Developer Tokens).")
    print(f"[ws] {args.ws}: name={wsinfo[1]!r} owner={wsinfo[2]} objects={wsinfo[4]}")

    all_objs = u.list_ws_objects(args.ws, include_metadata=False)
    genomes = [o for o in all_objs.values() if o[2].startswith(GENOME_TYPE_PREFIX)]
    recent = sorted((o for o in genomes if parse_ws_date(o[3]) >= cutoff), key=lambda o: o[0])
    print(f"[ws] {len(all_objs)} objects; {len(genomes)} genomes; "
          f"{len(recent)} saved since {cutoff:%Y-%m-%d %H:%M} UTC (last {args.days}d)")
    for o in recent:
        print(f"      {o[3]}  objid={o[0]:<6} v{o[4]:<3} {o[1][:60]:<60} {human(o[9])}")

    if args.dry_run:
        print("[dry-run] stopping before download.")
        return
    if not recent:
        print("[done] nothing to download.")
        return
    if args.limit:
        recent = recent[: args.limit]
        print(f"[limit] processing first {len(recent)}")

    dirs = {k: os.path.join(args.out, k) for k in ("genomes", "proteins", "contigs")}
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)

    manifest = []
    ok = skipped = failed = 0
    t0 = time.time()
    for i, o in enumerate(recent, 1):
        objid, name, otype, save_date, ver, size = o[0], o[1], o[2], o[3], o[4], o[9]
        ref = f"{args.ws}/{objid}/{ver}"
        stem = sanitize(name)
        tag = f"[{i}/{len(recent)}]"
        gpath = os.path.join(dirs["genomes"], stem + ".json")
        ppath = os.path.join(dirs["proteins"], stem + ".faa")
        entry = {"objid": objid, "name": name, "type": otype, "version": ver,
                 "ref": ref, "save_date": save_date, "ws_object_bytes": size}
        try:
            # Everything downstream needs the object, but on a resumed run the
            # JSON on disk is the same bytes -- reuse it instead of re-pulling
            # ~11 MB over the wire.
            have_genome = os.path.exists(gpath) and os.path.getsize(gpath) > 0
            need_prot = not args.skip_proteins and not os.path.exists(ppath)
            cpath_glob = [p for p in (os.path.join(dirs["contigs"], stem + ext)
                                      for ext in (".fna", ".fna.gz"))
                          if os.path.exists(p)]
            need_contigs = not args.skip_contigs and not cpath_glob

            # Already-complete genomes still get their metadata read back off the
            # local JSON, so manifest/summary are identical whether this was a
            # fresh run or a resume.
            complete = have_genome and not need_prot and not need_contigs
            if complete:
                print(f"{tag} have {stem[:60]}")

            if have_genome:
                with open(gpath) as fh:
                    data = json.load(fh)
            else:
                print(f"{tag} get {stem[:60]} ({human(size)})")
                obj = u.get_object(ref)
                if obj is None:
                    raise RuntimeError("get_object returned None")
                data = obj["data"]
                tmp = gpath + ".part"
                with open(tmp, "w") as fh:
                    json.dump(data, fh)
                os.replace(tmp, gpath)

            entry.update(
                genome_json=os.path.basename(gpath),
                genome_json_bytes=os.path.getsize(gpath),
                scientific_name=data.get("scientific_name"),
                taxonomy=data.get("taxonomy"),
                domain=data.get("domain"),
                dna_size=data.get("dna_size"),
                num_contigs=data.get("num_contigs"),
                gc_content=data.get("gc_content"),
                n_features=len(data.get("features") or []),
                assembly_ref=data.get("assembly_ref"),
            )

            if not args.skip_proteins:
                if os.path.exists(ppath):
                    entry["proteins"] = os.path.basename(ppath)
                else:
                    prots = genome_proteins(data)
                    if prots:
                        write_fasta(ppath, prots)
                        entry.update(proteins=os.path.basename(ppath), n_proteins=len(prots))
                    else:
                        entry["proteins_warning"] = "no protein_translation found"
                if "n_proteins" not in entry and entry.get("proteins"):
                    entry["n_proteins"] = sum(
                        1 for line in open(ppath) if line.startswith(">"))

            if not args.skip_contigs:
                if cpath_glob:
                    entry.update(contigs=os.path.basename(cpath_glob[0]),
                                 contigs_bytes=os.path.getsize(cpath_glob[0]))
                else:
                    aref = data.get("assembly_ref") or data.get("contigset_ref")
                    if not aref:
                        entry["contigs_warning"] = "genome has no assembly_ref"
                    else:
                        # Ref-path (genome;assembly): the assembly usually lives
                        # in another workspace and is only readable *through*
                        # the genome that references it.
                        adata = u.get_object(f"{ref};{aref}")["data"]
                        handle = adata.get("fasta_handle_ref") or adata.get("handle_ref")
                        if not handle:
                            entry["contigs_warning"] = "assembly has no fasta handle"
                        else:
                            dest = os.path.join(dirs["contigs"], stem + ".fna")
                            got = u.download_blob_file(handle, dest)
                            if not got:
                                raise RuntimeError(f"blob download failed for handle {handle}")
                            with open(dest, "rb") as fh:
                                magic = fh.read(2)
                            if magic == GZIP_MAGIC:
                                os.replace(dest, dest + ".gz")
                                dest += ".gz"
                            entry.update(contigs=os.path.basename(dest),
                                         contigs_bytes=os.path.getsize(dest),
                                         assembly_dna_size=adata.get("dna_size"),
                                         assembly_num_contigs=adata.get("num_contigs"))

            entry["status"] = "skipped" if complete else "downloaded"
            manifest.append(entry)
            if complete:
                skipped += 1
            else:
                ok += 1
                print(f"{tag} ok   {stem[:55]:<55} "
                      f"{entry.get('n_proteins', 0):>5} prot  "
                      f"{human(entry.get('contigs_bytes', 0)):>8} contigs  "
                      f"{(entry.get('taxonomy') or '')[:45]}")

        except Exception as e:
            print(f"{tag} FAILED {name}: {type(e).__name__}: {e}", file=sys.stderr)
            entry.update(status="failed", error=f"{type(e).__name__}: {e}")
            manifest.append(entry)
            failed += 1

    os.makedirs(args.out, exist_ok=True)
    man_path = os.path.join(args.out, "manifest.json")
    with open(man_path, "w") as fh:
        json.dump({"workspace": args.ws, "workspace_name": wsinfo[1],
                   "narrative": f"https://narrative.kbase.us/narrative/{args.ws}",
                   "window_days": args.days, "cutoff_utc": cutoff.isoformat(),
                   "generated_utc": datetime.now(timezone.utc).isoformat(),
                   "n_selected": len(recent), "downloaded": ok, "skipped": skipped,
                   "failed": failed, "genomes": manifest}, fh, indent=2)

    tsv_path = os.path.join(args.out, "genomes_summary.tsv")
    cols = ["name", "objid", "ref", "save_date", "taxonomy", "domain", "dna_size",
            "num_contigs", "gc_content", "n_features", "n_proteins", "assembly_ref", "status"]
    with open(tsv_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for e in manifest:
            fh.write("\t".join(str(e.get(c, "")) for c in cols) + "\n")

    print(f"\n[done] downloaded={ok} skipped={skipped} failed={failed} "
          f"in {time.time() - t0:.0f}s -> {args.out}")
    print(f"[manifest] {man_path}\n[summary]  {tsv_path}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
