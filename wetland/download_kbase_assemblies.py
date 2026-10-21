#!/usr/bin/env python
"""Download every assembly in a KBase narrative/workspace to a local folder.

Targets narrative https://narrative.kbase.us/narrative/261700 by default
(workspace id 261700). Uses cobrakbase's KBaseAPI for the Workspace + Handle
services and streams each assembly's FASTA straight from Shock/Blobstore.

Auth (first that resolves wins):
    --token <TOKEN>          explicit
    $KB_AUTH_TOKEN           environment
    ~/.kbase/token           file (cobrakbase default)

Get a fresh token from https://narrative.kbase.us -> Account -> "Developer
Tokens" (or copy the login token from a logged-in narrative session) and write
it to ~/.kbase/token, then re-run. Tokens expire, so "Invalid token" just means
you need a new one.

Usage:
    python download_kbase_assemblies.py                 # ws 261700 -> ./kbase_assemblies_261700
    python download_kbase_assemblies.py --dry-run       # inventory only, no downloads
    python download_kbase_assemblies.py --ws 261700 --out /some/dir

Re-running is safe: files already present with the expected byte size are
skipped, so an interrupted run resumes where it left off.
"""
import argparse
import json
import os
import re
import sys
import time

import requests

from cobrakbase import KBaseAPI

SHOCK_URL = "https://kbase.us/services/shock-api"
# Modern assemblies are KBaseGenomeAnnotations.Assembly; ContigSet is the legacy
# type (sequences sometimes stored inline rather than in Shock).
ASSEMBLY_TYPE_KEYS = ("Assembly", "ContigSet")

GZIP_MAGIC = b"\x1f\x8b"


def sanitize(name: str) -> str:
    """Turn a KBase object name into a safe-ish filename stem."""
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._-")
    return stem or "assembly"


def human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.1f}{unit}" if unit != "B" else f"{n}B"
        n /= 1024.0


def shock_headers(token: str) -> dict:
    return {"Authorization": "OAuth " + token}


def shock_node_meta(token: str, node_id: str) -> dict:
    r = requests.get(
        f"{SHOCK_URL}/node/{node_id}", headers=shock_headers(token), timeout=120
    )
    r.raise_for_status()
    return r.json()["data"]


def download_shock_node(token: str, node_id: str, dest: str) -> int:
    """Stream a Shock node's raw bytes to dest (atomic via .part). Returns bytes."""
    tmp = dest + ".part"
    written = 0
    with requests.get(
        f"{SHOCK_URL}/node/{node_id}?download_raw",
        headers=shock_headers(token),
        stream=True,
        timeout=1800,
    ) as r:
        r.raise_for_status()
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(1 << 20):  # 1 MiB
                if chunk:
                    written += fh.write(chunk)
    os.replace(tmp, dest)
    return written


def resolve_shock_node(api: KBaseAPI, data: dict):
    """Find the Shock node id backing an assembly object's FASTA.

    Handles the field variants seen across Assembly / ContigSet objects.
    Returns (node_id, handle_id_or_None) or (None, None) if not handle-backed.
    """
    # KBaseGenomeAnnotations.Assembly: handle service id in fasta_handle_ref
    href = data.get("fasta_handle_ref") or data.get("handle_ref")
    if href:
        handles = api.hs.hids_to_handles([href])
        return handles[0]["id"], href

    # Some objects embed the resolved handle directly
    fhi = data.get("fasta_handle_info") or {}
    handle = fhi.get("handle") or fhi
    if isinstance(handle, dict):
        if handle.get("id"):  # already a shock node id
            return handle["id"], handle.get("hid")
        if handle.get("hid"):
            handles = api.hs.hids_to_handles([handle["hid"]])
            return handles[0]["id"], handle["hid"]
    return None, None


def write_contigset_inline(data: dict, dest: str):
    """Fallback for legacy ContigSet objects that carry sequences inline.

    Returns bytes written, or None if the object has no usable inline sequence.
    """
    contigs = data.get("contigs") or []
    if not contigs or not any(c.get("sequence") for c in contigs):
        return None
    written = 0
    tmp = dest + ".part"
    with open(tmp, "w") as fh:
        for c in contigs:
            seq = c.get("sequence") or ""
            if not seq:
                continue
            cid = c.get("id") or c.get("name") or "contig"
            desc = c.get("description") or ""
            header = f">{cid} {desc}".rstrip()
            block = header + "\n" + "\n".join(
                seq[i : i + 70] for i in range(0, len(seq), 70)
            ) + "\n"
            written += fh.write(block)
    os.replace(tmp, dest)
    return written


def is_assembly_type(otype: str) -> bool:
    return any(k in otype for k in ASSEMBLY_TYPE_KEYS)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ws", type=int, default=261700, help="workspace id (default 261700)")
    default_out = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "kbase_assemblies_261700"
    )
    ap.add_argument("--out", default=default_out, help=f"output dir (default {default_out})")
    ap.add_argument("--token", default=None, help="KBase auth token (else env/file)")
    ap.add_argument("--limit", type=int, default=None,
                    help="only process the first N assemblies (for testing)")
    ap.add_argument("--dry-run", action="store_true", help="inventory only; no downloads")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)

    try:
        api = KBaseAPI(token=args.token)
    except Exception as e:
        sys.exit(f"[auth] could not build KBase client: {e}\n"
                 f"       Put a fresh token in ~/.kbase/token or pass --token.")
    token = api._token

    # Fail fast with a clear message if the token is stale.
    try:
        wsinfo = api.get_workspace_info(args.ws)
    except Exception as e:
        sys.exit(f"[auth] workspace {args.ws} not readable: {e}\n"
                 f"       Most likely an expired token — refresh ~/.kbase/token "
                 f"from https://narrative.kbase.us (Account -> Developer Tokens).")
    print(f"[ws] {args.ws}: name={wsinfo[1]!r} owner={wsinfo[2]} objects={wsinfo[4]}")

    all_objs = api.list_objects(args.ws)
    assemblies = [o for o in all_objs if is_assembly_type(o[2])]
    print(f"[ws] {len(all_objs)} objects total; {len(assemblies)} assembly-like")
    for o in assemblies:
        print(f"      objid={o[0]} v{o[4]:<3} {o[1][:55]:<55} {o[2]}")

    if args.dry_run:
        print("[dry-run] stopping before download.")
        return

    if args.limit:
        assemblies = assemblies[: args.limit]
        print(f"[limit] processing first {len(assemblies)} assemblies")

    manifest = []
    ok = skipped = failed = 0
    for i, o in enumerate(assemblies, 1):
        objid, name, otype, ver = o[0], o[1], o[2], o[4]
        ref = f"{args.ws}/{objid}/{ver}"
        stem = sanitize(name)
        try:
            # Pass the string ref (not the int objid) — cobrakbase's
            # process_workspace_identifiers calls .split(";") on it, which
            # crashes on an int. The ref also pins the exact version.
            data = api.get_object(ref)
            if data is None:
                raise RuntimeError("get_object returned None")

            node_id, handle_id = resolve_shock_node(api, data)
            entry = {
                "objid": objid, "name": name, "type": otype, "version": ver,
                "ref": ref, "shock_node": node_id, "handle": handle_id,
                "num_contigs": data.get("num_contigs"),
                "dna_size": data.get("dna_size"),
                "gc_content": data.get("gc_content"),
            }

            if node_id:
                meta = shock_node_meta(token, node_id)
                fmeta = meta.get("file") or {}
                exp_size = fmeta.get("size") or 0
                node_name = fmeta.get("name") or ""
                ext = ".fasta.gz" if node_name.endswith(".gz") else ".fasta"
                dest = os.path.join(args.out, stem + ext)

                if os.path.exists(dest) and exp_size and os.path.getsize(dest) == exp_size:
                    print(f"[{i}/{len(assemblies)}] skip (exists) {os.path.basename(dest)} "
                          f"({human(exp_size)})")
                    entry.update(file=os.path.basename(dest), bytes=exp_size, status="skipped")
                    manifest.append(entry); skipped += 1
                    continue

                print(f"[{i}/{len(assemblies)}] get {name} -> {os.path.basename(dest)} "
                      f"({human(exp_size) if exp_size else '?'})")
                got = download_shock_node(token, node_id, dest)
                # gzip sanity: rename if magic bytes disagree with chosen extension
                with open(dest, "rb") as fh:
                    magic = fh.read(2)
                if magic == GZIP_MAGIC and not dest.endswith(".gz"):
                    os.replace(dest, dest + ".gz"); dest += ".gz"
                if exp_size and got != exp_size:
                    raise RuntimeError(f"size mismatch: got {got} expected {exp_size}")
                entry.update(file=os.path.basename(dest), bytes=got, status="downloaded")
                manifest.append(entry); ok += 1

            else:
                # No Shock handle — try inline ContigSet sequences.
                dest = os.path.join(args.out, stem + ".fasta")
                got = write_contigset_inline(data, dest)
                if got is None:
                    raise RuntimeError(
                        f"no fasta handle and no inline sequences; top-level keys="
                        f"{sorted(data.keys())[:20]}")
                print(f"[{i}/{len(assemblies)}] rebuilt inline {os.path.basename(dest)} "
                      f"({human(got)})")
                entry.update(file=os.path.basename(dest), bytes=got, status="inline")
                manifest.append(entry); ok += 1

        except Exception as e:
            print(f"[{i}/{len(assemblies)}] FAILED {name}: {e}", file=sys.stderr)
            manifest.append({"objid": objid, "name": name, "type": otype,
                             "version": ver, "ref": ref, "status": "failed",
                             "error": str(e)})
            failed += 1

    man_path = os.path.join(args.out, "manifest.json")
    with open(man_path, "w") as fh:
        json.dump({"workspace": args.ws, "workspace_name": wsinfo[1],
                   "n_assemblies": len(assemblies), "downloaded": ok,
                   "skipped": skipped, "failed": failed, "objects": manifest},
                  fh, indent=2)
    print(f"\n[done] downloaded={ok} skipped={skipped} failed={failed} "
          f"-> {args.out}\n[manifest] {man_path}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
