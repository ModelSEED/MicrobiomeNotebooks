#!/usr/bin/env python
"""Download a KBaseGenomeAnnotations.Assembly's FASTA from this machine and verify it (read-only on KBase).

Path: Workspace get_objects2 (/fasta_handle_ref, /md5, /dna_size, /num_contigs) -> HandleService hids_to_handles
-> Blobstore node via the PUBLIC endpoint https://kbase.us/services/shock-api/node/<id>?download_raw
(the handle's own url, https://services.kbase.us/services/shock-api = prod-dock01, is unreachable from here).
Checks: md5 of the downloaded bytes == handle remote_md5 == node checksum; gzip is unpacked if present;
FASTA ids (first word) / lengths / per-contig md5 / AssemblyUtil's assembly md5 == the Assembly object.
Usage: python -I fetch_assembly.py <ref or ref;path> <out.fna> [--no-contig-check]
"""
import gzip, hashlib, http.client, json, os, shutil, socket, sys, time, urllib.error, urllib.request
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kbrpc import rpc, token

BLOB_PUBLIC = "https://kbase.us/services/shock-api"
NET_ERRORS = (urllib.error.URLError, TimeoutError, socket.timeout, ConnectionError, http.client.HTTPException)
WAITS = (10, 30, 60, 120, 300, 600)                  # ~19 min of retries for a Blobstore outage


def _retry(fn, what):
    """fn() again after a transport error or a 5xx; a 4xx is raised at once."""
    for wait in WAITS + (None,):
        try:
            return fn()
        except NET_ERRORS as e:
            if (isinstance(e, urllib.error.HTTPError) and e.code < 500) or wait is None:
                raise
            print(f"{what}: {type(e).__name__}: {str(e)[:200]}; retrying in {wait} s", file=sys.stderr, flush=True)
            time.sleep(wait)


def fetch(ref, out, contig_check=True):
    t0 = time.time()
    inc = ["/fasta_handle_ref", "/md5", "/dna_size", "/num_contigs", "/assembly_id", "/type"]
    if contig_check:
        inc += ["/contigs/*/contig_id", "/contigs/*/length", "/contigs/*/md5"]
    o = rpc("ws", "Workspace.get_objects2", [{"objects": [{"ref": ref, "included": inc}] if ";" not in ref else
                                              [{"ref": ref.split(";")[-1], "to_obj_ref_path": ref.split(";")[:-1],
                                                "included": inc}]}])[0]["data"][0]
    d, info = o["data"], o["info"]
    hid = d["fasta_handle_ref"]
    h = rpc("hs", "AbstractHandle.hids_to_handles", [[hid]])[0][0]
    node = h["id"]
    meta = _retry(lambda: json.load(urllib.request.urlopen(urllib.request.Request(
        f"{BLOB_PUBLIC}/node/{node}", headers={"Authorization": "OAuth " + token()}), timeout=120))["data"],
        "node metadata")
    tmp = out + ".part"

    def download():                                  # restarts the .part file from scratch on every try
        md5, n = hashlib.md5(), 0
        req = urllib.request.Request(f"{BLOB_PUBLIC}/node/{node}?download_raw",
                                     headers={"Authorization": "OAuth " + token()})
        with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as fh:
            want = r.headers.get("Content-Length") or meta.get("file", {}).get("size")
            while True:
                b = r.read(1 << 22)
                if not b:
                    break
                n += len(b)
                md5.update(b)
                fh.write(b)
        if want is not None and n != int(want):     # read(amt) returns b"" on an early close instead of raising
            raise http.client.IncompleteRead(b"", int(want) - n)          # -> _retry downloads again
        return md5.hexdigest()

    try:
        got = _retry(download, "download")
    except Exception:
        if os.path.exists(tmp):                      # do not leave a partial (possibly GB-sized) file behind
            os.remove(tmp)
        raise
    rep = {"ref": f"{info[6]}/{info[0]}/{info[4]}", "name": info[1], "type": info[2], "hid": hid, "node": node,
           "handle_url": h.get("url"), "handle_file_name": h.get("file_name"), "remote_md5": h.get("remote_md5"),
           "node_file": meta.get("file", {}).get("name"), "node_size": meta.get("file", {}).get("size"),
           "node_md5": (meta.get("file", {}).get("checksum") or {}).get("md5"), "download_md5": got,
           "bytes": os.path.getsize(tmp), "seconds": round(time.time() - t0, 1)}
    assert got == rep["node_md5"], ("node checksum mismatch", rep)
    if rep["remote_md5"]:
        assert got == rep["remote_md5"], ("handle remote_md5 mismatch", rep)
    with open(tmp, "rb") as fh:
        rep["gzip"] = fh.read(2) == b"\x1f\x8b"
    if rep["gzip"]:
        with gzip.open(tmp, "rb") as src, open(tmp + ".u", "wb") as dst:
            shutil.copyfileobj(src, dst, 1 << 22)
        os.replace(tmp + ".u", tmp)
    # parse: first word of the header = AssemblyUtil's contig_id (Bio.SeqIO record.id)
    seqs, order, desc = {}, [], 0
    cur, buf = None, []
    with open(tmp) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf).upper()
                hdr = line[1:].rstrip("\n")
                cur = hdr.split()[0] if hdr.split() else ""
                desc += hdr != cur
                assert cur not in seqs, ("duplicate id", cur)
                order.append(cur)
                buf = []
            else:
                buf.append(line.strip())
        if cur is not None:
            seqs[cur] = "".join(buf).upper()
    cm = {c: hashlib.md5(s.encode()).hexdigest() for c, s in seqs.items()}
    rep.update(n_fasta=len(seqs), bp_fasta=sum(map(len, seqs.values())), headers_with_description=desc,
               assembly_md5_recomputed=hashlib.md5(",".join(sorted(cm.values())).encode()).hexdigest(),
               obj_md5=d.get("md5"), obj_num_contigs=d.get("num_contigs"), obj_dna_size=d.get("dna_size"),
               first_ids=order[:3])
    assert rep["n_fasta"] == d["num_contigs"] and rep["bp_fasta"] == d["dna_size"], rep
    assert rep["assembly_md5_recomputed"] == d["md5"], rep
    if contig_check:
        oc = d["contigs"]
        assert set(oc) == set(seqs), "contig id sets differ"
        bad = [c for c in seqs if oc[c]["length"] != len(seqs[c]) or oc[c]["md5"] != cm[c] or oc[c]["contig_id"] != c]
        assert not bad, ("contig mismatches", bad[:5])
        rep["contig_check"] = f"{len(seqs)} ids, lengths and md5s identical to the Assembly object"
    if rep["headers_with_description"]:          # DAS Tool's prodigal_parallel/pullseq want bare ids
        with open(tmp) as src, open(tmp + ".ids", "w") as dst:
            for line in src:
                dst.write(">" + line[1:].split()[0] + "\n" if line.startswith(">") else line)
        os.replace(tmp + ".ids", tmp)
    os.replace(tmp, out)
    rep["seconds"] = round(time.time() - t0, 1)
    return rep


if __name__ == "__main__":
    print(json.dumps(fetch(sys.argv[1], sys.argv[2], "--no-contig-check" not in sys.argv), indent=1))
