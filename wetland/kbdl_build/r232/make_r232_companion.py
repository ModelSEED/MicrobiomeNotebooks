#!/usr/bin/env python3
"""Build the KBDL SkaniMetadata companion for the GTDB R232 skani database (GTDB-Tk R232 package, consolidated layout).

Same format as the r220 companion (951d52…) and admin._parse_gtdb_lineage: key = companion_key(skani Ref_file) =
'<accession>_genomic.fna.gz'; value = {genus (no g__ prefix), species ('s__…'), lineage, ncbi_taxid: null}.
Usage: python -I make_r232_companion.py <release232>/taxonomy/gtdb_taxonomy.tsv out.json
(2026-11-15 build: 199,923 entries, sha256 fd647052c321a51c…)"""
import hashlib, json, sys
entries = {}
for line in open(sys.argv[1]):
    acc, lin = line.rstrip("\n").split("\t")
    acc = acc[3:] if acc[:3] in ("RS_", "GB_") else acc
    d = dict(x.split("__", 1) for x in lin.split(";") if "__" in x)
    entries[f"{acc}_genomic.fna.gz"] = {"genus": d.get("g") or None, "lineage": lin, "ncbi_taxid": None,
                                        "species": ("s__" + d["s"]) if d.get("s") else None}
payload = json.dumps(entries, indent=2, sort_keys=True).encode()
open(sys.argv[2], "wb").write(payload)
print(len(entries), "entries; sha256", hashlib.sha256(payload).hexdigest())
