#!/usr/bin/env python
"""Stream-scan prFBA's single-line PGFam cache JSON ({"<gid>": ["PGF_...", ...], ...})
without loading it: regex the keys chunk by chunk, record each key's value span, and
derive the PGFam count from the span (each entry '"PGF_dddddddd", ' is 16 bytes).

Read-only on the cache; writes one small TSV (gid, n_pgfams, n_exact) + a JSON summary.
Usage: scan_pgfam_cache.py CACHE_JSON OUT_TSV OUT_SUMMARY_JSON
"""
import json
import os
import re
import sys
import time

CACHE, OUT_TSV, OUT_SUM = sys.argv[1:4]
KEY = re.compile(rb'"(\d+\.\d+)": \[')
CHUNK = 32 << 20
OVERLAP = 64

t0 = time.time()
size = os.path.getsize(CACHE)
keys = []            # (gid, abs_offset_of_value_start)  value start = index just after '['
seen_offsets = set()
with open(CACHE, "rb") as fh:
    pos = 0          # absolute offset of buf[0]
    buf = b""
    while True:
        data = fh.read(CHUNK)
        if not data:
            break
        buf = buf + data
        for m in KEY.finditer(buf):
            off = pos + m.end()
            if off in seen_offsets:
                continue
            seen_offsets.add(off)
            keys.append((m.group(1).decode(), off))
        keep = min(OVERLAP, len(buf))
        pos += len(buf) - keep
        buf = buf[-keep:]
        # bound the dedupe set to the overlap window
        seen_offsets = {o for o in seen_offsets if o >= pos}

keys.sort(key=lambda x: x[1])
# value span: from value start to the ']' before the next key (or end of file).
# Each value is a flat list of fixed-width strings, so the closing ']' sits right
# before '", "<next key>' -> ', "' is 3 bytes; the final value is followed by ']}'.
rows = []
n_empty = 0
n_inexact = 0
with open(CACHE, "rb") as fh:
    for i, (gid, off) in enumerate(keys):
        end = keys[i + 1][1] - len(f'"{keys[i + 1][0]}": [') - 2 if i + 1 < len(keys) else size - 1
        # end now points just past ']' for list i (', ' separator = 2 bytes)
        fh.seek(off)
        first = fh.read(1)
        if first == b"]":
            n, exact = 0, True
            n_empty += 1
        else:
            inside = (end - 1) - off            # bytes between '[' and ']'
            n_f = (inside + 2) / 16.0
            n = int(round(n_f))
            exact = abs(n_f - n) < 1e-9
            if not exact:
                n_inexact += 1
        rows.append((gid, n, int(exact)))

with open(OUT_TSV, "w") as out:
    out.write("bvbrc_genome_id\tn_pgfams\tcount_exact\n")
    for gid, n, ex in rows:
        out.write(f"{gid}\t{n}\t{ex}\n")
dups = len(rows) - len({r[0] for r in rows})
summary = dict(cache=CACHE, bytes=size, mtime=time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(os.path.getmtime(CACHE))),
               n_keys=len(rows), n_unique_keys=len({r[0] for r in rows}), n_duplicate_keys=dups,
               n_empty=n_empty, n_nonempty=len(rows) - n_empty, n_count_inexact=n_inexact,
               scan_seconds=round(time.time() - t0, 1))
json.dump(summary, open(OUT_SUM, "w"), indent=1)
print(json.dumps(summary, indent=1))
