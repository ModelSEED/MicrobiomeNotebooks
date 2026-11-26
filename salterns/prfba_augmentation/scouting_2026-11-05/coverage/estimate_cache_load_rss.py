#!/usr/bin/env python
"""Estimate the RSS cost of prFBA's whole-file json.load of the PGFam cache
(prefetch_gene_families.load_cache L90-93; select_references.bvbrc_gene_provider L705,
which also converts every list to a set) by parsing a contiguous slice of the first
N keys and extrapolating per PGFam entry. Read-only on the cache.
Usage: estimate_cache_load_rss.py CACHE_JSON INDEX_TSV N_KEYS
"""
import csv
import gc
import json
import sys

CACHE, INDEX, NK = sys.argv[1], sys.argv[2], int(sys.argv[3])


def rss():
    with open("/proc/self/statm") as fh:
        return int(fh.read().split()[1]) * 4096


rows = list(csv.DictReader(open(INDEX), delimiter="\t"))
total_entries = sum(int(r["n_pgfams"]) for r in rows)
sel = rows[:NK]
n_entries = sum(int(r["n_pgfams"]) for r in sel)
# the first NK keys are contiguous from byte 1; slice up to the ']' closing key NK-1
nbytes = 1 + sum(len(f'"{r["bvbrc_genome_id"]}": [') + max(16 * int(r["n_pgfams"]) - 2, 0) + 1 + 2 for r in sel) - 2
with open(CACHE, "rb") as fh:
    raw = fh.read(nbytes)
text = raw.decode() + "}"
del raw
gc.collect()
r0 = rss()
d = json.loads(text)
gc.collect()
r1 = rss()
assert len(d) == NK and sum(len(v) for v in d.values()) == n_entries, "slice mis-cut"
s = {g: set(v) for g, v in d.items()}
gc.collect()
r2 = rss()
del d
gc.collect()
r3 = rss()
per_list = (r1 - r0) / n_entries
per_both = (r2 - r0) / n_entries
per_set = (r3 - r0) / n_entries
out = dict(sample_keys=NK, sample_entries=n_entries, sample_json_bytes=len(text),
           bytes_per_entry_lists=round(per_list, 1), bytes_per_entry_lists_plus_sets=round(per_both, 1),
           bytes_per_entry_sets_only=round(per_set, 1), total_entries=total_entries,
           est_GB_json_load=round(per_list * total_entries / 1e9, 1),
           est_GB_peak_load_then_set_convert=round(per_both * total_entries / 1e9, 1),
           est_GB_sets_resident=round(per_set * total_entries / 1e9, 1))
print(json.dumps(out, indent=1))
