#!/usr/bin/env python
"""Compare the critic scan (critic/scan/hits.tsv) with the data part (files_plan.tsv) and the
docs part (occurrences.tsv). Writes critic/scan/missed_files.tsv and critic/scan/variant_hits.tsv."""
import csv, os, sys
from collections import defaultdict

PLAN = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
S = os.path.join(PLAN, 'critic', 'scan')
csv.field_size_limit(sys.maxsize)
hits = list(csv.DictReader(open(os.path.join(S, 'hits.tsv')), delimiter='\t', quoting=csv.QUOTE_NONE))
data_files = {r['file']: r for r in csv.DictReader(open(os.path.join(PLAN, 'data', 'files_plan.tsv')), delimiter='\t')}
docs_files = {r['file'] for r in csv.DictReader(open(os.path.join(PLAN, 'docs', 'occurrences.tsv')), delimiter='\t')}
binfiles = {r['file'] for r in csv.DictReader(open(os.path.join(PLAN, 'data', 'binary_for_docs.tsv')), delimiter='\t')} \
    if os.path.exists(os.path.join(PLAN, 'data', 'binary_for_docs.tsv')) else set()

per_file = defaultdict(lambda: defaultdict(int))
ctx = {}
for h in hits:
    if h['id_kind'] != 'old':
        continue
    per_file[(h['status'], h['file'])][h['variant']] += int(h['count'])
    ctx.setdefault((h['status'], h['file']), (h['id'], h['context'][:220]))

with open(os.path.join(S, 'missed_files.tsv'), 'w') as f:
    f.write('status\tfile\tin_data_plan\tin_docs_occ\tin_data_binary_list\tvariants\texample_id\texample_context\n')
    for (st, fn), v in sorted(per_file.items()):
        ind, ino, inb = fn in data_files, fn in docs_files, fn in binfiles
        if ind or ino:
            continue
        f.write(f"{st}\t{fn}\t{ind}\t{ino}\t{inb}\t{dict(v)}\t{ctx[(st, fn)][0]}\t{ctx[(st, fn)][1]}\n")

with open(os.path.join(S, 'variant_hits.tsv'), 'w') as f:
    f.write('status\tfile\tid\tvariant\tcount\tdata_action\tcontext\n')
    for h in hits:
        if h['variant'] in ('strict',):
            continue
        act = data_files.get(h['file'], {}).get('action', '(not in data plan)')
        f.write(f"{h['status']}\t{h['file']}\t{h['id']}\t{h['variant']}\t{h['count']}\t{act}\t{h['context'][:260]}\n")

n_tr = sum(1 for (st, fn) in per_file if st == 'tracked')
n_un = sum(1 for (st, fn) in per_file if st == 'untracked')
miss_tr = sum(1 for (st, fn) in per_file if st == 'tracked' and fn not in data_files and fn not in docs_files)
print('files with hits: tracked', n_tr, 'untracked', n_un, '| tracked not in data/docs plans', miss_tr)
