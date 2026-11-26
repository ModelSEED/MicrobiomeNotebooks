#!/usr/bin/env python
"""Offline QC of the RAST fetch (no KBase calls): completeness, duplicates, per-set stats, code sanity.

fetch_rast.py keys its manifest by bin, so a bin present twice (two sets or twice in one set) would silently
overwrite; this script cross-checks the manifest against rast_inputs.json (the submitted assembly refs per set),
the rename log, calls/summary.tsv, the per-bin TSVs and the fetch log lines.
Usage: ~/Documents/py_venv/bin/python check_fetch.py [--log fetch_rast_full.log]
"""
import argparse, collections, csv, glob, json, os, re, statistics as st, sys
csv.field_size_limit(sys.maxsize)
HERE = '/home/freiburger/Documents/EmilyKin/community_modeling/analysis/r13/annotations/rast'
EK = '/home/freiburger/Documents/EmilyKin'
ap = argparse.ArgumentParser(); ap.add_argument('--log', default=f'{HERE}/fetch_rast_full.log'); a = ap.parse_args()
G = f'{HERE}/genomes'
man = list(csv.DictReader(open(f'{G}/_manifest.tsv'), delimiter='\t'))
ren = json.load(open(f'{EK}/kbase_mag_rename_log.json'))
summ = {r['bin']: r for r in csv.DictReader(open(f'{HERE}/../calls/summary.tsv'), delimiter='\t')}
inp = json.load(open(f'{EK}/community_modeling/plan_r13/kbase_rast/rast_inputs.json'))
jobs = json.load(open(f'{EK}/community_modeling/plan_r13/kbase_rast/rast_jobs.json'))
want = {'emilykin_MAGs_RAST_code25': inp['code25']}
for i, refs in enumerate(inp['code11_batches']):
    want[f'emilykin_MAGs_RAST_code11_b{i+1}'] = refs
P = []
def dup(xs): return [k for k, v in collections.Counter(xs).items() if v > 1]
print(f'manifest rows {len(man)}')
for col in ('bin', 'genome_ref', 'assembly_ref', 'iterativeID'):
    d = dup([r[col] for r in man])
    print(f'  distinct {col}: {len({r[col] for r in man})}  duplicates: {d[:5]}')
    if d: P.append(f'duplicate {col}: {d}')
bins = {r['bin'] for r in man}
for nm, ref in (('rename log', set(ren)), ('calls/summary.tsv', set(summ))):
    miss, extra = sorted(ref - bins), sorted(bins - ref)
    print(f'  vs {nm}: missing {len(miss)} {miss[:10]}  extra {len(extra)} {extra[:10]}')
    if miss or extra: P.append(f'{nm}: missing {miss} extra {extra}')
# assembly_ref -> bin consistent with rename log
for r in man:
    if ren.get(r['bin'], {}).get('upa') != r['assembly_ref']:
        P.append(f'{r["bin"]}: manifest assembly_ref {r["assembly_ref"]} != rename upa {ren.get(r["bin"],{}).get("upa")}')
    if not r['genome_name'].endswith('_' + r['iterativeID']) or not r['genome_name'].startswith(r['genomeset'] + '_'):
        P.append(f'{r["bin"]}: genome_name {r["genome_name"]} vs set {r["genomeset"]} / iid {r["iterativeID"]}')
# per set: submitted assembly refs == fetched
print('per set: submitted vs fetched (assembly refs)')
for s, refs in want.items():
    got = sorted(r['assembly_ref'] for r in man if r['genomeset'] == s)
    miss, extra = sorted(set(refs) - set(got)), sorted(set(got) - set(refs))
    ok = sorted(refs) == got
    print(f'  {s}: submitted {len(refs)} (jobs n={jobs[s]["n"]}) fetched {len(got)}  exact_match {ok}  missing {miss} extra {extra}')
    if not ok: P.append(f'{s}: set mismatch missing {miss} extra {extra}')
# TSV files
tsvs = sorted(os.path.basename(p)[:-4] for p in glob.glob(f'{G}/*.tsv') if not os.path.basename(p).startswith('_'))
tmps = glob.glob(f'{G}/.*.tmp')
print(f'TSV files {len(tsvs)}; leftover tmp {len(tmps)}; TSV bins not in manifest {sorted(set(tsvs)-bins)}; manifest bins without TSV {sorted(bins-set(tsvs))}')
if tmps or set(tsvs) != bins: P.append('tsv set != manifest set or tmp leftovers')
# per-genome stats from TSVs
per = {}
for r in man:
    b = r['bin']; n = 0; prot_len = []; stops_int = 0; stops_end = 0; types = collections.Counter(); multi = 0
    with open(f'{G}/{b}.tsv') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            n += 1; types[row['type']] += 1
            p = row['protein_translation']
            if p:
                prot_len.append(len(p.rstrip('*')))
                stops_int += p.rstrip('*').count('*'); stops_end += p.endswith('*')
            if row['n_loc'] != '1': multi += 1
    if n != int(r['n_features']): P.append(f'{b}: TSV rows {n} != manifest n_features {r["n_features"]}')
    per[b] = dict(n=n, mean_aa=st.mean(prot_len) if prot_len else 0, stops_int=stops_int, stops_end=stops_end,
                  multi=multi, types=types)
# per set summary
print('per set summary (n_features, n_cds median [min-max]; function share = func/cds pooled and per-genome median; hypothetical; no call)')
allr = []
for s in want:
    R = [r for r in man if r['genomeset'] == s]; allr += R
    if not R: print(f'  {s}: none'); continue
    nf = [int(r['n_features']) for r in R]; nc = [int(r['n_cds']) for r in R]
    fn = sum(int(r['n_cds_function']) for r in R); hy = sum(int(r['n_cds_hypothetical']) for r in R); no = sum(int(r['n_cds_nofunction']) for r in R)
    tot = sum(nc); share = [int(r['n_cds_function']) / int(r['n_cds']) for r in R if int(r['n_cds'])]
    codes = collections.Counter(r['genetic_code'] for r in R)
    ma = [per[r['bin']]['mean_aa'] for r in R]
    print(f'  {s}: genomes {len(R)} codes {dict(codes)} | features med {st.median(nf):.0f} [{min(nf)}-{max(nf)}] | CDS med {st.median(nc):.0f} [{min(nc)}-{max(nc)}] total {tot}'
          f' | function {fn/tot:.3f} (per-genome med {st.median(share):.3f}, {min(share):.3f}-{max(share):.3f}) hypo {hy/tot:.3f} nocall {no/tot:.3f}'
          f' | mean aa per genome {min(ma):.0f}-{max(ma):.0f} | internal stops {sum(per[r["bin"]]["stops_int"] for r in R)} | multi-loc {sum(per[r["bin"]]["multi"] for r in R)}')
R = allr; tot = sum(int(r['n_cds']) for r in R)
print(f'ALL: genomes {len(R)} features {sum(int(r["n_features"]) for r in R)} CDS {tot} function {sum(int(r["n_cds_function"]) for r in R)/tot:.3f}'
      f' hypo {sum(int(r["n_cds_hypothetical"]) for r in R)/tot:.3f} nocall {sum(int(r["n_cds_nofunction"]) for r in R)/tot:.3f}')
# by phylum
byp = collections.defaultdict(lambda: [0, 0, 0])
for r in R:
    ph = summ.get(r['bin'], {}).get('gtdb_phylum', '?'); byp[ph][0] += 1; byp[ph][1] += int(r['n_cds']); byp[ph][2] += int(r['n_cds_function'])
print('function share by phylum:', {k: (v[0], round(v[2] / v[1], 3)) for k, v in sorted(byp.items(), key=lambda x: -x[1][0])})
# CDS vs pyrodigal genes
ratio = sorted((int(r['n_cds']) / int(summ[r['bin']]['genes']), r['bin']) for r in R)
print(f'RAST CDS / pyrodigal genes: min {ratio[0][0]:.3f} ({ratio[0][1]}) median {st.median([x for x,_ in ratio]):.3f} max {ratio[-1][0]:.3f} ({ratio[-1][1]})')
share = sorted((int(r['n_cds_function']) / int(r['n_cds']), r['bin'], summ[r['bin']]['gtdb_phylum']) for r in R)
print('lowest function share:', [(b, ph, round(x, 3)) for x, b, ph in share[:6]], 'highest:', [(b, ph, round(x, 3)) for x, b, ph in share[-3:]])
ma = sorted((v['mean_aa'], b) for b, v in per.items())
print(f'mean protein aa per genome: min {ma[0][0]:.0f} ({ma[0][1]}) max {ma[-1][0]:.0f} ({ma[-1][1]}); genomes with internal stops: {[b for b,v in per.items() if v["stops_int"]]}')
ty = collections.Counter(); [ty.update(v['types']) for v in per.values()]
print('feature types:', dict(ty))
# fetch log: each fetched genome line printed once
if os.path.exists(a.log):
    L = [l.split('\t') for l in open(a.log) if l.startswith('  ') and '\t' in l]
    lb = [x[0].strip() for x in L]
    print(f'fetch log: {len(lb)} genome lines, distinct {len(set(lb))}, duplicates {dup(lb)}; PROBLEM lines:')
    for l in open(a.log):
        if l.startswith('PROBLEM') or 'Traceback' in l or 'Error' in l: print('   ', l.rstrip())
print(f'QC problems: {len(P)}')
for p in P: print('  QC', p)
sys.exit(1 if P else 0)
