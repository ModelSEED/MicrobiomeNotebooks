#!/usr/bin/env python
"""Join reads metadata to assembly QC, classify each assembly, and diagnose the
failed ones (data-limited vs processing failure) using read depth + read QC."""
import re, json
import pandas as pd, numpy as np
from cobrakbase import KBaseAPI

WS = 261700
srr_re = re.compile(r'(SRR\d+)')
api = KBaseAPI()
objs = api.list_objects(WS, include_metadata=True)

# full reads metadata per SRR (keep richest = max total_bases)
reads = {}
for o in objs:
    if 'PairedEndLibrary' not in o[2]:
        continue
    m = srr_re.search(o[1])
    if not m:
        continue
    md = o[10] or {}
    def f(k):
        v = md.get(k)
        try: return float(v)
        except (TypeError, ValueError): return None
    rec = {'reads_name': o[1], 'total_bases': f('total_bases'), 'read_count': f('read_count'),
           'read_len': f('read_length_mean'), 'qual_mean': f('qual_mean'),
           'n_dup': f('number_of_duplicates'), 'reads_gc': f('gc_content'),
           'phred': md.get('phred_type')}
    srr = m.group(1)
    if srr not in reads or (rec['total_bases'] or 0) > (reads[srr]['total_bases'] or 0):
        reads[srr] = rec

df = pd.read_csv('assembly_qc.csv')
# make re-runnable: strip any derived/merged columns from a prior run
_drop = ['input_bases','input_reads','read_len','qual_mean','n_dup','reads_gc',
         'asm_per_input','dup_frac','tier','diagnosis','total_bases','read_count',
         'asm_Mb','input_Gb','reads_M','dup%','ratio']
df = df.loc[:, ~df.columns.duplicated()]
df = df.drop(columns=[c for c in _drop if c in df.columns])
for k in ['total_bases','read_count','read_len','qual_mean','n_dup','reads_gc']:
    df[k] = df['sample'].map(lambda s: reads.get(s, {}).get(k))
df.rename(columns={'total_bases':'input_bases','read_count':'input_reads'}, inplace=True)
df['asm_per_input'] = df['total_bp']/df['input_bases']
df['dup_frac'] = df['n_dup']/df['input_reads']

# ---- classification ----
def classify(r):
    bp = r['total_bp']
    ratio = r['asm_per_input']
    if bp < 2_000_000:
        return 'FAILED'
    if bp < 20_000_000 or (pd.notna(ratio) and ratio < 0.003):
        return 'POOR'
    return 'OK'
df['tier'] = df.apply(classify, axis=1)

def reason(r):
    if r['tier'] == 'OK':
        return ''
    ig = (r['input_bases'] or 0)/1e9
    rr = r['asm_per_input']
    bits = []
    if pd.notna(rr):
        bits.append(f"assembled {r['total_bp']/1e6:.2f}Mb from {ig:.2f}Gb input "
                    f"(ratio {rr:.4f} vs cohort median 0.012)")
    if r['input_bases'] and r['input_bases'] < 3e8:
        bits.append("LOW INPUT depth -> possibly data-limited")
    else:
        bits.append("ample input -> under-assembly = process/QC failure")
    if pd.notna(r['dup_frac']) and r['dup_frac'] > 0.5:
        bits.append(f"HIGH duplication {r['dup_frac']*100:.0f}%")
    if pd.notna(r['qual_mean']) and r['qual_mean'] < 25:
        bits.append(f"low read qual {r['qual_mean']:.0f}")
    return "; ".join(bits)
df['diagnosis'] = df.apply(reason, axis=1)

df.sort_values('total_bp', inplace=True)
df.to_csv('assembly_qc.csv', index=False)
flagged = df[df.tier != 'OK'].copy()
flagged.to_csv('assembly_flagged.csv', index=False)

print("=== tier counts ===")
print(df.tier.value_counts().reindex(['OK','POOR','FAILED']).to_string())
print(f"\n=== OK cohort reference (n={ (df.tier=='OK').sum() }) ===")
ok = df[df.tier=='OK']
print(f"  total_bp  median {ok.total_bp.median()/1e6:.0f} Mb  (5-95%: {ok.total_bp.quantile(.05)/1e6:.0f}-{ok.total_bp.quantile(.95)/1e6:.0f} Mb)")
print(f"  n50       median {ok.n50.median():.0f}")
print(f"  asm/input median {ok.asm_per_input.median():.4f}")
print(f"  dup_frac  median {ok.dup_frac.median():.2f}  qual_mean median {ok.qual_mean.median():.1f}")

print(f"\n=== FAILED + POOR (n={len(flagged)}) — sorted by size ===")
cols = ['sample','tier','total_bp','n50','input_bases','input_reads','dup_frac','qual_mean','asm_per_input']
show = flagged[cols].copy()
show['asm_Mb'] = (show.total_bp/1e6).round(3)
show['input_Gb'] = (show.input_bases/1e9).round(2)
show['reads_M'] = (show.input_reads/1e6).round(1)
show['dup%'] = (show.dup_frac*100).round(0)
show['ratio'] = show.asm_per_input.round(4)
print(show[['sample','tier','asm_Mb','n50','input_Gb','reads_M','dup%','qual_mean','ratio']].to_string(index=False))
