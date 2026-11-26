#!/usr/bin/env python3
"""Step 7 (supplement, poplar, after 03_aggregate): the five panel KOs whose KOfam threshold is '-' (never called by
the KOfam/DRAM rule): K16967 dmoA, K23996 mpnH/phnY*, K05915 phnZ, K18049 mpnS, K12905 hepD.

Rule (read_markers v1 'mpnS/hepD rule', 04_build_db.py special_label, extended to the five KOs): a gene with NO
passing KO (over all 27,981 profiles) gets the no-threshold KO if that hit has full E <= 1e-30 (-Z 12,055,118) and
HMM coverage >= 0.5 (min hmm_from .. max hmm_to over the hit's domains / profile length, as kofam postprocess.py
cov_of). If several of the five qualify, the lowest E wins. `is_top_scoring` = the KO is also the gene's top-scoring
profile over all 27,981 (v1 additionally required this); both variants are aggregated locally (05_panel_tables.py).
Two-pass inputs: hits and coverage from pass 1 (pass1/hits.tsv.gz, pass1/hmm/*.domtblout; every gene was screened
with these five profiles), 'no passing KO' and the top-scoring profile from pass 2 (every hit with E <= 1e-30 to one
of the five made its gene a candidate, 08_select_candidates.py).
Writes nothr/rep_calls.tsv, results/nothreshold_panel_by_sample_long.tsv.
usage: 07_nothreshold_panel.py <W>
"""
import gzip
import os
import sys
import pandas as pd

W = sys.argv[1]
os.chdir(W)
KOS = ["K16967", "K23996", "K05915", "K18049", "K12905"]
os.makedirs("nothr", exist_ok=True)
H = pd.read_csv("pass1/hits.tsv.gz", sep="\t", dtype={"rep": str, "ko": str})
H = H[H.ko.isin(KOS)]
calls = pd.concat([pd.read_csv(f"pass2/results/chunks/{f}", sep="\t", dtype=str, keep_default_na=False,
                               usecols=["rep", "best_ko", "top_ko"])
                   for f in sorted(os.listdir("pass2/results/chunks")) if f.endswith(".calls.tsv.gz")]).set_index("rep")
cand = H[H.full_evalue <= 1e-30].join(calls, on="rep")
assert cand.best_ko.notna().all(), "a no-threshold candidate was not searched in pass 2"
n_e30 = len(cand)
cand = cand[cand.best_ko == ""].copy()
need = set(zip(cand.rep, cand.ko))
cov = {}
for f in sorted(os.listdir("pass1/hmm")):
    if not f.endswith(".domtblout"):
        continue
    for line in open(f"pass1/hmm/{f}"):
        if line[0] == "#":
            continue
        x = line.split(None, 22)
        key = (x[0], x[3])
        if key not in need:
            continue
        leng, hf, ht = int(x[5]), int(x[15]), int(x[16])
        a = cov.get(key)
        cov[key] = [hf, ht, leng] if a is None else [min(a[0], hf), max(a[1], ht), leng]
cand["hmm_cov"] = [((cov[k][1] - cov[k][0] + 1) / cov[k][2]) if k in cov else 0.0 for k in zip(cand.rep, cand.ko)]
q = cand[cand.hmm_cov >= 0.5].sort_values(["rep", "full_evalue", "full_score"], ascending=[True, True, False])
q = q.drop_duplicates("rep").copy()
q["is_top_scoring"] = q.ko == q.top_ko
q.to_csv("nothr/rep_calls.tsv", sep="\t", index=False)
qi = q.set_index("rep")
rows = []
for s in sorted(f[:-7] for f in os.listdir("genes_ko") if f.endswith(".tsv.gz")):
    g = pd.read_csv(f"genes_ko/{s}.tsv.gz", sep="\t", usecols=["gene", "rep", "depth", "status", "domain"],
                    keep_default_na=False)
    g = g[g.rep.isin(qi.index)].join(qi[["ko", "is_top_scoring"]], on="rep")
    for (ko, top, st, dom), d in g.groupby(["ko", "is_top_scoring", "status", "domain"]):
        rows.append(dict(sample=s, ko=ko, is_top_scoring=top, status=st, domain=dom, n_genes=len(d),
                         sum_depth=d.depth.sum()))
pd.DataFrame(rows, columns=["sample", "ko", "is_top_scoring", "status", "domain", "n_genes", "sum_depth"]).to_csv(
    "results/nothreshold_panel_by_sample_long.tsv", sep="\t", index=False, float_format="%.6g")
print("hits of the five KOs", len(H), "with E <= 1e-30", n_e30, "of which no passing KO", len(cand),
      "called (cov >= 0.5)", len(q), q.ko.value_counts().to_dict(), "top-scoring", int(q.is_top_scoring.sum()))
