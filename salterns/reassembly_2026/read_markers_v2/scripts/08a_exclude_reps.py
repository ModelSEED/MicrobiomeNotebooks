#!/usr/bin/env python
"""Step 8a: marker representatives longer than their family's 99th length percentile are excluded at classification
time (hits to them are ignored). In the validation, background reads from non-ribosomal genes hit such over-long
'ribosomal protein' references (fusions / mis-predicted genes, e.g. a 263-aa 'uS11', a 344-aa KEGG 'uS9') at 50-70 %
identity. Writes refs/excluded_reps.tsv."""
import os
import pandas as pd
os.chdir("/scratch1/afreiburger/salterns_2026/read_markers_v2")
ln, cur = {}, None
for l in open("db/survey.faa"):
    if l.startswith(">"):
        cur = l[1:].strip()
    else:
        ln[cur] = ln.get(cur, 0) + len(l.strip())
c = pd.read_csv("refs/clu/clusters99.tsv", sep="\t", header=None, names=["label", "rep", "member"]).drop_duplicates("rep")
c["len"] = c.rep.map(ln)
MARKERS = ["uL2", "uL11", "uS3", "uS11", "uS12", "uS9", "uL4", "uL16"]
out = []
for lab, g in c[c.label.isin(MARKERS)].groupby("label"):
    q = g.len.quantile(0.99)
    x = g[g.len > q]
    out.append(x.assign(p99=q))
    print(lab, "p99", q, "excluded", len(x))
pd.concat(out)[["rep", "label", "len", "p99"]].to_csv("refs/excluded_reps.tsv", sep="\t", index=False)

# MpnS / HEPD references from the UniProt name queries ('methylphosphonate synthase' / gene mpnS; 'hydroxyethylphosphonate
# dioxygenase' / gene hepD) are kept only if they meet the same rule as the salterns proteins: top-scoring panel KO =
# K18049 or K12905, full E <= 1e-30, HMM coverage >= 0.5. 'gene:hepD' also matched heptaprenyl diphosphate synthase
# component genes, and two 'MpnS' entries have no K18049 hit at all.
import sys
sys.argv = ["x"]
exec(open("scripts/04_build_db_v2.py").read().split("if __name__")[0])
h = parse_hmm(["hmm_new/ext.tbl"], ["hmm_new/ext.domtbl"])
m = pd.read_csv("refs/labeled_meta.tsv", sep="\t", dtype=str).fillna("")
bad = []
for lab in ("mpnS", "hepD"):
    for pid in m.pid[(m.label == lab) & (m.source != "salterns")]:
        hits = h.get(pid, {})
        top = max(hits, key=lambda k: hits[k][1]) if hits else None
        ok = top in ("K18049", "K12905") and hits[top][0] <= 1e-30 and hits[top][3] >= 0.5
        if not ok:
            bad.append(pid)
reps = set(c.rep)
bad_reps = [p for p in bad if p in reps]
ex = pd.read_csv("refs/excluded_reps.tsv", sep="\t")
lab_of = dict(zip(c.rep, c.label))
ex = pd.concat([ex, pd.DataFrame({"rep": bad_reps, "label": [lab_of[p] for p in bad_reps], "len": [ln[p] for p in bad_reps],
                                  "p99": "mpnS/hepD rule failed"})])
ex.to_csv("refs/excluded_reps.tsv", sep="\t", index=False)
print("mpnS/hepD name-query refs failing the rule:", len(bad), "of which representatives:", len(bad_reps))
