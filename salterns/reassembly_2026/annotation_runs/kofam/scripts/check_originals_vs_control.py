#!/usr/bin/env python
"""Validate the originals run (originals_markers.py: KOfam 2026-10 markers, local Prodigal) against the earlier
control test on the same 38 original MAGs (annotation_choice/kofam_control_test: KOfam 2026-05, >= 2,500 bp).
Both use the same local Prodigal build and best-KO-among-markers, so every marker-count difference should come
from a profile/threshold change between the releases. For each differing (MAG, marker) cell it reports the
2026-05 and 2026-10 thresholds and the 2026-10 scores of the genes involved.
Writes results/originals_310/control_test_comparison.tsv. Run: ~/Documents/py_venv/bin/python -I <this>
"""
import re
from pathlib import Path

import pandas as pd

K = Path(__file__).resolve().parent.parent
S = K.parents[2]
CT = K.parents[1] / "annotation_choice/kofam_control_test"
PAT = re.compile(r"Salt_Pond_MetaG_?(.+?)_MG_DASTool_bins[._](concoct_out|metabat|maxbin)\.(\d+)")


def key(n):
    m = PAT.search(str(n))
    return f"{m[1]}|{m[2]}|{int(m[3])}" if m else None


mk = pd.read_csv(K / "meta/marker_ko_table.tsv", sep="\t")
sym2ko = dict(zip(mk.gene, mk.ko))
ct = pd.read_csv(CT / "marker_calls_by_mag.tsv", sep="\t")
ct = ct[(ct.min_contig == 2500) & ~ct.mag.str.startswith("new_")].set_index("mag").drop(columns="min_contig")
gr = pd.read_csv(S / "correlational/r232/guilds_r232.csv")
lab2key = {r.iterativeID_r232: key(r.mag) for r in gr.itertuples()}
h = pd.read_csv(K / "results/originals_310/marker_hits.tsv.gz", sep="\t", dtype={"threshold": str})
h["k"] = h.mag.map(key)
best = h[h.is_best_ko == 1]
cnt = best.groupby(["k", "marker"]).size().unstack(fill_value=0)
rows = []
for lab, r in ct.iterrows():
    k = lab2key[lab]
    for sym, v in r.items():
        new = int(cnt.loc[k, sym]) if (k in cnt.index and sym in cnt.columns) else 0
        if new != int(v):
            ko = sym2ko[sym]
            t = mk.set_index("gene").loc[sym]
            sc = h[(h.k == k) & (h.ko == ko)].sort_values("score", ascending=False)
            rows.append(dict(label=lab, marker=sym, ko=ko, n_2026_05=int(v), n_2026_10=new,
                             thr_2026_05=f"{t['kofam_threshold_2026-05']} {t['kofam_score_type_2026-05']}",
                             thr_2026_10=f"{t['kofam_threshold_2026-10']} {t['kofam_score_type_2026-10']}",
                             threshold_changed=(str(t["kofam_threshold_2026-05"]) != str(t["kofam_threshold_2026-10"]))
                             or (t["kofam_score_type_2026-05"] != t["kofam_score_type_2026-10"]),
                             scores_2026_10=";".join(f"{s:.1f}{'*' if b else ''}" for s, b in zip(sc.score, sc.is_best_ko))))
d = pd.DataFrame(rows)
d.to_csv(K / "results/originals_310/control_test_comparison.tsv", sep="\t", index=False)
n_cells = ct.size
print(f"{len(ct)} control MAGs x {ct.shape[1]} markers = {n_cells} cells; differing {len(d)}")
if len(d):
    print(d.groupby(["marker", "threshold_changed"]).agg(n=("label", "size"), up=("n_2026_10", lambda x: int((x > d.loc[x.index, 'n_2026_05']).sum()))).to_string())
    print("differences at markers whose threshold did NOT change:")
    print(d[~d.threshold_changed].to_string())
