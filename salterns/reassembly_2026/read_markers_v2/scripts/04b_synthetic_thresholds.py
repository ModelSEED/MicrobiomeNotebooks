#!/usr/bin/env python
"""Step 4b (v2): minimum full scores for the counted KOs without a KOfam threshold that are new in v2 (dmoA K16967,
phnY K23996, phnZ K05915), from the external hmmsearch (hmm_new/ext.tbl, -Z 20000000):
  own   = full scores of the KO's own KEGG prokaryote genes on its profile
  other = full scores on that profile of every other KEGG gene of the v2 panel (targets and decoys)
  min_full_score = max(other) + 1, i.e. no KEGG gene of another panel KO would pass (KOfam-like separation);
  the fraction of own genes kept at that score is reported. mpnS / hepD keep the v1 rule (no score floor).
Writes meta/synthetic_thresholds.tsv."""
import os
import numpy as np
import pandas as pd

os.chdir("/scratch1/afreiburger/salterns_2026/read_markers_v2")
KOS = {"K16967": "dmoA", "K23996": "phnY", "K05915": "phnZ"}
# characterised enzymes in the reviewed UniProt sets: DmoA of Hyphomicrobium sulfonivorans (E9JFX9; KEGG K16967 links
# to no genome gene, only GenBank ADU77278 = the same protein)
UP_OWN = {"K16967": ["up|rev_PF00296|E9JFX9"]}
rows = []
for l in open("hmm_new/ext.tbl"):
    if l.startswith("#"):
        continue
    f = l.split()
    if f[2] in KOS and f[0].startswith("kegg|"):
        rows.append((f[0], f[0].split("|")[1], f[2], float(f[4]), float(f[5])))
    elif f[2] in KOS and f[0].startswith("up|rev_"):        # reviewed UniProt: 'own' if annotated as that enzyme
        rows.append((f[0], "uniprot_reviewed", f[2], float(f[4]), float(f[5])))
h = pd.DataFrame(rows, columns=["pid", "gene_ko", "profile", "E", "score"])
meta = pd.read_csv("refs/ext_meta.tsv", sep="\t", dtype=str).set_index("pid")
h["domain"] = h.pid.map(meta.domain)
h = h[h.domain != "Eukaryota"]
out = []
for ko, lab in KOS.items():
    x = h[h.profile == ko]
    up_own = x[(x.gene_ko == "uniprot_reviewed") & x.pid.isin(UP_OWN.get(ko, []))]
    own = pd.concat([x[x.gene_ko == ko].score, up_own.score])
    oth = x[(x.gene_ko != ko) & (x.gene_ko != "uniprot_reviewed")].sort_values("score", ascending=False)
    n_own_total = int((meta.set == ko).sum() - ((meta.set == ko) & (meta.domain == "Eukaryota")).sum()) + len(up_own)
    omax = float(oth.score.iloc[0]) if len(oth) else 0.0
    thr = omax + 1.0
    out.append(dict(ko=ko, label=lab, min_full_score=round(thr, 1), n_own_kegg=n_own_total, n_own_hit=len(own),
                    own_min=own.min(), own_p05=own.quantile(0.05), own_median=own.median(),
                    own_kept_frac=(own >= thr).sum() / n_own_total if n_own_total else np.nan,
                    n_other_hits=len(oth), other_max=omax, other_top_ko=oth.gene_ko.iloc[0] if len(oth) else "",
                    other_top5=";".join(f"{k}:{s:.0f}" for k, s in zip(oth.gene_ko.head(5), oth.score.head(5)))))
o = pd.DataFrame(out)
o.to_csv("meta/synthetic_thresholds.tsv", sep="\t", index=False, float_format="%.1f")
pd.set_option("display.width", 250)
print(o.to_string())
