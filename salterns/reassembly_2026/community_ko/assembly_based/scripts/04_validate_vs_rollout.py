#!/usr/bin/env python3
"""Step 4 (check): genes on kept-MAG contigs vs the kofam_rollout MAG track.

The 935 rollout MAGs are DAS Tool bins of the same per-sample assemblies, and both tracks ran Prodigal 2.6.3
-p meta -g 11 (per-contig models), so every rollout MAG gene '<sample>__<bin>|<contig>_<n>' should exist here as
'<sample>|<contig>_<n>' with the identical protein, and get the identical KOfam call (best KO, score, all passing
KOs; only E-values differ, because -Z differs: 2,732,576 there, 12,055,118 here). Two-pass design: the panel/marker
label is compared for every gene, the complete call for the genes searched against the full DB in pass 2.
Writes results/validation_vs_kofam_rollout.txt and results/validation_vs_kofam_rollout_diffs.tsv.
usage: 04_validate_vs_rollout.py <W>
"""
import gzip
import hashlib
import os
import sys
import numpy as np
import pandas as pd

W = sys.argv[1]
S26 = "/scratch1/afreiburger/salterns_2026"
KR, GD = f"{S26}/kofam_rollout", f"{S26}/read_markers/genes"
os.chdir(W)


def mapname(g):
    mag, rest = g.split("|", 1)
    return mag.split("__", 1)[0] + "|" + rest


def read_faa(p):
    name, seq = None, []
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)


# 1. protein identity
md5_r = {}
for f in sorted(os.listdir(f"{KR}/genes")):
    if f.endswith(".faa"):
        for n, s in read_faa(f"{KR}/genes/{f}"):
            md5_r[mapname(n)] = hashlib.md5(s.encode()).digest()
n_r = len(md5_r)
found = same = 0
samples = sorted(f[:-7] for f in os.listdir("meta/genes") if f.endswith(".tsv.gz"))
for s in samples:
    for n, seq in read_faa(f"{GD}/{s}/{s}.faa"):
        h = md5_r.get(n)
        if h is not None:
            found += 1
            same += h == hashlib.md5(seq.encode()).digest()
print("rollout genes", n_r, "found by name", found, "identical protein", same, flush=True)
del md5_r

# 2. calls
old = pd.read_csv(f"{KR}/results/all_genes_ko.tsv.gz", sep="\t", dtype=str, keep_default_na=False,
                  usecols=["gene", "best_ko", "best_score", "all_passing"])
old["gene"] = old.gene.map(mapname)
old = old.set_index("gene")
ga = pd.read_csv(f"{KR}/genes_all.tsv.gz", sep="\t", dtype=str, usecols=["gene", "mag"])
ga["gene"] = ga.gene.map(mapname)
new = []
for s in samples:
    g = pd.read_csv(f"genes_ko/{s}.tsv.gz", sep="\t", dtype=str, keep_default_na=False,
                    usecols=["gene", "searched_full_db", "ko", "score", "all_passing", "status", "mag_id"])
    new.append(g[g.status == "kept_mag"])
new = pd.concat(new, ignore_index=True).set_index("gene")
P148 = {l.strip() for l in open("meta/pass1_kos.txt") if l.strip()}
cmp_ = ga.set_index("gene").join(new, how="outer").join(old[["best_ko", "best_score"]], how="left")
cmp_["in_rollout"] = cmp_.mag.notna()
cmp_["in_assembly_kept"] = cmp_.status.notna()
for c in ["ko", "score", "all_passing", "best_ko", "best_score"]:
    cmp_[c] = cmp_[c].fillna("")
cmp_["all_passing_old"] = old.all_passing.reindex(cmp_.index).fillna("")
both = cmp_[cmp_.in_rollout & cmp_.in_assembly_kept].copy()
# two-pass design: panel/marker label (best KO if one of the 148 panel + marker KOs, else none) for every gene;
# the complete call (best KO over all 27,981 profiles, score, all passing KOs) for the genes searched in pass 2
both["label_old"] = np.where(both.best_ko.isin(P148), both.best_ko, "")
both["label_new"] = np.where(both.ko.isin(P148), both.ko, "")
both["passes148_old"] = [any(x.split(":")[0] in P148 for x in a.split(";") if x) for a in both.all_passing_old]
srch = both[both.searched_full_db == "1"]
same_label = both.label_old == both.label_new
same_ko = srch.ko == srch.best_ko
same_score = same_ko & (srch.score == srch.best_score)
same_all = srch.all_passing == srch.all_passing_old
missed = both[both.passes148_old & (both.searched_full_db != "1")]
mag_ok = (both.mag == both.mag_id).all()
diffs = pd.concat([both[~same_label], srch[~(same_ko & same_score & same_all)]]).drop_duplicates()
diffs.to_csv("results/validation_vs_kofam_rollout_diffs.tsv", sep="\t")
with open("results/validation_vs_kofam_rollout.txt", "w") as fo:
    fo.write(f"rollout_genes\t{n_r}\nfound_by_name_in_assembly_genes\t{found}\nidentical_protein\t{same}\n"
             f"rollout_genes_on_kept_mag_contigs_here\t{len(both)}\n"
             f"kept_mag_genes_here_not_in_rollout\t{int((~cmp_.in_rollout & cmp_.in_assembly_kept).sum())}\n"
             f"rollout_genes_not_on_kept_mag_contigs_here\t{int((cmp_.in_rollout & ~cmp_.in_assembly_kept).sum())}\n"
             f"same_mag\t{mag_ok}\n"
             f"genes_with_panel_or_marker_ko_rollout\t{int((both.label_old != '').sum())}\t"
             f"here\t{int((both.label_new != '').sum())}\n"
             f"same_panel_marker_label_all_genes\t{int(same_label.sum())}\tof\t{len(both)}\tfraction\t{same_label.mean():.6f}\n"
             f"rollout_genes_passing_a_148_KO_not_searched_here\t{len(missed)}\n"
             f"genes_searched_full_db\t{len(srch)}\n"
             f"searched: same_best_ko\t{int(same_ko.sum())}\tfraction\t{same_ko.mean():.6f}\n"
             f"searched: same_best_ko_and_score\t{int(same_score.sum())}\tfraction\t{same_score.mean():.6f}\n"
             f"searched: same_all_passing\t{int(same_all.sum())}\tfraction\t{same_all.mean():.6f}\n"
             f"genes_differing\t{len(diffs)}\n")
print(open("results/validation_vs_kofam_rollout.txt").read())
