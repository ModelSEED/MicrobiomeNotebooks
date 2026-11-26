#!/usr/bin/env python
"""Step 6 (local): assembly-based KO profile vs the v1 read survey (read_markers).

For every read-survey label that is one panel KO with a KOfam threshold, or one marker family (uL4, uL16 = two KOs
each):
  read_copies_per_GE      read_markers/results/per_sample_KO.tsv copies_per_GE (reads -> depth / read GE)
  asm_copies_per_GE       this track: summed contig depth of the KO's genes / read-based GE depth
                          (markers: family depth incl. the 0.5 x threshold rescue, results/ge_assembly.tsv)
  assembled_fraction      asm / read
  v1_asm_copies_per_GE    read_markers/results/assembly_depth_complement.tsv, summed over contig status (v1:
                          62-profile panel search, best KO among panel KOs, score >= threshold)
  ratio_vs_v1             asm / v1_asm
Writes results/validation_vs_read_survey.tsv (per sample x label) and results/validation_vs_read_survey_summary.tsv
(per label x sample group: medians / ranges, Spearman rho across samples).
usage: 06_validate_vs_reads.py   (run from assembly_based/)
"""
import os
import warnings
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings("ignore", message="An input array is constant")   # markers: ~1 copy per GE everywhere

RES = __import__("os").environ.get("RES", "results")   # results dir (test runs: another dir)
RM = "../../read_markers"
CORE12 = [f"{a}_{b}_{d}" for a in ("R1", "R2") for b in "ABC" for d in ("D1", "D2")]


def group(s):
    if s in CORE12:
        return "unrestored_core"
    if "H2O" in s:
        return "water"
    if s.split("_")[0] in ("SF2", "R2A") and s.split("_")[-1] in ("D1", "D2"):
        return "other_sediment_core"
    return "R2_restored_other"


kp = pd.read_csv(f"{RM}/meta/ko_panel.tsv", sep="\t", dtype=str)
kp = kp[kp.role.isin(["target", "marker", "decoy_reported"])]
lab2ko = kp.groupby("gene").ko.apply(list).to_dict()
lab2role = kp.groupby("gene").role.first().to_dict()
rs = pd.read_csv(f"{RM}/results/per_sample_KO.tsv", sep="\t")
rs = rs[rs.label.isin(lab2ko)]
v1 = pd.read_csv(f"{RM}/results/assembly_depth_complement.tsv", sep="\t")
v1 = v1.groupby(["sample", "label"]).agg(v1_n_genes=("n_genes", "sum"), v1_sum_depth=("sum_depth", "sum"),
                                        v1_asm_copies_per_GE=("assembled_copies_per_GE", "sum"))
K = pd.read_csv(f"{RES}/ko_by_sample.tsv.gz", sep="\t").set_index(["ko", "sample"])
# two-pass design: ko_by_sample holds the 148 panel + marker KOs; read labels whose KO is outside them (grdB K10672,
# pmoA K10944, mmoX K16157) or has no KOfam threshold (mpnS K18049) are not compared
P148 = set(K.index.get_level_values(0))
thr = K.reset_index().drop_duplicates("ko").set_index("ko").kofam_threshold.astype(str)
skipped = sorted(l for l, ks in lab2ko.items() if not all(k in P148 and thr.get(k, "-") != "-" for k in ks))
rs = rs[~rs.label.isin(skipped)]
print("labels not compared:", skipped)
ge = pd.read_csv(f"{RES}/ge_assembly.tsv", sep="\t", index_col=0)
assert len(ge) == 39 or os.environ.get("TEST"), len(ge)
rs = rs[rs["sample"].isin(ge.index)]
resc = pd.read_csv(f"{RES}/marker_rescue_check.tsv", sep="\t").set_index(["sample", "family"])
rows = []
for r in rs.itertuples():
    kos = lab2ko[r.label]
    if lab2role[r.label] == "marker":
        q = resc.loc[(r.sample, r.label)]
        n = q.n_strict + q.n_rescue_v1panel
        d = ge.at[r.sample, f"depth_{r.label}"]
    else:
        n = sum(K.n_genes.get((k, r.sample), 0) for k in kos)
        d = sum(K.sum_depth.get((k, r.sample), 0.0) for k in kos)
    a = d / ge.at[r.sample, "GE_read_mean"]
    v = v1.loc[(r.sample, r.label)] if (r.sample, r.label) in v1.index else None
    rows.append(dict(sample=r.sample, sample_group=group(r.sample), label=r.label, role=lab2role[r.label],
                     kos="+".join(kos), read_count=r.reads, read_copies_per_GE=r.copies_per_GE, asm_n_genes=n,
                     asm_sum_depth=d, asm_copies_per_GE=a,
                     assembled_fraction=a / r.copies_per_GE if r.copies_per_GE > 0 else np.nan,
                     v1_n_genes=v.v1_n_genes if v is not None else 0,
                     v1_asm_copies_per_GE=v.v1_asm_copies_per_GE if v is not None else 0.0))
T = pd.DataFrame(rows)
T["ratio_vs_v1"] = T.asm_copies_per_GE / T.v1_asm_copies_per_GE.replace(0, np.nan)
T["genes_minus_v1"] = T.asm_n_genes - T.v1_n_genes
T.to_csv(f"{RES}/validation_vs_read_survey.tsv", sep="\t", index=False, float_format="%.6g")
S = []
for (lab, grp), d in pd.concat([T, T.assign(sample_group="all_39")]).groupby(["label", "sample_group"]):
    ok = d[d.read_count >= 20]                       # assembled fraction only where the read count is >= 20
    rho = spearmanr(d.read_copies_per_GE, d.asm_copies_per_GE)[0] if len(d) >= 5 else np.nan
    S.append(dict(label=lab, role=d.role.iloc[0], sample_group=grp, n_samples=len(d),
                  n_samples_reads_ge20=len(ok),
                  read_copies_per_GE_median=d.read_copies_per_GE.median(),
                  asm_copies_per_GE_median=d.asm_copies_per_GE.median(),
                  assembled_fraction_pooled=d.asm_copies_per_GE.sum() / d.read_copies_per_GE.sum()
                  if d.read_copies_per_GE.sum() > 0 else np.nan,
                  assembled_fraction_median=ok.assembled_fraction.median(),
                  assembled_fraction_min=ok.assembled_fraction.min(),
                  assembled_fraction_max=ok.assembled_fraction.max(),
                  spearman_rho_asm_vs_read=rho,
                  ratio_vs_v1_pooled=d.asm_copies_per_GE.sum() / d.v1_asm_copies_per_GE.sum()
                  if d.v1_asm_copies_per_GE.sum() > 0 else np.nan,
                  ratio_vs_v1_min=d.ratio_vs_v1.min(), ratio_vs_v1_max=d.ratio_vs_v1.max(),
                  genes_here=int(d.asm_n_genes.sum()), genes_v1=int(d.v1_n_genes.sum())))
pd.DataFrame(S).to_csv(f"{RES}/validation_vs_read_survey_summary.tsv", sep="\t", index=False, float_format="%.4g")
print(pd.DataFrame(S).query("sample_group == 'all_39'")[["label", "assembled_fraction_pooled",
                                                          "assembled_fraction_min", "assembled_fraction_max",
                                                          "spearman_rho_asm_vs_read", "ratio_vs_v1_pooled",
                                                          "genes_here", "genes_v1"]].to_string(index=False))
