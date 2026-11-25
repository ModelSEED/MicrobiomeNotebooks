#!/usr/bin/env python
"""Step 10: tables from the per-sample read classifications, the re-mapping of the classified reads and the
assembly depth tables.

Definitions (per sample s; reads = single mates, both mates counted, as in qc_read_totals.tsv):
  L_nt(label)   = 3 x median aa length of the label's complete references (refs/label_lengths.tsv)
  r_s           = mean QC read length of s (rollout reads/<s>/read_counts.json, qtrim result_bases / result_reads)
  L_eff         = L_nt + r_s - 2 * 3 * LMIN + 1   (read start positions giving >= LMIN aligned aa)
  depth         = reads x r_s / L_eff             (summed coverage of all copies of the gene family)
  GE depth      = mean over the 8 marker families of their depth (median reported too)
  copies_per_GE = depth / GE depth
  RPKM          = reads / (L_nt / 1e3) / (total QC reads / 1e6)
MAG side (final 585 species; CoverM read counts, >= 95 % id, >= 75 % aligned):
  genomes_per_GE(species) = (CoverM reads / MAG size bp) / (GE depth / r_s)   i.e. both in reads per bp
Outputs in results/ (see README).
"""
import os, json, gzip, glob, sys
from collections import defaultdict
import numpy as np
import pandas as pd

W = "/scratch1/afreiburger/salterns_2026/read_markers"
R = "/scratch1/afreiburger/salterns_2026/rollout"
os.chdir(W)
OUT = "results"
cut = pd.read_csv("cutoffs.tsv", sep="\t", dtype=str)
LMIN = int(dict(zip(cut.param, cut.value))["LMIN"])
MARKERS = ["uL2", "uL11", "uS3", "uS11", "uS12", "uS9", "uL4", "uL16"]
samples = os.environ["SAMPLES"].split() if os.environ.get("SAMPLES") else sorted(os.listdir(f"{R}/reads"))   # SAMPLES: test subsets
tot = pd.read_csv("meta/qc_read_totals.tsv", sep="\t", header=None, index_col=0)[1]
rlen = {}
for s in samples:
    q = json.load(open(f"{R}/reads/{s}/read_counts.json"))["qtrim"]
    assert q["result_reads"] == tot[s], s
    rlen[s] = q["result_bases"] / q["result_reads"]
LL = pd.read_csv("refs/label_lengths.tsv", sep="\t").set_index("label")
for lab, aa in pd.read_csv("meta/length_overrides.tsv", sep="\t").set_index("label").median_aa.items():
    LL.at[lab, "median_aa"] = aa                      # see meta/length_overrides.tsv
GROUP = {"mcrA": "methanogenesis", "mcrB": "methanogenesis", "mcrG": "methanogenesis",
         "mtaB": "methylotrophic_MT", "mtaA": "methylotrophic_MT", "mtbB": "methylotrophic_MT", "mtmB": "methylotrophic_MT",
         "mtbA": "methylotrophic_MT", "mttB": "mttB_family", "mtgB": "mttB_family", "mttB_family_other": "mttB_family",
         "phnG": "CP_lyase", "phnH": "CP_lyase", "phnI": "CP_lyase", "phnJ": "CP_lyase", "phnK": "CP_lyase",
         "phnL": "CP_lyase", "phnM": "CP_lyase", "grdH": "betaine_reductase", "grdI": "betaine_reductase",
         "grdB": "glycine_reductase(decoy class)", "cutC": "TMA_production", "cutD": "TMA_production",
         "cntA": "TMA_production", "cntB": "TMA_production", "mpnS": "MePn_synthesis", "hepD": "HEPD(decoy class)",
         "pmoA": "methanotrophy", "mmoX": "methanotrophy", "amoA": "ammonia_oxidation(decoy class)",
         "cuMMO_other": "other_CuMMO(decoy class)"}
GROUP.update({m: "marker" for m in MARKERS})

# ---------------------------------------------------------------- per sample x label
rows = []
for s in samples:
    c = pd.read_csv(f"results/per_sample/{s}.counts.tsv", sep="\t")
    for _, r in c.iterrows():
        rows.append(dict(sample=s, label=r.label, role=r.role, reads=int(r.reads), ambiguous_best=int(r.ambiguous_best),
                         low_identity_best=int(r.low_identity_best)))
df = pd.DataFrame(rows)
labels = sorted(set(df.label[df.role.isin(["target", "marker", "decoy_reported"])]) | set(MARKERS))
full = pd.MultiIndex.from_product([samples, labels], names=["sample", "label"])
df = df.set_index(["sample", "label"]).reindex(full).reset_index()
df["role"] = df.label.map(LL.role)
df[["reads", "ambiguous_best", "low_identity_best"]] = df[["reads", "ambiguous_best", "low_identity_best"]].fillna(0).astype(int)
df["group"] = df.label.map(GROUP).fillna("")
df["L_nt"] = df.label.map(LL.median_aa) * 3
df["read_len"] = df["sample"].map(rlen)
df["L_eff"] = df.L_nt + df.read_len - 6 * LMIN + 1
df["total_qc_reads"] = df["sample"].map(tot)
df["depth"] = df.reads * df.read_len / df.L_eff
df["RPKM"] = df.reads / (df.L_nt / 1e3) / (df.total_qc_reads / 1e6)
mk = df[df.label.isin(MARKERS)].pivot(index="sample", columns="label", values="depth")[MARKERS]
ge = pd.DataFrame({"GE_depth_mean": mk.mean(axis=1), "GE_depth_median": mk.median(axis=1),
                   "GE_depth_cv": mk.std(axis=1) / mk.mean(axis=1)})
for m in MARKERS:
    ge[f"depth_{m}"] = mk[m]
ge["read_len"] = pd.Series(rlen)
ge["total_qc_reads"] = tot
ge.to_csv(f"{OUT}/genome_equivalents.tsv", sep="\t")
df["copies_per_GE"] = df.depth / df["sample"].map(ge.GE_depth_mean)
df["copies_per_GE_mediannorm"] = df.depth / df["sample"].map(ge.GE_depth_median)
# merged reporting classes (labels that the hold-out test could not separate reliably at read level)
MERGED = {"MttB_superfamily(mttB+mtgB+other)": ["mttB", "mtgB", "mttB_family_other"],
          "GrdB_family(grdH+grdB)": ["grdH", "grdB"],
          "CuMMO_all(pmoA+amoA+other)": ["pmoA", "amoA", "cuMMO_other"]}
add = []
for name, mem in MERGED.items():
    sub = df[df.label.isin(mem)]
    for smp, g in sub.groupby("sample"):
        add.append(dict(sample=smp, label=name, role="merged", reads=g.reads.sum(), group="merged",
                        depth=g.depth.sum(), copies_per_GE=g.copies_per_GE.sum(),
                        copies_per_GE_mediannorm=g.copies_per_GE_mediannorm.sum(), RPKM=g.RPKM.sum(),
                        total_qc_reads=g.total_qc_reads.iloc[0], read_len=g.read_len.iloc[0]))
df = pd.concat([df, pd.DataFrame(add)], ignore_index=True)
df.to_csv(f"{OUT}/per_sample_KO.tsv", sep="\t", index=False, float_format="%.6g")

# ---------------------------------------------------------------- MAG side per guild
g = pd.read_csv("meta/guilds_final585.csv")
rc = pd.read_csv("meta/final585_read_count.tsv", sep="\t", index_col=0)
ra = pd.read_csv("meta/final585_relabund_24sediment.tsv", sep="\t", index_col=0)
fd = pd.read_csv("meta/final_dereplicated.tsv", sep="\t")
fd["mag"] = fd.genome.str.split("|").str[1]
size = dict(zip(fd.mag, fd.size_bp))
GUILD_GENES = {"methanogen": ["mcrA"], "methylamine_user": ["mtbB", "mtmB"], "mepn_consumer": ["phnJ"],
               "betaine_reducer": ["grdI", "grdH"], "tma_producer": ["cutC", "cntA"],
               "tma_candidate_one_sided": ["cutC", "cntA"]}
# read-based guild proxy: methanogen mcrA; methylamine_user mean(mtbB, mtmB); mepn_consumer phnJ;
# betaine_reducer grdI (grdH is reported but half of its held-out reads went to grdB); TMA: cutC + cntA
PROXY = {"methanogen": ("mean", ["mcrA"]), "methylamine_user": ("mean", ["mtbB", "mtmB"]), "mepn_consumer": ("mean", ["phnJ"]),
         "betaine_reducer": ("mean", ["grdI"]), "tma_producer": ("sum", ["cutC", "cntA"]),
         "tma_candidate_one_sided": ("sum", ["cutC", "cntA"])}
out = []
for s in samples:
    ge_rpb = ge.at[s, "GE_depth_mean"] / rlen[s]       # marker reads per bp
    for guild, genes in GUILD_GENES.items():
        sp = g[g[guild] == True].mag.tolist()
        reads = rc.loc[sp, s].sum() if sp else 0
        gpe = sum(rc.at[m, s] / size[m] for m in sp) / ge_rpb if sp else 0.0
        d = dict(sample=s, guild=guild, n_species=len(sp), mag_reads=int(reads), mag_read_fraction=reads / tot[s],
                 mag_genomes_per_GE=gpe)
        for gene in genes:
            r = df[(df["sample"] == s) & (df.label == gene)].iloc[0]
            d[f"{gene}_reads"] = r.reads
            d[f"{gene}_copies_per_GE"] = r.copies_per_GE
            d[f"{gene}_read_fraction"] = r.reads / tot[s]
        how, pg = PROXY[guild]
        vals = [d[f"{x}_copies_per_GE"] for x in pg]
        d["read_based_proxy"] = f"{how}({','.join(pg)})"
        d["read_based_copies_per_GE"] = np.mean(vals) if how == "mean" else sum(vals)
        d["fold_read_vs_MAG"] = d["read_based_copies_per_GE"] / gpe if gpe > 0 else np.nan
        out.append(d)
pd.DataFrame(out).to_csv(f"{OUT}/per_guild_read_vs_MAG.tsv", sep="\t", index=False, float_format="%.6g")
print("tables written")
