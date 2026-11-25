#!/usr/bin/env python
"""Step 10c: the key tables.
  results/key_unrestored_cores.tsv   12 samples of the six unrestored cores (R1_A/B/C, R2_A/B/C; D1 and D2) x guild:
                                     MAG read fraction, MAG genomes per GE, read-based copies per GE, fold, and the
                                     binned/unbinned/unassembled split of the guild's proxy-gene reads
  results/key_unrestored_cores_summary.tsv  mean / median / min / max over the 12 samples
  results/mcrA_taxonomy_by_sample.tsv       mcrA reads per sample by genus (best hit), identity bin and placement
  results/mcrA_taxonomy_unrestored_pooled.tsv  the same pooled over the 12 samples
  results/CP_lyase_taxonomy_unrestored_pooled.tsv  phnJ (and phnGHIJKLM) reads pooled over the 12 samples by genus
  results/split_per_guild.tsv               binned / unbinned / unassembled read fractions of each guild proxy gene
                                            (24 sediment samples and 15 waters, pooled per site group)
"""
import os
import numpy as np
import pandas as pd

os.chdir("/scratch1/afreiburger/salterns_2026/read_markers")
KEY = [f"{p}_{c}_{d}" for p in ("R1", "R2") for c in "ABC" for d in ("D1", "D2")]
gv = pd.read_csv("results/per_guild_read_vs_MAG.tsv", sep="\t")
sp = pd.read_csv("results/binned_unbinned_split.tsv", sep="\t")
tx = pd.read_csv("results/read_taxonomy.tsv", sep="\t", dtype={"genus": str}).fillna("")
att = pd.read_csv("results/mag_attributable_reads.tsv", sep="\t")
ko = pd.read_csv("results/per_sample_KO.tsv", sep="\t")
PROXY = {"methanogen": ["mcrA"], "methylamine_user": ["mtbB", "mtmB"], "mepn_consumer": ["phnJ"],
         "betaine_reducer": ["grdI"], "tma_producer": ["cutC", "cntA"], "tma_candidate_one_sided": ["cutC", "cntA"]}

rows = []
for _, r in gv[gv["sample"].isin(KEY)].iterrows():
    d = r.to_dict()
    genes = PROXY[r.guild]
    s = sp[(sp["sample"] == r["sample"]) & sp.label.isin(genes)]
    n = s.reads.sum()
    for k in ["kept_mag", "other_bin", "unbinned", "unassembled"]:
        d[f"proxy_reads_frac_{k}"] = s[k].sum() / n if n else np.nan
    d["proxy_reads"] = n
    a = att[(att["sample"] == r["sample"]) & (att.guild == r.guild) & att.label.isin(genes)].reads_attributable.sum()
    d["proxy_reads_attributable_to_guild_MAGs"] = a
    d["proxy_frac_attributable_to_guild_MAGs"] = a / n if n else np.nan
    rows.append(d)
k = pd.DataFrame(rows)
if len(k):
    mg = tx[(tx.label == "mcrA") & tx["sample"].isin(KEY)]
    meth = mg[mg.genus == "Methanolobus"].groupby("sample").reads.sum()
    tot = mg.groupby("sample").reads.sum()
    k["mcrA_frac_genus_Methanolobus_any_identity"] = [meth.get(s, 0) / tot.get(s, np.nan) if g == "methanogen" else np.nan
                                                      for s, g in zip(k["sample"], k.guild)]
k.to_csv("results/key_unrestored_cores.tsv", sep="\t", index=False, float_format="%.5g")
num = [c for c in k.columns if c not in ("sample", "guild", "read_based_proxy") and k[c].dtype != object]
summ = k.groupby("guild")[num].agg(["mean", "median", "min", "max"])
summ.columns = [f"{a}_{b}" for a, b in summ.columns]
summ.to_csv("results/key_unrestored_cores_summary.tsv", sep="\t", float_format="%.5g")

# mcrA taxonomy
m = tx[tx.label == "mcrA"]
m.groupby(["sample", "domain", "phylum", "family", "genus", "genus_source", "ident_bin", "status"]).reads.sum() \
    .reset_index().to_csv("results/mcrA_taxonomy_by_sample.tsv", sep="\t", index=False)
mk = m[m["sample"].isin(KEY)]
pool = mk.groupby(["family", "genus"]).reads.sum().sort_values(ascending=False).reset_index()
pool["frac"] = pool.reads / pool.reads.sum()
by_id = mk.pivot_table(index=["family", "genus"], columns="ident_bin", values="reads", aggfunc="sum", fill_value=0)
by_st = mk.pivot_table(index=["family", "genus"], columns="status", values="reads", aggfunc="sum", fill_value=0)
pool = pool.merge(by_id.reset_index(), on=["family", "genus"]).merge(by_st.reset_index(), on=["family", "genus"], suffixes=("", "_st"))
pool.to_csv("results/mcrA_taxonomy_unrestored_pooled.tsv", sep="\t", index=False, float_format="%.4g")
for lab_set, name in ((["phnJ"], "phnJ"), (["phnG", "phnH", "phnI", "phnJ", "phnK", "phnL", "phnM"], "phnGHIJKLM")):
    c = tx[tx.label.isin(lab_set) & tx["sample"].isin(KEY)]
    p = c.groupby(["phylum", "family", "genus"]).reads.sum().sort_values(ascending=False).reset_index()
    p["frac"] = p.reads / p.reads.sum()
    p.to_csv(f"results/CP_lyase_taxonomy_unrestored_pooled_{name}.tsv", sep="\t", index=False, float_format="%.4g")

# binned / unbinned split per guild proxy gene, by habitat
sp["habitat"] = np.where(sp["sample"].str.contains("H2O"), "water", "sediment")
sp["unrestored_core"] = sp["sample"].isin(KEY)
rows = []
labs = sorted(set(sum(PROXY.values(), [])) | {"mcrB", "mcrG", "phnG", "phnH", "phnI", "phnK", "phnL", "phnM", "grdH", "grdB",
                                               "mttB", "mtgB", "mtaB", "mtaA", "mtbA", "cutD", "cntB", "mpnS", "pmoA", "mmoX",
                                               "uL2", "uL11", "uS3", "uS11", "uS12", "uS9", "uL4", "uL16"})
for grp, g in [("unrestored_cores_12", sp[sp.unrestored_core]), ("all_sediment_24", sp[sp.habitat == "sediment"]),
               ("water_15", sp[sp.habitat == "water"])]:
    for lab in labs:
        x = g[g.label == lab]
        n = x.reads.sum()
        d = dict(group=grp, label=lab, reads=n)
        for c in ["kept_mag", "other_bin", "unbinned", "unassembled"]:
            d[f"frac_{c}"] = x[c].sum() / n if n else np.nan
        rows.append(d)
pd.DataFrame(rows).to_csv("results/split_per_gene_by_habitat.tsv", sep="\t", index=False, float_format="%.4g")
print(k[["sample", "guild", "mag_genomes_per_GE", "read_based_copies_per_GE", "fold_read_vs_MAG"]].to_string())

# assembly-depth complement vs the read survey: assembled copies per GE (KOfam-called genes x own-library contig depth)
# by contig status, as a fraction of the read-based copies per GE
ad = pd.read_csv("results/assembly_depth_complement.tsv", sep="\t")
p = ad.pivot_table(index=["sample", "label"], columns="status", values="assembled_copies_per_GE", aggfunc="sum", fill_value=0)
p.columns = [f"asm_copies_per_GE_{c}" for c in p.columns]
p = p.reset_index()
rb = ko[["sample", "label", "copies_per_GE", "reads"]].rename(columns={"copies_per_GE": "read_copies_per_GE", "reads": "read_count"})
x = rb.merge(p, on=["sample", "label"], how="left").fillna(0)
cs = [c for c in x.columns if c.startswith("asm_copies_per_GE_")]
x["asm_copies_per_GE_total"] = x[cs].sum(axis=1)
x["asm_over_read_ratio"] = np.where(x.read_copies_per_GE > 0, x.asm_copies_per_GE_total / x.read_copies_per_GE, np.nan)
x.to_csv("results/assembly_vs_reads.tsv", sep="\t", index=False, float_format="%.5g")
