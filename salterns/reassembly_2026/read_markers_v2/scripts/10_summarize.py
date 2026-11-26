#!/usr/bin/env python
"""Step 10 (v2): v1 10_summarize.py with the v2 labels, and the output contract of 2026-11-25 (same file names and v1
columns, new columns appended at the end).

Definitions (unchanged from v1; per sample s; reads = single mates, both mates counted):
  L_nt = 3 x median aa length of the label's complete references (refs/label_lengths.tsv; v1 overrides
         meta/length_overrides.tsv); r_s = mean QC read length; L_eff = L_nt + r_s - 6 x LMIN + 1;
  depth = reads x r_s / L_eff; GE depth = mean depth of the 8 marker families; copies_per_GE = depth / GE depth;
  RPKM = reads / (L_nt / 1e3) / (total QC reads / 1e6).
Appended to per_sample_KO.tsv (v2):
  ko                                 the KEGG KO(s) of a counted target / marker label (empty for decoy and merged rows)
  reads_archaeal, copies_per_GE_archaeal   reads whose best-hit representative is archaeal (refs/rep_info.tsv.gz
                                     'domain': GTDB R232 for salterns MAG and KEGG members, UniProt lineage, or transferred);
                                     copies_per_GE_archaeal = copies_per_GE x reads_archaeal / reads
  reads_methanogen_lineage, copies_per_GE_methanogen_lineage   same for best hits in a methanogen lineage
                                     (05b_rep_info_v2.py 'archaeal_lineage' == methanogen; ANME excluded)
  reads_bacterial, reads_eukaryotic, reads_domain_unknown
Merged reporting classes: the v1 three plus the v2 merges decided by the validation (MERGED_V2 below).
Outputs in results/: per_sample_KO.tsv, genome_equivalents.tsv, per_guild_read_vs_MAG.tsv.
"""
import os, json, gzip, glob, sys
from collections import defaultdict
import numpy as np
import pandas as pd

W = "/scratch1/afreiburger/salterns_2026/read_markers_v2"
V1 = "/scratch1/afreiburger/salterns_2026/read_markers"
R = "/scratch1/afreiburger/salterns_2026/rollout"
os.chdir(W)
OUT = "results"
cut = pd.read_csv("cutoffs.tsv", sep="\t", dtype=str)
LMIN = int(dict(zip(cut.param, cut.value))["LMIN"])
MARKERS = ["uL2", "uL11", "uS3", "uS11", "uS12", "uS9", "uL4", "uL16"]
samples = os.environ["SAMPLES"].split() if os.environ.get("SAMPLES") else sorted(os.listdir(f"{R}/reads"))
tot = pd.read_csv(f"{V1}/meta/qc_read_totals.tsv", sep="\t", header=None, index_col=0)[1]
rlen = {}
for s in samples:
    q = json.load(open(f"{R}/reads/{s}/read_counts.json"))["qtrim"]
    assert q["result_reads"] == tot[s], s
    rlen[s] = q["result_bases"] / q["result_reads"]
LL = pd.read_csv("refs/label_lengths.tsv", sep="\t").set_index("label")
for lab, aa in pd.read_csv(f"{V1}/meta/length_overrides.tsv", sep="\t").set_index("label").median_aa.items():
    LL.at[lab, "median_aa"] = aa
panel = pd.read_csv("meta/ko_panel_v2.tsv", sep="\t", dtype=str).fillna("")
KO_OF = panel[panel.role.isin(["target", "marker"])].groupby("gene").ko.agg(lambda x: ",".join(sorted(x))).to_dict()
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
for g, grp in zip(panel.gene, panel.group):          # v2 labels: the panel's compound|process
    GROUP.setdefault(g, grp)
GROUP["hepD"] = "phosphonates|2-hydroxyethylphosphonate dioxygenase"   # v2: counted target (v1: decoy class)

# ---------------------------------------------------------------- domain of the best hit, per read
ri = pd.read_csv("refs/rep_info.tsv.gz", sep="\t", dtype=str, usecols=["rep", "domain", "archaeal_lineage"]).fillna("").set_index("rep")
dom_rows = []
for s in samples:
    rd = pd.read_csv(f"results/per_sample/{s}.reads.tsv.gz", sep="\t", dtype={"read": str, "mate": str}, usecols=["label", "rep"])
    rd = rd.join(ri, on="rep")
    rd["dom"] = np.select([rd.domain == "Archaea", rd.domain == "Bacteria", rd.domain.isin(["Eukaryota", "Eukaryotes"])],
                          ["archaeal", "bacterial", "eukaryotic"], "domain_unknown")
    c = rd.groupby(["label", "dom"]).size().unstack(fill_value=0)
    meth = rd[rd.archaeal_lineage == "methanogen"].groupby("label").size()
    for lab, r in c.iterrows():
        dom_rows.append(dict(sample=s, label=lab, reads_archaeal=int(r.get("archaeal", 0)), reads_bacterial=int(r.get("bacterial", 0)),
                             reads_eukaryotic=int(r.get("eukaryotic", 0)), reads_domain_unknown=int(r.get("domain_unknown", 0)),
                             reads_methanogen_lineage=int(meth.get(lab, 0))))
    print("domain", s, flush=True)
dom = pd.DataFrame(dom_rows)

# ---------------------------------------------------------------- per sample x label
rows = []
for s in samples:
    c = pd.read_csv(f"results/per_sample/{s}.counts.tsv", sep="\t")
    for _, r in c.iterrows():
        rows.append(dict(sample=s, label=r.label, role=r.role, reads=int(r.reads), ambiguous_best=int(r.ambiguous_best),
                         low_identity_best=int(r.low_identity_best)))
df = pd.DataFrame(rows)
labels = sorted(set(df.label[df.role.isin(["target", "marker", "decoy_reported"])]) | set(MARKERS)
                | set(LL.index[LL.role.isin(["target", "marker", "decoy_reported"])]))
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
# counted panel targets without any reference in the DB (AOX K17066: KEGG lists only eukaryote genes, which form
# AOX_euk, and no salterns protein passes the profile): explicit zero rows, so that every panel KO has a row
missing = sorted(set(panel.gene[panel.role == "target"]) - set(df.label))
if missing:
    z = pd.DataFrame([dict(sample=s, label=l, role="target", reads=0, ambiguous_best=0, low_identity_best=0,
                           group=GROUP.get(l, "") + " (no prokaryotic reference in the DB)", read_len=rlen[s],
                           total_qc_reads=tot[s], depth=0.0, RPKM=0.0, copies_per_GE=0.0, copies_per_GE_mediannorm=0.0)
                      for s in samples for l in missing])
    df = pd.concat([df, z], ignore_index=True)
    print("targets without references (zero rows):", missing)
# v2 appended columns
df["ko"] = [KO_OF.get(l, "") if r in ("target", "marker") else "" for l, r in zip(df.label, df.role)]
df = df.merge(dom, on=["sample", "label"], how="left")
for c in ["reads_archaeal", "reads_bacterial", "reads_eukaryotic", "reads_domain_unknown", "reads_methanogen_lineage"]:
    df[c] = df[c].fillna(0).astype(int)
assert (df.reads_archaeal + df.reads_bacterial + df.reads_eukaryotic + df.reads_domain_unknown == df.reads).all()
frac = lambda x: np.where(df.reads > 0, x / df.reads.where(df.reads > 0, 1), 0.0)
df["copies_per_GE_archaeal"] = df.copies_per_GE * frac(df.reads_archaeal)
df["copies_per_GE_methanogen_lineage"] = df.copies_per_GE * frac(df.reads_methanogen_lineage)
# merged reporting classes
MERGED = {"MttB_superfamily(mttB+mtgB+other)": ["mttB", "mtgB", "mttB_family_other"],
          "GrdB_family(grdH+grdB)": ["grdH", "grdB"],
          "CuMMO_all(pmoA+amoA+other)": ["pmoA", "amoA", "cuMMO_other"]}
MERGED_V2 = json.load(open("meta/merged_v2.json")) if os.path.exists("meta/merged_v2.json") else {}
MERGED.update(MERGED_V2)
add = []
for name, mem in MERGED.items():
    sub = df[df.label.isin(mem)]
    for smp, g in sub.groupby("sample"):
        add.append(dict(sample=smp, label=name, role="merged", reads=g.reads.sum(), group="merged",
                        depth=g.depth.sum(), copies_per_GE=g.copies_per_GE.sum(),
                        copies_per_GE_mediannorm=g.copies_per_GE_mediannorm.sum(), RPKM=g.RPKM.sum(),
                        total_qc_reads=g.total_qc_reads.iloc[0], read_len=g.read_len.iloc[0], ko="",
                        reads_archaeal=g.reads_archaeal.sum(), copies_per_GE_archaeal=g.copies_per_GE_archaeal.sum(),
                        reads_methanogen_lineage=g.reads_methanogen_lineage.sum(),
                        copies_per_GE_methanogen_lineage=g.copies_per_GE_methanogen_lineage.sum(),
                        reads_bacterial=g.reads_bacterial.sum(), reads_eukaryotic=g.reads_eukaryotic.sum(),
                        reads_domain_unknown=g.reads_domain_unknown.sum()))
df = pd.concat([df, pd.DataFrame(add)], ignore_index=True)
V1COLS = ["sample", "label", "role", "reads", "ambiguous_best", "low_identity_best", "group", "L_nt", "read_len", "L_eff",
          "total_qc_reads", "depth", "RPKM", "copies_per_GE", "copies_per_GE_mediannorm"]
NEW = ["ko", "reads_archaeal", "copies_per_GE_archaeal", "reads_methanogen_lineage", "copies_per_GE_methanogen_lineage",
       "reads_bacterial", "reads_eukaryotic", "reads_domain_unknown"]
df[V1COLS + NEW].to_csv(f"{OUT}/per_sample_KO.tsv", sep="\t", index=False, float_format="%.6g")

# ---------------------------------------------------------------- MAG side per guild (unchanged from v1)
g = pd.read_csv(f"{V1}/meta/guilds_final585.csv")
rc = pd.read_csv(f"{V1}/meta/final585_read_count.tsv", sep="\t", index_col=0)
fd = pd.read_csv(f"{V1}/meta/final_dereplicated.tsv", sep="\t")
fd["mag"] = fd.genome.str.split("|").str[1]
size = dict(zip(fd.mag, fd.size_bp))
GUILD_GENES = {"methanogen": ["mcrA"], "methylamine_user": ["mtbB", "mtmB"], "mepn_consumer": ["phnJ"],
               "betaine_reducer": ["grdI", "grdH"], "tma_producer": ["cutC", "cntA"],
               "tma_candidate_one_sided": ["cutC", "cntA"]}
PROXY = {"methanogen": ("mean", ["mcrA"]), "methylamine_user": ("mean", ["mtbB", "mtmB"]), "mepn_consumer": ("mean", ["phnJ"]),
         "betaine_reducer": ("mean", ["grdI"]), "tma_producer": ("sum", ["cutC", "cntA"]),
         "tma_candidate_one_sided": ("sum", ["cutC", "cntA"])}
out = []
for s in samples:
    ge_rpb = ge.at[s, "GE_depth_mean"] / rlen[s]
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
