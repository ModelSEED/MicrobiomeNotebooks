#!/usr/bin/env python
"""Step 10b: (1) binned / unbinned / unassembled split of the classified reads, (2) assembly-based complement from the
depth tables, (3) genus of mcrA and C-P lyase reads, (4) MAG attribution of the reads.

(1) Each classified read (results/per_sample/<s>.reads.tsv.gz) is placed with map/<s>.hits.tsv.gz (re-mapping of the
    classified reads to the sample's own panel-carrying contigs, >= 95 % identity):
      kept_mag  = contig of a kept rollout MAG of that sample (rollout/<s>/keep/bins.tsv keep == yes)
      other_bin = contig of a DAS Tool bin that failed the KEEP rule
      unbinned  = contig >= 1 kb in no DAS Tool bin
      unassembled = no >= 95 % alignment to a contig >= 1 kb of the sample's own assembly
(2) Assembly complement: every labelled salterns gene of the sample (KOfam call, refs/labeled_meta.tsv) x its contig's
    own-library mean depth (rollout/<s>/depth/<s>.d1000.txt, jgi_summarize_bam_contig_depths, >= 97 % identity reads;
    binning_test/ for R1_A_D1 and SF2_C_D1), summed per label and contig status, divided by the read-based GE depth
    (results/genome_equivalents.tsv) -> assembled copies per GE, comparable to the read-based copies per GE.
(3) Genus of a read = genus of its best-hit representative (refs/rep_info.tsv.gz: binned salterns members -> GTDB R232;
    external members -> KEGG/UniProt organism; salterns-only unbinned clusters -> transferred from their best external
    hit). Identity bins of the read's best hit: >= 95, 85-95, 70-85, < 70 %.
(4) A read is MAG-attributable to a species cluster if its best-hit cluster contains a protein of a MAG of that
    species (rollout or co-assembly) and its identity is >= 95 %.
"""
import os, gzip
from collections import defaultdict
import numpy as np
import pandas as pd

W = "/scratch1/afreiburger/salterns_2026/read_markers"
R = "/scratch1/afreiburger/salterns_2026/rollout"
BT = "/scratch1/afreiburger/salterns_2026/binning_test"
os.chdir(W)
samples = os.environ["SAMPLES"].split() if os.environ.get("SAMPLES") else sorted(os.listdir(f"{R}/reads"))   # SAMPLES: test subsets
ri = pd.read_csv("refs/rep_info.tsv.gz", sep="\t", dtype=str).fillna("").set_index("rep")
ge = pd.read_csv("results/genome_equivalents.tsv", sep="\t", index_col=0)
kmeta = pd.read_csv("refs/labeled_meta.tsv", sep="\t", dtype=str).fillna("")
st = pd.read_csv("meta/contig_status.tsv.gz", sep="\t", dtype=str).fillna("")
g585 = pd.read_csv("meta/guilds_final585.csv")
m3 = pd.read_csv("meta/mags_three_way.tsv", sep="\t", dtype=str).fillna("")
GUILDS = ["methanogen", "methylamine_user", "mepn_consumer", "betaine_reducer", "tma_producer", "tma_candidate_one_sided"]
guild_clusters = {gd: set(g585[g585[gd] == True].species_cluster) for gd in GUILDS}
TAX_LABELS = {"mcrA", "mcrB", "mcrG", "phnJ", "phnG", "phnH", "phnI", "phnK", "phnL", "phnM", "pmoA", "mmoX", "amoA",
              "mttB", "mtgB", "mttB_family_other", "mtbB", "mtmB", "grdH", "grdI", "cutC", "cntA", "mpnS"}


def ibin(p):
    return ">=95" if p >= 95 else ("85-95" if p >= 85 else ("70-85" if p >= 70 else "<70"))


split_rows, tax_rows, attr_rows, asm_rows, mapcheck = [], [], [], [], []
for s in samples:
    rd = pd.read_csv(f"results/per_sample/{s}.reads.tsv.gz", sep="\t", dtype={"read": str, "mate": str})
    hp = f"map/{s}.hits.tsv.gz"
    if os.path.exists(hp):
        h = pd.read_csv(hp, sep="\t", header=None, names=["read", "mate", "contig", "pos", "mapq", "ident"],
                        dtype={"read": str, "mate": str})
        h = h[h.ident >= 0.95]
        sts = st[st["sample"] == s].set_index("contig")
        h["status"] = h.contig.map(sts.status)
        assert h.status.notna().all()
        rd = rd.merge(h[["read", "mate", "status"]], on=["read", "mate"], how="left")
        rd["status"] = rd.status.fillna("unassembled")
    else:
        rd["status"] = "not_mapped_yet"
    sp = rd.groupby(["label", "status"]).size().unstack(fill_value=0)
    for lab, r in sp.iterrows():
        d = dict(sample=s, label=lab, reads=int(r.sum()))
        for k in ["kept_mag", "other_bin", "unbinned", "unassembled"]:
            d[k] = int(r.get(k, 0))
        split_rows.append(d)
    # taxonomy + attribution
    rd = rd.join(ri[["genus", "family", "phylum", "domain", "genus_source", "species_clusters"]], on="rep")
    rd["ident_bin"] = rd.pident.map(ibin)
    t = rd[rd.label.isin(TAX_LABELS)]
    for key, gg in t.groupby(["label", "domain", "phylum", "family", "genus", "genus_source", "ident_bin", "status"]):
        tax_rows.append(dict(zip(["label", "domain", "phylum", "family", "genus", "genus_source", "ident_bin", "status"], key),
                             sample=s, reads=len(gg)))
    hi = rd[rd.pident >= 95]
    for gd, cl in guild_clusters.items():
        att = hi[hi.species_clusters.map(lambda x: bool(set(x.split(",")) & cl) if x else False)]
        for lab, n in att.groupby("label").size().items():
            attr_rows.append(dict(sample=s, guild=gd, label=lab, reads_attributable=int(n)))
    # assembly complement
    dp = f"{R}/{s}/depth/{s}.d1000.txt"
    if not os.path.exists(dp):
        dp = f"{BT}/{s}/depth/{s}.d1000.txt"
    dep = pd.read_csv(dp, sep="\t", usecols=[0, 2], names=["contig", "depth"], header=0).set_index("contig").depth
    km = kmeta[(kmeta.source == "salterns") & (kmeta["set"] == s) & kmeta.role.isin(["target", "marker", "decoy_reported"])].copy()
    km["contig"] = km.pid.str.split("|", n=1).str[1].str.rsplit("_", n=1).str[0]
    km["depth"] = km.contig.map(dep)
    assert km.depth.notna().all(), s
    for (lab, cs), gg in km.groupby(["label", "contig_status"]):
        asm_rows.append(dict(sample=s, label=lab, status=cs, n_genes=len(gg), sum_depth=gg.depth.sum(),
                             assembled_copies_per_GE=gg.depth.sum() / ge.at[s, "GE_depth_mean"]))
    print(s, flush=True)

sp = pd.DataFrame(split_rows)
for k in ["kept_mag", "other_bin", "unbinned", "unassembled"]:
    sp[f"frac_{k}"] = sp[k] / sp.reads
sp.to_csv("results/binned_unbinned_split.tsv", sep="\t", index=False, float_format="%.5g")
pd.DataFrame(tax_rows).to_csv("results/read_taxonomy.tsv", sep="\t", index=False)
pd.DataFrame(attr_rows).to_csv("results/mag_attributable_reads.tsv", sep="\t", index=False)
pd.DataFrame(asm_rows).to_csv("results/assembly_depth_complement.tsv", sep="\t", index=False, float_format="%.6g")
print("done")
