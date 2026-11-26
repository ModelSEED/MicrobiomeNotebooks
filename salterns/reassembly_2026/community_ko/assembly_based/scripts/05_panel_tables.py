#!/usr/bin/env python
"""Step 5 (local): panel tables for the KOs of ../ko_panel_methylated.tsv (138 KOs).

results/panel_status_split.tsv  per panel KO x sample: assembled depth on kept-MAG / other-bin / unbinned contigs
    (depth-weighted shares), with the binned part split by domain (kept MAG: GTDB R232 domain of the MAG; other
    DAS Tool bin: domain from its CheckM1 lineage, mapped on the kept bins, results/checkm1_lineage_domain_map.tsv;
    unbinned: no domain). Copies per read-based GE for the total and for the archaeal / bacterial binned parts.
results/panel_mag_attribution.tsv   per panel KO x species of the final 585 that carries it (KOfam best-KO calls of
    the MAG tracks: kofam_rollout for rollout MAGs, the first-pass kofam track for co-assembly MAGs, kofam_coasm38
    for the 38 co-assembly MAGs missing there), with GTDB taxonomy and the mean read fraction in the 12
    unrestored-core samples (final585_relabund_24sediment.tsv).
results/panel_mag_attribution_summary.tsv   per panel KO: number of carrier species (rep / any MAG), summed mean
    read fraction of the carriers, top carriers.
usage: 05_panel_tables.py   (run from assembly_based/)
"""
import numpy as np
import pandas as pd

RES = __import__("os").environ.get("RES", "results")   # results dir (test runs: another dir)
RA = "../.."                                    # reassembly_2026
SAL = "../../.."                                # salterns
panel = pd.read_csv("../ko_panel_methylated.tsv", sep="\t", dtype=str).fillna("")
assert panel.ko.is_unique
info = pd.read_csv(f"{RA}/annotation_runs/kofam_rollout/results/hmm_info.tsv", sep="\t", dtype=str).set_index("ko")
assert panel.ko.isin(info.index).all()
panel["kofam_threshold"] = panel.ko.map(info.threshold)
PCOLS = ["ko", "gene", "compound", "process", "status", "kofam_threshold"]
CORE12 = [f"{a}_{b}_{d}" for a in ("R1", "R2") for b in "ABC" for d in ("D1", "D2")]


def group(s):
    if s in CORE12:
        return "unrestored_core"
    if "H2O" in s:
        return "water"
    if s.split("_")[0] in ("SF2", "R2A") and s.split("_")[-1] in ("D1", "D2"):
        return "other_sediment_core"
    return "R2_restored_other"


# ------------------------------------------------------------------ status / domain split
K = pd.read_csv(f"{RES}/ko_by_sample.tsv.gz", sep="\t")
SD = pd.read_csv(f"{RES}/ko_status_domain_by_sample.tsv.gz", sep="\t")
ge = pd.read_csv(f"{RES}/ge_assembly.tsv", sep="\t", index_col=0)
samples = sorted(ge.index)
cats = ["kept_mag_Archaea", "kept_mag_Bacteria", "other_bin_Archaea", "other_bin_Bacteria", "other_bin_unknown",
        "unbinned"]
assert set(SD.status_domain) <= set(cats), set(SD.status_domain) - set(cats)
w = SD[SD.ko.isin(panel.ko)].pivot_table(index=["ko", "sample"], columns="status_domain", values="sum_depth",
                                         aggfunc="sum", fill_value=0).reindex(columns=cats, fill_value=0)
full = pd.MultiIndex.from_product([panel.ko, samples], names=["ko", "sample"])
w = w.reindex(full, fill_value=0)
w.columns = [f"depth_{c}" for c in cats]
w = w.reset_index()
kk = K.set_index(["ko", "sample"])[["n_genes", "sum_depth"]].reindex(full, fill_value=0).reset_index()
P = panel[PCOLS].merge(kk, on="ko").merge(w, on=["ko", "sample"])
P.insert(P.columns.get_loc("sample") + 1, "sample_group", P["sample"].map(group))
P["depth_kept_mag"] = P.depth_kept_mag_Archaea + P.depth_kept_mag_Bacteria
P["depth_other_bin"] = P.depth_other_bin_Archaea + P.depth_other_bin_Bacteria + P.depth_other_bin_unknown
assert np.allclose(P.depth_kept_mag + P.depth_other_bin + P.depth_unbinned, P.sum_depth, rtol=1e-4, atol=1e-4)   # both tables are written with 6 significant digits
tot = P.sum_depth.replace(0, np.nan)
for c in ["kept_mag", "other_bin", "unbinned"]:
    P[f"share_{c}"] = P[f"depth_{c}"] / tot
binned = (P.depth_kept_mag + P.depth_other_bin).replace(0, np.nan)
P["binned_share_archaea"] = (P.depth_kept_mag_Archaea + P.depth_other_bin_Archaea) / binned
P["binned_share_bacteria"] = (P.depth_kept_mag_Bacteria + P.depth_other_bin_Bacteria) / binned
P["binned_share_domain_unknown"] = P.depth_other_bin_unknown / binned
geR = P["sample"].map(ge.GE_read_mean)
P["copies_per_GE_readGE"] = P.sum_depth / geR
P["copies_per_GE_readGE_binned_archaea"] = (P.depth_kept_mag_Archaea + P.depth_other_bin_Archaea) / geR
P["copies_per_GE_readGE_binned_bacteria"] = (P.depth_kept_mag_Bacteria + P.depth_other_bin_Bacteria) / geR
P["copies_per_GE_readGE_kept_mag_archaea"] = P.depth_kept_mag_Archaea / geR
P["note"] = np.where(P.kofam_threshold == "-", "KOfam threshold '-': never called by the KOfam rule", "")
order = PCOLS + ["sample", "sample_group", "n_genes", "sum_depth", "copies_per_GE_readGE", "depth_kept_mag",
                 "depth_other_bin", "depth_unbinned", "share_kept_mag", "share_other_bin", "share_unbinned"] + \
        [f"depth_{c}" for c in cats[:5]] + ["binned_share_archaea", "binned_share_bacteria",
                                            "binned_share_domain_unknown", "copies_per_GE_readGE_binned_archaea",
                                            "copies_per_GE_readGE_binned_bacteria",
                                            "copies_per_GE_readGE_kept_mag_archaea", "note"]
P[order].to_csv(f"{RES}/panel_status_split.tsv", sep="\t", index=False, float_format="%.6g")
print("panel_status_split", P.shape)

# ------------------------------------------------------------------ MAG attribution
g = pd.read_csv(f"{SAL}/correlational/final585/guilds_final585.csv")
m3 = pd.read_csv(f"{RA}/comparison_rollout/mags_three_way.tsv", sep="\t", dtype=str)
mem = m3[m3.set.isin(["rollout", "coassembly"])].copy()
mem["mag"] = mem.genome.str.split("|", n=1).str[1]
assert len(mem) == 1455 and set(mem.species_cluster) == set(g.species_cluster)


def norm_first_pass(n):         # 'concoct.bin.002.fasta__R1_sediment' -> 'R1_sediment__concoct.bin.002'
    b, unit = n.split(".fasta__")
    return f"{unit}__{b}"


kr = pd.read_csv(f"{RA}/annotation_runs/kofam_rollout/results/all_mag_ko_counts.tsv", sep="\t")
kf = pd.read_csv(f"{RA}/annotation_runs/kofam/results/all_mag_ko_counts.tsv", sep="\t")
kf["mag"] = kf.mag.map(norm_first_pass)
kc = pd.read_csv(f"{RA}/annotation_runs/kofam_coasm38/results/all_mag_ko_counts.tsv", sep="\t")
coas = set(mem.mag[mem.set == "coassembly"])
rol = set(mem.mag[mem.set == "rollout"])
assert rol <= set(kr.mag), len(rol - set(kr.mag))
src = pd.concat([kr[kr.mag.isin(rol)].assign(source="kofam_rollout"),
                 kf[kf.mag.isin(coas)].assign(source="kofam_first_pass"),
                 kc[kc.mag.isin(coas) & ~kc.mag.isin(set(kf.mag))].assign(source="kofam_coasm38")])
have = set(src.mag)
missing = sorted(set(mem.mag) - have)
print("MAGs with KOfam calls", len(have & set(mem.mag)), "of", len(mem), "missing", len(missing), missing[:5])
src = src[src.ko.isin(panel.ko) & (src.n_best > 0)]
src = src.merge(mem[["mag", "species_cluster"]], on="mag")
rel = pd.read_csv(f"{RA}/abundance_final/final585_relabund_24sediment.tsv", sep="\t", index_col=0)
g["mean_read_fraction_12_unrestored_core_samples"] = rel.loc[g.mag, CORE12].mean(axis=1).values
g = g.set_index("species_cluster")
rows = []
for (ko, sc), d in src.groupby(["ko", "species_cluster"]):
    r = g.loc[sc]
    rep_n = int(d.loc[d.mag == r.mag, "n_best"].sum())
    rows.append(dict(ko=ko, species_cluster=sc, rep_mag=r.mag, rep_set=r.set, rep_n_genes=rep_n,
                     rep_carries=rep_n > 0, n_mags_in_species=int(r.n_mags_in_species),
                     n_mags_with_ko=d.mag.nunique(), domain=r.domain, phylum=r.phylum, family=r.family,
                     genus=r.genus if isinstance(r.genus, str) else "", gtdb_species=r.species if isinstance(r.species, str) else "",
                     completeness=r.completeness, contamination=r.contamination,
                     mean_read_fraction_12_unrestored_core_samples=r.mean_read_fraction_12_unrestored_core_samples,
                     mean_rel_abund_12cores_guilds_table=r.mean_rel_abund_12cores))
A = pd.DataFrame(rows)
A = panel[["ko", "gene", "compound", "process", "status", "kofam_threshold"]].merge(A, on="ko", how="left")
A = A.sort_values(["ko", "mean_read_fraction_12_unrestored_core_samples"], ascending=[True, False])
A.to_csv(f"{RES}/panel_mag_attribution.tsv", sep="\t", index=False, float_format="%.6g")
S = []
for _, p in panel.iterrows():
    d = A[(A.ko == p.ko) & A.species_cluster.notna()]
    dr = d[d.rep_carries == True]
    top = dr.head(5)
    S.append(dict(ko=p.ko, gene=p.gene, compound=p.compound, process=p.process, kofam_threshold=p.kofam_threshold,
                  n_species_rep_carries=len(dr), n_species_any_mag_carries=len(d),
                  n_archaeal_species_rep=int((dr.domain == "Archaea").sum()),
                  sum_mean_read_fraction_12_unrestored_rep_carriers=dr.mean_read_fraction_12_unrestored_core_samples.sum(),
                  top5_rep_carriers=";".join(f"{x.genus or x.family}({x.mean_read_fraction_12_unrestored_core_samples:.2e})"
                                              for x in top.itertuples()),
                  note="KOfam threshold '-': never called in the MAG tracks" if p.kofam_threshold == "-" else ""))
pd.DataFrame(S).to_csv(f"{RES}/panel_mag_attribution_summary.tsv", sep="\t", index=False, float_format="%.6g")
with open(f"{RES}/panel_mag_attribution_check.txt", "w") as fo:
    fo.write(f"final585 MAGs (rollout + coassembly)\t{len(mem)}\nwith KOfam calls\t{len(have & set(mem.mag))}\n"
             f"missing\t{len(missing)}\t{','.join(missing)}\n"
             f"mean read fraction: computed here = simple mean of the 12 unrestored-core samples (R1/R2 A-C D1+D2) of "
             f"final585_relabund_24sediment.tsv; guilds_final585.csv 'mean_rel_abund_12cores' = mean over all 12 cores "
             f"(D1/D2 merged DNA-weighted, all four sites), kept for reference\n")
print("panel_mag_attribution", A.shape)

# ------------------------------------------------------------------ supplement: the five no-threshold panel KOs
NOTHR = ["K16967", "K23996", "K05915", "K18049", "K12905"]
assert set(panel.ko[panel.kofam_threshold == "-"]) == set(NOTHR)
mrows = []
for track, fn in (("kofam_rollout", None), ("kofam_first_pass", norm_first_pass), ("kofam_coasm38", None)):
    t = track.replace("_first_pass", "")
    h = pd.read_csv(f"{RA}/annotation_runs/{t}/results/all_no_threshold_hits.tsv.gz", sep="\t", dtype=str,
                    keep_default_na=False)
    h = h[h.ko.isin(NOTHR)]
    h = h[(h.full_evalue.astype(float) <= 1e-30) & (h.hmm_cov.astype(float) >= 0.5) & (h.gene_best_ko == "")]
    if fn:
        h["mag"] = h.mag.map(fn)
    h = h[h.mag.isin(set(src_mags := (rol if t == "kofam_rollout" else coas)))]
    if t == "kofam_coasm38":
        h = h[~h.mag.isin(set(kf.mag))]
    for (mag, ko), d in h.groupby(["mag", "ko"]):
        mrows.append(dict(ko=ko, mag=mag, source=track, n_genes=len(d)))
M = pd.DataFrame(mrows, columns=["ko", "mag", "source", "n_genes"]).merge(mem[["mag", "species_cluster"]], on="mag")
M["is_species_rep"] = M.mag.isin(set(g.mag))
M = M.join(g[["mag", "domain", "phylum", "genus", "mean_read_fraction_12_unrestored_core_samples"]].rename(
    columns={"mag": "rep_mag"}), on="species_cluster")
M.to_csv(f"{RES}/nothreshold_panel_mag_attribution.tsv", sep="\t", index=False, float_format="%.6g")
print("no-threshold KOs, MAG side (no passing KO, E <= 1e-30, HMM cov >= 0.5):", M.groupby("ko").species_cluster.nunique().to_dict())
import os
if os.path.exists(f"{RES}/nothreshold_panel_by_sample_long.tsv"):
    NL = pd.read_csv(f"{RES}/nothreshold_panel_by_sample_long.tsv", sep="\t", keep_default_na=False)
    out = []
    for variant, sub in (("rule_E1e-30_cov0.5", NL), ("rule_plus_top_scoring_v1", NL[NL.is_top_scoring.astype(str) == "True"])):
        for ko in NOTHR:
            for s in samples:
                d = sub[(sub.ko == ko) & (sub["sample"] == s)]
                tot_d = d.sum_depth.sum()
                out.append(dict(ko=ko, gene=panel.set_index("ko").gene[ko], variant=variant, sample=s,
                                sample_group=group(s), n_genes=int(d.n_genes.sum()), sum_depth=tot_d,
                                copies_per_GE_readGE=tot_d / ge.at[s, "GE_read_mean"],
                                copies_per_GE_assembly=tot_d / ge.at[s, "GE_asm_mean"],
                                share_kept_mag=d[d.status == "kept_mag"].sum_depth.sum() / tot_d if tot_d else np.nan,
                                share_other_bin=d[d.status == "other_bin"].sum_depth.sum() / tot_d if tot_d else np.nan,
                                share_unbinned=d[d.status == "unbinned"].sum_depth.sum() / tot_d if tot_d else np.nan))
    pd.DataFrame(out).to_csv(f"{RES}/nothreshold_panel_by_sample.tsv", sep="\t", index=False, float_format="%.6g")
    print("nothreshold_panel_by_sample written")
