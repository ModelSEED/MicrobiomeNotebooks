#!/usr/bin/env python3
"""Step 3: per-gene KO calls x contig depth -> community KO profile per sample (two-pass design).

Two-pass design (2026-11-25, on the coordinator's request; the full 27,981-profile search of all 9.76 M unique
proteins was stopped after 4 of 521 chunks): pass 1 screened every protein with the 148 panel + marker profiles,
pass 2 searched the candidates (08_select_candidates.py) with all 27,981 profiles. The best-KO calls are exact for
the 138 panel KOs and the 10 marker KOs (and the marker rescue); for any other KO they are partial (only genes that
were candidates were searched against it), so the KO outputs are restricted to the 148 KOs.

Inputs (poplar):
  pass2/results/chunks/chunk_*.calls.tsv.gz  representative-level calls of the candidates (02_postprocess_chunk.py)
  pass2/candidates.tsv.gz                    the candidates (searched against the full DB)
  meta/pass1_kos.txt, meta/ko_panel_methylated.tsv   the 148 KOs (138 panel + 10 marker KOs)
  meta/genes/<s>.tsv.gz                      every gene: contig, partial, aa_len, representative (01_make_chunks.py)
  rollout/<s>/depth/<s>.d1000.txt            own-library mean depth per contig (jgi_summarize_bam_contig_depths,
                                             column totalAvgDepth; single BAM); binning_test/ for R1_A_D1, SF2_C_D1
  read_markers/meta/contig_status.tsv.gz     contig -> kept_mag / other_bin / unbinned, bin, mag_id
  rollout/results/mags_kept.tsv              kept MAG -> GTDB R232 domain
  rollout/<s>/keep/bins.tsv                  other DAS Tool bins -> CheckM1 lineage -> domain (map learnt on kept bins)
  read_markers/results/genome_equivalents.tsv   read-based GE depth (GE_depth_mean)
Outputs:
  genes_ko/<s>.tsv.gz                        per gene (all genes): searched_full_db, KO call (best passing KO over all
                                             27,981 profiles, candidates only), score, contig, depth, status, domain
  results/ko_by_sample.tsv.gz                the 148 panel + marker KOs x 39 samples
  results/ko_status_domain_by_sample.tsv.gz  per sample x KO (148) x contig status x binned domain: n_genes, sum_depth
  results/ge_assembly.tsv                    assembly-internal GE: 8 marker families, mean, median, CV; read GE
  results/marker_rescue_check.tsv            marker genes per sample and family: strict, rescue_all, rescue_panel
  results/checkm1_lineage_domain_map.tsv     lineage -> domain map used for other_bin contigs
  results/aggregate_qc.txt
usage: 03_aggregate.py <W>
Test mode (code-path check on a partial search): AGG_TEST_OUT=<dir> SAMPLES="s1 s2" -> tolerates unfinished chunks,
writes everything under <dir> instead of W.
"""
import os
import sys
import time
import numpy as np
import pandas as pd

W = sys.argv[1]
S26 = "/scratch1/afreiburger/salterns_2026"
R, BT, RM = f"{S26}/rollout", f"{S26}/binning_test", f"{S26}/read_markers"
os.chdir(W)
t0 = time.time()
TEST = os.environ.get("AGG_TEST_OUT", "")
O = TEST if TEST else "."
samples = sorted(f[:-7] for f in os.listdir("meta/genes") if f.endswith(".tsv.gz"))
assert len(samples) == 39
if TEST:
    samples = os.environ["SAMPLES"].split()
    os.makedirs(f"{O}/results", exist_ok=True)
FAMILIES = {"uL2": ["K02886"], "uL11": ["K02867"], "uS3": ["K02982"], "uS9": ["K02996"], "uS11": ["K02948"],
            "uS12": ["K02950"], "uL4": ["K02926", "K02930"], "uL16": ["K02878", "K02866"]}
KO2FAM = {k: f for f, ks in FAMILIES.items() for k in ks}

# ---------------------------------------------------------------- representative calls
P2 = os.environ.get("P2DIR", "pass2")
calls = []
chunks = sorted(f for f in os.listdir(f"{P2}/results/chunks") if f.endswith(".calls.tsv.gz"))
nch = sum(1 for f in os.listdir(f"{P2}/chunks") if f.startswith("chunk_") and f.endswith(".faa"))
assert TEST or len(chunks) == nch, (len(chunks), nch)
for f in chunks:
    calls.append(pd.read_csv(f"{P2}/results/chunks/{f}", sep="\t", dtype=str, keep_default_na=False))
calls = pd.concat(calls, ignore_index=True).set_index("rep")
assert calls.index.is_unique
cand = pd.read_csv(f"{P2}/candidates.tsv.gz", sep="\t", dtype=str).set_index("rep")
assert TEST or set(calls.index) <= set(cand.index)
P148 = [l.strip() for l in open("meta/pass1_kos.txt") if l.strip()]
panel = pd.read_csv("meta/ko_panel_methylated.tsv", sep="\t", dtype=str)
assert len(P148) == 148 and set(panel.ko) | set(KO2FAM) == set(P148)
print(f"rep rows {len(calls)} ({time.time() - t0:.0f} s)", flush=True)
calls = calls[["best_ko", "best_score", "n_passing", "all_passing", "top_ko", "top_full_score", "rescue_all",
               "rescue_panel"]]

# ---------------------------------------------------------------- domains
mk = pd.read_csv(f"{R}/results/mags_kept.tsv", sep="\t", dtype=str)
mk["domain"] = mk.gtdb_classification.str.extract(r"d__(\w+)")[0]
assert mk.domain.isin(["Archaea", "Bacteria"]).all()
dom_mag = dict(zip(mk.mag_id, mk.domain))
bins = pd.concat([pd.read_csv(f"{R}/{s}/keep/bins.tsv", sep="\t", dtype=str) for s in samples], ignore_index=True)
kb = bins[bins.keep == "yes"].assign(domain=lambda d: d.mag_id.map(dom_mag))
assert kb.domain.notna().all()
lm = kb.groupby(["checkm1_lineage", "domain"]).size().unstack(fill_value=0)
for c in ("Archaea", "Bacteria"):
    if c not in lm:
        lm[c] = 0
lm["n"] = lm.Archaea + lm.Bacteria
# CheckM1 lineage (marker-set lineage) -> domain by name: archaeal taxa -> Archaea, root -> unknown, any other
# (bacterial) taxon -> Bacteria. Checked on the 935 kept bins against their GTDB R232 domain (column 'by_name').
ARCH = ("k__Archaea", "p__Euryarchaeota", "p__Crenarchaeota", "p__Thaumarchaeota", "c__Halobacteria",
        "c__Methanomicrobia", "c__Methanobacteria", "c__Thermoplasmata", "c__Thermococci", "o__Methanosarcinales",
        "o__Halobacteriales", "f__Halobacteriaceae", "f__Methanosarcinaceae")


def lineage_domain(x):
    if not isinstance(x, str) or x.startswith("root") or not x.strip():
        return "unknown"
    return "Archaea" if x.split(" ")[0] in ARCH else "Bacteria"


lm["by_name"] = [lineage_domain(x) for x in lm.index]
lm["agree"] = np.where(lm.by_name == "Archaea", lm.Archaea, np.where(lm.by_name == "Bacteria", lm.Bacteria, 0))
n_agree, n_kept = int(lm.agree.sum()), int(lm.n.sum())
ob = bins[bins.keep != "yes"].copy()
ob["domain"] = ob.checkm1_lineage.map(lineage_domain)
dom_otherbin = {(r.sample, r.bin): r.domain for r in ob.itertuples()}
lm.to_csv(f"{O}/results/checkm1_lineage_domain_map.tsv", sep="\t")
ob_lin = ob.groupby(["checkm1_lineage", "domain"]).size().rename("n_other_bins").reset_index()
ob_lin.to_csv(f"{O}/results/checkm1_lineage_domain_other_bins.tsv", sep="\t", index=False)
print(f"lineage->domain by name agrees with GTDB for {n_agree}/{n_kept} kept bins; other bins by domain",
      ob.domain.value_counts().to_dict(), flush=True)

# ---------------------------------------------------------------- per sample
status = pd.read_csv(f"{RM}/meta/contig_status.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
ge = pd.read_csv(f"{RM}/results/genome_equivalents.tsv", sep="\t", index_col=0)
os.makedirs(f"{O}/genes_ko", exist_ok=True)
agg, sd_rows, ge_rows, resc_rows, qc = [], [], [], [], []
for s in samples:
    g = pd.read_csv(f"meta/genes/{s}.tsv.gz", sep="\t", dtype={"partial": str})
    n_total = len(g)
    g = g.join(calls, on="rep")
    g["searched_full_db"] = g.rep.isin(cand.index).astype(int)
    for c in ["best_ko", "best_score", "n_passing", "all_passing", "top_ko", "top_full_score", "rescue_all",
              "rescue_panel"]:
        g[c] = g[c].fillna("")
    dp = f"{R}/{s}/depth/{s}.d1000.txt"
    if not os.path.exists(dp):
        dp = f"{BT}/{s}/depth/{s}.d1000.txt"
    dep = pd.read_csv(dp, sep="\t")
    assert list(dep.columns[:3]) == ["contigName", "contigLen", "totalAvgDepth"] and dep.shape[1] == 5, dp
    dep = dep.set_index("contigName").totalAvgDepth
    g["depth"] = g.contig.map(dep)
    assert g.depth.notna().all(), (s, g.depth.isna().sum())
    st = status[status["sample"] == s].set_index("contig")
    g["status"] = g.contig.map(st.status)
    assert g.status.notna().all(), s
    g["bin"] = g.contig.map(st.bin)
    g["mag_id"] = g.contig.map(st.mag_id)
    dom = np.where(g.status == "kept_mag", g.mag_id.map(dom_mag),
                   np.where(g.status == "other_bin",
                            [dom_otherbin.get((s, b), "unknown") for b in g.bin], ""))
    g["domain"] = dom
    assert not ((g.status == "kept_mag") & g.domain.isna()).any()
    g.insert(1, "sample", s)
    g[["gene", "sample", "contig", "partial", "aa_len", "rep", "searched_full_db", "best_ko", "best_score", "n_passing", "all_passing",
       "rescue_all", "rescue_panel", "top_ko", "top_full_score", "depth", "status", "mag_id", "domain"]].rename(
        columns={"best_ko": "ko", "best_score": "score", "rescue_all": "marker_rescue_allprofile",
                 "rescue_panel": "marker_rescue_v1"}).to_csv(
        f"{O}/genes_ko/{s}.tsv.gz.part", sep="\t", index=False, compression={"method": "gzip", "compresslevel": 4})
    os.replace(f"{O}/genes_ko/{s}.tsv.gz.part", f"{O}/genes_ko/{s}.tsv.gz")
    called_any = g[g.best_ko != ""]
    called = called_any[called_any.best_ko.isin(P148)]
    a = called.groupby("best_ko").agg(n_genes=("gene", "size"), sum_depth=("depth", "sum"))
    a["sample"] = s
    a["n_genes_total_sample"] = n_total
    agg.append(a.reset_index().rename(columns={"best_ko": "ko"}))
    cat = called.status + np.where(called.domain != "", "_" + called.domain, "")
    b = called.assign(cat=cat).groupby(["best_ko", "cat"]).agg(n_genes=("gene", "size"), sum_depth=("depth", "sum"))
    b = b.reset_index().rename(columns={"best_ko": "ko", "cat": "status_domain"})
    b.insert(0, "sample", s)
    sd_rows.append(b)
    # markers: strict best-KO + rescue; v1 reading (top over the 62 v1 panel KOs) for the GE, all-profile for comparison
    fam_strict = called[called.best_ko.isin(KO2FAM)].assign(fam=lambda d: d.best_ko.map(KO2FAM))
    r_all = g[(g.best_ko == "") & (g.rescue_all != "")].assign(fam=lambda d: d.rescue_all.map(KO2FAM))
    r_pan = g[(g.best_ko == "") & (g.rescue_panel != "")].assign(fam=lambda d: d.rescue_panel.map(KO2FAM))
    # GE uses v1's reading of the rescue (top-scoring profile among the 62 v1 panel KOs): with the top over all
    # 27,981 profiles, archaeal uS9 / uS11 / uL2 genes whose top profile is the eukaryote-named archaeal homolog KO
    # (K02960 S16e, K02955 S14e, K02938 L8e; scores below those KOs' thresholds) are lost, e.g. in R1_A_H2O
    # (Haloquadratum-dominated). The all-profile reading is kept as depth_allprofile_<fam> for comparison.
    row = {"sample": s}
    for fam in FAMILIES:
        ds = fam_strict.loc[fam_strict.fam == fam, "depth"].sum()
        da = r_all.loc[r_all.fam == fam, "depth"].sum()
        dpn = r_pan.loc[r_pan.fam == fam, "depth"].sum()
        row[f"depth_{fam}"] = ds + dpn
        row[f"depth_allprofile_{fam}"] = ds + da
        resc_rows.append(dict(sample=s, family=fam,
                              n_strict=int((fam_strict.fam == fam).sum()), depth_strict=ds,
                              n_rescue_all=int((r_all.fam == fam).sum()), depth_rescue_all=da,
                              n_rescue_v1panel=int((r_pan.fam == fam).sum()), depth_rescue_v1panel=dpn))
    ge_rows.append(row)
    qc.append(dict(sample=s, genes=n_total, genes_searched_full_db=int(g.searched_full_db.sum()),
                   genes_with_any_ko_among_searched=len(called_any), genes_with_panel_or_marker_ko=len(called),
                   genes_marker_rescue_only=len(r_pan), genes_marker_rescue_only_allprofile=len(r_all), sum_depth_all_genes=g.depth.sum(),
                   genes_kept_mag=int((g.status == "kept_mag").sum()),
                   genes_other_bin=int((g.status == "other_bin").sum()),
                   genes_unbinned=int((g.status == "unbinned").sum())))
    print(s, n_total, len(called), f"{time.time() - t0:.0f} s", flush=True)

gdf = pd.DataFrame(ge_rows).set_index("sample")
fam_cols = [f"depth_{f}" for f in FAMILIES]
gdf["GE_asm_mean"] = gdf[fam_cols].mean(axis=1)
gdf["GE_asm_median"] = gdf[fam_cols].median(axis=1)
gdf["GE_asm_cv"] = gdf[fam_cols].std(axis=1) / gdf[fam_cols].mean(axis=1)
fam_cols_all = [f"depth_allprofile_{f}" for f in FAMILIES]
gdf["GE_asm_mean_allprofile"] = gdf[fam_cols_all].mean(axis=1)
gdf["GE_asm_cv_allprofile"] = gdf[fam_cols_all].std(axis=1) / gdf[fam_cols_all].mean(axis=1)
gdf["GE_read_mean"] = ge.GE_depth_mean.reindex(gdf.index)
gdf["GE_read_cv"] = ge.GE_depth_cv.reindex(gdf.index)
gdf["GE_asm_over_read"] = gdf.GE_asm_mean / gdf.GE_read_mean
assert gdf.GE_read_mean.notna().all()
gdf.to_csv(f"{O}/results/ge_assembly.tsv", sep="\t", float_format="%.6g")
pd.DataFrame(resc_rows).to_csv(f"{O}/results/marker_rescue_check.tsv", sep="\t", index=False, float_format="%.6g")

A = pd.concat(agg, ignore_index=True)
kos = P148
full = pd.MultiIndex.from_product([kos, samples], names=["ko", "sample"])
K = A.set_index(["ko", "sample"])[["n_genes", "sum_depth"]].reindex(full, fill_value=0).reset_index()
ntot = pd.DataFrame(qc).set_index("sample").genes
K["copies_per_GE_assembly"] = K.sum_depth / K["sample"].map(gdf.GE_asm_mean)
K["copies_per_GE_readGE"] = K.sum_depth / K["sample"].map(gdf.GE_read_mean)
K["gene_fraction"] = K.n_genes / K["sample"].map(ntot)
info = pd.read_csv("db/hmm_info.tsv", sep="\t", dtype=str).set_index("ko")
K.insert(1, "definition", K.ko.map(info.definition))
gene_of = dict(zip(panel.ko, panel.gene))
gene_of.update({k: f for k, f in KO2FAM.items()})
K["gene"] = K.ko.map(gene_of)
K["role"] = np.where(K.ko.isin(KO2FAM), "marker", "panel")
K["kofam_threshold"] = K.ko.map(info.threshold)
K = K[["sample", "ko", "gene", "role", "kofam_threshold", "definition", "n_genes", "sum_depth", "copies_per_GE_assembly", "copies_per_GE_readGE",
       "gene_fraction"]]
K.to_csv(f"{O}/results/ko_by_sample.tsv.gz", sep="\t", index=False, float_format="%.6g")
pd.concat(sd_rows, ignore_index=True).to_csv(f"{O}/results/ko_status_domain_by_sample.tsv.gz", sep="\t", index=False,
                                             float_format="%.6g")
q = pd.DataFrame(qc)
q.to_csv(f"{O}/results/aggregate_qc.tsv", sep="\t", index=False, float_format="%.6g")
with open(f"{O}/results/aggregate_qc.txt", "w") as fo:
    fo.write(f"samples\t{len(samples)}\ngenes\t{q.genes.sum()}\n"
             f"genes_searched_full_db\t{q.genes_searched_full_db.sum()}\n"
             f"genes_with_any_ko_among_searched\t{q.genes_with_any_ko_among_searched.sum()}\n"
             f"genes_with_panel_or_marker_ko\t{q.genes_with_panel_or_marker_ko.sum()}\t"
             f"fraction\t{q.genes_with_panel_or_marker_ko.sum() / q.genes.sum():.4f}\n"
             f"genes_marker_rescue_only\t{q.genes_marker_rescue_only.sum()}\n"
             f"panel_marker_kos_called_in_any_sample\t{int((K.groupby('ko').n_genes.sum() > 0).sum())}\t"
             f"of\t{len(kos)}\n"
             f"checkm1_lineage_domain_agrees_with_gtdb_kept_bins\t{n_agree}\t{n_kept}\n"
             f"GE_asm_cv_range\t{gdf.GE_asm_cv.min():.4f}\t{gdf.GE_asm_cv.max():.4f}\t"
             f"median\t{gdf.GE_asm_cv.median():.4f}\n"
             f"GE_asm_cv_range_allprofile_rescue\t{gdf.GE_asm_cv_allprofile.min():.4f}\t"
             f"{gdf.GE_asm_cv_allprofile.max():.4f}\tmedian\t{gdf.GE_asm_cv_allprofile.median():.4f}\n"
             f"genes_marker_rescue_only_allprofile\t{q.genes_marker_rescue_only_allprofile.sum()}\nseconds\t{time.time() - t0:.0f}\n")
print(open(f"{O}/results/aggregate_qc.txt").read())
