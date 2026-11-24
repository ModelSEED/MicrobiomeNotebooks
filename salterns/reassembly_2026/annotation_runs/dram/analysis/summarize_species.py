"""Per-species summary of the DRAM 1.5 run on the 558 best-per-species MAGs, for the manuscript.

Run from annotation_runs/dram/ with ~/Documents/py_venv/bin/python -I analysis/summarize_species.py
(after analysis/validate_vs_kofam.py, which shows DRAM's per-gene KO calls equal the KOfam track's).

Inputs
  results/merged/annotations.tsv.gz        DRAM gene table (source of truth for counts)
  results/distill/product.tsv              DRAM product table (module completeness / function presence)
  results/distill/genome_stats.tsv         tRNA / rRNA / quality
  results/distill_forms/function_heatmap_form.20261111.tsv   DRAM's product function definitions
  ../annotation_choice/custom_distillate_salterns.dram1.tsv  the 'Salterns methane' sheet definition (65 KOs)
  ../kofam/results/guilds_1324.tsv, guilds_species_clusters_new.tsv  KOfam guilds incl. the Sec read-through rescue
  ../kofam/results/per_mag/<MAG>.kofam.tsv.gz  per-gene KOfam scores (for the scores behind product calls)
  best_per_species_558.tsv                 species cluster of each MAG

Outputs (results/)
  species_558_dram_summary.tsv             one row per species cluster (wide, tidy: one value per cell)
  species_558_salterns_methane_long.tsv    long: species x 65 sheet KOs (count, gene ids)
  summary/salterns_methane_by_ko.tsv       per KO: species, genes, top families
  summary/product_functions_of_interest.tsv  per product function: species lit, contributing KOs, top families
  summary/guilds_dram_vs_kofam_558.tsv     guild totals: DRAM (no rescue) vs KOfam (with Sec rescue), best MAG
  summary/notable_product_calls.tsv        the methylotrophy / monooxygenase product calls listed per species
  summary/k14083_mttB_by_family.tsv
"""
import os
import re

import warnings

import numpy as np
import pandas as pd

warnings.simplefilter("ignore", pd.errors.PerformanceWarning)

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = os.path.dirname(HERE)
KOF = os.path.join(RUNS, "kofam")
RES = os.path.join(HERE, "results")
SUM = os.path.join(RES, "summary")
os.makedirs(SUM, exist_ok=True)


def sym_col(sym):
    return re.sub(r"[^A-Za-z0-9]+", "_", sym).strip("_")


best = pd.read_csv(os.path.join(HERE, "best_per_species_558.tsv"), sep="\t").set_index("name")
MAGS = list(best.index)

sheet_def = pd.read_csv(os.path.join(RUNS, "..", "annotation_choice", "custom_distillate_salterns.dram1.tsv"), sep="\t")
assert len(sheet_def) == 65
KO2SYM = dict(zip(sheet_def.gene_id, sheet_def.subheader))
KO2HDR = dict(zip(sheet_def.gene_id, sheet_def.header))

print("loading annotations ...", flush=True)
ann = pd.read_csv(os.path.join(RES, "merged", "annotations.tsv.gz"), sep="\t", index_col=0, low_memory=False,
                  usecols=lambda c: c in ("Unnamed: 0", "fasta", "scaffold", "ko_id", "kegg_hit", "pfam_hits",
                                          "cazy_best_hit", "peptidase_family") or c == "")
ann.index.name = "gene"
assert ann.fasta.nunique() == 558

# ------------------------------------------------------------------ taxonomy & quality
tax = pd.read_csv(os.path.join(HERE, "gtdb_taxonomy_558.tsv"), sep="\t").set_index("user_genome").loc[MAGS]
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
taxdf = tax.classification.str.split(";", expand=True)
taxdf.columns = RANKS
taxdf = taxdf.apply(lambda s: s.str.replace(r"^[dpcofgs]__", "", regex=True))
gs = pd.read_csv(os.path.join(RES, "distill", "genome_stats.tsv"), sep="\t").set_index("genome").loc[MAGS]

# ------------------------------------------------------------------ per-MAG gene statistics
g = ann.groupby("fasta")
stats = pd.DataFrame({
    "n_genes": g.size(),
    "n_genes_with_ko": g.ko_id.count(),
    "n_genes_with_pfam": g.pfam_hits.count(),
    "n_genes_with_cazy": g.cazy_best_hit.count(),
    "n_genes_with_merops": g.peptidase_family.count(),
}).loc[MAGS]
stats["frac_genes_with_ko"] = (stats.n_genes_with_ko / stats.n_genes).round(4)

# ------------------------------------------------------------------ 'Salterns methane' KO counts from annotations
sub = ann[ann.ko_id.isin(KO2SYM)]
cnt = sub.groupby(["fasta", "ko_id"]).size().unstack(fill_value=0).reindex(index=MAGS, columns=list(KO2SYM),
                                                                          fill_value=0)
# the sheet itself must agree (validated in validate_vs_kofam.py; re-asserted here)
sheet = pd.read_excel(os.path.join(RES, "distill", "metabolism_summary.xlsx"), sheet_name="Salterns methane") \
          .set_index("gene_id")[MAGS].astype(int)
assert (sheet.T.reindex(columns=list(KO2SYM)) == cnt).all().all()
C = {KO2SYM[k]: cnt[k] for k in KO2SYM}  # symbol -> Series over MAGS

# long table with gene ids
gene_ids = sub.reset_index().groupby(["fasta", "ko_id"]).gene.apply(lambda s: ";".join(sorted(s)))
long = cnt.stack().rename("n_genes").reset_index()
long.columns = ["mag", "ko", "n_genes"]
long["gene"] = long.ko.map(KO2SYM)
long["header"] = long.ko.map(KO2HDR)
long["species_cluster"] = long.mag.map(best.species_cluster)
long["gtdb_species_or_genus"] = long.mag.map(taxdf.species.where(taxdf.species != "", taxdf.genus + " sp."))
long["dram_gene_ids"] = [gene_ids.get((m, k), "") for m, k in zip(long.mag, long.ko)]
long = long[["species_cluster", "mag", "gtdb_species_or_genus", "header", "gene", "ko", "n_genes", "dram_gene_ids"]]
long.to_csv(os.path.join(RES, "species_558_salterns_methane_long.tsv"), sep="\t", index=False)

# ------------------------------------------------------------------ DRAM guild calls (KOfam-track rules, no Sec rescue)
phn7 = ["phnG", "phnH", "phnI", "phnJ", "phnK", "phnL", "phnM"]
mt = C["mttB"] + C["mtbB"] + C["mtmB"]
cut = (C["cutC"] > 0) & (C["cutD"] > 0)
cnt_ = (C["cntA/yeaW"] > 0) & (C["cntB/yeaX"] > 0)
dg = pd.DataFrame(index=MAGS)
dg["dram_methanogen"] = (C["mcrA"] > 0).astype(int)
dg["dram_mepn_consumer"] = (C["phnJ"] > 0).astype(int)
dg["dram_phnGHIJKLM_n_distinct"] = sum((C[s] > 0).astype(int) for s in phn7)
dg["dram_betaine_reducer"] = ((C["grdH"] > 0) | (C["grdI"] > 0)).astype(int)
dg["dram_methylamine_user"] = ((mt > 0) & (C["mcrA"] > 0)).astype(int)
dg["dram_mttB_family_without_mcrA"] = (((mt + C["mtgB"]) > 0) & (C["mcrA"] == 0)).astype(int)
dg["dram_tma_producer"] = (cut | cnt_).astype(int)
dg["dram_tma_candidate_one_sided"] = ((~(cut | cnt_)) & ((C["cutC"] + C["cntA/yeaW"]) > 0)).astype(int)
dg["dram_pylBCDS_n_distinct"] = sum((C[s] > 0).astype(int) for s in ["pylB", "pylC", "pylD", "pylS"])

# ------------------------------------------------------------------ KOfam guilds (with Sec rescue) for the same MAG
g1324 = pd.read_csv(os.path.join(KOF, "results", "guilds_1324.tsv"), sep="\t").set_index("mag").loc[MAGS]
assert (g1324.species_cluster == best.species_cluster).all()
kcols = ["methanogen", "mepn_consumer", "phnGHIJKLM_n_distinct", "phnGHIJKLM_n_contigs", "betaine_reducer",
         "betaine_reducer_evidence", "betaine_reducer_without_sec_rescue", "n_grdH", "n_grdH_sec_readthrough",
         "n_grdH_prodigal_only", "n_grdI", "n_grdB", "n_grdA", "methylamine_user", "mttB_family_without_mcrA",
         "tma_producer", "tma_candidate_one_sided", "methanogen_short_maxbin2_only", "mepn_consumer_short_maxbin2_only",
         "betaine_reducer_short_maxbin2_only", "methylamine_user_short_maxbin2_only",
         "marker_genes_on_short_maxbin2_contigs"]
kg = g1324[kcols].add_prefix("kofam_")
kg["kofam_betaine_reducer_evidence"] = kg.kofam_betaine_reducer_evidence.fillna("")
kg["kofam_marker_genes_on_short_maxbin2_contigs"] = kg.kofam_marker_genes_on_short_maxbin2_contigs.fillna("")

gsc = pd.read_csv(os.path.join(KOF, "results", "guilds_species_clusters_new.tsv"), sep="\t").set_index("species_cluster")
gsc = gsc.loc[best.species_cluster]
assert (gsc.best_mag.values == np.array(MAGS)).all(), "KOfam best_mag != DRAM MAG"
ccols = ["n_mags"] + [f"{x}_{y}" for x in ["methanogen", "mepn_consumer", "betaine_reducer", "methylamine_user",
                                            "tma_candidate_one_sided", "mttB_family_without_mcrA"]
                      for y in ["any_mag", "n_mags"]]
cl = gsc[ccols].add_prefix("cluster_")
cl.index = MAGS

# ------------------------------------------------------------------ product table, functions of interest
prod = pd.read_csv(os.path.join(RES, "distill", "product.tsv"), sep="\t").set_index("genome").loc[MAGS]
C1 = ["Methanogenesis, CO2 => methane", "Acetyl-CoA pathway, CO2 => acetyl-CoA",
      "Reductive acetyl-CoA pathway (Wood-Ljungdahl pathway)", "Other Reductases: TMAO reductase"]
POI = [c for c in prod.columns if c.startswith(("Methanogenesis and methanotrophy:", "Nitrogen metabolism:",
                                                "Sulfur metabolism:", "CAZy:"))]
POI = C1 + POI
prod_poi = prod[POI].astype(float)
prod_poi.columns = ["product: " + c for c in POI]

# reproduce DRAM's non-CAZy function calls from ko_id/EC to get the contributing ids (DRAM 1.5 make_functional_df:
# a row is present if any of its ids is present, a function if all its rows are)
fh = pd.read_csv(os.path.join(RES, "distill_forms", "function_heatmap_form.20261111.tsv"), sep="\t").fillna("")
fh = fh.apply(lambda s: s.str.strip() if s.dtype == object else s)
ec = ann.kegg_hit.dropna().str.findall(r"\[EC:([^\]]+)\]").explode().dropna().str.split().explode()
ids = pd.concat([ann.ko_id.dropna().rename("id").to_frame().assign(fasta=ann.fasta),
                 ("EC:" + ec).rename("id").to_frame().join(ann.fasta)])
ids_by = ids.groupby("fasta").id.apply(set).reindex(MAGS)
contrib = {}
repro_mismatch = {}
for fname, frame in fh.groupby("function_name", sort=False):
    cat = frame.category.iloc[0]
    col = f"{cat}: {fname}"
    if col not in prod.columns or cat == "CAZy":
        continue
    rows = [set(i.strip() for i in r.split(",")) for r in frame.function_ids]
    pres, present_ids = [], []
    for m in MAGS:
        s = ids_by[m]
        hits = [s & r for r in rows]
        pres.append(all(len(h) > 0 for h in hits))
        present_ids.append(",".join(sorted(set().union(*hits))))
    pres = pd.Series(pres, index=MAGS)
    repro_mismatch[col] = int((pres != prod[col].astype(bool)).sum())
    contrib[col] = pd.Series(present_ids, index=MAGS)
assert sum(repro_mismatch.values()) == 0, {k: v for k, v in repro_mismatch.items() if v}

# ------------------------------------------------------------------ assemble the per-species table
out = pd.DataFrame(index=MAGS)
out.index.name = "mag"
out["species_cluster"] = best.species_cluster
out["sample"] = best["sample"]
out["binner"] = [m.split(".")[0] for m in MAGS]
for r in RANKS:
    out[r] = taxdf[r]
out["completeness"] = best.completeness.round(2)
out["contamination"] = best.contamination.round(2)
out["size_bp"] = best.size_bp
out["n_contigs"] = best.n_contigs
out["n50"] = best.n50
out["tRNA_count"] = gs["tRNA count"]
for r in ["5S rRNA", "16S rRNA", "23S rRNA"]:
    out[r.replace(" ", "_") + "_present"] = gs[r].notna().astype(int)
out["dram_assembly_quality"] = gs["assembly quality"]
out = out.join(stats)
for ko, sym in KO2SYM.items():
    out[f"n_{sym_col(sym)}"] = cnt[ko]
out = out.join(dg).join(kg).join(cl)
# manuscript calls: DRAM KO layer + the KOfam Sec read-through rescue for grdH
out["betaine_reducer_final"] = out.kofam_betaine_reducer
out["grdH_final"] = out.kofam_n_grdH
assert (out.dram_betaine_reducer == out.kofam_betaine_reducer_without_sec_rescue).all()
assert (out.dram_betaine_reducer <= out.kofam_betaine_reducer).all()
for gname in ["methanogen", "mepn_consumer", "methylamine_user", "mttB_family_without_mcrA", "tma_producer",
              "tma_candidate_one_sided"]:
    assert (out[f"dram_{gname}"] == out[f"kofam_{gname}"]).all(), gname
out["product_TMA_to_DMA_lit"] = prod["Methanogenesis and methanotrophy: trimethylamine => dimethylamine"].astype(int)
out["k14083_caveat"] = np.where((out.product_TMA_to_DMA_lit == 1) & (out.n_mcrA == 0),
                                "K14083 (mttB-family) without mcrA: not methylamine methanogenesis", "")
out = out.join(prod_poi)
for col, s in contrib.items():
    if col.startswith(("Methanogenesis and methanotrophy:", "Nitrogen metabolism:", "Sulfur metabolism:",
                       "Other Reductases: TMAO")):
        out["ids: " + col] = s
out.reset_index().set_index("species_cluster").reset_index().to_csv(
    os.path.join(RES, "species_558_dram_summary.tsv"), sep="\t", index=False)

# slim manuscript table: one row per species, the methane-cycle columns only
bt = ["proV", "proW", "proX", "opuA", "opuBD", "opuC", "betT", "opuD", "BCCT"]
core = pd.DataFrame(index=out.index)
for c in ["species_cluster", "domain", "phylum", "class", "order", "family", "genus", "species", "completeness",
          "contamination", "n_genes", "cluster_n_mags"]:
    core[c] = out[c]
core["n_phnJ"] = out.n_phnJ
core["phnGHIJKLM_n_distinct"] = out.dram_phnGHIJKLM_n_distinct
core["phnGHIJKLM_n_contigs"] = out.kofam_phnGHIJKLM_n_contigs
core["phnCDE_n_distinct"] = sum((C[s] > 0).astype(int) for s in ["phnC", "phnD", "phnE"])
core["mepn_consumer"] = out.dram_mepn_consumer
core["mepn_consumer_needs_short_maxbin2_contig"] = out.kofam_mepn_consumer_short_maxbin2_only
core["cluster_mepn_consumer_n_mags"] = out.cluster_mepn_consumer_n_mags
core["n_grdI"] = out.n_grdI
core["n_grdH_dram"] = out.n_grdH
core["n_grdH_with_sec_rescue"] = out.kofam_n_grdH
core["n_grdH_sec_readthrough"] = out.kofam_n_grdH_sec_readthrough
core["grdACD_n_distinct"] = sum((C[s] > 0).astype(int) for s in ["grdA", "grdC", "grdD"])
core["betaine_reducer_dram"] = out.dram_betaine_reducer
core["betaine_reducer_final"] = out.betaine_reducer_final
core["betaine_reducer_evidence"] = out.kofam_betaine_reducer_evidence
core["cluster_betaine_reducer_n_mags"] = out.cluster_betaine_reducer_n_mags
core["n_mtgB"] = out.n_mtgB
core["betaine_transport_n_distinct"] = sum((C[s] > 0).astype(int) for s in bt)
for s_ in ["cutC", "cutD", "cntA/yeaW", "cntB/yeaX", "torA", "torZ", "dmsA", "tmm", "dmd-tmd", "tdm"]:
    core[f"n_{sym_col(s_)}"] = C[s_]
core["tma_producer"] = out.dram_tma_producer
core["tma_candidate_one_sided"] = out.dram_tma_candidate_one_sided
for s_ in ["mcrA", "mcrB", "mcrG", "mttB", "mtbB", "mtmB", "mtbA", "mtaB", "mtaC", "mtaA"]:
    core[f"n_{s_}"] = C[s_]
core["pylBCDS_n_distinct"] = out.dram_pylBCDS_n_distinct
core["methanogen"] = out.dram_methanogen
core["methylamine_user"] = out.dram_methylamine_user
core["mttB_family_without_mcrA"] = out.dram_mttB_family_without_mcrA
core["product_TMA_to_DMA_lit"] = out.product_TMA_to_DMA_lit
core["k14083_caveat"] = out.k14083_caveat
core["product_CO2_to_methane_module"] = prod["Methanogenesis, CO2 => methane"].astype(float).round(3)
for short, col in [("product_methanol_to_methane", "Methanogenesis and methanotrophy: methanol => methane"),
                   ("product_DMA_to_MMA", "Methanogenesis and methanotrophy: dimethylamine => monomethylamine"),
                   ("product_MMA_to_ammonia", "Methanogenesis and methanotrophy: monomethylamine => ammonia"),
                   ("product_dsrAB", "Sulfur metabolism: dissimilatory sulfate reduction (and oxidation) sulfate => sulfide"),
                   ("product_sox", "Sulfur metabolism: Thiosulfate oxidation by SOX complex, thiosulfate => sulfate"),
                   ("product_amo", "Nitrogen metabolism: Bacterial/Archaeal ammonia oxidation"),
                   ("product_nifHD", "Nitrogen metabolism: nitrogen => ammonia"),
                   ("product_nosZ", "Nitrogen metabolism: nitrous oxide => nitrogen")]:
    core[short] = prod[col].astype(int)
core.reset_index().set_index("species_cluster").reset_index().to_csv(
    os.path.join(RES, "species_558_methane_guilds.tsv"), sep="\t", index=False)

# ------------------------------------------------------------------ summaries
def top_taxa(mask, rank="family", n=6):
    vc = out.loc[mask, rank].replace("", "(unassigned)").value_counts()
    return "; ".join(f"{k} {v}" for k, v in vc.head(n).items()) + ("; ..." if len(vc) > n else "")


rows = []
for ko, sym in KO2SYM.items():
    m = cnt[ko] > 0
    rows.append(dict(header=KO2HDR[ko], gene=sym, ko=ko, species=int(m.sum()), genes=int(cnt[ko].sum()),
                     archaea=int((m & (out.domain == "Archaea")).sum()),
                     top_families=top_taxa(m) if m.any() else ""))
pd.DataFrame(rows).to_csv(os.path.join(SUM, "salterns_methane_by_ko.tsv"), sep="\t", index=False)

rows = []
for c in POI:
    v = prod[c].astype(float)
    lit = v > 0
    r = dict(function=c, species_nonzero=int(lit.sum()), species_complete=int((v >= 1).sum()),
             top_families_nonzero=top_taxa(lit) if lit.any() else "")
    if c in contrib:
        ids_ = contrib[c][prod[c].astype(bool)].str.split(",").explode()
        r["contributing_ids_in_lit_species"] = "; ".join(f"{k} {v}" for k, v in ids_.value_counts().items())
    rows.append(r)
pd.DataFrame(rows).to_csv(os.path.join(SUM, "product_functions_of_interest.tsv"), sep="\t", index=False)

# guild totals
GU = ["methanogen", "mepn_consumer", "betaine_reducer", "methylamine_user", "tma_producer", "tma_candidate_one_sided",
      "mttB_family_without_mcrA"]
rows = []
for gname in GU:
    rows.append(dict(guild=gname, dram_best_mag=int(out[f"dram_{gname}"].sum()),
                     kofam_best_mag=int(out[f"kofam_{gname}"].sum()),
                     kofam_any_mag_in_cluster=int(out[f"cluster_{gname}_any_mag"].sum())
                     if f"cluster_{gname}_any_mag" in out else None))
pd.DataFrame(rows).to_csv(os.path.join(SUM, "guilds_dram_vs_kofam_558.tsv"), sep="\t", index=False)

# per-gene KOfam scores behind the notable product calls
NOTE = {"Methanogenesis and methanotrophy: Key functional gene": ["K00399", "K00401", "K00402"],
        "Methanogenesis and methanotrophy: dimethylamine => monomethylamine": ["K16178"],
        "Methanogenesis and methanotrophy: monomethylamine => ammonia": ["K16176"],
        "Methanogenesis and methanotrophy: methanol => methane": ["K04480"],
        "Methanogenesis and methanotrophy: methane => methanol, with oxygen (pmo)": ["K10944", "K10945", "K10946"],
        "Nitrogen metabolism: ammonia => nitrite": ["K10944", "K10945", "K10946"],
        "Methanogenesis and methanotrophy: putative but not defining CO2 => methane": ["K00200"],
        "Nitrogen metabolism: Bacterial (aerobic-specific) ammonia oxidation": ["K10535"],
        "Methanogenesis, CO2 => methane": []}
need = set(k for v in NOTE.values() for k in v)
rows = []
for c, kos in NOTE.items():
    v = prod[c].astype(float)
    sel = v >= 1 if c == "Methanogenesis, CO2 => methane" else v > 0
    for m in v[sel].index:
        if kos:
            pm = pd.read_csv(os.path.join(KOF, "results", "per_mag", f"{m}.kofam.tsv.gz"), sep="\t", dtype=str)
            pm = pm[pm.kofam_id.isin(kos)]
            detail = "; ".join(f"{k} {s} (thr {t}, {a} aa, contig {cl_} bp)" for k, s, t, a, cl_ in
                               zip(pm.kofam_id, pm.kofam_score, pm.kofam_threshold, pm.aa_len, pm.contig_len))
        else:
            detail = "module completeness %.3f" % v[m]
        rows.append(dict(function=c, mag=m, species_cluster=out.at[m, "species_cluster"], phylum=out.at[m, "phylum"],
                         family=out.at[m, "family"], genus=out.at[m, "genus"], species=out.at[m, "species"],
                         completeness=out.at[m, "completeness"], n_mcrA=out.at[m, "n_mcrA"],
                         kofam_genes_scores=detail))
pd.DataFrame(rows).to_csv(os.path.join(SUM, "notable_product_calls.tsv"), sep="\t", index=False)

# K14083 caveat by family
m14 = out.product_TMA_to_DMA_lit == 1
k14 = out[m14].assign(family=lambda d: d.family.replace("", "(unassigned)")).groupby(["phylum", "family"]).agg(
    species=("n_mttB", "size"), mttB_genes=("n_mttB", "sum"), with_mcrA=("n_mcrA", lambda s: int((s > 0).sum())),
    with_pyl=("dram_pylBCDS_n_distinct", lambda s: int((s > 0).sum()))).sort_values("species", ascending=False)
k14.to_csv(os.path.join(SUM, "k14083_mttB_by_family.tsv"), sep="\t")

print("species table:", out.shape)
print("product reproduction mismatches:", sum(repro_mismatch.values()))
print(pd.read_csv(os.path.join(SUM, "guilds_dram_vs_kofam_558.tsv"), sep="\t").to_string(index=False))
