#!/usr/bin/env python
"""Metagenome sections on the 2026 data (2026-11-26): writes r232_review_edits_2026-11-26.json from
r232_review_edits_2026-11-25.json (run make_synthesis_edits_2026-11-25.py first). As tracked changes:
  - Methods P28 (metagenome processing points to the re-assembly) and P30 (the 138-KO panel, Table S1);
  - Methods P50 (statistics: Kruskal-Wallis + BH over the panel, Dunn, Ward; Figure S4 on methanogen-lineage reads;
    Mantel on both community definitions);
  - Results P57 (Figure 1 on the read survey), P62 (Figure S4 and the Mantel test), P67 (Figure S6, isolated genera);
  - Discussion P81 (Bray-Curtis ranges), P93 (cultivated genera, Figure S6);
  - captions: Figure 1 (P115), Figure S4 (P140), Figure S6 (P144, as in the Google Doc since 2026-11-25), Table S1 (P169),
    Table S4 (P175);
  - "figures_2026_11_26": the Figure 1, S4 and S6 images that make_r232_review_docx.py places.
The read-survey Methods paragraph, the new assembly-based paragraph and the Data availability folders are in
make_r232_review_docx.py. Numbers come from reassembly_2026/community_ko/results/ (panel_stats.py, figure1.py,
figureS4.py, figureS6.py, metabolome_taxonomy.py); the wording's preconditions are asserted.
Run: ~/Documents/py_venv/bin/python make_metagenome_edits_2026-11-26.py
"""
import copy
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
S = HERE.parent
CK = S / "reassembly_2026/community_ko"
RES = CK / "results"
sys.path.insert(0, str(CK))
import profiles as PR  # noqa: E402
import stats_lib as SL  # noqa: E402

DATE = "2026-11-26T00:00:00Z"
SRC = ("reassembly_2026/community_ko/ (panel_stats.py, figure1.py, figureS4.py, figureS6.py, metabolome_taxonomy.py, "
       "table_s1.py); reassembly_2026/read_markers_v2/README.md; reassembly_2026/community_ko/assembly_based/README.md")

base = json.load(open(HERE / "r232_review_edits_2026-11-25.json"))
P = PR.panel().set_index("ko")
st = pd.read_csv(RES / "panel_stats_read.tsv", sep="\t").set_index("ko")
st["gene"] = P.gene
f1 = json.load(open(RES / "figure1_meta.json"))
s4 = json.load(open(RES / "figureS4_meta.json"))
rho = pd.read_csv(RES / "figureS4_rho.tsv", sep="\t")
s6 = pd.read_csv(RES / "figureS6_stats.tsv", sep="\t")
s6s = pd.read_csv(RES / "figureS6_read_sensitivity.tsv", sep="\t").set_index("genus")
mt = json.load(open(RES / "metabolome_taxonomy.json"))
g = {r.gene: r for _, r in st.iterrows()}


def m(gene, ty):                                   # mean copies per 100 GE
    return g[gene][f"mean_per100GE_{ty}"]


def letters(gene):
    r = g[gene]
    return "/".join(str(r[f"letter_{t}"]) for t in ("unrestored", "restored", "reference")) if r.q < 0.05 else "ns"


def f(x, nd=1):
    return f"{x:.{nd}f}"


# ---- preconditions ---------------------------------------------------------------------------------------------------
n_fig1, n_sig = f1["n_detected"], f1["n_q05"]
assert f1["n_panel"] == 108 and n_fig1 == 107 and f1["not_detected"] == ["K17066"], f1["not_detected"]
assert f1["ward_k3_matches_restoration"], "Ward k = 3 no longer matches restoration; rewrite the clustering sentence"
for x in ("phnG", "phnH", "phnI", "phnJ", "phnK", "phnL", "phnM", "grdI", "mcrA", "mcrG", "dmmA", "dmmC"):
    assert letters(x) == "a/b/b", (x, letters(x))
for x in ("mgdA", "mgdB", "mgdC", "mgdD", "dmdA", "dddL", "MTO", "dmmB"):
    assert letters(x) == "a/b/a" and m(x, "unrestored") > m(x, "reference"), (x, letters(x))
for x in ("tmd", "mauA", "mauB", "gmaS", "gbcA", "gbcB", "CMO", "fdh"):
    assert letters(x).startswith("b/a/a") and m(x, "reference") == max(m(x, t) for t in ("unrestored", "restored", "reference")), x
for x in ("cntA", "mxaF", "xoxF", "gfa", "frmA", "dgcA", "dmdB", "dmdD"):
    assert letters(x) == "b/b/a", (x, letters(x))
for x in ("cutC", "mttB", "mtbB", "mtmB", "mtaB"):
    assert letters(x) == "b/a/b", (x, letters(x))
assert letters("bhmt") == "b/a/a" and m("bhmt", "restored") > max(m("bhmt", "unrestored"), m("bhmt", "reference"))
# methanogen-lineage share of the methyltransferases by site type (24 sediment samples)
t = PR.read_profile()
t = t[t["sample"].isin(PR.SED)].copy()
t["type"] = t["sample"].map(SL.site_type)
sh = t.groupby(["ko", "type"])[["reads", "reads_methanogen_lineage"]].sum()
sh = (sh.reads_methanogen_lineage / sh.reads).unstack()
sh.index = [P.loc[k, "gene"] for k in sh.index]
mtx = ["mtbB", "mtmB", "mtaB"]
res_lo, res_hi = sh.loc[mtx, "restored"].min(), sh.loc[mtx, "restored"].max()
unr_lo, unr_hi = sh.loc[mtx, "unrestored"].min(), sh.loc[mtx, "unrestored"].max()
assert res_hi < 0.10 and unr_lo > 0.2, (res_lo, res_hi, unr_lo, unr_hi)
# Figure S4
def sig(gene, sign):
    d = rho[(rho.gene == gene) & (rho.q < 0.05)]
    return set(d[d.rho > 0].metabolite) if sign > 0 else set(d[d.rho < 0].metabolite)
for gname in ("mcrA", "mcrB", "mcrG"):
    assert {"Methylamine", "Acetate", "Trimethylamine"} <= sig(gname, +1), gname
    assert {"Fructose", "Proline", "2-Oxoisocaproate", "Sucrose"} <= sig(gname, -1), gname
assert "Trehalose" in sig("mcrB", +1) and "Trehalose" not in sig("mcrA", +1) | sig("mcrG", +1)
assert "Glucose" in sig("mcrB", -1) and "Glucose" in sig("mcrG", -1) and "Glucose" not in sig("mcrA", -1)
mA = rho[(rho.gene == "mcrA") & (rho.q < 0.05)]
pos_lo, pos_hi = mA[mA.rho > 0].rho.min(), mA[mA.rho > 0].rho.max()
neg_lo, neg_hi = mA[mA.rho < 0].rho.max(), mA[mA.rho < 0].rho.min()
assert {"Trimethylamine", "Trehalose"} <= sig("phnJ", +1) and {"Proline", "Sucrose", "Glucose", "Fructose"} <= sig("phnJ", -1)
assert s4["n_samples"] == 19 and s4["n_metabolites"] == 40 and s4["n_tests"] == 1960, s4
assert sorted(x.split(" (")[0] for x in s4["missing_genes"]) == ["eutD", "pta", "pta_K13788"], s4["missing_genes"]
# Mantel and Bray-Curtis
mg, rg = mt["MAG species (585)"], mt["read marker genera"]
assert round(mg["mantel_r"], 2) == 0.78 and round(rg["mantel_r"], 2) == 0.70 and mg["p"] < 0.001 and rg["p"] < 0.001
bc_m = mg["metab_bc_where_taxo_ge_083"]; bc_t = mg["taxo_bc_ge_083"]
# Figure S6
def L6(row, genus):
    r = s6[(s6.row == row) & (s6.genus == genus)].iloc[0]
    return r, (r.letter_R1, r.letter_R2, r.letter_SF2, r.letter_R2A)
for gen in ("Halomonas", "Salegentibacter"):
    for row in ("Metagenome (reads)", "MAGs"):
        r, l = L6(row, gen)
        assert r.p < 0.05 and set(l[0]) & set(l[1]) and not (set(l[0]) & set(l[2])) and not (set(l[0]) & set(l[3])), (gen, row, l)
r, l = L6("Metagenome (reads)", "Idiomarina"); assert r.p > 0.05
r, l = L6("Metagenome (reads)", "Marinilabilia"); assert r.p < 0.05 and l == ("a", "a", "a", "b"), l
r, l = L6("Metagenome (reads)", "Clostridium"); assert r.p < 0.05 and l[0] == l[1] == "b" and l[2] == l[3] == "a", l
assert s6s.loc["Clostridium", "frac_best_hit_ge95"] < 0.1

FACTS = dict(fig1=dict(n_detected=n_fig1, n_q05=n_sig), methyltransferase_methanogen_share=dict(
    restored=[res_lo, res_hi], unrestored=[unr_lo, unr_hi]), mcrA_rho_pos=[pos_lo, pos_hi], mcrA_rho_neg=[neg_lo, neg_hi],
    mantel=dict(mag=mg["mantel_r"], reads=rg["mantel_r"]), bc=dict(metab=bc_m, taxo=bc_t))
NEW = {}


def put(key, paragraph, kind, old, new, rationale, **kw):
    NEW[key] = dict({"paragraph": paragraph, "kind": kind, "old": old, "new": new, "rationale": rationale,
                     "sources": SRC, "origin": "metagenome_2026_11_26", "date": DATE}, **kw)


def pct(x):
    return f"{x * 100:.0f}"


# ---- Methods -----------------------------------------------------------------------------------------------------------
put("p28_metagenomes", 28, "replace", "amplicons and metagenomes were processed as detailed by",
    "amplicons were processed as detailed by", "the metagenomes are re-processed (re-assembly, Methods below)")
put("p28_metagenomes_2", 28, "insert_before", "Sediments were sent to the UC Davis",
    "The 39 metagenomes, 26 of sediment and 13 of water (Table S4), were re-assembled from the raw reads as described "
    "below. ", "points to the re-assembly and Table S4")
put("p30_panel", 30, "replace_paragraph", "",
    "To build on the findings of Zhou et al. (2022), who reported differences in the abundance of some genes involved in "
    "methylotrophic methanogenesis and methylphosphonate cycling, we compiled a panel of 138 KEGG orthologs (KOs) "
    "involved in the cycling of methylated compounds and in methanogenesis (Table S1), using BioCyc, KEGG and the "
    "literature (49, 50). The panel extends the original list of 39 KOs: 34 were kept, six that KOfam assigns to "
    "unrelated functions were replaced by the KOs of the intended enzymes (Table S1, sheet ‘Original list’), and 68 were "
    "added. The added KOs cover the oxidation of trimethylamine N-oxide (TMAO), dimethylamine (DMA) and methylamines; "
    "the production of trimethylamine (TMA) from choline, carnitine and glycine betaine; the corrinoid methyltransferases "
    "of methylotrophic methanogens and bacteria and methyl-coenzyme M reductase; glycine betaine, dimethylglycine and "
    "sarcosine metabolism; methanol oxidation; the metabolism of dimethylsulfide (DMS), methanethiol and "
    "methanesulfonate; the cleavage and demethylation of dimethylsulfoniopropionate (DMSP); and the C–P lyase and "
    "oxidative pathways of phosphonate degradation. Thirty further KOs of the acetoclastic and hydrogenotrophic pathways "
    "of methanogenesis were analyzed for Figure S4. The abundance of each KO in each metagenome was measured both in all "
    "reads and in the assembled genes, as copies per genome equivalent (see below).",
    "the IMG/MEGAHIT/DESeq2 profile is replaced by the read- and assembly-based profiles of the 2026 data (the pipeline "
    "and versions are in the MAG and read-survey Methods paragraphs below)")
put("p50_stats", 50, "replace",
    "Differences in gene abundances among the three habitat types were tested with Kruskal-Wallis tests.",
    "Differences in gene abundances (copies per genome equivalent) among the three habitat types (unrestored R1 and R2, "
    "n = 12; restored SF2, n = 6; reference R2A, n = 6) were tested with Kruskal–Wallis tests, with Benjamini–Hochberg "
    "false discovery rate (FDR) correction across all panel KOs detected in the 24 sediment metagenomes, followed by "
    "Dunn’s tests with Holm correction; the metagenomes were clustered with Ward’s minimum-variance criterion on the "
    "per-gene z-scores. The abundances of the cultured genera (Figure S6) were compared among the four sites in the same "
    "way.", "statistics of Figure 1, Figure S6 and Table S1")
put("p50_s4", 50, "insert_before", "Only metabolites present in at least three samples",
    "For the correlations with gene abundance (Figure S4), only reads whose best hit was a methanogen were counted for "
    "the genes of the methanogenesis pathways. ", "Figure S4 counts methanogen-lineage reads")
put("p50_mantel", 50, "replace", "Correlations between Bray-Curtis dissimilarities were calculated with Mantel tests.",
    "Correlations between Bray–Curtis dissimilarities were calculated with Mantel tests (Pearson’s r, 9,999 "
    "permutations) over the 19 samples with NMR data, using all 50 metabolites and, for taxonomic composition, either the "
    "species of the re-assembled MAG set or the genera of the ribosomal-protein reads.", "Mantel inputs")

# ---- Results -----------------------------------------------------------------------------------------------------------
put("p57_figure1", 57, "replace_paragraph", "",
    f"Of the 108 methylated-compound KOs in the panel (Table S1), {n_fig1} were detected in the reads of the 24 sediment "
    "metagenomes; only the alcohol oxidase AOX, which has no prokaryotic reference, was absent. Hierarchical clustering "
    "of the metagenomes by these KOs separated the three sample types exactly (unrestored salterns, restored saltern, "
    f"reference wetland; Figure 1), and the abundance of {n_sig} of the {n_fig1} KOs differed among them "
    "(Kruskal–Wallis test, P~FDR~ < 0.05). The unrestored salterns had the most C–P lyase genes (*phnGHIJKLM*; *phnJ* "
    f"{f(m('phnJ', 'unrestored'))} copies per 100 genome equivalents, against {f(m('phnJ', 'restored'))} in the restored "
    f"saltern and {f(m('phnJ', 'reference'))} in the reference wetland), betaine reductase (*grdI* "
    f"{f(m('grdI', 'unrestored'))}, against {f(m('grdI', 'restored'))} and {f(m('grdI', 'reference'))}) and "
    f"methyl-coenzyme M reductase (*mcrA* {f(m('mcrA', 'unrestored'), 2)}, against {f(m('mcrA', 'restored'), 2)} and "
    "< 0.01), and the highest abundance of the genes for DMA oxidation by DMA monooxygenase (*dmmABC*), for "
    "methylamine oxidation through N-methylglutamate (*mgdABCD*), for DMSP demethylation and cleavage (*dmdA*, *dddL*) "
    "and for methanethiol oxidation (*MTO*). The reference wetland had the most genes for the aerobic oxidation of TMA "
    f"and methylamines (*tmd* {f(m('tmd', 'reference'))}, against {f(m('tmd', 'unrestored'))} in the unrestored "
    "salterns; *mauAB*; *gmaS*), carnitine monooxygenase (*cntA*), glycine betaine monooxygenase (*gbcAB*), choline "
    "monooxygenase (*CMO*), methanol dehydrogenases (*mxaF*, *xoxF*), and several genes of formaldehyde, sarcosine "
    "and methylthiopropanoate (MMPA) metabolism (*gfa*, *frmA*, *dgcA*, *dmdB*, *dmdD*). The restored saltern had the "
    f"most choline TMA-lyase (*cutC* {f(m('cutC', 'restored'))}, against {f(m('cutC', 'unrestored'))} and "
    f"{f(m('cutC', 'reference'))}) and betaine-homocysteine methyltransferase (*bhmt*), and the most corrinoid "
    "methyltransferases of TMA, DMA, methylamine and methanol (*mttB*, *mtbB*, *mtmB*, *mtaB*; *mtbB* "
    f"{f(m('mtbB', 'restored'))}, against {f(m('mtbB', 'unrestored'))} and {f(m('mtbB', 'reference'))}). These "
    "methyltransferases occur in methylotrophic methanogens and in bacteria: only "
    f"{pct(res_lo)}–{pct(res_hi)}% of the *mtbB*, *mtmB* and *mtaB* reads of the restored saltern came from methanogens "
    "(most of its *mtbB* and *mtmB* reads were of Desulfobacterota and Chloroflexota), against "
    f"{pct(unr_lo)}–{pct(unr_hi)}% in the unrestored salterns, where *Methanolobus* carried the largest share of the "
    "*mtmB* reads. In sum, the microbiomes at each site differed in their capacity to metabolize each of these "
    "substrates, both in terms of gene relative abundances and pathways utilized.",
    "Figure 1 on the read survey of the 2026 panel; the gene lists of the old paragraph (folD, purU, ftfL, soxBD, sdaA, "
    "gbsB as unrestored; fae, SARDH, sfnG, msuD as reference; ddhA as restored) are not supported by the new profile "
    "(Table S1)")
put("p62_figS4", 62, "replace",
    "TMA, trehalose, acetate, and MMA were significantly positively correlated with mcrABG while sucrose, "
    "2-oxoisocaproate, fructose, glucose, and proline were the most negatively correlated with mcrA (Figure S4).",
    "Counting only the reads of methanogens, MMA, acetate and TMA were significantly positively correlated with *mcrA*, "
    f"*mcrB* and *mcrG* (Spearman ρ = {f(pos_lo, 2)}–{f(pos_hi, 2)} for *mcrA*), and trehalose with *mcrB*, while "
    "fructose, proline, 2-oxoisocaproate and sucrose were negatively correlated with all three "
    f"(ρ = {f(neg_lo, 2).replace('-', '−')} to {f(neg_hi, 2).replace('-', '−')} for *mcrA*), and glucose with *mcrB* and "
    "*mcrG* (Figure S4). The C–P lyase genes shared most of these correlations: *phnJ*, for example, was positively "
    "correlated with TMA and trehalose and negatively with proline, sucrose, glucose and fructose.",
    "Figure S4 on the read survey (methanogen-lineage reads)")
put("p62_mantel", 62, "replace", "(Mantel test, r = 0.78, p < 0.001)",
    f"(Mantel test with the species of the re-assembled MAG set, r = {f(mg['mantel_r'], 2)}, and with the genera of the "
    f"ribosomal-protein reads, r = {f(rg['mantel_r'], 2)}; p < 0.001 for both)", "Mantel on the 2026 community data")
put("p67_figS6", 67, "replace",
    "All of these isolated genera were significantly more abundant in unrestored salterns than the restored saltern and "
    "reference wetland, demonstrating a preference for high salt concentrations (Figure S6).",
    "*Halomonas* and *Salegentibacter* were significantly more abundant in the unrestored salterns than in the restored "
    "saltern and the reference wetland, in the 16S rRNA gene data, the metagenomic reads and the MAGs, consistent with a "
    "preference for high salt concentrations, whereas *Idiomarina* was rare at all sites and differed among them only in "
    "the 16S rRNA gene data (Figure S6).", "Figure S6 with the read and MAG rows")

# ---- Discussion --------------------------------------------------------------------------------------------------------
put("p81_bc", 81, "replace",
    "(Bray-Curtis index = 0.31 to 0.85) in even the most dissimilar taxonomic communities (Bray-Curtis index = 0.83 to "
    "~0.96)",
    f"(Bray–Curtis index = {f(bc_m[0], 2)} to {f(bc_m[1], 2)}) in even the most dissimilar taxonomic communities "
    f"(Bray–Curtis index = {f(bc_t[0], 2)} to {f(bc_t[1], 2)} between the species compositions of the re-assembled MAG "
    "set)", "Bray-Curtis ranges on the 2026 community data")
put("p93_genera", 93, "replace",
    "The cultivated genera were all present in both 16S rRNA gene and metagenomic sequencing datasets and were more "
    "abundant in the unrestored salterns (Figure S6), demonstrating their preference for this environment.",
    "All five cultivated genera were detected in both the 16S rRNA gene and the metagenomic data, but in the metagenomic "
    "reads only *Halomonas* and *Salegentibacter*, and *Marinilabilia* relative to the reference wetland, were more "
    "abundant in the unrestored salterns (Figure S6); *Idiomarina* was rare at all sites, and the reads assigned to "
    "*Clostridium*, most of them distant from any reference genome, were more abundant in the restored saltern and the "
    "reference wetland.", "Figure S6 with the read and MAG rows")

# ---- captions ----------------------------------------------------------------------------------------------------------
CAP1 = ("Figure 1. Relative abundance of genes involved in the cycling of methylated compounds in the 24 sediment "
        "metagenomes: copies per genome equivalent in all reads, z-scored per gene. Samples (columns) are clustered with "
        "Ward’s criterion, which separates them by restoration type; genes (rows) are grouped by compound, in two "
        "blocks. Rows are annotated by whether the gene’s abundance differs among the three restoration types (KW: "
        "Kruskal–Wallis test, Benjamini–Hochberg FDR < 0.05 over all panel KOs), the restoration type with the highest "
        "abundance (Highest; shown when FDR < 0.05; pairwise Dunn’s tests in Table S1), the compound, the type of "
        "process, and the validation of the read survey (Validation: flagged = recall < 0.80 or precision < 0.90 on "
        "simulated reads, or not separable in short reads from a related gene; not testable = too few reference genes). "
        "The corrinoid methyltransferases (methyl transfer) occur in methylotrophic methanogens and in bacteria. The "
        "alcohol oxidase AOX (K17066) was not detected. TMA = trimethylamine, TMAO = trimethylamine N-oxide, DMA = "
        "dimethylamine, MMA = monomethylamine, GB = glycine betaine, DMG = dimethylglycine, DMS = dimethyl sulfide, "
        "MMPA = methylthiopropanoate, DMSP = dimethylsulfoniopropionate.")
CAPS4 = ("Figure S4. Spearman correlations between metabolite concentrations (NMR) and the abundance of methanogenesis "
         "and C–P lyase genes in the 19 sediment samples with NMR data, for the 40 metabolites present in at least three "
         "samples. Gene abundance is copies per genome equivalent in all reads; for the methanogenesis pathways only "
         "reads whose best hit is a methanogen (an archaeon of an order, family or genus named Methano-, excluding "
         "anaerobic methanotrophs) are counted, so that bacterial and haloarchaeal homologs are excluded; C–P lyase "
         "counts all reads. Three phosphotransacetylase genes (K13788, K15024 and *eutD*) had no methanogen reads and "
         "are not shown. Dots mark correlations with Benjamini–Hochberg FDR < 0.05 over all 1,960 tests. Metabolites "
         "are clustered by complete linkage of their correlation profiles.")
put("cap_fig1", 115, "replace_paragraph", "", CAP1, "Figure 1 on the read survey")
put("cap_figS4", 140, "replace_paragraph", "", CAPS4, "Figure S4 on the read survey")
put("cap_figS6_a", 144, "replace", " and metagenomic data (bottom row)",
    ", metagenomic reads (middle row) and metagenome-assembled genomes (MAGs; bottom row) of the 24 sediment samples",
    "Figure S6 rows (as suggested in the Google Doc on 2026-11-25)")
put("cap_figS6_b", 144, "replace",
    "Metagenomic abundances were computed by the IMG/M pipeline by using lastal to map contigs to the IMG reference "
    "database. ", "", "Figure S6: the IMG/lastal row is replaced")
put("cap_figS6_c", 144, "append", "",
    " Metagenomic abundance is the percentage of the reads of eight single-copy "
    "ribosomal-protein genes whose best DIAMOND hit is a reference genome of the genus (GTDB R232); *Halomonas* includes "
    "*Vreelandella*, which GTDB separates from it, and *Clostridium* includes the GTDB genera with *Clostridium* "
    "suffixes (e.g. Clostridium_B). Only 8% of the reads assigned to *Clostridium* had a best hit at least 95% identical, "
    "so that assignment is uncertain. MAG abundance is the percentage of all reads that mapped to the genus’s species in "
    "the re-assembled MAG set; no species of *Idiomarina*, *Clostridium* or *Marinilabilia* was recovered. In the "
    "metagenomic and MAG rows, letters show Dunn’s tests (Holm-adjusted) after a significant Kruskal–Wallis test across "
    "the four sites (p < 0.05); sites sharing a letter do not differ.", "Figure S6 definitions (as in the Google Doc)")
put("cap_tableS1", 169, "replace_paragraph", "",
    "Table S1. The 138 KEGG orthologs involved in the cycling of methylated compounds and in methanogenesis that were "
    "analyzed in this study (Excel file), with their abundance in the reads and in the assembled genes of the 24 sediment "
    "metagenomes by restoration type, Kruskal–Wallis and Dunn’s tests, the validation of the read survey, and the "
    "species of the re-assembled MAG set that carry them. A second sheet lists the 39 KOs of the original list and how "
    "each was carried into the panel.", "Table S1 is the 2026 panel")
put("cap_tableS4", 175, "replace_paragraph", "",
    "Table S4. The 39 metagenomes of this study (26 sediment and 13 water samples): NCBI BioProject, BioSample and SRA "
    "accessions, GOLD and IMG/M identifiers, read counts, and the MAGs kept from each re-assembly (Excel file).",
    "Table S4 covers all 39 metagenomes")

out = copy.deepcopy(base)
out["edits"] = base["edits"] + list(NEW.values())
out["figures_2026_11_26"] = {
    "Figure 1": {"image": "reassembly_2026/community_ko/figures/figure1_2026.png", "caption_paragraph": 115,
                 "image_inline_in_caption": True},
    "Figure S4": {"image": "reassembly_2026/community_ko/figures/figureS4_2026.png", "caption_paragraph": 140,
                  "image_inline_in_caption": True},
    "Figure S6": {"image": "reassembly_2026/community_ko/figures/figureS6_2026.png", "caption_paragraph": 144,
                  "image_paragraph": 143}}
out["note"] = ("r232_review_edits_2026-11-25.json plus the 2026-11-26 metagenome sections "
               "(make_metagenome_edits_2026-11-26.py).")
out["facts_2026_11_26"] = json.loads(json.dumps(FACTS, default=float))
json.dump(out, open(HERE / "r232_review_edits_2026-11-26.json", "w"), indent=1, ensure_ascii=False)
print(f"wrote r232_review_edits_2026-11-26.json: {len(out['edits'])} edits ({len(NEW)} new)")
print(NEW["p57_figure1"]["new"]); print(); print(NEW["p62_figS4"]["new"])
