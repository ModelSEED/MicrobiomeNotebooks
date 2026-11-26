# Community gene profiles of the methylated-compound panel (2026)

This folder replaces the original IMG/MEGAHIT/DESeq2 functional profile of the 24 sediment metagenomes. It holds the 138-KO panel and two independent abundance profiles of it, built on the 2026 data, plus the figures and tables made from them: Figure 1, Figure S4, Figure S6, Table S1 and Table S4. Run every script with `~/Documents/py_venv/bin/python` from this folder.

## Inputs

**The panel.** `ko_panel_methylated.tsv` lists 138 KOs:
- 34 of the 39 KOs in the original Table S1, kept, 6 of them relabeled;
- 6 KOs that replace originals KOfam assigns to unrelated functions: K02400 FlhA → K00121 frmA; K16259 → K14028/K14029 mxaFI; K03781 catalase → K17066 AOX; K00219 → K00148 fdhA; K08967 → K20034 dmdB;
- 68 added KOs;
- 30 KOs of acetoclastic and hydrogenotrophic methanogenesis, used for Figure S4.

`table_s1_original_39.tsv` holds the original list.

**Read-based profile.** The read survey named in `../read_survey_current.json`, currently `../read_markers_v2/`:
- DIAMOND blastx of all 2.08 billion reads against the panel references and decoys, validated on held-out genes;
- copies per genome equivalent, from eight ribosomal-protein genes;
- archaeal and methanogen-lineage reads;
- merged classes for genes that short reads cannot separate.

**Assembly-based profile.** `assembly_based/`:
- KOfam on all 12.06 million genes of the 39 second-pass assemblies, using a two-pass hmmsearch whose results are identical to a full search;
- copies per genome equivalent from contig depth.

See `assembly_based/README.md` for details.

## Scripts (in run order)

| Script | Output |
|---|---|
| `make_table_s4.py` | `../../Table S4 - metagenomes, accessions and re-assembly (2026).xlsx` |
| `panel_stats.py --source read` / `--source assembly` | `results/panel_stats_{read,assembly}.tsv`: Kruskal–Wallis across unrestored (n = 12), restored (n = 6) and reference (n = 6), BH FDR over all detected panel KOs, Dunn/Holm letters |
| `figure1.py` (`--source assembly` for the assembly version) | `figures/figure1_2026.png`, `results/figure1_stats.tsv`, `results/figure1_meta.json` |
| `metabolome_taxonomy.py` | `results/metabolome_taxonomy.json` (Mantel and Bray–Curtis ranges, metabolome vs MAG species and read genera) |
| `figureS4.py` | `figures/figureS4_2026.png`, `results/figureS4_rho.tsv`, `results/figureS4_meta.json` |
| `figureS6.py` | `figures/figureS6_2026.png`, `results/figureS6_{values,stats,read_sensitivity}.tsv` |
| `table_s1.py` | `../../Table S1 - methylated-compound gene panel (2026).xlsx` |

Shared modules:
- `profiles.py`: panel annotations, profile loaders, validation status;
- `stats_lib.py`: Kruskal–Wallis, Dunn, Holm, compact letters, BH;
- `heatmap_lib.py`: pheatmap-style drawing.

## Results (read survey v2, 2026-11-26)

**Figure 1:**
- 107 of the 108 methylated-compound KOs were detected in the 24 sediment metagenomes. AOX has no prokaryotic reference.
- 100 of them differ among restoration types (FDR < 0.05).
- Ward clustering separates the samples exactly by restoration type.

**By site type** (copies per 100 genome equivalents):
- **Unrestored:** highest in C–P lyase (phnJ 12.7, against 1.1 restored and 4.4 reference), grdI (6.5), mcrA (0.28), dmmABC and mgdABCD.
- **Reference:** highest in tmd, mauAB, gmaS, cntA, gbcAB, CMO and the methanol dehydrogenases.
- **Restored:** highest in cutC, bhmt and the corrinoid methyltransferases. Only 3–8% of those methyltransferase reads (mtbB, mtmB, mtaB) are from methanogens there, mostly Desulfobacterota and Chloroflexota, against 21–52% in the unrestored cores.

**Figure S4** (methanogen-lineage reads; 19 samples, 40 metabolites):
- mcrA is positively correlated with methylamine, acetate and TMA, and negatively with fructose, proline, 2-oxoisocaproate and sucrose (FDR < 0.05).
- C–P lyase shares most of these correlates.

**Mantel test** (metabolome vs community, Bray–Curtis on square-root-transformed data, all 50 metabolites, 19 samples): r = 0.78 with the MAG species and 0.70 with the read genera (p = 0.0001).

**Figure S6:**
- Halomonas (with GTDB Vreelandella) and Salegentibacter are higher in R1 and R2 than in SF2 and R2A, in both the read and the MAG rows.
- Marinilabilia is higher in R1 and R2 than in R2A.
- Idiomarina is rare and does not differ among sites.
- Reads assigned to Clostridium (8% with a best hit ≥ 95% identical) are higher in SF2 and R2A.

**Assembly vs reads:** where genes assemble, the two profiles agree (ρ 0.74–0.95 for phnJ, grdI, cntA and mtmB). But the assemblies hold only 0.3–0.5 of the read-based copies in the unrestored cores, and less elsewhere, so the figures use the read-based profile.

The text that reports these results is in `../../clauded_insights/make_metagenome_edits_2026-11-26.py`.
