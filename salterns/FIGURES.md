# Salterns manuscript: figures and tables, and the scripts that make them

Paths are relative to `salterns/`. Run every Python script with `~/Documents/py_venv/bin/python`. "Original" means the figure in the Google Doc was made by the authors outside this repository and was not regenerated in the 2026 revision.

The read-survey version behind every read-based number is named in `reassembly_2026/read_survey_current.json`, and all downstream scripts read it.

The review document `Metabolic Manuscript - R232 review (tracked changes).docx` is built by `clauded_insights/make_r232_review_docx.py`, which places the images listed below. The Google Doc suggestions come from `clauded_insights/gdoc_sync/tools/`.

## Main figures

| Figure | Script (run in this order) | Output |
|---|---|---|
| 1 | `reassembly_2026/community_ko/figure1.py` (needs the v2 read survey; statistics in `community_ko/stats_lib.py`) | `reassembly_2026/community_ko/figures/figure1_2026.png` |
| 2, 3, 4, 6 | Original (NMR heatmap, metabolite–methane correlations, substrate additions, isolate growth) | — |
| 5 | Original (KBase SpeciesTree, narrative 139158); only the caption is revised | — |
| 7 | `reassembly_2026/modelling_2026/simulations/per_gram_yields.py` → `specific_rates.py` → `figure7_pergram.py` | `reassembly_2026/modelling_2026/simulations/figure7_pergram.png` |
| 8 | `correlational/final585/rerun_correlational_final585.py --parts network`, which calls `correlational/louvain_cooccurrence_network.py`; the layout is deterministic since 2026-11-25 | `correlational/final585/within_condition_louvain_network_final585.png` |
| 9 | `correlational/final585/rerun_correlational_final585.py --parts guilds,corr` | `correlational/final585/guild_methane_raw_final585.png` |
| 10 | `reassembly_2026/modelling_2026/synthesis/pathway_summary.py` → `figure10_pathways.py` | `reassembly_2026/modelling_2026/synthesis/figure10_pathways.png` |

## Supplementary figures

| Figure | Script | Output |
|---|---|---|
| S1, S2, S3, S5, S7, S9 | Original (site photo, MetaboAnalyst ORA, NMR PCA, total gas, soil DNA vs methane, microcosm rates) | — |
| S4 | `reassembly_2026/community_ko/figureS4.py` (needs the v2 read survey; metabolite table from `community_ko/metabolome_taxonomy.py`) | `reassembly_2026/community_ko/figures/figureS4_2026.png` |
| S6 | `reassembly_2026/community_ko/figureS6.py`; the 16S row is cropped from the original figure (`community_ko/figures/figureS6_original_2021.png`), the read row comes from `read_markers_v2/results/marker_genus_profile.tsv` (`read_markers_v2/scripts/11c_marker_genus_profile.py`, run on poplar) | `reassembly_2026/community_ko/figures/figureS6_2026.png` |
| S8 | `correlational/final585/rerun_correlational_final585.py --parts guilds,corr` | `correlational/final585/guild_methane_dnanorm_final585.png` |
| S10 | `reassembly_2026/modelling_2026/figureS10/figureS10_final.py` | `reassembly_2026/modelling_2026/figureS10/figureS10_final.png` |
| S11 | `rerun_correlational_final585.py --parts heatmap,network` → `correlational/final585/figureS11_final585.py` | `correlational/final585/figureS11_final585.png` |

## Tables

| Table | Script | Output |
|---|---|---|
| 3 | `rerun_correlational_final585.py --parts table3` | `correlational/final585/table3_final585.json`, set into the review by `make_r232_review_docx.py` |
| 4, 5, S6, S7 | `clauded_insights/make_r232_review_docx.py` (from `correlational/final585/` and `reassembly_2026/manuscript_numbers_*.json`) | in the review docx |
| S1 | `reassembly_2026/community_ko/` (in progress: the KO panel `ko_panel_methylated.tsv` with read- and assembly-based statistics) | — |
| S4 | `reassembly_2026/community_ko/make_table_s4.py` | `Table S4 - metagenomes, accessions and re-assembly (2026).xlsx` |
| S5 | `correlational/r232/make_table_s5.py` (rebuild; `--compare` checks it against the committed workbook) | `Table S5 - MAG taxonomy, quality and guilds (GTDB R232).xlsx` |
| S8 | `reassembly_2026/make_table_s8.py` | `Table S8 - re-assembled MAG set taxonomy, quality and guilds (GTDB R232).xlsx` |

## Reproducibility check (2026-11-25)

Each script above was run again under a Python audit hook that recorded every file it opened.
- Every input is tracked in git. The exceptions are the poplar-side intermediates of the read surveys, which are summarized in `reassembly_2026/read_markers*/README.md`.
- The PNG outputs came out byte-identical. SVGs and xlsx files differ only in their embedded dates.
- Figure 8 was the one exception, because its layout depended on Python's per-process string hashing. Its network script now orders modules and subgraphs deterministically; the network, modules and statistics are unchanged.
