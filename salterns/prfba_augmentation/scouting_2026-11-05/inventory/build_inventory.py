#!/usr/bin/env python
"""Writes inventory.json + inventory.csv: every salterns R232 analysis, whether it depends on MAG gene content, and
how to point it at another annotation source. Line numbers refer to the working tree on 2026-11-04 (branch
clauded-insights-figure7). Paths are relative to /home/freiburger/Documents/MicrobiomeNotebooks/salterns unless absolute.
"""
import json
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parent
SCR = str(OUT)
W = f"{SCR}/rerun_with_alt_annotations.py"

A = []


def add(**k):
    A.append(k)


add(id="guilds", script="correlational/r232/rerun_correlational_r232.py --parts guilds (build_guilds, l.106-146)",
    annotation_dependent="yes (defines every guild)",
    inputs=["jplfaria/Cliff MAGs/RASTified_MAG_json/<mag>.contigs__.RAST.json (RAST l.47, path l.110)",
            "correlational/iterativeIDs_r232.json (l.73), iterativeIDs.json (l.74)", "data/Saltern_phylogeny_r232.json (l.75)",
            "data/Cliff_310MAG_relabund.tsv (l.76), data/SaltPondsDNA.csv (l.77), data/Cliff_Sample_Metadata_BGC_NMR.csv (l.78)",
            "MAG_QC/checkm/summary/CheckM_summary_table.tsv (l.118)", "correlational/analysis_2_marimo.py 2020 sets 'mepe' l.163-190, 'betaine' l.194-208 (read by marimo_set l.99-102)"],
    outputs=["correlational/r232/guilds_r232.csv (l.137)"],
    runtime="measured 11 s total with the regex parser (import 1.7 s, read 310 RAST JSONs 6.6 s, build 1.9 s); 58 s to read with json.load",
    swap_point="RAST constant l.47 + file pattern l.110; HERE (l.44) must be redirected or the CSV is written into the repo",
    one_line="import rerun_correlational_r232 as rr; rr.RAST = Path(ALT_DIR); rr.HERE = Path(OUT); g = rr.build_guilds()   # ALT_DIR holds <mag>.contigs__.RAST.json with KBase-genome 'functions' lists",
    wrapper=f"{W} --out OUT --annotations ALT_DIR_or_JSON [--add-functions EXTRA.json] [--parser json|trunc] --parts guilds",
    hard_coded_breakpoints=["assert 310 MAGs all with r232 IDs (l.80) only"],
    verified="identity run (original RAST, both parsers) reproduces guilds_r232.csv byte-for-byte",
    manuscript=["Table 4 members", "Table S5 guild flags + annotation columns", "every downstream item"])
add(id="guild_corr", script="rerun_correlational_r232.py --parts corr (guild_correlations l.156-183, per_mag_correlations l.186-201, plot_guilds l.204-247)",
    annotation_dependent="yes, through guild membership (per-MAG rho/q are abundance-only; guild flags merged in at l.199)",
    inputs=["guilds_r232.csv frame g", "abundance/CH4 frames built at import (l.86-96)"],
    outputs=["guild_methane_correlations_r232.csv (l.176)", "methanogen_archaea_by_core_r232.csv (l.182, sums g.methanogen and g.archaea members)",
             "mag_methane_correlations_r232.csv (l.200)", "guild_methane_raw_r232.png, guild_methane_dnanorm_r232.png (l.245-246)"],
    runtime="measured 159 s (4,000 random-set draws N_NULL l.50, 8 guilds x 2 bases x 2 core sets)",
    swap_point="takes g from build_guilds; outputs to rr.HERE",
    one_line="corr = rr.guild_correlations(g); rr.per_mag_correlations(g); [rr.plot_guilds(g, corr, b) for b in rr.AB]",
    wrapper=f"{W} --out OUT ... --parts guilds,corr",
    hard_coded_breakpoints=["none in code; but one RNG stream default_rng(0) (l.157) is shared by all guilds in GUILDS order (l.149-152): verified with the synthetic test (phnJ 22 -> 24 MAGs) that null_frac_ge_rho shifts for every guild, including unchanged ones (e.g. betaine_reducer DNA-normalised 0.202 -> 0.1875, archaea raw 0.3155 -> 0.2858), while their rho values are unchanged"],
    verified="identity run: all four CSVs byte-identical to the committed files",
    manuscript=["Table 5 (make_r232_review_docx.py l.619-632)", "Figure 9 = guild_methane_raw_r232.png (l.658)", "Figure S8 = guild_methane_dnanorm_r232.png (l.668)",
                "Table S6 = methanogen_archaea_by_core_r232.csv (l.434, 704-711)", "Table S5 per-MAG rho/q", "Results P78 numbers (l.498-516, partly literal)"])
add(id="heatmap_4site", script="rerun_correlational_r232.py --parts heatmap (heatmap l.276-371)",
    annotation_dependent="cosmetic only: label colours/boxes and legend from taxon_guilds (l.258-265, 324-336, 356-361); statistics abundance-only",
    inputs=["relabund + DNA", "guilds g", "correlational/iterativeID_levels_r232.json (l.326)"],
    outputs=["abundance_correlatons_one_triangle_reduced_defined_r232.png (l.366)", "significantly_connected_organisms_r232.json (l.298, guild-independent)"],
    runtime="not run; ~50 s inferred from mtimes (guilds_r232.csv 02:31:15 -> heatmap png 02:32:04, 2026-10-21), 90x100 in figure at 300 dpi (l.307, 367)",
    swap_point="same rr.RAST/rr.HERE monkeypatch; then rr.heatmap(g)", one_line="rr.heatmap(g)", wrapper="not wrapped (not manuscript-facing)",
    hard_coded_breakpoints=[], verified="not run", manuscript=["not used by make_r232_review_docx.py"])
add(id="within_condition_network", script="rerun_correlational_r232.py --parts network (within_condition l.374-487, network l.490-498 -> ../louvain_cooccurrence_network.py)",
    annotation_dependent="no (build_guilds is not called for --parts network, l.505)",
    inputs=["relabund, DNA, iterativeIDs_r232.json, Saltern_phylogeny_r232.json, averaged_normalized_MAG_abundances.csv"],
    outputs=["all_candidate_pairs_r232.csv, direct_interactions_r232.csv, glasso_direct_edges_r232.csv, C_{marginal,within,between}_r232.csv, within_condition_meta_r232.json, marginal_vs_within_r232.png, within_condition_louvain_network_r232.png, within_condition_louvain_modules_r232.json"],
    runtime="not run; ~2.2 min estimated: same pipeline logged [131.7 s] for the baseline in prefilter_sensitivity/within_run.log, Louvain ~7 s (C_within 02:29:42 -> network png 02:29:49)",
    swap_point="none needed", one_line="", wrapper="not needed", hard_coded_breakpoints=[], verified="not run (annotation-independent)",
    manuscript=["Figure 8 image (make_r232_review_docx.py l.640)", "Figure S11 inputs", "Figure 8 caption names the guild taxa retained in the network as literals (l.652-654)"])
add(id="figureS11", script="correlational/r232/figureS11_r232.py",
    annotation_dependent="highlights only (colour, box, legend counts) from guilds_r232.csv (l.35-43); matrix/asterisks from annotation-independent files",
    inputs=["C_within_r232.csv (l.33)", "glasso_direct_edges_r232.csv (l.34)", "guilds_r232.csv (l.35)"],
    outputs=["figureS11_r232.png (l.111)", "figureS11_r232.order.json (l.112-114)"], runtime="measured 4.8 s",
    swap_point="HERE l.28 (inputs and outputs)", one_line="run with HERE patched to a folder holding the alt guilds_r232.csv + links to the two static CSVs",
    wrapper=f"{W} --out OUT ... --parts guilds,s11",
    hard_coded_breakpoints=["assert the three highlighted guilds share no taxon (l.45)"],
    verified="identity run: PNG and order.json byte-identical", manuscript=["Figure S11 (make_r232_review_docx.py l.676-692; caption numbers 245/169/13/186 are annotation-independent)"])
add(id="membership_2020_vs_r232", script="correlational/r232/make_guild_membership_comparison.py",
    annotation_dependent="yes: membership from guilds_r232.csv (l.117) + extra RAST strings (EXTRA l.48-58) read with the string-aware parser rast_functions (l.73-77, used l.225-227)",
    inputs=["guilds_r232.csv (l.117)", "RAST JSON (RAST l.45)", "MAG_QC/comparison_vs_2020.csv (l.122)", "correlational/iterativeIDs_r232_vs_2020.csv (l.123)",
            "data/Cliff_310MAG_relabund.tsv (l.125)", "correlational/saltern_fluxes.json, saltern_fluxes_betaine.json (l.131-138)", "models/MePn_consumers/*.json (l.140)",
            "correlational/models/*_gapfilled_betaine.json + *_betaine.json (l.141, 150-155)", "correlational/analysis(2).ipynb recorded output (l.143-147)",
            "correlational/betaine_grdH_r232/fva_tma_by_core_r232.csv (l.158-164)"],
    outputs=["guild_membership_2020_vs_r232.csv (l.261)", "guild_membership_2020_vs_r232.md (l.518)"], runtime="measured 1.8 s",
    swap_point="HERE l.42 (reads guilds_r232.csv, writes outputs), COR l.43, RAST l.45",
    one_line="run with HERE -> alt folder and COR -> repo correlational/",
    wrapper=f"{W} --out OUT ... --parts guilds,membership",
    hard_coded_breakpoints=["assert rows per guild == {MePn 27, betaine 15, methanogen 1, TMA 4} (l.265)", "assert 20 MePn / 15 betaine modelled (l.266-267)",
                            "assert no 'added' member (l.268)", "assert exactly one Halanaerobium with FVA min 0 (l.371)", "fva_tma file must cover exactly the 15 betaine MAGs (l.160-161)",
                            "md prose literals, e.g. 'Both gene-based sets are strict subsets' (l.391), 'No model file (3)' (l.451), 'Zero growth (4)' (l.456)"],
    verified="identity run: CSV and MD byte-identical; synthetic test adding the PhnJ string to JAWOCM01.1 stops at l.265 with AssertionError {'MePn consumer': 28, ...}",
    manuscript=["Table 4 last column summarises it but is typed as literals (make_r232_review_docx.py l.605-606)", "cited in R232 review notes l.94 and Google Doc comments review l.115"])
add(id="archaeal_mepn_tma_report", script="correlational/r232/make_archaeal_mepn_tma_report.py",
    annotation_dependent="yes, gene-level: RULES (l.67-148), REF_FAMS (l.163-168), X_FLAGS (l.171-172) on each feature's ' / '-joined functions (l.222); targets from guilds_r232.csv (l.183, 207-210)",
    inputs=["guilds_r232.csv (l.183)", "RAST JSON incl. location + protein_translation (RAST l.49, load_genome l.213-233)", "MAG_contigs/<mag>.contigs__.fna (l.237, 289)",
            "relabund, DNA, CheckM, GTDB-Tk summaries (l.185-196)", "saltern_fluxes*.json (l.197-201)", "models/MePn_consumers (l.202), correlational/models (l.435-445)"],
    outputs=["archaeal_mepn_and_tma_candidates.md (l.54)", "archaeal_mepn_and_tma_candidates_genes.csv (l.55)"], runtime="measured 106 s",
    swap_point="HERE l.46 / COR l.47 swap the targets only; the gene inventory needs real coordinates on MAG contigs, which prFBA-merged features lack (create_synthetic_genome gives location [[<id>.contig, 1, '+', len]], codiffusion_bioreactor/bvbrc_to_kbase_genome.py l.1157), so contig_seqs(mag)[contig] (l.314) would KeyError",
    one_line="not a swap: keep RAST=original and only point HERE at the alt guilds_r232.csv",
    wrapper=f"{W} --out OUT ... --parts guilds,archaeal  (targets from the alt guilds, genes from the original RAST)",
    hard_coded_breakpoints=["next(i for i in A if i.startswith('Natrinema'/'Haloquadratum'/'Halorubrum')) l.747", "next(... 'CAZFJG01'/'J153'/'UBA11579'/'T1Sed10') l.769-772",
                            "synteny lookups need POS_G/phnJ/phnM/POS_X/phnO-like rows (l.750-754)", "asserts l.792, 794, 1209, 1299", "narrative is hand-written around the current members"],
    verified="identity run: genes CSV byte-identical; MD differs only in the generation date (l.814)",
    manuscript=["not in the manuscript; Google Doc comments review l.138 suggests it as a Discussion source"])
add(id="prefilter_abundance", script="correlational/r232/prefilter_sensitivity_abundance.py",
    annotation_dependent="partly: guild roles label key organisms and drive the guild-sum sensitivity (l.91-99, 189-196, 534-590, 652-688); filter statistics abundance-only",
    inputs=["imports rerun_correlational_r232 (l.51-52)", "guilds_r232.csv (l.71)", "mag_methane_correlations_r232.csv (l.457)", "guild_methane_correlations_r232.csv (l.532)"],
    outputs=["prefilter_sensitivity/abund_* (OUT l.61; mkdir at import l.62)"], runtime="not run, not documented (unverified)",
    swap_point="HERE l.50 (inputs, OUT, and sys.path for the rr import)", one_line="run patched with HERE -> alt folder that already holds the corr outputs",
    wrapper=f"{W} --out OUT ... --parts guilds,corr,prefilter_abund (implemented, not executed)",
    hard_coded_breakpoints=["abundance-only asserts (l.224-225, 263, 353-354, 433-434, 497) unaffected", "consistency asserts vs the corr outputs (l.462, 569) hold if those are regenerated together"],
    verified="not run", manuscript=["conclusions quoted as literals in the review prose (make_r232_review_docx.py l.536-538)"])
add(id="prefilter_within", script="correlational/r232/prefilter_sensitivity_within.py",
    annotation_dependent="labels only (ROLES/KEY_TAXA from guilds_r232.csv l.78, 113-125); all statistics abundance-only",
    inputs=["guilds_r232.csv (l.78)", "relabund, DNA, iterativeIDs (l.74-77)", "committed r232 network files for verification (l.526-541)"],
    outputs=["prefilter_sensitivity/within_* (OUT l.55, mkdir at import l.56)"], runtime="not run; docstring: ~15 min (14 settings) + ~20 min (--null 20) + --variance",
    swap_point="HERE l.52 (guilds + OUT), COR l.53", one_line="patch HERE/COR", wrapper="not wrapped", hard_coded_breakpoints=[], verified="not run",
    manuscript=["review prose literals (e.g. 'Relaxing the prevalence filter gave no edge to any guild taxon', make_r232_review_docx.py l.536-538)"])
add(id="betaine_grdH_rerun", script="correlational/betaine_grdH_r232/rerun_betaine_grdH_r232.py",
    annotation_dependent="yes, but only as a subset filter and diagnostics: grdh9 = the 15 stored betaine models whose MAG has betaine_reducer (l.266-269); GRD_PATTERNS (l.198-203, re.I) count grd genes in the RAST JSON (l.232-238). The FBA models themselves are fixed 2020 files",
    inputs=["guilds_r232.csv (GUILDS l.65)", "RAST JSON (l.61)", "correlational/models/<mag>.contigs__.RAST{,_gapfilled,_gapfilled_betaine,_betaine}.json (l.97, 211-214)",
            "correlational/saltern_fluxes_betaine.json (l.64)", "inputs/averaged_normalized_abundances_unrestored_22f4629.csv (l.62)", "data/averaged_normalized_abundances_unrestored.csv (l.63)", "data/Saltern_phylogeny.json (l.66)"],
    outputs=["fluxes_grdH_r232.json (l.374)", "reproduction_check_r232.csv (l.306)", "fva_tma_by_core_r232.csv (l.340)", "diagnostics_grdH_r232.csv (l.342)"],
    runtime="not run here; README l.18 says ~75 s (cobra 0.29.1, GLPK; Gurobi FVA if licensed, l.278-285; Gurobi thread count not checked)",
    swap_point="GUILDS l.65; HERE l.57 (outputs; AB_AS_RUN l.62 derives from HERE)", one_line="patch GUILDS and HERE (and AB_AS_RUN back to the repo path)",
    wrapper="not wrapped (FBA; would need a thread cap)",
    hard_coded_breakpoints=["assert len(grdh9) == 9 (l.269)", "a MAG that gains grdH outside the 15 stored betaine models has no model to simulate"],
    verified="not run", manuscript=["feeds fva regimes to make_guild_membership_comparison.py and modelling_numbers_r232.py; betaine text P74-P76 via edits JSON"])
add(id="betaine_stats", script="correlational/betaine_grdH_r232/betaine_stats_grdH_r232.py",
    annotation_dependent="indirect (grdH9 list from fluxes_grdH_r232.json meta, l.139)",
    inputs=["fluxes_grdH_r232.json (l.138)", "saltern_fluxes*.json (l.136-137)", "clauded_insights/figure7_stats.json (l.31, asserts l.160-162)"],
    outputs=["stats_grdH_r232.json (l.245)"], runtime="README l.18: a few seconds (not run)", swap_point="HERE l.27", one_line="", wrapper="not wrapped",
    hard_coded_breakpoints=["asserts reproduction of figure7_stats.json (l.160, 162)"], verified="not run",
    manuscript=["superseded basis (betaine_grdH_r232/README.md l.47); read by modelling_numbers_r232.py"])
add(id="modelling_consistent_basis", script="clauded_insights/modelling_consistent_basis_r232.py",
    annotation_dependent="yes, as subsets of the STORED 2020 fluxes: phnJ17 and grdH9 (l.92-100); 'MePn_0.07_20' (all 20 models) is guild-independent",
    inputs=["correlational/r232/guilds_r232.csv (l.48)", "correlational/saltern_fluxes.json, _betaine.json, _addition.json (l.89-91)", "betaine_grdH_r232/inputs/...22f4629.csv (l.44)", "data/* (l.42-45)"],
    outputs=["clauded_insights/modelling_consistent_basis_r232.json (l.127)"], runtime="measured ~1 s (modelling part 2.9 s incl. Figure 7 and numbers)",
    swap_point="GUILDS l.48, HERE l.34, S l.35", one_line="patch GUILDS path and HERE/S",
    wrapper=f"{W} --out OUT ... --parts guilds,modelling",
    hard_coded_breakpoints=["assert phnJ17 has 17 and nonbinding3 has 3 (l.104)", "assert grdH9 has 9 (l.105)", "a gained phnJ/grdH MAG without a stored flux is silently absent"],
    verified="identity run: JSON byte-identical", manuscript=["Figure 7 (betaine series = betaine_grdH9)", "edits JSON numbers (P36, P69-P77, P129)"])
add(id="figure7", script="clauded_insights/figure7_r232.py", annotation_dependent="indirect (betaine_grdH9 series)",
    inputs=["clauded_insights/modelling_consistent_basis_r232.json (l.30)"], outputs=["figure7_r232.png/.svg (l.68-69)", "figure7_r232_stats.json (l.70-71)"],
    runtime="measured < 2 s", swap_point="HERE l.26", one_line="", wrapper=f"{W} --parts guilds,modelling",
    hard_coded_breakpoints=["labels '20 models' / '9 grdH carriers' are literals (l.32-33)"],
    verified="identity run: stats numbers identical; the committed figure7_r232_stats.json 'basis' text predates the 2026-11-04 docstring edit of modelling_consistent_basis_r232.py; PNG differs in 565 of 4.0M pixels (label placement)",
    manuscript=["Figure 7 (make_r232_review_docx.py l.637)"])
add(id="modelling_numbers", script="clauded_insights/modelling_numbers_r232.py", annotation_dependent="yes, subsets of stored fluxes (phnJ17, grdH9) + genus/family shares from guilds (l.116-129, 176-179)",
    inputs=["guilds_r232.csv (l.116)", "stored fluxes", "betaine_grdH_r232/fva_tma_by_core_r232.csv, stats_grdH_r232.json", "clauded_insights/figure7_stats.json (l.245)"],
    outputs=["clauded_insights/modelling_numbers_r232.json (OUT l.50)"], runtime="measured ~1 s", swap_point="G l.116, HERE l.47, S l.48",
    one_line="patch G path and HERE/S", wrapper=f"{W} --parts guilds,modelling",
    hard_coded_breakpoints=["assert 20/17/15/9 (l.129)", "curated-not-modelled MAG lists typed in (l.163-166)"], verified="identity run: JSON byte-identical",
    manuscript=["traces the draft R2 values; superseded basis for model-vs-measured"])
add(id="review_docx", script="clauded_insights/make_r232_review_docx.py",
    annotation_dependent="Tables 4/5/S6 and Figures 7/9/S8/S11 are read from the r232 outputs; much guild prose is literal",
    inputs=["correlational/r232/guild_methane_correlations_r232.csv (l.433)", "methanogen_archaea_by_core_r232.csv (l.434)", "guilds_r232.csv (l.435-437)", "PNGs (l.637-676)", "r232_review_edits_2026-10-26.json (l.44)"],
    outputs=["Metabolic Manuscript - R232 review (tracked changes).docx in the repo root (OUT l.42)"], runtime="not run",
    swap_point="R l.41 (input folder) and OUT l.42 (else it overwrites the repo docx)", one_line="patch R and OUT", wrapper="not wrapped",
    hard_coded_breakpoints=["Table 4 last column and TMA members typed in (l.605-606, 612)", "Table 5 row labels carry counts 22/27/9/15 (l.622-624)",
                            "Results prose literals: 'Twenty-two MAGs carried phnJ ... 11 of the 12 Roseovarius', 'Nine MAGs carried grdH', 'five of the 27' (l.488-494), '254 of the 310', '24-46%' (l.507-509), '99% of random 22-MAG sets' (l.511, 672)",
                            "Figure 8 caption guild taxa (l.652-654)", "edits JSON items 14-16, 29-32, 36 (e.g. 'Three of the 20 (DAINOF01.1, Thiolinea.1, CAZRZN01.1) lack phnJ', 'the nine MAGs carrying grdH')"],
    verified="not run", manuscript=["the manuscript itself"])
add(id="table_S5", script="(no generator in the repo; rebuilt by " + f"{SCR}/make_table_s5.py)",
    annotation_dependent="yes (guild flags + phnJ/grdH/mcrA/cntA/cutC-cluster annotation columns)",
    inputs=["guilds_r232.csv", "mag_methane_correlations_r232.csv", "data/Saltern_phylogeny_r232.json", "MAG_QC/gtdbtk_local_r232/gtdbtk.{ar53,bac120}.summary.tsv", "data/Cliff_310MAG_relabund.tsv"],
    outputs=["Table S5 - MAG taxonomy, quality and guilds (GTDB R232).xlsx (repo root; created by openpyxl 2026-10-21 17:15:52)"], runtime="measured < 5 s",
    swap_point="--run DIR pointing at an alt run folder", one_line=f"~/Documents/py_venv/bin/python -B {SCR}/make_table_s5.py --run OUT --out OUT/TableS5.xlsx",
    wrapper="make_table_s5.py", hard_coded_breakpoints=[],
    verified="rebuilt from the committed inputs: all 32 columns equal the committed workbook (floats within 1e-12), same row order (archaea first, then phylum, then iterativeID)",
    manuscript=["Table S5"])
add(id="gene_annotation_presence", script="clauded_insights/gene_annotation_presence.py",
    annotation_dependent="yes: substring PATTERNS (l.32-45) on the lower-cased, concatenated output of the truncating 'functions' regex (l.56-62) over all 450 RAST files",
    inputs=["RAST dir (l.26, glob *.json l.58)", "data/Cliff_310MAG_relabund.tsv (l.66)", "clauded_insights/read_abundance_by_reaction_totMAGs.csv (l.79)"],
    outputs=["reaction_presence_model_vs_geneannot.csv, read_abundance_by_reaction_totMAGs_geneannot.csv, read_abundance_by_reaction_totMAGs_byType_geneannot.csv (clauded_insights)"],
    runtime="not run", swap_point="RAST l.26 (writes into clauded_insights via HERE l.25)", one_line="patch RAST and HERE", wrapper="not wrapped",
    hard_coded_breakpoints=["norm() joins names (l.47-53); no collisions among the 450 files (checked)"], verified="not run",
    manuscript=["none found (not referenced by the review builder, edits JSON, review notes or comments review)"])
add(id="read_abundance_by_reaction", script="clauded_insights/read_abundance_by_reaction.py (+ reaction_abundance_by_sampletype.py, reaction_abundance_arbitration.md)",
    annotation_dependent="yes: same substring approach (RAST l.63, PATTERNS l.64-77, parser l.78-84)",
    inputs=["RAST dir", "data/Cliff_310MAG_relabund.tsv", "~/Documents/ModelSEEDDatabase/Objects/Reactions.json.gz (l.27)"],
    outputs=["clauded_insights/read_abundance_by_reaction.csv (l.127), read_abundance_by_reaction_totMAGs.csv (l.134) -> reaction_abundance_by_sampletype.py -> _byType.csv, reaction_abundance_vs_sheet2.csv"],
    runtime="not run", swap_point="RAST l.63; outputs hard-wired to REPO/clauded_insights (l.127, 134)", one_line="patch RAST and the two output paths", wrapper="not wrapped",
    hard_coded_breakpoints=[], verified="not run", manuscript=["none found"])
add(id="annotation_independent", script="MAG_QC/compare_vs_2020.py; correlational/make_iterativeIDs_r232.py; correlational/louvain_cooccurrence_network.py; within-condition/GLasso/maxT; 4-site heatmap statistics; Table 3 (make_r232_review_docx.table3 l.396-418, table3_2020_corrected.py, abundances.ipynb)",
    annotation_dependent="no (CheckM, GTDB-Tk, taxonomy labels, abundances, methane)", inputs=[], outputs=["MAG_QC/comparison_vs_2020.csv", "iterativeIDs_r232.json", "network files", "Table 3"],
    runtime="", swap_point="none", one_line="", wrapper="", hard_coded_breakpoints=[], verified="code read only", manuscript=["Table 3", "Figure 8 image", "MAG QC/taxonomy text"])
add(id="fba_2020_models", script="correlational/analysis(1).ipynb @22f4629, analysis(2).ipynb, analysis_2_marimo.py, betaine_gapfilling.ipynb; models in correlational/models/ and models/MePn_consumers/",
    annotation_dependent="yes through the RAST-derived ModelSEED models (GPRs use RAST CDS ids), but frozen: re-running with augmented annotations means rebuilding/gapfilling the models, not swapping an input",
    inputs=["RAST genomes -> MSBuilder models"], outputs=["saltern_fluxes.json, saltern_fluxes_betaine.json, saltern_fluxes_addition.json, fits*.json"],
    runtime="not assessed", swap_point="model rebuild", one_line="", wrapper="", hard_coded_breakpoints=["the 20 MePn / 15 betaine model sets are fixed lists of files"],
    verified="not run", manuscript=["Figure 7 inputs (stored fluxes)"])

RUN_ORDER = [
    "1 guilds (rerun_correlational_r232 build_guilds via the wrapper; use --parser json for augmented sources)",
    "2 corr (guild_correlations, per_mag_correlations, plot_guilds)",
    "3 if betaine_reducer changed: betaine_grdH_r232/rerun_betaine_grdH_r232.py (GUILDS l.65, HERE l.57) then betaine_stats_grdH_r232.py; "
    "otherwise modelling_numbers_r232.py keeps the old grdH flags from fva_tma_by_core_r232.csv (written l.319, filtered l.219) and its assert l.242 fails",
    "4 membership (make_guild_membership_comparison.py), modelling_consistent_basis_r232.py -> figure7_r232.py -> modelling_numbers_r232.py",
    "5 figureS11_r232.py, make_archaeal_mepn_tma_report.py (targets only), make_table_s5.py (scratch), prefilter_sensitivity_* (labels only)",
    "6 make_r232_review_docx.py with R (l.41) and OUT (l.42) patched; then hand-edit the literal prose/Table 4 cells and the edits JSON",
]
CAVEATS = [
    "make_guild_membership_comparison.py writes its CSV (l.261) before the count assert (l.265): a failed augmented run leaves a CSV but no MD",
    "with --parser trunc, strings added via --add-functions bypass the ']' truncation applied to original strings; use --parser json for augmented runs (identical guild columns on the original data)",
    "the guild regexes were checked only against the 2020 RAST/SEED vocabulary (parsing/guild_pattern_matched_strings.csv) and BV-BRC product facets (bvbrc_vocab/); a DRAM/KO-based augmentation would need new rules",
]
json.dump({"analyses": A, "run_order": RUN_ORDER, "caveats": CAVEATS}, open(OUT / "inventory.json", "w"), indent=1)
pd.DataFrame([{k: ("; ".join(v) if isinstance(v, list) else v) for k, v in a.items()} for a in A]).to_csv(OUT / "inventory.csv", index=False)
print("wrote", OUT / "inventory.json", "and inventory.csv:", len(A), "analyses")
