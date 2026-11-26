# Data stage: R220 -> R232 taxonomy and MAG iterativeIDs in ~/Documents/EmilyKin (data/PLAN.md steps A and B)

Started 2026-11-11 16:41 CST. Repo `~/Documents/EmilyKin`, branch `iterativeid-nomenclature`, HEAD 9e087c3.
Nothing is committed or pushed by this stage (the orchestrator commits). Apply scripts live in
`rename_plan/apply_scripts/data/`; scratch outputs in `rename_plan/apply_scratch/data/`.

## 0. Pre-flight

- `git status --porcelain`: 0 lines; HEAD 9e087c3.
- `apply_log/asv_md5_before.tsv`: 1,443 of 1,443 files match before any change.
- `free -g` before the skani run: 251 GB total, 137 GB available.

## 1. D2 network freeze (done FIRST, before the new map was written)

Script: `apply_scripts/data/d2_freeze.py` (asserts HEAD 9e087c3 and an unmodified map and guild table).

| file (new, LF) | content | rows |
|---|---|---|
| `network/mag_linked_asvs_r220_frozen.tsv` | `asv_iterativeID`, `MAG`: the network MAG link at 9e087c3 = abundances.csv columns (935) whose ID is a value of `mag_iterativeID_old_to_new.json` | 140 |
| `network/mag_guild_iterative_ID_stale_frozen.tsv` | `MAG`, `stale_iterative_ID`: the guild table's `iterative_ID` column at 9e087c3 (6a04588-era map), used only by the root `spieceasi_modules_flanking.py` stale-keying audit, because D4 rewrites that column | 291 |

Checks (all passed):
- 140 pairs; for all 140, `asv_to_mag_mapping.json` `best_mag` == the bin (critic C1).
- Every node universe a consumer can see is a subset of the 935 abundances.csv columns (aggregate 935, all-phases 140,
  phases I-V 371/128/870/284/271, the 15 SPIEC-EASI `out/asv/*.csv` with 89 / 53 / 52 nodes), so membership in the
  frozen 140 equals membership in the HEAD map for every node any consumer evaluates.

Consumers changed to read the frozen list (code only; names stay ASV IDs; the MAG's current ID is looked up by bin):

| script | before | after | added output fields |
|---|---|---|---|
| `scripts/render_cooccurrence_grid.py` | `load_abundances`: values of the map | new `MAG_LINK_FROZEN` + `load_mag_linked_asvs()`; `load_abundances` uses it (same column order) | none |
| `scripts/render_cooccurrence_grid_all_phases.py` | values of the map | `rcg.load_mag_linked_asvs()` | none |
| `scripts/build_module_mag_membership.py` | `iid2mag` = inverted map | `iid2mag` = frozen list; `_bin2iid` = map, labels only | each member gets `mag_iterativeID` (current MAG ID); `_description` / `name_map` text rewritten to state the frozen rule |
| `scripts/render_module_denitrification_genes.py` | `iid2mag` = inverted map (critic C1 silent guard) | frozen list | none (not re-run: gate verdict "do NOT use") |
| `scripts/render_module_pathway_completeness.py` | same | frozen list; docstring | none (not re-run: gate verdict) |
| `network/compositional_networks/render_spieceasi_mag_network.py` | `MAG_IDS`/`NEW_TO_BIN` from the map | frozen list; `BIN_TO_MAG_ID` = map | `spieceasi_mag_nodes.csv`: `MAG_iterativeID` |
| `network/compositional_networks/spieceasi_modules_flanking.py` | same | frozen list; `mag_iid()` | edge lists: `source_MAG_iterativeID`, `target_MAG_iterativeID`; `spieceasi_module_membership.csv`: `MAG_iterativeID` |
| `spieceasi_modules_flanking.py` (root, writes `tablesforpaper/`) | `id_to_bin` = inverted map; stale audit read the guild table's `iterative_ID` | `id_to_bin` = frozen list; stale audit reads `mag_guild_iterative_ID_stale_frozen.tsv`; docstring section "GTDB R232 rename" | `modules_membership_89node.csv`: `MAG_iterativeID`; `guild_flanking_*.csv`: `guild_node_MAG_iterativeID`; `pao_gao_direct_edges.csv`: `source_MAG_iterativeID`, `target_MAG_iterativeID` |

`render_gao_pao_networks.py` and `render_cooccurrence_grid_per_phase.py` need no edit: they reach the MAG filter only
through `rcg.load_abundances()` (and the gate showed the committed GAO/PAO figures are the default all-ASV mode, C14).

D2 acceptance test (run in the main repo after the new map was written; logs in `apply_scratch/data/d2_accept/`,
programmatic comparison in `acceptance.txt`): **PASSED**.
- `render_cooccurrence_grid_all_phases.py`: both JSONs and both PNGs byte-identical to HEAD (140 nodes).
- `render_spieceasi_mag_network.py`: PNGs identical; `spieceasi_mag_nodes.csv` = HEAD + column `MAG_iterativeID`.
- compositional `spieceasi_modules_flanking.py`: 4 edge lists and the membership CSV = HEAD + the added columns.
  `spieceasi_flanking_by_module.csv` and `spieceasi_modules_summary.json` differed only in the order of tied members
  (PYTHONHASHSEED, gate C6i); identical as sets and they carry no MAG field, so both were restored from HEAD.
- root `spieceasi_modules_flanking.py`: `modules_membership_89node.csv`, `guild_flanking_all.csv`, `guild_flanking_nonMAG.csv`,
  `pao_gao_direct_edges.csv` = HEAD + the added columns; `pao_gao_shared_driver_summary.json` = HEAD + 28 added
  lines (the 7 direct edges' `source/target_MAG_iterativeID`); the other 5 tables byte-identical;
  `module_by_carbon_source.csv` differed in one float's last digit only (gate tier B) and was restored from HEAD.
- `build_module_mag_membership.py`: `louvain_module_mag_membership.json` = HEAD after removing the added
  `mag_iterativeID` fields, except `_description` / `name_map` (rewritten on purpose) and 6 `gamma_rule` strings
  copied from the current network summaries (gate C6e: the committed file had older wording).
- `render_module_flanking_figures.py`: PNGs identical. The 4 untracked `figures/spieceasi_*.pdf` were rewritten.
- Not re-run (gate verdicts): module panels C6f/C6g, GAO/PAO networks (default mode never reads the map), aggregate and
  per-phase networks (all-ASV; no MAG fields), `render_spieceasi_docx.py`.

## 2. A1-A5: GTDB R232 taxonomy (by bin)

**A1** `gtdbtk.bac120.summary.tsv` <- prep `results/gtdbtk.r232.summary.tsv` (md5 a4b37ad9…), converted LF -> CRLF so the
repo file keeps its CRLF terminator (rule: preserve line endings). Same 291 bins, same row order, same 20 columns.
`git diff --numstat` 279/279 lines (130 lineages; ANI/AF/RED/reference fields of the re-run; the best-of-both
hybrid is gone, as the plan intends).

**A2-A5 and B3**: new repo script `scripts/apply_r232_taxonomy.py` (reads the R232 lineage from the repo's own
`gtdbtk.bac120.summary.tsv`, joins on the bin, copies rank names verbatim, `s__` empty -> ""). Its writer is the new
`scripts/_binkeyed.py` (`patch_table`: unchanged lines byte-identical, CRLF/LF kept per file, refuses mixed EOL /
multi-line quoted fields; every target was round-trip tested first: 31 files, 0 lines not reproduced). Idempotent
(second run: 0 changes). Run logs: `apply_scratch/data/apply_r232_taxonomy.run{1,2}.log`.

| step | file | rows changed / rows | cells changed | EOL |
|---|---|---|---|---|
| A2 | `mag_quality_summary.tsv` | 130 / 291 | Phylum 22, Class 23, Order 19, Family 24, Genus 60, Species 43 | CRLF |
| A3 | `mag_abundance_summary.tsv` | 130 / 291 | same as A2 (`iterativeIDs` ASV column untouched) | CRLF |
| A3 | `mag_abundance_summary.json` | 291 records | 148 taxonomy cells (6-rank schema kept; `json.dumps` form verified before writing) | LF |
| A4 | `leftover_reads_analysis/mag_taxonomy.json` | 291 | 148 cells (blank rank stays `null`, as the file stores it; `indent=0` form verified) | LF |
| A4 | `leftover_reads_analysis/mag_genera.json` | — | 172 -> 183 genera (-26, +37) | LF |
| A5 | `mimag_hq_classification.tsv` | 98 / 291 | Phylum 22, Genus 60, Species 43 | CRLF |

The counts include rank gains/losses (Family +1/-1, Genus +11/-1, Species +26/-1); the plan's 22/49/17 count renames
only. Check: `mag_quality_summary.tsv` equals `data/mag_lineage_r220_r232.tsv` r232_* for 291 x 7 cells (0 mismatches).

**C25 leftover reads.** `aggregate_report.py` (gate tier A) re-run after A4: all 5 outputs byte-identical (no
orphan-ASV `genus_recovered` flag changes). `parse_paf_16s.py` not runnable (PAF from a pipe, never stored), so
`out/CAN_*.ont.json` and the labels in `leftover_taxa.csv` / `REPORT.md` keep R220 genus names (Holophaga, JACRCT01,
JAIEMT01, UBA5200, UBA919). Frozen with a note: new file `leftover_reads_analysis/R232_NOTE.md`.

## 3. B1 + D1: the R232 map

- Backups of the three map files: `apply_scratch/data/backups/`.
- `scripts/recompute_mag_ids.py` run after A3 (log `recompute_mag_ids.r232.log`): cases 138/70/83, 0 duplicates,
  0 numeric-root `_m`, 0 inherited IDs that are not ASV IDs, 51 changed (50 at the root).
- **GATE (b1_gate.log): PASSED.** New json == `data/mag_id_map_r220_r232.json` 291/291; == simulation
  `results/iterativeid_simulation_r232.tsv` 291/291; `assignment_case` == `case_r232` 291/291; `_new_ids.json`
  mag_new/mag_case consistent 291/291 and its 3 ASV keys identical to HEAD; pre-c796fdf "old" CSV columns preserved
  291/291. numstat: csv 51/51 (CRLF kept), json 1/1, `_new_ids.json` 1/1.
- **D1 (new, LF):** `mag_iterativeID_r220_to_r232.csv` (header `bin,old_id,new_id,case_old,case_new`, 291 rows, 51
  changed) and `mag_iterativeID_r220_to_r232.json` (`{bin: {old_id, new_id, case_old, case_new}}`). old = HEAD map/case,
  new = gated R232 map/case.

## 4. B2: apply_mag_ids extension (EOL-preserving)

`scripts/apply_mag_ids.py`: `set_mag_col` now calls `_binkeyed.patch_table` (C11 fix: the old `csv.writer` wrote CRLF into
every file); docstring records D4 and the exclusions. Added the 12 non-r13 ADD targets of `apply_mag_ids_extension.tsv`
(fig1_annotation.tsv: `label` and `iterativeID` of MAG tips only). **Not added:** the 3 r13 ADD rows (`taxonomy.tsv`,
`panel/mag_panel.tsv`, `rast/genomes/_manifest.tsv`): B9 waits for the KBase stage. Log `apply_mag_ids.r232.log`;
re-run idempotent (0 changes).

| file | ID cells changed (rows) | not in map | EOL |
|---|---|---|---|
| `mag_abundance_summary.tsv` | 51 | 0 | CRLF |
| `mag_abundance_by_day_intersection_with_ids.csv` | 51 | 0 | CRLF |
| `gene_ab_figure/data/taxonomy_labels.tsv`, `offline_bundle/…` | 2 each | 0 | CRLF |
| `gene_ab_figure/data/gene_copy_per_mag.tsv` | 2 | TOTAL row | CRLF |
| `ko_gene_abundance/ko_rpkm_combined.tsv` | 50 | CAN_4_bin.225 (276-era) | CRLF |
| `ko_gene_abundance/ko_rpkm_per_sample.tsv` | 2,100 | 42 rows, CAN_4_bin.225 | CRLF |
| `ko_copy_number_matrix_nosZ.csv` | 51 | 0 | CRLF |
| `clade_classify/out/allbins_nosz_master.tsv` | 27 | CAN_4_bin.225 | CRLF |
| `Table_S7_MAG_functional_classification.csv` (ADD) | 82 (51 R232 + 31 stale, D4) | 0 | LF |
| `mag_guild_classification.csv` (ADD) | 82 | 0 | LF |
| `mag_percell_abundance.csv` (ADD) | 51 | 0 | LF |
| `ko_gene_abundance/ko_rpkm_combined_ext.tsv` (ADD) | 86 | CAN_4_bin.225 | CRLF |
| `ko_gene_abundance/ko_rpkm_per_sample_ext.tsv` (ADD) | 3,870 | 45 rows | CRLF |
| `gene_ab_figure/data/rpkm_long_MAG_gene_sample_phase.tsv` (ADD) | 19,350 | 225 rows | LF |
| `gene_ab_figure/data/rpkm_long_MAG_gene_sample_291_preview.tsv` (ADD) | 7,135; 1,115 never-filled rows (blank ID and genus) left blank | 0 | CRLF |
| `gene_ab_figure/data/novel_bacteria_summary.tsv` (ADD) | 1 (stale Dokdonella.14 -> Dokdonella.8, D4) | 0 | CRLF |
| `fig1_tree/fig1_annotation.tsv` (ADD) | 16 (8 MAG tips x label+ID) | 61 reference tips | CRLF |
| `fig1_tree/fig1_mags_table.tsv` (ADD) | 8 | 0 | LF |
| `phase_rpkm_chord/mag_category_assignments.csv` (ADD) | 10 (incl. pre-c796fdf `Root_m.N`, C19) | 0 | CRLF |
| `phase_rpkm_chord/phase_mag_gene_rpkm_long.csv` (ADD) | 185 | 0 | CRLF |
| `mag_abundance_summary.json` | 291 set (51 differ) | — | LF |

**C18 / D1 gate (`b2_gate.log`): PASSED** for all 22 ID columns: every row whose bin is in the map carries the R232 ID
(0 wrong), one ID per bin, no ID shared by two bins; the only blanks are the 1,115 never-filled preview rows. Every
`git diff --numstat` of these files has added == deleted lines (no EOL rewrite).

C12: `ko_rpkm_combined_hq.tsv` / `ko_rpkm_per_sample_hq.tsv` have ID columns that are blank in every row (ffe668c era);
**left blank** (consumers use the map). Only their genus was patched (B3).

## 5. B3: genus / species of derived MAG tables (same script, `apply_r232_taxonomy.py`)

Every cell equalled the R220 value before the patch (profiled first; blank-genus preview rows: see below).

| file | rows changed | cells |
|---|---|---|
| `ko_gene_abundance/ko_rpkm_combined.tsv` / `_ext` / `_hq` | 57 / 57 / 35 | genus |
| `ko_gene_abundance/ko_rpkm_per_sample.tsv` / `_ext` / `_hq` | 2,394 / 2,565 / 591 | genus |
| `gene_ab_figure/data/rpkm_long_MAG_gene_sample_phase.tsv` / `_phase_hq` | 12,825 / 3,469 | genus |
| `gene_ab_figure/data/rpkm_long_MAG_gene_sample_291_preview.tsv` | 4,315 (incl. 480 genus gains in filled rows); 1,115 never-filled rows left blank | genus |
| `gene_ab_figure/data/taxonomy_labels.tsv` and `offline_bundle/…` | 5 each | Genus 2, Species 5 |
| `gene_ab_figure/data/gene_copy_per_mag.tsv` | 2 | genus |
| `gene_ab_figure/data/nosz_clades.tsv`, `nosz_cladeII_subclades.tsv`, `clade_classify/out/nosz_clades_updated.tsv` | 1 each | genus |
| `clade_classify/out/allbins_nosz_master.tsv` | 22 | genus |
| `gene_ab_figure/data/novel_bacteria_summary.tsv` | 0 (phylum..genus unchanged in R232) | — |
| `phase_rpkm_chord/mag_category_assignments.csv` / `phase_mag_gene_rpkm_long.csv` | 2 / 50 | Genus |

Validation: `gene_ab_figure/08_gene_copy_per_mag.py` (gate: .tsv tier A) re-run afterwards reproduces the bin-patched
`gene_copy_per_mag.tsv` byte for byte; its `.md` (stale at HEAD) now carries the R232 IDs/genera (6 rows; also fixes
the D4 drift Dokdonella.14 -> .8, Ca_Competibacter.139 -> .85).

## 6. B4-B8: regenerations

- **B4 `04_root_aggregate.py`: NOT run** (gate: tier C; committed root tables are 924fa18-era, 150 vs 146 roots; rows
  keyed by root, so no bin-keyed patch is possible). `ko_rpkm_root_{combined,per_sample}.tsv` left as they are
  (stale before R232). Needs the orchestrator's approval to re-run.
- **B5 `scripts/build_percell_abundance.py`** (gate tier A) re-run: `mag_percell_abundance.csv` differs from HEAD only in
  iterativeID 51, Phylum 22, Class 23, Order 19, Family 24, Genus 60, Species 43; all numeric columns identical;
  "REGRESSION PASSED" vs the submission workbook; deterministic (2 runs). numstat 131/131 (LF kept).
- **B6 builder workbook `MAG_target_gene_filled.xlsx` only** (needed by B7 and written by B8's S6 sheet; `build_xlsx_291.py`
  NOT used, gate tier C). Script `apply_scripts/data/patch_builder_workbook.py` (openpyxl, formulas kept; HEAD copy in
  `apply_scratch/data/backups/MAG_target_gene_filled.HEAD.xlsx`); change list `apply_scratch/data/builder_workbook_changes.tsv`
  (1,052 cells). S1 and S1b: ID 82, Phylum 22, Class 23, Order 19, Family 24, Genus 60, `Taxonomy ` 130 each; S2a ID 54;
  S3 ID 82 + Taxonomy 130; S4 ID 42 + genus 22; S5 ID 2. Cross-check vs docs `table_cells_plan.tsv` for this file
  (S1-S5): 1,052 / 1,052 cells agree (the 4 "(blank in R232)" cells written as empty). **The submission workbook
  `MAG_target_gene_filled_UPDATED (1).xlsx`, `cladeexcel.xlsx`, `abundance and coverage (1).xlsx` were NOT touched
  (docs stage).** Zip parts unchanged; S2b and KO_order_mapping unchanged.
- **B8a S6 sheet** `scripts/build_table_s6_ani.py` after the C7 fix (`GTDB_RELEASE = "R232"`; `REF_DB` "GTDB R232
  (GTDB-Tk v2.8.0, classify_wf)"; the skani note and `--skani` help say R232; docstring no longer says the summary
  lives on the cluster). Run with no arguments (auto-loads the repo R232 summary) and **without `--skani`** (method
  parity with HEAD, gate). 291 local-filled, 291 GTDB-Tk-filled, 0 via skani, 169 unresolved (HEAD 194). Calls:
  121 same species (HEAD 96), 1 same species AF<0.60, 169 unresolved. Changed cells: Reference database 291, Notes 137,
  ID 82, GTDB taxonomy 130, ranks 22/23/19/24/60/43, reference name 43, accession 36, ANI 71, AF 74, call 27,
  interpretation 40. The 276-era `.pre_s6.bak.xlsx` (gitignored) was left as it is.
- **B8b D3 skani search** (`apply_scratch/data/skani_r232/run_skani_r232.sh`, log `skani_search.log`): skani 0.4.0
  `search --min-af 0 -n 3 -t 24` vs `~/Documents/gtdbtk_data/release232/skani/database` (199,923 reps), nice 10,
  systemd scope MemoryMax 60G; free memory checked first (137 GB available); 3 min 38 s wall, 41.8 GB max RSS;
  727 hits for 268 MAGs (23 with no hit). Committed to the repo with its run script:
  `reprocess_291/skani_novelty_r232/{mag_vs_gtdb_r232_skani.tsv,run_skani_r232.sh}` (328 KB, not gitignored).
- **B8c standalone `Table_S6_skani_novelty.tsv`**: `reprocess_291/build_skani_novelty_291.py` given environment
  overrides (`SKANI_TSV`, `GTDB_TAX`, `GTDBTK_SUMMARY`, `MAG_DIR`, `OUT`; defaults unchanged = the poplar R220 paths)
  and gzip support for the taxonomy; logic unchanged (gate logic check 291/291). Reference taxonomy =
  `~/Documents/gtdbtk_data/release232/taxonomy/gtdb_taxonomy.tsv` (the GTDB-Tk release file, as the R220 run used its
  release220 copy). Result (CRLF kept, 165/165 lines): **>=95: 121, 85-95: 76, <85: 28 (23 no hit), low_AF_unresolved:
  66** (R220: 96 / 70 / 48 / 77). Checks: `mag_gtdb_taxonomy` == repo R232 summary 291/291; every best-hit reference
  taxonomy == `gtdb_downloads/r232/bac120_taxonomy_r232.tsv.gz` (268/268). coasm_bin.476 -> GCF_032229525.1
  *Accumulibacter* sp032229525, 98.68 % ANI, AF_query 72.88 % (>=95; R220: low_AF_unresolved).
  **Recipe artifact to flag:** coasm_bin.535 is `low_AF_unresolved` because `-n 3` keeps the three highest-ANI hits,
  all tiny-AF (100 % ANI on 0.44 % AF); GTDB-Tk R232 assigns it *Novosphingobium* sp017302615 (98.4 % ANI, AF 0.90).
  It is the only MAG with a GTDB-Tk species that the recipe does not call >=95.
- **B7 `scripts/build_table_s7_classification.py`** after the D6 fix, reading the patched S1/S3: Table S7 CSV differs from
  HEAD only in iterative_ID 82, Phylum 22, Family 24, Genus 60; `mag_guild_classification.csv` only in iterative_ID 82,
  Genus 60; every classification column (tiers, guild, storage and denitrification categories) is identical, so D6 and
  R232 change no call. S7 sheet: the same 4 columns. LF kept (124/124, 103/103).

## 7. Code fixes (`code_fixes.tsv`, critic C7/C8, D6)

| file | change |
|---|---|
| `scripts/render_fig1_tree.py` | root on Phylum in {Patescibacteria, Patescibacteriota}; `GTDB2PALETTE` + Patescibacteriota->Patescibacteria, Myxococcota_A->Myxococcota; **D3:** reference tips relabelled from the GTDB-Tk release232 `taxonomy/gtdb_taxonomy.tsv` (new `--ref-taxonomy`, `none` = old behaviour); GB_GCA_017302355.1 (absent from R232) keeps its R220 lineage and is labelled "GCA_017302355.1 (not in GTDB R232)"; title "… GTDB r220 reference genomes (GTDB R232 names)"; docstring. The `how` string is unchanged. |
| `fig1_tree/compare_trees.py` | same rooting test (both names); help text |
| `community_modeling/analysis/r13/annotations/build_taxonomy.py` | `DIDERM` + Myxococcota_A; `FLAG['Patescibacteriota'] = FLAG['Patescibacteria']` (not run: B9) |
| `…/panel/tests/test_build_panel.py` | fixture `p__Patescibacteria` -> `p__Patescibacteriota`; `pytest -p no:cacheprovider panel/tests`: 15 passed |
| `scripts/render_percell_composition.py` | `DOCX_COLORS` keys follow the same MAGs with the same colours: SCMT01->Lysobacter, Azospira->Proximibacter, Giesbergeria->Simplicispira; docstring. Checked with the module's own `load()/build()`: the 20 bands are the same entities (3 renamed, Unclassified reordered), 0 fallback colours; JJ008 keeps CAN_5_bin.112 (coasm_bin.492 -> CADEDV01, pooled in Others). OLB5's label becomes "(Ignavibacteria (Bacteroidota))" because its class is now Ignavibacteria. Not rendered (C5). |
| `scripts/build_table_s7_classification.py` | D6: Tier-1 GAO regex `Competibacter|Contendobacter|Contendibacter`; docstring note |
| `scripts/build_table_s6_ani.py` | C7 (above) |
| `reprocess_291/build_skani_novelty_291.py` | env overrides + gzip taxonomy (above) |
| `clade_classify/render_nosz_top5_tree.py` | docstring: CAN_3_bin.221 / Tepidaquicellales.1_m (Ignavibacteria.21 before R232) |
| `fig1_tree/score_ref_candidates.py`, `scripts/build_fig1_alignment.py`, `fig1_tree/run_fig1_raxml.sh` | comments: g__Azospira_A (r220) = g__Proximibacter (R232), Azospira_A restricta = Proximibacter restrictus; **hazard notes**: both scripts read the repo summary / Table S6 (now R232) but match against r220 references, so re-running them now would not reproduce the committed Figure-1 alignment (use the 9e087c3 R220 files) |
| left on purpose | `calls/call_genes.py:6`, `panel/validate_panel.py:80` ("five Patescibacteria MAGs": optional wording, logic by bin); `reprocess_291/run_tree62.sh` (historical server script); the `leave` rows of `code_fixes.tsv` |

Test render of Fig 1 to scratch only (`apply_scratch/data/fig1_test/`): "root: Patescibacteria clade (7 tips)" (gate
string), labels unique, 9 reference labels change (8 R232 names + the GCA_017302355.1 note), 3 reference phyla.
`compare_trees.py` with the R232 outgroup: both trees rooted on "Patescibacteria clade (7 tips)". The committed figure
and annotation are left for C1.

## 8. Final checks (`apply_scratch/data/final_checks.txt`)

- `asv_md5_before.tsv`: 1,443 / 1,443 files unchanged (incl. the C13 all-phases network JSONs).
- 72 tracked files modified + 9 new files; **no line-terminator change** in any modified file (HEAD vs now).
- `git diff --numstat`: every data table has added == deleted lines; the unequal ones are code, `mag_genera.json`
  (172 -> 183 entries) and the two JSONs that gained the D2 fields.
- D5: no log, run record, `*.pre_*.bak.*`, payload/verify JSON, r12 snapshot or track2_guards file modified.
  `community_modeling/plan_r13/integrator/AMENDMENTS.sha256` (run from `community_modeling/`): 14/14 OK;
  `track2_guards_sha256_before.txt`: 7/7 OK.
- Both apply scripts re-run on the final tree: working tree unchanged (idempotent).
- Bytecode written by the checks was deleted; no other new ignored files except the 4 regenerated untracked
  `network/compositional_networks/figures/spieceasi_*.pdf`.
- Nothing committed, nothing pushed, nothing sent to KBase.

## 9. Files changed (`git diff --numstat`, added / deleted lines)

| file | + | - |
|---|---|---|
| `MAG_target_gene_filled.xlsx` | - | - |
| `Table_S6_skani_novelty.tsv` | 165 | 165 |
| `Table_S7_MAG_functional_classification.csv` | 124 | 124 |
| `clade_classify/out/allbins_nosz_master.tsv` | 29 | 29 |
| `clade_classify/out/nosz_clades_updated.tsv` | 1 | 1 |
| `clade_classify/render_nosz_top5_tree.py` | 2 | 1 |
| `community_modeling/analysis/r13/annotations/build_taxonomy.py` | 3 | 1 |
| `community_modeling/analysis/r13/annotations/panel/tests/test_build_panel.py` | 1 | 1 |
| `fig1_tree/compare_trees.py` | 4 | 2 |
| `fig1_tree/fig1_annotation.tsv` | 8 | 8 |
| `fig1_tree/fig1_mags_table.tsv` | 8 | 8 |
| `fig1_tree/run_fig1_raxml.sh` | 1 | 1 |
| `fig1_tree/score_ref_candidates.py` | 4 | 1 |
| `gene_ab_figure/data/gene_copy_per_mag.md` | 6 | 6 |
| `gene_ab_figure/data/gene_copy_per_mag.tsv` | 3 | 3 |
| `gene_ab_figure/data/nosz_cladeII_subclades.tsv` | 1 | 1 |
| `gene_ab_figure/data/nosz_clades.tsv` | 1 | 1 |
| `gene_ab_figure/data/novel_bacteria_summary.tsv` | 1 | 1 |
| `gene_ab_figure/data/rpkm_long_MAG_gene_sample_291_preview.tsv` | 8020 | 8020 |
| `gene_ab_figure/data/rpkm_long_MAG_gene_sample_phase.tsv` | 22950 | 22950 |
| `gene_ab_figure/data/rpkm_long_MAG_gene_sample_phase_hq.tsv` | 3469 | 3469 |
| `gene_ab_figure/data/taxonomy_labels.tsv` | 6 | 6 |
| `gene_ab_figure/offline_bundle/data/taxonomy_labels.tsv` | 6 | 6 |
| `gtdbtk.bac120.summary.tsv` | 279 | 279 |
| `ko_copy_number_matrix_nosZ.csv` | 51 | 51 |
| `ko_gene_abundance/ko_rpkm_combined.tsv` | 70 | 70 |
| `ko_gene_abundance/ko_rpkm_combined_ext.tsv` | 102 | 102 |
| `ko_gene_abundance/ko_rpkm_combined_hq.tsv` | 35 | 35 |
| `ko_gene_abundance/ko_rpkm_per_sample.tsv` | 2940 | 2940 |
| `ko_gene_abundance/ko_rpkm_per_sample_ext.tsv` | 4590 | 4590 |
| `ko_gene_abundance/ko_rpkm_per_sample_hq.tsv` | 591 | 591 |
| `leftover_reads_analysis/mag_genera.json` | 37 | 26 |
| `leftover_reads_analysis/mag_taxonomy.json` | 148 | 148 |
| `mag_abundance_by_day_intersection_with_ids.csv` | 51 | 51 |
| `mag_abundance_summary.json` | 1 | 1 |
| `mag_abundance_summary.tsv` | 131 | 131 |
| `mag_guild_classification.csv` | 103 | 103 |
| `mag_iterativeID_old_to_new.csv` | 51 | 51 |
| `mag_iterativeID_old_to_new.json` | 1 | 1 |
| `mag_percell_abundance.csv` | 131 | 131 |
| `mag_quality_summary.tsv` | 130 | 130 |
| `mimag_hq_classification.tsv` | 98 | 98 |
| `network/compositional_networks/out/asv/spieceasi_glasso_edgelist.csv` | 718 | 718 |
| `network/compositional_networks/out/asv/spieceasi_glasso_edgelist_MAGonly.csv` | 267 | 267 |
| `network/compositional_networks/out/asv/spieceasi_mag_nodes.csv` | 54 | 54 |
| `network/compositional_networks/out/asv/spieceasi_mb_edgelist.csv` | 402 | 402 |
| `network/compositional_networks/out/asv/spieceasi_mb_edgelist_MAGonly.csv` | 160 | 160 |
| `network/compositional_networks/out/asv/spieceasi_module_membership.csv` | 90 | 90 |
| `network/compositional_networks/render_spieceasi_mag_network.py` | 12 | 7 |
| `network/compositional_networks/spieceasi_modules_flanking.py` | 15 | 2 |
| `network/louvain_module_mag_membership.json` | 1290 | 649 |
| `phase_rpkm_chord/mag_category_assignments.csv` | 11 | 11 |
| `phase_rpkm_chord/phase_mag_gene_rpkm_long.csv` | 195 | 195 |
| `reprocess_291/build_skani_novelty_291.py` | 13 | 7 |
| `scripts/_new_ids.json` | 1 | 1 |
| `scripts/apply_mag_ids.py` | 36 | 20 |
| `scripts/build_fig1_alignment.py` | 5 | 1 |
| `scripts/build_module_mag_membership.py` | 22 | 13 |
| `scripts/build_table_s6_ani.py` | 6 | 6 |
| `scripts/build_table_s7_classification.py` | 3 | 1 |
| `scripts/render_cooccurrence_grid.py` | 17 | 3 |
| `scripts/render_cooccurrence_grid_all_phases.py` | 3 | 2 |
| `scripts/render_fig1_tree.py` | 43 | 12 |
| `scripts/render_module_denitrification_genes.py` | 4 | 1 |
| `scripts/render_module_pathway_completeness.py` | 7 | 3 |
| `scripts/render_percell_composition.py` | 9 | 5 |
| `spieceasi_modules_flanking.py` | 32 | 10 |
| `tablesforpaper/guild_flanking_all.csv` | 872 | 872 |
| `tablesforpaper/guild_flanking_nonMAG.csv` | 340 | 340 |
| `tablesforpaper/modules_membership_89node.csv` | 90 | 90 |
| `tablesforpaper/pao_gao_direct_edges.csv` | 15 | 15 |
| `tablesforpaper/pao_gao_shared_driver_summary.json` | 28 | 0 |

New files: `mag_iterativeID_r220_to_r232.csv`, `mag_iterativeID_r220_to_r232.json`, `network/mag_linked_asvs_r220_frozen.tsv`,
`network/mag_guild_iterative_ID_stale_frozen.tsv`, `scripts/_binkeyed.py`, `scripts/apply_r232_taxonomy.py`,
`reprocess_291/skani_novelty_r232/mag_vs_gtdb_r232_skani.tsv`, `reprocess_291/skani_novelty_r232/run_skani_r232.sh`,
`leftover_reads_analysis/R232_NOTE.md`.

## 10. Not done here / hand-offs

- **B9 (r13 tables)**: waits for the KBase stage (rename log). Code fixes are in. `PANEL_FROZEN.md` note belongs to B9.
- **B4 root tables**: stale, not regenerated (needs approval).
- **Docs stage**: submission + collaborator workbooks (the builder workbook is DONE here: do not patch it again);
  `render_spieceasi_docx.py` text still names `mag_iterativeID_old_to_new.json` as the MAG-link source (the counts,
  53 nodes, are unchanged under D2); `novel_bacteria_summary.{tsv,md}` novelty wording for coasm_bin.260 (gains
  s__JAAEKA01 sp049294015); whether builder S6 should use `--skani` (open decision; parity kept).
- **Figures stage**: C1 (Fig 1 code ready and test-rendered), C2 (08 already re-run here), C3, C4 (patched CSVs; 01 not
  usable), C5 (percell colours fixed), C22 PDFs. `fig1_annotation.tsv` carries R232 IDs but still R220 taxonomy until C1.
- **Scientific changes surfaced by the data (for the D7 list):** Table S6 buckets 121/76/28/66 (was 96/70/48/77);
  builder S6 same-species calls 121 (was 96); coasm_bin.476 no longer novel (sp032229525, 98.68 %); guild / Table S7
  calls unchanged; R232 has 15 phylum names (Myxococcota_A separate); percell legend genera renamed.

## 11. Late checks and housekeeping

- `mag_abundance_summary.json` vs HEAD, per record: only `taxonomy.{Phylum,Class,Order,Family,Genus}` (22/23/19/24/60)
  and `mag_iterativeID` (51) differ; the ASV `iterativeIDs` lists and abundance dicts are untouched. The `.tsv` differs
  only in Phylum..Species and `mag_iterativeID`.
- Readers of the outputs that gained fields (`render_module_flanking_figures.py`, `render_spieceasi_docx.py`, the root
  `spieceasi_modules_flanking.py`): all read by column / key name; none indexes by position or asserts a key set.
- **PANEL_FROZEN.md:** A1 changed `gtdbtk.bac120.summary.tsv`, so the input hash that `panel/PANEL_FROZEN.md` cites for it
  is stale from this stage on (no code checks it). The dated note belongs to B9 (do not edit the freeze-2 hashes).
- Scratch artefacts to keep until the R232 apply is verified (`rename_plan/apply_scratch/data/`): `skani_r232/`
  (raw search, run script, log), `backups/` (HEAD copies of the map files, builder workbook, S7/guild CSVs, percell,
  gene_copy), `fig1_test/` (scratch Fig 1 render + annotation + compare_trees TSV), `d2_accept/` (run logs,
  `acceptance.txt`), `builder_workbook_changes.tsv`, the per-step `*.log` files, `b1_gate.log`, `b2_gate.log`,
  `final_checks.txt`, `numstat_table.md`. Apply scripts: `rename_plan/apply_scripts/data/{d2_freeze.py,patch_builder_workbook.py}`.

Finished 2026-11-12 13:40 CST (the run spanned the date change; it is one apply, not two).
