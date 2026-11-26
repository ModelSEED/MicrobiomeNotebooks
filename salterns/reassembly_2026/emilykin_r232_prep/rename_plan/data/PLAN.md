# EmilyKin R220 → R232: rename plan for tracked data files and scripts (data part)

Prepared 2026-11-11, read-only. Nothing in `~/Documents/EmilyKin` (branch `iterativeid-nomenclature`, clean at 9e087c3)
or in KBase was changed. Everything here is in `rename_plan/data/`. The scripts that produced it are in `scripts/`
(section 10). The manuscript, SI documents and KBase are other parts. This part only lists their files and the
dependencies on them (section 9).

## 1. Bottom line

- **The maps are verified.** The naming rule, re-run here with the R220 taxonomy, reproduces the committed
  `mag_iterativeID_old_to_new.json` for 291 of 291 MAGs. The committed IDs also agree 291/291 with
  `mag_abundance_summary.{json,tsv}`, the intersection CSV and the committed `assignment_case`.
- **What R232 changes.** 51 MAG iterativeIDs change, 50 of them at the root. There are no duplicate IDs.
  130 of the 291 lineages change; 161 are identical.
- **The rename is not a string replacement.** It must be a bin-keyed join:
  - 21 old IDs are the IDs of ASVs that keep them.
  - 4 new IDs already exist as ASV IDs: `Simplicispira.1`, `Holophagaceae.9`, `Taibaiella.5` and `Hoeflea.1`.
  - One pair of IDs is a chain: `Chitinophagaceae.1_m` belongs to coasm_bin.117 today and to CAN_5_bin.149 after
    the rename.
- **Most occurrences of the old ID strings are ASV context.** Of 117,114 occurrences in 1,391 tracked files:
  - 19,314 are MAG context. This covers all 11,743 occurrences of the 30 MAG-only `_m` IDs, plus 7,571 occurrences of
    the 21 inherited IDs that sit in MAG columns.
  - 97,725 are ASV context.
  - 75 are ambiguous.
  - 1,068 of the files are community-modelling run records. Their member names are ASV synthetic genomes and they
    stay as they are.
- **One change affects results, not just labels: the network MAG link.**
  - Several analyses call an ASV node a MAG when its ID string equals a MAG's iterativeID (`NETWORK_MAG_ONLY`):
    the all-phases network, the GAO/PAO networks, the SPIEC-EASI MAG display, Louvain MAG membership, the module
    denitrification and pathway panels, and the `tablesforpaper/` flanking tables.
  - Under R232, 17 ASV nodes lose their MAG and 4 gain one.
  - The MAG-only universe goes from 140 to 127 ASVs, and the SPIEC-EASI MAG display from 53 to 45 nodes.
  - This needs a decision: D2 (section 6).
- **Community modelling needs no model re-run.**
  - Members are the 16S-ASV synthetic genomes (community_modeling/README.md line 4), named by ASV iterativeID.
  - MAG IDs appear only in the r13 annotation tables (keyed by bin) and in KBase RAST object names.
  - What needs updating: the r13 `taxonomy.tsv`, `mag_panel.tsv` and `calls/summary.tsv` (rebuild or patch), and
    two phylum lists in `build_taxonomy.py`.
- **What the task assumed about the taxonomy sources was partly wrong.**
  - `iterativeID_taxonomy.json`, `iterativeID_phylums.json`, `iterativeID_color_map.json` and `iterativeID_levels.json`
    are keyed by the 3,945 ASV IDs and hold MiDAS taxonomy.
  - `Genus_color_map.json`, `Species_color_map.json`, `Phylum_color_map.json` and `phylum_colors.json` are keyed by
    48–49 MiDAS phylum names (`Patescibacteria`, `Proteobacteria`, …).
  - None of these is a MAG source. Leave them unchanged.
  - The MAG figures map GTDB names onto these palettes in code (`render_fig1_tree.py` `GTDB2PALETTE`), and that map
    needs the new phylum names.

## 2. Authoritative maps (built and checked here)

| file | content |
|---|---|
| `mag_id_map_r220_r232.tsv` | all 291 bins: `old_id` (committed R220), `new_id`, `case_old` / `case_new`, `changed`, `root_changed`, `old_is_asv_id`, `new_is_asv_id`, chain columns, R220 / R232 genus |
| `mag_id_map_r220_r232.json` | `{bin: new_id}` for 291 bins (the R232 equivalent of `mag_iterativeID_old_to_new.json` / `_new_ids.json` `mag_new`) |
| `mag_id_changes_r220_r232.csv` | the 51 changes only, with chain notes |
| `mag_lineage_r220_r232.tsv` | all 291: R220 lineage (repo `gtdbtk.bac120.summary.tsv`), R232 lineage (`results/gtdbtk.r232.summary.tsv`), per-rank R220/R232 names, changed ranks, R232 method/ANI/AF/RED, and whether the repo's own R220 copies agree |
| `taxon_name_changes_r220_r232.tsv` | 109 distinct (rank, R220 name → R232 name) changes over the 291 MAGs, and whether the R220 name is also an ASV (MiDAS) taxon name in `taxonomy.json` |
| `substring_hazards.json` | old IDs that are string prefixes of other live IDs (only `Rhodobacteraceae.1` → `.10`…`.19`) |
| `build_maps.log` | the verification output |

**Checks** (`build_maps.log`):
- The simulation, re-run in `build_maps.py`, gives the same IDs as `results/iterativeid_simulation_r232.tsv` for
  291 of 291.
- **Old → new cases** (51): asv_inherit → gtdb_genus 13, asv_inherit → gtdb_lowest 8, gtdb_genus → gtdb_genus 13,
  gtdb_lowest → gtdb_lowest 10, gtdb_genus → asv_inherit 4, gtdb_genus → gtdb_lowest 2, gtdb_lowest → gtdb_genus 1.
- **R220 taxonomy copies.**
  - `mag_quality_summary.tsv` equals the gtdbtk summary for 291 of 291.
  - `mag_abundance_summary.{json,tsv}` match 290/291. coasm_bin.221 has a blank Genus there but `1-14-0-10-56-10` in
    gtdbtk/quality; the R232 genus is CAZOZK01.
  - Being numeric, the genus does not affect naming.
  - The json holds 6 ranks (no Species); the tsv and the quality table hold 7.
- **Taxon-name collisions.** R220 names that change and are *also* MiDAS names of ASVs:
  - `Patescibacteria`, `Myxococcota`, `Verrucomicrobiae`, `ABY1`, `Leptospirae`, `SJA-28`, `Flavobacteriales`,
    `Hydrogenedentiales`, `Burkholderiaceae`, `Bosea`.
  - Never replace these as bare strings. GTDB-prefixed forms (`p__Patescibacteria`) or bin-keyed cells only.
- **Display variants.** The only `Ca. ` display form in the tracked text and code is in
  `network/compositional_networks/render_spieceasi_docx.py:402` (`Ca. Competibacter.32`), and there it names an ASV
  node.

## 3. How every occurrence was classified

Each tracked file was scanned: 2,896 candidates out of 13,425 tracked files, including `.gz` text and the text of
xlsx/docx/pptx/pdf. A boundary-aware tokenizer was used, so `Thauera.3` does not match `Thauera.34`, `X.5` does not match
`X.5_m` or `X.5_ma`, and embedded KBase names (`emilykin_MAGs_RAST_code11_b1_Ca_Falkowbacteria.7`) are resolved. Coverage
was cross-checked against `git grep -F`. The only files it adds are three where `Rhodobacteraceae.1` is a prefix of
`.10+`, which are correctly excluded.

The rules, in priority order:

1. **The 30 `_m` old IDs are always MAG context.** Minted IDs are never ASV IDs.
2. **The 21 inherited old IDs (`X.N`, also live ASV IDs) are classified per column or JSON path, by its other
   tokens.**
   - MAG-only IDs and no ASV-only IDs make the column MAG.
   - ASV-only IDs and no MAG-only IDs make it ASV.
   - Both make it ambiguous.
   - In tables with a bin column, a column is MAG when at least 50 % of its rows equal the row bin's current ID,
     unless the header is a member, node, ASV or partner column.
   - The header rule matters because of a trap: in `d5_1_all_members_mag_proxy.tsv`, `member=Azospira.3` (the ASV)
     sits beside `best_mag=CAN_4_bin.208`, whose MAG ID was inherited from that ASV. The member column stays.
3. **Area rules** override the column result:
   - **ASV data** (always ASV): abundances, correlations, `network_method_comparison`, `dbRDA`, `leftover_reads_analysis`
     `asv_*`, `iterativeID_*.json`, `nonzero_per_day.json`, the ASV co-occurrence membership files.
   - **Community modelling** (members are ASVs): run provenance and pipeline code.
   - **Network MAG link**: ASV node labels stay, but MAG membership changes (D2).
   - **KBase object-name records** (`kbase_kofam_push/`, `rast/genomes/`): MAG.
4. **Text and code files are classified line by line.**
   - A line is MAG context if it has a bin, a MAG-only ID, `MAG`, `bin` or `GTDB` in it.
   - It is ASV context if it has an ASV-only ID, `ASV`, `member`, `16S` or `MiDAS` in it.
   - With both or neither, the line is ambiguous.

**Results.** Across the 1,391 in-scope files:

| action bucket (`files_plan.tsv` `action`) | files |
|---|---|
| leave: community_modeling run provenance (r12/r13 logs, runs, escher_maps, plan_*, analysis r12/data_inconsistencies; members are ASVs) | 1,068 |
| leave: ASV context | 124 |
| KBase-part dependent (r13 RAST/KOfam object-name records, `kbase_mag_rename_log.json`) | 60 |
| regenerate: network MAG link (D2) | 46 |
| data: taxonomy patch / apply IDs / re-run generator | 37 |
| leave: backup / superseded / historical note (D5) / code whose strings are ASV member names or comments | 23 |
| text: MAG sentences / R220 taxon names (manual) | 20 |
| code fix (section 4(d)) | 12 |
| documents part (`MANUSCRIPT_methods_results_discussion.md`) | 1 |

- The 75 ambiguous occurrences are all in 3 text files:
  - `clade_classify/NOSZ_TOP5_NOTES.md` (1).
  - `community_modeling/analysis/r13/annotations/DESIGN_ANSWERS.md` (9).
  - `design_answers/compute_answers.py` (2). Most of these are r13 member names, which are ASVs.
- They need a human read.

## 4. Ordered regeneration plan

The full plan is in `regeneration_order.tsv`: 24 steps, with inputs, outputs, runtime and whether data missing from the repo is
needed.

**Principle (as in c796fdf).**
- Make value-preserving relabels and bin-keyed patches to data tables, then re-render the figures.
- Do not re-run the one-shot or poplar-bound builders:
  - `build_mag_quality_table.py`, `build_mag_abundance_291.py` and `compose_static.py`.
  - `03_aggregate_ko_rpkm.py`: the committed 275-MAG RPKM tables come from the old mapping.
  - `build_mimag_hq.py` and `build_skani_novelty_291.py`.

### (a) Taxonomy source of truth (steps A1–A5)

Write one new script (`apply_r232_taxonomy.py`, not written here). It joins on the bin and copies the R232 strings
verbatim from `mag_lineage_r220_r232.tsv`.

| file | columns | cells that change |
|---|---|---|
| `gtdbtk.bac120.summary.tsv` | replace by `results/gtdbtk.r232.summary.tsv` (byte copy; this also removes the best-of-both ANI/AF/RED hybrid) | 130 lineages, 99 pplacer, 30 closest-genome |
| `mag_quality_summary.tsv` | Domain..Species | P 22, C 23, O 19, F 22, G 49, S 17 |
| `mag_abundance_summary.tsv` | Domain..Species (the `iterativeIDs` column is ASV: leave it) | P 22, C 23, O 19, F 22, G 48, S 17 |
| `mag_abundance_summary.json` | `taxonomy{Domain..Genus}` (keep the 6-rank schema; `iterativeIDs` is ASV) | P 22, C 23, O 19, F 22, G 48 |
| `leftover_reads_analysis/mag_taxonomy.json` (+ `mag_genera.json` list) | phylum..genus | P 22, C 23, O 19, F 22, G 49 |
| `mimag_hq_classification.tsv` | Phylum, Genus, Species | 22 / 49 / 17 |

**Not MAG taxonomy, so leave:**
- `mag_abundance_by_day*.csv`: the Kingdom..Species there is the mapped ASV's MiDAS taxonomy.
- The `iterativeID_*.json` files and the colour maps (section 1).
- The `mag_skani_hits/` R226 run and `server_handoff_out/tree62` (frozen).

### (b) Apply the new MAG IDs (steps B1–B3)

1. **Run `scripts/recompute_mag_ids.py` after A3.** It reads `mag_abundance_summary.json`.
   - **Gate:** its new json must equal `mag_id_map_r220_r232.json` for 291 of 291 (51 changes, 0 duplicates).
   - Its csv keeps the **pre-c796fdf** "old" columns, not the R220 IDs. Commit the R220 → R232 pair separately, for
     example as `mag_iterativeID_r220_to_r232.csv` from `mag_id_changes_r220_r232.csv`.
2. **Extend `scripts/apply_mag_ids.py`** with every `ADD` row of `apply_mag_ids_extension.tsv`. All of them are
   bin-join, whole-cell updates.

   | file | MAG column | ID column | note |
   |---|---|---|---|
   | `Table_S7_MAG_functional_classification.csv` | MAG | iterative_ID | 33 stale (pre-408cb39) IDs are also fixed |
   | `mag_guild_classification.csv` | MAG | iterative_ID | 33 stale |
   | `mag_percell_abundance.csv` | MAG | iterativeID | or re-run `build_percell_abundance.py` (B5) |
   | `ko_gene_abundance/ko_rpkm_combined_ext.tsv` | MAG | mag_iterativeID | 42 stale |
   | `ko_gene_abundance/ko_rpkm_per_sample_ext.tsv` | mag | iterativeID | 1,890 stale rows |
   | `gene_ab_figure/data/rpkm_long_MAG_gene_sample_phase.tsv` | mag | iterativeID | 9,450 stale rows |
   | `gene_ab_figure/data/rpkm_long_MAG_gene_sample_291_preview.tsv` | mag | iterativeID | 3,860 stale; 1,115 blank |
   | `gene_ab_figure/data/novel_bacteria_summary.tsv` | bin | iterativeID | no changed MAG (no-op) |
   | `fig1_tree/fig1_annotation.tsv` | tip | label, iterativeID | or re-render (C1) |
   | `fig1_tree/fig1_mags_table.tsv` | MAG | iterativeID | or `select_fig1_mags.py` |
   | `phase_rpkm_chord/{mag_category_assignments,phase_mag_gene_rpkm_long}.csv` | MAG | MAG_iteration_ID | stale IDs; or re-run 01 |
   | `community_modeling/analysis/r13/annotations/{taxonomy.tsv, panel/mag_panel.tsv, rast/genomes/_manifest.tsv}` | bin | iterativeID | or rebuild (B9); `_manifest` `genome_name` follows KBase |

   - The 9 files already listed in apply_mag_ids.py stay as they are. So does the
     `mag_abundance_summary.json` `mag_iterativeID` (all MAG context).
   - Do **not** add `member`, `rep_iterativeID`, `node`, `guild_node`, `source`/`target` or `asv_*` columns. They hold ASV
     IDs that equal a MAG ID only where the MAG inherited them.
   - Skip the `*.bak.*`, `superseded_*`, `test62_*` and `server_handoff_out/` copies.
3. **Genus/species patch by bin** of the derived MAG tables (taxonomy_plan rows with action "patch"):
   - `ko_rpkm_{combined,per_sample}{,_ext,_hq}.tsv` and `rpkm_long_*.tsv` genus.
   - `gene_ab_figure` `taxonomy_labels.tsv` (Genus, Species) and `gene_copy_per_mag.tsv` genus.
   - `nosz_clades*.tsv`, `nosz_cladeII_subclades.tsv` and `allbins_nosz_master.tsv` genus.
   - Keep each file's suffix convention.

### (c) Scripts to re-run (steps B4–C7; all local, from seconds to minutes; nothing needs poplar)

- **B4 `ko_gene_abundance/04_root_aggregate.py`.** Its rows are keyed by MAG-ID roots. The Ignavibacteria, Rhodobacter,
  Azospira, Flavipsychrobacter, Chiayiivirga, Sphingobacteriales and Paceibacteria rows go, and new roots appear.
- **B5 `scripts/build_percell_abundance.py`.** It uses the local `counts_291_mag_CAN_*.json`, the untracked
  `meta/results/01_nanoplot_raw` (present locally) and `mag_metagenome_abundance.csv`.
- **B6** `MAG_target_gene_filled*.xlsx`: owned by the documents part.
  - Run `reprocess_291/build_xlsx_291.py` with `EMILYKIN_REPO=~/Documents/EmilyKin` (it defaults to `/scratch1`), or
    patch the cells by bin.
  - B7, C2 (`12_novel_bacteria_summary.py`), C5 and `render_fig1_tree.py` read it.
- **B7 `scripts/build_table_s7_classification.py`.** It writes the Table S7 CSV and `mag_guild_classification.csv`
  from the xlsx.
- **B8 Table S6.**
  - Rebuild `mag_gtdb_taxonomy` from R232. That changes 129 lineages.
  - Rebuild `skani_best_taxonomy` against R232 references from the prep `skani_preview/`. The builder has poplar
    paths, so port it locally.
- **B9 r13 annotations.**
  - Fix `build_taxonomy.py` first: `FLAG` gets `Patescibacteriota` and `DIDERM` gets `Myxococcota_A`
    (coasm_bin.258). Without the fix, the CPR flag drops silently for 21 MAGs and one diderm call is lost.
  - Then re-run `build_taxonomy.py` and `panel/build_panel.py`, and patch `calls/summary.tsv`.
  - `CODE25` is keyed by bin. All 5 MAGs are still BD1-5 / Absconditabacterales in R232, so it is unaffected.
- **C1 Fig 1.** No RAxML re-run is needed: the tree tips are bins.
  - Fix the rooting in `render_fig1_tree.py` and `fig1_tree/compare_trees.py` (`Patescibacteria` or
    `Patescibacteriota`). With R232 names they would fall back silently to midpoint rooting.
  - Add `Patescibacteriota` and `Myxococcota_A` to `GTDB2PALETTE`.
  - Re-run `select_fig1_mags.py` and `render_fig1_tree.py --annotation-out`.
  - The 61 reference tips carry R220 names from `~/Documents/emilykin_tools/gtdbtk_r220_markers`. In R232, 21 of
    them change lineage, and GB_GCA_017302355.1 is no longer in GTDB (D3).
- **C2 gene_ab_figure:** `08` → `04` → `05_bundle_offline` → `06`, `09`, `10`, `12`, `13`. 2 of the 21 focal MAGs change
  ID.
- **C3 clade_classify:** `render_nosz_tree_fig.py`, `render_nosz_top5_tree.py`, `render_organism_tree.py` (after B7).
- **C4 `phase_rpkm_chord`:** `01` → `03` → `04`. This folder still carries pre-408cb39 IDs.
- **C5 other MAG figures:**
  - `render_chord_rpkm_function_guild.py`, `render_n2o_sink_producer.py`, `render_denitrification_heatmaps.py`.
  - `render_percell_composition.py`: fix its colour keys first.
  - `ppk1_classify/render_ppk1_tree_fig.py`.
- **C6 network (per D2).** Run `render_cooccurrence_grid_all_phases.py` and `render_gao_pao_networks.py` with
  `NETWORK_MAG_ONLY=1`, then the rest:
  - `build_module_mag_membership.py`.
  - `render_module_denitrification_genes.py` and `render_module_pathway_completeness.py`.
  - In `network/compositional_networks`: `render_spieceasi_mag_network.py`, `spieceasi_modules_flanking.py` (and the
    root copy, which writes `tablesforpaper/`), `render_module_flanking_figures.py` and `render_spieceasi_docx.py`.
    The last has hard-coded counts ("**53** of those that are MAGs") and edge lists.
  - Runtime: minutes. The all-phases Spearman, the Louvain gamma scan and the 16–27 MB grid PNGs take the time.
- **C7 community modelling: nothing to run.**
- **D1 verify.**
  - Join every MAG ID column on its bin against `mag_id_map_r220_r232.json`.
  - Check the ASV files are md5-identical before and after.
  - **Do not verify by "old-string count → 0".** The 21 old IDs remain valid ASV IDs, and the 4 new IDs already
    occur as ASV IDs: 20,184 times in 604 tracked files.

### (d) Hard-coded strings in scripts (`code_fixes.tsv`, 73 rows)

12 code files have an entry that is not "leave".
- **5 need a logic or key change:** `render_fig1_tree.py`, `fig1_tree/compare_trees.py`, `build_taxonomy.py`,
  `render_percell_composition.py`, and the `test_build_panel.py` fixture.
- **7 are comments, docstrings or a historical script:** `render_nosz_top5_tree.py`, `call_genes.py`,
  `validate_panel.py`, `score_ref_candidates.py`, `build_fig1_alignment.py`, `run_fig1_raxml.sh` and `run_tree62.sh`.
- **`render_spieceasi_docx.py` depends on D2.** It holds ASV node names, plus the hard-coded MAG-node count 53 and
  edge lists that change under D2a.
- **No tracked file is named by an old MAG ID.** `git ls-files | grep -F -f old_ids` returns 0 files. The ID-named
  `*_gf.xml` model files are named after ASV members.

| script | string | fix |
|---|---|---|
| `scripts/render_fig1_tree.py:137–150` | `r["Phylum"] == "Patescibacteria"` | accept `Patescibacteriota`; `GTDB2PALETTE` add `Patescibacteriota→Patescibacteria`, `Myxococcota_A→Myxococcota` (palette keys are MiDAS names) |
| `fig1_tree/compare_trees.py:96,103,114,143` | same rooting rule | same |
| `community_modeling/analysis/r13/annotations/build_taxonomy.py:43,48` | `DIDERM {'Myxococcota',…}`, `FLAG['Patescibacteria']` | add `Myxococcota_A`, `Patescibacteriota` (or strip GTDB `_X` suffixes before lookup) |
| `…/annotations/panel/tests/test_build_panel.py:106` | `p__Patescibacteria` lineage fixture | update the fixture |
| `scripts/render_percell_composition.py:15,72–78` | `DOCX_COLORS` genus keys `Giesbergeria`, `SCMT01`, `JJ008`, (`Azospira` for `Azospira_A`) | `Simplicispira`; `Lysobacter` (`Lysobacter_I` after suffix strip); JJ008 keeps 1 of 2 MAGs (other → CADEDV01); add `Proximibacter`; the top-20 band set may change |
| `fig1_tree/score_ref_candidates.py:3,7`, `scripts/build_fig1_alignment.py:623–624`, `fig1_tree/run_fig1_raxml.sh:5` | `g__Azospira_A` (CAN_4_bin.208) | comments/defaults: now `g__Proximibacter` |
| `clade_classify/render_nosz_top5_tree.py:29` | `CAN_3_bin.221 / Ignavibacteria.21` | docstring: `Tepidaquicellales.1_m` |
| `…/r13/annotations/verify_design/{verify_design.py:16, v7_storage_gpr_reach.py:11}` | `"Ignavibacteria.31": ["CAN_3_bin.221"]` | **leave**: `Ignavibacteria.31` is the ASV member; the bin is its MAG |
| `…/calls/call_genes.py:6`, `…/panel/validate_panel.py:80` | "five Patescibacteria MAGs" (text) | optional wording only; the logic is keyed by bin |
| `reprocess_291/run_tree62.sh:7` | `--outgroup_taxon p__Patescibacteria` | historical server script; leave unless re-run on R232 |

**Leave these.** They hold ASV, NCBI or MiDAS names:
- `community_modeling/member_reference_phyla.py`, `ek_pipeline.py:731`, `clade_classify/partA_reclassify.py:38`.
- The `render_spieceasi_docx.py` phylum mixes of ASV modules.
- `compute_answers.py`, `v8_a52_shares.py` and `curate_templates.py:649`: member names.

## 5. The directories named in the task

**community_modeling**
- Members and model IDs are ASV iterativeIDs: `Ca_Competibacter.109`, `Thauera.34`, `Azospira.3`, and so on.
  - The r12 and r13 configs, runs, logs, Escher maps, `member_colors.json`, `template_corrections.csv` and the
    `analysis/r12_*` / `data_inconsistencies/` tables all use them.
  - These 1,068 files hold 37,713 occurrences, and they are ASV context. The aggregated rows in `columns_plan.tsv`
    have one row per subtree.
- The MAG-keyed artefacts are under `analysis/r13/annotations/`:
  - `taxonomy.tsv`, `panel/mag_panel.tsv` and `calls/summary.tsv`: taxonomy and iterativeID.
  - `rast/genomes/_manifest.tsv` and `kbase_kofam_push/{placement_by_bin.tsv, verify/*.json}`: KBase genome names
    `emilykin_MAGs_RAST_code11_b<k>_<iterativeID>`.
  - Track-2 feature IDs embed the same names (`plan_r13/track2/GENOME_OBJECT_SPEC.md`).
  - `fetch_rast.py` checks the genome-name suffix against `kbase_mag_rename_log.json`. So the KBase rename and the
    new log must land together.
- **Re-running r12 or r13 is not needed and not proportionate.** The r12 75-day run took 4.4–4.8 days wall. A re-run
  changes nothing but the KBase feature-name strings in the MAG-backed Track-2 objects. The one affected Track-2 bin
  is CAN_3_bin.221, which backs the ASV member Ignavibacteria.31 and is renamed Ignavibacteria.21 → Tepidaquicellales.1_m.
- **Rebuild only the r13 annotation tables** (B9), after the `FLAG`/`DIDERM` fix. Their envelope and template calls
  must come out unchanged.

**network/** (and `tablesforpaper/`, `module_pathway_completeness/`)
- **The node IDs are ASV IDs and never change.**
- **The MAG link changes with the string identity.**

  | network | nodes | MAG-linked, R220 | MAG-linked, R232 |
  |---|---|---|---|
  | aggregate | 935 | 140 | 127 |
  | per-phase I | 371 | 60 | 52 |
  | per-phase II | 128 | 11 | 8 |
  | per-phase III | 870 | 137 | 125 |
  | per-phase IV | 284 | 92 | 84 |
  | per-phase V | 271 | 87 | 78 |
  | all-phases (MAG-only universe) | 140 | 140 | 127 |

  - **Aggregate and per-phase networks:** all ASVs. The edges and Louvain modules are unchanged; only the MAG
    annotation changes.
  - **All-phases network:** its edges, BH-FDR and modules must be recomputed.
  - **SPIEC-EASI:** the MAG display goes from 53 to 45 nodes (9 drop, 1 joins). The SPIEC-EASI fit on 89 ASVs is
    unchanged.
- **Dropped ASVs:**
  - Azospira.3, Bosea.3, Ca_Competibacter.32, Ca_Falkowbacteria.7, Gammaproteobacteria.7, Hydrogenedensaceae.8.
  - Ignavibacteria.5/.7/.13/.31, Myxococcaceae.11, Rhodobacter.5/.13, Rhodobacteraceae.1.
  - Sphingobacteriales.9/.13/.16.
- **Added ASVs:** Hoeflea.1, Holophagaceae.9, Simplicispira.1, Taibaiella.5.

**gene_ab_figure/**
- `data/taxonomy_labels.tsv`, `gene_copy_per_mag.tsv` and the `offline_bundle` copy are already in
  apply_mag_ids.py. 2 MAGs change; Genus and Species need patching.
- Add the `rpkm_long_*` files, and patch their genus.
- Re-run 04–13 (C2). `12_novel_bacteria_summary.py` reads Table S1 from the xlsx, so it waits for B6.

**clade_classify/**
- `out/allbins_nosz_master.tsv` is already listed: 27 changed MAGs, 22 genus cells change.
- Patch the genus in `nosz_clades_updated.tsv`.
- The trees use bin or locus tips, so no tree is re-run. Re-render the 3 figures.

**fig1_tree/**
- The annotation and MAG table are covered in C1. 8 of the 65 tip MAGs change ID, and 33 MAG taxonomy cells change.
- The tree files have bin and accession tips, so they do not change.
- Fix the rooting string.
- Reference-tip names: D3.

**dbRDA/**
- ASV/iterativeID-root based: 553 ASV roots, MiDAS genus.
- None of the 51 old IDs occurs there. **No action.**

**ko_gene_abundance/**
- `ko_rpkm_combined.tsv` and `ko_rpkm_per_sample.tsv` are already listed.
- Add the `_ext` pair to apply_mag_ids.py; the `_hq` pair needs only a genus patch (it has no ID column).
- Re-run `04_root_aggregate.py`. Do not re-run `03`/`07`: they need the old count JSONs, and the tables are a
  275-MAG basis.

**asv_mag_mapping/**
- `mag_16s.fasta`, the SAM file and `asvs.fasta` are keyed by bin and by ASV hash. `asv_to_mag_mapping.json` is keyed
  by hash and bin, with tiers that do not depend on taxonomy. **No action.**
- The asv_inherit naming rule reads this mapping. So the new IDs depend on it, but it is not changed.

## 6. Decisions needed (before the run)

- **D1. Commit an explicit R220 → R232 map.** `recompute_mag_ids.py` overwrites the map with the R232 IDs and keeps
  only the ancient "old" columns. Recommendation: commit `mag_id_changes_r220_r232.csv` alongside it.
- **D2. Network MAG link (`NETWORK_MAG_ONLY` and its users).** Choose one:
  - **(a)** Keep the string-identity rule. The MAG-only results change, not just their labels (section 5): all-phases
    network 140 → 127 nodes, new edges and modules, GAO/PAO figures, SPIEC-EASI 53 → 45, module panels,
    `tablesforpaper/`. The numbers in the SPIEC-EASI docx and the manuscript follow.
  - **(b)** Decouple. Freeze the R220-era list of 140 MAG-representative ASVs as an explicit file, or derive it from
    `asv_to_mag_mapping.json`. Then only the MAG labels change, through the bin.
  - Recommendation: (b), unless the paper should report the R232-consistent network.
  - Caveat on (b): for 17 nodes, it keeps a MAG link whose stated basis is taxonomic agreement at the deepest
    non-numeric rank, and R232 no longer supports that agreement.
    - So the frozen list needs its rule written down: either "R220-era MAG-representative ASVs" or "ASV → best_mag
      tier from `asv_to_mag_mapping.json`".
  - Whichever option is chosen, all seven consumers must move together: `render_cooccurrence_grid*.py`,
    `render_gao_pao_networks.py`, `build_module_mag_membership.py`, `render_module_denitrification_genes.py`,
    `render_module_pathway_completeness.py`, `render_spieceasi_mag_network.py`, and both
    `spieceasi_modules_flanking.py`.
- **D3. Reference-genome taxonomy:**
  - Table S6 `skani_best_taxonomy`, `closest_genome_taxonomy`, and the Fig 1 reference tips (21 of 61 change; 1 is
    no longer in GTDB).
  - Options: relabel from the R232 taxonomy files (`gtdb_downloads/r232`), or leave with a note.
  - The species-level novelty buckets do change: 122 species-level MAGs in R232, against 96 in R220.
- **D4. Stale IDs.** Table S7, `mag_guild_classification.csv`, `ko_rpkm_*_ext`, `rpkm_long_*` and `phase_rpkm_chord/`
  still hold pre-408cb39 IDs. A bin-keyed apply overwrites them with R232 IDs, which is presumably wanted.
- **D5. Scope of text edits.** Most `.md` notes are historical records, as are the community-modelling plans and
  reports. Edit only living documents: READMEs, METHODS, and the files listed with action "edit text" or "check R220
  taxon names".

## 7. Risks and gotchas

- **Order matters.** The sequence is: patch the taxonomy (A3), run `recompute_mag_ids.py` (B1) and check it against
  this map, then apply (B2).
  - If B2 runs from the old `_new_ids.json`, nothing changes.
  - If B1 runs before A3, it reproduces the R220 IDs.
- **No sequential string replacement anywhere.**
  - The chain `Chitinophagaceae.1_m` (coasm_bin.117 → `.2_m`; CAN_5_bin.149 takes `.1_m`) breaks it.
  - So does `Simplicispira.1` against `Simplicispira.1_m`, two different MAGs after the rename.
  - So does the `Rhodobacteraceae.1` prefix of `.10`–`.19`.
- **`Patescibacteria`, `Myxococcota`, `Burkholderiaceae`, `Bosea` and `Flavobacteriales` are live MiDAS ASV names.** A
  global replacement would corrupt `taxonomy.json`, `iterativeID_*.json`, the palettes and every ASV figure.
- **The 4 new asv_inherit IDs make the MAG and the ASV share a label.** The repo already does this for 138 other MAGs,
  but the D2 consequences above follow from it.
- **`gtdbtk.bac120.summary.tsv` loses the best-of-both hybrid.** ANI/AF/RED now come from the committed sequences, and
  scripts that read `closest_genome_*` see R232 references. This is intended.
- **Do not re-run poplar-bound or one-shot builders** (section 4 principle).

## 8. Machine-readable outputs (all TSV, one header row)

| file | rows | one row per |
|---|---|---|
| `files_plan.tsv` | 1,391 | in-scope tracked file: area, `generated` (curated / heuristic / unknown), `generated_by`, readers, occurrence counts (`occ_old_minted`, `occ_old_inherit`, `occ_MAG`, `occ_ASV`, `occ_ambiguous`, `occ_new_id_already_ASV`), `mag_taxonomy_cells_change`, `r220_name_hits`, `code_strings_to_fix`, `action` |
| `columns_plan.tsv` | 778 | (file, column / JSON key path / text part): context, occurrence split, old IDs found, basis of the call, `generated_by`, `action`; community_modeling provenance aggregated to 74 subtree rows |
| `taxonomy_plan.tsv` | 215 | (file, MAG taxonomy column): rank, GTDB vs MiDAS vs reference-DB, cells matching R220, cells that change in R232, action |
| `apply_mag_ids_extension.tsv` | 41 | (file, delimiter, MAG column, ID column): current / stale / blank counts and status (`ADD`, already listed, skip, KBase) |
| `code_fixes.tsv` | 73 | (code file, line, string): kind, context, fix |
| `binary_for_docs.tsv` | 319 | xlsx/docx/pdf/png: old-ID and R220-name hits from the OOXML/PDF text, plus the heuristic writer script, for the documents part |
| `regeneration_order.tsv` | 24 | ordered step: script, inputs, outputs, data not in repo, runtime, depends_on |

- **`generated_by` is curated** for about 50 key files (section 4).
- **For the rest, it is a heuristic:** a tracked script names the file's basename near a write call. A basename shared
  by several tracked files (for example `README.md`) is not attributed.

## 9. Hand-offs to the other parts

**Documents part**
- `binary_for_docs.tsv`. Old IDs occur in:
  - `MAG_target_gene_filled*.xlsx`: up to 175 `_m` + 124 inherited per workbook.
  - `cladeexcel*.xlsx`, `abundance and coverage (1).xlsx`.
  - `SPIEC-EASI_network_analysis.docx`: ASV node names, D2.
  - `papert3_070626 (2)_correction.docx`.
- R220 names occur in all the manuscript and SI docx files: `Giesbergeria`, `Azospira_A`, `Patescibacteria`/`Myxococcota`
  in MAG sentences only, `Accumulibacter phosphatis_G` / `sp017302455` (Wastewater manuscript-SI).
- The 271 PNGs are listed with their heuristic generator.
- The xlsx is upstream of B7, C2 and C5.

**KBase part**
- `kbase_mag_rename_log.json` (`new_name`) is consumed by `build_taxonomy.py`, `build_panel.py` and `fetch_rast.py`.
- The r13 records with object names (`_manifest.tsv` `genome_name`, `placement_by_bin.tsv`, `kbase_kofam_push/verify/*.json`,
  `plan_r13/kbase_rast/rast_jobs.json`) are provenance. Refresh them only if the r13 fetch is re-run after the KBase
  rename.

## 10. How to regenerate this folder (read-only on the repo; uses ~/Documents/py_venv)

```
P=rename_plan/data/scripts; PY=~/Documents/py_venv/bin/python
$PY -I $P/build_maps.py        # maps + verification            (<10 s)
$PY -I $P/scan_repo.py         # token scan of all tracked files (~4.5 min) -> scan/*_raw.tsv
$PY -I $P/gen_index.py         # basename -> writer/reader scripts (~2 min) -> scan/generators.tsv
$PY -I $P/bin_id_columns.py    # bin-keyed ID / taxonomy columns (~40 s) -> scan/bin_keyed_*.tsv
$PY -I $P/code_hardcoded.py    # hard-coded strings in code (2 s) -> scan/code_hardcoded_raw.tsv
$PY -I $P/assemble_plan.py     # final TSVs (5 s)
$PY -I $P/write_regen_order.py # regeneration_order.tsv
```
