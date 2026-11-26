# EmilyKin R220 → R232 rename plan: completeness and ambiguity critique

Prepared 2026-11-11. Read-only: `~/Documents/EmilyKin` is still clean at 9e087c3, and nothing was written to KBase.
This reviews the three part plans (`data/PLAN.md`, `docs/PLAN.md`, `kbase/PLAN.md`). The supporting scan and scripts
are in `critic/`.

## 1. Verdict

The three plans are sound on the main points:
- every rename is keyed by bin or object ID, never done as a string replacement;
- the 21 old IDs that are also ASV IDs are handled correctly;
- the 4 new IDs that are already ASV IDs are handled correctly;
- the Chitinophagaceae.1_m chain is handled;
- KBase rename order and rollback are correct.

They are not yet executable as written:
- **Two items change results and need a decision first.** These are the network decision D2 and the Table S6 input.
- **Two steps have duplicate owners.** Table S6 and Table S7 are each rebuilt by both the data part and the docs
  part.
- **Several steps are not runnable locally or are ordered wrongly.** See section 3.
- **The apply tool would rewrite line endings.** Every LF file it touches would be converted to CRLF.

Sections 2 to 6 list the corrections. Section 7 gives the independent scan, and section 8 the re-check of classifications.

## 2. Required corrections: decisions and results

**C1. D2, option (b) "derive from `asv_to_mag_mapping.json`", is not label-only.**
- Of the 935 network ASV nodes, 689 have a non-weak `best_mag` (226 distinct bins); 140 are linked by ID identity.
- So option (b) defined that way changes the MAG-only network far more than D2a does (140 → 689 nodes, not 140 → 127).
- **Only a frozen HEAD list keeps the results unchanged:** the 140 (ASV ID, bin) pairs. All 140 have
  `best_mag` equal to the bin that holds the ID (checked).
- Write the list before the new map is written (for example `network/mag_representative_asvs_r220.tsv`), and point
  every consumer at it.
- The consumers include the two inverse maps as well:
  - `scripts/render_module_denitrification_genes.py:158` and `render_module_pathway_completeness.py:100` (`iid2mag`);
  - `spieceasi_modules_flanking.py:41-42` (`MAG_IDS`, `NEW_TO_BIN`);
  - `render_cooccurrence_grid.py:107`;
  - `render_cooccurrence_grid_all_phases.py:80`;
  - `build_module_mag_membership.py:32`;
  - `render_spieceasi_mag_network.py:36`;
  - the root copy `spieceasi_modules_flanking.py:62`.
- Under the R232 map, `iid2mag` silently loses the 17 dropped nodes and links the 4 new ones. It is silent because
  both call sites guard the lookup with `if m in iid2mag` (`render_module_pathway_completeness.py:108`,
  `render_module_denitrification_genes.py:187`), so nothing raises.
- Both plans recommend "(b)". Restate it as "(b-frozen)" and record the rule: "R220-era MAG-representative ASVs".

**C2. D2 label sub-decision (new). After R232, the MAG-network node label and the MAG's ID differ for 17 nodes.**
- Today the node label (the ASV ID) equals the MAG ID, so no one had to choose between them.
- After R232 they differ for 17 nodes. Examples: node `Azospira.3` is bin CAN_4_bin.208, which is now
  `Proximibacter.2_m`; node `Ca_Competibacter.32` is coasm_bin.365, now `Contendibacter.1_m`.
- `network/louvain_module_mag_membership.json` stores `{"iterativeID": <node ID>, "MAG": <bin>}`. After R232 the
  meaning of that `iterativeID` field is ambiguous.
- Decide for the MAG-network figures and tables which label to use:
  - the ASV ID (both plans' default), with the MAG ID in a legend or table; or
  - the MAG ID.
- Then rename the `iterativeID` field to `asv_iterativeID` (or add `mag_iterativeID`). The SPIEC-EASI docx
  "gene-content PAO-like MAG (Azospira.3 …)" sentences are this same issue.

**C3. Cross-part count mismatch.**
- docs `PLAN.md` R1 says "21 ASVs drop out, 4 join". Within the network universe the correct figure is **17 drop,
  4 join** (140 → 127), as the data plan says. I verified this independently.
- Ignavibacteria.21, Ignavibacteria.25, Paracoccus.4 and Mesorhizobium.7 are not network nodes.
- SPIEC-EASI goes 53 → 45: 44 nodes are kept and 1 joins (verified).

**C4. Table S6 input: the prep `skani_preview` is not adequate.**
- `skani_preview/skani_search_r232.tsv` was run with skani defaults (min-AF 15 %, no `-n`). Only 219 of 291 MAGs
  have any hit.
- Table S6 was built from `skani search --min-af 0 -n 3` (`reprocess_291/run_skani.sh`). The builder labels a best
  hit with AF < 15 % as `low_AF_unresolved`, and a MAG with no hit as `<85`.
- With the preview, the 72 hitless MAGs would be wrongly called `<85`.
- Fix: run a new local search before S6 is rebuilt, with
  `-d ~/Documents/gtdbtk_data/release232/skani/database --min-af 0 -n 3`. The database is present locally.
- Use `gtdb_downloads/r232/bac120_taxonomy_r232.tsv.gz` as the reference taxonomy.
- Data B8 must be changed to this. Docs step 1 already hedges; resolve it to "run the new search".

**C5. Record the `Contendobacter` spelling as a decision. It changes no current result.**
- R232 names coasm_bin.365's genus `Contendibacter`. The repo spells it `Ca_Contendobacter` in three places:
  - `gao_pao_categories.json` (GAOs);
  - `scripts/_gao_pao_label_helper.py`;
  - the Table S7 regex `Competibacter|Contendobacter` (`build_table_s7_classification.py:141`).
- **No current consumer changes:**
  - The label helper and the GAO/PAO networks apply these lists to ASV (MiDAS) genera. The GTDB spelling does not
    matter there.
  - Table S7 calls coasm_bin.365 canonical GAO through its family rule (docs R5).
  - `phase_rpkm_chord/01_build_long_table.py` `guild_of()` has no family rule, but coasm_bin.365 is not one of its
    19 focal MAGs.
  - No genus-to-guild call among the 291 MAGs changes from R220 to R232 (checked).
- **Why record it anyway.** If the spelling is later aligned to GTDB, any future genus-only rule on MAG genera would
  flip coasm_bin.365 to "Canonical GAO" without anyone noticing. Note the decision in the docs part's R5; do not let it
  happen as a silent "typo fix".

## 3. Required corrections: ordering, ownership and runnability

**C6. One owner for Table S6 and one for Table S7.**
- Data B8 and docs steps 1–3 both rebuild S6, with different inputs.
- Data B7 and docs step 5 both run `build_table_s7_classification.py`.
- Assign each to the docs part. The data part keeps only its downstream dependency edges.

**C7. `scripts/build_table_s6_ani.py` writes R220 provenance.**
- Line 40 sets `REF_DB = "GTDB R220 (GTDB-Tk v2.7.2, classify_wf)"`. Line 307 writes the note text "skani vs GTDB r220
  species reps", and line 173 has R220 help text.
- Docs step 3 runs this builder, so it would write R220 strings into the builder workbook's S6 sheet.
- Patch it first. Add it to `code_fixes.tsv`; it is missing there.

**C8. C1 (Fig 1 / SI Fig S4) has missing dependencies and a hidden R220 lookup.**
- `render_fig1_tree.py:73` reads `Table_S6_skani_novelty.tsv`: the red-circle markers, ANI and ANI_class. So C1
  depends on S6 (C4, C6), not only on A1, B1 and B6.
- Reference-tip names come from `RES/sp_clusters_r220.tsv` (line 86), with `RES = ~/Documents/emilykin_tools/
  gtdbtk_r220_markers`. The code then asserts that every reference is in R220.
- So D3's option "relabel from R232" is a code change:
  - point the code at `gtdb_downloads/r232/sp_clusters_r232.tsv`;
  - handle GB_GCA_017302355.1, which is absent from R232 and would trip the assert;
  - re-check the `labels unique` assert.
- **Rooting fix.** Accept *both* `Patescibacteria` and `Patescibacteriota`. Under a relabel-MAGs-only D3, the
  reference tips stay `Patescibacteria`, so a single-name test gives a partial outgroup and a silent midpoint fallback.
- Gate: the printed `how` string must still say "Patescibacteria clade (N tips)", as at HEAD.

**C9. Scripts listed as plain re-runs have hard-coded server paths.**

| step | script | line | what it hard-codes |
|---|---|---|---|
| C4 | `phase_rpkm_chord/01_build_long_table.py` | 31 | `REPO = "/scratch1/…/EmilyKin"`. Data says "data_not_in_repo: no". It needs a local REPO patch; its inputs `gene_ab_figure/data/counts_CAN_*.json` and `read_counts_raw.tsv` are local, but untracked |
| C2 | `gene_ab_figure/05_bundle_offline.py` | 16 | `WORK = /scratch1/…` |
| docs step 2 | `reprocess_291/build_skani_novelty_291.py` | 6–10 | all `/scratch1` paths |
| docs step 7 | `reprocess_291/compute_manuscript_stats_291.py` | 5, 155 | `/scratch1` paths |

- For `05_bundle_offline.py`: drop it from C2, because `apply_mag_ids.py` already patches
  `offline_bundle/data/taxonomy_labels.tsv`. Or patch WORK.
- `build_xlsx_291.py` does honour `EMILYKIN_REPO` (line 35), and all its inputs exist locally (checked).

**C10. B9 (r13 tables) depends on the KBase rename log.**
- `build_taxonomy.py:38` and `panel/build_panel.py:71` read bin → iterativeID from `kbase_mag_rename_log.json`, not
  from `mag_iterativeID_old_to_new.json`.
- `rast/fetch_rast.py:166` checks each Genome object name's suffix against that log.
- Sequence:
  1. B1 gate (final map).
  2. KBase step (b) (renames).
  3. Update the log's `new_name`, in the same change set as step (b).
  4. B9.
- The r13 rebuild runs offline, but a RAST re-fetch fails if the log and KBase disagree.

**C11. Line terminators: `apply_mag_ids.set_mag_col` writes CRLF.**
- It uses `csv.writer` with the default `lineterminator='\r\n'`.
- The 9 files it handles today are already CRLF, which is how that happened.
- The proposed ADD targets are LF:
  - `Table_S7_MAG_functional_classification.csv`, `mag_guild_classification.csv`, `mag_percell_abundance.csv`;
  - `gene_ab_figure/data/rpkm_long_*.tsv`, `fig1_tree/fig1_mags_table.tsv`;
  - r13 `taxonomy.tsv`, `panel/mag_panel.tsv` and `rast/genomes/_manifest.tsv`.
- Extending the tool as it is would rewrite every line, and leave `\r` on the last field for readers that split
  lines themselves.
- Fix: detect and preserve each file's terminator and quoting.
- Gate: `git diff --numstat` shows only the rows of changed bins.
- A1 (byte copy of the R232 summary) also turns `gtdbtk.bac120.summary.tsv` from CRLF into LF. That is harmless for
  csv readers, but say so in the commit.

**C12. A wrong fact in the data plan: `ko_rpkm_combined_hq.tsv` and `ko_rpkm_per_sample_hq.tsv` do have ID columns.**
- The columns are `mag_iterativeID` and `iterativeID`. They are blank in all rows (ffe668c era).
- The data plan says "no ID column; genus patch only".
- The consumers (`render_chord_rpkm_function_guild.py`, `build_xlsx_291.py`) use the map, so blank is safe. Fill
  them by bin or leave them blank, but state which.

**C13. A wrong fact in the data plan: the all-phases network is listed as "leave".**
- `files_plan.tsv` has `network/network_module_membership_p_value_FDR_allphases.json` as "leave (ASV context)". It is
  the MAG-only network (`network_summary_p_value_FDR_allphases.json` gives n_nodes = 140).
- Under D2a both files change, and the summary file is missing from the plan. Move both to "regenerate per D2".

**C14. GAO/PAO network mode is unverified.**
- C6 says to run `render_gao_pao_networks.py` with `NETWORK_MAG_ONLY=1`.
- The script's docstring says its aggregate graph uses "all ASVs", and its documented run line sets no environment.
- Before attributing any change to D2, establish which mode produced the committed PNGs: re-run at HEAD and compare.

**C15. Add a baseline gate for every regenerate step.**
- Before changing anything, re-run each table or figure script at HEAD into a scratch copy and diff it against the
  committed output.
- Only a script that reproduces HEAD can be trusted to show the R232 effect.
- c796fdf left 15 figures with no generator; any figure that fails the baseline goes to the same manual list.

**C16. B5 has an input missing from the plan.**
- `build_percell_abundance.py` also reads `MAG_target_gene_filled_UPDATED (1).xlsx`, for its numeric regression check.
- That does not change the order, because the check is numeric only. List the input.
- The other untracked inputs are present: `meta/results/01_nanoplot_raw/CAN_{1..5}/NanoStats.txt` and
  `ko_gene_abundance/counts_291_mag_CAN_*.json` (checked).

**C17. All 2,910 formula cells in the submission workbook's S1 sheet have no cached values.**
- This is already the case at HEAD, from an earlier openpyxl save.
- Any script that reads that sheet with `data_only=True` gets `None` for those cells.
- The docs part's surgical edit will not make this worse, but do not add a reader of those columns. Or recalculate with
  LibreOffice headless after editing.

## 4. Required additions: chains and collisions

**C18. There is more than one chain once the stale-era files count.**
- The data plan says "one chain". That is true for HEAD → R232 only.
- The 6a04588-era files hold some IDs that are now held by other bins:
  - Files affected: all workbooks, the Table S7 CSV, `mag_guild_classification.csv` and
    `abundance and coverage (1).xlsx`.
  - `Chitinophagaceae.2_m` = coasm_bin.61 (HEAD `.14`). After R232 it is coasm_bin.117's ID.
  - `Chitinophagaceae.14` = coasm_bin.409 (HEAD `.18`).
  - At 924fa18, `Ignavibacteria.21` and `.31` are swapped.
  - In the pre-c796fdf columns, `Flavipsychrobacter.1_m` was coasm_bin.213.
- The bin-keyed join handles all of these. The docs `table_cells_plan.tsv` already has the coasm_bin.61 and
  coasm_bin.409 drift cells (checked).
- **Add a D1 gate:** after the apply, each MAG ID column is unique within its file and equals the R232 map for every
  bin.

**C19. An ID encoding the tokenizers do not look for: `Root_m.N`.**
- The pre-c796fdf form is `Root_m.N` (`Giesbergeria_m.1`, `Holophaga_m.1`, `Luteimonas_D_m.1`,
  `Zambryskibacteraceae_m.1`, `Flavipsychrobacter_m.1`, `Giesbergeria_m.2`).
- It occurs in:
  - `phase_rpkm_chord/{mag_category_assignments,phase_mag_gene_rpkm_long}.csv` and their PNG labels;
  - `MAG_target_gene_empty_tables.xlsx`.
- Both are covered, by the bin join plus C4 regeneration and by docs "do not edit". But no verification grep for
  `X.N_m` will ever see this form. The D1 bin-keyed gate is the only check that catches it.

**C20. Prefix hazards (confirmed, none new).**
- Old IDs: only `Rhodobacteraceae.1` is a prefix of a live ID (`.10`–`.19`).
- New IDs: only `Simplicispira.1` is a prefix (of `Simplicispira.1_m`).
- No other old or new ID is a digit-continuation prefix of any ASV or MAG ID.

**C21. KBase (the plan passes; two small additions).**
- Re-simulated: no collision in forward order or in rollback order.
- Object IDs are contiguous 1–891, so no deleted object reserves a name.
- The rename log matches the inventory for 291 of 291 bins (by objid, name and HEAD ID).
- The swap chain is ordered: Assemblies #11 before #12, Genomes #62 before #63.
- Additions:
  - **Pre-flight.** Assert that none of the 102 target names exists in 262319 (the four ASV-shared names
    `Simplicispira.1`, `Holophagaceae.9`, `Taibaiella.5` and `Hoeflea.1` included), except `Chitinophagaceae.1_m`
    before #11. Also assert the name of every object is either its old or its new name, so a partial run can resume by
    objid.
  - **Map provenance.** Build `renames.tsv` from the committed R232 map that passes the B1 gate, not from the
    simulation TSV. Run step (b) after B1 and before B9 (C10).

## 5. Ordering summary (merged; replaces the per-part orders where they conflict)

**Step 0. Decide.**
- D2: (a), or (b-frozen) plus the label rule (C1, C2).
- D3: relabel the references (needs the code change in C8) or leave them.
- The Contendobacter spelling (C5).
- Freeze the network ASV → bin list (C1).

**Step 1. Baseline gate (C15).**

**Step 2. Taxonomy and map (data A1–A5, B1).**
- Commit the R220 → R232 map (D1).
- Gate: 291 of 291.

**Step 3. KBase steps (a) and (b), then the `kbase_mag_rename_log.json` update (C10, C21).**

**Step 4. Data B2 and B3.**
- Use the newline-preserving apply (C11).
- Gate: per-file uniqueness (C18).

**Step 5. Data B4 and B5, then B9.**

**Step 6. Tables and workbooks (docs part).**
1. The new skani search (C4).
2. S6, with the builder patched (C7).
3. The builder workbook (B6).
4. S7 and the guild CSV (C6).
5. The submission workbook and the collaborator workbooks: cell edits keyed by bin.

**Step 7. Figures.**
- C1, with C8 applied.
- C2, without `05_bundle_offline`.
- C3.
- C4, with the path patched.
- C5.
- C6 (per D2).

**Step 8. Text (docs part), then KBase step (c), the AssemblySet re-save.**

**Step 9. Verification (D1).** Use the bin-keyed gate. Also check that the ASV files are md5-identical before and after.

## 6. Smaller corrections

**C22. Untracked PDF siblings of regenerated figures still carry old IDs** (the plans list only tracked files).
- The PDFs:
  - `fig1_tree/fig1_tree.pdf`, `denitrification_MAG_heatmap.pdf`, `clade_classify/out/nosz_top5_tree.pdf`,
    `clade_classify/organism_tree/organism_tree.pdf`;
  - the 5 `network/compositional_networks/figures/spieceasi_*.pdf`;
  - `network/{N2O_sink_producer,chord_*}.pdf`, `figures/percell_composition.pdf`.
- The same scripts rewrite them. Add them to the output lists, and to the D1 check via `pdftotext`, so that a stale PDF
  is not inserted into a DOCX.

**C23. Frozen KBase feature IDs: list the two large untracked r13 tables as "leave".**
- `rast/crosswalk.tsv` (255 MB; e.g. 3,540 `…_Giesbergeria.1_m.CDS.*`) and `rast/concordance_pairs.tsv` (186 MB) embed
  the old KBase names in `rast_feature_id`.
- So does every `rast/genomes/<bin>.tsv` (e.g. 3,648 in CAN_5_bin.40).
- All of them carry a `bin` column. They fit the KBase plan's "internal IDs frozen; join on bin" rule.
- List them explicitly as "leave (frozen KBase feature IDs)" so that no one string-replaces them.

**C24. The r13 Track-2 backing of member Ignavibacteria.31 is CAN_3_bin.221 on purpose, not by mapping.**
- `panel/member_mag_choice.tsv` records it: "DESIGN_R13 3.1(c): not best_mag CAN_5_bin.167".
- The ASV's own `best_mag` is CAN_5_bin.167, the bin that holds the MAG ID `Ignavibacteria.31` at HEAD. At 924fa18,
  CAN_3_bin.221 held that ID.
- The data plan's "leave" is right. But its rationale ("the bin is its MAG") should cite the design choice, because
  `DESIGN_ANSWERS.md` "Ignavibacteria.31 (CAN_3_bin.221)" reads like a 924fa18-era MAG label.

**C25. Leftover-reads outputs downstream of A4 are not in the plan.**
- Data A4 patches `leftover_reads_analysis/mag_taxonomy.json` and `mag_genera.json`.
- Those files are inputs to `parse_paf_16s.py:39` and `aggregate_report.py:19`. Their outputs label reads as
  "recovered_genus_relative" by MAG genus, and still carry R220 GTDB genus names:
  - `leftover_taxa.csv`, which has Holophaga, UBA919, JACRCT01, JAIEMT01, UBA2386 and UBA5200;
  - `out/CAN_*.ont.json`;
  - `leftover_summary.csv`, `leftover_analysis.json`, `REPORT.md` and `COLLABORATOR_SUMMARY.md`.
- Decide between two options:
  - re-run `aggregate_report.py`, and `parse_paf_16s.py` if its PAF inputs exist locally; or
  - freeze the outputs and mark them as "R220 genus labels".
- Patching only the input JSONs leaves the inputs and outputs inconsistent.

## 7. Independent scan (critic/scan/)

**Method** (`critic/scripts/critic_scan.py`).
- Files scanned:
  - every tracked file, at 9e087c3;
  - untracked and ignored files, except the trees listed under "Coverage caveats";
  - 39,123 files in all.
- Containers opened:
  - `.gz` files;
  - all XML parts of xlsx, docx and pptx: charts, comments and embedded OOXML, plus docx/pptx paragraphs with their
    runs joined;
  - PDF, through `pdftotext`;
  - FASTA header lines only.
- What was searched for:
  - the 51 old IDs and the 4 ASV-shared new IDs, in strict form;
  - `strict_embedded`: the ID preceded by `_`, as inside KBase names;
  - `ca_display`: `Ca. X.N`, `Candidatus X.N` and a bare `X.N`;
  - `sep`: root followed by `_`, `-` or a space, then N;
  - `escaped`: markdown/LaTeX `\_`.
- `compare_coverage.py` compares the hits with `data/files_plan.tsv` and `docs/occurrences.tsv`.

**Results**

| check | result |
|---|---|
| tracked files with an old ID | 1,326 (119,271 occurrences) |
| tracked files missed by both plans | 9, all `*.pre_*.bak.xlsx` (docs "do not edit"; they are in data `binary_for_docs.tsv`). **No live file is missed.** |
| new ASV-shared IDs | 605 tracked files, 20,190 occurrences (data: 604 / 20,184; the small difference is from docx joining and PDF extraction) |
| `sep` variants (`Ignavibacteria_5`, `Giesbergeria_1_m`, …) | 0 |
| escaped variants | 0 |
| docx/pptx IDs visible only after joining runs | 0 |
| `ca_display` | 2: `Ca. Competibacter.32` in `SPIEC-EASI_network_analysis.docx` and `render_spieceasi_docx.py`. Both are ASV nodes, and both plans already have them |
| `strict_embedded` (KBase names) | 267 rows. The 58 tracked files are all in the data plan (KBase-part dependent). The rest are the untracked `kbase_kofam_push/payloads` and `rast/genomes` files (C23) |
| JSON files keyed by an old ID | all ASV context (correlations, `asv_iterativeID_old_to_new.json`, r13 member reports), all in the data plan as "leave" |
| a file or path named by an old ID | 0, tracked or untracked |
| untracked hits | the regenerated PDFs (C22); r13 `_freeze1`, `_prefix_full`, `probe/`, `_pre_finalize` and `review/` snapshots (historical, leave); `escher_maps/prototype*` (ASV members, leave); `Experimental Data Inconsistencies.pdf` (a historical report, leave); 10 ignored `*.bak.xlsx` |
| `Root_m.N` (pre-c796fdf form) | found by a separate `git grep` only in `phase_rpkm_chord/*.csv` and `MAG_target_gene_empty_tables.xlsx` (C19) |
| IDs split across columns (root + number) | the only `root`/`suffix` tables are ASV-root tables: community-modelling `member_day.csv`, `q3_*_by_root_phase.csv` and `compositional_networks/out/root_collapsed/`. The only MAG-root tables are `ko_rpkm_root_*` (B4, regenerated) |
| R220 GTDB-only taxon names (42 that are not also MiDAS names) in code | 21 `.py/.sh/.ipynb` files; all are in `code_fixes.tsv` or `files_plan.tsv`. `data_processing.ipynb` hits are MiDAS ASV outputs |
| the same names in tracked data/text files | 114 files; 67 are in a plan and 47 are not. Of the 47, 38 are reference or ASV context, or frozen runs: `bvbrc_*` (14), `synthetic_genomes*` (18), the `mag_skani_hits/` R226 run (3), the `server_handoff_out/` R220 skani raw (1) and r13 `guild_exceptions.tsv` with its `.pre_full` copy (2; JARSSO01 in a family-rule note). 3 are outputs of regeneration steps that `files_plan.tsv` does not list: `ko_rpkm_root_{combined,per_sample}.tsv` (B4) and `leftover_reads_analysis/mag_genera.json` (A4). **The remaining 6 are a real gap: the leftover-reads outputs (C25)** |

**Coverage caveats** (`skipped.tsv`).
- Untracked `community_modeling/{runs,_pinned,plan_evidence,analysis/data_inconsistencies}/` and
  `metabolic_complementarity/profiles/`, 17,223 files, were not scanned. They are ASV run trees.
- 17 untracked files over 60 MB were not scanned: minimap2 indexes, kofam hits, `bvbrc_*` ASV JSON, and
  `rast/crosswalk.tsv` / `concordance_pairs.tsv`, which were checked by hand under C23.
- PNG pixels cannot be scanned; their generators are covered by the plans.

## 8. Re-check of the MAG/ASV context calls (24 occurrences)

| # | file / location | ID | part's call | critic verdict |
|---|---|---|---|---|
| 1 | SPIEC-EASI docx para 115 (edge list) | Azospira.3 | ASV-keep (docs) | agree: SPIEC-EASI node = ASV |
| 2 | SPIEC-EASI docx para 118 | Ignavibacteria.5 | ASV-keep | agree |
| 3 | SPIEC-EASI docx para 127 | Ignavibacteria.13, Ca. Competibacter.32 | ASV-keep | agree |
| 4 | SPIEC-EASI docx paras 153 and 165 | Azospira.3 — Ca. Competibacter.32/.51 | ASV-keep | agree |
| 5 | SPIEC-EASI docx para 140 "gene-content PAO-like MAG (Azospira.3 …)" | Azospira.3 | ambiguous-author | agree. It is resolved by the C2 label rule |
| 6 | SPIEC-EASI docx para 140 "Ca. Competibacter (.26, .32, .51)" | .32 | ambiguous-author | agree. An ID regex cannot see an abbreviated list |
| 7 | SPIEC-EASI docx para 357 | Azospira.3 | ambiguous-author | agree |
| 8 | `papert3_070626 (2)_correction.docx` Table 1 row 4 | Ignavibacteria.21 (CAN_3_bin.221) | MAG-rename → Tepidaquicellales.1_m | agree (row keyed by bin) |
| 9 | `MAG_target_gene_filled_UPDATED (1).xlsx` S1!B9 | Ca_Falkowbacteria.7 (CAN_1_bin.195) | MAG-rename → Patescibacterium.1_m | agree |
| 10 | `cladeexcel.xlsx` S4!B14 | Rhodobacter.13 (CAN_3_bin.10) | MAG-rename → Pseudogemmobacter.1_m | agree |
| 11 | `abundance and coverage (1).xlsx` Sheet1!B7 | Ca_Falkowbacteria.7 | MAG-rename | agree |
| 12 | `Table_S7_MAG_functional_classification.csv` row 7 | Ca_Falkowbacteria.7 | MAG-rename (bin column) | agree |
| 13 | UPDATED (1) S1!B284 (6a04588 era) | Chitinophagaceae.2_m (coasm_bin.61) | drift fix → Chitinophagaceae.14 | agree. This is the era chain of C18 |
| 14 | r13 `design_answers/d5_1_all_members_mag_proxy.tsv` `member` | Azospira.3 (best_mag CAN_4_bin.208) | ASV (data) | agree: r13 members are ASV synthetic genomes |
| 15 | r13 `verify_design.py:16` `MEMB["Ignavibacteria.31"] = ["CAN_3_bin.221"]` | Ignavibacteria.31 | leave, ASV member (data) | agree. The rationale needs the C24 wording |
| 16 | r13 `taxonomy.tsv` rows 36 and 81 | Ignavibacteria.21 / .31 next to their bins | MAG (data) | agree: rebuild in B9 |
| 17 | r13 `kbase_kofam_push/verify/*.json` `output_name` | 51 old IDs | MAG, KBase record, leave as provenance | agree |
| 18 | r13 `diagnostics/rast_only_markers.tsv` `feature` | Rhodobacter.5, Hydrogenedensaceae.8 … (embedded) | MAG (data) | agree: frozen KBase feature IDs, leave |
| 19 | `mag_abundance_by_day_intersection.csv` Kingdom..Species | (MiDAS taxonomy; coasm_bin.365 = Ca_Competibacter / midas_s_34217) | ASV taxonomy, leave (data) | agree, checked |
| 20 | `network/louvain_module_mag_membership.json` `{"iterativeID": "Azospira.3", "MAG": "CAN_4_bin.208"}` | Azospira.3 | regenerate per D2 | partly agree: the field's meaning becomes ambiguous (C2) |
| 21 | `tablesforpaper/guild_flanking_all.csv` `partner` / `partner_is_MAG` | Ignavibacteria.13 | D2 | agree: `partner_is_MAG` is computed by ID identity |
| 22 | `network/network_module_membership_p_value_FDR_allphases.json` | (17 dropped nodes) | leave, ASV context (data) | **disagree**: MAG-only network, regenerate per D2 (C13) |
| 23 | `render_spieceasi_docx.py:402` "Ca. Competibacter.32" | Ca_Competibacter.32 | ASV node (data) | agree |
| 24 | `escher_maps/prototype/*` (untracked) member labels | Azospira.3, Ca_Competibacter.32 | not listed (community-modelling rule: ASV) | agree: ASV members, leave |

**Confirmations from independent checks.**
- The simulated cases are 155/62/74 → 138/70/83 (R220 → R232). 51 IDs change and none collide.
- The MAG-only universe goes 140 → 127.
- SPIEC-EASI goes 53 → 44 kept + 1 joined.
- The 21 old IDs and 4 new IDs that are ASV IDs each have the ASV's `best_mag` equal to the MAG that holds the ID
  (25 of 25).
- No genus-based guild (`gao_pao_categories.json` `norm()` rule) changes between R220 and R232 for any of the 291
  MAGs.
- `mag_abundance_by_day*.csv` taxonomy is MiDAS.
- KBase: the rename log equals the inventory for 291 of 291, the object IDs are contiguous, and there is no collision
  forward or in rollback.

## 9. Files

- `critic/scripts/critic_scan.py` and `critic/scripts/compare_coverage.py` are the scanners. Each is read-only on the
  repo; run them with `~/Documents/py_venv/bin/python -I`.
- The scan outputs are in `critic/scan/`:
  - `hits.tsv`: file × ID × variant, with count and context;
  - `missed_files.tsv`;
  - `variant_hits.tsv`;
  - `skipped.tsv`;
  - `path_hits.tsv` (empty).
