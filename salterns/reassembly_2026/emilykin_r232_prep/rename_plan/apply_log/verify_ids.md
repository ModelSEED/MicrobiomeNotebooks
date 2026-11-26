# verify_ids: independent check of bin-keyed MAG iterativeIDs (GTDB R232 apply)

Run 2026-11-12, after the gate, data, kbase, figures and docs stages. I checked the actual files in `~/Documents/EmilyKin`
(uncommitted working tree on `iterativeid-nomenclature`, HEAD 9e087c3) and KBase ws 262319 (read-only calls only:
`list_objects`, `get_objects2`). I did not use the stage reports as evidence.

**Result: no ID error found, so no file was changed by this stage.** Nothing in the repo or in KBase was written.
The open items below need a decision; none is a mechanical fix.

The scripts and outputs are in `rename_plan/apply_scratch/verify_ids/`. They are listed at the end.

## Reference map (independent source)

- The reference is `results/iterativeid_simulation_r232.tsv`. Its `committed_id` column is the R220 ID and
  `simulated_r232_id` is the R232 ID.
- 291 bins. 51 IDs change, 50 of them at the root. No duplicates in the old or the new IDs.
- Cases: R232 138/70/83, R220 155/62/74.
- `committed_id` equals HEAD `mag_iterativeID_old_to_new.json` for 291/291.

## 1. Map files (D1): all pass

| file | check | result |
|---|---|---|
| `mag_iterativeID_r220_to_r232.csv` / `.json` | bin, old_id, new_id, case_old, case_new equal the simulation | 291/291 |
| `mag_iterativeID_old_to_new.json` | equals the R232 map | 291/291, no extra keys |
| `mag_iterativeID_old_to_new.csv` | `new_mag_iterativeID` and `assignment_case` equal R232 | 291/291 |
| `mag_iterativeID_old_to_new.csv` | `old_mag_iterativeID` and `old_intersection_id` (pre-c796fdf columns) equal HEAD | 291/291 |
| `scripts/_new_ids.json` | `mag_new` and `mag_case` equal R232 | 291/291 |
| `scripts/_new_ids.json` | `asv_old2new`, `asv_hash2new`, `asv_new2hash` byte-equal to HEAD | 3,945 each |
| `rename_plan/data/mag_id_map_r220_r232.json` | equals the simulation | 291/291 |

## 2. Every MAG ID column, joined on bin

**Method (`v03_scan.py`).**
- Scope: every tracked and untracked (non-ignored) csv/tsv/txt/json/xlsx file that contains a bin name. That is
  611 files.
- A bin column is any column in which at least half the non-blank values are bin names.
- Every other column that holds an ID from any map era is joined row by row on its best-matching bin column. The ID
  sources are: R220, R232, the `alias_holders.tsv` aliases, the 924fa18 / 6a04588 / ffe668c IDs, and the
  pre-c796fdf columns.
- JSON is walked the same way: dicts keyed by bin, and records that carry a bin value.
- `v07_json_multi.py` checks JSON records that carry two bins: `source_bin` / `source_MAG_iterativeID`, and the same for
  `target`.

**Living MAG ID columns: 192 columns in 48 files, all equal to the R232 map for every row, with no stale and no
foreign ID.** The full list is in `verified_mag_id_columns.tsv`. It covers:
- every ADD and already-listed row of `data/apply_mag_ids_extension.tsv`;
- the two workbooks the docs stage edited, sheets S1/S1b/S2a/S3/S4/S5/S6/S7;
- `cladeexcel.xlsx` S4 and `abundance and coverage (1).xlsx` Sheet1;
- Table S7 and `mag_guild_classification.csv` (D4: the 33 stale older-map IDs are gone);
- the r13 `taxonomy.tsv`, `mag_panel.tsv` and `_manifest.tsv`;
- `phase_rpkm_chord/*.csv` (the C19 `Root_m.N` IDs are gone);
- `fig1_tree` tables, gene_ab_figure tables, ko_rpkm tables, `mag_percell_abundance.csv`;
- the network outputs' `*_MAG_iterativeID` fields, `louvain_module_mag_membership.json` `mag_iterativeID`
  (641 members), and `pao_gao_shared_driver_summary.json` (28 records);
- `kbase_mag_rename_log.json` `new_name`.

**Uniqueness (D1/C18 gate).** In every living MAG ID column each bin carries one ID, and no ID is on two bins of the
291. There is one exception, which predates R232 and is listed under open items: `CAN_4_bin.225`.

**Other workbook cells.** I scanned every string cell of the 4 living workbooks, including notes and comments.
- The only old-ID string is `Chitinophagaceae.1_m` on CAN_5_bin.149. That is its correct R232 ID (the swap chain).
- No "ID next to bin" text pair disagrees with R232.

**Flagged columns, all expected (`scan_tracked.tsv`).**
- **ASV node and member columns: leave.** `source`, `target`, `node`, `name`, `asv_iterativeID` and `guild_node` in
  the network tables, the Louvain `iterativeID` field, r12/r13 `member` and `rep_iterativeID`, and `asv_maps.json`
  `iid`. Under D2 they hold ASV IDs.
- **The map files' own old columns.** `old_id`, `name_r220`, `kbase_r232/rename_mags_kbase_r232_log.tsv` `old_name`,
  and `network/mag_guild_iterative_ID_stale_frozen.tsv`.
- **Superseded or do-not-edit files.** These are `*.bak.*`, `fig1_tree/superseded_*`, `test62_*`,
  `server_handoff_out/`, `MAG_target_gene_empty_tables.xlsx`, `MAG_target_gene_filled (1).xlsx` and
  `MAG_target_gene_filled_UPDATED.xlsx`.
- **Historical data-inconsistency and r12 records.**

## 3. Network MAG link (D2)

- **`network/mag_linked_asvs_r220_frozen.tsv`** holds 140 (ASV, bin) pairs.
  - Each ASV equals the bin's R220 MAG ID (140/140), and all 140 are ASV IDs.
  - The set equals the node set of the committed all-phases network (140 = 140).
- **All 9 consumers read the frozen list for membership and the map only for the label.** I checked this by grep.
  - The 9 are: `rcg.load_abundances` (and through it the all-phases, per-phase and GAO/PAO scripts),
    `build_module_mag_membership`, the two module panels, `render_spieceasi_mag_network`, both
    `spieceasi_modules_flanking`, and `render_spieceasi_docx`.
  - No consumer reads `mag_iterativeID_old_to_new.json` for membership.
- **Network tables: (node, bin) equals the frozen pair and the MAG ID field equals R232[bin] for every row.**
  - Files: the 6 `out/asv/spieceasi_*.csv`, the 5 `tablesforpaper` outputs and the Louvain JSON.
  - `partner_is_MAG` equals membership of the frozen list.
- **Each of these files equals HEAD once the added ID fields are dropped.** The Louvain JSON also differs in 6
  `gamma_rule` provenance strings and `_description` / `name_map`, which the gate predicted.

## 4. Old IDs left in MAG contexts

The 51 old IDs split into three groups:
- 29 old MAG-only `_m` IDs that are no longer any bin's ID;
- 1 reused ID, `Chitinophagaceae.1_m`;
- 21 inherited IDs.

**Tracked and untracked text** (`git grep`, boundary-checked). An old MAG-only `_m` ID is left in only these files.
- **Historical records (D5).** The `kbase_kofam_push` placement and verify records, `rast_inputs.json`,
  `fetch_rast_full.log` and the D54/D59 data-inconsistency records.
- **Frozen KBase feature IDs (C23).** `rast_only_markers.tsv` and `rast_only_marker_cds.tsv`.
- **Superseded files.** `fig1_tree/superseded_*`, `test62_*`, `server_handoff_out/` and `ko_copy_number_matrix_nosZ.pre_bb.bak.csv`.
- **The map files and the KBase r232 logs.**
- **The pre-c796fdf columns of `mag_iterativeID_old_to_new.csv`**, unchanged from HEAD.
- **Old-to-new change notes** in `fig1_tree/README.md` and in a docstring of `scripts/rename_mags_kbase_r232.py`.

**Critic scanner re-run over 39k tracked, untracked and ignored files** (`critic_rescan/`, same code with the output sent
to scratch).
- 143 files still hold an old MAG-only `_m` ID. Every one is in a category above, or is:
  - a backup;
  - an r13 snapshot (`_freeze1`, `_prefix_full`, `probe`, `rast/_probe`, `review/`);
  - a KOfam push payload;
  - a `rast/genomes/<bin>.tsv` feature table.
- In the two edited DOCX the hits are only in `word/comments.xml` (explanations) and in `<w:delText>` (tracked
  deletions).

**DOCX accepted view** (`v06_docx.py`: `w:del` / `w:moveFrom` removed).
- **Main DOCX.** Table 1 has 53 rows with a bin, and all 53 have the R232 ID. No old MAG-only or inherited ID is in the
  accepted text.
- **SI DOCX.** No old ID.
- **SPIEC-EASI DOCX.** The inherited IDs appear only as ASV node labels (D2): `Azospira.3` ×5,
  `Ignavibacteria.13` ×1 and `Ignavibacteria.5` ×1. Both MAG mentions read "Azospira.3 [MAG CAN_4_bin.208,
  Proximibacter.2_m under GTDB R232]".
- **Comments.** All 25 "bin old → new" pairs in the main DOCX comments equal R232.

**"ID next to bin" in text** (`v05_textpairs.py`, `v05b_tabpairs.py`). Every mismatch is one of these:
- an r13/r12 member (ASV) name beside its backing bin (C24);
- a historical record;
- a network ASV node label beside its bin (D2).

None is a living MAG label.

**Inherited old IDs in edited text and code.**
- Only change notes ("X before GTDB R232", "X → Y") and ASV node labels (the SPIEC-EASI docx generator).
- Every "old → new" pair written in edited notes equals the map (12 pairs).
- The `paper3/R232_change_list_for_Emily.md` bin/ID statements agree with the map. Examples: coasm_bin.260 =
  Ca_Sarcinithrix.1, CAN_5_bin.70 = Azonexus.12, coasm_bin.157 = Steroidobacterales.1_m.

**PDFs** (`pdftotext`, 45 PDFs).
- The regenerated PDFs carry no old MAG ID.
- The SPIEC-EASI PDFs keep inherited IDs as ASV node labels (D2).
- `denitrification_MAG_heatmap.pdf` is stale; see open items.

**Prefix traps.** I counted each trap ID per changed file against HEAD:
- `Rhodobacteraceae.10`–`.19`, `Simplicispira.1` versus `.1_m`, and the 4 new IDs that are also ASV IDs.
- The count changes only in bin-keyed MAG columns.
- `Rhodobacteraceae.1X` never changed.

## 5. ASV-keyed files and historical records

- **ASV baseline.** `asv_md5_before.tsv` matches for 1,443/1,443 files, so no ASV-keyed baseline file changed. The
  files that D2 changed are not in that baseline. They are ASV-node-keyed and gained bin and R232 MAG ID fields only:
  - `network/compositional_networks/out/asv/` (6 files);
  - `network/louvain_module_mag_membership.json`;
  - `tablesforpaper/{guild_flanking_all, guild_flanking_nonMAG, modules_membership_89node, pao_gao_direct_edges}.csv`
    and `pao_gao_shared_driver_summary.json`.
- **Tracked files.** 109 changed against `tracked_md5_head.txt`, and none is missing. No changed file is a log, run
  record, `*.pre_*.bak.*`, push payload, verify JSON, KBase upload log, r12 snapshot or track2_guards file.
  - `kbase_mag_rename_log.json` changed on purpose (C10).
  - `kbase_mag_upload_log.json` is unchanged.
- **Files newer than the gate stamp** (`newer_files.txt`), apart from the files git reports as changed or new:
  - The tracked ones are rewrites that match `tracked_md5_head.txt`, so their content is unchanged.
  - The ignored ones are regenerated figure outputs whose content I did not check. I text-scanned the PDFs with
    `pdftotext` (section 4); the PNGs cannot be scanned.
  - No r13 snapshot, payload or log is among them.
- **New backups.** The 6 new `*.pre_r232.bak.*` files are byte-identical to the HEAD versions of their originals.
- **AMENDMENTS.sha256**, run from `community_modeling/`: 14/14 OK.
- **track2_guards sha256**: 7/7 OK, covering AMENDMENT_1–4, CRITERIA.md, run_plan.json and acceptance.json.
- **Line endings.** No changed text file changed its line endings.
- **Row counts.** In `git diff --numstat` every csv/tsv has added equal to deleted.

## 6. KBase ws 262319 (read-only, joined on bin through `data.assembly_id`)

- **Objects.** 892 objects with no duplicate names: 291 Assemblies, 291 Genomes, 8 GenomeSets, 1 AssemblySet, 300
  reports and the narrative.
- **Assemblies.** 291/291 are named with the R232 ID of their own bin. All are at v2, and each carries
  `std_lineages.gtdb` equal to the R232 summary lineage (source_ver 232.0).
- **Genomes.** 291/291 are named `<genomeset>_<R232 ID>`, joined through `assembly_ref` → Assembly objid → bin.
  - The script's regex is `emilykin_MAGs_RAST_code\d+(?:_b\d+)?_(.+)`, so the 5 code25 Genomes are included.
  - The first run used a code11-only regex and flagged those 5. The re-run with the regex above found 0 bad.
- **AssemblySet 585** is at v3. Its 291 item labels equal R232[bin], and every ref is at v2.
- **GenomeSets.** The 8 sets have elements keyed by ref and contain no old ID string.
- **Names.** No object is named with an old MAG-only ID.
- **Rename log** (`kbase_r232/rename_mags_kbase_r232_log.tsv`). 51 Assembly and 51 Genome renames, all on the 51
  changed bins. Old name = R220[bin] and new name = R232[bin] for all 102.

## Open items (need a decision; not changed)

1. **Pre-existing, not caused by R232: `CAN_4_bin.225` shares an ID with `CAN_4_bin.64`.**
   - CAN_4_bin.225 is a 276-era bin that is not among the 291 (SERVER_TASKS_291.md: "must be dropped everywhere").
   - It is still labelled `Accumulibacter.1_m`, which is CAN_4_bin.64's ID, in 6 tables:
     - `clade_classify/out/allbins_nosz_master.tsv` (1 row);
     - `gene_ab_figure/data/rpkm_long_MAG_gene_sample_phase.tsv` (225 rows);
     - `ko_gene_abundance/ko_rpkm_combined{,_ext}.tsv` (1 row each);
     - `ko_rpkm_per_sample.tsv` (42 rows) and `ko_rpkm_per_sample_ext.tsv` (45 rows).
   - These rows are identical at HEAD. Dropping or blanking them changes data, so it is a judgment call.
2. **Figures still showing old MAG IDs in a MAG context, not re-rendered per the gate verdicts:**
   - the `phase_rpkm_chord` chord and heatmap PNGs (their CSVs are bin-patched);
   - `denitrification_MAG_heatmap.png` and its untracked `.pdf` (Giesbergeria.1_m, Gammaproteobacteria.1_m,
     Azospira.3, Ca_Competibacter.32, Ignavibacteria.21/.5);
   - the ID-labelled Figure 4 image inside `papert3_070626 (2)_correction.docx`, which has a bin-keyed relabel comment
     and no generator.
3. **`ko_gene_abundance/ko_rpkm_root_{combined,per_sample}.tsv` (B4)** are keyed by MAG-ID roots from the 924fa18 era
   and were not regenerated (gate tier C). They need approval to re-run.
4. **KBase Genomes still pin Assembly v1** in `assembly_ref`, which carries no lineage. Step (d) was deferred.
5. **`community_modeling/plan_r13/review_feasibility/FEASIBILITY_REVIEW.md:90`** says "Ignavibacteria.21 for
   CAN_3_bin.221 … are MAG iterativeIDs". It is a review record (D5), so I left it.

## Files (rename_plan/apply_scratch/verify_ids/)

- `ref_map.json` and `id_sets.json`: the reference map and the old/new ID groups.
- `v01_maps.py`: the map files.
- `v02_md5.py`: the ASV md5 baseline.
- `v03_scan.py`: the bin-join column scan. Outputs `scan_tracked.tsv` and `verified_mag_id_columns.tsv`.
- `v04_kbase.py`: the read-only KBase check.
- `v05_textpairs.py` and `v05b_tabpairs.py`: "ID near bin" pairs. Outputs `textpairs*.tsv`.
- `v06_docx.py`: the DOCX accepted view.
- `v07_json_multi.py`: JSON records with two bins.
- `critic_scan_rerun.py` and `critic_rescan/`: the critic scanner re-run.
- `changed_tracked.txt` and `newer_files.txt`: the changed-file lists.
