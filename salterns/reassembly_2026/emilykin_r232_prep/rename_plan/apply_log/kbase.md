# KBase stage: GTDB R232 taxonomy and the 51 MAG renames in ws 262319 (2026-11-12)

Writes go to KBase ws 262319 only. Nothing is committed or pushed in git. Token: `~/.kbase/token` (never printed).
Step (d), the Genome re-save, is deferred and was not done: the r13 KOfam annotation-ontology push lives in Genome v2
(GenomeSet 262319/891/1).

## Step 0. Plan re-run against the committed map (read-only)

- `apply_scripts/kbase/build_kbase_plan_committed.py` is a copy of `kbase/build_kbase_plan.py` with two changes:
  - it reads `~/Documents/EmilyKin/mag_iterativeID_r220_to_r232.json` (written by the data stage) instead of
    `results/iterativeid_simulation_r232.tsv`;
  - it writes to `apply_scratch/kbase/plan_rerun/`.
- Run 13:44 CST (read-only whitelist; 891 objects, ws modified 2026-11-09T21:51:02Z, max objid 891), `problems: 0`.
- **Result: identical to the plan in every row.**
  - `assemblies_lineage.tsv`, `renames.tsv`, `sets.tsv`, `genomes_taxonomy.tsv`, `inventory_262319.tsv` and
    `no_action.tsv` are byte-identical to `kbase/`.
  - `preflight_snapshot.json` is equal in every key except `built` (the timestamp).
- Direct check: the committed map equals the simulation for 291/291 bins (old, new, case_old, case_new); 51 changed;
  the rename log's `new_name` equals the map's `old_id` for 291/291.

## Order used inside KBase: (b) → (a) → (c) → (e)

(b) ran first so that its pre-flight could assert `preflight_snapshot.json` verbatim (every object still v1/v2 as
planned). (e) creates a new object, so every max-objid-891 assertion ran before it. (a) looks objects up by objid, so
it does not depend on the names.

## Step (b). Rename 51 Assemblies and 51 Genomes (2026-11-12 19:52–19:53 UTC)

**Script:** new `~/Documents/EmilyKin/scripts/rename_mags_kbase_r232.py` (subcommands `check`, `apply [--apply]
[--limit N]`, `verify`).
- Inputs, all in the repo: `mag_iterativeID_r220_to_r232.json` (old/new ID), `kbase_mag_rename_log.json` (Assembly
  objid from `upa`; the version part is ignored), and r13 `rast/genomes/_manifest.tsv` (Genome objid from
  `genome_ref`, name prefix from `genomeset`).
- It computes the order itself (target held by another renamed object → that object first; Assemblies 1–51, then
  Genomes 52–102 in the same bin order).
- Pre-flight (read-only), run before every apply:
  - Assembly type and `data.assembly_id == bin` (51/51); Genome type and `assembly_ref` objid == the bin's Assembly
    objid (51/51).
  - Every object's current name is its old or new name (C21, resumable by objid).
  - `--snapshot kbase/preflight_snapshot.json`: names, versions and checksums of all 591 snapshot objects unchanged
    (a renamed object may carry its new name), max objid 891.
  - `--plan kbase/renames.tsv`: computed (order, kind, bin, objid, name, new_name) equal to the plan, 102/102.
  - Simulated renames in order against all 891 names (hidden and deleted included): no collision; the end state
    has unique names.
- Each `rename_object` result is checked: same objid, new name, **same version and checksum**. A failed call is
  reconciled by objid and logged.

**Run.**
1. `check --plan … --snapshot …`: 0 problems, 0 of 102 done.
2. Canary `apply --apply --limit 1`: order 1, objid 10 `Ca_Falkowbacteria.7` → `Patescibacterium.1_m`. Checked by
   objid (v1, chsum unchanged) and by name (new name resolves, old name does not).
3. `apply --apply` for the other 101 (pre-flight re-run: 1 of 102 done, 0 problems). Orders 11 then 12 ran as
   planned: objid 195 `Chitinophagaceae.1_m` → `.2_m` before objid 150 `Flavipsychrobacter.1_m` →
   `Chitinophagaceae.1_m`; Genomes 744 before 785 likewise.
4. `verify`: **582/582 (objid, name) pairs match the R232 map** (291 Assemblies, 291 Genomes); 0 duplicate names in the
   workspace.

**Record:** `~/Documents/EmilyKin/kbase_r232/rename_mags_kbase_r232_log.tsv`, 102 rows, all `renamed`. Its
(order, objid, old name, new name, version) equal `kbase/renames.tsv` row for row.

**Not changed, by design (PLAN §3b):** Genome `data.id`, `source_id` and every feature/CDS/mRNA ID keep the old ID
string (e.g. `…code11_b5_Chitinophagaceae.1_m.CDS.*` inside coasm_bin.117's Genome, now named `…_Chitinophagaceae.2_m`).
Join on bin or objid, never on those strings. Rollback: `rename_object` in reverse order with the old names.

## Step (a). R232 lineage onto the 291 Assemblies, v1 → v2 (20:01–20:05 UTC)

**Script:** new `~/Documents/EmilyKin/scripts/writeback_gtdbtk_r232_kbase.py`, an adapted copy of
`salterns/reassembly_2026/local/writeback_gtdbtk_r232.py` (the salterns script is unchanged). Subcommands `apply`,
`verify`, `sets` (step c) and `report` (step e). Adaptations:
- `WSID = 262319`; `OUT = ~/gtdbtk_kbase/emilykin_r232/gtdbtk_out`; `BATCH = ~/gtdbtk_kbase/emilykin_r232/batchfile.tsv`.
  `run_info()` checks `gtdbtk.json` is 2.8.0 / r232, all steps completed, over this batchfile.
- `summaries()` reads `classify/gtdbtk.bac120.summary.tsv` (the ar53 file is optional; there is none here).
- `mags()` replaced: all 291 MAGs (no CheckM gate). Classified MD5 = KBase-style MD5 of the batchfile FASTA
  (md5 of `','.join(sorted(md5(upper(contig))))`) plus the contig count. Bin sets of batchfile, summary, R232 map,
  rename log and the R220 summary must be equal (291).
- `infos_by_name()` replaced by lookup **by objid** (from `kbase_mag_rename_log.json` `upa`). Before saving, each
  object must have: Assembly type; name ∈ {old ID, new ID}; `data.assembly_id == bin`; `data.md5` and info `MD5` ==
  classified MD5; `num_contigs` / `N Contigs` == FASTA contig count; latest version v1. The saved version must have
  the same objid and name, version +1, the GTDB_* metadata, and unchanged MD5 / Size / N Contigs.
- `CLASSIFIED_BY` is byte-identical to the salterns string, so `carries()` is unchanged.
- Provenance as salterns, with the description "… over the 291 dereplicated EmilyKin MAGs (best-of-both catalogue,
  git 6541627) …". `custom` holds the GTDB-Tk details plus `classified_assembly_version` (v1), `bin`,
  `iterativeID_r232`, `iterativeID_r220`, `previous_gtdb_r220` (from `git show 9e087c3:gtdbtk.bac120.summary.tsv`,
  which equals the prep `r220_vs_r232_per_mag.tsv` for 291/291) and its source note.
- Log in the repo: `~/Documents/EmilyKin/kbase_r232/writeback_gtdbtk_r232_log.tsv` (not the salterns `writeback/`).

**Run.**
1. Dry run: 291 to write, 0 skipped (all 291 MD5 / contig counts equal the classified FASTA, all at v1), 14 batches.
2. Canaries `--limit 2`: CAN_1_bin.107 → 262319/3/2, CAN_1_bin.148 → 262319/7/2. Checked: v2 data equal to v1 in every
   field except the added `std_lineages`; metadata `GTDB_lineage` / `GTDB_source_ver 232.0` / `GTDB_taxon_id` /
   `GTDB_classified_by` present; MD5, Size, N Contigs unchanged.
3. `apply --apply`: the other 289 in 2 single saves plus 12 batches of ≤25, no stop.
4. `verify`: **291/291** carry the R232 metadata; 291/291 `std_lineages.gtdb` + `assembly_id` + `md5` read back
   exactly; 291/291 MD5 and contig count equal the classified FASTA; 291/291 named with the R232 ID; all at v2.
5. Log vs plan `assemblies_lineage.tsv`: 291 rows; lineage, new ref (`objid/expected_new_ver`) and MD5 equal for 291/291.

## Step (c). AssemblySet 262319/585 re-saved, v2 → v3 (20:08 UTC)

`writeback_gtdbtk_r232_kbase.py sets --plan kbase/sets.tsv --apply`.
- Each item → the latest version carrying the R232 lineage (`carries()` by objid; all v2). Label = the object's current
  name, which must equal the R232 ID (so (b) had to be complete). Items sorted by label (Python sort;
  `alphaI_cluster.1` last).
- Checked against the plan before saving: items equal `sets.tsv` for 291/291 (old index, new index, label, ref).
- Description: "291 dereplicated MAGs from the EmilyKin denitrification bioreactor 291-MAG catalogue, named by
  iterativeID after the GTDB R232 update (R220 -> R232 map: mag_iterativeID_r220_to_r232.json; bin -> object:
  kbase_mag_rename_log.json). GTDB R232 (GTDB-Tk 2.8.0) lineages in std_lineages.gtdb of each Assembly version listed."
- User metadata: `n_assemblies` 291, `naming` "iterativeID (GTDB R232, 2026-11)". (`description` and `item_count` are
  extracted by the type, `@metadata ws`.) Provenance `input_ws_objects` = `[262319/585/2]`.
- Result 262319/585/3: 291 items, 51 labels changed, all refs v2; read back equal to the plan.
- GenomeSets: untouched (refs, no labels; 891 pins Genome v2 and only needs re-pointing if step (d) is ever done).
- Record: `~/Documents/EmilyKin/kbase_r232/sets_saved.tsv`.

## Step (e). Report GTDBTk_R232_local_report_emilykin = 262319/892/1 (20:10 UTC)

`writeback_gtdbtk_r232_kbase.py report --apply` (staging under `$TMPDIR` = `apply_scratch/kbase/`; a first dry run
had staged inside the repo's `kbase_r232/`, which I deleted and moved out of the repo in the script).
- 12 files (10.0 MB) + HTML: `emilykin_MAGs_gtdb_r232.tsv` (291 rows: bin, Assembly v2 ref, RAST Genome ref, R232 / R220
  ID, cases, R232 / R220 lineage, GTDB-Tk details), `emilykin_MAGs_iterativeID_r220_to_r232_changes.tsv` (51 rows:
  objids, old → new ID, old → new Genome name, lineages), the GTDB-Tk summary, 8 trees, and a zip of `gtdbtk_out`
  (22 files; intermediates and the 342 MB public reference MSA omitted).
- HTML: the 51 renames, lowest named rank, phyla, and the per-MAG table.
- `objects_created` = 291 Assembly v2 + set 585/3 (292). All 13 handles resolve. Record:
  `~/Documents/EmilyKin/kbase_r232/report_saved.json`.

## Final KBase state (read-only check, after (e))

Against `kbase/inventory_262319.tsv`: 892 objects, max objid 892; 291 Assemblies at v2 with their planned names; set
585 at v3; 291 Genomes, 8 GenomeSets, 299 reports and the narrative unchanged (version and checksum; 51 Genomes with
their new names); the only new object is 892, the report. Step (d) not done (deferred).

## Repo side (~/Documents/EmilyKin; nothing committed)

Every table edit below is keyed by bin (`scripts/_binkeyed.patch_table`, which keeps each file's terminator and rewrites
only changed lines) or is a builder re-run gated against HEAD. No string replacement.

| file | what | how | numstat (add/del) | before → after |
|---|---|---|---|---|
| `kbase_mag_rename_log.json` | `new_name` = R232 ID; `upa` = Assembly v2 (the lineage-carrying version); history fields `name_r220` (the R220-era name) and `upa_r220` (the v1 ref the log held; v1 has no lineage and is what the RAST Genomes' `assembly_ref` pins) for all 291 | `json.dumps(indent=1, sort_keys=True)`, no trailing newline (round-trip of the HEAD bytes proved first); v2 refs from the writeback log | 924/342 | new_name changed 51; upa 291 (v1 → v2); LF; 24,014 → 43,214 bytes |
| r13 `rast/genomes/_manifest.tsv` | `iterativeID` and `genome_name` (= `genomeset` + `_` + R232 ID) for the 51; `genome_ref` / `assembly_ref` unchanged (a rename keeps versions) | patch_table on `bin`; asserted the old values first | 51/51 | 51 rows; other columns unchanged 291/291; LF |
| r13 `rast/check_fetch.py` | line 37: compare the rename log's `upa` with the manifest's `assembly_ref` by ws/objid, not as a full string (with `upa` → v2 the old check flags all 291; the plan's "reads only the objid" was wrong for this line) | Edit, comment added | 3/2 | post-check below |
| r13 `taxonomy.tsv` | B9a hybrid (gate verdict tier C): `build_taxonomy.py` (FLAG/DIDERM fix already in from the data stage) run with `--out` to scratch, then `kbase_rast_genome_ref`, `rast_genomeset`, `envelope_evidence_rast` restored to HEAD by bin, then installed | script + restore | 291/291 | see gate below |
| r13 `panel/mag_panel.tsv` | B9b: `build_panel.py --out <scratch>` (gate: tier A at HEAD), only `mag_panel.tsv` differed, installed | copy | 291/291 | see gate below |
| r13 `calls/summary.tsv` | B9e: `gtdb_phylum` / `gtdb_order` from the R232 summary (`--reflag` does not refresh them; flags do not depend on taxonomy) | patch_table on `bin` | 37/37 | phylum 22, order 19 cells; other columns unchanged |
| r13 `panel/PANEL_FROZEN.md` | dated note under the title + section "Freeze 2a: GTDB R232 relabel (2026-11-12)" at the end: why, the rebuild check, new hashes of `mag_panel.tsv` and `tests/test_build_panel.py` (backtick row format), changed-input hashes. Freeze-2 tables untouched | append | 39/0 | see below |
| `scripts/rename_mags_kbase.py` | docstring warning: do not run with a new map; use the two new scripts | insert | 5/0 | parses |
| new `scripts/rename_mags_kbase_r232.py` | step (b) script | – | new | – |
| new `scripts/writeback_gtdbtk_r232_kbase.py` | steps (a), (c), (e) script | – | new | – |
| new `kbase_r232/` | run records: `rename_mags_kbase_r232_log.tsv` (102), `writeback_gtdbtk_r232_log.tsv` (291), `sets_saved.tsv` (1), `report_saved.json` | written by the scripts | new | – |

**B9 gates.**
- `taxonomy.tsv` vs HEAD by bin: differing columns are only `kbase_assembly_ref` (291, v1 → v2, from the log), `iterativeID`
  (51) and the GTDB-derived columns (`gtdb_classification` 130, phylum 22, class 23, order 19, family 24, genus 60,
  species 43, classification method / note / closest-genome / placement / RED / msa / translation-table / warnings, all
  equal to the R232 summary). **Equal to HEAD 291/291:** `template_class`, `gram_stain_inferred`, `envelope_lineage`,
  `gram_flag`, `lps_expected`, `genetic_code_used`, `gram_rule`, `envelope_evidence_kofam(_basis)` and the 3 restored
  RAST columns. Phylum renames: Patescibacteria → Patescibacteriota 21, Myxococcota → Myxococcota_A 1; both keep their
  HEAD envelope class through the FLAG/DIDERM fix. Run without `--kofam`, as HEAD was (all 291 `NA`).
- `build_panel.py` (1 min 35 s, 2.38 GB RSS): `mag_panel_loci.tsv`, `member_panel.tsv`, `member_panel_long.tsv`
  byte-identical; `build_full.log` identical; `build_report.json` differs only in `inputs.out` (scratch path) → committed
  file kept. `mag_panel.tsv` differs only in `iterativeID` 51, `kbase_ref` 291 (→ v2), `gtdb_phylum` 22, `gtdb_class`
  23, `gtdb_family` 24, `gtdb_genus` 60; every marker / count / tier / class / flag column identical.
- `validate_panel.py`: output identical to the committed `validate_full.txt` (18 tests, PASS 16, WARN 2, FAIL 0); not
  written over.
- `pytest panel/tests/` (TMPDIR in scratch, no bytecode, no cache): 15/15 pass.
- `check_fetch.py` (fixed), output to scratch: `QC problems: 0`; vs the committed `check_fetch.out` only the two phylum
  label lines differ (Patescibacteriota, Myxococcota_A). The HEAD version of the script, on the same files, reports 291
  "!= rename upa" problems, so the fix was needed. `check_fetch.out` is left as the historical record of the fetch.
- PANEL_FROZEN: running only `compute_answers.py`'s regex line over the new file, the 7 files it checks all match their
  current sha256, and all 16 backticked hash rows match the files on disk. `compute_answers.py` itself was not run (it
  writes outputs).
- D1 gate on the r13 ID columns: `taxonomy.tsv`, `mag_panel.tsv` and `_manifest.tsv` `iterativeID` equal the R232 map for
  291/291 and are unique; manifest `genome_name` = `genomeset_<R232 ID>` 291/291; rename log `new_name` 291/291 R232,
  unique, `upa` all v2, `upa_r220` all v1.

**Left unchanged on purpose.**
- `kbase_kofam_push/placement_by_bin.tsv`, `payloads/`, `verify/`, `push_log.tsv`, `plan_r13/kbase_rast/rast_jobs.json`,
  `rast_inputs.json`, `kbase_mag_upload_log.json`, `kbase_mag_upload.out`, `check_fetch.out`, `fetch_rast_full.log`:
  records of what was pushed or fetched (D5); they keep the old names. `rast/genomes/<bin>.tsv`, `crosswalk.tsv`,
  `concordance_pairs.tsv`: frozen KBase feature IDs (C23).
- `scripts/apply_mag_ids.py`: the data stage documented that the r13 tables are rebuilt by their own builders after the
  KBase rename and are not listed there; I kept it that way.
- `panel/build_report.json`, `validate_full.txt`, `build_full.log`: equal to a rebuild (apart from the scratch path).
- C24: the r13 Track-2 backing of ASV member `Ignavibacteria.31` by CAN_3_bin.221 (`member_mag_choice.tsv`,
  `verify_design.py`) is a design choice, not a MAG ID; untouched.
- C25 (leftover-reads outputs) was handled by the data stage (`leftover_reads_analysis/R232_NOTE.md`).
- The plan builder `kbase/build_kbase_plan.py` now fails its own assertion (`rename log != committed ids`) because the
  log holds the R232 names. That is expected; it described the pre-rename state.

## Final checks

- ASV-keyed files: 1,443/1,443 md5-identical to `asv_md5_before.tsv`.
- `community_modeling/plan_r13/integrator/AMENDMENTS.sha256` (run from `community_modeling/`): 14/14 OK.
  `track2_guards_sha256_before.txt`: 7/7 OK.
- Line endings: none of my 8 modified tracked files changed terminator (all LF, as at HEAD). Data tables have
  added = deleted in numstat; the rename log (fields added) and the two text files (append-only) are the expected
  exceptions.
- Bytecode: my heredoc imports of `_binkeyed` wrote `scripts/__pycache__/_binkeyed.cpython-312.pyc` (gitignored);
  deleted. No other `.pyc` newer than the stage start.
- **Concurrent writer seen.** Between 13:49 and 13:51 CST, files I did not touch changed in the repo:
  `scripts/render_fig1_tree.py`, `fig1_tree/*` (annotation, MAG table, tree PNG/PDF, `compare_trees_changed_splits.tsv`,
  `fig1_mags.txt`), `clade_classify/organism_tree/*`, `clade_classify/out/nosz_top5_tree.*`,
  `figures/percell_composition.*`, `gene_ab_figure/*` figures and `gene_abundance_over_time.tsv`, `ppk1_classify/*`.
  Presumably the figures stage running in parallel. Not mine; not inspected or reverted.
- Nothing committed or pushed.
- **Freeze-2 hash pinning.** `grep` over `community_modeling/analysis/r13` and `plan_r13` for the freeze-2 sha256 of
  `mag_panel.tsv` (ad22880f…) and `tests/test_build_panel.py` (952ad10f…): no file other than `PANEL_FROZEN.md` holds
  them. `git grep PANEL_FROZEN` hits only prose (`ANNOTATIONS.md`, `DESIGN_ANSWERS.md`, `FINALIZE_REVIEW_LOG.md`,
  `crosswalk_NOTES.md`, `panel_NOTES.md`, `plan_r13/track2/GENOME_OBJECT_SPEC.md` and its two `_as_of_AMENDMENT_*`
  copies) and `compute_answers.py` (handled by freeze 2a). No Track-2 code or JSON pins the superseded hashes.

## Scripts and scratch of this stage (outside the repo)

- `apply_scripts/kbase/build_kbase_plan_committed.py`: the step-0 plan re-run (read-only).
- `apply_scripts/kbase/probe_readonly.py`: read-only probe of set 585, one Assembly and the Assembly type spec.
- `apply_scratch/kbase/`: `plan_rerun/` (step-0 outputs), `rename_check_1.log`, `rename_apply.log`,
  `rename_verify.log`, `writeback_dryrun.log`, `writeback_apply.log`, `writeback_verify.log`, `report_apply.log`,
  `emilykin_r232_report_stage_*/` (the files uploaded in the report), `final_state_check.log`, and `b9/` (taxonomy
  HEAD / run / hybrid copies and gate log, `panel_out/` rebuild, validator and check_fetch outputs, pytest tmp).
