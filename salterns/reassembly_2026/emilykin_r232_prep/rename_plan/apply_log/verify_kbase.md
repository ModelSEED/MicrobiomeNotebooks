# Independent verification, lens: KBase (2026-11-12, 15:05–15:25 CST)

**Verdict: PASS. Nothing needed fixing.** 57 checks pass with 0 failures: 43 on KBase, 14 on the repo, plus the
cross-cutting end checks below. 13 NOTEs are recorded; none is a defect.

**Scope of what I did:**
- KBase: read-only calls only. The whitelist in `kbro.py` allows `get_workspace_info`, `list_objects`,
  `get_object_info3`, `get_objects2`, `list_referencing_objects`, `get_object_history` and
  `AbstractHandle.hids_to_handles`, plus blob-node record GETs.
- Repo: I edited no file in `~/Documents/EmilyKin`.
- Hygiene: everything ran under nice 10 with `python -I -B`, and nothing was committed.

## Independence: where the expected values came from

I took no expected value from a file the kbase stage wrote.

| expected value | source |
|---|---|
| R232 ID per bin | prep `results/iterativeid_simulation_r232.tsv` (not the committed map, not `renames.tsv`) |
| R232 lineage per bin | prep `results/gtdbtk.r232.summary.tsv`. Columns 1–2 are byte-equal to `~/gtdbtk_kbase/emilykin_r232/gtdbtk_out/classify/` |
| R220 lineage | prep `results/r220_vs_r232_per_mag.tsv` |
| bin → Assembly objid | `git show 9e087c3:kbase_mag_rename_log.json` |
| Genome → bin | KBase itself: Genome `data.assembly_ref` → Assembly objid → `data.assembly_id` |
| pre-run state | `kbase/inventory_262319.tsv` (891 objects), `preflight_snapshot.json`, `no_action.tsv`, all built 2026-11-11 before any write |
| classified MD5 | recomputed here from the 291 batchfile FASTAs with the KBase recipe |

`renames.tsv`, `sets.tsv` and `kbase_r232/*` were then compared with that truth as consistency checks.

## KBase checks (`apply_scratch/verify_kbase/results_kbase.tsv`, log `verify_kbase_run3.log`)

**Census of ws 262319** (`list_objects` with hidden, deleted and all versions included):
- Max objid is 892, objids 1–892 are contiguous, and there are 0 deleted objects.
- The only object not in the inventory is 892, `GTDBTk_R232_local_report_emilykin` (Report, v1).
- For all 891 inventory objects, the inventory version still exists with the same checksum, save date and type.
- Versions:
  - every Assembly has exactly v1 and v2;
  - set 585 has exactly v1, v2 and v3;
  - no other object has a new version. That covers all 291 Genomes (still v2), 8 GenomeSets (v1), 299 reports and the
    narrative (v5).
- All 591 `preflight_snapshot.json` (version, checksum) pairs are unchanged.
- The 291 Assembly v2 were saved 20:01:34–20:07:10 UTC, set 585/3 at 20:08:26 and the report at 20:10:39.
- All 892 names are distinct.

**Assemblies, all 291 (v1 and v2 data fetched in full):**
- `data.assembly_id` equals the bin in both versions.
- v1 has no `std_lineages`.
- **v2 data equals v1 in every field except the added `std_lineages`** (the stage checked this only on its 2 canaries).
- `std_lineages` holds only the key `gtdb`, which equals `{lineage: prep R232, source_ver: "232.0", taxon_id: lowest named rank}`.
- MD5 agrees everywhere: v1 data, v2 data, v1 metadata, v2 metadata and the recomputed FASTA MD5.
- The contig count is the same in v1 and v2 and equals the FASTA count.
- The type is still Assembly-5.1.
- v2 metadata equals v1 metadata plus exactly the 4 GTDB_* keys, with the expected values.
- Provenance: `input_ws_objects` is the v1 ref, and `custom` holds bin, R232 ID, R220 ID and `previous_gtdb_r220`, which
  equals the prep R220 lineage.
- The name (on both versions) equals the sim R232 ID.
- 130 of the 291 lineages differ from R220, as the plan says.

**Renames (by objid):**
- **Assemblies:** all 291 carry the sim R232 ID. Exactly 51 objects were renamed, and they are exactly the 51 changed bins.
- **Genomes:**
  - all 291 map one-to-one to the bins through `assembly_ref`, and every `assembly_ref` pins Assembly v1;
  - every name equals the pre-run prefix plus the R232 ID; 51 were renamed;
  - `data.id` keeps the pre-run name (frozen, by design);
  - the v2 checksums are unchanged.
- No other object was renamed.
- **Swap chain:**
  - objid 195 (coasm_bin.117) is `Chitinophagaceae.2_m` and objid 150 (CAN_5_bin.149) is `Chitinophagaceae.1_m`;
  - Genome 744 ends in `_Chitinophagaceae.2_m` and Genome 785 in `_Chitinophagaceae.1_m`;
  - the v1 of each carries the same name;
  - the log shows order 11 before 12 and order 62 before 63.
- **Name resolution:**
  - 582 of 582 current names resolve to the expected objid;
  - 50 retired Assembly names and 50 retired Genome names no longer resolve. The one reused old name is
    `Chitinophagaceae.1_m`.
  - The prefix trap behaves: `Simplicispira.1` → 169 (CAN_5_bin.40) and `Simplicispira.1_m` → 308 (coasm_bin.249).

**AssemblySet 585 v3:**
- It has 291 items covering all 291 Assemblies.
- Each ref is `262319/<objid>/2`, and each label equals the R232 ID and the live object name.
- The labels are unique and in Python sort order, from `Absconditicoccaceae.1_m` to `alphaI_cluster.1`.
- Compared with v2, 51 labels changed, and v2's items are still pinned to v1.
- The name is `emilykin_MAGs_291`. The metadata has `naming = iterativeID (GTDB R232, 2026-11)` and
  `n_assemblies = 291`, provenance input is `262319/585/2`, and the description names
  `mag_iterativeID_r220_to_r232.json`.
- `sets.tsv` (index_new, ref_new, label_new) equals the live set for 291 of 291 items.

**Genomes and GenomeSets untouched:**
- The 8 GenomeSets are still v1 with unchanged checksums.
- The elements of 891 are exactly the 291 Genome v2 refs.
- **All 291 Genome v2 still carry the ontology events [SSO RAST 1.9.1, KO r13-kofam-2].**

**Report 262319/892/1:**
- `objects_created` holds exactly the 291 Assembly v2 refs plus 585/3, with no others.
- It has 13 links (12 files and 1 HTML). Every handle resolves through `hids_to_handles` to its node, and every node
  GET returns 200 with a non-zero size.
- The md5 of each blob equals its local source: 12 files, plus `index.html`, which is `output.zip`.
- The uploaded `emilykin_MAGs_gtdb_r232.tsv` has 291 rows, with ID, lineage and v2 ref matching the truth for every row.
  `..._changes.tsv` has 51 rows, exactly the changed bins, with correct old and new IDs.
- `list_referencing_objects` on the 291 Assembly v2: each is referenced only by 585/3 and 892/1, and nothing outside
  262319 refers to them.

**Other workspaces:**
- All 49 ID-named objects in `no_action.tsv` (31 in 261948, 18 in 233212) keep their name, version and type, and all
  were last saved before 2026-11-11.
- Workspace modification dates: 261948 2026-11-04, 233212 2025-11-11, 275568 2026-11-05, 221042 2025-10-25. None falls
  inside the stage's write window.

## Repo checks (`apply_scratch/verify_kbase/results_repo.tsv`, log `verify_repo_run1.log`)

**`kbase_mag_rename_log.json`:**
- HEAD and the working copy are both exactly `json.dumps(indent=1, sort_keys=True)`, LF, with no trailing newline.
- For 291 of 291 bins:
  - the keys are exactly `name_r220`, `new_name`, `upa` and `upa_r220`;
  - `name_r220` and `upa_r220` equal HEAD;
  - `new_name` equals the sim R232 ID and the live name;
  - `upa` equals the live latest Assembly ref (v2, same objid).

**`kbase_r232/` run records.** The folder holds exactly 4 files.
- The rename log has 102 rows, orders 1–102. Objid, bin, kind, old name, new name and version match live KBase and the
  inventory, and the objids are exactly the 102 renamed objects.
- The writeback log has 291 rows. `old_ref` equals the HEAD upa, `new_ref` the live v2, the lineage the prep R232 and the
  MD5 the recomputed value.
- `sets_saved.tsv` reads 585/2 → 585/3, 291 items, 51 labels changed.
- `report_saved.json` matches the live report.

**r13 `rast/genomes/_manifest.tsv`:**
- For 291 of 291 bins, `genome_ref` and `assembly_ref` are byte-equal to HEAD, and both map to the bin through live KBase.
- `iterativeID` equals the sim ID, and `genome_name` equals the live Genome name and `<genomeset>_<R232 ID>`.
- No other column changed, and the column order is the same.

**r13 tables:**
- `taxonomy.tsv`: `kbase_assembly_ref` equals the live v2 ref, and `kbase_rast_genome_ref` equals HEAD and is the bin's
  Genome. `iterativeID` and `gtdb_classification` equal the truth (291 of 291).
- `mag_panel.tsv`: `kbase_ref` equals the live v2 ref and `iterativeID` the sim ID (291 of 291).

**D5 KBase records:**
- The tracked records are unchanged against HEAD: `kbase_mag_upload_log.json`/`.out`, `kbase_kofam_push/`,
  `plan_r13/kbase_rast/`, `rast/check_fetch.out`, `fetch_rast_full.log` and `rast/genomes/<bin>.tsv`.
- Under `kbase_kofam_push/`, `plan_r13/kbase_rast/` and `rast/`, no file has been written since 13:40 CST apart from
  `_manifest.tsv` and `check_fetch.py`. That includes ignored `payloads/` and `__pycache__/`.

**Hygiene:**
- The KBase token appears in none of 819 files: the modified and untracked repo files, `apply_scratch/` and `apply_log/`.
- The 2 new scripts, `rename_mags_kbase.py` and `check_fetch.py` all parse.
- No report staging, zip or bytecode appears in `git status`.
- No tracked non-JSON file cites a MAG object as `262319/<name>`.

**Readers of the rename log's `upa`:**
- `fetch_rast.py` keys on objid only, so it is unaffected.
- `check_fetch.py` now compares ws/objid. The fix is correct; the HEAD version would flag all 291.
- `build_taxonomy.py` and `build_panel.py` copy `upa` into `kbase_assembly_ref` / `kbase_ref`. No code reads those
  columns (`git grep`).
- `validate_panel.py` uses only the keys, and `rename_mags_kbase.py` only writes the log.

## Cross-cutting end checks

- `PANEL_FROZEN.md`: all 16 hash rows that `compute_answers.py` reads (the backticked rows; for duplicate keys the last
  row wins, i.e. freeze 2a) match the files on disk, including after the docs stage's later r13 note edits.
- `AMENDMENTS.sha256`, run from `community_modeling/`: 14 of 14 OK. The track2_guards baseline: 7 of 7 OK.
- The ASV-keyed md5 baseline: 1,443 of 1,443 files unchanged.

## Fixes made

None. I found no mechanical error in KBase or in the repo's KBase-dependent files.

## Notes for the orchestrator (judgment, not defects)

1. **The report is visible in the narrative data panel.** Report 892 is visible, while all 299 older (app) reports are
   hidden. The salterns script it was adapted from does not hide its report either, so this is consistent. Hiding it is
   cosmetic and needs a KBase write, so it is Andrew's call.
2. **Step (d) was deferred (documented).** Every Genome's `assembly_ref` pins Assembly v1, which has no lineage, and
   Genome `data.id` and feature IDs keep the old strings. Join on bin or objid. A consumer that wants the lineage must use
   Assembly v2 or set 585/3.
3. **`upa` moved to v2 (documented choice).** The plan left open whether to update `upa` or add an `upa_r232` field. The
   stage updated it, which changed all 291 rows of `taxonomy.tsv` `kbase_assembly_ref` and `mag_panel.tsv` `kbase_ref`
   (freeze 2a records this). Nothing breaks.
4. **Small wrong fact in `kbase/PLAN.md` §6.** It says the manifest `genome_ref` is v2. It is v1 for all 291 bins (the
   RAST-only Genome version), both at HEAD and now. The KOfam ontology lives in Genome v2 (GenomeSet 891). Nothing was
   changed, so this has no effect on the apply. It matters only to anyone who resolves the manifest ref expecting KOfam
   terms.
5. **`kbase/build_kbase_plan.py` now fails its pre-rename assertion**, as the stage reported. This is expected; do not
   re-run it.

## Concurrent writer seen

`paper3/R232_change_list_for_Emily.md` was modified at 15:14:25 CST, during my run. It is the only repo file newer
than my start (15:04 CST). I did not touch it. It is probably the docs verifier (`apply_log/verify_docs.md` was
written at 15:17). `git status` shows 127 lines, the same as at the end of the docs stage.

## Files

- Scripts: `apply_scripts/verify_kbase/kbro.py` (read-only RPC), `verify_kbase.py` (KBase), `verify_repo.py` (repo vs
  KBase).
- Outputs: `apply_scratch/verify_kbase/`, which holds `results_kbase.tsv`, `results_repo.tsv`, `cache.json` (live
  snapshot: latest infos, set 585 items, report links and blob md5s, classified MD5s), `verify_kbase_run{1,2,3}.log` and
  `verify_repo_run1.log`. Runs 1 and 2 stopped on a blob-host timeout and on a bug in my script; run 3 is the complete
  run.
