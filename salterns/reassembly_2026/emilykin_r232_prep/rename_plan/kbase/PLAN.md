# EmilyKin GTDB R232 update: KBase change plan

**Built:** 2026-11-11, read-only, by `build_kbase_plan.py` in this folder.
- Only `get_workspace_info`, `list_objects`, `get_object_info3`, `get_objects2` and `list_referencing_objects` were called.
- Nothing was saved, renamed or deleted in KBase, and nothing in `~/Documents/EmilyKin` was changed.

**Inputs:**
- `../../results/iterativeid_simulation_r232.tsv`: 51 MAG iterativeIDs change.
- `../../results/gtdbtk.r232.summary.tsv`: GTDB-Tk 2.8.0 / R232, 291 of 291 classified. It is byte-identical in its first
  two columns to `~/gtdbtk_kbase/emilykin_r232/gtdbtk_out/gtdbtk.bac120.summary.tsv`.
- `~/Documents/EmilyKin/kbase_mag_rename_log.json`.

**If the final ID map differs from the simulation:** for example, if the repo part changes a tie-break, re-run
`build_kbase_plan.py` against it before any write. Every table here is keyed on the simulation's `committed_id` →
`simulated_r232_id`.

## 1. Scope: one workspace, not three

| ws | what it is | MAG-named objects | verdict |
|---|---|---|---|
| **262319** `freiburgermsu:narrative_1781059579957` "denitrification_bioreactor" | 291 MAG Assemblies, 291 RAST Genomes, 1 AssemblySet, 8 GenomeSets, 299 hidden reports, narrative. Last modified 2026-11-09T21:51:02Z; max objid 891 | 582 | **all changes happen here** |
| 275568 `narrative_1793811858679` "EmilyKin prFBA synthetic genomes (depth 100)" | 3,755 synthetic Genomes named by **ASV md5**, plus GenomeSet `emilykin_synthetic_genomes_prFBA_d100` whose element keys are those md5s | 0 | ASV context; no action |
| 278957 | **no such workspace** | – | no action |

**Why 275568 and 278957 looked relevant:** every repo "reference" to them is a numeric coincidence.
- 275568 appears as a cache file size in the `r12_snapshot_*.json` files and as `-0.275568` in
  `network/cooccurrence_edges_p_value_FDR.csv`.
- 278957 appears as `0.278957` in the data-inconsistency notes, as a gene coordinate in `rast/genomes/coasm_bin.53.tsv`,
  `calls/genes.tsv` and `crosswalk.tsv`, and as a correlation in the network CSV.
- No repo file contains `275568/` or `278957/` as an object reference. The ~630 repo files that cite `262319/…` are the
  real dependency: the KOfam push payloads and verify files, the RAST manifest, `rast_inputs.json`, the rename and upload
  logs, and the design notes.

**Other workspaces checked** (`no_action.tsv`):
- **261948 "codiffusion bioreactor"** holds 31 synthetic Genomes whose names equal EmilyKin MAG IDs.
  - Six of them equal **changed** IDs: old `Rhodobacter.5`, `Paracoccus.4`, `Bosea.3`, `Mesorhizobium.7`,
    `Rhodobacteraceae.1`, and new `Simplicispira.1`.
  - Their `asv` metadata md5s are not in EmilyKin's `iterativeIDs.json`, and their Source is "Synthetic". This is another
    project's ASV namespace that happens to use the same strings.
  - **Never run a name-based sweep across workspaces.**
- **233212 "ASV_genomes"** has 18 ASV-named Assemblies, Genomes and models (`Bosea.1`, `Hydrogenophaga.1`,
  `Azoarcus.2`). This is ASV context, probably EmilyKin's own 2025 ASV synthetic genomes: Source=Sludge, `.ASV.mdl`
  models, and the 262319 narrative's batch-import cell creates an `ASV_genomes` set. All of these IDs are unchanged, so
  there is nothing to do.
- **221042 "Sludge reactor microbiome"** has none.
- **`list_referencing_objects`** over all 873 Assembly/Genome versions in 262319 finds **no referencing object outside
  262319**: no models, FBAs or sets elsewhere. Inside 262319, the references come from the Genomes (`assembly_ref` and
  provenance), the sets and the reports.

## 2. What 262319 holds now (`inventory_262319.tsv`, 891 rows)

| objects | type | versions | name | taxonomy now |
|---|---|---|---|---|
| 291 MAG Assemblies (objid 3–577) | `KBaseGenomeAnnotations.Assembly-5.1` | all **v1** | current (R220-era) iterativeID, e.g. `Rhodobacter.5`. `data.assembly_id` is the **bin name** (291/291) | **none**: no `std_lineages`, no GTDB_* metadata; user metadata empty |
| 291 RAST Genomes (objid 586–888) | `KBaseGenomes.Genome-11.1` | all **v2**: v1 is RAST, v2 is the r13 KOfam ontology event from `push_kofam.py` | `emilykin_MAGs_RAST_<code25 \| code11_bN>_<iterativeID>` | `taxonomy` "Unconfirmed Organism", `scientific_name` "Unknown species", `taxon_ref` 12570/1/1, no `std_lineages` or `taxon_assignments` |
| AssemblySet 585 `emilykin_MAGs_291` | `KBaseSets.AssemblySet-1.2` | v2 | – | items `{ref: 262319/<objid>/1, label: <iterativeID>}`, sorted by label |
| 7 RAST GenomeSets (598, 829, 833, 848, 854, 885, 889) | `KBaseSearch.GenomeSet-2.1` | v1 | – | elements keyed **by ref**, pinned to Genome v1; no ID labels |
| GenomeSet 891 `emilykin_MAGs_RAST_KOfam` | `KBaseSearch.GenomeSet-2.1` | v1 | – | elements keyed by ref, pinned to Genome v2 |
| 299 hidden reports | `KBaseReport.Report-2.0` | v1 | upload / RAST / batch | historical; reference objects by UPA |
| narrative 1 | `KBaseNarrative.Narrative-4.0` | v5 | – | 4 cells: welcome, an old sludge-genome upload loop, a bulk-import cell and an `ASV_genomes` batch import. **No cell names a MAG object.** The MAG upload, rename, RAST and KOfam work was all done by API scripts, not by cells |

**No object in KBase carries a GTDB lineage at all, from any release.** Step (a) adds one for the first time; it does not
replace an R220 lineage.

**Identity checks done** (`assemblies_lineage.tsv`):
- For 291 of 291 MAGs, the KBase Assembly MD5 equals the KBase-style MD5 of the FASTA that GTDB-Tk R232 classified, and
  the contig count matches. KBase-style MD5 = md5 of `','.join(sorted(md5(upper(contig))))`, verified on CAN_1_bin.107.
  The classified FASTA is `~/gtdbtk_kbase/emilykin_r232/genomes/<bin>.fa` from the batchfile; its file md5 equals
  `local_catalogue_md5.txt` and `staged_md5.txt`.
- `data.assembly_id == bin` holds for 291 of 291.
- Every object name equals the rename log's `new_name`, and the rename log's objids are correct.

## 3. Change plan, in order

### Step 0. Pre-flight (read-only; repeat immediately before any write)
- Check that the token is valid.
- Re-list 262319 with hidden objects included and assert that:
  - the names, versions and checksums of the 291 Assemblies, 291 Genomes and 9 sets equal `preflight_snapshot.json`;
  - the maximum objid is 891.
- Other sessions wrote RAST and KOfam objects here until 2026-11-09, so a later writer must not be overwritten silently.
- Assert that none of the 102 target names exists, except the one name that the chain frees (step b).

### Step (a). R232 lineage onto the 291 Assemblies: v1 → v2 (`assemblies_lineage.tsv`)
- **What changes:** 130 lineages differ from R220 and 161 are identical. All 291 get the field, because none has it now.
- **Method:** the same as `salterns/reassembly_2026/local/writeback_gtdbtk_r232.py`, `apply`.
  - Get the full object (Assemblies are tiny: 3.9 MB in total).
  - Set `data.std_lineages.gtdb = {lineage, source_ver: "232.0", taxon_id: <lowest named rank, e.g. g__JAEXTR01>}`.
  - Keep user metadata and add `GTDB_classified_by = "local GTDB-Tk 2.8.0, GTDB R232 package"`.
  - Keep the same type: `KBaseGenomeAnnotations.Assembly-5.1` already defines `std_lineages`, as checked with
    `get_type_info`.
  - Save by **objid**. The Workspace extracts `GTDB_lineage`, `GTDB_source_ver` and `GTDB_taxon_id` from the field.
  - Order: 2 single-object canaries, then batches of 25.
- **Adaptations the salterns script needs** (write a copy, e.g. `writeback_gtdbtk_r232_emilykin.py`; do not edit the
  salterns one):
  - `WSID = 262319`.
  - `OUT = ~/gtdbtk_kbase/emilykin_r232/gtdbtk_out`. The summary is `OUT/gtdbtk.bac120.summary.tsv`. There is no ar53
    file and the `classify/` path differs.
  - `BATCH = ~/gtdbtk_kbase/emilykin_r232/batchfile.tsv`.
  - **Replace `mags()`.** It reads the salterns CheckM `selected_*.tsv` files. EmilyKin has no CheckM gate: all 291
    dereplicated MAGs are in. The "classified MD5" is computed from the batchfile FASTA with the KBase recipe above.
  - **Replace `infos_by_name()` with lookup by objid**, taken from `kbase_mag_rename_log.json` `upa`. Then assert
    `data.assembly_id == bin` before saving. Names change in step (b), so name lookup would make (a) depend on whether (b)
    ran. With objids, (a) and (b) can run in either order.
  - Provenance: the same structure, with the description "291 dereplicated EmilyKin MAGs (best-of-both catalogue,
    git 6541627)". Put `classification_method`, `closest_genome_*` and `red_value` in `custom`, and keep the R220
    lineage there too as `previous_gtdb_r220`.
  - `carries()`, the version-jump check, the reconciliation after a failed save and the TSV log are reused unchanged.
    Write the log to `rename_plan/kbase/…` or `EmilyKin/…`, not to the salterns `writeback/` folder.
- **Expected:** 291 Assemblies at v2, MD5 / Size / N Contigs unchanged, GTDB_* metadata present. Then run `verify`.
- **Not needed:** the KBase `kb_gtdbtk` app. Running it would recompute the classification on KBase's own (older)
  GTDB-Tk and data.

### Step (b). Rename 51 Assemblies and their 51 Genomes (`renames.tsv`, run in the `order` column)
- **Mechanics:** `Workspace.rename_object({obj: {wsid: 262319, objid}, new_name})`.
  - It renames every version of the object, keeps its objid and versions, and **creates no new version**.
  - Every UPA reference stays valid, both in KBase (sets, Genome `assembly_ref`, reports) and in the repo (`rast_inputs.json`,
    `_manifest.tsv`, `push_log.tsv`, `genomeset.json`).
  - Only references *by name* change meaning.
- **Order:** orders 1–51 are the Assemblies and 52–102 the Genomes, in the same bin order.
- **The one swap chain:**
  - `coasm_bin.117` `Chitinophagaceae.1_m` → `Chitinophagaceae.2_m` (Assembly 195, Genome 744) must run **before**
    `CAN_5_bin.149` `Flavipsychrobacter.1_m` → `Chitinophagaceae.1_m` (Assembly 150, Genome 785).
  - Both Genomes are in batch `code11_b5`, so the Genome names chain in the same way.
  - There is no cycle: `Chitinophagaceae.2_m` is free. No temporary name is needed; the `order` column encodes the
    dependency. The builder simulated all 102 renames in order against all 891 names, including hidden objects, and
    found no collision.
- **New names are valid** Workspace names (`^[\w|.-]+$`, not all digits).
  - Four MAGs newly take an ASV-shared ID (`Simplicispira.1`, `Holophagaceae.9`, `Taibaiella.5`, `Hoeflea.1`).
  - `Simplicispira.1` and `Simplicispira.1_m` will both exist, as different MAGs. Exact names never collide, but any
    prefix search (`get_names_by_prefix`, grep) will conflate them.
- **Do NOT re-run `EmilyKin/scripts/rename_mags_kbase.py` with a new map.**
  - It decides "already renamed" by `if iid in current` on names alone, and assumes objects still carry bin names.
  - With the current state, `CAN_5_bin.149`'s target `Chitinophagaceae.1_m` is "present" because it is coasm_bin.117's
    object, so the script would be marked done and never renamed.
  - Every other MAG would be "missing", because its current name is not its bin name, and the script would exit.
  - It also re-saves the set from v1 refs it looked up by name.
  - **Fix:** add an objid-keyed mode that reads the old → new map plus objids from `kbase_mag_rename_log.json`, renames
    in this `order`, and verifies `(objid, name)` pairs with `get_object_info3` by objid. That mode should also rename
    the Genome by its objid (from `renames.tsv`, or Genome `assembly_ref` → Assembly objid).
- **Rollback:** run `rename_object` in reverse `order` with the old names. This is lossless.

**What a Genome rename does not change (accepted, by design):**
- The Genome `data.id`, `data.source_id` (metadata "Source ID") and **every feature, CDS and mRNA ID** keep the old ID
  string, for example `emilykin_MAGs_RAST_code11_b5_Chitinophagaceae.1_m.CDS.10` inside CAN_5_bin.149's renamed Genome
  after step (b).
- Rewriting them would require re-saving 3.3 GB of Genomes.
- It would also break the KOfam ontology-term keys (`push_kofam.py`, which records `ftrs_found` against these IDs), the
  repo's `kbase_kofam_push/payloads|verify/*.json` and `rast/genomes/*.tsv`, and the crosswalks.
- Treat the internal IDs as frozen provenance identifiers. The object name and the set label are the user-facing IDs.
- **Consequence:** an internal ID such as `…_Chitinophagaceae.1_m.CDS.*` names a *different* MAG
  from the object `Chitinophagaceae.1_m` after the swap. Any consumer must join on bin or objid, never on the ID string
  inside feature IDs.
- `rast/fetch_rast.py:166` checks that the Genome **object name** ends with the rename log's iterativeID. Renaming the
  Genomes keeps that check passing once the rename log is updated.

### Step (c). Re-save AssemblySet 585 once, after (a) and (b): v2 → v3 (`sets.tsv`)
- **Items:** all 291 refs are re-pointed from `262319/<objid>/1` to the latest Assembly version that carries the R232
  lineage, as the salterns `sets` step and kb_gtdbtk do. That is expected to be `/2`, but a retried batch could make it
  `/3`, so resolve it at write time with `get_object_info3` plus `carries()`, as salterns `set_moves()` does. 51 labels
  are set to the new IDs.
- **Order:** items are re-sorted by new label (Python sort, as `rename_mags_kbase.py` did; `alphaI_cluster.1` sorts last).
  `item_index_old` and `item_index_new` are both recorded.
- **Description:** the current one names `mag_iterativeID_old_to_new.json`. Rewrite it to name the R232 map that the
  repo part writes, plus "GTDB R232 (GTDB-Tk 2.8.0) lineages in std_lineages.gtdb".
- **Metadata:** keep `n_assemblies` = 291; set `naming` = `iterativeID (GTDB R232, 2026-11)`.
- **Provenance:** `input_ws_objects` = `[262319/585/2]`.
- **GenomeSets:** nothing to do now. Their element keys are refs, so a rename changes no label.
  - The 7 RAST GenomeSets are app outputs pinned to Genome v1. Leave them as provenance.
  - GenomeSet 891 (KOfam) pins Genome v2. It needs re-pointing only if step (d) creates Genome v3.

### Step (d). Genome taxonomy: optional, recommended **deferred** (`genomes_taxonomy.tsv`)
- **Fields `Genome-11.1` supports:** `taxonomy` (string, metadata "Taxonomy"); `std_lineages.gtdb` (metadata GTDB_*,
  the same struct as on Assemblies); `taxon_assignments` (`{"gtdb": <taxon_id>}`); `taxon_ref` (currently 12570/1/1,
  "Unknown"; leave it).
- **Proposal, per MAG:**
  - `taxonomy` = the R232 names without prefixes, joined by `"; "`.
  - `std_lineages.gtdb` as for the Assembly.
  - `taxon_assignments.gtdb` = taxon_id.
  - Keep `scientific_name`, `id`, `source_id` and all features.
- **Cost:**
  - 291 full get-and-save operations, 3.32 GB in total (largest 28 MB). The 1 GB `get_objects2` cap applies even with
    `included` paths, so fetch in batches of ≤ 25.
  - Then re-point GenomeSet 891 to v3, which makes it v2.
  - The repo's KOfam verify logs pin Genome v2. That is harmless, but those files become historical.
- **Why defer:** the authoritative GTDB lineage is on the Assembly, the Genomes are RAST derivatives, and nothing reads
  Genome taxonomy today.
- **Caveat if deferred:** each Genome's `assembly_ref` pins Assembly **v1**, which has no lineage. A consumer that walks
  Genome → `assembly_ref` sees no GTDB lineage unless it resolves the latest Assembly version.

### Step (e). Optional report
- A `KBaseReport.Report`, as the salterns `report` step makes: the per-MAG R232 lineage, the 51-row rename map and the
  GTDB-Tk outputs. Its `objects_created` would list the 291 Assembly v2s and set 585/3.
- Name it e.g. `GTDBTk_R232_local_report_emilykin`. No object in 262319 has this name today (`inventory_262319.tsv`).

## 4. The 51 renames (full columns in `renames.tsv`)

| # | bin | Assembly objid | old ID | new ID | case | Genome objid | Genome name prefix |
|---|---|---|---|---|---|---|---|
| 1 | CAN_1_bin.195 | 10 | Ca_Falkowbacteria.7 | Patescibacterium.1_m | asv_inherit → gtdb_genus | 727 | emilykin_MAGs_RAST_code11_b1_ |
| 2 | CAN_1_bin.2 | 11 | Ignavibacteria.5 | Tepidaquicellaceae.1_m | asv_inherit → gtdb_lowest | 709 | emilykin_MAGs_RAST_code11_b1_ |
| 3 | CAN_1_bin.48 | 12 | Paceibacteria.1_m | Minisyncoccia.1_m | gtdb_lowest → gtdb_lowest | 595 | emilykin_MAGs_RAST_code11_b1_ |
| 4 | CAN_2_bin.173 | 31 | Burkholderiaceae.1_m | Burkholderiaceae_C.1_m | gtdb_lowest → gtdb_lowest | 754 | emilykin_MAGs_RAST_code11_b6_ |
| 5 | CAN_3_bin.10 | 49 | Rhodobacter.13 | Pseudogemmobacter.1_m | asv_inherit → gtdb_genus | 663 | emilykin_MAGs_RAST_code11_b2_ |
| 6 | CAN_3_bin.208 | 66 | Azospira_A.1_m | Proximibacter.1_m | gtdb_genus → gtdb_genus | 610 | emilykin_MAGs_RAST_code11_b1_ |
| 7 | CAN_3_bin.221 | 68 | Ignavibacteria.21 | Tepidaquicellales.1_m | asv_inherit → gtdb_lowest | 636 | emilykin_MAGs_RAST_code11_b3_ |
| 8 | CAN_3_bin.27 | 70 | Burkholderiaceae.2_m | Burkholderiaceae_C.2_m | gtdb_lowest → gtdb_lowest | 874 | emilykin_MAGs_RAST_code11_b6_ |
| 9 | CAN_4_bin.197 | 96 | Ignavibacteria.13 | Tepidaquicellaceae.2_m | asv_inherit → gtdb_lowest | 850 | emilykin_MAGs_RAST_code11_b4_ |
| 10 | CAN_4_bin.208 | 104 | Azospira.3 | Proximibacter.2_m | asv_inherit → gtdb_genus | 632 | emilykin_MAGs_RAST_code11_b2_ |
| **11** | **coasm_bin.117** | 195 | **Chitinophagaceae.1_m** | **Chitinophagaceae.2_m** | gtdb_lowest → gtdb_lowest | 744 | emilykin_MAGs_RAST_code11_b5_ |
| **12** | **CAN_5_bin.149** | 150 | Flavipsychrobacter.1_m | **Chitinophagaceae.1_m** (after #11) | gtdb_genus → gtdb_lowest | 785 | emilykin_MAGs_RAST_code11_b5_ |
| 13 | CAN_5_bin.150 | 152 | Paceibacteria.2_m | Minisyncoccia.2_m | gtdb_lowest → gtdb_lowest | 745 | emilykin_MAGs_RAST_code11_b1_ |
| 14 | CAN_5_bin.167 | 153 | Ignavibacteria.31 | Tepidaquicellales.2_m | asv_inherit → gtdb_lowest | 611 | emilykin_MAGs_RAST_code11_b2_ |
| 15 | CAN_5_bin.40 | 169 | Giesbergeria.1_m | Simplicispira.1 | gtdb_genus → asv_inherit | 816 | emilykin_MAGs_RAST_code11_b2_ |
| 16 | coasm_bin.105 | 188 | Rhodobacter.5 | Pseudogemmobacter.2_m | asv_inherit → gtdb_genus | 725 | emilykin_MAGs_RAST_code11_b4_ |
| 17 | coasm_bin.114 | 193 | Sphingobacteriales.16 | Aurantibacillus.1_m | asv_inherit → gtdb_genus | 742 | emilykin_MAGs_RAST_code11_b3_ |
| 18 | coasm_bin.121 | 206 | Sphingobacteriales.13 | Aurantibacillus.2_m | asv_inherit → gtdb_genus | 720 | emilykin_MAGs_RAST_code11_b1_ |
| 19 | coasm_bin.134 | 212 | Lysobacter_B.1_m | Solilutibacter.1_m | gtdb_genus → gtdb_genus | 857 | emilykin_MAGs_RAST_code11_b6_ |
| 20 | coasm_bin.135 | 214 | Flavipsychrobacter.2_m | Nemorincola.1_m | gtdb_genus → gtdb_genus | 799 | emilykin_MAGs_RAST_code11_b2_ |
| 21 | coasm_bin.147 | 227 | Holophaga.1_m | Holophagaceae.1_m | gtdb_genus → gtdb_lowest | 688 | emilykin_MAGs_RAST_code11_b6_ |
| 22 | coasm_bin.157 | 233 | Gammaproteobacteria.1_m | Steroidobacterales.1_m | gtdb_lowest → gtdb_lowest | 667 | emilykin_MAGs_RAST_code11_b5_ |
| 23 | coasm_bin.213 | 269 | Flavipsychrobacter.3_m | Polluticoccus.1_m | gtdb_genus → gtdb_genus | 641 | emilykin_MAGs_RAST_code11_b5_ |
| 24 | coasm_bin.215 | 272 | Chiayiivirga.1_m | Denitratimonas.1_m | gtdb_genus → gtdb_genus | 867 | emilykin_MAGs_RAST_code11_b3_ |
| 25 | coasm_bin.233 | 292 | Gammaproteobacteria.7 | Steroidobacterales.2_m | asv_inherit → gtdb_lowest | 623 | emilykin_MAGs_RAST_code11_b2_ |
| 26 | coasm_bin.237 | 293 | Holophaga.2_m | Holophagaceae.9 | gtdb_genus → asv_inherit | 847 | emilykin_MAGs_RAST_code11_b3_ |
| 27 | coasm_bin.249 | 308 | Giesbergeria.2_m | Simplicispira.1_m | gtdb_genus → gtdb_genus | 838 | emilykin_MAGs_RAST_code11_b5_ |
| 28 | coasm_bin.25 | 307 | Paceibacteria.3_m | Phycocordibacterales.1_m | gtdb_lowest → gtdb_lowest | 630 | emilykin_MAGs_RAST_code11_b4_ |
| 29 | coasm_bin.268 | 329 | Flavipsychrobacter.4_m | Polluticoccus.2_m | gtdb_genus → gtdb_genus | 866 | emilykin_MAGs_RAST_code11_b6_ |
| 30 | coasm_bin.279 | 336 | Xinfangfangia.1_m | Pseudogemmobacter.3_m | gtdb_genus → gtdb_genus | 888 | emilykin_MAGs_RAST_code11_b6_ |
| 31 | coasm_bin.286 | 347 | Sphingobacteriales.9 | Aurantibacillus.3_m | asv_inherit → gtdb_genus | 606 | emilykin_MAGs_RAST_code11_b4_ |
| 32 | coasm_bin.305 | 353 | Hydrogenedensaceae.8 | Hydrogenedentales.1_m | asv_inherit → gtdb_lowest | 835 | emilykin_MAGs_RAST_code11_b3_ |
| 33 | coasm_bin.31 | 355 | Rhodobacteraceae.1 | Pseudogemmobacter.4_m | asv_inherit → gtdb_genus | 626 | emilykin_MAGs_RAST_code11_b4_ |
| 34 | coasm_bin.312 | 356 | Flavobacteriales.1_m | Flavobacteriales_F.1_m | gtdb_lowest → gtdb_lowest | 707 | emilykin_MAGs_RAST_code11_b5_ |
| 35 | coasm_bin.314 | 364 | Spirosomataceae.1_m | Ravibacter.1_m | gtdb_lowest → gtdb_genus | 605 | emilykin_MAGs_RAST_code11_b1_ |
| 36 | coasm_bin.322 | 371 | Flavipsychrobacter.5_m | Nemorincola.2_m | gtdb_genus → gtdb_genus | 862 | emilykin_MAGs_RAST_code11_b6_ |
| 37 | coasm_bin.337 | 384 | Ignavibacteria.7 | Tepidaquicellaceae.3_m | asv_inherit → gtdb_lowest | 661 | emilykin_MAGs_RAST_code11_b5_ |
| 38 | coasm_bin.338 | 387 | Paracoccus.4 | Pseudogemmobacter.5_m | asv_inherit → gtdb_genus | 864 | emilykin_MAGs_RAST_code11_b6_ |
| 39 | coasm_bin.348 | 391 | Myxococcaceae.11 | Xiheimicrobium.1_m | asv_inherit → gtdb_genus | 749 | emilykin_MAGs_RAST_code11_b4_ |
| 40 | coasm_bin.365 | 401 | Ca_Competibacter.32 | Contendibacter.1_m | asv_inherit → gtdb_genus | 681 | emilykin_MAGs_RAST_code11_b2_ |
| 41 | coasm_bin.367 | 406 | Luteimonas_D.1_m | Luteimonas.1_m | gtdb_genus → gtdb_genus | 822 | emilykin_MAGs_RAST_code11_b4_ |
| 42 | coasm_bin.395 | 424 | Chiayiivirga.2_m | Denitratimonas.2_m | gtdb_genus → gtdb_genus | 882 | emilykin_MAGs_RAST_code11_b6_ |
| 43 | coasm_bin.399 | 426 | Bosea.3 | Allobosea.1_m | asv_inherit → gtdb_genus | 706 | emilykin_MAGs_RAST_code11_b1_ |
| 44 | coasm_bin.40 | 428 | Flavipsychrobacter.6_m | Taibaiella.5 | gtdb_genus → asv_inherit | 756 | emilykin_MAGs_RAST_code11_b2_ |
| 45 | coasm_bin.420 | 444 | Mesorhizobium.7 | Arvimicrobium.1_m | asv_inherit → gtdb_genus | 708 | emilykin_MAGs_RAST_code11_b4_ |
| 46 | coasm_bin.466 | 460 | Ignavibacteria.25 | Tepidaquicellales.3_m | asv_inherit → gtdb_lowest | 671 | emilykin_MAGs_RAST_code11_b3_ |
| 47 | coasm_bin.471 | 470 | Lysobacter_B.2_m | Solilutibacter.2_m | gtdb_genus → gtdb_genus | 825 | emilykin_MAGs_RAST_code11_b3_ |
| 48 | coasm_bin.479 | 472 | Zambryskibacteraceae.1_m | Phycocordibacterales.2_m | gtdb_lowest → gtdb_lowest | 796 | emilykin_MAGs_RAST_code11_b4_ |
| 49 | coasm_bin.532 | 513 | Paceibacteria.4_m | Phycocordibacterales.3_m | gtdb_lowest → gtdb_lowest | 665 | emilykin_MAGs_RAST_code11_b5_ |
| 50 | coasm_bin.540 | 525 | Pseudohoeflea.1_m | Hoeflea.1 | gtdb_genus → asv_inherit | 672 | emilykin_MAGs_RAST_code11_b4_ |
| 51 | coasm_bin.558 | 544 | Chiayiivirga.3_m | Denitratimonas.3_m | gtdb_genus → gtdb_genus | 755 | emilykin_MAGs_RAST_code11_b4_ |

**Genome names:** each is the prefix + new ID, e.g. objid 785 becomes
`emilykin_MAGs_RAST_code11_b5_Chitinophagaceae.1_m` after objid 744 has become `…_Chitinophagaceae.2_m`.

**ASV context:**
- 21 old IDs were `asv_inherit`, i.e. identical to an EmilyKin ASV ID, such as `Rhodobacter.5`, `Paracoccus.4`,
  `Bosea.3` and `Mesorhizobium.7`. Those ASVs keep their IDs.
- No ASV object lives in 262319, so in KBase every one of these strings is MAG context, and only these objids are
  renamed.
- The same strings in 261948 and 233212 are other ASV namespaces and must stay untouched (§1).

## 5. Risks

1. **Swap chain.** If #12 runs before #11, the Workspace rejects it with a name collision. A non-objid script could
   instead "skip as done", which is the `rename_mags_kbase.py` bug. Enforce `order` and verify by objid.
2. **Name-based references change meaning.**
   - After (b), `262319/Chitinophagaceae.1_m` resolves to CAN_5_bin.149, not coasm_bin.117.
   - A name such as `262319/Rhodobacter.5` stops resolving at all.
   - The repo's KBase scripts reference MAG objects by UPA or objid only: `fetch_rast.py` (by `upa` objid),
     `push_kofam.py`, `make_genomeset.py`, `rast_submit.py`. Sets are looked up by set name, and set names do not change.
   - No narrative cell names a MAG. The KBase data panel will simply show the new names.
   - Any note, slide or message that cites a KBase object *name* is stale. Cite bin names or UPAs.
3. **Concurrent writers.** Other sessions have written RAST and KOfam objects into 262319, most recently 2026-11-09.
   Step 0 must re-assert `preflight_snapshot.json`. The salterns-style per-batch version-jump check stops on drift.
4. **Internal Genome IDs keep the old IDs** (step b). Joins must use bin or objid. Do not "fix" them with a re-save.
5. **Version pinning.**
   - After (a), Assembly v1 still exists and is still what the Genomes' `assembly_ref`, `rast_inputs.json` and
     `_manifest.tsv assembly_ref` point to. Nothing breaks.
   - Lineage-aware consumers must use the latest version, or set 585 v3.
   - Workspace versions cannot be deleted. To roll back (a), save a v3 that copies v1.
6. **Map finality.** The plan uses the simulation's 51-row map. If the repo part's final map differs, regenerate this
   plan. A partial apply under a different map would leave KBase and the repo inconsistent.
7. **Prefix and substring hazards.** `Simplicispira.1` and `Simplicispira.1_m`, `Holophagaceae.1_m` and
   `Holophagaceae.9`, `Pseudogemmobacter.1_m`–`.5_m`. KBase names are exact, but greps and prefix searches over names
   or over the repo need word boundaries and the `_m` suffix.
8. **KBase search and indexing** may show old names or metadata for a while after (a) and (b). This is cosmetic.

## 6. Repo-side consequences (for the repo part; not done here)

- **`kbase_mag_rename_log.json`:**
  - Set `new_name` for the 51.
  - Update `upa` to `/2` for all 291 after (a), or add a field `upa_r232`.
  - `fetch_rast.py` and `check_fetch.py` read only the objid, so either works.
  - Keep the R220-era map as history.
- **`scripts/rename_mags_kbase.py`:** add the objid-keyed, ordered, Genome-aware mode (§3b), or write a new
  `rename_mags_kbase_r232.py`. Do not run the current script on a new map.
- **`community_modeling/analysis/r13/annotations/rast/genomes/_manifest.tsv`:**
  - Columns `iterativeID` and `genome_name` for the 51 rows.
  - `genome_ref` and `assembly_ref` are unchanged: v2 and v1, and a rename keeps versions.
- **Historical artefacts that should stay as they are**, recording what was pushed:
  - `kbase_kofam_push/payloads/*.json` and `verify/*.json`, and `push_log.tsv`. Their `genome_name` and feature IDs
    hold the old IDs.
  - `kbase_mag_upload_log.json`, `kbase_mag_upload.out`.
  - The `r12_snapshot_*.json` files.
- **Docs that mention `emilykin_MAGs_291` or the RAST sets by name** (`plan_r13/maps/mags.md`, `DESIGN_R13*.md`,
  `ANNOTATIONS.md`, `panel_NOTES.md`): set names do not change. Only prose that quotes a MAG object name needs the new
  ID, and only where it means the MAG.

## 7. Files here

| file | rows | content |
|---|---|---|
| `build_kbase_plan.py` | – | builds everything below (read-only KBase calls, whitelist enforced). Re-run: `~/Documents/py_venv/bin/python -I build_kbase_plan.py` (~3 min) |
| `inventory_262319.tsv` | 891 | every object, including hidden: context, planned action, new name, checksum |
| `assemblies_lineage.tsv` | 291 | step (a): current (empty) GTDB, R220 and R232 lineage, taxon_id, MD5 identity check |
| `renames.tsv` | 102 | step (b): `order`, kind, bin, old/new ID, case, chain note, `name_check` |
| `sets.tsv` | 299 | step (c): 291 AssemblySet items (ref v1 → latest R232 version, expected v2; label old→new; index old/new) and 8 GenomeSets (no action) |
| `genomes_taxonomy.tsv` | 291 | step (d), optional: current and proposed `taxonomy` / `std_lineages.gtdb` / `taxon_assignments` |
| `no_action.tsv` | 55 | 275568, 278957, the ID-named ASV objects in 261948 (other project) and 233212 (EmilyKin 2025 ASV genomes) with evidence, and the `list_referencing_objects` result |
| `preflight_snapshot.json` | – | name, version and checksum of the 291 + 291 + 9 objects, and the rename order, for step 0 |

All TSVs start with `ws, objid, name, type, ver, action, new_name`.
