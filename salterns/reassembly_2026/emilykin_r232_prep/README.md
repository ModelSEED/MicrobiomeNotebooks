# EmilyKin MAGs: re-classify with GTDB R232? (desk assessment + ready-to-run plan)

Prepared 2026-11-11 on frankensteinsmonster while poplar was unreachable. Nothing here writes to KBase, poplar or the
EmilyKin repo. The GTDB-Tk run itself has **not** been started (plan in section 6).

## 1. Bottom line

- **Consistency.** The salterns MAGs (310 original + 1,324 re-assembly) and all 12,611 wetland MAGs are on GTDB-Tk 2.8.0 /
  R232. EmilyKin is the only project still on R220.
- **The R220 run used an unsupported pairing.** GTDB-Tk v2.7.2 was run with R220 data
  (`~/Documents/EmilyKin/MANUSCRIPT_methods_results_discussion.md` line 45: every run logged "you are not using the
  reference data intended for this release: r232"). The GTDB-Tk docs support R220 only with GTDB-Tk 2.4.0–2.6.1, and
  R232 needs 2.7.0 or later (`gtdb_downloads/gtdbtk_docs_installing_index.rst`).
- **The repo summary is a hybrid.** In `gtdbtk.bac120.summary.tsv`, the 134 best-of-both reverted MAGs carry the
  classification of their own sequence, but their ANI, AF and RED come from the diffcov genomes. One run on the
  committed sequences fixes both problems.
- **Species-level result, measured.** Section 3 explains the method: a skani search against the R232 GTDB-Tk database,
  then GTDB-Tk 2.8.0's species rule. Results for the 291 MAGs:
  - **122 get a GTDB species, against 97 in R220.**
  - 26 gain a species and 1 loses one. CAN_5_bin.70 loses its species because its R220 representative,
    GCA_017302355.1, left GTDB in R226.
  - Of the 96 that keep a species, 80 keep the same name and 16 get a new name.
- **The focal EBPR taxa change.**

  | MAG | ID | R220 | R232 |
  |---|---|---|---|
  | coasm_bin.185 | *Ca.* Accumulibacter.8 | s__Accumulibacter sp000585075 | **s__Accumulibacter regalis** |
  | coasm_bin.250 | *Ca.* Accumulibacter.2 | s__Accumulibacter phosphatis_G | **s__Accumulibacter delftensis** |
  | coasm_bin.347 | *Ca.* Accumulibacter.3 | s__Accumulibacter sp017302455 | **s__Accumulibacter phosphatis** |
  | coasm_bin.476 | Accumulibacter.2_m | g__Accumulibacter (no species) | **s__Accumulibacter sp032229525** |
  | CAN_4_bin.72 | Competibacter_A.1_m | g__Competibacter_A (no species) | **s__Competibacter_A sp037045535** |
  | coasm_bin.365 | *Ca.* Competibacter.32 | g__JARSSO01 | g__Contendibacter (name map; to be confirmed by the run) |

  - For coasm_bin.347, the R220 representative was merged into the *A. phosphatis* cluster.
  - **Unchanged:** CAN_1_bin.98 (*A.* sp014879755), CAN_4_bin.236 (*Thauera* sp017997795), and the three
    *Competibacter phosphatis* MAGs (CAN_2_bin.85, CAN_4_bin.192, CAN_5_bin.106).
- **Ranks above species.** These are a projection: the R220 taxon is followed to the R232 name of the same reference
  genomes. For the 291 MAGs:
  - 176 (60 %) have an identical full lineage.
  - 197 (68 %) have an identical lineage down to genus.
  - 232 (80 %) keep the same phylum and genus names.
  - The biggest single change is p__Patescibacteria → **p__Patescibacteriota** (21 MAGs, renamed in R226).
- **Empirical prior.** On the wetland MAGs, R214 → R232 left 974 of 2,449 lineages identical (40 %) and 1,531 with the
  same phylum and genus (63 %). That spans three release steps; EmilyKin spans two (R220 → R226 → R232), so the higher
  stability projected here is plausible.
- **Main risk: identifiers.** Feeding R232 taxonomy into the current EmilyKin naming rule (`scripts/recompute_mag_ids.py`)
  would change **48 of 291 MAG iterativeIDs**. Each of the 48 changes its root, not only its number (section 5).
  - **Recommendation:** freeze the iterativeIDs, which are labels, and update only the taxonomy columns and text.
  - Several scripts match `Phylum == "Patescibacteria"` exactly; they need the new name.
- **Cost.** About 2–2.5 h wall time, with roughly 1.7 h of it at ~130–140 GB RAM in a 170 GB memory-capped scope
  (section 6).

## 2. What changed between R220, R226 and R232 (sources in gtdb_downloads/)

| | R220 (Apr 2024) | R226 | R232 |
|---|---|---|---|
| genomes | 596,859 | 732,475 | 901,341 |
| species clusters | 113,104 | 143,614 | 199,923 |
| bacterial phyla / genera | 196 / 23,112 | 183 / 27,326 | 172 / 34,834 |
| archaeal phyla / genera | 21 / 1,847 | 21 / 2,079 | 25 / 2,669 |
| GTDB-Tk versions supported | 2.4.0–2.6.1 | 2.4.1–2.6.1 | **2.7.0+** (2.8.0 used here) |

Counts are from `tables/release_stats.tsv` (computed from the taxonomy files). Release totals and dates come from
`gtdb_downloads/r*/RELEASE_NOTES.txt` and `VERSION.txt`. The GTDB-Tk compatibility row is from the GTDB-Tk installation
docs.

**GTDB release notes:**
- **R220:**
  - ANI moved from FastANI to skani for species clustering.
  - CheckM2 values were shown but not yet used.
  - 15 taxon spellings were announced for the next release, for example Leptospirae → Leptospiria and
    Hydrogenedentiales → Hydrogenedentales.
- **R226:**
  - QC now requires completeness ≥50 %, contamination <5 % and quality ≥50 under **both CheckM1 and CheckM2**. This is why
    some R220 representatives left GTDB, including EmilyKin's GCA_017302355.1.
  - Mash is no longer a clustering prefilter.
  - 20 % of heterogeneous archaeal MSA columns were pruned.
- **R232:**
  - WitChi was applied to remove compositional bias from both MSAs before tree inference.
  - The 16S extraction method changed.
  - The R202–R226 resource paper was published: NAR doi 10.1093/nar/gkaf1040.

**GTDB-Tk 2.4 → 2.8** (`gtdb_downloads/gtdbtk_changelog.rst`, `gtdbtk_docs_announcements.rst`):
- **2.4.0:** skani replaces FastANI. The summary columns `fastani_*` became `closest_genome_*`.
- **2.5.0:** skani only; the Mash prefilter was removed.
- **2.6.0:** skani 0.3.1 and pplacer 1.1.alpha19 are pinned. The 100-genome cap in dense genera was removed, so all
  representatives of a genus are now compared.
- **2.7.0:**
  - Uses a pre-sketched skani database of all R232 representatives, searched once at the start.
  - The package shrank from 198 GB to 98 GB.
  - `--skip_ani_screen` was replaced by `--place_species`.
  - Requires R232.
- **2.8.0 (2026-10-13):**
  - `classify` and `classify_wf` now use the same species rule: take the closest reference that passes the AF threshold,
    then check only that reference against its own radius.
  - Some genomes previously assigned by the ANI screen are now placed in the tree. GTDB says this affects only fuzzy
    species in condensed genera (4 of 10,000 GlobDB genomes).
  - AF rounding was harmonised and the `other_related_references` field was restored.
- **Unchanged across these versions:**
  - A species is assigned when ANI ≥ the representative's circumscription radius (95 % by default; higher for fuzzy
    species).
  - AF ≥ `--min_af` (default 0.5), where AF = max(AF_ref, AF_query).
  - The bac120 and ar53 marker sets.
  - **Data and code are version-locked.** R232 data cannot run under GTDB-Tk < 2.7.0, and R220 is not supported by 2.7.x.
    The poplar `bin_qc` environment therefore cannot do this run as-is; the local venv (`~/Documents/gtdbtk_venv`,
    2.8.0 + `~/Documents/gtdbtk_data/release232`) can.

**New species representatives:**
- 92,333 of the 199,923 R232 representatives are genomes that were not in R220 at all.
- 57,971 were not in R226 either.
- 1,068 were in R220 but not as representatives.

**Genome by genome, all of GTDB** (`gtdb_release_diff.py` → `tables/`). The R220 → R232 fate of the 113,104 R220
species representatives:
- 106,522 are still representatives.
- 2,945 were absorbed into another cluster.
- 3,637 left GTDB. 3,381 R220 representatives were already absent in R226.
- Of the 109,467 still present, **11,286 (10.3 %) changed species name**. In most of these the genus part of the name
  changed; only 1,478 change at species rank itself.
- **First rank that changed** (R220 representatives still in R232):

  | First changed rank | Representatives |
  |---|---|
  | none | 65,478 |
  | phylum, suffix only | 16,904 (mostly Bacillota_A/B/C… → Bacillota, merged in R226) |
  | phylum, new name | 6,541 |
  | class | 3,214 |
  | order | 4,245 |
  | family | 5,175 |
  | genus | 6,432 |
  | species | 1,478 |

- **Genus status** (`taxon_map_r220_r232.tsv`): 22,828 of 24,549 R220 genera keep their name (93 %); 1,271 were
  renamed, 371 split, and 79 changed only their suffix.
- **Phylum renames that matter beyond EmilyKin:**
  - Patescibacteria → Patescibacteriota.
  - Myxococcota split, mostly to Myxococcota_A (the wetland's main phylum "move").
  - Bdellovibrionota split.
  - Bacillota_A/B/C/… merged into Bacillota.
  - Desulfobacterota_B/E/G merged into Desulfobacterota.
  - Nanoarchaeota → Nanobdellota, Deferribacterota → Chrysiogenota, Dependentiae → Babelota.
- R220 → R226 and R226 → R232 are tabulated separately (`rep_fate_summary_r220_r226.tsv`,
  `rep_fate_summary_r226_r232.tsv`).

## 3. EmilyKin projection (`project_emilykin_r232.py` → `tables/emilykin_r232_projection.tsv`)

**Inputs:**
- The 291 committed best-of-both MAG FASTAs: `~/Documents/EmilyKin/dereplicated_genomes`, git 6541627, md5s frozen in
  `local_catalogue_md5.txt`. The files are git-tracked and clean at HEAD. The names match the R220 summary 291/291.
- The current R220 summary (the hybrid).

**Evidence 1, skani preview:** `skani_preview/skani_search_r232.tsv`.
- This is the exact search GTDB-Tk 2.8.0 runs first: `skani search -t 24 -d release232/skani/database --ql … --short-header`,
  skani 0.4.0 from the GTDB-Tk venv.
- Run under nice 10 in a 60 GB scope: 2 min 44 s, 42 GB peak RSS, 1,063 hits.
- The hits were filtered with GTDB-Tk 2.8.0's rule, copied from `gtdbtk/ani_screen.py`: AF = max(AF_ref, AF_query)/100
  ≥ 0.5; take the closest AF-passing reference; assign its species if ANI ≥ its radius in `radii/gtdb_radii.tsv`.
- That is what classify_wf will do, so the species column should come out the same.

**Evidence 2, name map:** each R220 taxon is followed to the dominant R232 name of the same genomes.
- This applies to ranks above species, and to the 169 MAGs with no R232 species.
- It is approximate. pplacer and RED can place a MAG differently, and 13 MAGs sit in R220 taxa that R232 split
  (dominant share < 0.9), for example g__Flavipsychrobacter (13 R232 targets) and g__Holophaga.

**Results:**

| species status (R220 → R232) | MAGs |
|---|---|
| kept, same name | 80 |
| kept, renamed (genus transfer or placeholder → Latin; e.g. Giesbergeria → Simplicispira, Lysobacter_B → Solilutibacter, Mesorhizobium_I → Arvimicrobium) | 16 |
| **gained** (e.g. *Zoogloea* sp041917255, *Brevundimonas* sp036732845, *Dokdonella* sp964390215, OLB9 sp964501555) | 26 |
| lost (CAN_5_bin.70, reference left GTDB) | 1 |
| no species either release | 168 |

| first changed rank (projected) | MAGs |
|---|---|
| identical | 176 |
| species gained / lost / renamed | 17 / 1 / 3 |
| genus: renamed / placeholder→Latin / suffix only | 20 / 8 / 2 |
| family | 15 |
| order | 15 |
| class (Verrucomicrobiae→Verrucomicrobiia 9, Paceibacteria→Minisyncoccia 5, JAEDAM01→Gracilibacteria 5, …) | 12 |
| phylum (Patescibacteria→Patescibacteriota 21, Myxococcota→Myxococcota_A 1) | 22 |

- The class row counts MAGs whose *first* change is at class; the examples count every MAG with that rename.
- **Trend.** Species-level hits rose from 96 (R220, Table S6 skani) to 105 (R226, the unadopted `mag_skani_hits/` run)
  to 122 (R232, GTDB-Tk rule). The R220 and R226 buckets used slightly different AF rules.
- **Where the projection is weakest.** Six MAGs without a species have their closest AF-passing R232 reference in a
  different genus from the projected one. Check these first after the run.

## 4. Empirical prior from the other projects
- **Wetland** (memory `wetland-gtdbtk-local.md`; KBase GTDB-Tk 2.3.2 / R214 vs local 2.8.0 / R232):
  - 974 of 2,449 lineages identical; 1,531 with the same phylum and genus.
  - Most phylum "moves" were Myxococcota → Myxococcota_A.
- **Salterns** (`salterns/MAG_QC/comparison_vs_2020.csv`; 2020 GTDB-Tk 1.x vs R232):
  - Only 42 of 310 lineages identical; 148 changed first at phylum (renamed phyla).
  - Named species rose from 19 to 176. Far more release steps lie in between, so treat this as an upper bound on churn.
- **EmilyKin (projected):** 176/291 identical and 232/291 with the same phylum and genus, in line with fewer steps.

## 5. Downstream EmilyKin artefacts that use the GTDB taxonomy (all in `~/Documents/EmilyKin` unless noted)

**Taxonomy tables:**
- `gtdbtk.bac120.summary.tsv` (copied from poplar `processed/mag/gtdbtk_291/classify/`).
- `mag_abundance_summary.{json,tsv}` (Domain…Species columns; the authoritative MAG taxonomy).
- `mag_quality_summary.tsv`.
- `iterativeID_taxonomy.json`, `iterativeID_phylums.json`, `Genus_color_map.json`.

**Naming:**
- `scripts/recompute_mag_ids.py` → `apply_mag_ids.py` (map `mag_iterativeID_old_to_new.{csv,json}`, 408cb39).
- `scripts/assign_mag_iterative_ids.py`, `scripts/regen_iterative_ids.py`.
- On poplar: `processed/mag/abundance/assign_mag_iterativeids.py`.
- These generators are one-shot. Re-running them writes an identity old → new map.
- **Simulation** (`simulate_iterativeids_r232.py`, read-only):
  - It reproduces the committed map 291/291 from the R220 taxonomy.
  - With the R232 projection, **48 IDs change**. 21 asv_inherit MAGs lose their inherit because the GTDB name no longer
    agrees with the MiDAS ASV name; the other 27 are GTDB-derived IDs whose root is renamed.
  - Examples: Ca_Competibacter.32 → Contendibacter.1_m; Ignavibacteria.N → Tepidaquicellales/-aceae;
    Paceibacteria.N_m → Minisyncoccia/Phycocordibacterales; Flavipsychrobacter.N_m → Nemorincola.N_m;
    Giesbergeria → Simplicispira; Rhodobacter.13/.5 and Paracoccus.4 → Pseudogemmobacter.
  - Full list: `tables/iterativeid_simulation_r232.tsv`.
  - The IDs are also KBase object names in ws 262319 (`kbase_mag_rename_log.json`, AssemblySet 262319/585 and the RAST /
    KOfam GenomeSets used by the r13 remodelling), so renaming them is expensive.
  - **Recommendation:** keep the IDs and update only the taxonomy columns.

**Supplementary tables and manuscript:**
- `scripts/build_table_s6_ani.py` and `reprocess_291/build_skani_novelty_291.py` (Table S6 novelty buckets, AF 0.60
  rule).
- `scripts/build_mag_quality_table.py` (re-running it regresses best-of-both; see memory).
- `reprocess_291/compute_manuscript_stats_291.py`, `build_xlsx_291.py`.
- `reprocess_291/bestofboth/compose_static.py` (it created the hybrid).
- Manuscript text: §2 GTDB-Tk paragraph, line 45, and the phylum counts on line 123.

**Figures and trees:**
- `scripts/render_fig1_tree.py:137` and `fig1_tree/compare_trees.py:96` root on `r["Phylum"] == "Patescibacteria"`.
  With R232's "Patescibacteriota" they would silently fall back to midpoint rooting.
- Also: `scripts/select_fig1_mags.py`, `fig1_tree/score_ref_candidates.py`.
- `gene_ab_figure/10_…`, `12_…`, `13_…` (novel-bacteria module).
- `clade_classify/render_nosz_tree_fig.py`, `render_nosz_top5_tree.py`.
- `network/compositional_networks/render_spieceasi_docx.py`.

**Community modelling:**
- `community_modeling/analysis/r13/annotations/build_taxonomy.py` keys its envelope rules on phylum names: `FLAG`
  ['Patescibacteria'] and `DIDERM` {'Myxococcota', …}.
- R232's Patescibacteriota / Myxococcota_A would drop the CPR flag and the diderm call for those MAGs.
- `panel/build_panel.py`, `calls/call_genes.py`, `verify_authority/verify_calls.py`.
- The genetic-code-25 set is hard-coded by MAG name (`CODE25`), so it is unaffected.

## 6. Ready-to-run plan (do not run without the user's go-ahead)

All paths are relative to this folder. The scripts are in `run/`.

| step | what | notes |
|---|---|---|
| `run/00_poplar_inventory.sh` | read-only ssh to poplar. Saves `config.yaml`, `reprocess_291/env.sh`, the `mag/` listing, both gtdbtk.log heads (version and data warning), the **server-native** R220 summaries, poplar's MAG md5s (compared with `local_catalogue_md5.txt`), and `assign_mag_iterativeids.py` → `poplar_inventory/` | needs the user's tmux + Duo login (memory `poplar-annotation-envs`). It is optional for the run itself. Per `EmilyKin/reprocess_291/env.sh`, the R220 run used `--genome_dir /scratch1/afreiburger/emilykin/EmilyKin/dereplicated_genomes` (no batchfile), data `processed/refdata/gtdbtk/release220`, env `.snakemake_envs/8c703f8a…_` |
| `run/10_stage_genomes.sh` | copies the 291 FASTAs to `~/gtdbtk_kbase/emilykin_r232/genomes`, checks them against the frozen md5s, and writes `batchfile.tsv` | the default `SRC=local` is the committed best-of-both catalogue, the directory the KBase 262319 upload script read (the upload log has no md5s to match). `SRC=poplar` rsyncs poplar's copy; use it only if 00 shows 0 differing files. Dry-tested 2026-11-11 (291/291 md5 OK) |
| `run/20_launch_gtdbtk.sh` | `systemd-run --user --scope --unit=gtdbtk-emilykin-r232 -p MemoryMax=170G -p MemorySwapMax=0 nice -n 10 bash run/run_gtdbtk_emilykin.sh …` + `../local/mem_guard.sh` (UNIT, GUARD_DIR) | the same launch as `pipeline.py gtdbtk_stage`. It refuses to start if MemAvailable < 150 GB (`WAIT=1` re-checks every 10 min), if another `gtdbtk classify_wf` is running, or if the scope exists. classify_wf runs with `--cpus 24 --pplacer_cpus 8 --min_perc_aa 10`, GTDB-Tk 2.8.0, R232 (checked via metadata.txt) |
| `run/30_collect_and_diff.sh` | copies the summaries to `results/`, then runs `diff_r220_r232.py` (per-MAG change classes; scores the desk prediction) and `simulate_iterativeids_r232.py --r232-summary` (read-only ID impact) | also diffs against the server-native R220 summary if 00 fetched it. Writes nothing into the EmilyKin repo. Dry-tested on the projection |

**Run time, from the salterns R232 run on 2026-11-11** (1,324 genomes, 2 h 21 min at 48 CPUs):
- skani screen: ~3 min (measured for these 291).
- Prodigal and markers on ~169 genomes: ~5–10 min at 24 CPUs.
- Backbone placement: ~5–7 min.
- Class-level trees: by their projected classes, the 169 genomes needing placement fall in 7 of the 8 R232 class-level
  trees (tree 4 untouched).
  Each tree costs ~10–20 min of mostly fixed pplacer load time, ~1.7 h in total.
- No archaea, so no archaeal tree.
- **Total ≈ 2–2.5 h.**

**Memory:** the salterns run peaked at 138 GB scope anon / 132.7 GB pplacer RSS with 8 pplacer threads. Hence the
170 GB cap and the 150 GB start gate. If memory is tight, `PPLACER_CPUS=1` or GTDB-Tk's `--scratch_dir` are slower but
use less.

**After the run:**
1. Read `results/r220_vs_r232_summary.tsv` and check the 6 flagged MAGs and the 13 split-taxon MAGs.
2. Decide whether IDs stay frozen (recommended).
3. Then regenerate the taxonomy columns. Use a new script that patches only the Domain…Species columns, joined on
   the bin name; do not re-run the one-shot ID generators.
4. Rebuild Table S6/S7 text, the manuscript GTDB paragraph (GTDB-Tk 2.8.0 / R232; drops the version-mismatch caveat)
   and the Fig 1 rooting string.
5. If wanted, write R232 lineages to the KBase 262319 Assemblies the way the wetland/salterns writeback did
   (`std_lineages.gtdb`, not meta).

## 7. Open items that need poplar
- Confirm what actually ran for `gtdbtk_291`: the GTDB-Tk version in the log header and the data path.
- Confirm whether poplar's `EmilyKin/dereplicated_genomes` is the all-overwrite 291 or the best-of-both.
- Fetch the server-native (non-hybrid) R220 summary.
- Compare `processed/mag/abundance/assign_mag_iterativeids.py` with the local `scripts/recompute_mag_ids.py` rule.

## 8. Files here
- `gtdb_downloads/`: GTDB R220/R226/R232 `{ar53,bac120}_taxonomy_rNNN.tsv.gz` (md5 = GTDB `MD5SUM.txt`),
  `sp_clusters_rNNN.tsv`, release notes, METHODS, VERSION, plus the GTDB-Tk changelog and docs (raw GitHub, master,
  2026-11-11).
- `gtdb_release_diff.py` → `tables/release_stats.tsv`, `rep_fate_*` (the per-representative tables were gzipped after writing), `taxon_map_*`.
- `project_emilykin_r232.py` → `tables/emilykin_r232_projection{,_summary}.tsv`.
- `simulate_iterativeids_r232.py` → `tables/iterativeid_simulation_r232.tsv`.
- `diff_r220_r232.py` (used by run/30).
- `skani_preview/`: query list, skani log, raw hits.
- `local_catalogue_md5.txt`: md5 of the 291 committed MAG FASTAs.
- `run/`: the plan scripts (section 6).
