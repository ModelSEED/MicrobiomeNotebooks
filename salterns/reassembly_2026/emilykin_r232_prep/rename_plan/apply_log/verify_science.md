# Independent verification, science lens (GTDB R232 update of ~/Documents/EmilyKin)

Run 2026-11-12, 15:00-15:40 CST, on the uncommitted working tree left by the gate, data, KBase, figures and docs stages
(HEAD 9e087c3, branch `iterativeid-nomenclature`, 127 `git status` lines at start and at end). The stage reports were not
trusted. Every check below reads the repo files, the prep `results/gtdbtk.r232.summary.tsv`, the GTDB R232 downloads,
the skani database or KBase directly.

- My scripts: `rename_plan/apply_scripts/verify_science/`.
- Their outputs: `rename_plan/apply_scratch/verify_science/`.
- Everything ran with `~/Documents/py_venv/bin/python -I`. skani ran under `nice 10` with 24 threads or fewer, `free -g`
  was checked first (about 173 GB available), and peak memory was 42 GB.
- KBase was only read (`get_object_info3`, `get_workspace_info`). Nothing was committed or pushed.

## 1. Verdict

The science checks pass. The one exception is a statement about the skani `-n 3` method (section 3.1). It is a
judgment item and is reported rather than fixed.

| check | result | evidence |
|---|---|---|
| Lineages by bin in every taxonomy-bearing file | 0 mismatches across 30 files and 80 file/column pairs, the 2 JSON files and `mag_genera.json` | `check_taxonomy.out` |
| Workbook taxonomy and ID cells | 0 mismatches in the builder and submission workbooks (S1, S1b, S2a, S3, S4, S5, S6, S7), cladeexcel S4 and abundance Sheet1 | `check_workbooks.py` (console) |
| KBase ws 262319 | All 291 Assemblies are at v2, carry the R232 name, and their `GTDB_lineage` metadata equals the R232 classification. All 102 renamed objects carry the new name, and the Genomes keep their version. Set 585 is at v3. Max objid is 892 | `kbase_readonly.out` |
| Table S6, re-derived from the skani output | Re-derived buckets 121 / 76 / 28 / 66 equal the table, with 0 cell mismatches out of 291 × 6 columns. All best-hit taxonomies are in R232 | `check_table_s6.out` |
| Table S6 skani run reproduced | A full re-run with the same settings gives 727 of 727 rows identical, and the best hit and bucket are identical for 291 of 291 | `skani_full/`, section 2.2 |
| Table S7 / guild CSV | The only changed columns are `iterative_ID` (82), Genus (60), Family (24) and Phylum (22). No classification column changed | `check_table_s7.py` |
| D2 network freeze | The frozen list equals the HEAD map ∩ the 935 abundance columns, 140 of 140. All network tables equal HEAD row for row once the added fields are dropped, and those fields equal the R232 map by bin | `check_network.py` |
| Manuscript numbers | Every recomputed taxonomy number equals the text (section 2.4) | `manuscript_stats.out` |
| D7 claim edits | Accurate (section 2.5) | |
| Figures against HEAD renders | Every change has an R232 cause or is a D4-type stale-label fix, apart from one figure (section 3.2) | `figs/` |
| D1 sweep | 44 file/column pairs plus the JSON: the ID equals the R232 map by bin, and no ID is shared by two bins | `d1_sweep.out` |
| Line endings | No line-ending change in any of the 109 modified tracked files | `d1_sweep.out` |
| ASV-keyed md5 | 1,443 of 1,443 unchanged | |
| AMENDMENTS.sha256 | 14 of 14 OK, run from `community_modeling/` | |

## 2. Checks in detail

### 2.1 Taxonomy by bin

`check_taxonomy.py` works per file and column:
- it learns the HEAD convention against R220 (raw GTDB name versus suffix-stripped);
- it then requires each working-tree cell to equal R232 under that same convention, for every row whose HEAD cell
  matched R220.

Results:
- **No mismatches.** Every file in `data/taxonomy_plan.tsv` marked patch or regenerate has 0 mismatches.
- **Files covered.** These include `gtdbtk.bac120.summary.tsv`, whose content equals the R232 summary and stays CRLF,
  `mag_quality_summary`, `mag_abundance_summary` (.tsv and .json), mimag, S6, S7, guild, percell, `gene_ab_figure/*`,
  the nosZ tables, `ko_rpkm_*`, `phase_rpkm_chord/*`, the r13 taxonomy/panel/calls tables, `fig1_annotation` and the
  leftover-reads JSONs.
- **The 291-preview table.** Its 13 bins with a blank genus are the never-filled preview rows, which were already blank
  at HEAD.
- **r13 `taxonomy.tsv`.** All of its GTDB-Tk fields equal the R232 summary for 291 of 291: method, note, closest
  genome and placement reference, ANI, AF, taxonomy, RED, msa and warnings.
- **Fig S4 references.** 60 of 60 labels, genera, species and phyla equal `bac120_taxonomy_r232`.
  GB_GCA_017302355.1 is absent from R232 and is labelled "(not in GTDB R232)", as D3 requires.

### 2.2 Table S6

`check_table_s6.py` applies the Table S6 rule as I wrote it independently of the builder:
- the best hit is the highest ANI among hits with AF_query ≥ 15 %;
- if there is none, the highest ANI overall gives low_AF_unresolved;
- with no hit, the MAG is `<85`.

Results:
- **Equality with the table.** The re-derivation equals the committed table in `ANI_class`, `skani_best_ref`, ANI,
  AF_query, AF_ref, `skani_best_taxonomy` and `mag_gtdb_taxonomy` for 291 of 291, and the sentence agrees with the
  class for 291 of 291.
- **Reference taxonomy.** The `bac120_taxonomy_r232.tsv.gz` copy and the GTDB-Tk `gtdb_taxonomy.tsv` copy agree for
  every best reference.
- **Agreement with GTDB-Tk.** Bucket ≥95 and a GTDB-Tk R232 species agree for 290 of 291. The exception is
  coasm_bin.535, which is low_AF in S6 but *Novosphingobium* sp017302615 by GTDB-Tk.
- **Bucket transitions.** From HEAD to R232, 95 MAGs stay ≥95 and 26 move into ≥95. The full transition table is in
  `check_table_s6.out`.
- **Reproducibility.** A full re-run of the stage's exact command (skani 0.4.0, the release232 database, `--min-af 0
  -n 3`, `-t 24`) took 2 min 54 s at 41.8 GB and gave 727 rows identical to
  `reprocess_291/skani_novelty_r232/mag_vs_gtdb_r232_skani.tsv`.
  - An earlier 5-MAG spot run (`-t 16`) differed only in one lower-ranked tiny-AF hit, the third hit of
    coasm_bin.476. Its best hits and buckets were identical.
  - The tail of the `-n 3` list therefore depends on the query batch. The full-batch run that S6 uses is
    reproducible.
- **Workbooks.** In the submission workbook, the S6_ANI skani block (U..AA) equals the TSV for 291 of 291. The
  builder workbook has no skani block. Columns B..T are GTDB-Tk only: there are 169 "Unresolved by ANI" rows, and
  "Reference database" reads "GTDB R232 (GTDB-Tk v2.8.0, classify_wf)", which shows the C7 fix is in.

### 2.3 Networks (D2)

- **The frozen list.** Taking the HEAD `mag_iterativeID_old_to_new.json` values present among the 935 `abundances.csv`
  columns reproduces `network/mag_linked_asvs_r220_frozen.tsv` exactly: 140 of 140 (ASV, bin) pairs. Under the R232
  map the link would be 127 (17 drop and 4 join: Hoeflea.1, Holophagaceae.9, Simplicispira.1, Taibaiella.5), which is
  why the freeze matters.
- **SPIEC-EASI and paper tables.** For the 4 edge lists and 2 node or membership files in `out/asv/` and the 4
  `tablesforpaper` CSVs, the rows equal HEAD row for row once the added columns are dropped. Row counts: 717 / 266 /
  401 / 159 edges, 53 MAG nodes, 89 module nodes, 871 / 339 flanking rows and 14 PAO-GAO edges. Every added
  `*_MAG_iterativeID` cell, 4,555 in all, equals the R232 map by bin.
- **Louvain membership JSON.** `network/louvain_module_mag_membership.json` has 641 member records with the same
  (node, bin) set as HEAD, and every `mag_iterativeID` is correct. The other differences from HEAD are text only:
  `_description`, `name_map` and 6 `gamma_rule` strings. The gamma values, module sizes and members are unchanged; the
  gamma_rule wording was refreshed to match the network summaries.
- **Untouched outputs.** No other file under `network/` changed (all-phases n_nodes 140, n_edges 5,978), and
  `module_pathway_completeness/` did not change.

### 2.4 Manuscript numbers (`MANUSCRIPT_methods_results_discussion.md`, working DOCX, SI DOCX, change list)

These were recomputed from the R232 summary (`manuscript_stats.py`) and from the repo tables:

- **Methods (R232).** ani_screen 122, topology+ANI 97, topology 51, RED 21. R220 was 97 / 145 / 49.
- **Phyla.** There are 15:
  - Pseudomonadota 131, Bacteroidota 62, Chloroflexota 36, Patescibacteriota 21;
  - Verrucomicrobiota 10, Planctomycetota 10, Acidobacteriota 8, Spirochaetota 4;
  - Myxococcota 2, Bdellovibrionota 2;
  - 1 each for Myxococcota_A, Gemmatimonadota, Hydrogenedentota, Armatimonadota and Elusimicrobiota.
- **Rank resolution.** Order 291, family 290, genus 283 (183 genera) and species 122 (169 unnamed). R220 was
  genus 273 / 172 genera and species 97. Species: 26 gained and 1 lost (CAN_5_bin.70).
- **Cases.** 138 / 70 / 83 (R220 155 / 62 / 74), with 51 IDs changed (50 at the root). The map equals the simulation,
  and `mag_iterativeID_old_to_new.{csv,json}` and `_new_ids.json` `mag_new` equal the map for 291 of 291.
- **Table S6.** 121 / 76 / 28 / 66. Of the 28 `<85`, 23 have no hit. ≥95 is 121/291 = 41.6 %, as "two-fifths (42 %)"
  says.
- **Accumulibacter.** There are six:
  - 185 *A. regalis* 97.58 / 0.838;
  - 250 *A. delftensis* 98.60;
  - 347 *A. phosphatis* 98.14;
  - 476 sp032229525 98.68 / AF 0.837;
  - CAN_1_bin.98 sp014879755 99.48;
  - CAN_4_bin.64 sp003538495 96.05.

  GCF_032229525.1 is absent from `bac120_taxonomy_r220` and present as a species representative in R226 and R232.
  That confirms "added in R226".
- **Competibacter.** *Competibacter* + *Competibacter_A* make 8 MAGs and *Accumulibacter* 6. The percell genus sums are
  unchanged: Competibacter 29.54 / 30.32 / 30.14 / 14.63 / 10.49, Accumulibacter 10.89 / 6.63 / 7.24 / 31.12 / 36.51,
  Thauera peaks at 11.30 in CAN_2 and Azonexus at 6.64 in CAN_4.
- **Fig 1 bands.** The recomputed top-20 order equals the manuscript list, ending "... *Proximibacter*,
  *Simplicispira*, Unclassified". Unclassified goes from 1.89 to 0.57, and JJ008 from 1.08 to 1.04.
- **Other species named in the text.** CAN_2_bin.85, CAN_4_bin.192 and CAN_5_bin.106 are *C. phosphatis*,
  CAN_4_bin.236 is *Thauera* sp017997795, CAN_3_bin.203 is JAEULV01, CAN_3_bin.221 is OLB5, and the nosZ II.D genera
  are as stated.
- **Table 1, PDF.** `papert3_071626.pdf` pp. 10-12 gives 54 bins. Of these, 28 have a lineage change, which equals the
  change list's §1.1 set exactly. All 28 "what changes" entries and R232 lineages are correct.
- **Table 1, working DOCX.** In the accepted view, 53 of 53 rows have ID and genus (or f__family) equal to R232 by bin.
- **Stale-name scan.** I scanned the .md and the accepted text of both DOCX files for the R220-only GTDB names and the
  51 old MAG IDs. The only hits are deliberate r220 history (`Azospira_A` in the r220 reference-selection wording,
  "Patescibacteria in r220") and `Bosea.1`. `Bosea.1` is the unchanged ASV-inherited ID of CAN_4_bin.220, whose genus
  JAIUNP01 is correct.

### 2.5 D7 claim edits

These were checked against the data:
- **Main DOCX para 655 (the coasm_bin.476 reversal).** "all six belong to GTDB R232 species, and five are ... ≥97 %
  ANI, including coasm_bin.476 (... sp032229525, 98.68 %)" is true by both GTDB-Tk and S6 skani. The ≥97 % set is 185,
  250, 347, 476 and CAN_1_bin.98, and CAN_4_bin.64 is 96.05. The Word comment gives the correct evidence.
- **Matching clades.** "Matching these clades" (SI para 208, md line 157) holds: IA regalis, IC delftensis, IIA
  phosphatis. At r220 the species of coasm_bin.250 was *phosphatis_G*.
- **coasm_bin.260.** It is no longer species-novel (JAAEKA01 sp049294015, 98.33 %, ani_screen). The
  `novel_bacteria_summary.tsv` RED, method and novelty columns equal the R232 summary for 7 of 7.
- **Fig S4 comment.** "10 MAGs change ANI class" holds; the 10 bins listed equal the `fig1_annotation` diff.

## 3. Findings that need a decision (not changed)

### 3.1 The skani `-n 3` artefact affects 17 MAGs, not one

Two documents state that the `-n 3` recipe artefact touches a single MAG:
- `MANUSCRIPT_methods_results_discussion.md` line 103, added by the docs stage: "... is called `low_AF_unresolved`
  even when a full-length hit exists; **this affects one MAG**, coasm_bin.535 ...";
- the change list `paper3/R232_change_list_for_Emily.md` §3 item 5: "**One** Table S6 call is an artefact".

That is true only at species level. I ran an unrestricted skani search: same database, `--min-af 0`, no `-n`, on the 66
low_AF_unresolved MAGs (`skani_lowaf_non/`, `check_n3_artifact.py`, `n3_artifact.out`). 17 of the 66 have a hit with
AF_query ≥ 15 % that `-n 3` pushed out behind three tiny high-identity fragments. They re-bucket as follows:
- **1 to ≥95:** coasm_bin.535.
- **12 to 85-95:**
  - CAN_3_bin.124, CAN_3_bin.43, CAN_4_bin.208, CAN_4_bin.25;
  - coasm_bin.2, coasm_bin.208, coasm_bin.256, coasm_bin.324;
  - coasm_bin.335, coasm_bin.558, coasm_bin.563, coasm_bin.75.
- **4 to <85:** CAN_1_bin.177, CAN_3_bin.208, coasm_bin.315, coasm_bin.395.

| | ≥95 | 85-95 | <85 | low_AF |
|---|---|---|---|---|
| Table S6 as built (`-n 3`) | 121 | 76 | 28 | 66 |
| without the truncation | 122 | 88 | 32 | 49 |

- **Why this is not an R232 error.** The recipe matches R220, which also used `-n 3` and had 77 low_AF. The bucket
  counts in the text are therefore method-consistent.
- **What is affected.** The quantifier is wrong, and so is the main-text sentence "a further 66 could not be
  confidently placed because <15 % of the genome aligned to any reference" (working DOCX para 655 and PDF p13
  l359-360). That sentence is false for 17 of the 66.
- **Options for Andrew and Emily:**
  - (a) Keep the recipe and change the md sentence to: "this changes the species-level call of one MAG,
    coasm_bin.535 (...); in all, 17 of the 66 `low_AF_unresolved` MAGs have an alignment of ≥ 15 % to a
    representative outside their three highest-ANI hits (1 ≥95 %, 12 at 85-95 %, 4 <85 %)". Add the same note to
    change-list item 5.
  - (b) Rebuild S6 without `-n 3`. That changes the method and needs S6 rebuilt on both releases for comparison.

### 3.2 A figure changed for a non-R232 reason

- **Figure 13.** `gene_ab_figure/novel_metabolic_potential.png` differs from HEAD only in one x-axis label,
  `Dokdonella.14` → `Dokdonella.8`, for CAN_1_bin.210. The bbox is 1372,3961-1579,4140 and 0.03 % of pixels differ.
  Dokdonella.14 is a pre-R232 stale older-map ID. The re-render equals the gate's HEAD re-run, so the change is a
  D4-type label fix. D4 names only the workbooks, Table S7 and the guild CSVs, so this is outside D4's literal scope.
  The new label is the correct current ID. Keep it, or restore it from git if figures must change only for R232.
- **Other figures with D4-type fixes.** Figures 04 and 06 (`Ca. Competibacter.139` → `.85`, `Dokdonella.14` → `.8`)
  and the organism tree (`Ca_Accumulibacter.6`/`.8`, `Ca_Competibacter.49`, `Ca_Leptovillus.11`,
  `Ca_Sarcinithrix.2`) carry the same kind of fix alongside their R232 changes.
- **Organism tree canvas.** The canvas goes from 3842×2495 to 3843×2479 because the committed PNG was stale (gate
  finding). The topology and support dots are identical (visual check, `figs/sbs_organism.png`).
- **Every other changed figure has an R232 cause**, by pixel bbox and visual side-by-side (`figs/`):
  - **Fig S4:** the only diffs are in labels, red circles, title and legend. Branches are identical
    (`figs/fig1_diff_overlay.png`).
  - **Main Fig 1:** the band renames and the Unclassified and JJ008 changes.
  - **Fig S5:** one label changes, and the tree is horizontally compressed about 3 % to fit the longer label
    `Tepidaquicellales.1_m`. The topology is identical.
  - **Figures 04 and 06:** changes only in the label band.
  - **Figure 10:** the novel set changes.
- **The DOCX images match.** The embedded `image_r232_*` files are md5-identical to the regenerated `fig1_tree.png`,
  `nosz_top5_tree.png` and `percell_composition.png`.

### 3.3 Other open items (confirmed, already raised by the stages)

- **Novel-bacteria module.** It is split between sets:
  - figure 10 uses the R232 set (+CAN_5_bin.70 Azonexus.12, −coasm_bin.260);
  - figure 13 and Tables N2-N3 keep the r220 seven;
  - figure 13's subtitle "no GTDB species assignment" is now false for coasm_bin.260.
  - CAN_5_bin.70 is "species-novel" in R232 only because its r220 reference GCA_017302355.1 (98.67 % ANI) left GTDB.
    Biologically it is not novel against known genomes.
- **Stale renders.** These were skipped by gate verdict and still carry R220 or older labels: the nosZ reference tree,
  the phase_rpkm_chord chords and heatmaps, `denitrification_MAG_heatmap` and the legacy module panels. I did not
  verify them against R232.
- **KBase.** The Genome `assembly_ref` still pins Assembly v1, which has no lineage. The lineage is on v2 and in set
  585/3.
- **B4 `ko_rpkm_root_*` tables.** They still use the 924fa18-era roots (tier C, not regenerated).
- **Builder S6, columns B..T.** These are GTDB-Tk only (no `--skani`): 169 rows read "No close ANI hit" while the skani
  block shows hits, for example CAN_1_bin.107 at 86.64 %. This is pre-existing and part of the open `--skani` decision.

## 4. Fix made (clear, mechanical)

| file | what | how | before → after |
|---|---|---|---|
| `reprocess_291/skani_novelty_r232/run_skani_r232.sh` (new, untracked, data stage) | The header comment said "727 hits for 269 of the 291 MAGs". The hit table has 268 distinct query MAGs, 23 with none, and both my re-derivation and the full re-run confirm it | One-line Edit: "(727 hits for 268 of the 291 MAGs; 23 have no hit)" | md5 2ce11c6a… → 27b34e5e…. LF kept (0 CR). Line count unchanged (18). No other repo file was edited |

After the fix, ASV md5 is 1,443 of 1,443, AMENDMENTS.sha256 is 14 of 14 OK, and `git status` still has 127 lines (the
file is untracked, so the count is unchanged).

## 5. Scratch outputs (`apply_scratch/verify_science/`, 19 MB)

- `check_taxonomy.out`, `check_table_s6.out`, `manuscript_stats.out`, `changed_columns.out`, `d1_sweep.out`;
- `kbase_readonly.out`, `n3_artifact.out`, `stale_names_scan.txt`;
- `main_docx_changes.txt`, `SI-paper_changes.txt`, `SPIEC-EA_changes.txt`, `table1_pdf_bins.txt`;
- `skani_spot/`, `skani_full/` (the full re-run), `skani_lowaf_non/` (the unrestricted search);
- `figs/` (HEAD and current PNGs, side-by-sides, the Fig S4 diff overlay).
