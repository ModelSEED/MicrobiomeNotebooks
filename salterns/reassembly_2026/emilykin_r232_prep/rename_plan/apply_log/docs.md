# Docs stage: workbooks, manuscript text, tracked-change DOCX and Emily's change list (docs/PLAN.md)

Started 2026-11-12 14:20 CST in `~/Documents/EmilyKin` (branch `iterativeid-nomenclature`, HEAD 9e087c3), on top of
the uncommitted data, KBase and figures stages (`git status` 101 lines at start). Nothing is committed or pushed by
this stage and nothing is sent to KBase. Apply scripts: `rename_plan/apply_scripts/docs/`; scratch and evidence:
`rename_plan/apply_scratch/docs/`.

## 0. Pre-flight

- Baselines recorded before any edit (`apply_scratch/docs/before/`): `git status --porcelain` (101 lines),
  `git diff --numstat`, the 47 `__pycache__` directories, and md5 + CRLF count + tracked flag of every file this stage
  might touch (`files_md5_eol.tsv`). Every .md target is LF (0 CRLF lines); `novel_bacteria_summary.tsv` is CRLF.
- R232 numbers re-derived from the current repo files, not from the plan: `compute_manuscript_stats_291.py` copied to
  scratch with `R="."` and the case column fixed to `assignment_case` (the repo script is unchanged), run read-only in
  the repo (`apply_scratch/docs/manuscript_stats_r232.txt`):
  15 phyla (Patescibacteriota 21, Myxococcota 2, Myxococcota_A 1), order 291, family 290, genus 283 (183 genera),
  species 122 (169 unnamed); methods 51 topology / 97 topology+ANI / 122 ani_screen / 21 RED; cases 138/70/83;
  Table S6 >=95 121 / 85-95 76 / <85 28 (23 no hit + 5 divergent; R220: 40 + 8) / low_AF 66; guilds unchanged
  (62/13/77/117/22); nosZ 37/88/1; ppk1 6.
- GTDB-Tk command and versions from the prep `results/gtdbtk.log`: GTDB-Tk v2.8.0, `classify_wf --batchfile ...
  --cpus 24 --pplacer_cpus 8 --min_perc_aa 10 --scratch_dir ...`, reference data r232, skani v0.4.0 ANI screen
  (122 genomes classified by the ANI pre-screen).
- Accumulibacter values (repo `gtdbtk.bac120.summary.tsv` closest_genome_ani / af, and `Table_S6_skani_novelty.tsv`
  skani ANI / AF_query): coasm_bin.185 *A. regalis* 97.58 / 0.838; coasm_bin.250 *A. delftensis* 98.60 / 0.851;
  coasm_bin.347 *A. phosphatis* 98.14 / 0.827; CAN_1_bin.98 sp014879755 99.48 / 0.963; coasm_bin.476 sp032229525
  98.68 / 0.837 (skani AF_query 72.88 %); CAN_4_bin.64 sp003538495 96.05 / 0.743.
- Small inconsistency noticed in a data-stage file (not edited here): the header comment of
  `reprocess_291/skani_novelty_r232/run_skani_r232.sh` says "727 hits for 269 of the 291 MAGs"; the hit table has
  268 distinct query MAGs (23 with no hit), as the data log says.

## A. Workbooks (bin-keyed cell edits)

Script `apply_scripts/docs/patch_docs_workbooks.py`; gate `apply_scripts/docs/gate_docs_workbooks.py`; change list
`apply_scratch/docs/workbook_changes.tsv` (2,805 cells); logs `wb_verify_head.log`, `wb_apply.log`, `wb_gate.log`,
`wb_idempotent.log`.

- Backups (byte copies, new untracked files, not gitignored): `MAG_target_gene_filled_UPDATED (1).pre_r232.bak.xlsx`,
  `cladeexcel.pre_r232.bak.xlsx`, `abundance and coverage (1).pre_r232.bak.xlsx`.
- Rule: rows keyed by the bin column; ID -> `mag_iterativeID_old_to_new.json` (R232; overwrites the 6a04588-era drift
  IDs, D4); rank / lineage -> repo `gtdbtk.bac120.summary.tsv` (R232); a blank R232 rank is written as an empty cell
  (as the sheets store blanks). The builder workbook `MAG_target_gene_filled.xlsx` (data stage B6) was not touched.
- **S6_ANI method.** The submission sheet = builder sheet columns A..T + the 7 `Table_S6_skani_novelty.tsv` columns
  (U..AA). Verified at HEAD (`--verify-head`, HEAD builder copy + HEAD TSV): the copy rule changes 0 of the S6 cells.
  So S6 B..T are copied by bin from the rebuilt builder sheet (no `--skani`, method parity, data stage B8a) and U..AA
  from the rebuilt TSV (D3 skani search), with numbers typed as the sheet stores them.

| file / sheet | ID | Phylum | Class | Order | Family | Genus | Species | lineage | other |
|---|---|---|---|---|---|---|---|---|---|
| UPDATED (1) S1_MAG quality | 82 | 22 | 23 | 19 | 24 | 60 | – | 130 (`Taxonomy `) | |
| UPDATED (1) S2a_RPKM per MAG | 54 | | | | | | | | |
| UPDATED (1) S3_gene copies | 82 | | | | | | | 130 | |
| UPDATED (1) S4_nosZ_genecopy_type | 42 | | | | | 22 (`genus`) | | | |
| UPDATED (1) S5 _ppk_ACC | 2 | | | | | | | | |
| UPDATED (1) S6_ANI | 82 | 22 | 23 | 19 | 24 | 60 | 43 | 130 | reference name 43, accession 36, ANI 71, AF 74, call 27, interpretation 40, Notes 137, Reference database 291; skani block: best_ref 93, best_taxonomy 142, ANI 93, AF_query 92, AF_ref 93, ANI_class 54, Novelty_interpretation 90 |
| UPDATED (1) S7_functional | 82 | 22 | | | 24 | 60 | | | |
| cladeexcel.xlsx S4 | 42 | | | | | 22 | | | |
| abundance and coverage (1).xlsx Sheet1 | 82 | | | | | | | | |

ID 82 = 51 R232 changes + 31 stale-era cells (D4); S2a 54 = 31 + 23; S4 42 = 27 + 15.

**Gate (PASSED):**
- all 1,429 non-skani `table_cells_plan.tsv` cells of these three files are changed to the planned value (0 missing,
  0 wrong) and no non-S6 cell outside the plan changed; 414 of the 1,779 S6 recompute cells are not listed in the plan
  (blank cells that R232 fills, e.g. a closest reference for a newly resolved MAG);
- the full cell diff against each backup equals the change list (2,659 / 64 / 82 cells); no style change; S1 keeps
  its 2,910 formulas (still without cached values, critic C17, unchanged); zip part list and defined names unchanged;
- D1: every ID column (7 sheets + cladeexcel S4 + abundance Sheet1) equals the R232 map for every bin, no duplicates;
- idempotent: a dry re-run finds 0 cells to change (so S6 now equals builder B..T + TSV U..AA for 291/291).
- A second apply from the backups gave identical sheet XML; only `docProps/core.xml` (openpyxl save timestamp)
  differs between runs.
- Housekeeping: the first dry-run overwrote the change list with an empty file (script bug, fixed: dry runs no longer
  write it); the workbooks were restored from the backups and re-patched to regenerate it.

## B. `MANUSCRIPT_methods_results_discussion.md` (direct edits)

Script `apply_scripts/docs/edit_manuscript_md.py` (22 exact, unique sentence-level replacements, each asserted to
occur once; no global replace); log `apply_scratch/docs/md_edit.log`. LF kept (0 CRLF before and after);
`git diff --numstat` 15/15 lines (long paragraphs, one line each).

| line | change |
|---|---|
| 45 | GTDB-Tk v2.7.2 / r220 -> **v2.8.0 classify_wf against GTDB R232** (command from the R232 `gtdbtk.log`: `--batchfile --cpus 24 --pplacer_cpus 8 --min_perc_aa 10 --scratch_dir`; release232 data; skani v0.4.0 ANI pre-screen; KBase report 262319/892/1; the map file `mag_iterativeID_r220_to_r232.csv`); method counts 145/97/49 -> **122 ani_screen / 97 topology+ANI / 51 topology / 21 RED**; the version-mismatch caveat **deleted** |
| 61 | 155/62/74 -> **138/70/83** on R232; "51 MAG identifiers changed (50 at the root), no ASV identifier changed" |
| 87, 153 (x2) | CAN_5_bin.40 *Giesbergeria* -> *Simplicispira* (JJ008 at line 153 is CAN_5_bin.112 and stays) |
| 103 | skani v0.3.2 / r220 database -> **skani v0.4.0 / GTDB-Tk release232 database (199,923 reps)**, `reprocess_291/skani_novelty_r232/`, R232 taxonomy; added the `-n 3` caveat (coasm_bin.535 is low_AF_unresolved although GTDB-Tk gives *Novosphingobium* sp017302615, 98.4 %, AF 0.90) |
| 107 | CAN_4_bin.208 "*Azospira_A*" -> "*Proximibacter* (r220 *Azospira_A*)"; reference choice now reads "r220 *Azospira_A* representative ... *Azospira_A restricta* (GCF_016858125.1; *Proximibacter restrictus* in GTDB R232)"; "61 GTDB r220 reference genomes (labelled with their GTDB R232 names; GCA_017302355.1 ... keeps its r220 name)". The r220 tree inputs (mask, reps MSA, v2.7.2 identify/align re-implementation) are kept: the tree was relabelled, not rebuilt (D3) |
| 123 | r220 -> R232; **15 phyla** (14 if Myxococcota_A is merged); Patescibacteriota (21); Myxococcota (2) + Myxococcota_A (1); family 290; genus **283 (183 genera)**; species **122 (169 unnamed)**; RED **21**; S6 **121 / 76 / 28 (23 no hit) / 66**; "only about a third ... conspecific" -> "about two-fifths (42 %)" |
| 141 | 155/62/74 -> 138/70/83 (r220 values in parentheses) |
| 147 | coasm_bin.185 *Accumulibacter* sp000585075 -> ***Accumulibacter regalis*** |
| 157 | added (D7): GTDB R232 species agree with the ppk1 clades (regalis 97.58 %, delftensis 98.60 %, phosphatis 98.14 %) and coasm_bin.476 = *Accumulibacter* sp032229525 (98.68 %, AF 0.84; species added in R226), not species-novel |
| 181 | Discussion caveat "Finally, the GTDB taxonomy is r220 although the GTDB-Tk binary (v2.7.2) natively expects r232 ..." **deleted** |
| 191 | Fig 1 band list replaced by the regenerated figure's 20 bands (figures log `C5_render_percell_composition.log`): ... *Dokdonella*, *Lysobacter* (GTDB *Lysobacter_I*), JJ008, SpSt-398, *Desulfobacillus*, *Zoogloea*, *Proximibacter*, *Simplicispira*, Unclassified (Unclassified moved from 8th to 20th; SCMT01/Azospira/Giesbergeria renamed) |
| 225 | Table S5 legend: coasm_bin.476 genome in GTDB R232 species *Accumulibacter* sp032229525 (98.68 %) |
| 229 | Fig S4 legend: references "labelled with their GTDB R232 names"; rooted on the **Patescibacteriota** clade; "no skani hit to any GTDB R232 representative"; references labelled with their GTDB R232 species |
| 231 | Table S6 legend: GTDB R232 / skani v0.4.0; **121 / 76 / 28 (23 no hit) / 66** (r220 96/70/48/77 in parentheses) |

Remaining r220 / v2.7.2 / Azospira_A / Patescibacteria strings in the file are deliberate (history and comparison
values: lines 45, 61, 107, 123, 141, 229, 231). Line 77's `*Ca.* Contendobacter` is the MiDAS genus roster of
`gao_pao_categories.json` (ASV level), not a MAG genus: unchanged (D6 only touched the Table S7 regex).

## C. Tracked changes in the two working DOCX files

Helper `apply_scripts/docs/docx_track.py` (own surgical OOXML code; the docx skill's `comment.py` for comments and
its `validate.py` / `soffice.py` / `accept_changes.py` for checks). Author **"Claude (GTDB R232 update)"**, initials CR.
- A match is located in the paragraph's visible text; only the runs at its two ends are split (rPr cloned); matched
  runs go into `w:del` (`w:t` -> `w:delText`), the replacement into a `w:ins` after them. No paragraph is rebuilt, so
  the Oct-19 tracked changes of "Claude (A. Freiburger)" are preserved. A run inside such an existing `w:ins` (all of
  Table 1) is isolated by splitting that `w:ins` (author/date kept) and the `w:del` is nested inside it.
- Images: the drawing run goes into `w:del`, a clone pointing at a new media part (`word/media/image_r232_<k>.png`,
  new relationship, unique docPr id, cy rescaled to the new aspect) into `w:ins`.
- Driver `apply_scripts/docs/run_docx.sh`: byte-copy backup `<name>.pre_r232.bak.docx` (new untracked file), unzip the
  backup into scratch, edit, rezip, then `validate.py --original <backup>` (XSD, new errors only) and
  `validate.py --original <backup> --author "Claude (GTDB R232 update)"` (redlining). Only a validated scratch file
  was copied into the repo. Logs: `apply_scratch/docs/{main,si,spieceasi}_{edit,validate,redline}.log`; change
  records `apply_scratch/docs/{main,si,spieceasi}_docx_changes.json`.
- **Pre-existing defect fixed in both `_correction` files:** `word/settings.xml` carried `<w:trackChanges/>` (written
  by the Oct-19 helper; not an OOXML settings element; the baseline XSD check of the untouched files fails on it). It
  was replaced by a schema-valid `<w:trackRevisions/>` in CT_Settings order. No content change; it also turns Track
  Changes on when Emily opens the file. (SPIEC-EASI's pre-existing `w:zoom` defect, from python-docx, was left.)

### C1. `papert3_070626 (2)_correction.docx` (52 tracked changes, 6 comments; XSD and redlining PASSED; 897 -> 897 paragraphs)

| where | change | comment |
|---|---|---|
| para 45 | "v2.7.2" -> "v2.8.0 (classify_wf)"; "release r220" -> "release R232" | R232 rerun; optional R232 citation (plan's Parks et al. 2025 NAR doi, flagged "please verify") |
| para 57 | "14 bacterial phyla" -> 15; "Patescibacteria (21)" -> "Patescibacteriota (21)"; "Myxococcota (3)" -> "Myxococcota (2), Myxococcota_A (1)" | 15 vs 14 is the author's convention choice (coasm_bin.258) |
| para 59 (Table 1 legend) | "GTDB-Tk (r220) genus" -> "(R232)" | lists all 16 Table 1 cell changes and notes the old IDs remain ASV names |
| Table 1 | 8 iterativeID + 8 GTDB-genus cells, located by the bin in column 1 (`table1_plan.tsv`, each current value asserted); Latin genus names italic, alphanumeric placeholders roman (the table's convention). **D1 gate on all 53 rows: ID == R232 map and genus (or f__family) == R232, 0 mismatches** | – |
| Figure 1 (image1) | replaced (tracked) by the regenerated `figures/percell_composition.png` (4543x2271; the embedded image was a 1543x770 export of the same figure, same aspect) | legend changes, JJ008 = CAN_5_bin.112 only, Unclassified 1.89 -> 0.57 % |
| para 656 | "48" -> "28", "77" -> "66"; **D7**: "three of these were closely related ... likely represents a novel Accumulibacter lineage." -> "all six belong to GTDB R232 species, and five are closely related to their species representative (>=97% ANI), including the near-complete coasm_bin.476 (...; *Accumulibacter* sp032229525, 98.68% ANI)." | R232 evidence: S6 rebuild numbers, coasm_bin.476 GTDB-Tk 98.68 % / AF 0.84 to GCF_032229525.1 (added in R226), skani AF_query 72.9 %, all six ANIs, ppk1 statement kept |
| Figure 4 caption (para 779) | no text change | Emily's 2026-10-21 Figure 4 (PDF p. 24) labels the x-axis **by bin** and needs nothing; the older ID-labelled image here would need 17 relabels, listed bin-keyed (`apply_scratch/docs/fig4_label_map.tsv`: all 47 visible labels resolved uniquely through the workbook's 6a04588 map and each bin's 5-class matches its panel; 7 R232 changes + 10 older-map drifts) |
| para 657 | unchanged: 4/5/41/118 verified (Table S7 classifications unchanged) | – |

Images 2/3/5 (no IDs) and the legacy module panels image6/image7 (ASV node labels; not in Emily's master; D2) were
not touched.

### C2. `SI-paper3 (3)_correction.docx` (22 tracked changes, 4 comments; XSD and redlining PASSED; 478 -> 478 paragraphs)

| where | change | comment |
|---|---|---|
| para 209 (ppk1) | **D7** sentence inserted after "... coasm_bin.476 to clade ID (100/100).": R232 species *A. regalis* / *A. delftensis* / *A. phosphatis* match the clades; coasm_bin.476 = *Accumulibacter* sp032229525 (98.68 %) | evidence + flags (not changed) the pre-existing "candidate-novel" wording for CAN_1_bin.98, a known GTDB species (99.48 %) in both releases |
| para 212 (SI-6) | *Giesbergeria* -> *Simplicispira* (CAN_5_bin.40); ", labelled with their GTDB R232 names," after "closest GTDB r220 reference genomes" | – |
| para 216 (SI-7) | " classified with GTDB R232" inserted; 155/62/74 -> **138/70/83** | 51 IDs changed (50 root), no ASV ID changed, map file |
| para 252 (Fig S4 legend) | "(labelled with their GTDB R232 names; GCA_017302355.1 is not in R232)"; Patescibacteria -> **Patescibacteriota** clade; CAN_4_bin.208 "(Proximibacter in GTDB R232)", "r220 Azospira_A representative" x2, "; Proximibacter restrictus in R232" after GCF_016858125.1; skani "GTDB r220 representative" -> "R232". "61 GTDB r220 ..." kept (relabel, D3) | – |
| Figure S4 (image5, md5 = HEAD `fig1_tree.png`) | replaced by the regenerated `fig1_tree/fig1_tree.png` | 8 tip relabels, R232 reference names, GCA_017302355.1 note, 10 red-circle class changes |
| Figure S5 (image4, md5 = HEAD `nosz_top5_tree.png`) | replaced by the regenerated `clade_classify/out/nosz_top5_tree.png` (5311x2416) | only CAN_3_bin.221 -> Tepidaquicellales.1_m |

### C3. Checks on both files

- Rendered with LibreOffice (`soffice --convert-to pdf`, 32 and 22 pages) and inspected: Table 1 pages, the Fig 1
  page (old image struck through, new image below), para 656, SI pages with Fig S4/S5 and the legends.
- Accepted view (`accept_changes.py`, then text of every `w:t`): no "v2.7.2", "release r220", "14 bacterial phyla",
  "Patescibacteria (", "Patescibacteria clade", "GTDB-Tk (r220)", Giesbergeria, SCMT01, JARSSO01, UBA919,
  "likely represents a novel", "155 ASV-inherited", Ignavibacteria.21, Ca_Falkowbacteria.7, Paceibacteria.1_m remain;
  the 3 remaining "Azospira_A" in the SI are the deliberate r220 reference-selection wording.
- The repo files are byte-identical to the validated scratch outputs.

## D. `SPIEC-EASI_network_analysis.docx` (12 tracked changes, 2 comments; XSD and redlining PASSED)

Node labels stay ASV IDs (D2). python-docx paragraph numbers below; the plan's 140/357 are its all-paragraph indices.

| para | change |
|---|---|
| 53 (plan 140) | "gene-content PAO-like MAG (Azospira.3," -> "... PAO-like MAG-linked ASV node (Azospira.3 [MAG CAN_4_bin.208, Proximibacter.2_m under GTDB R232],"; "(.26, .32, .51)" -> "(.26, .32 [MAG coasm_bin.365, Contendibacter.1_m under GTDB R232], .51)"; comment (Pseudoxanthomonas.23, Zoogloea.2, .26, .51 keep their IDs, checked by bin) |
| 72 (plan 357) | "gene-content PAO-like MAGs (Azospira.3," -> "... MAG-linked ASV nodes (Azospira.3 [MAG CAN_4_bin.208, Proximibacter.2_m under GTDB R232]," |
| 14 | the link rule now reads "equalled a MAG iterativeID under the GTDB r220-based naming ...; these 140 ASV–MAG pairs are frozen ..." and the source file `network/mag_linked_asvs_r220_frozen.tsv` replaces `mag_iterativeID_old_to_new.json` (this sentence stated the rule D2 replaced; comment explains D2) |
| 16 | guild join: "the current naming map (`mag_iterativeID_old_to_new.json`), not through the guild table's own iterativeID column, which carries an older naming map." -> "the frozen ASV–MAG list (`network/mag_linked_asvs_r220_frozen.tsv`), never through an iterativeID string." (the root `spieceasi_modules_flanking.py` joins by the frozen list since the data stage) |

Generator kept consistent: `network/compositional_networks/render_spieceasi_docx.py` got the same wording (new
`mag_note()` reads the frozen list + current map, so the bracket notes are data-driven; numstat 27/9). It was NOT run
into the repo (that would overwrite the tracked document). Check: run with its output path argument into scratch
(`apply_scratch/docs/spieceasi_regen_check.docx`); its text equals the accepted view of the tracked document for
370 / 370 paragraphs, and no repo file was written by that run.

## E. Living notes and the novel-bacteria module (D5)

Script `apply_scripts/docs/edit_living_notes.py` (exact unique substrings; bin-keyed cell patches via the repo's
`scripts/_binkeyed.patch_table`); log `apply_scratch/docs/living_notes.log`. Line endings kept (all LF except the
CRLF `novel_bacteria_summary.tsv`).

| file | change | numstat |
|---|---|---|
| `fig1_tree/README.md` | CAN_4_bin.208 = Proximibacter.2_m (Azospira.3 before); rooting sentence notes Patescibacteriota and the subtitle rule; reference labels = R232 names since 2026-11-12; new section 7 "GTDB R232 relabel" (tip/reference relabels, red circles, root name, hazard note on the alignment scripts). r220 names in sections 2–6 kept as build history | 24/4 |
| `clade_classify/NOSZ_TOP5_NOTES.md` | PNG 5311x2416 since the relabel (was 5309); CAN_3_bin.221 = Tepidaquicellales.1_m (3 places) | 6/5 |
| `clade_classify/CLADE_I_II_DISCREPANCY.md` | CAN_5_bin.40 *Simplicispira* (r220 *Giesbergeria*), 3 places | 3/3 |
| `gene_ab_figure/METHODS_nosz_clades.md` | same, 2 places | 2/2 |
| `clade_classify/out/summary.md` (output of `classify_nosz.py`, not re-run) | bin-keyed genus cell CAN_5_bin.40 -> Simplicispira (the only genus change among its bin rows) | 1/1 |
| `phase_rpkm_chord/README.md` | "Giesbergeria.1" -> Simplicispira.1 = CAN_5_bin.40, noting the un-re-rendered chord figures still say Giesbergeria_m.1 | 2/1 |
| r13 `crosswalk_NOTES.md` | dated R232 note (phylum rename, 51 renames by objid with versions/checksums unchanged, Assemblies v2 with R232 lineage, feature IDs keep old names: join on bin, B9 tables, PANEL_FROZEN freeze 2a); Patescibacteria -> Patescibacteriota (3); Genome-name sentence notes the renames | 6/4 |
| r13 `kofam_NOTES.md` | Patescibacteriota (2; coasm_bin.479 order UBA9973 -> Phycocordibacterales) | 2/2 |
| r13 `panel_NOTES.md` | "JARSSO01" -> "Contendibacter coasm_bin.365 (r220 JARSSO01)" | 1/1 |
| `gene_ab_figure/data/novel_bacteria_summary.tsv` | bin-keyed RED / method / novelty for the 7 rows from the R232 summary (coasm_bin.260: "not species-novel in GTDB R232 (s__JAAEKA01 sp049294015; ANI 98.33%)", ani_screen, RED N/A; the other six "species-novel (s__; no GTDB R232 species; topology and ANI)" with R232 RED) | 7/7 |
| `gene_ab_figure/data/novel_bacteria_summary.md` | Table N1 ID + RED by bin, N2/N3 IDs via N1's ID->bin (Dokdonella.14 -> .8, D4); R232 note (6 of 7 remain species-novel; CAN_5_bin.70 lost its species; figure 10 now R232 set vs figure 13 / N2-N3 r220 set: author decision) | 11/9 |
| `novel_bacteria_supplementary.md` (276-era narrative) | dated R232 note under the title only; body kept as written | 2/0 |

**Not edited, on purpose:**
- `community_modeling/analysis/r13/annotations/ANNOTATIONS.md` and `plan_r13/track2/GENOME_OBJECT_SPEC.md`: both are
  in the ASV md5 baseline (`asv_md5_before.tsv`) and contain no MAG iterativeID. The only MAG-meaning string is "5
  Patescibacteria" (the r13 gene-calling run) in ANNOTATIONS.md. The R232 note was put into crosswalk_NOTES.md.
- `DESIGN_ANSWERS.md`, `plan_r13/TRACK2_DECISIONS_2026-11-09.md`: their IDs (Ignavibacteria.31, Ignavibacteria.5,
  ...) are r13 member (ASV) names (critic C24); no MAG object name is cited. `FINALIZE_REVIEW_LOG.md`, `validate_full.txt`,
  `check_fetch.out`, `fig1_raxml_run.log`, `SERVER_HANDOFF_*`, `REVIEW_papert3_071626.md`, `manuscript_stats_*.txt`:
  logs / historical records (D5).
- `bvbrc_embedding_hits/REPORT_emilykin_hits.md`, `synthetic_genomes_d100/abstained_asvs.md`: "Patescibacteria" there
  is a BV-BRC reference / MiDAS ASV name. `dbRDA/...scores.md`: Bosea is an ASV root.
- `MANUSCRIPT_methods_results_discussion.md` line 77 "*Ca.* Contendobacter": MiDAS genus roster (ASV level).
- The do-not-edit list of docs/PLAN.md section 3 (re-downloads, superseded drafts, `MAG_target_gene_filled (1).xlsx`,
  `MAG_target_gene_filled_UPDATED.xlsx`, `MAG_target_gene_empty_tables.xlsx`, `abundance and coverage.xlsx`, etc.)
  and every `*.pre_*.bak.*`: untouched. The builder workbook `MAG_target_gene_filled.xlsx` (data stage) untouched.
- Brief deck / program overview: ASV context; the optional slide-13 "five Patescibacteria" note was not changed
  (generated deck; would need `figure_data.json` / `make_brief_pptx.py`).

## F. Change list for Emily

`~/Documents/EmilyKin/paper3/R232_change_list_for_Emily.md` (new, untracked, not ignored), written by
`apply_scripts/docs/build_change_list.py`. Section 1 is keyed by page / manuscript line of `paper3/papert3_071626.pdf`
(14 items: p7 l229–230, p8–9 l307–310, p9 l321–322, p10 l345, Table 1, p13 Fig 1, p13 l359–360, p15 l389–391,
l408–412, l416–418, l422, p28 l711–712, abstract l17–18, p24 Fig 4 = no change). The 28 Table 1 lineage rows are
generated by bin from HEAD (r220) and repo (R232) summaries; the 28 equal `table1_plan.tsv` and are exactly the rows
of the 54-MAG Table 1 set whose lineage changes. Section 2 covers the SI by section/paragraph and the workbook;
sections 3–5 list the scientific (D7) changes, name traps and open choices. Master facts checked: the closed genomes
still span five phyla (13/6/4/1/1); master Figure 4 labels by bin; master Figure 1 is the same 1543x770 export as the
working DOCX image.

## G. Final checks (`apply_scratch/docs/`)

- ASV-keyed md5 (`apply_log/asv_md5_before.tsv`): **1,443 / 1,443 unchanged** (`asv_md5_end.txt`).
- `community_modeling/plan_r13/integrator/AMENDMENTS.sha256` (run from `community_modeling/`): **14 / 14 OK**;
  r13 track2_guards (`track2_guards_sha256_before.txt`): **7 / 7 OK**.
- D5: no log, run record, `*.pre_*.bak.*` (existing), payload/verify JSON, upload log, r12 snapshot or track2_guards
  file touched. New backups only: the six `*.pre_r232.bak.{xlsx,docx}` (byte copies, mtime preserved).
- Line endings: no change in any edited text file (`eol_numstat_end.txt`); the CRLF table keeps CRLF (7/7 numstat).
- `git status`: 101 -> 127 lines, all attributable to this stage (`git_status_end.txt`): 21 modified tracked files
  (3 xlsx, 3 docx, `MANUSCRIPT_methods_results_discussion.md`, `render_spieceasi_docx.py`, 12 notes / module files,
  plus `novel_bacteria_summary.tsv`, already modified by the data stage) and 7 new untracked files (6 backups +
  `paper3/R232_change_list_for_Emily.md`). `find -newer <stage start>` lists only these files.
- No `__pycache__` created (47 directories before and after; all scripts run with `-B` / `PYTHONDONTWRITEBYTECODE`).
- Nothing committed, nothing pushed, nothing sent to KBase. No one-shot ID generator or poplar-bound builder was run;
  `render_spieceasi_docx.py` was run only into scratch.

## H. Hand-offs / open items

- **Author decisions (also in Emily's list):** 15 vs 14 phyla; whether to cite the R232 release (the plan's Parks et
  al. 2025 NAR doi was not verified here); Fig S4 relabel (done) vs rebuild; novel-bacteria module set (r220 seven vs
  R232 set; figure 10 vs 13 disagree); CAN_1_bin.98 "candidate-novel" ppk1 wording in the SI (pre-existing, flagged).
- **Not regenerated, still showing older labels (figures-stage skips):** the nosZ reference tree, the phase_rpkm_chord
  figures, `render_denitrification_heatmaps` outputs, the module panels. None is embedded in Emily's master; the
  working DOCX's older Figure 4 got a comment instead.
- **Builder S6 `--skani` question (data stage):** the submission S6 A..T follows the builder (no `--skani`, 169
  unresolved); U..AA carry the D3 skani values.
- Emily's newest SI is not on disk; SI items are keyed by section/paragraph of `SI-paper3 (3)_correction.docx`.
- `run_skani_r232.sh` header comment says 269 MAGs with hits; the table has 268 (data-stage file; not edited).

Finished 2026-11-12 ~15:25 CST.
