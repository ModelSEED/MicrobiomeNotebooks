# Independent verification, docs lens (verify_docs)

Run 2026-11-12, 14:58–15:20 CST, in `~/Documents/EmilyKin` (branch `iterativeid-nomenclature`, HEAD 9e087c3) on top of the
uncommitted data, KBase, figures and docs stages. I did not rely on the stage reports: every expected value below was
re-derived from primary sources. These were the committed R232 map (`mag_iterativeID_old_to_new.json`, equal to
`results/iterativeid_simulation_r232.tsv` for 291/291), the repo `gtdbtk.bac120.summary.tsv` (equal to the prep
`results/gtdbtk.r232.summary.tsv` apart from CR), HEAD's R220 summary from `git show`, the rebuilt `Table_S6_skani_novelty.tsv`,
the builder `MAG_target_gene_filled.xlsx`, and the PDF `paper3/papert3_071626.pdf`.
Scripts and evidence: `rename_plan/apply_scratch/verify_docs/`. Nothing was committed or pushed and nothing was sent to KBase.

**Result:** all checks pass apart from one mechanical error in Emily's change list, which I fixed (section 6). The items that
need a judgment call are in section 8.

## 1. Backups

- All six `*.pre_r232.bak.*` files (3 DOCX, 3 XLSX) are byte-identical to `git show HEAD:<file>`, so every
  diff below is measured against HEAD.

## 2. DOCX files (`docx_verify.py`, own code; logs `main_verify.log`, `si_verify.log`, `spiec_verify.log`)

Method:
- On the current `document.xml`, every `w:ins` by "Claude (GTDB R232 update)" was removed and every `w:del` by that
  author was unwrapped, along with the new comment markers.
- The result was compared with the backup as a stream of characters and objects. Each character carried its text, its
  run properties (rsids ignored) and the chain of other authors' `w:ins`/`w:del` wrappers (author and date). The
  stream also covered paragraph, table, row, cell and section properties, bookmarks, and drawings (`r:embed` and extent).
- "Accepted view" means the same author's changes accepted. I diffed it against the backup paragraph by paragraph.

| check | main `papert3_070626 (2)_correction` | SI `SI-paper3 (3)_correction` | `SPIEC-EASI_network_analysis` |
|---|---|---|---|
| python-docx load | OK (219 body paragraphs, 2 tables) | OK (221, 3) | OK (86, 7) |
| new-author revisions | 26 ins + 26 del = 52, one date 2026-11-12T20:54:43Z | 13 + 9 = 22, 20:47:00Z | 6 + 6 = 12, 20:48:02Z |
| other authors' revisions | "Claude (A. Freiburger)": del 43 and ins@rPr 594 unchanged; ins 695 → 697 (split around nested deletions) | "Claude (A. Freiburger)" ins 52 → 64 (split), del 38 and rPr 5 unchanged; "Andrew Freiburger" 3/3 unchanged | none before or after |
| **reject check (stream equality)** | **PASS** (67,099 items) | **PASS** (43,553) | **PASS** (31,888) |
| nesting | 0 ins-in-ins, 0 ins-in-del, 0 del-in-del; 17 del-in-(old)ins (valid: text inserted by A, deleted by B); 0 `w:t` inside `w:del`; 0 `delText` outside `w:del` | same; 6 del-in-ins | none nested |
| revision ids | new ids unique; duplicate ids 9000–9003 are **pre-existing** in the backup | new unique; duplicates 9004–9006 pre-existing | unique |
| paragraphs | 897 → 897 | 478 → 478 | 370 → 370 |
| comments | 6 new (author and initials CR correct), each with start, end and reference markers, anchored on the edited paragraph (p44, 56, 58, 653 image, 655, 778 Fig 4 caption); the old comment is unchanged | 4 new, on p208, 215, 250 (Fig S4 image), 253 (Fig S5 image) | 2 new, on p14 and p139 |
| skill `validate.py --original` (XSD) and `--author` (redlining), re-run on the **repo** files | PASSED / PASSED | PASSED / PASSED | PASSED / PASSED |
| zip parts changed | document.xml, rels, comments*.xml, settings.xml; added media/image_r232_1.png | same; added image_r232_1/2.png | document.xml, rels, Content_Types; comments*.xml added |

- The reject check covers `document.xml` only. `word/settings.xml` in both `_correction` files also changed:
  `<w:trackChanges/>` became `<w:trackRevisions/>`, in schema order. That edit cannot be tracked, and "reject all" does
  not undo it. It is benign and schema-correct (the backup's element is not a CT_Settings child), and it turns Track
  Changes on when the file is opened.
- python-docx counts one fewer `inline_shapes` in main and SI (6 vs 7, 3 vs 4). Its XPath only sees drawings in runs
  that are direct children of `w:p`, and both the deleted and the inserted image runs now sit inside `w:del`/`w:ins`.
  This is expected.

### 2.1 Planned edits against the accepted view (docs/PLAN.md §4.2/§4.3, every non-ASV-keep row of `occurrences.tsv`)

- **Main.**
  - p44: v2.7.2 → v2.8.0 (classify_wf); r220 → R232.
  - p56: 15 phyla, Patescibacteriota (21), Myxococcota (2) + Myxococcota_A (1). The counts recomputed from the summary
    sum to 291.
  - p58: Table 1 legend, (r220) → (R232).
  - Table 1: 16 cells. Independently, **53/53 rows** have ID equal to the map and genus (or f__family) equal to R232 by bin.
  - p653: Figure 1 image.
  - p655: 48 → 28, 77 → 66, and the D7 Accumulibacter sentence. Recomputed: the <85 class is 28 with 23 no-hit;
    low_AF is 66; five of six are ≥97 % (CAN_4_bin.64 is 96.05); coasm_bin.476 is sp032229525 at 98.68 %.
  - para 657 (4/5/41/118): unchanged, correctly.
  - Fig 4: comment only (see section 7).
  - No other paragraph changed.
  - Accepted-view scan: no r220, R220, v2.7.2, "Patescibacteria", Azospira, Giesbergeria, SCMT01, JARSSO01, UBA919,
    "aligned poorly", sp000585075, Flavipsychrobacter, Paceibacteria, Ignavibacteria.N, Ca_Falkowbacteria or "14
    bacterial" remains.
- **SI.**
  - p208: D7 sentence. The R232 species of 185/250/347/476 and the ANI values match the summary.
  - p211: Simplicispira, plus "labelled with their GTDB R232 names".
  - p215: 138/70/83. Matches `case_new` in `mag_iterativeID_r220_to_r232.json`.
  - p250 and p253: Fig S4 and Fig S5 images.
  - p251: Fig S4 legend.
  - The remaining "r220"/"Azospira_A" strings are the deliberate D3 tree-history wording.
- **SPIEC-EASI.**
  - p14 and p16: the D2 frozen-list rule.
  - p139 and p356: "MAG-linked ASV node(s)" with the bracketed bin and R232 ID.
  - The ASV node labels (Azospira.3, Ca. Competibacter.32, Ignavibacteria.5/.13, and SILVA Patescibacteria/Myxococcota)
    are kept, as D2 requires.

### 2.2 Comment facts (checked against data)

- **Fig 4 comment (main).** The 17-label relabel map is correct 17/17:
  - every old label equals the backup workbook's S7 `MAG iteration ID` for that bin;
  - every new label equals the R232 map;
  - every "R232" or "older map" tag is right.
- **Fig S4 comment (SI).** Checked against `fig1_annotation.tsv`, HEAD vs current:
  - exactly 10 MAGs change `ANI_class`, and the named examples are right;
  - the reference relabels are right (Contendibacter sp018262135, Lysobacter_I sp019510125, Patescibacterium
    sp002344425, JAFIGV01 sp937871155, Niabella defluvii, A. regalis, A. delftensis, Proximibacter restrictus);
  - GCA_017302355.1 has the "not in GTDB R232" label.
- **"GCF_032229525.1 added to GTDB in R226".** It is absent from `bac120_taxonomy_r220.tsv.gz` and present in r226 and
  r232 as s__Accumulibacter sp032229525.
- **Accumulibacter species and ANI/AF values in all comments** equal the R232 summary: 185 regalis 97.58/0.838; 250
  delftensis 98.60/0.851; 347 phosphatis 98.14/0.827; 98 sp014879755 99.48/0.963; 476 sp032229525 98.68/0.837;
  64 sp003538495 96.05/0.743.

### 2.3 Images

| file | old image (now inside `w:del`) | new image (inside `w:ins`) | extent / aspect |
|---|---|---|---|
| main Fig 1 | media/image1.png (1543×770). Visually the same figure as the HEAD `percell_composition.png`, in a different export | image_r232_1.png, md5 equal to the repo `figures/percell_composition.png` (a0728e7d…) | 6409715×3204152, aspect 2.0004 = pixel aspect |
| SI Fig S4 | image5.png, md5 equal to HEAD `fig1_tree.png` (inside the old author's `w:ins`) | image_r232_1.png, md5 equal to the repo `fig1_tree/fig1_tree.png` (6094e7b0…) | square, 1.0000 |
| SI Fig S5 | image4.png, md5 equal to HEAD `nosz_top5_tree.png` | image_r232_2.png, md5 equal to the repo `clade_classify/out/nosz_top5_tree.png` (ab2b1b2d…) | aspect 2.1983 = pixel aspect |

- No stale PDF was embedded. `docPr` ids are unique.
- LibreOffice renders of copies (`render/main_fig1_pages.png`, `render/si_figS4_pages.png`): the old image is struck
  through and the new image appears below it.
- The new Fig 1 legend reads Lysobacter, Proximibacter, Simplicispira, UBA5532 (Patescibacteriota), with Unclassified
  last.

## 3. Workbooks (`wb_verify.py`; log `wb_verify.log`; cell list `wb_changes_independent.tsv`)

Method:
- Both files of each pair were loaded with openpyxl (`data_only=False`), and every coordinate of every sheet was
  compared for value, type and style.
- Expected values came from the sources:
  - ID → map;
  - rank and lineage → R232 summary, with each column's format first proved against HEAD R220 for 291/291 (blank rank =
    empty cell);
  - S6 B..T → builder sheet by bin;
  - S6 U..AA → TSV by bin.
- S6 HEAD parity: the backup S6 equals the HEAD builder B..T plus the HEAD TSV U..AA, with 0 mismatching cells. This
  confirms the copy rule independently.

| file | changed cells | explained | wrong | unexplained | missing | style diffs | formulas | structure (merged, freeze, filter, CF, DV, widths, heights, names, sheets) |
|---|---|---|---|---|---|---|---|---|
| `MAG_target_gene_filled_UPDATED (1).xlsx` | **2,659** (S1 360, S2a 54, S3 212, S4 64, S5 2, S6 1,779, S7 188; S2b, KO_order_mapping and Notes_per_cell 0) | 2,659 | 0 | 0 | 0 | 0 | 2,910 → 2,910 (raw `<f>` 2,910 in both) | identical |
| `cladeexcel.xlsx` | **64** (ID 42, genus 22) | 64 | 0 | 0 | 0 | 0 | 0 | identical |
| `abundance and coverage (1).xlsx` | **82** (ID; Sheet2 0) | 82 | 0 | 0 | 0 | 0 | 0 | identical |

- D1: every ID column (S1, S2a, S3, S4, S5, S6, S7, cladeexcel S4, abundance Sheet1) equals the map for every bin row,
  with no duplicates.
- No title, notes or non-bin row changed.
- The zip part lists are identical. The changed parts are only the edited sheets plus `docProps/core.xml`.
- S6 checked directly against R232: lineage and Phylum..Species have 0 mismatches over 291 rows, and `Reference
  database` reads "GTDB R232 (GTDB-Tk v2.8.0, classify_wf)" in 291/291 rows.
- String scan of every cell of all three workbooks, for R220-only names, "r220", "v2.7.2" and the 51 old-only IDs
  (with boundaries): 0 stale hits. The only matches are the R232-valid family/order names GWB1-50-10 and UBA2400, and
  the formula `=AVERAGE(N220:R220)` (a cell reference).
- The builder workbook was not touched by the docs stage; its S6 is the data stage's.

## 4. `MANUSCRIPT_methods_results_discussion.md`

- `git diff --numstat` 15/15, LF.
- Every changed number was recomputed:
  - methods 122/97/51/21;
  - 15 phyla with the per-phylum counts;
  - genus 283 (183 genera); species 122 (169 unnamed); RED 21;
  - S6 121/76/28 (23 no hit)/66;
  - cases 138/70/83;
  - 122/291 = 42 %;
  - release232 skani reps 189,801 + 10,122 = 199,923.
- The Fig 1 band list equals the regenerated legend order.
- r220 and Patescibacteria remain only as explicit history or comparison values.

## 5. Change list for Emily (`paper3/R232_change_list_for_Emily.md`) against the PDF (40 pp.)

- All 14 section-1 items were checked with `pdftotext -layout`.
  - Each quote is present on the stated page at the stated margin lines: p7 l229–230, p8–9 l307–310, p9 l321–322,
    p10 l345, p13 l359–360, p15 l389–391, l408–412, l416–418, l422, p28 l711–712, p1 l17–18.
  - For l416–418, the rest of that sentence (l418–420, "best presented as reactor-specific genome evidence") still
    holds.
  - Figure 4 is the image on p24 (caption p25). I looked at it: the x-axis is labelled by bin, so "no change" is right.
  - "Unchanged on purpose" is right: p9 l329–332 have 8/6 MAGs; p13 l368 has 4/5/41/118.
- Every non-`no-change` row of `occurrences.tsv` for the PDF (40 rows: 28 Table 1 lineages + 12 text/figure rows, including the p9 l322 Patescibacteria row)
  maps to an item.
- Full-text scan of all 40 pages for R220-only names, r220, v2.7.2, Azospira and the 51 old-only IDs (Table 1
  lineages excluded): the only hits are p7 l229–230 and p8/p9 Patescibacteria. Items 1–3 cover them.
- **Table 1 (pp. 10–12).**
  - All 54 rows were parsed (bin and full lineage, rejoined across wrapped lines), and 54/54 equal HEAD r220.
  - Exactly **28** change in R232, and they are exactly the 28 rows of §1.1.
  - Each listed R232 lineage equals the repo summary.
  - The "what changes" rank-diff column is correct for 28/28.
  - The PDF-only bin (54 vs the DOCX's 53) is coasm_bin.347. It is in the list: sp017302455 → A. phosphatis.
- Numbers in sections 3–4, recomputed:
  - species 97 → 122 (26 gained, lost CAN_5_bin.70);
  - genus 273 → 283 (172 → 183 genera);
  - RED 49 → 21;
  - S6 96/70/48/77 → 121/76/28/66;
  - coasm_bin.535 *Novosphingobium* sp017302615 98.4 %/0.903;
  - coasm_bin.260 JAAEKA01 sp049294015 98.33 %;
  - Unclassified Fig 1 band 1.89 → 0.57 %. Recomputed: the blank-genus five-phase mean of `per_cell_rel_abun` goes
    from 2.14 % (18 MAGs) to 0.63 % (8 MAGs), times the ~0.88 recruitment rate. 11 MAGs gained a genus and 1
    (CAN_5_bin.110) lost one.
  - The figure-10/13 split matches the scripts (10 = denitrification bubble, R232 set; 13 = metabolic potential, r220 set).
- SI section labels were checked against the SI DOCX headings. **One error:** see section 6.

## 6. Fix applied (mechanical)

| file | what | how | before → after |
|---|---|---|---|
| `paper3/R232_change_list_for_Emily.md` (untracked, LF) | Section 2, first row: the ppk1 paragraph was labelled "SI-3". In `SI-paper3 (3)_correction.docx` that paragraph (python-docx p208) is under the heading "SI-5. Genome and ppk1-based reference-tree construction for Ca. Accumulibacter clade" (p204). SI-3 (p174) is "GAO/PAO/denitrifier classification logic". The other two labels (SI-6 p209, SI-7 p212) are right | one exact Edit: `\| SI-3, ppk1 paragraph (` → `\| SI-5, ppk1 paragraph (` (1 occurrence); copy before the edit in `apply_scratch/verify_docs/R232_change_list_for_Emily.before_verify.md` | md5 541b5204… → 085039d5…; 125 lines before and after; 1 line differs; 0 CR before and after |

The generator `apply_scripts/docs/build_change_list.py` (line 84) still says "SI-3". I did not edit the docs stage's
script. If the list is regenerated, apply the same one-word change there.

## 7. Final gates (re-run now, after all stages)

- ASV-keyed md5 (`asv_md5_before.tsv`): **1,443 / 1,443 unchanged**.
- `community_modeling/plan_r13/integrator/AMENDMENTS.sha256`, run from `community_modeling/`: **14 / 14 OK**.
- track2_guards (`track2_guards_sha256_before.txt`): **7 / 7 OK**.
- `git status --porcelain` is identical (127 lines) to the docs stage's `git_status_end.txt`.
- The only repo file newer than that snapshot is my fix to the change list.
- 47 `__pycache__` directories, as before.
- Nothing was committed or pushed.

## 8. Judgment items (reported, not changed)

1. **settings.xml** in both `_correction` DOCX (`trackChanges` → `trackRevisions`) is an untracked edit outside
   `document.xml`. It is benign; see section 2.
2. **Figure 4 in the working main DOCX.** The stale ID-labelled image was kept, with a comment holding a correct
   17-label bin-keyed map. Emily's master Fig 4 is labelled by bin and needs nothing.
3. The data-stage file `reprocess_291/skani_novelty_r232/run_skani_r232.sh` has a header saying "269 of the 291 MAGs".
   The hit table has 268 distinct `Query_name` values. This is outside the docs lens; I did not edit it.
4. The novel-bacteria module mixes sets: figure 10 uses the R232 set, while figure 13 and Tables N2–N3 use the r220
   set. This is an author decision and is listed for Emily.
5. Open author choices are carried correctly in the list: 15 vs 14 phyla; the R232 citation (the doi was not
   verified); Fig S4 relabel vs rebuild; the CAN_1_bin.98 "candidate-novel" ppk1 wording.
6. The DOCX files were validated (XSD and redlining) and rendered in LibreOffice. They were **not** opened in Microsoft
   Word.
7. A cosmetic point, not a defect. In main p56 and change-list item 2, "Myxococcota (2), Myxococcota_A (1),
   Bdellovibrionota (2)" breaks the list's descending-count order. It belongs to the open 15-vs-14 phyla choice, so it
   is Emily's call.
