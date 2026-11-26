# EmilyKin R232 rename plan: manuscript, SI, tables and figures (part "docs")

Prepared 2026-11-11, read-only. Nothing in `~/Documents/EmilyKin` (clean at 9e087c3, branch `iterativeid-nomenclature`)
or in KBase was modified. Everything below is a plan. The scripts that produced the tables are in `work/`.

Files in this folder:

| file | what it is |
|---|---|
| `occurrences.tsv` | **The main deliverable.** 3,710 rows: file, doc_status, location, bin, quote, context, current, proposed_replacement, decision, note |
| `id_change_map.tsv` | The 51 changed MAG IDs, with ASV collisions, swap partners and the 924fa18 / 6a04588 / ffe668c aliases |
| `alias_holders.tsv` | Every alias string and the bins that held it in each map era (HEAD, 408cb39, ffe668c, 6a04588, 924fa18, R232) |
| `table1_plan.tsv` | Main-text Table 1, bin-keyed: 16 cells in the working DOCX and 28 lineage rows in Emily's Oct-21 PDF |
| `table_cells_plan.tsv` | Every workbook / S6 / S7 cell that changes, bin-keyed, with the map era of the current value (14,577 rows, including superseded workbooks) |
| `table_sheet_summary.tsv` | Change counts per sheet, and which ID map each sheet uses |
| `work/` | Scan scripts and intermediates: segment extraction of 95 documents, regex scans, embedded-figure extraction, Fig S4 reference check |

`decision` values in `occurrences.tsv`:

| decision | meaning |
|---|---|
| MAG-rename | MAG ID cell or text, keyed by bin; replace with the R232 ID |
| MAG-rename (drift fix) | the cell holds the 6a04588-map ID; the bin-keyed rewrite fixes it too (the HEAD ID is unchanged by R232) |
| MAG-taxon-update | GTDB rank, lineage or species name of a MAG |
| version-update | "GTDB-Tk v2.7.2 / r220" statements and the version-mismatch caveat |
| number-update | counts derived from GTDB or IDs (phyla, case counts) |
| recompute-table | R220-derived skani / ANI / novelty values: rebuild, do not hand-edit |
| claim-revise | scientific statement that R232 contradicts (not a rename) |
| regenerate-figure | figure whose labels or bands carry IDs or GTDB names |
| ASV-keep | the string names an ASV: **do not change** |
| ambiguous-author | needs the author's call (listed in section 8) |
| no-change | checked; correct under R232 |
| superseded-no-edit | backup or superseded download; leave untouched |

---

## 1. Bottom line

1. **51 MAG iterativeIDs change, not 48.** 50 of them change root. The README §1/§5 figure of 48 came from the desk
   projection; the real R232 run (`results/iterativeid_simulation_r232.tsv`) supersedes it.
   - **21 are asv_inherit → minted** (e.g. Rhodobacter.5 → Pseudogemmobacter.2_m). The old string stays valid as the
     ASV's name, so it must not be touched in ASV contexts.
   - **4 are minted → asv_inherit.** Simplicispira.1, Holophagaceae.9, Taibaiella.5 and Hoeflea.1 are now shared with an
     ASV.
   - The other 26 are minted → minted with a new root (Paceibacteria.N_m → Minisyncoccia / Phycocordibacterales …).
2. **Every replacement must be keyed by MAG bin.** No document can be updated by plain string replacement:
   - **Swap chains.** `Chitinophagaceae.1_m` moves from coasm_bin.117 to CAN_5_bin.149, and coasm_bin.117 takes
     `Chitinophagaceae.2_m`. That string was coasm_bin.61's ID in the 6a04588 map.
   - **Era collisions.** `Ignavibacteria.21` and `.31` swapped between CAN_3_bin.221 and CAN_5_bin.167 (924fa18 vs later).
     `Ca_Competibacter.109` is coasm_bin.365 at 924fa18, CAN_4_bin.116 in the workbooks (6a04588) and CAN_2_bin.85 at
     HEAD. `Ca_Accumulibacter.8` is coasm_bin.250 in the workbooks but coasm_bin.185 in the manuscript.
   - **The workbooks use the 6a04588 map**, not HEAD. Each full sheet has 51 R232 changes plus 31 drift cells, so a
     HEAD→R232 string map misses rows such as coasm_bin.114 "Bacteroidia.2_m" and coasm_bin.466 "Ignavibacteria.27".
3. **The current master is Emily's Word file, which is not on disk.** Its newest seen state is `paper3/papert3_071626.pdf`
   (Word export 2026-10-21; gitignored; no DOCX). A Google Drive search found no paper-3 document; the only hit was the
   June "Paper 2" folder with `paper draft-PAOs.docx`.
   - The newest editable files are `papert3_070626 (2)_correction.docx` and `SI-paper3 (3)_correction.docx`. They hold
     the tracked changes of 2026-10-19 (commits 2b52feb, dce179e).
   - **Edit route:** tracked changes on a `_correction` copy, plus a page/line change list for Emily's master.
     **Do not regenerate the DOCX from the .md.** The DOCX has its own numbering, captions and collaborator text.
4. **The biggest manuscript change is not a rename.** Emily's §2.2 and Conclusions say that coasm_bin.476 "aligned poorly
   to any reference and likely represents a novel Accumulibacter lineage".
   - Under R232 it is **s__Accumulibacter sp032229525**: 98.68 % ANI, AF 0.837 to GCF_032229525.1, a genome that entered
     GTDB in R226.
   - All six Accumulibacter MAGs now have a GTDB species, and five are ≥97 % ANI to it (not three).
   - The R232 species of coasm_bin.185, coasm_bin.250 and coasm_bin.347 are *A. regalis*, *A. delftensis* and
     *A. phosphatis*. They now agree with the ppk1 clades IA, IC and IIA.
5. **Cross-part risk R1 (section 9): network membership is defined by ID equality.**
   - The network scripts' `NETWORK_MAG_ONLY` convention counts an ASV as "MAG-associated" if its ID is a value of
     `mag_iterativeID_old_to_new.json`.
   - Renaming MAGs therefore silently changes the MAG-node set of the SPIEC-EASI and co-occurrence networks: 21 ASVs drop
     out and 4 join, among them Taibaiella.5, currently a "flanking" node in SPIEC-EASI module 1.
   - Freeze or redefine this by bin before any network figure is re-rendered.

---

## 2. The ID change set (inputs)

- **Source.** `results/iterativeid_simulation_r232.tsv` (repo rule `scripts/recompute_mag_ids.py`). It reproduces the
  committed map 291/291 from R220.
- **Cases.** R220: 155 inherit / 62 genus / 74 lowest. R232: **138 / 70 / 83**. No duplicate IDs in the R232 map.
- **Full list.** `id_change_map.tsv`. Columns: old_is_asv_id (21 True), new_is_asv_id (4 True), swap partners, and the
  aliases from earlier map eras.
- **ASV-ID collisions that matter for text:**
  - Old IDs that remain ASV names, so leave ASV uses alone: Ca_Falkowbacteria.7, Ignavibacteria.5/.13/.21/.31/.7/.25,
    Rhodobacter.13/.5, Azospira.3, Sphingobacteriales.16/.13/.9, Gammaproteobacteria.7, Hydrogenedensaceae.8,
    Rhodobacteraceae.1, Paracoccus.4, Myxococcaceae.11, Ca_Competibacter.32, Bosea.3, Mesorhizobium.7.
  - New IDs that already exist as ASV names: Simplicispira.1, Holophagaceae.9, Taibaiella.5, Hoeflea.1.
    - Any document that already mentions these ASVs keeps them.
    - After the rename the same string also names a MAG, which is by design (asv_inherit).
- **Substring and suffix hazards.**
  - Boundaries: `X.1` vs `X.1_m` vs `X.13`, e.g. Simplicispira.1 vs Simplicispira.1_m. Use the lookarounds
    `(?<![A-Za-z0-9_])ID(?![0-9A-Za-z_])`.
  - Prefix forms: `Ca_`, `Ca.`, `*Ca.*` and `Candidatus` (e.g. "Ca. Competibacter.32").
  - Abbreviated lists escape any ID regex. Examples in SPIEC-EASI_network_analysis.docx: "Ca. Competibacter
    (.26, .32, .51)" and ".26/.32/.49/.51/.109".
  - In prose with no bin key, use a two-pass replace (old → unique placeholder → new) so that a swap chain cannot cascade.
- **Partial genus renames.** A genus can change for some MAGs only. Examples: JJ008 (coasm_bin.492 → CADEDV01;
  CAN_5_bin.112 stays JJ008), Flavobacteriales (3 of 7 → _E/_F), Myxococcota (1 of 3 → Myxococcota_A). A name replace
  without a bin is wrong for these.

---

## 3. Which document is the master, and how each should be edited

| document | status | edit route |
|---|---|---|
| `paper3/papert3_071626.pdf` | Emily's newest content (2026-10-21, 40 pp). It partially merges the Oct-19 fixes (memory `paper3-071626-review-findings`). No DOCX on disk | **Send a change list** (all `paper3/...pdf` rows of `occurrences.tsv`, keyed by page and manuscript line). Ask Andrew/Emily for the current DOCX, copy it to `<name>_correction.docx` and apply the same rows as tracked changes |
| `papert3_070626 (2)_correction.docx` | newest editable main text; tracked changes by "Claude (A. Freiburger)" | tracked changes in place, after a `.pre_r232.bak.docx` backup. Use the existing helper (memory `manuscript-docx-integration`) |
| `SI-paper3 (3)_correction.docx` | newest editable SI | as above |
| `MANUSCRIPT_methods_results_discussion.md` | Claude's git-tracked long-form source of truth. Not the submission text | direct edits (section 4.4) |
| `MAG_target_gene_filled_UPDATED (1).xlsx` | **submission workbook** (S1–S7, Notes_per_cell). Hand-edited per-cell, s16_equiv and 5-class columns with no generator | **bin-keyed surgical cell edits** with a `.pre_r232.bak.xlsx` backup. Never run `build_xlsx_291.py` or `build_mag_quality_table.py` against it |
| `MAG_target_gene_filled.xlsx` | builder workbook. `render_denitrification_heatmaps.py` and `build_table_s7_classification.py` read it | rebuild through the builders after the inputs are R232 (section 6), or apply the same bin-keyed cell plan |
| `Table_S6_skani_novelty.tsv`, `Table_S7_MAG_functional_classification.csv` | standalone S6/S7 | regenerate (section 6) |
| `cladeexcel.xlsx` (S4 source), `abundance and coverage (1).xlsx` | collaborator workbooks with ID columns | bin-keyed ID / genus cell edits |
| `SPIEC-EASI_network_analysis.docx` | supporting write-up (SI-8 source) | no rename (ASV labels); 3 author items (section 8) |
| `novel_bacteria_supplementary.md` | 276-era SI module, partly stale | update only if kept (section 7) |
| brief pptx / program_overview html | community-modelling presentation | no rename (ASV context); one optional speaker-note edit |

**Do not edit:**
- every `*.pre_*.bak.*`;
- the `(1)/(2)/(3)` re-downloads (`papert3_070626*.docx` except `(2)_correction`, `SI-paper3*.docx` except
  `(3)_correction`);
- `CAN_denitrification_manuscript.docx`, `wastewater manuscript_draft1.docx`, `Wastewater manuscript-SI.docx`;
- `MAG_target_gene_filled (1).xlsx` (276-era), `MAG_target_gene_filled_UPDATED.xlsx`, `MAG_target_gene_empty_tables.xlsx`,
  `abundance and coverage.xlsx` (no ID column);
- `TableS7PAOGAOdenitrifierlogic.docx`, `PAOGAOdenitrifierlogic.docx`, `Andrewgraphs.pptx`,
  `Metagenomic_Pipeline_Overview*.pptx`;
- the historical hand-off notes `SERVER_HANDOFF_S6_ani_novelty.md` and `MANUSCRIPT_EDITS_MIMAG_rRNA_tRNA.md`.

They are counted as `superseded-no-edit`.

**DOCX mechanics:**
- Use python-docx with the repo's `track_replace` helper. Python-docx `.text` hides `<w:ins>`, so verify through the raw
  XML or `soffice --headless --convert-to pdf`.
- **Gotcha:** a second `track_replace` in the same paragraph can wipe the first if the multi-run fallback triggers.
  Consolidate the target paragraph to a single run first.
- Table 1 cells are single-paragraph cells: edit only the iterativeID and GTDB genus cells listed in `table1_plan.tsv`.

---

## 4. Per-document change list (summary; every row is in `occurrences.tsv`)

### 4.1 Emily's master (`paper3/papert3_071626.pdf`)

| where | current | proposed | decision |
|---|---|---|---|
| p7 l229–230 | "GTDB-Tk v2.7.2 against … release r220 (50, 51)" | "GTDB-Tk v2.8.0 (classify_wf) against … release R232"; add the R232 citation, Parks et al. 2025 NAR doi 10.1093/nar/gkaf1040 | version-update |
| p8 l307–310 | "14 bacterial phyla … Patescibacteria (21) … Myxococcota (3)" | "15 … Patescibacteriota (21) … Myxococcota (2), Myxococcota_A (1)". Alternatively keep 14 by merging _A, which is an author decision | number-update |
| p9 l322 | "Patescibacteria (0.67–1.17 Mb)" | Patescibacteriota | MAG-taxon-update |
| p9 l330/332 | "eight MAGs, five high-quality" / "six MAGs, four high-quality" | no change: the GTDB genus counts are identical in R232 | no-change |
| p10 Table 1 | 28 of 54 lineage strings change (`table1_plan.tsv`) | R232 lineages; caption adds "(GTDB R232)" | MAG-taxon-update |
| p13 Fig 1 | R220 genus bands | re-rendered figure (section 5) | regenerate-figure |
| p13 l359–361 | "48 … ANI <85% … a further 77 …" | recompute from the R232 skani rebuild of Table S6 | recompute-table |
| p15 l389–391 | "three … ≥97% ANI … coasm_bin.476 … aligned poorly … novel Accumulibacter lineage" | revise the claim: five are ≥97 % (185 97.58, 250 98.60, 347 98.14, 98 99.48, 476 98.68; CAN_4_bin.64 96.05). coasm_bin.476 = *Accumulibacter* sp032229525 | claim-revise |
| p15 l408–412 | "falls near Ca. A. delftensis … sits near proposed Ca. A. regalis" | optionally say the R232 GTDB species now agree (regalis / delftensis / phosphatis) | claim-revise (author) |
| p15 l416–420 | "strongest novelty point … aligned poorly to any reference" | keep the ppk1-level statement (nearest ppk1 BPBW2884/BPBW2899); drop the genome-level novelty | claim-revise |
| p15 l422 | "nearest GTDB references are largely metagenome-derived drafts" | re-check against the R232 representatives | claim-revise (author) |
| p28 l712 | "a novel near-complete Clade ID lineage MAG" | e.g. "a near-complete Clade ID MAG (GTDB R232 *Accumulibacter* sp032229525) whose ppk1 is nearest only to sparsely characterized references" | claim-revise |

**Table 1 in Emily's PDF** has no iterativeID column; it prints full lineages. 28 rows change:
- p__Patescibacteria → p__Patescibacteriota in 4 rows (CAN_1_bin.48, CAN_1_bin.195, CAN_1_bin.174, CAN_5_bin.90);
- c__JAEDAM01 → c__Gracilibacteria, c__Paceibacteria → c__Minisyncoccia, c__Kapaibacteriia → c__Ignavibacteria,
  c__Verrucomicrobiae → c__Verrucomicrobiia;
- o__SJA-28 → o__Tepidaquicellales; f__Burkholderiaceae → f__Burkholderiaceae_C (Ottowia, Hydrogenophaga, Inhella);
- o__Rhodospirillales_A → o__Rhodospirillales and f__UXAT02 → f__Defluviicoccaceae (coasm_bin.88; genus Defluviicoccus
  unchanged);
- genus renames SCMT01, Flavipsychrobacter, JARSSO01, Azospira_A, UBA919, and genus gains for CAN_1_bin.48, CAN_1_bin.2
  and CAN_5_bin.2;
- species gains (UBA2287, Leptovillus, OLB9, PGYR01, Dokdonella, PHOS-HE28) and species renames (Accumulibacter ×3);
- the Azonexus sp017302355 species loss for CAN_5_bin.70.

### 4.2 Main working DOCX (`papert3_070626 (2)_correction.docx`)

- **para 45.** "GTDB-Tk v2.7.2 … release r220" → "v2.8.0 … R232" (version-update).
- **para 57.** Phylum sentence, as in the PDF (number-update).
- **para 58.** Competibacter 8/4 and Accumulibacter 6/4: no change.
- **para 59 (Table 1 legend).** "GTDB-Tk (r220) genus" → "GTDB-Tk (R232) genus".
- **Table 1** (bin-keyed; 8 iterativeID cells + 8 GTDB-genus cells):

  | row | bin | iterativeID | GTDB genus |
  |---|---|---|---|
  | 4 | CAN_3_bin.221 | Ignavibacteria.21 → **Tepidaquicellales.1_m** | OLB5 (no change) |
  | 11 | CAN_5_bin.177 | unchanged | SCMT01 → **Lysobacter_I** |
  | 15 | CAN_4_bin.208 | Azospira.3 → **Proximibacter.2_m** | Azospira_A → **Proximibacter** |
  | 16 | coasm_bin.365 | Ca_Competibacter.32 → **Contendibacter.1_m** | JARSSO01 → **Contendibacter** |
  | 17 | CAN_1_bin.2 | Ignavibacteria.5 → **Tepidaquicellaceae.1_m** | f__B-1AR → **JAFDZO01** (genus gained) |
  | 18 | CAN_1_bin.195 | Ca_Falkowbacteria.7 → **Patescibacterium.1_m** | UBA919 → **Patescibacterium** |
  | 19 | CAN_5_bin.2 | unchanged | f__UBA2400 → **JBDHTQ01** (genus gained) |
  | 21 | CAN_1_bin.48 | Paceibacteria.1_m → **Minisyncoccia.1_m** | f__GWB1-50-10 → **CAYHLH01** (genus gained) |
  | 22 | coasm_bin.40 | Flavipsychrobacter.6_m → **Taibaiella.5** | Flavipsychrobacter → **JAFIGV01** |
  | 54 | coasm_bin.337 | Ignavibacteria.7 → **Tepidaquicellaceae.3_m** | UBA9002 (no change) |

  The working DOCX's Table 1 is the 53-MAG version. Emily's PDF has 54 MAGs with a lineage column and no ID column.
  Re-check which layout survives before editing.
- **para 656.** The same Accumulibacter novelty revision as the PDF, plus 48/77 from the S6 rebuild.
- **para 657.** Guild counts 4/5/41/118 are expected unchanged; confirm after the S7 rebuild.
- **Figures 1, 4 and the module panel (image6):** section 5.

### 4.3 SI working DOCX (`SI-paper3 (3)_correction.docx`)

- **para 209 (ppk1 text).** Optionally add that the R232 species agree with the ppk1 clades. Note coasm_bin.476 =
  *Accumulibacter* sp032229525. CAN_1_bin.98 is unchanged (sp014879755, 99.48 %).
- **para 212.** "Giesbergeria CAN_5_bin.40" → *Simplicispira*. "closest GTDB r220 reference genomes (Figure S4)" is an
  author item (section 8).
- **para 216 (iterativeID method).** "155 ASV-inherited, 62 GTDB-genus, and 74 GTDB-lowest" → **138 / 70 / 83**. Add
  "on GTDB R232".
- **para 252 (Fig S4 legend):**
  - "rooted on the Patescibacteria clade" → Patescibacteriota;
  - "Azospira_A representative … (Azospira_A restricta, GCF_016858125.1)" → *Proximibacter* (*Proximibacter
    restrictus*);
  - "best skani ANI to a GTDB r220 representative (Table S6)" → the release of the rebuilt S6;
  - "61 GTDB r220 species-representative reference genomes" is an author item (keep it if the tree is only relabelled).
- **para 255 (Fig S5 legend).** Text unchanged; the image is re-rendered (CAN_3_bin.221 label).
- **Tables S1/S10/S11** in the SI DOCX carry no MAG IDs or taxa: no change.

### 4.4 `MANUSCRIPT_methods_results_discussion.md` (direct edits)

| line | change |
|---|---|
| 45 | GTDB-Tk **v2.8.0** `classify_wf`, local run (`--batchfile … --cpus 24 --pplacer_cpus 8 --min_perc_aa 10 --scratch_dir …`, `release232` data, skani v0.4.0 ANI screen; verified in `results/gtdbtk.log`). Method counts 145 topology / 97 ANI / 49 RED → **51 topology, 97 topology+ANI, 122 ani_screen, 21 RED** (R232 labels). **Delete** the version-mismatch caveat ("As a caveat, GTDB-Tk v2.7.2 expects … reported as such throughout") |
| 61 | 155 / 62 / 74 → **138 / 70 / 83**; "51 MAG IDs changed; ASV IDs unchanged" |
| 87, 153 | *Giesbergeria* (CAN_5_bin.40) → *Simplicispira*. Line 153's **JJ008 stays**: it is CAN_5_bin.112 |
| 103 | skani v0.3.2 / GTDB r220 skani database → the R232 database and the skani version used in the S6 rebuild |
| 107 | Fig S4 methods: *Azospira_A* → *Proximibacter*, *Azospira_A restricta* → *Proximibacter restrictus*. "r220" for the tree inputs (references, bac120 mask, `bac120_msa_reps_r220.faa`, the GTDB-Tk v2.7.2 identify/align re-implementation) is factual tree history. Keep it unless the tree is rebuilt (section 8) |
| 123 | "GTDB-Tk (r220)" → R232. **15 phyla**, Patescibacteria → **Patescibacteriota (21)**, Myxococcota (3) → **Myxococcota (2) + Myxococcota_A (1)**. Rank resolution: family 290 (no change), genus **283 (183 distinct genera)**, species **122 (169 unnamed)**, RED novelty 49 → **21**. skani buckets 96/70/48/77 → from the rebuilt Table S6 |
| 141 | 155 / 62 / 74 → **138 / 70 / 83** |
| 147 | coasm_bin.185 (*Accumulibacter* sp000585075) → *Accumulibacter regalis* |
| 157, 225 | CAN_1_bin.98 unchanged; add coasm_bin.476 = *Accumulibacter* sp032229525 (98.68 %) wherever ID novelty is implied |
| 181 | **Delete** the Discussion caveat "Finally, the GTDB taxonomy is r220 although the GTDB-Tk binary (v2.7.2) natively expects r232 …" |
| 191 | Fig 1 legend genus list: replace with the list from the re-rendered figure. Expected: SCMT01 → Lysobacter_I, Giesbergeria → Simplicispira, Azospira → Proximibacter, UBA5532 phylum note → Patescibacteriota; JJ008 loses coasm_bin.492; order may change |
| 229 | Fig S4 legend: Patescibacteria → Patescibacteriota. "no skani hit to any GTDB r220 representative" follows the rebuilt S6. "61 GTDB r220 references" is an author item |
| 231 | Table S6: r220 → R232, and the recomputed buckets |

The .md contains none of the 51 old IDs as text (only Accumulibacter.8 and Competibacter.85, both unchanged). The 155/62/74
numbers and the GTDB statements are its only rename-related edits.

### 4.5 Workbooks (bin-keyed cells; `table_cells_plan.tsv`)

| file / sheet | ID cells (R232 + drift) | rank cells | lineage cells | other |
|---|---|---|---|---|
| `MAG_target_gene_filled_UPDATED (1).xlsx`, S1_MAG quality | 51 + 31 | 136 changed + 12 gained | 130 (`Taxonomy ` col) | |
| … S2a_RPKM per MAG | 31 + 23 | | | |
| … S3_gene copies | 51 + 31 | | 130 (`Taxonomy`) | |
| … S4_nosZ_genecopy_type | 27 + 15 | 22 (genus) | | |
| … S5 _ppk_ACC | 0 + 2 (coasm_bin.185 Ca_Accumulibacter.6 → .8; coasm_bin.250 .8 → .2) | | | |
| … S6_ANI | 51 + 31 | 153 + 38 gained (Species col) | 130 | 2,436 R220 skani/ref/provenance cells → **recompute** |
| … S7_functional | 51 + 31 | 94 + 12 gained | | |
| `MAG_target_gene_filled.xlsx` (S1, S1b, S2a, S3, S4, S5, S6, S7) | same pattern | | | 676 S6 cells → recompute |
| `cladeexcel.xlsx` S4 | 27 + 15 | 22 | | |
| `abundance and coverage (1).xlsx` Sheet1 | 51 + 31 | | | |
| `Table_S7_MAG_functional_classification.csv` | 51 + 31 | 106 | | |
| `Table_S6_skani_novelty.tsv` | (no ID column) | | 129 (+1 coasm_bin.221 genus mismatch) | 1,760 skani cells → recompute |

- **Rule.** For every row whose first column is a bin, set the ID column to the R232 ID from
  `results/iterativeid_simulation_r232.tsv` and the rank/lineage columns from `results/gtdbtk.r232.summary.tsv`.
- **Never look up by the current ID string.** Header rows and notes that say "GTDB R220 (GTDB-Tk v2.7.2, classify_wf)"
  (S6 `Reference database` column) become "GTDB R232 (GTDB-Tk v2.8.0, classify_wf)".
- **Formatting.** Keep formulas (`=AVERAGE(...)` in S1) and styles: load with openpyxl `data_only=False`, write values
  only into the planned cells.

### 4.6 Supporting documents

- **`SPIEC-EASI_network_analysis.docx`.** Every MAG-like label is an **ASV iterativeID**:
  `render_spieceasi_mag_network.py` says "label = ASV iterativeID (+MAG bin)". So these stay as they are: Azospira.3,
  Ignavibacteria.5/.13, Ca. Competibacter.32/.109, the abbreviated ".32" lists, and "Patescibacteria" / "Myxococcota"
  (SILVA names in the module and lineage tables). Three sentences call Azospira.3 a "MAG" (paras 140 and 357, plus the
  abbreviated list in para 140); see section 8.
- **`novel_bacteria_supplementary.md`.** None of its 7 MAGs changes ID in R232. coasm_bin.260 (Ca_Sarcinithrix.1) gains
  s__JAAEKA01 sp049294015, so "seven species-novel" becomes 6 of 7. The module is 276-era: Dokdonella.14 = CAN_1_bin.210,
  which is Dokdonella.8 at HEAD.
- **Brief deck and program overview.** "Ca_Competibacter.109 … top member on 64 of 75 days" is the 16S by-day ASV, and
  CAN_2_bin.85's MAG ID is unchanged anyway, so no edit. Slide 13 notes say "five Patescibacteria" (genetic-code-25 MAGs)
  → optionally Patescibacteriota. The deck is generated by `make_brief_pptx.py`; edit `figure_data.json` / the script,
  not the pptx.

---

## 5. Figures

| figure (where) | generator | what shows R220 / old IDs | action |
|---|---|---|---|
| Main Fig 1, per-cell composition (DOCX image1; PDF p13) | same generator as `figures/percell_composition.{png,pdf}`, judging by the legend format. The embedded image is **not** byte-identical to the committed PNG (no md5 match for any main-DOCX image), so regenerate rather than assume the repo file is what Emily embedded: `scripts/render_percell_composition.py` (genus from `mag_quality_summary.tsv`, suffix-stripped; "Ca_" prefix if an iterativeID root is Ca_<genus>) | legend: SCMT01, JJ008, Azospira, Giesbergeria, "UBA5532 (Patescibacteria)" | **regenerate** once `mag_quality_summary.tsv` carries R232. Bands may merge or split (JJ008 → JJ008 + CADEDV01; Azospira_A → Proximibacter) |
| Main Fig 2/S3 chords (image2; SI image3) | network/chord scripts | none | no change |
| Main Fig 3 heatmap + bars (image3) | (guild-grouped, unlabelled) | none | no change unless guilds change |
| **Main Fig 4**, denitrifier-class gene heatmap (image4) | Emily-made 5-class version; the closest repo generator is `scripts/render_denitrification_heatmaps.py` (4-facet), which reads the `MAG iteration ID` column of the builder workbook | x-axis **MAG iterativeIDs in the 6a04588 map** (visible: Ca_Falkowbacteria.7, Giesbergeria.*, Azospira.*, Gammaproteobacteria.*, Ignavibacteria.*, Lysobacter_B.*, Flavobacteriales.*) | **regenerate** after the workbook ID column is rewritten. Emily must re-render her version |
| Main Fig 5, class × phase heatmap (image5) | — | none | no change |
| DOCX image6 (`module_denitrification_genes_phaseIII`, pre-408cb39 render) and image7 (`module_pathway_completeness_phaseIII`) | `scripts/render_module_denitrification_genes.py` | image6 x-axis MAG IDs; both depend on MAG-node membership | **Not in Emily's Oct-21 master.** That PDF has main Figures 1–5 only and no co-occurrence-module or Louvain panel. They are legacy images left in the Oct-19 working DOCX: no regeneration is needed for the master. Regenerate only if they are reinstated, and only after risk R1 is resolved |
| SI Fig S1, ppk1 tree | `ppk1_classify/…` | bin names only | no change |
| SI Fig S2, framework diagram | hand-made | literature taxa | no change |
| **SI Fig S4**, 65-MAG genome tree (`fig1_tree/fig1_tree.png`) | `scripts/render_fig1_tree.py` (+ `fig1_annotation.tsv`) | 8 MAG tips change (Ca_Competibacter.32, Azospira.3, Flavipsychrobacter.6_m, Ignavibacteria.7/.5/.21, Ca_Falkowbacteria.7, Paceibacteria.1_m). 8 reference names change (JARSSO01 → Contendibacter, SCMT01 → Lysobacter_I, sp000585075 → regalis, phosphatis_G → delftensis, Azospira_A restricta → Proximibacter restrictus, Niabella → defluvii, Flavipsychrobacter → JAFIGV01, UBA919 → Patescibacterium). 1 reference (GCA_017302355.1, Azonexus sp017302355) left GTDB. Phylum colours / legend: Patescibacteria (7 tips) | **fix `render_fig1_tree.py:137` and `fig1_tree/compare_trees.py:96`** (`== "Patescibacteria"`; with R232 they silently fall back to midpoint rooting), regenerate `fig1_annotation.tsv` from R232, re-render. A full rebuild with R232 references is an author decision |
| **SI Fig S5**, nosZ top-5 (`clade_classify/out/nosz_top5_tree.png`) | `clade_classify/render_nosz_top5_tree.py` (reads `mag_iterativeID_old_to_new.json`) | rank-5 MAG CAN_3_bin.221 labelled Ignavibacteria.21 | **re-render** → Tepidaquicellales.1_m (genus OLB5 unchanged) |
| SPIEC-EASI figures (`spieceasi_mag_network.png`, `spieceasi_modules_flanking_89.png`) | `network/compositional_networks/render_*` | ASV labels; the MAG-node set is defined by ID equality | labels unchanged; **do not re-render before R1 is fixed** |
| Other repo PNGs (co-occurrence grids, chord RPKM, module panels, gene_ab_figure 10–13) | listed in README §5 | ID labels | owned by the figures part. The same R1 caveat applies to every `NETWORK_MAG_ONLY` grid |

---

## 6. Tables S6 and S7: how they are built and what to run

**Table S6 (ANI novelty).** It is R220 throughout: references, ANI/AF and bucket counts 96 / 70 / 48 / 77.

- **`reprocess_291/build_skani_novelty_291.py`** writes `Table_S6_skani_novelty.tsv` (9 columns).
  - Inputs, as hard-coded `/scratch1` paths: the skani search output, `gtdb_taxonomy_r220.tsv` and the gtdbtk summary.
  - Rule: best hit = highest ANI with AF_query ≥ 15 %; classes ≥95 / 85–95 / <85 / low_AF_unresolved.
  - **To do:** point it at an R232 skani search, the R232 taxonomy (`gtdb_downloads/r232/bac120_taxonomy_r232.tsv.gz`)
    and `results/gtdbtk.r232.summary.tsv`.
  - `skani_preview/skani_search_r232.tsv` already holds raw R232 hits, but it was run with GTDB-Tk's settings. Its
    `--min-af` and number of hits per query may differ from the S6 recipe (`--min-af 0 -n 3`), so re-run skani with the
    S6 settings against `release232/skani/database`, or confirm the preview has enough hits.
- **`scripts/build_table_s6_ani.py`** fills the `Table S6_ANI` sheet of `MAG_target_gene_filled.xlsx`, not the
  submission workbook.
  - Local stage: ID (from `mag_iterativeID_old_to_new.json`), lineage and ranks (from `mag_quality_summary.tsv`) and
    "Reference database".
  - ANI stage: GTDB-Tk summary + skani.
  - Novelty rule: ANI 95 % / AF 0.60.
  - Run it after the ID map and `mag_quality_summary.tsv` are R232; pass `--gtdbtk results/gtdbtk.r232.summary.tsv` and
    `--skani <R232 search>`.
- **Submission workbook S6_ANI** (27 columns, including the skani_* block): copy the rebuilt columns in, bin-keyed.

**Table S7 (functional classification).**
- **`scripts/build_table_s7_classification.py`** reads S1 and S3 of `MAG_target_gene_filled.xlsx`. It writes
  `Table_S7_MAG_functional_classification.csv`, the workbook sheet and `mag_guild_classification.csv` (the guild source of
  truth for figures).
- **Taxonomy-dependent tiers:**
  - Tier 1 GAO = `Competibacter|Contendobacter` genus OR family Competibacteraceae.
  - Tier 2 PAO/GAO genus lists.
- **Check:** coasm_bin.365's R232 genus is "**Contendibacter**", which the MiDAS-spelled regex `Contendobacter` does not
  match. It stays canonical GAO through family Competibacteraceae, so no guild change is expected.
  - No other changed genus or family hits the Tier 1/2 lists. This was verified over all 291 R220 → R232 lineages: the
    Tier-1/Tier-2 taxon flags do not change for any MAG. Defluviicoccus is a genus that is unchanged; only its family
    UXAT02 → Defluviicoccaceae changes.
  - Confirm by diffing the rebuilt CSV against HEAD. Only Genus/Family/Phylum and iterative_ID should change.
- **The S7 CSV / guild CSV currently carry the 6a04588 ID map** (31 drift rows). Rebuilding from an R232 S1 fixes this.

**`reprocess_291/build_xlsx_291.py`** rebuilds S1/S2a/S2b/S3/S4/S5/S6 data rows of the builder workbook from the ID map,
`mag_quality_summary.tsv` and the gtdbtk summary.
- Its header comment `R=/scratch1` path must be patched to the repo root (`REPO`).
- It has a surgical mode for S3 nosZ. **Do not run it against the submission workbook.**

**`reprocess_291/compute_manuscript_stats_291.py`** prints the manuscript ground-truth numbers. It reads
`gtdbtk.bac120.summary.tsv`, `mag_iterativeID_old_to_new.csv` and `Table_S6_skani_novelty.tsv`.
- Patch `R="/scratch1/..."` → `R="."` and guard the server-only `mag_16s.fasta` / Bakta glob.
- Re-run it after the R232 inputs land to confirm the section 4.4 numbers: phyla, case counts and skani buckets.

---

## 7. Scientific consequences beyond renaming (for Andrew and Emily)

1. **coasm_bin.476 is no longer novel at species level.**
   - R232: s__Accumulibacter sp032229525, 98.68 % ANI, AF 0.837. The reference GCF_032229525.1 was added in R226.
   - The R220 Table S6 row says "low_AF_unresolved … likely novel".
   - This hits the PDF (l389–391, 416–420, 712), DOCX para 656, the SI ppk1 text and possibly the abstract. The ppk1-level
     statement (nearest ppk1 BPBW2884/BPBW2899, clade ID) remains valid.
2. **Accumulibacter species names now match the ppk1 clades.** coasm_bin.185 *regalis* (IA), coasm_bin.250
   *delftensis* (IC), coasm_bin.347 *phosphatis* (IIA). This strengthens rather than contradicts the text.
3. **Species-level resolution rises from 97 to 122** (26 gained, 1 lost: CAN_5_bin.70). RED-novelty falls from 49 to 21.
   - The ".md" sentence "only about a third of the catalogue is confidently conspecific" and the S6 buckets must be
     recomputed.
   - The novel-bacteria SI module loses one of its seven (coasm_bin.260).
4. **Phylum count** becomes 15 if Myxococcota_A is counted separately (author convention).
5. **coasm_bin.365.** Its ID changes from Ca_Competibacter.32 to Contendibacter.1_m.
   - The MAG was never in the GTDB Competibacter genus (R220 JARSSO01), so the "eight Competibacter MAGs" and
     "five canonical GAO" counts are unchanged.
   - Text that groups ".32" with Competibacter by its ID would now read differently.

---

## 8. Ambiguous items (need an author decision)

1. **SPIEC-EASI "Azospira.3" called a MAG.** Paras 140 and 357: "gene-content PAO-like MAG (Azospira.3, …)", and the
   "(.26, .32, .51)" list. The node is the ASV. The MAG CAN_4_bin.208 becomes Proximibacter.2_m, and coasm_bin.365's MAG
   becomes Contendibacter.1_m.
   - Proposal: keep the ASV labels and write "Azospira.3 (ASV; MAG CAN_4_bin.208, Proximibacter.2_m)", or call them
     "nodes".
2. **Fig S4: relabel or rebuild?** The tree was inferred with 61 R220 references and an R220 bac120 mask.
   - Relabel only: update tip and reference names to R232 and keep "r220 reference genomes" in the captions (.md
     l107/229, SI para 212/252). One reference (GCA_017302355.1) has no R232 name, so label it "Azonexus sp017302355 (R220)".
   - Full rebuild with R232 representatives: change "r220" everywhere.
3. **Phylum count convention.** 15 phyla (Myxococcota_A separate) or 14 (suffix-merged). Fig 1 already strips genus
   polyphyly suffixes.
4. **Table 1 layout.** The working DOCX has 53 MAGs with iterativeID + GTDB genus columns. Emily's PDF has 54 MAGs with
   full lineages and no IDs. Edit whichever survives in her master.
5. **Accumulibacter species wording** (PDF l408–412, SI para 209): whether to cite the new GTDB species alongside the
   ppk1 clades.

---

## 9. Risks and cross-part dependencies

- **R1, network MAG-node membership by ID equality (highest).**
  - Where: `scripts/render_cooccurrence_grid.py` (`NETWORK_MAG_ONLY`), `render_cooccurrence_grid_per_phase.py`,
    `network/compositional_networks/render_spieceasi_mag_network.py`, `spieceasi_modules_flanking.py`.
  - Their rule: an ASV is a MAG node if its ID ∈ values(`mag_iterativeID_old_to_new.json`).
  - After the R232 map is written, 21 ASVs lose MAG status (Rhodobacter.5/.13, Azospira.3, Ignavibacteria.*,
    Ca_Competibacter.32, …) and 4 gain it (Simplicispira.1, Holophagaceae.9, Taibaiella.5, Hoeflea.1).
  - Effect: the SPIEC-EASI "53 MAG-resolved nodes", the module "(MAG+flank)" counts, the guild-edge tables and every
    MAG-only co-occurrence grid would change for a naming reason only.
  - **Fix before the new map lands:** define membership by bin, via `asv_to_mag_mapping.json` non-weak tiers or a frozen
    ASV→bin list taken from the HEAD map.
- **R2, workbook map era.** All workbook sheets, `Table_S7_*.csv`, `mag_guild_classification.csv` and Emily's
  `abundance and coverage (1).xlsx` use the 6a04588 map. A HEAD-based string map would corrupt the swap-chain IDs
  (Ca_Accumulibacter.6/.8/.2, Ca_Competibacter.48/.109/.49/.57/.85, Chitinophagaceae.2_m, Dokdonella.17/.21,
  Ca_Leptovillus.9/.11, Ca_Sarcinithrix.1/.2). Use bin-keyed edits only.
- **R3, one-shot generators.** `regen_iterative_ids.py`, `apply_new_ids*.py` and `assign_mag_iterative_ids.py` write an
  identity old→new map if re-run. Produce the R232 map from the simulation output. `apply_mag_ids.py`'s bin-keyed
  `set_mag_col` is the right pattern for data files; it reads `scripts/_new_ids.json`.
- **R4, hard-coded "Patescibacteria" in code.**
  - `render_fig1_tree.py:137`, `fig1_tree/compare_trees.py:96`: fall back to midpoint rooting.
  - `community_modeling/analysis/r13/annotations/build_taxonomy.py` (`FLAG`, `DIDERM` with Myxococcota): drops the CPR
    flag and the diderm call.
  - Owned by the code part, but Fig S4 depends on it.
- **R5, S7 regex spelling.** `Contendobacter` vs GTDB "Contendibacter". Harmless today through the family rule; add
  `Contendibacter` to be explicit.
- **R6, `mag_abundance_summary.tsv` vs gtdbtk summary.** coasm_bin.221's R220 genus `1-14-0-10-56-10` is blank in the
  summary (the only one of 291 lineages that differs). R232 maps it to CAZOZK01. Rebuild the taxonomy columns from the
  R232 summary rather than patching the old ones.
- **R7, Emily's partial merges.** Her next download may drop tracked fixes, as the Oct-21 PDF did. Re-diff every new
  download against `occurrences.tsv` before closing the task.
- **R8, KBase object names.** ws 262319 Assemblies, RAST/KOfam GenomeSets and r13 model inputs carry the old IDs (other
  part). Manuscript text never cites KBase object names, so documents do not depend on that order.

---

## 10. Execution order (after Andrew approves)

1. **Data first (other parts).**
   - R232 taxonomy into `mag_abundance_summary.{tsv,json}` and `mag_quality_summary.tsv`, bin-keyed.
   - The R232 ID map (`mag_iterativeID_old_to_new.{csv,json}` + `scripts/_new_ids.json`) from
     `results/iterativeid_simulation_r232.tsv`.
   - **Before that map is written, fix R1.**
2. **Table S6.** R232 skani search → patched `build_skani_novelty_291.py` → `Table_S6_skani_novelty.tsv`; then
   `build_table_s6_ani.py --gtdbtk results/gtdbtk.r232.summary.tsv --skani …` → builder workbook S6.
3. **Builder workbook.** Patched `build_xlsx_291.py` (S1–S5 rows) → `build_table_s7_classification.py` (S7 csv, sheet,
   guild csv). Diff against HEAD: only ID / taxonomy / S6 columns may change.
4. **Submission workbook.** Bin-keyed cell edits from `table_cells_plan.tsv`, rows where file = `MAG_target_gene_filled_UPDATED (1).xlsx`.
   - Copy S6 columns from step 2.
   - Leave per-cell / s16 / 5-class columns untouched.
   - Back up to `.pre_r232.bak.xlsx`.
   - Same for `cladeexcel.xlsx` and `abundance and coverage (1).xlsx`.
5. **Stats.** Patched `compute_manuscript_stats_291.py` → confirm the section 4.4 numbers.
6. **Figures.**
   - Fig 1: `render_percell_composition.py`.
   - Fig S4: fix the rooting string, regenerate `fig1_annotation.tsv`, run `render_fig1_tree.py`.
   - Fig S5: `render_nosz_top5_tree.py`.
   - Fig 4 and module panels: after the workbook (Fig 4) and after R1 (module panels). Emily re-renders her Fig 4 version.
7. **Text.**
   - `.md` edits (section 4.4).
   - Tracked changes in the two `_correction` DOCX files (sections 4.2–4.3), with images replaced from step 6.
   - PDF-keyed change list for Emily's master (section 4.1 + `table1_plan.tsv`).
   - Optional: SPIEC-EASI and novel-bacteria items after the author decisions in section 8.
