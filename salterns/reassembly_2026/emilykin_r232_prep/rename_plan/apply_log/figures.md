# Figures stage: data/PLAN.md step C (C1–C6) and docs/PLAN.md figure steps 8–11, with the C22 PDFs

Run on 2026-11-12, 13:44–14:05 CST, in `~/Documents/EmilyKin` (branch `iterativeid-nomenclature`, HEAD 9e087c3). It
started on top of the data stage's uncommitted tree: 80 `git status` lines.
- Nothing was committed or pushed, and nothing was sent to KBase.
- Every run went through `apply_scratch/figures/run.sh`, which sets:
  - `nice -n 10`;
  - thread variables at 24;
  - `MPLBACKEND=Agg` and `PYTHONDONTWRITEBYTECODE=1`;
  - `/usr/bin/time -v`;
  - the venv Python `~/Documents/py_venv/bin/python`.
- Per-step logs are in `apply_scratch/figures/logs/`. Pre-stage copies of every output this stage could touch are in
  `apply_scratch/figures/before/`, and the inspection crops are in `apply_scratch/figures/inspect/`.

**Renderer rule.** I used only the scripts the gate marked "use" (they reproduce HEAD, tier A) or "use (diff vs gate)"
(the script reproduces HEAD and only the committed file was stale).
- Each "diff vs gate" output was compared with the gate's own HEAD re-run in `gate_scratch/<step>/out/`. Every
  difference had to be an R232 change, or a D4 overwrite of a stale older-map ID.
- Every gate "do NOT use" or render-drift script was skipped and documented (section 5).

## 0. Pre-flight

- `apply_log/asv_md5_before.tsv`: 1,443 of 1,443 files match before any change (`logs/asv_md5_start.txt`).
- I snapshotted `git status --porcelain` (80 lines), `git diff --numstat`, the list of all 45 PDFs and the list of
  `__pycache__` directories. I recorded the md5 of the 25 figure outputs this stage could touch
  (`before/outputs_md5.tsv`) and copied each into `before/`.
- `free -g`: 161 GB available. Load average about 18 on 64 cores. No step needed more than 0.3 GB.

## 1. Code change (one, in the repo)

| file | change | numstat delta |
|---|---|---|
| `scripts/render_fig1_tree.py` | `root_tree()`: the figure subtitle now names the rooting clade by the phylum name its tips carry, `"/".join(sorted({Phylum of the outgroup tips}))`, instead of the literal "Patescibacteria". | +4 / −1 (file total vs HEAD now 47/13; the data stage left it at 43/12) |

The line is: `names = "/".join(sorted({r["Phylum"] for r in ann if r["Phylum"] in PATESCI})); how = f"{names} clade ({len(pates)} tips)"`.

- **Why.** With D3, all 7 outgroup tips are `Patescibacteriota` in R232 (4 MAGs and 3 references). The phylum legend
  said "Patescibacteriota (4 MAGs, 3 refs)", but the subtitle still said "rooted on Patescibacteria clade". The SI
  Fig S4 legend is being changed to Patescibacteriota (docs PLAN §4.3, para 252).
- **Effect on R220 inputs.** With R220 inputs the string is "Patescibacteria clade (N tips)", exactly as at HEAD.
- **Critic C8 gate, reinterpreted.** The gate exists to catch a silent fallback to midpoint rooting. It now reads
  "`<name> clade (7 tips)`, not `midpoint`". This run prints `root: Patescibacteriota clade (7 tips)`, with 116 of
  123 nodes marked, as at HEAD.
- `fig1_tree/compare_trees.py` was left as it is. It prints the literal "Patescibacteria clade" to stdout only, and
  its TSV holds bins and accessions.

## 2. Renders, in order

| step | script | outputs | result |
|---|---|---|---|
| C1a | `scripts/select_fig1_mags.py` | `fig1_tree/fig1_mags.txt`, `fig1_mags_table.tsv` | asserts pass (65 tips, the 3 expected new MAGs, 6 focal); **both byte-identical** to the data stage's bin-patched files |
| C1b | `scripts/render_fig1_tree.py --tree fig1_tree/fig1_mags.raxml.support --out fig1_tree/fig1_tree --annotation-out fig1_tree/fig1_annotation.tsv` | `fig1_annotation.tsv` (CRLF kept: 127/127 lines; numstat 53/53), `fig1_tree.png` (4125×4125, unchanged size), untracked `fig1_tree.pdf` | prints "references not in GTDB R232 …: ['GB_GCA_017302355.1']"; root `Patescibacteriota clade (7 tips)`; support ≥ 70: 116/123, as at HEAD |
| C1c | `fig1_tree/compare_trees.py --tsv fig1_tree/compare_trees_changed_splits.tsv` | TSV | **byte-identical**; both trees rooted on the 7-tip clade, RF = 2 as at HEAD |
| C3 | `clade_classify/render_nosz_top5_tree.py` (= SI Fig S5) | `out/nosz_top5_tree.png`, untracked `.pdf` | PDF token diff vs the gate HEAD re-run: exactly `Ignavibacteria.21 → Tepidaquicellales.1_m`. Canvas 5309 → 5311 px wide (longer label, tight bbox) |
| C5 | `scripts/render_percell_composition.py` (= main Fig 1) | `figures/percell_composition.png`, untracked `.pdf` | token diff vs the gate run: `SCMT01→Lysobacter`, `Azospira→Proximibacter`, `Giesbergeria→Simplicispira`, `(Patescibacteria)→(Patescibacteriota)`. Canvas unchanged. Per-phase totals 100 ×5. Band values: see 3.3 |
| C5 | `ppk1_classify/render_ppk1_tree_fig.py` | PNG and tracked PDF | PNG **byte-identical**. Tracked PDF: same size and text, only `CreationDate` differs, so **restored from HEAD** (`git checkout`) |
| C2 | `gene_ab_figure/09_gene_abundance_over_time.py` | TSV, PNG, tracked PDF | TSV and PNG **byte-identical**; tracked PDF metadata-only, so **restored from HEAD** |
| C2 | `gene_ab_figure/04_plot_heatmap.py` (diff vs gate) | `gene_abundance_figure.png` / tracked `.pdf` | vs the gate HEAD re-run, the differing pixels (0.058 %) all lie in the x-axis label band (bbox y 1579–1732). Visible changes: `Ignavibacteria.31 → Tepidaquicellales.1_m` (R232), `Giesbergeria.1_m → Simplicispira.1` (R232). Against the committed file it also fixes the D4 labels `Ca. Competibacter.139 → .85` and `Dokdonella.14 → .8`. Canvas unchanged |
| C2 | `gene_ab_figure/06_plot_copy_number.py` (diff vs gate) | `gene_copy_number_figure.png` / tracked `.pdf` | as for 04: 0.117 % of pixels, all in the label band (y 815–981); same label changes. Canvas unchanged |
| C2 | `gene_ab_figure/10_novel_denitrification_bubble.py` (diff vs gate) | `novel_denitrification_bubble.png`, untracked `.pdf` | **column set changes for an R232 reason** (3.4): `+Azonexus.12` (CAN_5_bin.70), `−Ca_Sarcinithrix.1` (coasm_bin.260). Still 7 columns; canvas unchanged; layout checked by eye |
| C2 | `gene_ab_figure/13_novel_metabolic_potential.py` (diff vs gate) | `novel_metabolic_potential.png`, untracked `.pdf` | **pixel-identical to the gate HEAD re-run.** Against the committed PNG only the D4 label `Dokdonella.14 → .8` differs (bbox 1372,3961–1579,4140). Its MAG set is fixed by `data/novel_functional_ko_long.tsv` (3.4) |
| C3 | `clade_classify/render_organism_tree.py` (diff vs gate) | `organism_tree/organism_tree.png`, untracked `.pdf` | canvas 3843×2479, equal to the gate re-run (the committed file was 3842×2495). R232 changes: `Giesbergeria · CAN_5_bin.40 [Giesbergeria.1_m] → Simplicispira · … [Simplicispira.1]`, `[Ignavibacteria.21] → [Tepidaquicellales.1_m]`. D4 fixes through the guild CSV: Ca_Accumulibacter.6→.2, Ca_Competibacter.49→.51, Ca_Leptovillus.11→.9, Ca_Sarcinithrix.2→.1. **All 20 labels parsed from the PDF: ID equals the R232 map, genus equals the R232 GTDB genus (0 mismatches)** |

Run times were 0.2–4.2 s each, and the peak RSS was 0.3 GB.

## 3. Verification detail

### 3.1 Fig S4 annotation (`logs/C1b_annotation_diff.txt`)

- **Shape.** Columns and tip order are identical to HEAD: 126 tips.
- **Rows changed.** 31 MAG rows and 22 reference rows. Every changed cell is explained:
  - **MAG tips (65).** `label` and `iterativeID` equal `mag_iterativeID_old_to_new.json` (R232) for 65/65. The 7
    ranks equal the repo R232 `gtdbtk.bac120.summary.tsv` for 65/65. `skani_best_ref`, `skani_ANI`,
    `skani_AF_query` and `ANI_class` equal the R232 `Table_S6_skani_novelty.tsv` for 65/65. 0 mismatches.
  - **Columns unchanged.** `per_cell_mean_pct`, `black_circle_focal`, `closest_ref_tip` and `membership`.
- **The 8 MAG tip relabels:**
  - coasm_bin.365: Ca_Competibacter.32 → Contendibacter.1_m
  - CAN_4_bin.208: Azospira.3 → Proximibacter.2_m
  - coasm_bin.40: Flavipsychrobacter.6_m → Taibaiella.5
  - coasm_bin.337: Ignavibacteria.7 → Tepidaquicellaceae.3_m
  - CAN_1_bin.2: Ignavibacteria.5 → Tepidaquicellaceae.1_m
  - CAN_3_bin.221: Ignavibacteria.21 → Tepidaquicellales.1_m
  - CAN_1_bin.195: Ca_Falkowbacteria.7 → Patescibacterium.1_m
  - CAN_1_bin.48: Paceibacteria.1_m → Minisyncoccia.1_m
- **Red-circle changes (10 MAGs).** They follow the R232 Table S6 rebuild (D3; buckets 96/70/48/77 → 121/76/28/66).
  They are not render drift.
  - coasm_bin.321 85–95 → ≥95 (86.39 → 99.70)
  - CAN_5_bin.70 ≥95 → 85–95 (98.67 → 89.29; its R220 best reference GCA_017302355.1 left GTDB)
  - CAN_5_bin.1 <85 → 85–95
  - coasm_bin.53 and coasm_bin.10: open ring → ≥95
  - coasm_bin.579 <85 → ≥95
  - coasm_bin.489 none → ≥95
  - CAN_1_bin.2 none → open ring (low_AF)
  - coasm_bin.521 85–95 → ≥95
  - CAN_1_bin.195 none → 85–95
- **References (61).** 60 lineages and labels equal the prep `gtdb_downloads/r232/bac120_taxonomy_r232.tsv.gz`
  (60/60), and all 60 are still species representatives in `sp_clusters_r232.tsv`.
  - **GB_GCA_017302355.1** is absent from R232. Its row equals HEAD apart from the label "GCA_017302355.1 (not in GTDB
    R232)", and it keeps its R220 lineage *Azonexus* sp017302355 (D3).
  - **9 label changes:** JARSSO01 sp018262135 → Contendibacter sp018262135; SCMT01 sp019510125 → Lysobacter_I
    sp019510125; Accumulibacter sp000585075 → Accumulibacter regalis; Accumulibacter phosphatis_G → Accumulibacter
    delftensis; Azospira_A restricta → Proximibacter restrictus; Azonexus sp017302355 → GCA_017302355.1 (not in GTDB
    R232); Niabella sp019316985 → Niabella defluvii; Flavipsychrobacter sp937871155 → JAFIGV01 sp937871155; UBA919
    sp002344425 → Patescibacterium sp002344425.
  - **3 phylum changes:** Patescibacteria → Patescibacteriota (GCA_002344425.1, GCA_021814805.1, GCA_002426885.1).
- **Figure text** (checked visually):
  - title: "bac120 ML tree: 65 MAGs + 61 GTDB r220 reference genomes (GTDB R232 names)";
  - subtitle: "… rooted on Patescibacteriota clade (7 tips)";
  - phylum legend: Pseudomonadota 32/29, Bacteroidota 22/22, Chloroflexota 4/4, Patescibacteriota 4/3, Acidobacteriota
    2/2, Verrucomicrobiota 1/1.

### 3.2 PDF token diffs against the gate HEAD re-runs

The files are `logs/pdf_token_diff_tierA.txt` and `logs/pdf_token_diff_vs_gate.txt`, as listed in section 2.
- The Fig S4 PDF text comes out fragmented, because its labels are rotated. That PDF was therefore verified through
  the annotation table (3.1) and by eye, not through the token diff.
- The diffs for 04, 06 and 10 are fragmented in the same way. 04 and 06 were checked by pixel bbox and by a
  side-by-side crop of the label band (`inspect/gene_ab_labels_committed_vs_R232.png`).
- The labels that overlap the caption in 04 are identical in the committed figure. That problem predates this stage.

### 3.3 Per-cell composition (main Fig 1): band values

The printed table equals the gate run except in these bands. Five-phase means, gate → now:
- **Renames:** SCMT01 → Lysobacter (1.078), Azospira → Proximibacter (0.661), Giesbergeria → Simplicispira (0.644).
  The values are identical; colours follow the same MAGs (data-stage DOCX_COLORS fix).
- **Unclassified 1.892 → 0.568.** 11 MAGs gain a GTDB genus in R232. It drops from 8th to 20th (the last band) in the band order.
- **JJ008 1.081 → 1.043.** coasm_bin.492 is CADEDV01 in R232 and is now pooled in Others.
- **Others 15.462 → 16.825.**
- The 20 bands are the same entities: 3 renamed, Unclassified reordered, no fallback colour.
- The OLB5 legend text "(Ignavibacteria (Bacteroidota))" is the same as at HEAD (a hard-coded parenthetical). The
  data-stage note that it changes was inaccurate, but harmless.

### 3.4 The novel-bacteria figures (R232 science, for D7 and docs)

- **Figure 10** defines "species-novel" as a blank GTDB species in `taxonomy_labels.tsv`. Under R232:
  - coasm_bin.260 (Ca_Sarcinithrix.1) gains *JAAEKA01* sp049294015, so it drops out;
  - CAN_5_bin.70 (Azonexus.12) loses *Azonexus* sp017302355, whose reference left GTDB, so it joins.
  - The printed novel list differs from the gate run by exactly those two lines.
- **Figure 13** takes its 7 MAGs from `data/novel_functional_ko_long.tsv` (`11_novel_functional_annotation.py`, not in
  the plan and not run). It still shows coasm_bin.260 and lacks CAN_5_bin.70.
- **So 10 and 13 now disagree on membership.** The subtitle of 13 ("species-novel MAGs (no GTDB species assignment)")
  is no longer true for coasm_bin.260.
- `data/novel_bacteria_summary.{tsv,md}` (bin-patched by the data stage; `12_…` fails at HEAD) still lists the R220
  set of 7.
- Hand-off: decide whether the "novel bacteria" module keeps its R220 set of 7 (with a note) or moves to the R232 set,
  which would need `11_` re-run. Either way the figure 13 subtitle and the summary text need a decision.

### 3.5 C22 / D1 PDF check (`logs/D1_pdf_check.tsv`)

All 45 PDFs were checked with `pdftotext -layout` and a boundary-aware match.
- **What was searched for.** The 30 `_m` old IDs (MAG-only), the 21 inherited old IDs (reported, but they are ASV
  labels in network PDFs), and R220-only GTDB names (Giesbergeria, SCMT01, JARSSO01, Azospira_A, UBA919,
  Flavipsychrobacter, sp000585075, phosphatis_G).
- **Regenerated in this stage, all clean.** No old `_m` ID, no old inherited ID and no R220-only name in:
  - `fig1_tree/fig1_tree.pdf`, `clade_classify/out/nosz_top5_tree.pdf`, `clade_classify/organism_tree/organism_tree.pdf`,
    `figures/percell_composition.pdf` (all untracked);
  - `gene_ab_figure/{gene_abundance_figure,gene_copy_number_figure}.pdf` (tracked);
  - `gene_ab_figure/novel_{denitrification_bubble,metabolic_potential}.pdf` (untracked).
- **`gene_abundance_over_time.pdf` and `ppk1_reference_tree.pdf`** were restored from HEAD and are clean. They show as
  "rewritten" only because `git checkout` reset their mtime.
- **SPIEC-EASI PDFs** (`spieceasi_mag_network{,_mb}`, `spieceasi_modules_{flanking_89,mag_display}`, rewritten by the
  data stage, and `spieceasi_mag_modules.pdf`, last written 2026-07-10). They hold 9 inherited old IDs (Azospira.3,
  Ca_Competibacter.32, …) as **ASV node labels**, which is correct under D2. No old `_m` IDs.
  `spieceasi_root_{modules,network}.pdf` are clean.
- **Still stale, untracked, NOT for insertion into a DOCX.** Each comes from a skipped renderer or is historical:
  - `denitrification_MAG_heatmap.pdf`: `Giesbergeria.1_m`, `Gammaproteobacteria.1_m`, `Azospira.3`,
    `Ca_Competibacter.32`, `Ignavibacteria.21/.5`, `Giesbergeria`;
  - `clade_classify/out/nosz_bb_reference_tree.pdf`: Giesbergeria, sp000585075;
  - `fig1_tree/superseded_ref030917225/fig1_tree.pdf` and `fig1_tree/test62_fig1_tree.pdf`: superseded test renders;
  - `paper3/papert3_071626.pdf`: Emily's master (docs stage);
  - `Experimental Data Inconsistencies.pdf`: historical.
  - `network/N2O_sink_producer.pdf` is a 276-era render (gate) but contains no changed ID. `network/chord_*.pdf` are
    clean.

### 3.6 C6 under D2 (frozen membership)

- **No network figure changed in this stage**, and none needed to:
  - The data stage re-ran C6b, e, h, i, j and k and proved their PNGs byte-identical to HEAD (D2 acceptance).
  - The aggregate (C6a) and per-phase (C6c) grids, and the GAO/PAO figures in their committed default mode (C6d),
    never read the MAG map or the frozen list. They call `rcg.load_abundances()` without `NETWORK_MAG_ONLY`.
  - So I did not re-run them. Their outputs are in `asv_md5_before.tsv` and still match (4.).
- Node labels remain ASV IDs. The 9 inherited old IDs visible in the SPIEC-EASI figures are those ASV labels.
- **C6f/C6g module panels** (`module_pathway_completeness/*`): skipped (gate: do NOT use, pre-BH-FDR content). Their
  column labels are the ASV node IDs of the frozen MAG-linked nodes, so under D2 they are not R232-stale.
- **C6l `render_spieceasi_docx.py`**: it writes a document, so it belongs to the docs stage. Not run.

## 4. End checks

| check | result |
|---|---|
| ASV-keyed md5 (`asv_md5_before.tsv`) | **1,443 / 1,443 unchanged** (`logs/asv_md5_end.txt`) |
| `community_modeling/plan_r13/integrator/AMENDMENTS.sha256` (run from `community_modeling/`) | **14/14 OK** |
| r13 `track2_guards` (`track2_guards_sha256_before.txt`) | **7/7 OK** |
| EOL | `fig1_annotation.tsv` CRLF 127/127 (as at HEAD); `fig1_mags_table.tsv` (LF), `compare_trees_changed_splits.tsv` (LF), `gene_abundance_over_time.tsv` (CRLF): byte-identical |
| `git diff --numstat` of `fig1_annotation.tsv` | 53/53 (added == deleted) |
| `__pycache__` | none created (same 47 directories as before) |
| D5 | no log, run record, `*.pre_*.bak.*`, payload/verify JSON, r12 snapshot or track2_guards file touched (`git diff --numstat -- '*.bak*' '*.log'`: 0) |
| `git status` | 80 → 92 lines: +10 modified by this stage, +2 untracked **not from this stage** (see below) |

**Files changed by this stage**
- **Tracked (12):** 10 newly modified; `render_fig1_tree.py` and `fig1_annotation.tsv` were already modified by the
  data stage.
  - `scripts/render_fig1_tree.py`;
  - `fig1_tree/fig1_annotation.tsv`, `fig1_tree/fig1_tree.png`;
  - `clade_classify/out/nosz_top5_tree.png`, `clade_classify/organism_tree/organism_tree.png`;
  - `figures/percell_composition.png`;
  - `gene_ab_figure/gene_abundance_figure.{png,pdf}`, `gene_ab_figure/gene_copy_number_figure.{png,pdf}`;
  - `gene_ab_figure/novel_denitrification_bubble.png`, `gene_ab_figure/novel_metabolic_potential.png`.
- **Untracked, rewritten (C22):** `fig1_tree/fig1_tree.pdf`, `clade_classify/out/nosz_top5_tree.pdf`,
  `clade_classify/organism_tree/organism_tree.pdf`, `figures/percell_composition.pdf`,
  `gene_ab_figure/novel_denitrification_bubble.pdf`, `gene_ab_figure/novel_metabolic_potential.pdf`.
- **Re-run, byte-identical or restored:** `fig1_mags.txt`, `fig1_mags_table.tsv`, `compare_trees_changed_splits.tsv`,
  `gene_abundance_over_time.{tsv,png,pdf}`, `ppk1_reference_tree.{png,pdf}`.
- The md5 of each, before and after, is in `apply_scratch/figures/outputs_md5_after.tsv`.

**Not from this stage.** `kbase_r232/` (`rename_mags_kbase_r232_log.tsv`) and `scripts/rename_mags_kbase_r232.py`
appeared in the repo during this run. They are written by the concurrent KBase stage. I did not touch them.

## 5. Skipped, with the reason (gate verdicts) — the manual/decision list

| renderer | outputs left as they are | gate verdict | what stays stale (visible) |
|---|---|---|---|
| `gene_ab_figure/05_bundle_offline.py` | `offline_bundle/`, `.zip` | do NOT use (older snapshot; C9) | none. `offline_bundle/data/taxonomy_labels.tsv` was bin-patched by the data stage |
| `gene_ab_figure/12_novel_bacteria_summary.py` | `data/novel_bacteria_summary.{tsv,md}` | fails at HEAD (S1 layout drift) | the summary text keeps the R220 novelty of coasm_bin.260 (3.4; docs) |
| `clade_classify/render_nosz_tree_fig.py` | `out/nosz_bb_reference_tree.png` (+ untracked pdf) | do NOT use (content: tree inputs changed in ffe668c) | tip labels "GTDB genus/species · bin": **Giesbergeria** (CAN_5_bin.40, now Simplicispira) and *Accumulibacter* sp000585075 (now *A. regalis*). These are the two R220-only names that the pdftotext check (3.5) found in its PDF. A re-render would also fold in the ffe668c tree change. Decision for Andrew |
| `phase_rpkm_chord/01_build_long_table.py` | category CSVs | do NOT use (values: local `counts_CAN_*.json` differ from poplar's) | — (`mag_category_assignments.csv` / `phase_mag_gene_rpkm_long.csv` were bin-patched by the data stage) |
| `phase_rpkm_chord/03_plot_chords.py`, `04_plot_heatmap.py` | 10 chord PNGs, 2 heatmap PNGs | render drift (0.5–2.4 % of pixels, ±1 px) | rendered from the HEAD CSV (`logs/skipped_outputs_stale_ids.txt`). **R232:** CAN_3_bin.221 `OLB5_m.1` (now Tepidaquicellales.1_m); CAN_5_bin.40 `Giesbergeria_m.1` / genus Giesbergeria (now Simplicispira.1 / Simplicispira); CAN_5_bin.177 genus SCMT01 (now Lysobacter_I). **Pre-R232 (D4):** Dokdonella.14, UBA5066_m.1, Pseudoxanthomonas_A_m.1, Denitratisoma.2, JAEULV01_m.1, Ca_Competibacter.139, JAAEKA01_m.2, Leptovillus_m.2. Re-rendering 03/04 from the patched CSVs would fix the labels at the cost of the drift. Decision for Andrew |
| `scripts/render_chord_rpkm_function_guild.py` | `network/chord_*.png` | styling drift | none (guilds and genes only; no IDs or GTDB names). No R232 action |
| `scripts/render_n2o_sink_producer.py` | `network/N2O_sink_producer.png` (+ 276-era untracked pdf) | render drift (canvas +2/+4 px) | **none.** Its only labels are the `mag_iterativeID` of MAGs whose GTDB genus matches Accumulibacter/Competibacter/Thauera. No such MAG changes ID, and no genus changes into or out of that match (coasm_bin.365 is JARSSO01 → Contendibacter, neither matched). Panel A counts come from KO data. So R232 changes nothing in this figure |
| `scripts/render_denitrification_heatmaps.py` (docs step 11; closest to main Fig 4) | `denitrification_MAG_heatmap.png`, `denitrification_group_abundance.png` (+ untracked pdfs) | do NOT use (reads S1 N–R as relative abundance, but they now hold mapped reads; output 28,523 px wide) | x-axis IDs in the 6a04588 map, including Giesbergeria.1_m, Gammaproteobacteria.1_m, Azospira.3, Ca_Competibacter.32, Ignavibacteria.21/.5. Emily re-renders her Fig 4; or decide on the code fix (% = 100·mapped/Σmapped) |
| `scripts/render_module_denitrification_genes.py`, `render_module_pathway_completeness.py` (C6f/C6g) | `module_pathway_completeness/*` | do NOT use (pre-BH-FDR content) | none for R232 under D2 (labels are ASV node IDs). Legacy DOCX images 6/7, not in Emily's Oct-21 master |
| `scripts/render_gao_pao_networks.py`, `render_cooccurrence_grid{,_per_phase}.py` (C6a/c/d) | network PNGs/JSONs | tier A, but independent of the map | — (not re-run, 3.6) |
| `network/compositional_networks/render_spieceasi_docx.py` (C6l) | `SPIEC-EASI_network_analysis.docx` | — | docs stage |

## 6. Hand-offs

- **Docs.**
  - Fig S4 legend: "rooted on the Patescibacteriota clade". The figure now says so. The title reads "61 GTDB r220
    reference genomes (GTDB R232 names)", which matches the "relabel only" option of docs §8 item 2.
  - Main Fig 1 legend: Lysobacter, Proximibacter, Simplicispira, UBA5532 (Patescibacteriota). JJ008 is now
    CAN_5_bin.112 only. "Unclassified" falls to 0.57 % (five-phase mean).
  - Fig S5: Tepidaquicellales.1_m.
  - Replace the embedded images from the regenerated PNGs; do not use any PDF in the stale list of 3.5.
  - The novel-bacteria module: 3.4.
  - Living notes next to the regenerated figures still quote the old subtitle, IDs or names ("Patescibacteria clade",
    Ignavibacteria.21, Giesbergeria): `fig1_tree/README.md`, `clade_classify/NOSZ_TOP5_NOTES.md`,
    `clade_classify/CLADE_I_II_DISCREPANCY.md` and `gene_ab_figure/METHODS_nosz_clades.md`. They are text edits for
    the docs stage; this stage did not edit them.
  - **Rollback path for the novel-bacteria module.** If the module is to stay at its R220 set until Andrew decides,
    restore `gene_ab_figure/novel_denitrification_bubble.png` from git and the untracked `.pdf` from
    `apply_scratch/figures/before/gene_ab_figure/`.
- **D7 / change list for Emily** (scientific, not renames):
  - The Fig S4 red-circle changes listed in 3.1 (10 MAGs; for example CAN_5_bin.70 ≥95 → 85–95 because its R220
    reference left GTDB).
  - The per-cell Unclassified band shrinks.
  - The novel set of figure 10 is now `+CAN_5_bin.70 −coasm_bin.260`.
- **Andrew's decisions:**
  - the C4 chord re-render (accept drift to get R232 labels?);
  - the nosZ reference tree re-render (folds in the ffe668c tree change);
  - the denitrification-heatmap code fix;
  - whether to re-run `11_novel_functional_annotation.py` for figure 13.
