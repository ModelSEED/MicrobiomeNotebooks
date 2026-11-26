# Gate stage: every builder and renderer re-run at HEAD 9e087c3 (critic C15 / correction 17)

Run 2026-11-11, 15:40–16:30 CST.
- **Nothing in `~/Documents/EmilyKin` was changed.** Every run happened in a detached worktree at
  `rename_plan/head_worktree` (9e087c3), which is removed at the end of this stage.
- Main-repo checks at the end:
  - `find -newer <stamp>` (outside `.git`) returns nothing;
  - `git status --porcelain` is unchanged (empty);
  - all 13,425 tracked files match `tracked_md5_head.txt`.
- Nothing was committed, and nothing was sent to KBase.

**Method.**
- Each script was run alone from a clean HEAD worktree, not chained, so an upstream script that does not reproduce
  cannot affect the verdict on the scripts downstream of it.
- Each run used `nice -n 10`, `/usr/bin/time -v`, thread variables set to 24, `MPLBACKEND=Agg`, and `strace` to find
  untracked inputs.
- Every file a run wrote was compared with its HEAD blob. A gitignored PDF had no HEAD blob, so it was compared with
  the main repo's untracked copy.
- After each run the worktree was restored (`git checkout`, then deleting the run's new files).

**Comparison tiers.**

| tier | meaning |
|---|---|
| A | byte-identical |
| B | equivalent: pixel-identical PNG, text-identical PDF, cell-identical XLSX, text/table/media-identical DOCX, parsed-equal JSON, or numerically identical CSV/TSV (≤1e-9) |
| C | differs (each case is described below) |
| N | new file with no committed counterpart |

- **Software:** Python 3.12.3, pandas 3.0.3, numpy 2.4.6, matplotlib 3.10.3, openpyxl 3.1.5, networkx 3.6.1,
  scipy 1.17.1, statsmodels 0.14.6, pycirclize 1.10.1, Biopython 1.85, Pillow 12.1.1, python-docx 1.2.0.
- **Runtime:** 54 harness runs (50 distinct steps plus 4 determinism repeats), 20.4 min wall in total, plus the
  16 s xlsx round trip. The longest were the per-phase networks (7.0 min, 3.1 GB), the aggregate network (5.0 min,
  1.6 GB) and `build_panel.py` (1.5 min, 2.3 GB).
- **Machine-readable results:**
  - `apply_log/gate_results.tsv`: one row per written file, with tier, comparison text, runtime and RSS.
  - `gate_scratch/<step>/`: `run.log`, `time.txt`, `strace.txt`, `result.json`, and in `out/` every output that was
    not byte-identical.

## 1. Baselines recorded before any change (main repo = HEAD)

| file (in `apply_log/`) | content |
|---|---|
| `asv_md5_before.tsv` | md5, bytes and EOL for 1,443 files. 1,215 are the `files_plan.tsv` rows whose action starts with `leave`. 227 are core ASV files with no old-ID hit, which `files_plan.tsv` does not list (abundances, `taxonomy.json`, `iterativeID*.json`, `relative_abundance.json`, `asv_to_mag_mapping.json`, `asvs.fasta`, the colour maps, `correlations/`, `dbRDA/`, `asv_mag_mapping/`, `network_method_comparison/`, the aggregate and per-phase network JSONs, `cooccurrence_edges_p_value_FDR.csv`, …). Column `asv_keyed`: 1,419 yes, 23 "leave, not ASV-keyed", and 1 flagged (`dbRDA/dbRDA_performance_iterativeID_scores.md`, whose plan action is a text check). The two all-phases network JSONs carry the C13 note: MAG-only network, results kept by D2-frozen, fields may be added |
| `tracked_md5_head.txt` | md5 of all 13,425 tracked files at HEAD. This is for the final "nothing else changed" check |
| `eol_before.tsv` | line terminator, line count and size of every `apply_mag_ids_extension.tsv` target (C11) |
| `track2_guards_sha256_before.txt` | sha256 of r13 `track2_guards` `AMENDMENT_1..4.md`, `CRITERIA.md`, `run_plan.json` and `acceptance.json`. `CRITERIA.md` = df0e8dc3…, which equals the value registered in `run_plan.json` |

- `community_modeling/plan_r13/integrator/AMENDMENTS.sha256`: **14/14 OK**.
  - It must be checked from `community_modeling/`, because its paths are relative to that directory. Run from the repo
    root, it reports 14 files as missing.
- **C11 correction.** The ADD targets are not all LF.
  - LF: Table S7, `mag_guild_classification.csv`, `mag_percell_abundance.csv`, `fig1_mags_table.tsv`,
    `rpkm_long_MAG_gene_sample_phase.tsv`, and r13 `taxonomy.tsv`, `mag_panel.tsv`, `_manifest.tsv`.
  - CRLF: `fig1_annotation.tsv`, `novel_bacteria_summary.tsv`, `rpkm_long_…_291_preview.tsv`, `ko_rpkm_combined_ext.tsv`,
    `ko_rpkm_per_sample_ext.tsv`, and both `phase_rpkm_chord/*.csv`.
  - The 9 targets already in `apply_mag_ids.py` are all CRLF, and `mag_abundance_summary.json` is LF.
  - The extension must keep each file's own terminator.

## 2. Verdict per script

**Key**
- "use": the script reproduces HEAD and may regenerate the committed output after the R232 inputs land.
- "use (diff vs gate)": its output differs from the committed file only because the committed file is stale. Judge the
  R232 effect by diffing the R232 output against the gate's own HEAD re-run in `gate_scratch/<step>/out/` (or against
  the committed file when tier A), not by assuming the committed file was current.
- "do NOT use": the script does not reproduce HEAD, so the fallback applies.

### Data part B1–B9

**B1 `scripts/recompute_mag_ids.py`.** 0.2 s.
- Tier A for 3 of 3 outputs (`mag_iterativeID_old_to_new.{csv,json}`, `scripts/_new_ids.json`). Cases 155/62/74,
  0 duplicates.
- Use as planned, with the gate against the simulation map.

**B2 `scripts/apply_mag_ids.py`.** 0.2 s.
- Tier A for 10 of 10. It is idempotent at HEAD, and its CRLF output equals the CRLF files.
- Use it for the 9 listed files. The ADD extension must preserve EOL (C11). Note `ko_rpkm_per_sample.tsv`: 42 rows are
  not in the map, because they are 276-era MAGs.

**B4 `ko_gene_abundance/04_root_aggregate.py`.** 0.2 s. **Tier C. Do NOT use.**
- What differs: 150 → 146 roots in `ko_rpkm_root_combined.tsv`; 2,907 → 2,835 rows in `_per_sample.tsv`.
- Why: the committed root tables were last built at 924fa18 (2026-06-27), with that era's IDs. Their input,
  `ko_rpkm_combined.tsv`, carries the HEAD IDs. Totals are identical (36,027.0096 RPKM, 276 MAGs). Roots only in the
  committed file: Chthoniobacteraceae, Flavitalea, Mesorhizobium_I, Pseudorhodobacter, Sandaracinus. Only in the
  re-run: Polyangiales.
- No bin-keyed patch is possible, because the rows are keyed by root. Leave the tables as they are and list them as
  stale. Or, with the orchestrator's explicit approval, re-run and separate the two effects by diffing the R232 output
  against `gate_scratch/B4_root_aggregate/out/`.

**B5 `scripts/build_percell_abundance.py`.** 2.1 s.
- Tier A. Workbook regression PASSED.
- Use it, or apply the ADD bin join.

**B6 `reprocess_291/build_xlsx_291.py`** (`EMILYKIN_REPO=<worktree>`). 1.5 s. **Tier C. Do NOT use.**
- S1: 1,779 cells.
  - The builder writes relative abundance into columns N–R and T. At HEAD those columns hold `mapped_reads in CAN_i`
    and `average mapped reads`.
  - It also reads `ko_gene_abundance/counts_291_CAN_*.json`, which has no `per_mag_mapped`, so it wrote 0 everywhere.
- ID cells rewritten from the 6a04588 map to HEAD: S2a 25, S3 33, S4 17, S5 2.
- S6: 1,585 cells cleared to the skeleton, by design.
- Fallback: a bin-keyed cell patch of `MAG_target_gene_filled.xlsx`, using the docs part's `table_cells_plan.tsv`
  rows for that file.

**B7 `scripts/build_table_s7_classification.py`.** 2.0 s.
- Table S7 CSV and the guild CSV: tier A. Builder workbook S7 sheet: tier B, cell-identical.
- Use it, after S1/S3 of the builder workbook are patched by bin. `iterative_ID` comes from S1 column B, which holds
  6a04588-era IDs, so B7 alone keeps the stale IDs. Apply D6 (the `Contendibacter` regex) first.

**B8 / docs step 2 `reprocess_291/build_skani_novelty_291.py`.** **Cannot be re-run at HEAD.**
- Its R220 291-MAG raw search (`/scratch1/…/skani_novelty_291/mag_vs_gtdb_skani.tsv`) and `gtdb_taxonomy_r220.tsv`
  exist only on poplar.
- Logic check instead, with the 276-era raw `server_handoff_out/mag_vs_gtdb_skani.tsv` and the R220 rep taxonomy
  built from `~/Documents/emilykin_tools/gtdbtk_r220_markers/sp_clusters_r220.tsv`:
  - 214 of the 238 MAGs that had a 276-era query get all 7 skani columns identical to the committed 291 Table S6.
  - The other 24 also differ from the committed 276 handoff table. That is consistent with their genome changing in
    the best-of-both catalogue.
  - `ANI_class` and the sentence re-derive from the committed ANI/AF for 291 of 291.
  - The committed `skani_best_taxonomy` equals the R220 rep taxonomy for 251 of 251.
  - `mag_gtdb_taxonomy` equals the repo summary for 290 of 291. The exception is coasm_bin.221, whose genus is blank
    in S6.
- Use it for the D3 rebuild against the new local R232 search. Patch `SK`, `TAXR`, `GTDBTK`, `MAGDIR` and `OUT`.
  `TAXR` is read with `open()`, so give it a decompressed `bac120_taxonomy_r232.tsv`, or patch the script to use
  `gzip.open`. The output is a new computation (D3), not a reproduction.

**Docs step 3 `scripts/build_table_s6_ani.py`** (no args; auto-loads the repo gtdbtk summary). 1.9 s. Tier C, ID cells
only.
- 33 cells in S6 column B differ: 6a04588-era IDs become HEAD IDs. That is the D4 drift, and overwriting it is wanted.
- Every other S6 cell reproduces: closest reference, ANI, AF, call, interpretation, notes and the "Reference database"
  string. 291 rows were filled from GTDB-Tk and 0 from skani, so the HEAD sheet was built without `--skani`.
- Use it, after the C7 patch (`REF_DB`, the line-307 note, the line-173 help text).
- **Method parity.** The HEAD sheet used the GTDB-Tk `closest_genome_*` columns only: 0 MAGs were resolved via skani
  and 194 remain "unresolved". Docs step 2 says to pass `--skani <R232 search>`. Doing so changes the method for those
  194 rows, which would be filled from skani. That is a method change, not a rename.
- For parity, run `--gtdbtk <R232 summary>` without `--skani`. Use `--skani` only if the orchestrator decides so, and
  say so in the change log.
- The skani search (D3) feeds the standalone `Table_S6_skani_novelty.tsv` through `build_skani_novelty_291.py`, as at
  HEAD.
- It writes `.pre_s6.bak.xlsx` only if no backup exists. The main repo already has a 276-era one (gitignored), so the
  apply must take its own backup.

**Docs step 7 `reprocess_291/compute_manuscript_stats_291.py`** (patch `R="."`). 0.1 s. Stdout only.
- Every printed number equals `MANUSCRIPT_methods_results_discussion.md`: 14 phyla, Patescibacteria 21,
  Myxococcota 3, 290 / 273 (172 genera) / 97 (194 unnamed), methods 145/97/49, skani 96/70/48/77, tiers
  1116/286/694/1703, 3,799 / 151, 244 MAGs, 353 16S genes in 247 MAGs, nosZ 37/88/1.
- Use it for the R232 number check, after two fixes:
  - It prints "naming case: (no case col)", because it looks for `mag_case` while the CSV column is
    `assignment_case`. The 155/62/74 was confirmed separately from the CSV. Fix the column name if the script is used
    to confirm 138/70/83.
  - The Bakta glob points at `/scratch1` and prints `{}` here.

**B9a r13 `build_taxonomy.py`** (patch `EK`; the 286 untracked `rast/genomes/*.tsv` copied in). 8.1 s. **Tier C. Do NOT
re-run blindly.**
- Only 3 columns differ, each for 286 bins: `kbase_rast_genome_ref`, `rast_genomeset` and `envelope_evidence_rast`.
- They are filled now because the code-11 RAST sets landed after `taxonomy.tsv` was built. ANNOTATIONS.md and the
  FINALIZE log record this as an open operator item.
- The other 33 columns are identical, including `iterativeID`, taxonomy, `envelope_lineage`, `gram_flag` and
  `lps_expected`.
- **Recommended fallback (hybrid).** Run the fixed `build_taxonomy.py`, so the script itself applies the
  FLAG/DIDERM logic to the R232 phyla. Then restore `kbase_rast_genome_ref`, `rast_genomeset` and
  `envelope_evidence_rast` to their HEAD values, by bin. This keeps exactly what the script computes and drops only the
  RAST fill. Gate: after the restore, those 3 columns equal HEAD for 291 of 291.
- **Other options:**
  - A hand patch of `iterativeID` and the taxonomy columns by bin. Not preferred: it would re-implement the
    envelope/flag logic by hand.
  - The orchestrator explicitly accepts the RAST fill as part of B9.

**B9b r13 `panel/build_panel.py`** (patch `EK`). 89.6 s, 2.3 GB.
- Tier A: `mag_panel.tsv`, `mag_panel_loci.tsv`, `member_panel.tsv` and `member_panel_long.tsv` are byte-identical.
  `build_report.json` is equal after normalising the absolute input paths.
- Use it.
- **Note:** `panel/PANEL_FROZEN.md` holds the freeze-2 sha256 of `mag_panel.tsv` and `build_report.json`, and the input
  hashes of `gtdbtk.bac120.summary.tsv` and `kbase_mag_rename_log.json`.
  - Only that .md cites them; no code checks them.
  - B9 changes these frozen files. Add a dated freeze note or addendum; do not edit the freeze-2 hashes (D5).

**B9c r13 `panel/validate_panel.py`** (patch `EK`). 1.2 s.
- 18 tests: 16 PASS, 2 WARN, 0 FAIL. This is the FINALIZE record.
- Use it as the B9 post-check.

**B9d `pytest panel/tests/`.** 2.1 s.
- 15 of 15 pass.
- Use it as the post-check after the `p__Patescibacteriota` fixture edit.

**B9e r13 `calls/call_genes.py --reflag`** (patch `REPO`). 0.1 s.
- `calls/summary.tsv` is tier A.
- `--reflag` recomputes only `flags`. It does not refresh `gtdb_phylum` / `gtdb_order`, and there is no
  "summary-only" mode.
- A full run re-calls 291 MAGs with pyrodigal and rewrites the frozen, untracked `genes.tsv` / `proteins.faa`.
- So patch `gtdb_phylum` / `gtdb_order` by bin.

### Data part C1–C6 and docs figure steps

**C1a `scripts/select_fig1_mags.py`.** 1.0 s.
- Tier A for 2 of 2.
- Use it.

**C1b `scripts/render_fig1_tree.py --tree … --out fig1_tree/fig1_tree --annotation-out …`.** 11.8 s.
- Tier A: annotation TSV and PNG. Tier B: the untracked PDF, text-identical.
- It prints "root: Patescibacteria clade (7 tips)", which is the C8 gate string.
- Use it, after the C8 code change.

**C1c `fig1_tree/compare_trees.py --tsv …`.** 0.4 s.
- Tier A. Both trees are rooted on "Patescibacteria clade (7 tips)".
- Use it.

**C2 `08_gene_copy_per_mag.py`.** 0.1 s.
- `.tsv` tier A. `.md` tier C on 4 lines: the committed `.md` still has the 6a04588 IDs (Dokdonella.14,
  Ignavibacteria.31, …). The `.tsv` was updated by `apply_mag_ids` and the `.md` was not.
- Use (diff vs gate).

**C2 `04_plot_heatmap.py` / `06_plot_copy_number.py` / `10_…bubble.py` / `13_…potential.py`.** 37 s / 4 s / 9 s / 8 s.
- Tier C, labels only. The committed figures carry pre-408cb39 labels: Ca. Competibacter.139 → .85,
  Ignavibacteria.31 → .21, Dokdonella.14 → .8.
- The differing pixels are confined to those x-axis labels (cropped and checked by eye). The PDF text diff also
  consists only of those tokens.
- Use (diff vs gate).

**C2 `09_gene_abundance_over_time.py`.** 7.2 s.
- TSV and PNG tier A. PDF tier B.
- Use it.

**C2 `12_novel_bacteria_summary.py`.** **Fails at HEAD** (rc 1).
- It raises `ValueError: Unknown format code 'f' … str` at line 101. It reads S1 of the builder workbook by fixed
  index (`s[7]`), and the S1 layout has drifted.
- Fallback: patch `novel_bacteria_summary.{tsv,md}` by bin. No MAG among its 7 changes ID. coasm_bin.260's species
  gain is a text edit.

**C2 `05_bundle_offline.py`** (patch `WORK`). **Tier C. Do NOT use** (drop it, per C9).
- The committed `offline_bundle/` is an older snapshot:
  - the `gene_rpkm*` rows say `norC/nosZ` where the current data say `norBC/nosZ`, with values;
  - `target_genes.bed` has 42 rows, not 73;
  - `04_plot_heatmap.py` is an older copy.
- Fallback: the existing `apply_mag_ids` bin join of `offline_bundle/data/taxonomy_labels.tsv`, and nothing else.

**C3 `render_nosz_top5_tree.py`.** 10.4 s.
- PNG tier A. PDF tier B.
- Use it.

**C3 `render_nosz_tree_fig.py`.** 7.8 s. **Tier C, content. Do NOT use.**
- The PNG goes from 6568×3088 to 6569×2776. Tips change (CAN_5_bin.70, coasm_bin.260 and coasm_bin.481 are absent in
  the re-run) and supports change (e.g. [64] → [21], [100] → [98]).
- The committed figure (9c35150, 2026-07-06) predates the ffe668c nosZ tree inputs (2026-07-08).
- Its tips are labelled "GTDB genus/species · bin", with no iterativeIDs. Giesbergeria will change under R232.
- Fallback: leave it, and put it on the manual/decision list. A re-render also folds in the ffe668c tree change.

**C3 `render_organism_tree.py`.** 6.2 s. Tier C, labels only.
- 9 ID labels differ. The size goes from 3842×2495 to 3843×2479.
- The script takes IDs from `mag_guild_classification.csv`, which holds 6a04588-era IDs. The committed PNG shows an
  even older set.
- Use (diff vs gate), after the guild CSV is bin-patched.

**C4 `phase_rpkm_chord/01_build_long_table.py`** (patch `REPO`). 0.1 s. **Tier C, values. Do NOT use.**
- `MAG_gene_RPKM` differs by a uniform per-sample factor (CAN_1–5: 0.858, 0.860, 0.858, 0.904, 0.906). So the local,
  untracked `gene_ab_figure/data/counts_CAN_*.json` `total_mapped` values are not the ones poplar used.
- Guild and N-role calls also differ (e.g. CAN_1_bin.210: "Strong putative PAO-like" → "Other").
- The IDs are stale, including the pre-c796fdf `Root_m.N` form (C19).
- Fallback: a bin-keyed patch of `MAG_iteration_ID` in `mag_category_assignments.csv` and
  `phase_mag_gene_rpkm_long.csv`. Both are CRLF; keep that. Leave the category CSVs, which have no IDs.

**C4 `03_plot_chords.py` / `04_plot_heatmap.py`** (render the committed CSVs). 15 s / 8 s. Tier C, render drift.
- Labels and legend are the same (checked side by side). 0.5–2.4 % of pixels differ in ribbons and cells, and widths
  move by ±1 px.
- These are 2026-06-30 renders and their PDFs were never committed.
- Content-faithful but not pixel-reproducible. Use them only if the orchestrator accepts render drift; otherwise
  leave the figures, which already carry pre-c796fdf labels.

**C5 `render_chord_rpkm_function_guild.py`.** 23 s. Tier C, styling.
- Committed PNG versus re-run: label size and circle radius differ. The PDFs are text-identical to the main-repo
  untracked PDFs. The re-run is deterministic (two runs are byte-equal).
- No R232 action is needed: the figures show genres and genes, with no IDs or GTDB names.

**C5 `render_n2o_sink_producer.py`.** 4.1 s. Tier C, render drift.
- Visually identical: labels, counts, n = 291. The canvas goes from 4050×1699 to 4052×1703. The re-run is
  deterministic (seeded).
- The main-repo untracked PDF is a stale 276-era render (C22).
- Use it, with a visual check; its labels are MAG IDs.

**C5 `render_denitrification_heatmaps.py`.** 19.5 s, 1.1 GB. **Tier C. Do NOT use.**
- It reads S1 columns N–R (`r[13:18]`) as relative abundance, but at HEAD they hold absolute mapped reads.
- Result: the colour scale reaches 1e6 %, and the MAG heatmap grows from 6,324 to 28,523 px wide.
- Fallback: leave it (Emily re-renders main Fig 4), or make a code fix (% = 100·mapped / Σmapped, as
  `build_table_s7` does) by an explicit decision.

**C5 `render_percell_composition.py`.** 3.6 s.
- PNG tier A. PDF tier B.
- Use it, after the DOCX_COLORS fix.

**C5 `ppk1_classify/render_ppk1_tree_fig.py`.** 10.9 s.
- PNG tier A. PDF tier B.
- Use it.

**C6a `render_cooccurrence_grid.py`** (aggregate, all ASVs). 301 s, 1.6 GB.
- Tier A for 5 of 5.
- Use it, for labels/fields only (the D2 acceptance test).

**C6b `render_cooccurrence_grid_all_phases.py`** (MAG-only, `NETWORK_MAG_ONLY` defaults to 1). 34 s.
- Tier A for 4 of 4 (140 nodes).
- Use it. **D2 acceptance test:** with the frozen list, the re-run must be byte-identical apart from any added fields.

**C6c `render_cooccurrence_grid_per_phase.py`.** 419 s, 3.1 GB.
- Tier A for 20 of 20.
- Use it, as for C6a.

**C6d `render_gao_pao_networks.py`, default mode (no environment).** 65 s.
- **Tier A for 10 of 10.**
- In this mode the script never reads the MAG map, so it is independent of R232 and D2. Do not re-run it.

**C6d, the same script with `NETWORK_MAG_ONLY=1`.** 38 s.
- Tier C: 9 PNGs differ.
- C14 is resolved: the committed figures are the all-ASV mode. The plan's "run with `NETWORK_MAG_ONLY=1`" is wrong.

**C6e `build_module_mag_membership.py`.** 0.1 s. Tier C, provenance text only.
- 6 `gamma_rule` description strings differ, copied from network summaries whose wording was revised later. All
  module and MAG data are identical.
- Use it. A re-run also refreshes those 6 strings.

**C6f `render_module_denitrification_genes.py`.** 11 s. **Tier C, content. Do NOT use.**
- Module member sets differ.
- The committed outputs (408cb39, 2026-07-08) predate the BH-FDR per-phase networks (8e1f66b, 2026-10-16).
- Fallback: leave these legacy DOCX images 6/7, which are not in Emily's Oct-21 master. Put them on the decision list.

**C6g `render_module_pathway_completeness.py`.** 9 s. **Tier C, content. Do NOT use.**
- Same cause as C6f: the completeness values differ.
- Fallback as C6f.

**C6h `compositional_networks/render_spieceasi_mag_network.py`.** 8 s.
- Tier A for 3 of 3. The 2 PDFs are tier B.
- Use it (D2 acceptance test).

**C6i `compositional_networks/spieceasi_modules_flanking.py`.** 3 s. Tier C, order only.
- Hub and flanking-member order varies from run to run (3 runs): members tied on degree come from a set, so the order
  depends on PYTHONHASHSEED.
- As sets the content is identical. The 5 edgelist and membership CSVs are tier A.
- Use it, but compare list fields as sets. Under D2-frozen nothing should change.

**C6j root `spieceasi_modules_flanking.py` → `tablesforpaper/`.** 3 s.
- Tier A for 9 of 10. `module_by_carbon_source.csv` is tier B (1 cell differs in float representation only).
- Use it.

**C6k `render_module_flanking_figures.py`.** 11 s.
- PNGs tier A. PDFs tier B.
- Use it.

**C6l `render_spieceasi_docx.py`.** 1.1 s.
- Tier B: paragraphs, tables and media are identical.
- Use it. Its hard-coded counts are a D2 matter.

**C25 `leftover_reads_analysis/aggregate_report.py`.** 0.1 s.
- Tier A for 5 of 5.
- It can be re-run after A4, but its genus labels come from `out/CAN_*.ont.json`.
  `parse_paf_16s.py` was **not run**: it reads minimap2 PAF from a pipe (`run_all_ont.sh`, a full re-mapping of 5
  Nanopore samples), and the PAF is not stored.
- C25 still needs a decision.

**Docs step 6: openpyxl load/save round trip** (`gate/xlsx_roundtrip.py`; the surgical-edit mechanism). 16 s.
- All four workbooks are **cell-identical**: `MAG_target_gene_filled_UPDATED (1).xlsx`, `MAG_target_gene_filled.xlsx`,
  `cladeexcel.xlsx`, `abundance and coverage (1).xlsx`.
- No zip parts are lost, and there is no change to merged cells, data validation, conditional formatting, column
  widths, freeze panes, filters or defined names.
- UPDATED (1) S1: 2,910 formulas with 0 cached values before and after. C17 is confirmed and does not get worse.

## 3. Findings for the orchestrator

1. **No network result needs a re-run for naming reasons.**
   - Aggregate, per-phase, all-phases, GAO/PAO (default mode) and SPIEC-EASI all reproduce byte-identically at HEAD.
     Louvain is seeded (42).
   - That makes a strong D2 test possible. After the consumers read the frozen 140-pair list, re-running C6b, C6h, C6j
     and C6k must give byte-identical outputs. Only fields added on purpose (bin, R232 MAG ID) may differ.
2. **The GAO/PAO figures are independent of D2** (C14).
3. **Several committed outputs were stale before R232.** Their scripts were never re-run after the 408cb39 ID update or
   later input changes. For "use (diff vs gate)" scripts, the R232 diff should be read against the gate re-run. The
   affected outputs:
   - B4 root tables;
   - gene_ab_figure 04/06/10/13 labels and the 08 `.md`;
   - the organism tree;
   - the module panels;
   - the nosZ reference tree;
   - the phase_rpkm_chord tables and figures;
   - the builder workbook ID columns (6a04588 map).
4. **The workbook S1 layout drift** (`rel_abun` became `mapped_reads`) breaks three scripts: `build_xlsx_291.py`,
   `12_novel_bacteria_summary.py` and `render_denitrification_heatmaps.py`. None of them may regenerate committed
   outputs. Use the fallbacks above.
5. **Do not regenerate** (fallback per section 2):
   - B4, B6, C2-05, C2-12, C3 nosZ reference tree, C4-01 (C4-03/04 only with accepted drift);
   - C5 denitrification heatmaps;
   - C6f/C6g module panels;
   - B9a unless the RAST fill is accepted.
6. **Untracked PDFs (C22).**
   - Stale in the main repo, by comparison with a HEAD re-render:
     - `network/N2O_sink_producer.pdf` (276-era);
     - `denitrification_*.pdf`;
     - `clade_classify/out/nosz_bb_reference_tree.pdf`;
     - `clade_classify/organism_tree/organism_tree.pdf`;
     - `gene_ab_figure/novel_{denitrification_bubble,metabolic_potential}.pdf` (Dokdonella.14).
   - Text-identical to the HEAD re-render: `fig1_tree.pdf`, `nosz_top5_tree.pdf`, `percell_composition.pdf`,
     `chord_*.pdf`, `ppk1_reference_tree.pdf` (tracked), and `spieceasi_*.pdf`.
7. **The B9 scripts read absolute paths.** `build_taxonomy.py`, `build_panel.py`, `validate_panel.py` and
   `call_genes.py` hard-code `/home/freiburger/Documents/EmilyKin`. Run in the main repo, they read the main repo, as
   intended. A future scratch-copy run must patch them, as done here.

## 4. Worktree-only patches (machine-specific paths; reverted after each run)

The patch files are in `rename_plan/gate/patches/`.

| patch | change |
|---|---|
| `build_skani_novelty_291.logiccheck.diff` | `SK` → `server_handoff_out/mag_vs_gtdb_skani.tsv` (276-era R220 raw); `TAXR` → `gate_scratch/inputs_external/gtdb_taxonomy_r220_from_sp_clusters.tsv` (built from `sp_clusters_r220.tsv`); `GTDBTK`, `MAGDIR`, `OUT` → repo-relative |
| `compute_manuscript_stats_291.R.diff` | `R` = `/scratch1/…` → `"."` |
| `build_taxonomy.EK.diff`, `build_panel.EK.diff`, `validate_panel.EK.diff`, `call_genes.REPO.diff` | `EK` / `REPO` = `/home/freiburger/Documents/EmilyKin` → the worktree, so that they do not read the main repo |
| `01_build_long_table.REPO.diff` | `REPO` = `/scratch1/…` → the script's parent directory |
| `05_bundle_offline.WORK.diff` | `WORK` = `/scratch1/…` → the script's directory |

- `build_xlsx_291.py` was run with `EMILYKIN_REPO=<worktree>`.
- No logic was patched.

## 5. Untracked inputs copied into the worktree

The files were copied, not symlinked, and are listed in `gate_scratch/copied_inputs.tsv`.

| purpose | files | size |
|---|---|---|
| `meta/results/01_nanoplot_raw/` (B5) | 90 | 14.2 MB |
| `ko_gene_abundance/counts_291_{,mag_}CAN_*.json` (B5/B6) | 10 | 2.4 MB |
| `gene_ab_figure/data/counts_CAN_*.json` (C4) | 5 | — |
| r13 `calls/genes.tsv`, `kofam/best_ko.tsv`, `kofam/hits.tsv`, `rast/crosswalk.tsv` (B9b) | 4 | 2.44 GB |
| r13 `rast/genomes/<bin>.tsv` (B9a) | 286 | 468 MB |

- `strace` found no other missing input. The one ENOENT on `MAG_target_gene_filled.pre_s6.bak.xlsx` is the backup
  existence test, not an input.

## 6. Not run

- **One-shot or poplar-bound generators**, per the rules: `regen_iterative_ids.py`, `apply_new_ids*.py`,
  `assign_mag_iterative_ids.py`, `build_mag_quality_table.py`, `build_mag_abundance_291.py`, `compose_static.py`,
  `03_aggregate_ko_rpkm.py`, `build_mimag_hq.py`.
- **`parse_paf_16s.py`** (C25), and a full `call_genes.py` run (B9e).
- **KBase scripts**: they belong to another stage.

## 7. Artifacts and teardown

- Harness: `rename_plan/gate/harness.py`. Usage: `run STEP [--cwd] [--env K=V] [--patch F] [--strace] -- CMD`.
- Supporting scripts in `rename_plan/gate/`: `asv_md5.py`, `copy_inputs.sh`, `xlsx_roundtrip.py`, `summarize.py`.
- Outputs: `rename_plan/gate_scratch/` (67 MB). HEAD re-runs that differ from the committed files are in
  `<step>/out/`, and side-by-side crops are in `_c*.png`.
- **Keep `gate_scratch/` until the R232 apply is verified.** The "use (diff vs gate)" verdicts compare R232 outputs
  with the HEAD re-runs at
  `/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan/gate_scratch/<step>/out/`.
  Re-running a step with `harness.py` deletes that step's folder first.
- The worktree is removed with `git worktree remove --force` followed by `git worktree prune`. The result is recorded
  below.

### Teardown record (2026-11-11 16:29 CST)

**Before removal:**
- worktree tracked changes: 0;
- main repo: `find -newer <stamp>` outside `.git` 0, `git status --porcelain` 0 lines, `md5sum -c tracked_md5_head.txt`
  OK;
- `~/Documents/emilykin_tools` and `~/Documents/gtdbtk_data`: 0 files newer than the stamp.

**Removal:**
- `git -C ~/Documents/EmilyKin worktree remove --force rename_plan/head_worktree` succeeded, then `git worktree prune`.
- `git worktree list` shows only `~/Documents/EmilyKin 9e087c3 [iterativeid-nomenclature]`.
- The `head_worktree` directory is gone.

**After removal:**
- HEAD is 9e087c3, on branch `iterativeid-nomenclature`.
- `git status --porcelain`: 0 lines.
- `md5sum -c tracked_md5_head.txt`: OK.
- `find -newer`: 0.
- `AMENDMENTS.sha256` (from `community_modeling/`): 14/14 OK.
- `asv_md5_before.tsv` re-check: 1,443 of 1,443 match.

**The worktree's untracked copies.** They were deleted with it: 395 files, 2.9 GB. Their sources in the main repo
were never touched.
