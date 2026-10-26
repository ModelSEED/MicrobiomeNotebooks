# Figure S11: what it takes to regenerate it with GTDB R232

Written 2026-10-26. This note changes no existing file. The prototypes and scripts named below live in the session scratch folder
`/tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns/a24fbb2e-9d33-5e05-aeea-a07d777ccde7/scratchpad/wf/figure-s11/`
(abbreviated `SCRATCH/` below). **That folder is per-session and temporary.** The prototype renderer `render_s11.py`, the
part worth keeping, is reproduced in full in the last appendix. Re-run it with
`~/Documents/py_venv/bin/python -I render_s11.py --salterns <salterns folder> --mode r232 --size 34 --label-font 12 --out <png>`.

## Short answer

- **The statistics are already done.** S11 is built from the within-condition correlation matrix plus the Graphical Lasso
  edges. Both already exist with R232 labels:
  - `C_within_r232.csv` (169 taxa);
  - `glasso_direct_edges_r232.csv` (326 edges, 186 with |partial r| > 0.05).
- **Only the renderer is missing.** `rerun_correlational_r232.py` never ported the plotting cell (the README lists
  `within_condition_triangle.png` as "not re-rendered").
- **The 2020 renderer still works.** Copied verbatim into scratch, it reproduces the manuscript image from the 2020 CSVs: the
  same 140 taxa, the same order, 138 asterisks, and the same 6 blue and 10 orange labels.
- **Pointed at the R232 inputs, it renders a 169-taxon version** with all guild taxa present, including the 13 archaeal taxa
  and *Methanolobus*.
- **Changes needed:**
  - a larger canvas or a smaller font;
  - guild highlight sets taken from `guilds_r232.csv` instead of the hard-coded 2020 MAG lists;
  - a rewritten caption, which is wrong even for the 2020 figure (see "Caption").
- **Effort:** about 2–3 h, plus the author's choice of which guild sets to highlight.

## Which image is S11

- **Location in the draft:** `Metabolic Manuscript.docx` P162 holds the picture and P163 the caption. The picture is
  `word/media/image10.png` (2006 × 2048 px), copied to `SCRATCH/FigureS11_2020_from_docx.png`.
- **What it shows:** its embedded title reads "Within-condition Pearson r on CLR residuals (asterisks: direct interaction in
  Graphical Lasso, |r_partial| > 0.05)". It is a lower-triangle heatmap of 140 taxa:
  - green asterisks mark Graphical Lasso edges;
  - 10 orange labels and boxed rows/columns mark the 2020 MePn set;
  - 6 blue ones mark the 2020 betaine set.
- **Source file:** `correlational/within_condition_triangle.png` (600 dpi, about 278 Mpx; PIL refuses to open it without
  raising its pixel limit). The docx copy is a downsample of it.

## The generating code

Everything is in `correlational/analysis_2_marimo.py`.

| Step | Lines | What it does |
|---|---|---|
| Cell `within_condition_correlation`, function `_wc_run` | 1926–2197 (def at 1959) | Parameters `EFFECT_THRESHOLD=0.7, FDR_ALPHA=0.05, MIN_COOCC=6, ZERO=1e-5`. |
| Depth merge | 1970–1982 | DNA-weighted D1/D2 merge to 12 cores (`data/Cliff_310MAG_relabund.tsv`, `data/SaltPondsDNA.csv`). |
| Presence filter | 1996 | Keeps a **MAG** with > 1e-5 in ≥ 6 of 12 cores. |
| Label, "defined IDs" filter, collapse | 1997–2004 | Relabels with `iterativeIDs.json`, **drops every MAG without an ID**, and sums MAGs into taxa (suffix dropped). |
| CLR | 2012–2016 | CLR with pseudocount = ½ × smallest non-zero value (3.705e-08 in both 2020 and R232). |
| Within-condition residuals | 2020–2024 | Subtracts the mean of each site (R1, R2, R2A, SF2) and computes `C_within` (Pearson). |
| Graphical Lasso | 2133–2144 | `GraphicalLassoCV(alphas=10, max_iter=200)` on the residuals gives `glasso_direct_edges.csv`. |
| Matrix output | 2149–2151 | Writes `C_marginal.csv`, `C_within.csv`, `C_between.csv`. |
| **Renderer**, cell `within_condition_triangle_heatmap` | **2200–2391** | Takes `wc_C_within` and `wc_glasso_df`. |
| ↳ Hard-coded 2020 guild MAG lists | 2237–2282 | `betaine_full` (15 MAGs) and `mepe_full` (27 MAGs), mapped to taxa by `_to_short` (2230–2235) through `iterativeIDs.json`. |
| ↳ Highlight colours | 2285–2286 | Betaine `#1f77b4` (blue), MePn `#e8671b` (orange). |
| ↳ Clustermap | 2288–2294 | `figsize=(28, 28)`, `coolwarm_r`, seaborn defaults: average linkage, Euclidean distance between matrix rows. Dendrograms hidden. |
| ↳ Labels | 2301–2322 | 16 pt; highlighted labels 19.2 pt bold. |
| ↳ Mask | 2324–2331 | Upper triangle masked. |
| ↳ Asterisks | 2338–2355 | GLasso pairs with \|partial r\| > 0.05. |
| ↳ Boxes | 2357–2366 | Boxed guild rows and columns. |
| ↳ Title | 2369–2371 | |
| ↳ Save | 2384 | `savefig("within_condition_triangle.png", dpi=600)`. |

Notes on the code:

- **No colour map is used.** The S11 cell has no `row_colors`/`col_colors`, so `iterativeID_color_map.json` and
  `Phylum_color_map.json` play no part in S11 (see "Colour maps" below).
- **The cell docstring is stale (2202–2208).** It says the asterisks mark |r_within| ≥ 0.7. The code marks GLasso edges.
- **A second, older within-condition cell (1746–1923) also writes `glasso_direct_edges.csv` and `direct_interactions.csv`.**
  - It uses a different presence filter: ≥ 3 of 12 cores.
  - It calls `DataFrame.groupby(..., axis=1)`, which pandas 3.0.3 (the venv) rejects: "unexpected keyword argument 'axis'"
    (tested).
  - The on-disk 2020 `glasso_direct_edges.csv` is the S11 version: 216 edges, 138 at |partial r| > 0.05, matching the marimo
    session log line "kept 138 of 216 edges". Re-running the notebook could overwrite it depending on cell order.
- **The R232 GLasso fit raised warnings.** The R232 run log (`network_r232.log`) shows scikit-learn RuntimeWarnings
  (overflow / invalid value in `dot`) during `GraphicalLassoCV`. The fit still completed (alpha = 0.5064, 326 edges). If the
  edges are recomputed, check that they reproduce.
- **The 2020 cell output is cached** in `__marimo__/session/analysis_2_marimo.py.json`:
  - 140 taxa, 9,730 pairs, 1,321 with |r_within| > 0.7;
  - 8,180 candidate pairs, 0 passing FDR;
  - maxT 95th-percentile |r| = 0.992;
  - 216 GLasso edges.

## Does it still run?

- **Not as a notebook.** `marimo` is not installed in `~/Documents/py_venv`.
- **The two cells behind S11 are plain Python and run under the venv** (pandas 3.0.3, seaborn 0.13.2, scikit-bio 0.7.2,
  scikit-learn 1.7.0).
  - `_wc_run` was already ported to `rerun_correlational_r232.py::within_condition()` (lines 374–487). Its README records
    that, with the 2020 IDs, the port reproduces the 2020 `C_within.csv` to 1e-16.
  - I re-derived the R232 `C_within` independently in scratch (`SCRATCH/absent_core_share.py`). It matches
    `C_within_r232.csv` to 1.1e-16.
  - The renderer, copied into `SCRATCH/render_s11.py` with only its inputs parameterised, reproduces the docx image from
    `C_within.csv` + `glasso_direct_edges.csv`. The output is `SCRATCH/S11_2020_reproduction.png`, which shows:
    - the same row order (top: Thiohalobacterales, Thiohalomonadales, UBA5704, GCA-001735895, JABDJL01, Halofilum; bottom:
      JAAYHV01, UBA12077, BM507, RZYO01);
    - 138 asterisks;
    - the same 6 betaine and 10 MePn labels.

## Inputs: 2020-specific vs. R232 equivalents

| Input | 2020 (used for S11) | R232 equivalent | Status |
|---|---|---|---|
| Abundance | `data/Cliff_310MAG_relabund.tsv` | same | shared |
| DNA weights | `data/SaltPondsDNA.csv` | same | shared |
| MAG → label | `correlational/iterativeIDs.json` (280 MAGs; **no archaea**) | `correlational/iterativeIDs_r232.json` (310 MAGs) | exists |
| Within-condition matrix | `correlational/C_within.csv` (140 taxa) | `correlational/r232/C_within_r232.csv` (169 taxa, 13 archaeal) | exists |
| GLasso edges | `correlational/glasso_direct_edges.csv` (216; 138 > 0.05) | `correlational/r232/glasso_direct_edges_r232.csv` (326; 186 > 0.05: 136 positive, 50 negative; CV alpha = 0.5064) | exists |
| Guild highlight sets | hard-coded `betaine_full` / `mepe_full` MAG lists in the cell | `correlational/r232/guilds_r232.csv` columns `mepn_consumer` (phnJ), `betaine_reducer` (grdH), `methanogen` (mcrA), `tma_producer`, `archaea`; plus `mepn_2020` / `betaine_2020` for the old sets | exists |
| Phylum / level per taxon | `iterativeID_phylums.json`, `iterativeID_levels.json` | `iterativeID_phylums_r232.json`, `iterativeID_levels_r232.json` | exist; not used by S11, but useful for italicising genera, as the R232 4-site heatmap does |
| Colour maps | `iterativeID_color_map.json`, `Phylum_color_map.json` | none | not used by S11 |

**S11 is built from `C_within`**, the Pearson correlation of within-site CLR residuals. It is not built from a raw-abundance
co-variation matrix. Neither `C_marginal` nor the 4-site Spearman matrix (`abundance_correlatons_one_triangle_*`) enters it.

## Changes needed for R232

### 1. Labels (63 shared, 77 only in 2020, 106 only in R232)

- **The 2020 figure has 140 taxon labels; the R232 one has 169.** Only 63 names are common to both. Some names vanish because
  R232 resolves 2020 family- or order-level labels to genera; others are new placeholder genus names.
- **2020 guild labels a reader will look for** (2020 label → R232 label of the same MAGs):
  - **MePn (2020 set):**
    - *Salinarimonas* → *Saliniramus*;
    - *Roseicyclus* → JAVJPJ01 + *Roseibacterium*;
    - IMCC34051 → DAINOF01;
    - *Maliponia* → CAZRZN01;
    - UBA6960 → *Gaopeijia*.
  - **Betaine (2020 set):**
    - Acidaminobacteraceae → DATFJG01;
    - Synergistales → CAZFJG01;
    - Halanaerobiales → CAZOON01 / DATFIM01 / M55B127.
  - ***Halomonas*:** the 2020 row summed *Halomonas*.1 and *Halomonas*.2. R232 splits them into *Halomonas* and
    *Vreelandella* (both rows present). This MAG split is separate from the isolate reclassification in the Methods.
- **The full 2020 → R232 mapping** is in the appendix.

### 2. Filter: parameters unchanged, but archaea now pass, and many guild taxa sit on the threshold

- **The caption's filter matches the code.** "Below 1e-5 in more than half of the 12 cores" is equivalent to keeping ≥ 6 of 12.
  The caption omits two details:
  - the filter is applied **per MAG, before** MAGs are summed into a taxon;
  - in 2020 every MAG without an iterativeID was then dropped.
- **2020:** 245 of 310 MAGs pass the presence filter. Of these, 222 carried a 2020 ID. The 23 archaeal MAGs that passed were
  dropped only for lacking an ID.
- **R232:** all 245 enter, giving 169 taxa, 13 of them archaeal.
- **Two guild MAGs fail the per-MAG filter, so their abundance is in no row:**
  - BM714.2, a grdH carrier detected in 5 of 12 cores (the BM714 row holds the 3 other BM714 MAGs);
  - *Thiolinea*.1, in the 2020 MePn set, 5 cores (absent from both the 2020 and the R232 figure).
- **Fifty of the 169 R232 taxa sit exactly at the threshold (detected in 6 cores).**
  - So are these guild MAGs: *Haloquadratum*.4, *Phaeodactylibacter*.1, the three *Halanaerobium*, BM714.1/.3/.4, both
    CAZFJG01, DATFIM01.1, UBA11579.1 and *Methanolobus*.1.
  - A filter one core stricter (≥ 7 of 12) would remove *Methanolobus*.1 and 5 of the 9 grdH MAGs: BM714.2 is already out,
    and BM714.4 and the three *Halanaerobium* would follow. It would also remove the phnJ MAGs *Haloquadratum*.4 and
    *Phaeodactylibacter*.1.
- **A caveat for the low-prevalence rows (computed in `SCRATCH/absent_core_share.py`).**
  - For a taxon absent from whole sites, the residual in those cores is the CLR of the pseudocount: it varies only with each
    core's geometric mean, not with the taxon.
  - Across the 110 taxa absent from at least one core, the median share of their residual sum of squares from absent cores is
    0.12. It exceeds 0.5 for 7 taxa.
  - ***Methanolobus*** (detected only in R1 and R2) **has 0.42.** DATFJG01 has 0.36; *Halanaerobium*, M55B127 and UBA8670 have
    0.19–0.23.
  - Inference: much of the *Methanolobus* row reflects the community composition of the cores where it is absent, not
    *Methanolobus* itself.

### 3. Guild highlight sets

The R232 sets come from `guilds_r232.csv`, collapsed to taxa. All are present in the R232 matrix.

| Set | MAGs | Taxa in R232 S11 |
|---|---|---|
| MePn consumers, phnJ | 22 | 10: *Albimonas*, *Desulfotignum*, *Haloquadratum*†, JAVJPJ01, *Marivita*, *Natrinema*†, *Phaeodactylibacter*, *Roseibacterium*, *Roseovarius*, *Saliniramus* |
| MePn, 2020 curated set | 27 | 14 of 15 (adds CAZRZN01, DAINOF01, *Gaopeijia*, *Halorubrum*†; *Thiolinea* filtered out) |
| Betaine reducers, grdH | 9 | 5: BM714, DATFJG01, *Halanaerobium*, M55B127, UBA8670 |
| Betaine, 2020 curated set | 15 | 7 (adds CAZFJG01, DATFIM01) |
| TMA-route candidates | 4 | 4: CAZFJG01, J153, T1Sed10-7, UBA11579 (weak evidence, no CutC; optional) |
| Methanogen, mcrA | 1 | *Methanolobus*† |

† = archaeal.

- **Highlighting is per taxon, so a box covers every MAG with that label.** For example, the *Halorubrum* row in the 2020-set
  variant sums 4 MAGs, of which only *Halorubrum*.3 has accessory *phn* genes.
- **No taxon falls in both the MePn and betaine sets.** CAZFJG01 is in both the 2020 betaine set and the TMA set. The renderer
  colours by first match, so the precedence order must be fixed explicitly.
- **GLasso edges (|partial r| > 0.05) touching the R232 guild taxa:**

| Taxon | Guild | Edges | Partners |
|---|---|---|---|
| DATFJG01 | grdH | 23 | |
| *Phaeodactylibacter* | phnJ | 4 | DT-110, CAZWDL01, DATFKA01, *Psychroflexus* |
| UBA11579 | TMA | 3 | |
| *Haloquadratum* | phnJ | 2 | *Halosegnis*, *Natronomonas* |
| M55B127 | grdH | 2 | |
| UBA8670 | grdH | 1 | DATFJG01 |
| J153 | TMA | 1 | |
| *Methanolobus*, *Halanaerobium*, BM714 and the other 8 phnJ taxa | | 0 | |

### 4. Archaea

- **13 archaeal taxa enter**, all Halobacteriota: DATFIX01, Haladaptataceae, *Haloarcula*, *Halohasta*, *Halonotius*,
  *Haloplanus*, *Haloquadratum*, *Halorubrum*, *Halosegnis*, Hikarchaeaceae, *Methanolobus*, *Natrinema*, *Natronomonas*.
- **In the prototype's clustering:**
  - eight haloarchaeal taxa form one block with GLasso asterisks (rows 108–118 of 169, counted from the top), consistent
    with the single haloarchaeal Louvain module in `within_condition_louvain_modules_r232.json`;
  - *Methanolobus* sits at row 7, next to *Desulfotignum* and *Halanaerobium*, the unrestored-only anaerobes.
- ***Natrinema* and *Haloquadratum* are phnJ carriers**, so with R232 guilds two MePn highlights are archaeal.
- **Marker:** the prototype appends "†" to archaeal labels rather than adding a fifth colour.

### 5. Colour maps (asked about; not used by S11)

S11 has no phylum colour bar in either version. If one is wanted, the 2020 JSON maps cannot be reused.

- **`Phylum_color_map.json` (14 keys, 2020 GTDB names) lacks 12 of the 21 phyla in R232 S11:** Actinomycetota, Atribacterota,
  Bacillota, Bipolaricaulota, Chloroflexota, Halobacteriota, Macinerneyibacteriota, Patescibacteriota, Pseudomonadota, SAR324,
  Synergistota, Thermotogota.
  - Across all 310 R232 MAGs, 18 phyla are missing. The 6 extra are Asgardarchaeota, Campylobacterota, Desulfobacterota_D,
    Fidelibacterota, Nanohalarchaeota and Thermoproteota.
  - The map was already incomplete for the 2020 figure's own phyla: Bipolaricaulota, Chloroflexota, Firmicutes_F,
    Mcinerneyibacteriota, Patescibacteria, SAR324, Synergistota and Thermotogota are missing.
- **`iterativeID_color_map.json` (94 taxon keys) lacks 126 of the 169 R232 S11 taxa** (it lacked 46 of the 140 2020 taxa).
  Missing:

  Algicoccus, Alkalispirochaeta, BM002, BM003, BM004, BM714, CAJXSD01, CAJXVY01, CAVTAV01, CAYPLZ01, CAYUNZ01, CAYVBJ01,
  CAZCJJ01, CAZFJG01, CAZLBJ01, CAZOON01, CAZRZN01, CAZUKV01, CAZVWW01, CAZWDL01, CAZYIP01, Christiangramia,
  Coraliomargarita_B, DAHVFK01, DAINOF01, DASQZJ01, DASQZK01, DASQZP01, DATFIM01, DATFIX01, DATFJC01, DATFJG01, DATFJT01,
  DATFJZ01, DATFKA01, DATFLG01, DATFLL01, DATFMH01, DATFML01, DATFMP01, DATFMT01, DATFMZ01, DATGCC01, Defluviilinea,
  Draconibacterium, Flagellimonas, GCA-001735895, Gaopeijia, Gilvirhabdus, Haladaptataceae, Halalkalibaculum, Halanaerobium,
  Haloarcula, Halohasta, Halonotius, Haloplanus, Haloquadratum, Halorubrum, Halosegnis, Hikarchaeaceae, Humimonas, JAACFI01,
  JAAXQD01, JAAYHV01, JABDJL01, JADGBM01, JAFDBW01, JAFMPW01, JAGXFP01, JAGXLC01, JAHEIH01, JAIPAH01, JAIPJO01, JAIVHP01,
  JAJRXK01, JAKDMP01, JAKVAV01, JAQYVZ01, JASBRS01, JASJAR01, JAUIBE01, JAVJPJ01, JAWLFI01, JAWLGC01, JAWOCM01, JAZFUB01,
  JBAUSV01, JBBDSH01, Lentisalinibacter, Lutibacter, M55B111, M55B112, M55B127, M55B154, M55B175, Marimicrobium,
  Methanolobus, Multivoradales, NaphS2, Natrinema, Natronomonas, Oleincola, Ornithinicoccus, Pelovirga, Phaeodactylibacter,
  Pleomorphochaeta, RI-121, Rhodohalobacter, Roseibacterium, SKVR01, SM23-33, SZUA-36, Salibaculum, Saliniramus,
  Salinivirga, Salsipaludibacter, Sunxiuqinia, Sva-07, UBA11574, UBA11579, UBA1847, UBA8670, UBA9300, VAYN01, VRUE01,
  Vreelandella.
- **Recommendation:** derive colours at render time from `guilds_r232.csv`, as `heatmap()` in `rerun_correlational_r232.py`
  (lines 300–305) already does. It uses turbo, with archaea in the first 15 % of the range. Do not extend the static JSONs.

### 6. Layout and fonts

- **Row height:** about `figsize × 0.92 × 72 / N` pt per row. 2020: 28 in / 140 taxa = 13.2 pt per row, already tighter than
  the 16 pt labels (adjacent labels touch in the 2020 image). R232 at 28 in: 169 taxa = 11.0 pt per row.
- **A direct port (28 in, 16/19.2 pt) overlaps badly** (`SCRATCH/S11_r232_direct_port.png`).
  - The worst collision is the run of adjacent highlighted rows 5–7: *Desulfotignum*, *Halanaerobium*, *Methanolobus*.
  - Their x-axis labels overprint each other.
  - Other adjacent highlighted pairs are JAVJPJ01–T1Sed10-7 (only if TMA is boxed) and *Roseibacterium*–BM714.
- **34 in with 12 pt labels restores 13.3 pt per row** (`SCRATCH/S11_r232_prototype.png`). Ordinary labels are legible.
  Highlighted labels at 14.4 pt still touch in the three-row run, so keep highlights at about 13 pt (1.1× not 1.2×) or bold
  only.
- **Save at dpi ≈ 300.** At 28 in, 600 dpi gives a file of about 278 Mpx that PIL rejects. Word only holds a downsample of
  about 2000 px anyway.
- **Legend:** the 2020 figure has none. Its colours are explained only in the caption. The prototype adds a guild/archaea
  legend in the empty upper triangle.

## Prototype result (scratch only)

| | 2020 (docx) | R232 prototype |
|---|---|---|
| Taxa | 140 (0 archaeal) | 169 (13 archaeal) |
| Pairs | 9,730 | 14,196 |
| Pairs with \|r_within\| ≥ 0.7 | 1,321 | 1,779 |
| FDR-significant | 0 of 8,180 candidates | 0 of 11,573 candidates |
| Asterisks (GLasso, \|partial r\| > 0.05) | 138 | 186 |
| Highlighted taxa | 10 MePn (2020) + 6 betaine (2020) | 10 phnJ + 5 grdH + 1 mcrA (+4 TMA in the variant) |

All prototype files are in `SCRATCH/`:

- `render_s11.py` (renderer; `--mode 2020|r232 --size --label-font --tma`);
- `S11_2020_reproduction.png`;
- `S11_r232_direct_port.png`;
- `S11_r232_prototype.png`;
- `S11_r232_prototype_with_TMA.png`;
- `guild_filter_check.py`, `label_colour_check.py`, `absent_core_share.py`.

The PNGs are rendered at 75 dpi for viewing. No colour gaps: S11 uses a continuous colormap plus 3–4 fixed guild colours.

## Caption (P163)

**Current text, for the record.** The python-docx paragraph dump shows "Pairs whose correlations exceed  are denoted…" with an
apparent blank. The blank is an OMML equation, `ρ=0.7`, which `paragraph.text` drops. The full sentence reads "Pairs whose
correlations exceed ρ=0.7 are denoted with a green asterisk inside their cell."

**Errors in the current caption** (they apply to the 2020 figure too):

1. **Asterisks.** They mark Graphical Lasso direct associations (|partial r| > 0.05; 138 pairs), not correlations above 0.7.
   1,321 pairs exceed |r| 0.7, and those are not marked. The figure's own title states the GLasso rule.
2. **Symbol.** "ρ" denotes Spearman, but the colours are Pearson r on within-site CLR residuals.
3. **Colour.** "red and blue": the MePn highlight is orange (`#e8671b`).
4. **Wording.** "subtracting per-condition average abundances to remove indirect relationships": subtracting site means
   removes co-variation caused by shared responses to the site. The step meant to remove indirect (conditional) relationships
   is the Graphical Lasso. (This is a judgement about wording, not a code error.)
5. **Missing facts:**
   - the filter is applied per MAG, and MAGs are then summed per taxon;
   - unbinned reads are excluded;
   - in 2020, archaea were excluded because they had no label;
   - no pair passes FDR, so the panel is descriptive.

**Proposed R232 caption** (numbers from the files above; edit the guild wording to whichever sets the paper adopts):

> **Figure S11.** Within-site co-variation of prevalent taxa. Relative abundances of the 310 dereplicated MAGs were merged
> across the two depths of each core, weighted by DNA concentration (12 cores). MAGs above 10⁻⁵ relative abundance in at
> least 6 of the 12 cores (245 MAGs) were summed by GTDB R232 taxon (169 taxa, 13 archaeal, marked †). They were then
> CLR-transformed (pseudocount = half the smallest non-zero abundance), and the mean of each site (R1, R2, SF2, R2A; three
> cores each) was subtracted, which removes co-variation caused by shared responses to site. Colours show Pearson
> correlations between these within-site residuals, with rows and columns ordered by hierarchical clustering (average
> linkage). Green asterisks mark the 186 pairs that are direct associations in the Graphical Lasso network (|partial r| >
> 0.05). No pair passed the permutation-based false discovery rate threshold (q < 0.05), so the heatmap is descriptive.
> Labels and boxed rows/columns mark methylphosphonate consumers (*phnJ*; orange), betaine reducers (*grdH*; blue), and the
> methanogen *Methanolobus* (*mcrA*; green).

**Minimal fix if the 2020 figure is kept.** Replace the asterisk sentence with "Green asterisks mark pairs that are direct
associations in the Graphical Lasso network (|partial r| > 0.05)", change "red" to "orange", and add "Archaea were not
included (they had no taxonomic label in the 2020 classification)."

**Already in the repo:** `make_r232_review_docx.py` (line 595) adds a reviewer comment on P163 saying S11 was not regenerated.
It does not edit P162/P163. In "Metabolic Manuscript - R232 review notes.md", item 10 and the "Not changed" list record
the same.

## Steps and effort

1. **Port the renderer** into `rerun_correlational_r232.py` as a new `--parts triangle`, about 1 h:
   - read `C_within_r232.csv`, `glasso_direct_edges_r232.csv` and `guilds_r232.csv`;
   - reuse `render_s11.py`'s logic;
   - assert there is no NaN in the matrix before clustering;
   - fix the guild colour precedence;
   - mark archaea;
   - italicise genera using `iterativeID_levels_r232.json`.

   This is a recommendation; the script was not edited.
2. **Decide the highlight sets** (author's call): phnJ/grdH/mcrA (gene-based, as in the rest of the R232 revision) or the 2020
   curated sets; TMA candidates boxed or not.
3. **Tune the layout**, about 0.5 h: 34 in, 12 pt labels, highlights ≤ 13 pt, dpi 300, legend in the upper triangle. Check the
   rows 5–7 run.
4. **Update the README** of `correlational/r232/` (output list) and the review notes, about 15 min.
5. **Replace the picture at P162 and the caption at P163** in the review docx builder (or by hand), and drop the "not
   regenerated" comment, about 0.5 h.
6. **Optional, about 0.5 h:** add the absent-core caveat for *Methanolobus* to the caption or text, or show the 2020-set
   variant beside it.

**Total: about 2–3 h.** No statistics need re-running; the slow parts of the pipeline (1,296-permutation null, maxT, GLasso CV)
are already done.

## Appendix: 2020 S11 labels whose MAGs carry other R232 names (83)

Computed over all MAGs that carry each 2020 label in `guilds_r232.csv`. Labels marked "(kept)" keep their name for some MAGs
and gain new names for others.

4572-130 → DATFJZ01; 70-9 → VAYN01, VRUE01; Acidaminobacteraceae → DATFJG01; Alphaproteobacteria → CAVTAV01;
BM103 → JAQYVZ01; BM507 → JAIOTG01, JAIPJO01; Balneolaceae → Balneolaceae (kept), CAYVBJ01, Fodinibius, Rhodohalobacter;
Balneolales → DATFMH01; Bipolaricaulia → M55B111; Cryomorpha → JAIVHP01, JAKVAV01; DT-110 → DATFKA01, DT-110 (kept);
DT-91 → JAHEIH01; Desulfohalobiaceae → CAYUNZ01; Desulfomicrobiaceae → JAWLGC01; Desulfurivibrionaceae → JADGBM01, JAGXFP01;
Flavobacteriaceae → Gilvirhabdus, JAFMPW01; GCA-2748055 → DASQZK01, M55B112; GWC2-71-9 → DATFMZ01; Gallionellaceae → RI-121;
Geothermobacteraceae → DATFMT01; Gramella → Christiangramia; HIMB11 → CAJXSD01, HIMB11 (kept); HLUCCA08 → Salibaculum;
Halanaerobiales → CAZOON01, DATFIM01, M55B127; Halioglobus → JAACFI01, Marimicrobium;
Halomonas → Halomonas (kept), Vreelandella; IMCC34051 → DAINOF01; JAACFB01 → UBA1847; JAAYYF01 → JASBRS01;
JABDQB01 → Humimonas; JABSQH01 → Sva-07; KS3-K002 → Humimonas; Krumholzibacteria → JAJRXK01; M0040 → JAFDBW01, M0040 (kept);
MWEM01 → CAZVWW01; Maliponia → CAZRZN01; Mcinerneyibacterium → CAZUKV01; Microbacteriaceae → CAJXVY01;
Muricauda → Flagellimonas; Nevskiales → CAZYIP01; Nitriliruptorales → Salsipaludibacter; Opitutales → CAZCJJ01;
Ornithinicoccus_A → Ornithinicoccus; PWYM01 → XN24; Paceibacterales → CAYPLZ01; RZYO01 → Pleomorphochaeta;
Rhodobacteraceae → JASJAR01; Rhodothermales → DATGCC01; Roseicyclus → JAVJPJ01, Roseibacterium; SAR324 → DATFMP01;
SCGC-AAA027-K21 → Algicoccus; SG8-39 → DAHVFK01; SHZL01 → SZUA-36; SLDV01 → Pelovirga; SLND01 → JAKDMP01, JBAUSV01;
SURF-3 → M55B175; SZUA-36 → JAUIBE01; SZUA-401 → JAWOCM01; Salinarimonas → Saliniramus;
Saprospiraceae → CAZWDL01, DATFJT01, JAWLFI01; Sediminispirochaetaceae → DSBN01; Synergistales → CAZFJG01;
Thiohalobacterales → DATFLG01, DATFML01; Thiohalomonadales → DATFLL01; UBA10428 → DASQZP01; UBA11386 → JAZFUB01;
UBA12077 → JAGXLC01, UBA12077 (kept); UBA12170 → JAIPAH01; UBA12270 → SM23-33; UBA12294 → Defluviilinea;
UBA4486 → JBBDSH01; UBA4802 → M55B154; UBA5335 → DASQZJ01; UBA5704 → CAZLBJ01, Multivoradales; UBA6264 → DATFJC01;
UBA6794 → Oleincola; UBA6960 → Gaopeijia; UBA7960 → SKVR01; VadinHA17 → UBA9300; WC007 → Draconibacterium;
WN38 → Coraliomargarita_B; Woeseiaceae → Lentisalinibacter; YR4-1 → Halalkalibaculum.

Some R232 names in this list (e.g. JAIOTG01) belong to MAGs that fail the presence filter, so they do not appear in the R232
figure.

## Appendix: prototype renderer `render_s11.py` (copied from the scratch folder)

Rendering logic is that of `analysis_2_marimo.py` lines 2200–2391; only the inputs, highlight sets, sizes and output path are parameterised. Seaborn `clustermap` defaults (checked in 0.13.2): `method='average'`, `metric='euclidean'`.

```python
"""Scratch prototype for manuscript Figure S11 (within-condition triangle heatmap).

Rendering logic copied from the `within_condition_triangle_heatmap` cell of
correlational/analysis_2_marimo.py (lines 2200-2391); only the inputs, highlight
sets, sizes and output path are parameterised. Nothing is written into the repo.

  mode 2020 : C_within.csv + glasso_direct_edges.csv, 2020 hand-curated guild lists
              (the betaine_full / mepe_full MAG names of the cell, mapped through
              iterativeIDs.json exactly as `_to_short` does)
  mode r232 : r232/C_within_r232.csv + r232/glasso_direct_edges_r232.csv, gene-based
              guilds from r232/guilds_r232.csv (phnJ, grdH, mcrA); archaea marked
"""
import argparse, json, re, sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns

ap = argparse.ArgumentParser()
ap.add_argument("--salterns", required=True)
ap.add_argument("--mode", choices=["2020", "r232"], required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--dpi", type=int, default=150)
ap.add_argument("--size", type=float, default=28.0, help="figure size in inches (square)")
ap.add_argument("--label-font", type=float, default=16.0)
ap.add_argument("--tma", action="store_true", help="r232: also box the 4 TMA-route candidate taxa")
a = ap.parse_args()

S = Path(a.salterns); COR = S / "correlational"

if a.mode == "2020":
    C = pd.read_csv(COR / "C_within.csv", index_col=0)
    gl = pd.read_csv(COR / "glasso_direct_edges.csv")
    ids = json.load(open(COR / "iterativeIDs.json"))
    src = (COR / "analysis_2_marimo.py").read_text()
    cell = src[src.index("def within_condition_triangle_heatmap"):src.index("def within_condition_network")]
    def mag_list(var):
        body = re.search(rf"{var} = \[(.*?)\]", cell, re.S)[1]
        return re.findall(r'"(Salt_Pond[^"]+)"', body)
    def to_short(full_id):              # verbatim logic of _to_short in the cell
        bare = full_id
        for suf in (".contigs__", ".contigs"):
            if bare.endswith(suf): bare = bare[: -len(suf)]; break
        iid = ids.get(bare, bare)
        return ".".join(iid.split(".")[:-1]) if "." in iid else iid
    groups = [  # (name, colour, taxa) - order = priority, as in the cell (betaine first)
        ("Betaine consumers (2020 curated set)", "#1f77b4", {to_short(x) for x in mag_list("betaine_full")}),
        ("MePn consumers (2020 curated set)", "#e8671b", {to_short(x) for x in mag_list("mepe_full")}),
    ]
    archaea = set()
    title = ("Within-condition Pearson r on CLR residuals  "
             r"(asterisks: direct interaction in Graphical Lasso, $|r_{\rm partial}| > 0.05$)")
else:
    C = pd.read_csv(COR / "r232/C_within_r232.csv", index_col=0)
    gl = pd.read_csv(COR / "r232/glasso_direct_edges_r232.csv")
    g = pd.read_csv(COR / "r232/guilds_r232.csv")
    g["taxon"] = g.iterativeID_r232.str.rsplit(".", n=1).str[0]
    tx = lambda col: set(g.loc[g[col], "taxon"])
    groups = [
        ("Methanogen (mcrA)", "#009E73", tx("methanogen")),
        ("Betaine reducers (grdH)", "#1f77b4", tx("betaine_reducer")),
        ("MePn consumers (phnJ)", "#e8671b", tx("mepn_consumer")),
    ]
    if a.tma:
        groups.append(("TMA-route candidates (cntA / cutC-cluster genes)", "#8e44ad", tx("tma_producer")))
    archaea = tx("archaea")
    title = ("Within-condition Pearson r on CLR residuals, GTDB r232 taxa  "
             r"(asterisks: Graphical Lasso direct interaction, $|r_{\rm partial}| > 0.05$)")

assert not C.isna().any().any(), "NaN in C_within - clustering would fail"
N = len(C)
LABEL_FONT = a.label_font
HIGHLIGHT_FONT = LABEL_FONT * 1.2
TITLE_FONT, STAR_FONT, CBAR_LABEL_FONT, CBAR_TICK_FONT = 30 * a.size / 28, 18 * a.label_font / 16, 18, 14
STAR_COLOR, BOX_LW = "#90EE90", 2.0

def group_of(t):
    for name, col, members in groups:
        if t in members:
            return name, col
    return None, None

def shown(t):   # archaea get a dagger suffix in the r232 prototype
    return f"{t} †" if t in archaea else t

cm = sns.clustermap(C, cmap="coolwarm_r", center=0, vmin=-1, vmax=1, figsize=(a.size, a.size),
                    dendrogram_ratio=(0.08, 0.08), xticklabels=True, yticklabels=True)
cm.ax_row_dendrogram.set_visible(False); cm.ax_col_dendrogram.set_visible(False)
ax = cm.ax_heatmap
ax.yaxis.set_ticks_position("left"); ax.yaxis.set_label_position("left")
for lbls, rot in ((ax.get_xticklabels(), 90), (ax.get_yticklabels(), 0)):
    for lbl in lbls:
        t = lbl.get_text(); _, col = group_of(t)
        lbl.set_rotation(rot)
        if col:
            lbl.set_color(col); lbl.set_fontweight("bold"); lbl.set_fontsize(HIGHLIGHT_FONT)
        else:
            lbl.set_fontsize(LABEL_FONT)
        lbl.set_text(shown(t))
    if rot == 90: ax.set_xticklabels(lbls)
    else: ax.set_yticklabels(lbls)

ro, co = cm.dendrogram_row.reordered_ind, cm.dendrogram_col.reordered_ind
reo = C.iloc[ro, co]
mask = np.triu(np.ones_like(reo, dtype=bool), k=1)
mesh = ax.collections[0]
arr = mesh.get_array().reshape(reo.shape); arr[mask] = np.nan; mesh.set_array(arr.ravel())

direct = set()
for _, r in gl.iterrows():
    if abs(r["partial_r"]) > 0.05:
        direct |= {(r["taxon_a"], r["taxon_b"]), (r["taxon_b"], r["taxon_a"])}
ox, oy = list(reo.columns), list(reo.index)
n_star = 0
for i in range(N):
    for j in range(N):
        if mask[i, j] or i == j: continue
        if (oy[i], ox[j]) in direct:
            ax.text(j + 0.5, i + 0.5, "*", ha="center", va="center", color=STAR_COLOR,
                    fontsize=STAR_FONT, fontweight="bold", clip_on=False); n_star += 1
for i, t in enumerate(oy):
    _, col = group_of(t)
    if not col: continue
    ax.add_patch(patches.Rectangle((0, i), i + 1, 1, fill=False, edgecolor=col, linewidth=BOX_LW, clip_on=False))
    ax.add_patch(patches.Rectangle((i, i), 1, N - i, fill=False, edgecolor=col, linewidth=BOX_LW, clip_on=False))
ax.set_xlabel(""); ax.set_ylabel("")
ax.set_title(title, fontsize=TITLE_FONT, pad=20)
hm = ax.get_position()
cax = cm.cax
cax.set_position([hm.x0 + hm.width * 0.78, hm.y0 + hm.height * 0.62, hm.width * 0.025, hm.height * 0.30])
cax.tick_params(labelsize=CBAR_TICK_FONT); cax.set_xlabel(""); cax.set_ylabel("")
cax.set_title("Pearson r", fontsize=CBAR_LABEL_FONT, pad=10)
present = [(n, c, m & set(C.index)) for n, c, m in groups]
handles = [Line2D([], [], ls="", marker="s", ms=14, mfc="none", mec=c, mew=2.5, label=f"{n}: {len(m)} taxa")
           for n, c, m in present if m]
if archaea:
    handles.append(Line2D([], [], ls="", marker="$†$", ms=14, color="black",
                          label=f"archaeal taxon: {len(archaea & set(C.index))} taxa"))
if a.mode == "r232":
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.55, 0.55), fontsize=18 * a.size / 28,
              frameon=True, title="Label / box colour", title_fontsize=20 * a.size / 28)
cm.figure.savefig(a.out, dpi=a.dpi, bbox_inches="tight")
pt_per_row = a.size * (1 - 0.08) * 72 / N
print(f"mode={a.mode} taxa={N} asterisks(lower triangle)={n_star} figsize={a.size} in dpi={a.dpi} "
      f"~{pt_per_row:.1f} pt per row, label font {LABEL_FONT} pt (highlight {HIGHLIGHT_FONT:.1f})")
for n, c, m in present:
    print(f"  {n}: {sorted(m)}")
print("  top rows:", oy[:6], " bottom rows:", oy[-4:])
json.dump({"order": oy}, open(Path(a.out).with_suffix(".order.json"), "w"))
```
