# prFBA augmentation: scouting results (2026-11-05)

These four read-mostly scouts prepared a fresh prFBA augmentation of the saltern MAGs. Per the user, no existing augmentation was used, and all augmentation is to be regenerated from prFBA.

The augmentation itself was postponed. The 39 metagenomes are being re-assembled and re-binned in KBase (`../../reassembly_2026/`), so it should run on the new MAGs.

The scripts here are prototypes copied out of the session scratchpad. Their paths still point at that scratchpad, so treat them as working notes rather than tools. `scout_results_full.json` holds every scout's complete structured output.

## 1. prFBA MAG arm (`mag_arm/`)

**What works unchanged**
- `mag_skani_hits.py search` (skani 0.3.2 against the GTDB r226 sketch, `-t 16`) works as shipped.
- `select_references.py --knobs-json ani_knobs.json` also works as shipped.
- Every MAG gets one reference, so p = 1 for every reference function.

**What does not work, or is missing**
- **Accession mapping:** GTDB accession → BV-BRC genome_id is not wired in. `load_acc_map` has no callers, and the conversion is listed under `integration_todo`. It has to be applied by hand.
- **The "fast" gene builder is broken in two ways:**
  - It sends no User-Agent, so BV-BRC returns 403.
  - It does not filter features by annotation source, so it mixes PATRIC/RASTtk and RefSeq CDSs and roughly doubles the function strings.
- **Synthetic-genome builder:** `build_synthetic_genomes*.py` does not run on a MAG selection unchanged. It uses the 16S thresholds and refuses estimated gene data.
- **No MAG annotation merge:** the MAG's own RAST annotation is never merged into the synthetic genome. A merge prototype is in `map_fetch_build.py`: MAG features at p = 1, reference-only functions at k/n, each with its provenance.

**Smoke test (3 MAGs):** none of the marker genes is truly absent from any of the three MAGs. Every apparent absence is a RefSeq-named duplicate of a gene the MAG already has, or a false positive on the name.

## 2. Coverage of the 310 old MAGs (`coverage/`)

**R226 skani run on all 310** (`skani_r226_310/`; 3:49 wall, about 37 GB RSS):
- **Resolved, 193 of 310:**
  - 125 at species level (ANI ≥ 95);
  - 65 at genus level (85–95);
  - 3 at family level.
- **Abstained, 117:**
  - 67 with too little alignment (low AF);
  - 48 with no hit;
  - 2 below ANI 83.
- **By group:**
  - MePn consumers: 22 of 22 resolved.
  - grdH carriers: 8 of 9.
  - TMA candidates: 0 of 4.
  - Archaea: 17 of 30.

**R232 as a ceiling:** 241 MAGs reach a tier. 48 of them are lost against r226, mostly because their best reference exists only in R232.

**Reference quality:**
- Of the 176 MAGs at ANI ≥ 95, 107 sit at ≥ 99.5.
- Their references are mostly other MAGs: 65 from the SMAG soil-MAG project and 3 deposits of this study's own MAGs.
- Only 10 references are isolates.
- An augmentation therefore often inherits the same incompleteness as the MAG.

**Fetching and resource traps:**
- The shared 5.1 GB PGFam cache needs 23–40 GB of RAM to load, and it is rewritten after every fetch. Use a small new cache instead.
- 1,800 cached empty entries are stale and need re-fetching.
- BV-BRC takes 0.7–1.6 s per genome with a browser User-Agent.

## 3. Analyses that depend on annotation (`inventory/`)

**Single source:** all gene content enters through `build_guilds()` in `correlational/r232/rerun_correlational_r232.py`, which reads `RASTified_MAG_json` and writes `guilds_r232.csv`.

**Swap seam:** `rerun_with_alt_annotations.py` swaps in another annotation source and re-runs the dependent scripts. Identity runs reproduce every output byte for byte.

**What breaks or drifts on a real swap**
- **Stops:**
  - hard-coded count asserts (membership l.265–268, l.371; S11 l.45; betaine l.269; consistent basis l.104–105; numbers l.129);
  - `next(...)` name lookups in the archaeal report.
- **Silent drift:**
  - the shared RNG (`default_rng(0)`, l.157) moves the random-set fractions even for unchanged guilds;
  - the docx prose and the Table 4 cells are literals.

**Vocabulary:** the guild regexes are specific only within the 2020 RAST vocabulary. On BV-BRC product names:
- phnJ also matches PhnG/H/K/L and C-P lyase regulators, and misses the RefSeq name "C-P-lyase PhnJ";
- grdH, cutC_cluster and cntA all become less specific.

**Table S5:** `make_table_s5.py` regenerates Table S5 exactly; the workbook previously had no generator.

## 4. Modelling with augmented annotations (`modelling/`)

**Verdict:** both arms can be rebuilt, gapfilled and simulated locally with one pipeline (`build_pipeline.py`):
- ModelSEEDpy;
- GramNeg template v7 plus Core v6;
- a hand-added MePn ABC transporter;
- ATP correction;
- gapfilling on lactate, then on the Pi-free medium;
- the original protocol from analysis(1) cell 27.

A model takes about 2–3.5 minutes and 1.3 GB.

**Key finding: the MePn methane of a ratio-bound model does not depend on its genome.**
- In 34 of 40 stored MAG × ratio cases, methane per unit biomass = r/(1 + r) × 0.660485, where 0.660485 is the template biomass's phosphorus content (error ≤ 1e-10).
- The three high-flux MAGs (*Roseovarius*.9, *Roseovarius*.11, *Albimonas*.1) are exactly the 3 of 24 original models that cannot grow without MePn. In *Roseovarius*.9, glycine cannot be made without MePn.
- Those three produce 72% of the simulated methane. A rebuilt *Roseovarius*.9 grows without MePn, so their MePn dependence may be a gapfilling artefact. **Verify this before relying on the "about half of the measured methane" statement.**

**Consequence for an original-vs-augmented comparison:** the two arms differ only where augmentation flips a regime: MePn obligacy, a methane sink, or infeasibility.

**Protocol cautions**
- Gapfill first on a medium without MePn or betaine. Mock arm b1 became MePn-obligate after gapfilling on the MePn-containing medium.
- prFBA (ModelSEEDpy-prfba worktree) cannot run the saltern objective unchanged. Use a lexicographic second stage instead.

**Security note:** `salterns/cliff_mags/mac/salterns/build_models.ipynb` (cell 1) and `gapfill.ipynb` (cell 86) contain a hard-coded KBase token.
- Both files are tracked in the public repository, since commit c708d10 (2026-01-13).
- It is not the token currently in `~/.kbase/token`.
- Its owner should revoke it.
