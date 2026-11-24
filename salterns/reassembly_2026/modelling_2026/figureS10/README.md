# Figure S10 for the 2026 guild models

The original Figure S10 ("Core Gap Filling - MePn consumers / BET consumers") showed, for each draft model of the 2020 set, how many reactions ModelSEEDpy's ATP correction (MSATPCorrection, core gap-filling) had to add for the draft to make ATP on each of 22 ATP test media. Its heatmap layout is reproduced here for the 2026 arm-A drafts of the 38 MePn-consumer and 18 betaine-reducer models.

- **`../pipeline/atp_media_scan.py`** rebuilds each draft exactly as `build_model.py` does before ATP correction:
  - RAST draft, plus ATPM;
  - for MePn consumers, also the cell-5 disables and the MePn module.

  It then runs `MSATPCorrection.evaluate_growth_media()` and records, per medium, the reactions added + reversed (0 = makes ATP already; null = no solution).
  - Ran on TF on 2026-11-23 (`~/tf_runs/salterns_modelling`, same code md5 and inputs as the build, Gurobi, 14 in parallel, about 16 s per model).
  - Output: `scan/<mag>__<guild>.json`, 56 files.
- **`figureS10_final.py`** writes:
  - `figureS10_final.png/.svg`, two clustered heatmaps with LLac.O2 (the simulation medium) boxed;
  - `figureS10_counts.tsv`;
  - `figureS10_summary.json`.

**Result:**
- **MePn consumers:** aerobic L-lactate is tied for the cheapest medium, at a median of 0 reactions.
- **Betaine reducers:** it ranks 16th of 22 (median 5.5). Anaerobic pyruvate and fumarate need none.
- **Manuscript:** P37's claim that lactate needed the fewest reactions "for all MePn and betaine consumers" is corrected, with a comment.
