# Community methane simulations of the 2026 guild models (final MAG set) and Figure 7

Run 2026-11-23. Models: `../models/` (38 MePn consumers, 18 betaine reducers; see its README).

## Scripts

| Script | Where | What |
|---|---|---|
| `simulate_guilds.py` | TF (`~/tf_runs/salterns_sim_final`, `~/venvs/salterns_modelling`, GLPK; models md5-checked against local) | **MePn:** the original cell-27 protocol (analysis(1).ipynb at 22f4629, as in `clauded_insights/mepn_obligacy_check/analysis/vcommon.py`): ratio constraint MePn ≥ r × Pi (r = 0.07, 0.33), cell-5 medium, bio1 fixed to the core abundance, lactate uptake minimised. **Betaine:** the original protocol: medium without MePn + betaine uptake ≤ 0.1, same objective, TMA export plus its FVA range. **Abundance:** a = mean over D1, D2 of relative abundance / DNA (final585 CoverM table). **Also:** the MePn models with the RAST methane monooxygenase (rxn00844) closed. Output `fluxes_final.json`, `simulate_guilds_TF.log` (327 s). |
| `simulate_guilds.py --reference ../models reference_k.json` | local (seconds) | Every MePn model at a fixed biomass of 1e-3 gives each model's CH4 per unit biomass, k. The ratio binds in 76/76 runs. |
| `consistent_basis_final.py` | local | Areal comparison with the measured CH4 of the six unrestored cores, on the consistent basis of `clauded_insights/modelling_consistent_basis_r232.py`. MePn: k × biomass. Betaine: the pipeline conversion (Σ 0.06 × TMA × volume-weighted DNA). Sensitivities: relab_only, mmo_blocked, C–P lyase gene-supported subset, betaine rescaled and FVA min/max. Output `consistent_basis_final.json`. |
| `figure7_final.py` | local | `figure7_final.png/.svg`, `figure7_final_stats.json` |

**Solver tolerance.** At core abundances of about 1e-4 and below, the MePn fluxes fall under GLPK's feasibility tolerance:
- one cell, `R2_C_H2O__metabat2.bin.005` in R1_C at r = 0.33, is "infeasible";
- some cells return 0.

The analysis therefore takes k from the reference run, and asserts that it equals k in each model's highest-abundance core (maximum relative difference 4e-12).

## Results

| Set | Model / measured, total (range among cores) | Spearman ρ (exact p, n = 6) | R² |
|---|---|---|---|
| MePn, 38 models, r = 0.07 | 0.145 (0.047–0.328) | +0.03 (1.0) | 0.02 |
| MePn, r = 0.33 | 0.554 (0.180–1.265) | +0.03 (1.0) | 0.03 |
| MePn, methane monooxygenase closed, 0.07 / 0.33 | 0.152 / 0.578 | +0.03 | 0.04 |
| MePn, C–P lyase gene-supported (29), 0.07 | 0.125 | −0.09 | 0.03 |
| MePn, relative abundance only (constant DNA), 0.07 | 0.147 | +0.77 (0.10) | |
| Betaine, 18 models | 2.56 (1.60–5.53) | +0.60 (0.24) | 0.08 |
| Betaine, TMA at FVA minimum / maximum | 2.49 / 4.52 | +0.66 / −0.09 | |
| Betaine, rescaled to biomass | 6.58 | −0.26 | 0.16 |

**MePn, k per unit biomass.** k is set by the ratio and the template biomass P content:
- GramNeg V7: 0.0432 (28 models), and one at 0.0435;
- Archaea V6: 0.0498 (5) and 0.0501 (1);
- GramPos V7: 0.0765 (2).

**MePn, *Roseovarius*.2.** As built it releases no methane at r = 0.07 (most of it is oxidised at 0.33), because its RAST-annotated methane monooxygenase oxidises it at the lactate optimum.

**Shares of the simulated MePn methane** (r = 0.07): Rhodobacteraceae 75 %, *Roseovarius* 49 %, *Roseovarius*.10 alone 16.5 %.

**Betaine:**
- TMA is at the 0.1 bound in 35 of 108 cells, never zero, and not unique in 43.
- Summed TMA varies 3.0-fold across cores, against 13.2-fold for summed abundance.
- Simulated betaine methane against volume-weighted DNA gives R² = 0.15, so the old "betaine methane follows DNA (R² = 0.97)" no longer holds.

## Compared with the 2026-11-08 models (310 original MAGs, three repaired)

| | Old models | New models |
|---|---|---|
| MePn, r = 0.07 | 0.190 of measured | 0.145 |
| MePn, r = 0.33 | 0.721 | 0.554 |
| Betaine | 9 grdH models, ρ = −0.09 | exceeds measured in all six cores |
| MePn-obligate models | — | none |
