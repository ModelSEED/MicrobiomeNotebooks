# Betaine (TMA) models restricted to the 9 *grdH* carriers (GTDB R232)

This folder answers the R232 review comment on Results P76: "the *Halanaerobium* MAGs carry *grdH* but make no TMA in the models, and 6 of the 15 betaine models lack *grdH*. Consider re-running the betaine models on the 9 *grdH* carriers."

Nothing outside this folder was modified. That covers `correlational/models/`, `saltern_fluxes_betaine.json`, `saltern_fluxes.json` and `clauded_insights/`.

| File | What it is |
|---|---|
| `rerun_betaine_grdH_r232.py` | Re-simulates all 15 betaine models with the original settings and checks the result against `saltern_fluxes_betaine.json`. Also runs FVA on TMA, diagnoses each model and re-simulates the 9 *grdH* carriers. |
| `betaine_stats_grdH_r232.py` | First reproduces `clauded_insights/figure7_stats.json`, then computes the per-core modelled betaine CH4 and its statistics for each MAG set on each basis. |
| `inputs/averaged_normalized_abundances_unrestored_22f4629.csv` | The abundance table the original simulation actually used (`git show 22f4629:salterns/correlational/averaged_normalized_abundances_unrestored.csv`). |
| `fluxes_grdH_r232.json` | Fluxes for each MAG and core. It holds the 9 existing values, GLPK and Gurobi re-simulations of all 15 models, the FVA min/max envelope for all 15, and a sensitivity run of the 9 on the current abundance file. |
| `stats_grdH_r232.json` | Every statistic quoted below, plus the modelled CH4 per core and the saturation diagnostics. |
| `reproduction_check_r232.csv` | Largest difference from the stored fluxes for each MAG and solver, and whether each stored value lies inside its FVA range. |
| `fva_tma_by_core_r232.csv` | For each of the 15 MAGs and 6 cores: stored TMA, re-simulated TMA, FVA min/max, regime, lactate cost of betaine reduction, and whether the model is feasible without O2. |
| `diagnostics_grdH_r232.csv` | For each MAG: taxonomy, *grdH* evidence, where `rxn17220_c0` came from and its GPR, betaine transporters, annotated *grd* genes and whether any of them is in the model, and the sensitivity checks. |

Run both scripts with `~/Documents/py_venv/bin/python`. They use cobra 0.29.1 with GLPK, and Gurobi for FVA when a licence is present. The simulation script takes about 75 s and the statistics script a few seconds.

## Summary

1. **The review comment is partly right.** It is correct that 6 of the 15 modelled MAGs lack *grdH*. The three *Halanaerobium* models behave differently. *Halanaerobium*.3 is forced to make TMA (FVA minimum = maximum, 0.0105–0.1 mmol gDW⁻¹ h⁻¹), and *Halanaerobium*.2 is forced to make none. *Halanaerobium*.1's TMA is not unique: anywhere from 0 to 0.039–0.1, depending on the core. Its stored 0.027–0.1 is only the vertex GLPK returns, and Gurobi returns 0 in every core. An equally optimal solution therefore has two of the three *Halanaerobium* making no TMA, so the comment is partly defensible. The four zero-TMA models in the stored fluxes (the "four MAGs" of P76) are *Halanaerobium*.2 and BM714.4, which carry *grdH*, and BM714.1 and BM714.3, which do not.
2. **No model links its betaine reductase to *grd* genes.** The reaction `rxn17220_c0` was hand-added to all 15 models with an empty GPR by `add_betaine_reductase_pathway()` in `betaine_gapfilling.ipynb` (cell 2); MSGapfill did not add it. In each of the 9 carriers, the annotated betaine-reductase component B genes (EC 1.21.4.4; 1–8 genes per genome, IDs in the diagnostics CSV) appear in no model GPR. The 6 MAGs without *grdH* carry glycine-reductase genes (component B EC 1.21.4.2, protein A, component C) but no betaine-specific component B. Four of those six still make TMA, only through the gene-less `rxn17220_c0`.
3. **No structural fix is supported by the annotations, so no model was changed.** Neither zero-TMA *grdH* model is missing a reaction or transporter.
   - **BM714.4:** the zero is one of many equally optimal solutions.
   - ***Halanaerobium*.2:** the zero is forced, because betaine reduction costs lactate in this oxygen-respiring model.

   The only annotation-supported change would be a GPR on `rxn17220_c0`, which does not change any flux.

   Re-simulating the 9 carriers with the original settings, with GLPK in the original call order, reproduces their stored fluxes (largest difference in TMA 4.9e-15), so set (c) below is numerically the same as set (b). The match depends on the solver path, because TMA is not unique in 20 of the 54 carrier MAG×core cells. Gurobi returns a different TMA for 4 of the 9 (*Halanaerobium*.1, BM714.2, BM714.4 and *M55B127*.2).
4. **Restricting the set to the 9 *grdH* carriers barely changes Figure 7.** The table below uses the bases of `figure7_corrected.py` with n = 6 unrestored cores. Pearson p is the p-value of r.

| Set | Basis | Pearson R² | Pearson r (p) | Spearman ρ (p) |
|---|---|---|---|---|
| (a) original 15, stored fluxes | raw areal | 0.097 | −0.31 (0.55) | −0.086 (0.87) |
| | DNA-normalised (submitted Figure 7a) | 0.295 | −0.54 (0.27) | −0.43 (0.40) |
| (b) 9 *grdH*, stored fluxes | raw areal | 0.099 | −0.31 (0.55) | −0.086 (0.87) |
| | DNA-normalised | 0.285 | −0.53 (0.28) | −0.43 (0.40) |
| (c) 9 *grdH*, re-simulated, same settings | raw areal | 0.099 | −0.31 (0.55) | −0.086 (0.87) |
| | DNA-normalised | 0.285 | −0.53 (0.28) | −0.43 (0.40) |
| (d) 9 *grdH*, every value at its FVA minimum | raw / DNA-norm | 0.085 / 0.151 | −0.29 / −0.39 | −0.20 / −0.54 |
| (e) 9 *grdH*, every value at its FVA maximum | raw / DNA-norm | 0.141 / 0.290 | −0.38 / −0.54 | −0.20 / −0.54 |
| MePn comparator (`saltern_fluxes.json`) | raw / DNA-norm | 0.370 / 0.886 | +0.61 (0.20) / +0.94 (0.005) | +0.60 (0.21) on both |

   The modelled betaine CH4 per core for set (b) is R1_A 4032, R1_B 3202, R1_C 1666, R2_A 1900, R2_B 2036 and R2_C 1313 µmol m⁻² d⁻¹. The full per-core values for each set are in `stats_grdH_r232.json`.

   **The model-vs-measured values in this table are superseded.** Every modelled value, and the per-core µmol m⁻² d⁻¹ values above, uses the areal conversion of `figure7_corrected.py` (per-core volume-weighted DNA, 0.7 g cm⁻³). That basis is superseded by a dimensionally consistent one (per-MAG flux per unit of model abundance × relative abundance × DNA × layer thickness), on which the MePn models give ρ = −0.20, R² = 0.08 and the nine grdH betaine models ρ = −0.09, R² = 0.10. Do not quote the old numbers (MePn R² = 0.37 / ρ = 0.60, betaine R² ≈ 0.29) as model-vs-measured results. See `clauded_insights/modelling_consistent_basis_r232.json`. What the table still supports is the within-basis comparison of sets (a)–(e).

   **The negative ρ on the DNA-normalised basis is built in.** There, x = CH4/DNA is dominated by 1/DNA (R2_B and R2_C are about 5300, against ≤ 589 in the other cores), while y tracks DNA density (R² = 0.973, see item 6). CH4/DNA against DNA density alone gives ρ = −0.77. The sign comes from the opposite DNA terms on the two axes and says nothing about betaine.

5. **The text's R² = 0.862 (P74) comes from a third basis.** It is not the R² = 0.30 that submitted Figure 7a shows (row a). It is `fits_betaine.json`, computed in `analysis(1).ipynb`. There, x is the DNA-normalised CH4, and y is the summed betaine flux multiplied by a constant (2 µg DNA g⁻¹ and 1.4 g cm⁻³) instead of each core's DNA. I reproduce it to 8 decimals. On that basis the 9-MAG set gives R² = 0.953 (ρ = 0.64, p = 0.17). That fit rests on R2_B and R2_C, the two cores with one very low-DNA depth (R2_B_D1 and R2_C_D2, 0.16 µg g⁻¹ soil; that is 2.0 and 1.99 µg mL⁻¹ in the extract), which inflate CH4/DNA. The number is also unstable across equally optimal flux solutions. For the 15 MAGs, with every value at its FVA minimum or maximum, it falls to 0.82 or 0.19. P74's 0.862 should therefore not be quoted.
6. **TMA saturates, so the betaine fit measures DNA density, not abundance.**
   - **Bound:** in the stored 9-MAG fluxes, 33 of 54 MAG×core values (61 %) sit at the 0.1 betaine uptake bound and 12 are zero.
   - **Little response to abundance:** across the 6 cores, summed TMA varies only 1.32-fold, while the summed DNA-normalised abundance of the same MAGs varies 14.5-fold.
   - **DNA dominates:** the unit conversion multiplies each core by its DNA density, which varies 3.85-fold. As a result, modelled betaine CH4 against DNA density gives R² = 0.973 (0.976 for the 15).
   - **MePn contrast:** MePn CH4 flux is exactly proportional to abundance, at one of three constants per MAG (0.04321, 0.46974, 0.55344 CH4 per unit abundance). Its summed flux varies 17-fold across cores, tracking abundance (ρ = 0.94), and its modelled CH4 against DNA density gives only R² = 0.30.

   The betaine and MePn fits are therefore not like-for-like. A low betaine R² follows from the 0.1 uptake cap and the DNA multiplier; it is not evidence against betaine as a methane source. Conversely, the MePn fit adds no information beyond MePn-consumer abundance times a constant.

## 1. How the original simulation was run

| Item | Value (source) |
|---|---|
| Producer | `analysis(1).ipynb` at commit 22f4629 (cells 37 and 40), which wrote `saltern_fluxes_betaine.json` and `fits_betaine.json`. The same code is in `analysis(2).ipynb` (cells 33 and 36) and `analysis_2_marimo.py`. In the marimo copy, `_add_atpm(model)` uses the undefined name `_model`, so that copy would not run as is. |
| Models | `models/<MAG>.contigs__.RAST_gapfilled_betaine.json`. `*_betaine.json` has the same reaction set in all 15. The pipeline was: RAST draft → `MSATPCorrection` + `MSGapfill` on the aerobic lactate medium (`*_gapfilled.json`) → `add_betaine_reductase_pathway()`, which added `rxn17220_c0`, `rxn07207_c0`, `EX_cpd00441_e0` and `rxn_TMA_transport_c0` to all 15 models. |
| ATPM | `ATPM_c0` was added with bounds (0, 1000) whenever absent, which it was in all 15. |
| Medium | The fixed lactate medium with MePn removed: lactate 100, Pi 100, **O2 1000**, NH3 1000, and so on. Betaine was added as `EX_cpd00540_e0` = 0.1. |
| The 0.1 bound | A betaine **uptake** bound, not a TMA export cap; TMA export is open to 1000. FVA never finds TMA above 0.1 in any MAG or core, because each betaine yields one TMA. |
| Biomass | `bio1` was fixed to the core's abundance (lower bound = upper bound = abundance). |
| Objective | `EX_cpd00159_e0` with bounds (−1000, 0), maximised (the cobra default). This minimises lactate uptake at fixed biomass. |
| Abundances | Each core's value is the **mean over the two depths** of relative abundance divided by DNA. This is the 22f4629 version of `averaged_normalized_abundances_unrestored.csv`. It matches every abundance logged in `sols_log_betaine.txt` (for example, *UBA8670*.1 in R1_C is 0.0004784199386761). The current `data/averaged_normalized_abundances_unrestored.csv`, which uses a DNA-weighted depth merge and was rewritten after 2026-04-22, does **not** match the log. |
| How cores differ | Only through the fixed biomass, that is, the abundance. The medium and bounds are the same in every core. |

**Reproduction.** On the 22f4629 abundances, GLPK with the original call order reproduces the lactate objective for all 15 MAGs (to within 1.3e-11) and TMA to within 1e-9 for 14 of 15. The call order matters: copy the model, call `slim_optimize()` on the stored medium, then set the medium. Gurobi reproduces TMA for 10 of 15. Each mismatch, and the solver dependence itself, comes from a MAG and core where TMA is not unique. Every stored TMA value lies inside the FVA range at the lactate optimum. The original Mac environment probably had no Gurobi, which would make GLPK the faithful solver; this is an inference, not something checked.

## 2. Diagnosis

**TMA regimes** at the lactate optimum (gurobi FVA, `fva_tma_by_core_r232.csv`):

| Regime | 9 *grdH* (54 cells) | 6 without *grdH* (36 cells) |
|---|---|---|
| Fixed at the 0.1 uptake bound | 23 | 11 |
| Fixed below the bound (the only abundance-responsive case) | 5 | 6 |
| Not unique (alternative optima; the stored value is the solver's pick) | 20 | 19 |
| Forced to zero | 6 (*Halanaerobium*.2) | 0 |

**The four zero-TMA models:**

| MAG | *grdH* | Why TMA = 0 |
|---|---|---|
| *Halanaerobium*.2 | yes | **Forced.** The FVA maximum is 0 in every core. Betaine reduction costs lactate here, +0.043 mmol lactate per mmol TMA in R1_A. When betaine reduction is forced, its acetyl-phosphate and oxidised thioredoxin are recycled through the draft glycine-reductase reaction `rxn11676_c0` (glycine-synthesis direction) and leave as glycine; they produce no ATP. Closing `rxn11676_c0` still gives a TMA maximum of 0 (checked in R1_A only). The model needs O2 in every core, so respiration beats a 1-ATP Stickland-type step for reductant. Nothing is missing: the model has `rxn17220_c0`, both betaine transporters with GPRs (ABC `rxn05181_c0` and H+ symport `rxn05579_c0`), and TMA transport and exchange. |
| BM714.4 | yes (beta subunit only, weaker evidence than the other 8) | **Not unique.** TMA can be anywhere between 0 and 0.002–0.1 depending on the core, and betaine reduction is neutral to the objective (marginal cost ≈ 0). The stored 0 is the solver's pick. |
| BM714.1 | no | Same as BM714.4: TMA anywhere between 0 and 0.0009–0.1. |
| BM714.3 | no | Same: TMA anywhere between 0 and 0.0006–0.1. |

*Halanaerobium*.1 (*grdH*) is not in this table only because its stored TMA is non-zero. It is just as non-unique: TMA anywhere between 0 and 0.039–0.1. The stored 0.027–0.1 is GLPK's pick, and Gurobi returns 0 in every core.

**Why 4 of the 6 MAGs without *grdH* make TMA:** in every model, TMA comes from `rxn17220_c0`, the hand-added gene-less betaine reductase.
- *UBA8670*.1 and *CAZFJG01*.2 sit at 0.1 in every core.
- *DATFIM01*.1 has non-unique TMA in all cores, so its stored 0.0014–0.059 is arbitrary; it is the one MAG that GLPK does not reproduce.
- *CAZFJG01*.1 has unique values (0.0047–0.055) that scale with abundance.

**GPR check for the 9 *grdH* models:** `rxn17220_c0` has an empty GPR in all 9, and no betaine-reductase component B gene appears in any reaction GPR. The candidate genes for an annotation-based GPR are listed in `diagnostics_grdH_r232.csv` (column `betaine_reductase_compB_genes`).

**Minor inconsistencies:**
- **H+ coefficient:** `rxn17220_c0` carries −1 H+ in 4 models (*UBA8670*.1, *UBA8670*.2, *M55B127*.1, *M55B127*.2) and 0 in the other 11, a by-product of the automatic H+ balancing in `add_betaine_reductase_pathway()`.
- **Oxygen:** 65 of 90 MAG×core cells are infeasible without O2. The models that stay feasible without O2 in all 6 cores are *Halanaerobium*.1, *Halanaerobium*.3, *CAZFJG01*.1 and BM714.4; BM714.2 does so in 1 core. In other words, an anaerobic, Stickland-type enzyme is being simulated inside aerobic models. Halanaerobiales are generally described as anaerobes; that is background knowledge, not checked against a source here.

## 3. What was "corrected"

- **The fix applied:** the betaine guild is now the 9 *grdH* carriers. The 6 MAGs whose TMA comes only from a reaction without gene support were dropped.
- **The 9 carriers:** they were re-simulated with the original settings (set c), which reproduces the stored values.
- **Zero-TMA *grdH* models:** neither *Halanaerobium*.2 nor BM714.4 was altered, because no annotation-supported reaction or transporter is missing.
- **Non-unique TMA:** the FVA envelope (sets d and e) is reported instead of choosing a value. Treat sets d and e as two corner cases, not bounds on R², because each MAG and core can sit anywhere in its own range.
- **Sensitivity on the current abundance file** (DNA-weighted depth merge, 9 *grdH* re-simulated): raw R² = 0.142 with ρ = −0.26, and DNA-normalised R² = 0.397 with ρ = −0.54. The FVA corners give raw 0.149–0.185 and DNA-normalised 0.296–0.388. All of these are in `stats_grdH_r232.json` under `SENS_*`.

## Open issues (not fixed here)

- **Both Figure 7 panels need one abundance basis.** Both flux files (`saltern_fluxes.json` for MePn and `saltern_fluxes_betaine.json`) were computed on the 22f4629 depth-mean abundances. The MePn CH4/abundance ratio is constant only on that basis. `figure7_corrected.py` panel (b), and the current `data/` files, use the DNA-weighted merge. So panel (a) and panel (b) of Figure 7 rest on different normalisations of the same abundances.
- **The betaine models are aerobic.** They should be re-gapfilled and simulated anaerobically before betaine reduction can be judged mechanistically.
- **`rxn17220_c0` needs a GPR** from the component B genes if the models are kept.
- **The fixed 0.1 betaine uptake stays an assumption.** The text cites ref. 54 for it; I did not verify it. Because most models saturate at it, the modelled betaine CH4 mostly reflects this assumption.
- **The P76 wording does not match the data.** "Six MAGs identified as Dethiosulfatibacteraceae and Acidaminobacteraceae produced the maximum methane" fits poorly. The six MAGs at 0.1 in all 6 cores are:
  - *UBA8670*.1 and *UBA8670*.2 (Dethiosulfatibacteraceae);
  - *DATFJG01*.1 (Acidaminobacteraceae);
  - BM714.2 (BM714);
  - *CAZFJG01*.2 (JAFGMB01 in R232; Synergistales in 2020);
  - *M55B127*.2 (M55B127).

  *M55B127*.1 reaches the bound in 5 of 6 cores. Of these, *UBA8670*.1 and *CAZFJG01*.2 have no *grdH*.
