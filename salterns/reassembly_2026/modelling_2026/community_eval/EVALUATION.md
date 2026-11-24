# Would the EmilyKin or codiffusion_bioreactor community-modelling pipeline reproduce the saltern methane better?

Evaluated 2026-11-23 by Claude, for Andrew.
- **Pipelines:** surveyed read-only by two agents, with file:line references in their reports. Nothing in either project was modified.
- **Ceiling test:** `input_signal_check.py` (→ `input_signal_check.tsv`).

## Short answer

No, not for the between-core differences, and only by assumption for the total.
- **The two pipelines' own track records.** Their agreement with measured gas fluxes comes from an *imposed* input flux, not from the community model. The parts the community model decides do not track sample-to-sample differences.
- **The salterns' own inputs.** None of the measured inputs ranks the six unrestored cores' methane, so no formulation built on them can be expected to.

## 1. codiffusion_bioreactor (prFBA + MSCommunity, methane bioreactor)

**What it is:**
- **prFBA:** probabilistic FBA. Each 16S ASV gets a synthetic genome merged from its best-hit BV-BRC references. The per-reaction prevalence weights are stored but not used in the simulations.
- **Community model:** an MSCommunity compartmentalised model with a shared e0. Abundances enter as biomass coefficients, and growth is fixed at 1/SRT.
- **Measured spec:** each day's measured H2 and CO2 consumption, medium supply and VFA effluent caps are imposed. CH4 is left free and reported as an FVA interval.

**Validation, 11 days, v7fo:**
- Measured CH4 lies inside [min, max] on 11/11 days.
- The **maximum ≈ H2/4**: the electron ceiling set by the imposed H2 consumption. Its ρ = 0.99 with measured CH4 follows from the input.
- The **minimum**, the part the community decides, has ρ = −0.13 (p = 0.71).
- The project's own manuscript says this "checks the stoichiometry of the models rather than validating them".

**Transfer to the salterns:**
- The salterns cores have no donor-consumption rate, growth rate or supply rates, only standing concentrations.
- Microcosm rates would re-introduce the near-tautology.
- **Compute depends on scope:**
  - All 585 species would be about twice the project's largest day, which took about 10 h, and one of its solves returned a false zero.
  - A guild community of the 56 modelled species would be small and cheap.

## 2. EmilyKin (fork of the codiffusion pipeline; denitrifying P-removal reactor, N2O)

**What it is:**
- the same MSCommunity structure;
- an export policy that closes CH4 and every C/H/O compound more reduced than acetate;
- curated, eQuilibrator-oriented templates, element closure and Gibbs-budget rows;
- round-13 member caps, maintenance and storage.

**Validation:**
- Round 12 (75 days): N gas is **0 at every reported vertex**, and the model takes up 0.355 of the measured nitrite.
- Round 13: the member caps, maintenance (a 32-fold spread, falsified by leave-one-out) and the staged-LP gate all **failed** their pre-registered criteria.
- Across days, model nitrite uptake vs measured N reduction gives ρ = 0.19. The *input* nitrite cap alone gives ρ = 0.74.
- The methane ancestor gave a CH4 interval of [0, ~50] against a measured 6–35: "the constraints did not pin the output".

**Transfer:** the same missing-input problem. The export policy would also have to be opened for CH4 and TMA.

## 3. The ceiling test: can anything measured rank the six cores?

Thirty candidate inputs were fixed before looking at the results:
- guild biomass and relative abundance;
- soil DNA;
- 11 NMR metabolites (methanogenic substrates, donors, osmolytes);
- P, S, N, salinity and TOC;
- three guild × substrate products.

Each was tested with Spearman's ρ against measured CH4, n = 6, exact permutation p (|ρ| ≥ 0.886 is needed for p < 0.05).

| Input | ρ | exact p |
|---|---|---|
| Porewater betaine | 0.83 | 0.058 |
| SO4, N,N-dimethylglycine, trehalose | 0.77 | 0.10 |
| MePn-consumer relative abundance | 0.66 | 0.18 |
| Soil DNA | −0.60 | 0.24 |
| Betaine-reducer biomass × betaine | −0.26 | 0.66 |
| **MePn-consumer areal biomass** (= what the MePn model computes) | **0.03** | **1.0** |
| MePn biomass × N/P | 0.03 | 1.0 |

- **No input reaches p < 0.05.** About 1.5 nominal hits would be expected by chance among 30.
- **The MePn model's ρ = 0.03 is the ρ of its input.** At a binding phosphate-source ratio, simulated CH4 = constant × MePn-consumer biomass. FBA adds nothing between cores, a point both surveys reached independently.
- **The areal conversion erases the weak abundance signal.** MePn-consumer *relative* abundance has ρ = 0.66, but soil DNA is anti-correlated with CH4 (ρ = −0.60). Multiplying by DNA therefore removes it. On the relative-abundance-only basis the models give ρ = 0.77 (p = 0.10, `consistent_basis_final.json` "relab_only").

**The DNA values of two depths deserve a check before anything is concluded about tracking.** R2_B_D1 and R2_C_D2 both read exactly 0.16 µg DNA g⁻¹ soil, against 1.6–4.8 for the other ten unrestored depths, and they belong to the two highest-methane cores.
- Every areal conversion multiplies by DNA, so these two values pull the simulated emission of R2_B and R2_C down.
- Sensitivity (not a result): replacing them with their core's other depth, or with the unrestored median, raises the MePn model's ρ from 0.03 to 0.37 (p = 0.50) and its total from 14 % to 16–17 % of measured. That helps, but it is still no tracking.
- If 0.16 is a detection-limit placeholder or an entry error, the areal model numbers in the manuscript move. This is a measurement question for the authors, not a modelling one.

## 4. What would and would not help

**Would not help** (added free parameters spend power that n = 6 does not have):
- member uptake caps (k·a·U);
- fitted maintenance or GAM;
- a growth rate per core (there is no SRT analogue in a sediment);
- trace-flux floors and effluent caps;
- prFBA weights (one reference per MAG gives p = 1 for every function);
- the staged or storage LPs;
- imposing any measured gas (circular).

**Worth considering:** these improve correctness or transparency, not the fit.
1. **Report CH4 as an FVA interval per core** (minimum and maximum at the protocol optimum) instead of a single vertex.
   - This matters for betaine: TMA is non-unique in 43 of 108 cells, and the total spans 2.5–4.5× measured.
2. **Betaine → TMA → methanogen cross-feeding** (a guild MSCommunity with a shared TMA pool) would be a stoichiometric sanity check, not an improvement.
   - Methylotrophic methanogenesis yields about 2.25 CH4 per TMA (4 TMA → 9 CH4 + 3 CO2 + 4 NH4+). The manuscript's 0.06 is the empirical community yield from the TMA-addition microcosms, so replacing it with a methanogen model would raise the betaine overshoot (already 2.6× measured) about 37-fold.
   - The only way the methanogen model could cap that is through its own fixed biomass (a qmax-type parameter). Simulated CH4 would then be a constant × *Methanolobus* biomass, whose ρ with measured CH4 is −0.43 (table above).
   - The *Methanolobus* model is also acetogenic as built and would need curation first (models/README.md, caveat 1).
3. **Element closure and ATP/electron-from-nothing gates** (both pipelines). The 2026 models already pass the free-ATP tests; closure rows add no parameters.
4. **Compare on all 12 cores (site gradient) as well as within the 6 unrestored ones.** Across sites, guild abundance and CH4 co-vary strongly (Figure 9), and any abundance-driven model reproduces that. It is fair to report, but it tests the site gradient rather than the mechanism.

**Recommendation:** do not port either pipeline to "regenerate" the measured methane. State plainly in the manuscript that:
- the models are a stoichiometric feasibility check (how much CH4 MePn and betaine could supply);
- between-core differences are not explained by any measured input at n = 6.

No community extension tested here is expected to improve between-core tracking. The most useful next steps are:
- checking the two 0.16 µg g⁻¹ DNA values;
- reporting FVA intervals.

The prfba_augmentation scouting flagged (2026-11-05) that three MePn-obligate MAGs produced 72 % of the simulated methane and that this should be checked before the totals are relied on. That flag referred to the 2020 models. **No 2026 model is MePn-obligate** (growth_triplet, 38/38), so it is resolved.

## 5. Can abundance-weighted methanogens, given their substrates, regenerate the measured methane? (2026-11-24)

Scripts: `methanogen_test.py` (→ `methanogen_test.json`, `methanogen_test_by_core.tsv`) and `bracket.py` (→ `bracket.tsv`).

### The manuscript's model magnitudes rest on an unstated growth rate of 1 h⁻¹

- **Where it comes from.** The areal conversion (`clauded_insights/modelling_consistent_basis_r232.py`, 2026-10-26, inherited by `simulations/consistent_basis_final.py`) multiplies two quantities:
  - k, the CH4 per g of biomass formed (flux ÷ bio1);
  - the standing biomass (g m⁻²).

  It then multiplies by 24 h d⁻¹ with no rate in between. That is dimensionally a daily flux only if every member replaces its biomass once an hour (μ = 1 h⁻¹, doubling in about 40 min).
- **The original pipeline** (bio1 fixed to abundance ÷ DNA) used a different, equally arbitrary rate.
- **Ranking is unaffected:** ρ, exact p and the ceiling test are all unchanged, because a uniform μ cancels.
- **Magnitudes are not:**
  - the MePn "14 % / 55 % of measured";
  - the Abstract's "one seventh to one half";
  - the betaine 2.6×, which uses the pipeline's conversion of each model's per-gDW TMA rate × the core's total biomass.

  All of these scale with this assumption.

### Specific rates each guild would need, at the pipeline's biomass constants (per core)

| Core | Measured CH4 (µmol m⁻² d⁻¹) | MePn μ needed, r = 0.07 (h⁻¹) | MePn μ needed, r = 0.33 | *Methanolobus* q needed (mmol CH4 g⁻¹ h⁻¹) | Methyl-substrate pool turnover needed (d) |
|---|---|---|---|---|---|
| R1_A | 737 | 3.0 | 0.79 | 2.3 | 5.3 |
| R1_B | 1347 | 5.2 | 1.4 | 21 | 2.2 |
| R1_C | 490 | 4.0 | 1.0 | 8.1 | 6.4 |
| R2_A | 1327 | 8.8 | 2.3 | 99 | 1.7 |
| R2_B | 1607 | 9.1 | 2.4 | 90 | 2.1 |
| R2_C | 1589 | 21 | 5.6 | 60 | 3.7 |

- **Methanogen yield:** in the methanogenic KOfam-arm builds of *Methanolobus* on methanol, 196–406 mmol CH4 per g biomass formed. This is about 5,000–9,000× the MePn consumers' 0.043.
- **Methanogen activity:** the builds' maximum specific rate is 50–64 mmol g⁻¹ h⁻¹, but that maximum is set by the 100-mmol uptake bound, not by physiology.
- **The 2026 build (arm A) makes no CH4.**
- **Substrate pools:** the standing porewater methanol and methylamines could supply the measured CH4 if turned over every 2–6 days (NMR concentrations × porosity 0.736 × methanogenesis stoichiometry). The betaine pool is 10–20× larger again.

### Bracket: simulated ÷ measured CH4, total of the six cores (`bracket.tsv`)

| Biomass factor f | μ (h⁻¹) | MePn r = 0.07 | MePn r = 0.33 | MePn, all P from MePn | Methanogen (arms C–D / B) |
|---|---|---|---|---|---|
| 1 | 0.001 | 0.0001 | 0.0006 | 0.002 | 0.014 / 0.029 |
| 1 | 0.01 | 0.001 | 0.006 | 0.022 | 0.14 / 0.29 |
| 1 | 0.05 | 0.007 | 0.028 | 0.11 | 0.7 / 1.4 |
| 10 | 0.01 | 0.014 | 0.055 | 0.22 | 1.4 / 2.9 |
| 10 | 0.05 | 0.072 | 0.28 | 1.1 | 7 / 14 |
| 1 | 1 (the manuscript's implicit rate) | 0.14 | 0.55 | 2.2 | 14 / 29 |

- **f** = 1 uses the pipeline's constants (DNA = 10 % of dry biomass, complete extraction). f = 10 corresponds to DNA at about 3 % of dry biomass and about 30 % extraction efficiency. These are illustrative values.
- **μ** of 0.001–0.05 h⁻¹ corresponds to sediment doubling times of about four weeks to 14 h.
- **MePn consumers** reach the measured flux only in the most favourable corner: MePn as essentially the sole P source, biomass about 10× the pipeline estimate, and doubling every 14 h. At the manuscript's ratios (0.07 and 0.33) they supply ≤ 7 % and ≤ 28 % at realistic growth.
- **The single methanogen species** reaches the measured flux under moderate assumptions:
  - doubling every 14 h at the pipeline's biomass;
  - or every 3 days at 10× biomass.

  The bracket counts growth-coupled CH4 only; maintenance methanogenesis would add to it.
- **Ranking does not improve.** Methanogen biomass vs measured CH4: ρ = −0.43 (p = 0.42). Weighted by the methyl-substrate CH4 potential: ρ = −0.60 (p = 0.24). The substrate potential alone: ρ = 0.26.

### Caveats and the check still open

- **Methanogen biomass may be underestimated.** It comes from one MAG species, 0.03–0.26 % of reads, in cores where 40–60 % of reads are unbinned.
- **That check is open.** The total *mcrA* abundance in the metagenomes would show whether unbinned methanogens add materially. The IMG KO table behind Figure 1 is not in the repository; the check needs that table from IMG, or an *mcrA* search of the unbinned contigs on poplar.
- **R2 cores:** the *Methanolobus* rates needed in R2_A and R2_B (90–99 mmol g⁻¹ h⁻¹) are high, so this matters most there.

### Proposed framing for the manuscript (not applied)

- Present model output as CH4 per unit of biomass formed.
- Present the comparison as the specific rate each guild would need to supply the measured flux.

This removes the hidden 1 h⁻¹. On that framing:
- growth-coupled MePn assimilation at the stated ratios cannot plausibly supply the measured methane;
- methylotrophic methanogenesis by the detected *Methanolobus* can, at moderate activity, from the standing methylated-substrate pools.

This bears directly on the manuscript's central claim and is for the authors to decide.
