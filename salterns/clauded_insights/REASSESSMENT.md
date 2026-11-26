# Re-assessment — what the earlier corrections got right, and the deeper issue they missed

A second, more skeptical pass (new tests in `reassessment_figure.py`) changes the
conclusion in one important way and corrects one of my earlier claims.

## Correction to my own earlier analysis

- **"DNA double-normalisation inflates the correlation" was over-stated.** It is
  subset-dependent, not a general artefact: for panel (b) with all 12 cores the *raw*
  fit (R² = 0.75) is actually **higher** than the DNA-normalised fit (0.64); the
  inflation only appears in the n = 6 subset. And for panel (a) the two axes have
  *opposite* DNA dependence, so the "shared-divisor spurious correlation" mechanism I
  invoked does not apply there. The real problem is **n and an environmental-gradient
  confound**, not normalisation per se.
- **The corrected Figure 7b (all 12 cores, ρ = 0.73/0.70, both significant) is more
  honest but still measures the wrong thing** — see below. My correction improved a
  fundamentally weak analysis rather than replacing it.

## The finding that reframes everything: one environmental gradient drives methane

Core-level correlations with CH₄ (n = 12, depth-averaged; `reassessment_figure.py`):

| Variable | Spearman ρ | p |
|---|---|---|
| **Sulfate (SO4-S)** | **+0.97** | <0.001 |
| Chloride | +0.94 | <0.001 |
| Salinity | +0.79 | 0.002 |
| **N:P ratio** | **+0.73** | 0.007 |
| MePn-consumer abundance | +0.73 | 0.007 |
| Betaine-consumer abundance | +0.70 | 0.011 |
| Olsen-P (phosphate) | −0.29 | 0.35 |

Methane is set by the **salinity / sulfate / hypersalinity axis** (unrestored salterns
high, restored/reference low) far more strongly than by any biological variable. Every
cross-sample quantity — metabolites, MAG abundances, modelled fluxes — sorts along this
one axis, so **any** cross-sample correlation with methane is largely a gradient echo.

## Two tests neither the manuscript nor the earlier pass ran

1. **Specificity (Fig `reassessment_gradient_confound.png` b).** The summed abundance of
   a *random* 20-MAG subset correlates with methane about as well as the MePn consumers:
   **37% of random subsets reach ρ ≥ 0.73** (betaine: 45%). → the consumer→methane
   correlation (Fig 7b / S8) is **not consumer-specific** and cannot, alone, support a
   MePn- or betaine-specific mechanism.
2. **Biomass confound (partial correlation).** Controlling for DNA density (ρ(DNA,CH₄) =
   0.52), the consumer signal weakens to ρ = 0.59 (MePn, p = 0.04) / 0.57 (betaine,
   p = 0.06). Part of the correlation is just "more biomass → more methane."

## What the evidence actually supports (honest hierarchy)

1. **Strongest — manipulative, gradient-free:** substrate additions (Fig 4: MePn-acid →
   more CH₄) and Halomonas culturing (Fig 6). These are the paper's real backbone.
2. **Moderate — mechanistic but gradient-confounded:** measured-metabolite correlations
   (`...c`): **TMA ρ = 0.86**, methylamine 0.72, betaine 0.70, DMA −0.61. Interpretable
   because these are actual substrates; TMA (the direct methanogen substrate) is the
   standout, and it is the *amine/betaine* route.
3. **Weakest — non-specific, gradient-confounded:** the metabolic-modelling and
   consumer-abundance correlations (Fig 7a/b, S8) — exactly the part this repo produces.

## Implications for the manuscript

- **De-emphasise the abundance correlation** (Fig 7b/S8) as evidence, or report the
  random-subset null so readers see it is not specific.
- **Recast the modelling** (Fig 7a) as *mechanistic sufficiency + attribution* ("MePn can
  produce enough methane; Roseovarius (77%)/Albimonas (11%) carry it" — family_contribution.png),
  not as validation-by-correlation.
- **Add a geochemistry–methane figure.** Sulfate ρ = 0.97, salinity 0.79, N:P 0.73 are
  strong, directly relevant (N:P supports the P-limitation → MePn story), and currently
  under-shown. Use partial correlation / variance partitioning to ask what explains
  methane *beyond* the salinity/sulfate gradient (probably little at the sample level —
  which is why the incubations matter).
- **Lead with the causal evidence** (Fig 4, 6) and the TMA metabolite result; treat
  cross-sample correlations as consistent-but-confounded context.

## Scope note

This directory is the **metabolic-modelling half** only. The substrate-addition (Fig 4,
Table 2), 16S (Fig S6), and IMG-KO gene-abundance (Fig 1) data/code are **not here** —
they live in the JGI/R analysis (Zenodo). The NMR concentrations and geochemistry *are*
here (in `data/Cliff_Sample_Metadata_BGC_NMR.csv`) and are under-analysed; the tests
above were run directly from that file.
