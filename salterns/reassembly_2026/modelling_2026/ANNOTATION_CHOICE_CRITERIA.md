# Annotation source for the 2026 salterns metabolic models: pre-registered comparison

Fixed on 2026-11-20, before any of the compared models was built. Andrew asked to "consider other annotation tools besides RAST if the resulting models would be sufficiently higher quality".

## Arms

All arms use the same pipeline on the 57 modelling species representatives (`../comparison_rollout/modelling_species.tsv`):
- the same combined curated templates (`templates/`, newest upstream plus all of Andrew's fork/PR curations);
- the same ModelSEEDpy (Andrew's fork, branch `sync-cshenry-main`);
- the same ATP correction, gap-filling media and minimum biomass.

| Arm | Gene → reaction evidence |
|---|---|
| **A: RAST** (baseline) | RAST functional roles (SSO), the templates' native vocabulary |
| **B: RAST + KOfam** | A, plus KOfam KOs (2026-10-26, DRAM rule; already computed) translated to ModelSEED reactions through the ModelSEED KO→reaction mapping |
| **C: RAST + KOfam + eggNOG** | B, plus eggNOG-mapper KOs, EC numbers and BiGG reactions translated to ModelSEED reactions |
| **D: RAST + KOfam + eggNOG + Bakta** | C, plus Bakta EC numbers and products |

Other annotations are attached to RAST features by coordinate overlap (≥ 80 % reciprocal overlap, same strand). An unmatched gene enters as an extra feature carrying only its own annotation.

## Metrics, per model

1. Gene-associated (non-gap-filled) reactions.
2. Gap-filled reactions needed to grow:
   - on the lactate + MePn, phosphate-free medium for MePn consumers;
   - on the lactate + betaine medium for betaine reducers;
   - with the existing simulation constraints.
3. Guild pathway support:
   - C–P lyase (*phnGHIJKLM*) for MePn consumers;
   - betaine reductase (*grdH/grdI*) for betaine reducers;
   - methanogenesis core for the methanogen.

   Counted as the fraction of the pathway's reactions that are gene-associated rather than gap-filled or hand-added.
4. Energy sanity: no ATP made from nothing (maximise ATPM on a medium without carbon or electron donor = 0), and the ATP-correction test pass rate.
5. Blocked and dead-end reaction fraction after gap-filling.
6. Spurious gene support: reactions whose only evidence is a single non-RAST source with weak support (EC-only, or eggNOG without KO), reported separately.

## Decision rule

Adopt the arm with the most added evidence that satisfies **all** of the following against arm A, on the median over the 57 models:
- (a) gene-associated reactions ≥ 10 % higher;
- (b) gap-filled reactions ≥ 20 % lower;
- (c) guild-pathway gene support not lower in any guild;
- (d) no model gains an ATP-from-nothing loop, and the ATP-correction pass rate is not lower;
- (e) the blocked/dead-end reaction fraction rises by no more than 2 percentage points.

If no merged arm meets (a)–(e), use arm A (RAST). Report every metric for every arm, whichever wins.
