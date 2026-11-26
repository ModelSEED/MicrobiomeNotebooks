# Numbers to reconcile — Figure 7 vs. the manuscript text

Every R² currently in the Results/Discussion/Abstract for the metabolic-modelling
section was computed on a **pre-DNA-weighting** analysis and matches **nothing** in
the current data or the submitted Figure 7. The submitted figure itself is on the
**DNA-double-normalised** basis (both axes divided by DNA density), which inflates
the correlations. The table gives, for each claim, the stale text value, the value
the submitted figure actually shows, and the honest corrected value (raw areal
units; Spearman ρ with p because n is tiny).

All corrected values are reproduced by `figure7_corrected.py` → `figure7_stats.json`.

## Panel (a) — modelled vs. measured methane  (n = 6 unrestored)

| Series | Text says | Submitted Fig 7 basis (DNA-norm) | **Corrected (raw areal units)** |
|---|---|---|---|
| MePn | R² = 0.978 (Disc.) / 0.966 (Res.) | R² = 0.89, ρ = 0.60 | **R² = 0.37, ρ = 0.60, p = 0.21 (NS)** |
| Betaine | R² = 0.862 | R² = 0.30, ρ = −0.43 | **R² = 0.10, ρ = −0.09, p = 0.87 (NS)** |

- The DNA normalisation alone moves MePn from R² = 0.37 → 0.89. The rank statistic
  (ρ = 0.60) is unchanged and, at n = 6, **not significant** (p = 0.21).
- On a consistent unit basis the model **under-predicts** measured emission (MePn
  points sit below the 1:1 line — the text already notes "predict half"), so the
  honest claim is "the model reproduces the *rank order / pattern*," not the magnitude.

## Panel (b) — methane vs. consumer abundance

| Series / set | Text says | Submitted Fig 7 (DNA-norm, n=6) | **Corrected raw, all 12 cores** |
|---|---|---|---|
| MePn consumers | R² = 0.974 (all) / 0.966 (unrestored) | R² = 0.56, ρ = 0.60 | **R² = 0.75, ρ = 0.73, p = 0.007** |
| Betaine consumers | R² = 0.852 (all) / 0.798 (unrestored) | R² = 0.19, ρ = 0.03 | **R² = 0.38, ρ = 0.70, p = 0.011** |

- The text's 0.974 / 0.966 / 0.852 / 0.798 reproduce **no** current computation — they are orphaned.
- **Key scientific correction:** with all 12 cores on raw units, betaine-consumer
  abundance (ρ = 0.70) tracks methane **about as well** as MePn-consumer abundance
  (ρ = 0.73), and *both* are significant. The "betaine barely correlates" story
  (ρ = 0.03) is an artefact of the n = 6 unrestored subset + DNA double-normalisation.
  MePn is still somewhat favoured (higher R²), but the separation is far narrower
  than the text claims.
- On the raw unrestored-only subset (n = 6), *neither* set is significant
  (MePn R² = 0.06, betaine R² = 0.17) — n = 6 simply cannot support these claims.

## Internal inconsistency to fix

- Discussion (Metabolic modelling): "MPn-consuming organisms were much more associated
  with methane emissions than betaine-consuming organisms (**R² = 0.798**)". The 0.798
  is the *old betaine* unrestored-abundance number — it is attached to the wrong
  organism set here.

## Other numbers (unchanged / still supported)

- "1.4 methane : 1 methylphosphonate" stoichiometry — recompute from the current
  fluxes and confirm; unaffected by normalisation.
- Ratio sensitivity (0.07 → R², 0.3 → R², unconstrained → R²): report these on the
  **corrected raw basis** and as a sensitivity curve, and lead with the **unconstrained**
  (parameter-free) fit so the result cannot be read as tuned to the answer.

## Recommended framing for the corrected text

1. Report ρ with p and n on every correlation; drop R² to 2 decimals (never 0.97828).
2. State the DNA-weighting/normalisation explicitly and show the raw-unit result as primary.
3. Soften "primary methanogenic substrate … MePn ≫ betaine" to reflect that betaine
   consumers also correlate with methane once all 12 cores are used; keep MePn as the
   *better-supported* (not exclusive) route, buttressed by the independent culturing
   (Fig 6) and substrate-addition (Fig 4) evidence, which do not depend on n = 6.
