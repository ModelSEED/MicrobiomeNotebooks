# MePn / Bet pathway abundance by sample type — my values vs. Sheet2, and arbitration

Files: `read_abundance_by_reaction_totMAGs_byType.csv` (my %s by type),
`reaction_abundance_vs_sheet2.csv` (side-by-side). Sheet2 of `cliff_mags.xlsx` is a
prior manual attempt at the same split.

> Update: reaction presence is now taken from the **RAST gene annotation** (the local
> RASTified genome JSONs), not the ModelSEED models. This recovered rxn00270 (48 MAGs)
> and rxn05011 (41 MAGs) that the base models never built, confirmed C–P lyase in 22 MAGs
> (matching the curated table), and confirmed methyl-CoM reductase in ~0–1 MAGs (archaea gap).

## The two methods measure different things
| | denominator | weighting | presence source |
|---|---|---|---|
| **Mine** | whole mapped community | **read abundance** of each MAG | RAST gene annotation (per-MAG functions) |
| **Sheet2** | MePn+Bet gene pool only (MePn%+Bet% = 100%) | **gene copy counts** (unweighted) | direct RAST/KO gene annotation |

So Sheet2's "MePn% / Bet%" is a *relative gene-inventory split*, not a community-abundance %.

## By sample type
**Mine — % of mapped community carried by possessing MAGs (pathway union):**
| | Unrestored | Restored SF2 | Reference R2A |
|---|---|---|---|
| MePn | 8.4 | 0.2 | 1.0 |
| Bet  | 20.1 | 1.7 | 2.7 |

**Sheet2 — % of the MePn+Bet gene pool:**
| | Unrestored | Restored | Reference |
|---|---|---|---|
| MePn% | 19.4 | 0.2 | 3.6 |
| Bet%  | 80.6 | 99.8 | 96.4 |

Both agree qualitatively: everything is far higher in unrestored; MePn is the minority
of the combined MePn+Bet pool.

## Where they differ, and the arbitration
1. **Sheet2's "Bet" is inflated by non-methanogenic betaine genes.** Its `betaine` total sums
   betaine *transport* (Bet/Pro ABC, Opu, Ota) and *biosynthesis* (gbsA/B) — osmolyte usage,
   which does **not** make methane. The methanogenic betaine route is only `bet_reductase → tma_corr → mcom`.
   Comparing pathway *signatures* (C–P lyase vs betaine reductase), abundance-weighted:
   **MePn rxn43120 = 7.3 % ≥ Bet rxn17220 = 2.8 %** in unrestored — i.e. MePn ≥ betaine-reduction,
   the opposite of Sheet2's 19:81. Sheet2's split conflates osmolyte betaine with methanogenic betaine.
2. **Gene-count vs read-abundance weighting.** Sheet2 counts gene copies (a rare MAG counts
   the same as an abundant one). Mine weights by read abundance, so the abundant MePn carriers
   (Roseovarius) raise MePn's share (29.6 % of MePn+Bet vs Sheet2's 19.4 %). For "which pathway
   the abundant organisms run" (i.e. what drives methane flux), abundance-weighting is correct.
3. **Methanogenesis (mcom / rxn03127): Sheet2 is more reliable here.** My model-based count finds
   it in **1** MAG because the dereplicated MAG set has essentially no archaea (the binning gap
   flagged earlier); Sheet2's gene annotation finds mcom more broadly. Neither is large, but for
   the terminal methanogenesis step trust the annotation, not my models.
4. My `rxn24612` (tma_corr) is annotated in 58 MAGs — a broad/promiscuous corrinoid
   methyltransferase that inflates my **Bet union**. Prefer the signature `rxn17220`.

## Recommended ("correct") values
- **For "% abundance of MAGs possessing each reaction"** (the stated goal): use my by-type file —
  abundance-weighted, whole-community denominator. These are the correct *kind* of value.
- **For the MePn-vs-betaine methane comparison**: use the *signature* reactions, abundance-weighted
  (MePn C–P lyase `rxn43120` vs betaine reductase `rxn17220`): **MePn ≈ 7.3 % ≥ Bet ≈ 2.8 %** in
  unrestored, ~0 in restored/reference. Do **not** use Sheet2's 19:81 (osmolyte-inflated) or my
  Bet *union* (rxn24612-inflated).
- **Methanogenesis**: defer to Sheet2's gene annotation; my model count is an archaeal-binning artefact.

Net: the two didn't disagree by error so much as by definition. On the metric the question asks for
(abundance of the pathway-carrying organisms) my by-type values are correct; Sheet2's headline
Bet-dominance is an artefact of counting osmolyte-betaine genes and not weighting by abundance.
