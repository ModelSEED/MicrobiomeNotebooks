# Gapfill Integration Bug Report

**Notebook:** `correlational/test_bet.ipynb`
**Model:** `Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs__.RAST`
**Date:** 2026-03-23

---

## Problem

After gapfilling the model for biomass (`bio1`) production, the post-gapfilling model failed to grow biomass in the same media it was gapfilled in. The gapfill solver found a valid solution (141 new reactions), and `test_gapfill_database` confirmed the gapfill database model could grow, yet the integrated model could not.

---

## Root Cause

The notebook used a custom `_integrate_solution` function to add gapfilled reactions to the model. This function only added reactions found in the GramNeg template, silently discarding 4 reactions that were part of the gapfill solution but are not template reactions:

| Reaction | Type | Direction | Identity |
|----------|------|-----------|----------|
| `DM_cpd02701_c0` | Demand | `>` (export) | S-adenosyl-4-methylthio-2-oxobutanoate |
| `EX_cpd00048_e0` | Exchange | `<` (uptake) | Sulfate |
| `EX_cpd00063_e0` | Exchange | `<` (uptake) | Calcium |
| `EX_cpd00099_e0` | Exchange | `<` (uptake) | Chloride |

The filtering logic responsible:

```python
for rxn_id, d in gapfill_res['new'].items():
    if rxn_id[:-1] in template_gramneg.reactions:  # EX_ and DM_ reactions fail this check
        gap_sol[rxn_id[:-1]] = get_reaction_constraints_from_direction(d)
```

The `DM_cpd02701_c0` demand reaction was the critical missing piece. cpd02701 (S-adenosyl-4-methylthio-2-oxobutanoate) is a dead-end byproduct of methionine salvage metabolism that is produced as a side-effect of biomass precursor synthesis. Without a demand reaction to drain it, flux through those pathways is blocked, preventing any biomass production.

---

## Fix Applied

Replaced the custom `_integrate_solution` function with the built-in `MSGapfill.integrate_gapfill_solution()` method. This method:

1. Copies ALL reactions from the gapfill model (template reactions, exchange reactions, and demand reactions) into the original model
2. Sets reaction bounds according to the gapfilled direction
3. Assigns gene associations from reaction scores when available
4. Tests each reaction to determine if it is actually needed, removing unneeded ones
5. Validates that the integrated model achieves the minimum objective

### Before (broken)

```python
# Custom function that only handles template reactions
def _integrate_solution(template, model, gap_fill_solution):
    added_reactions = []
    for rxn_id, (lb, ub) in gap_fill_solution.items():
        template_reaction = template.reactions.get_by_id(rxn_id)
        model_reaction = template_reaction.to_reaction(model)
        model_reaction.lower_bound = lb
        model_reaction.upper_bound = ub
        added_reactions.append(model_reaction)
    model.add_reactions(added_reactions)
    add_exchanges = MSBuilder.add_exchanges_to_model(model)
    return added_reactions, add_exchanges

gap_sol = {}
for rxn_id, d in gapfill_res['new'].items():
    if rxn_id[:-1] in template_gramneg.reactions:
        gap_sol[rxn_id[:-1]] = get_reaction_constraints_from_direction(d)

new_model = model.copy()
added_reactions, add_exchanges = _integrate_solution(template_gramneg, new_model, gap_sol)
```

### After (fixed)

```python
gapfill.integrate_gapfill_solution(gapfill_res)
```

---

## Verification Results

Full end-to-end execution with `venv_official` (modelseedpy 0.4.3):

| Step | Result |
|------|--------|
| Model loaded | 834 reactions |
| ATP correction | 2 test conditions (empty, Glc) |
| Gapfilling | 141 new reactions found |
| Integration | 834 -> 970 reactions (135 added, 6 filtered as unneeded) |
| **Biomass flux** | **0.2351 (optimal)** |

The 6 reactions filtered as unneeded by the built-in optimizer: `rxn00152_c0`, `rxn01640_c0`, `rxn09427_c0`, `rxn00375_c0`, `rxn01013_c0`, `rxn01303_c0`.

One media compound (`cpd25960`) has no corresponding exchange reaction in the model, but this did not affect growth.

---

## Recommendations

1. **Always use `MSGapfill.integrate_gapfill_solution()`** rather than custom integration code. It handles all reaction types and includes built-in validation.
2. **Be aware that gapfill solutions can include non-template reactions** (EX\_, DM\_, SK\_ prefixes). Any custom integration logic must account for these.
3. **Demand reactions for dead-end metabolites** are a common and essential part of gapfill solutions, particularly for biomass precursor pathways involving methionine salvage, cofactor biosynthesis, and lipid metabolism.
