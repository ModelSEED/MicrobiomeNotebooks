#!/usr/bin/env python3
"""
Diagnose why gap-filled model isn't growing on betaine media
"""

import cobra
from cobra.flux_analysis import find_blocked_reactions, find_essential_reactions

print("="*80)
print("DIAGNOSING GROWTH BOTTLENECK IN GAP-FILLED MODEL")
print("="*80)
print()

# Load gap-filled model
model_path = './models/Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs__.RAST.gapfilled.json'
print(f"Loading gap-filled model...")
model = cobra.io.load_json_model(model_path)

# Set up betaine media
media_betaine_minimal = {
    'EX_cpd00540_e0': 10,
    'EX_cpd00001_e0': 1000,
    'EX_cpd00067_e0': 1000,
    'EX_cpd00099_e0': 1000,
    'EX_cpd00205_e0': 1000,
    'EX_cpd00009_e0': 1000,
    'EX_cpd00048_e0': 1000,
    'EX_cpd00034_e0': 10,
    'EX_cpd00058_e0': 10,
    'EX_cpd00063_e0': 10,
    'EX_cpd00149_e0': 10,
    'EX_cpd00244_e0': 10,
    'EX_cpd00030_e0': 10,
    'EX_cpd10516_e0': 10,
}

test_media = {k: v for k, v in media_betaine_minimal.items() if k in model.reactions}
model.medium = test_media

if 'bio1' in model.reactions:
    model.objective = 'bio1'
else:
    model.objective = 'ATPM_c0'

print(f"Model: {model.id}")
print(f"Objective: {model.objective}")
print(f"Media components: {len(test_media)}")
print()

# Try to optimize
solution = model.optimize()
print(f"Initial optimization:")
print(f"  Status: {solution.status}")
print(f"  Objective value: {solution.objective_value}")
print()

if solution.objective_value < 0.001:
    print("Zero growth detected. Running diagnostics...")
    print()

    # Check 1: Can the model grow on ANY media?
    print("Check 1: Testing with complete media...")
    complete_media = {}
    for rxn in model.exchanges:
        if rxn.id.startswith('EX_') and rxn.lower_bound < 0:
            complete_media[rxn.id] = abs(rxn.lower_bound)

    model.medium = complete_media
    complete_solution = model.optimize()
    print(f"  Growth on complete media: {complete_solution.objective_value:.6f}")

    if complete_solution.objective_value > 0.001:
        print(f"  ✓ Model CAN grow - betaine media is limiting")
    else:
        print(f"  ✗ Model CANNOT grow - fundamental issue with biomass reaction")
    print()

    # Reset to betaine media
    model.medium = test_media

    # Check 2: What reactions are blocked?
    print("Check 2: Finding blocked reactions...")
    try:
        blocked = find_blocked_reactions(model)
        print(f"  Total blocked reactions: {len(blocked)}")

        # Check if betaine pathway is blocked
        betaine_rxns = ['rxn17220_c0', 'rxn00002_TMA', 'EX_cpd00441_e0', 'EX_cpd00540_e0']
        blocked_betaine = [r for r in betaine_rxns if r in blocked]
        if blocked_betaine:
            print(f"  ⚠ Betaine pathway reactions blocked: {blocked_betaine}")
        else:
            print(f"  ✓ Betaine pathway reactions NOT blocked")
    except Exception as e:
        print(f"  Error: {e}")
    print()

    # Check 3: Biomass reaction requirements
    print("Check 3: Analyzing biomass reaction...")
    if 'bio1' in model.reactions:
        biomass_rxn = model.reactions.bio1
        print(f"  Biomass reaction: {biomass_rxn.id}")
        print(f"  Metabolites required ({len(biomass_rxn.metabolites)}):")

        # Group by compartment
        compartments = {}
        for met, coef in biomass_rxn.metabolites.items():
            if coef < 0:  # Required (consumed)
                comp = met.compartment
                if comp not in compartments:
                    compartments[comp] = []
                compartments[comp].append((met.id, met.name, coef))

        for comp, mets in sorted(compartments.items()):
            print(f"\n    Compartment {comp}: {len(mets)} precursors")
            for met_id, met_name, coef in mets[:10]:  # Show first 10
                print(f"      {met_id}: {coef:.4f}")
            if len(mets) > 10:
                print(f"      ... and {len(mets)-10} more")
    print()

    # Check 4: Can betaine be taken up?
    print("Check 4: Testing betaine uptake in isolation...")
    test_model = model.copy()

    # Set betaine exchange as objective
    test_model.objective = 'EX_cpd00540_e0'
    test_model.objective_direction = 'min'  # Minimize = maximize uptake

    uptake_solution = test_model.optimize()
    print(f"  Max betaine uptake: {abs(uptake_solution.objective_value):.6f}")

    if abs(uptake_solution.objective_value) > 0.001:
        print(f"  ✓ Betaine CAN be taken up")
    else:
        print(f"  ✗ Betaine CANNOT be taken up - transport issue")
    print()

    # Check 5: Can betaine reductase carry flux?
    print("Check 5: Testing betaine reductase flux capacity...")
    test_model2 = model.copy()

    # Set betaine reductase as objective
    if 'rxn17220_c0' in test_model2.reactions:
        test_model2.objective = 'rxn17220_c0'
        test_model2.objective_direction = 'max'

        reductase_solution = test_model2.optimize()
        print(f"  Max reductase flux: {reductase_solution.objective_value:.6f}")

        if reductase_solution.objective_value > 0.001:
            print(f"  ✓ Betaine reductase CAN carry flux")

            # What's consumed and produced?
            print(f"\n  Key fluxes at max reductase activity:")
            for rxn_id in ['EX_cpd00540_e0', 'rxn17220_c0', 'EX_cpd00441_e0', 'EX_cpd00029_e0']:
                if rxn_id in test_model2.reactions:
                    flux = reductase_solution.fluxes.get(rxn_id, 0.0)
                    print(f"    {rxn_id}: {flux:.6f}")
        else:
            print(f"  ✗ Betaine reductase CANNOT carry flux - missing cofactors/electron carriers?")
    print()

    # Check 6: Try relaxing biomass constraint
    print("Check 6: Testing with relaxed biomass requirements...")
    test_model3 = model.copy()

    # Allow model to produce incomplete biomass
    if 'bio1' in test_model3.reactions:
        biomass_rxn = test_model3.reactions.bio1
        # Scale down all biomass coefficients
        original_coefs = dict(biomass_rxn.metabolites)

        for met in biomass_rxn.metabolites:
            coef = biomass_rxn.metabolites[met]
            if coef < 0:  # Required metabolite
                # Reduce requirement by 50%
                biomass_rxn.add_metabolites({met: -coef * 0.5})

        relaxed_solution = test_model3.optimize()
        print(f"  Growth with 50% reduced biomass requirements: {relaxed_solution.objective_value:.6f}")

        if relaxed_solution.objective_value > 0.001:
            print(f"  ⚠ Growth possible with relaxed biomass - missing biosynthesis pathways")
        else:
            print(f"  ✗ Still no growth - not a biomass precursor issue")

        # Restore original coefficients
        for met, coef in original_coefs.items():
            biomass_rxn.add_metabolites({met: -biomass_rxn.metabolites[met] + coef})
    print()

    # Check 7: Energy balance
    print("Check 7: Checking energy metabolism...")
    atp_producing = []
    for rxn in model.reactions:
        if 'cpd00002_c0' in [m.id for m in rxn.metabolites]:  # ATP
            atp_coef = rxn.metabolites.get(model.metabolites.cpd00002_c0, 0)
            if atp_coef > 0:  # Produces ATP
                atp_producing.append(rxn.id)

    print(f"  ATP-producing reactions: {len(atp_producing)}")

    # Check if any ATP production is possible
    if 'ATPM_c0' in model.reactions:
        test_model4 = model.copy()
        test_model4.objective = 'ATPM_c0'
        atp_solution = test_model4.optimize()
        print(f"  Max ATP production rate: {atp_solution.objective_value:.6f}")

        if atp_solution.objective_value > 0.001:
            print(f"  ✓ ATP production is possible")
        else:
            print(f"  ✗ Cannot produce ATP - energy metabolism broken")
    print()

print("="*80)
print("DIAGNOSIS COMPLETE")
print("="*80)
print()

print("Recommended actions based on diagnostics:")
print("  1. Check if model can grow on complete media")
print("  2. If betaine pathway is blocked, add missing electron carriers")
print("  3. If biomass precursors missing, add biosynthesis pathways")
print("  4. If energy metabolism broken, check ATP/cofactor pathways")
print("  5. Consider using ModelSEEDpy's automated gap-filling")
