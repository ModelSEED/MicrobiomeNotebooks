#!/usr/bin/env python3
"""
Gap-fill betaine metabolism pathway into saltern MAG models
Based on literature evidence for betaine reductase (GrdHI) pathway
"""

import cobra
from cobra import Reaction, Metabolite
import json
from copy import deepcopy
import pandas as pd

print("="*80)
print("BETAINE METABOLISM PATHWAY GAP-FILLING")
print("="*80)
print()

# Load the test model
model_id = 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs'
model_path = f'./models/{model_id}__.RAST.json'

print(f"Loading model: {model_id}")
print(f"  Taxonomy: Halanaerobiales > Halanaerobiaceae > Halanaerobium saccharolyticum")
print(f"  Literature confidence: HIGH")
print()

model = cobra.io.load_json_model(model_path)

print(f"Original model stats:")
print(f"  Reactions: {len(model.reactions)}")
print(f"  Metabolites: {len(model.metabolites)}")
print(f"  Genes: {len(model.genes)}")
print()

# Check what metabolites already exist
print("Checking for required metabolites...")
required_metabolites = {
    'cpd00540_c0': 'Glycine betaine (cytosol)',
    'cpd00441_c0': 'Trimethylamine (cytosol)',
    'cpd00441_e0': 'Trimethylamine (extracellular)',
    'cpd00029_c0': 'Acetate (cytosol)',
    'cpd00067_c0': 'Proton (cytosol)',
}

missing_metabolites = []
for met_id, met_name in required_metabolites.items():
    if met_id not in model.metabolites:
        missing_metabolites.append((met_id, met_name))
        print(f"  ✗ Missing: {met_id} ({met_name})")
    else:
        print(f"  ✓ Present: {met_id} ({met_name})")

print()

# Add missing metabolites
print("Adding missing metabolites...")
for met_id, met_name in missing_metabolites:
    compartment = met_id.split('_')[-1]
    base_id = met_id.replace(f'_{compartment}', '')

    new_met = Metabolite(
        id=met_id,
        name=met_name,
        compartment=compartment
    )
    print(f"  + Added: {met_id} ({met_name})")

print()

# Define betaine reductase reaction
# Biochemistry: Betaine + 2 Reduced_Ferredoxin → Trimethylamine + Acetate + 2 Oxidized_Ferredoxin
print("Adding betaine reductase reaction (rxn17220_c0)...")

betaine_reductase = Reaction('rxn17220_c0')
betaine_reductase.name = 'Glycine betaine reductase (GrdHI)'
betaine_reductase.subsystem = 'Betaine metabolism'
betaine_reductase.lower_bound = 0  # Irreversible
betaine_reductase.upper_bound = 1000

# Get or create metabolites
if 'cpd00540_c0' not in model.metabolites:
    betaine_c = Metabolite('cpd00540_c0', name='Glycine betaine', compartment='c0')
else:
    betaine_c = model.metabolites.cpd00540_c0

if 'cpd00441_c0' not in model.metabolites:
    tma_c = Metabolite('cpd00441_c0', name='Trimethylamine', compartment='c0')
else:
    tma_c = model.metabolites.cpd00441_c0

if 'cpd00029_c0' not in model.metabolites:
    acetate_c = Metabolite('cpd00029_c0', name='Acetate', compartment='c0')
else:
    acetate_c = model.metabolites.cpd00029_c0

if 'cpd00067_c0' not in model.metabolites:
    h_c = Metabolite('cpd00067_c0', name='H+', compartment='c0')
else:
    h_c = model.metabolites.cpd00067_c0

# Simplified reaction without explicit ferredoxin (assume it's coupled to other processes)
# Betaine + H+ → Trimethylamine + Acetate
betaine_reductase.add_metabolites({
    betaine_c: -1.0,
    h_c: -1.0,
    tma_c: 1.0,
    acetate_c: 1.0
})

model.add_reactions([betaine_reductase])
print(f"  ✓ Added: {betaine_reductase.id}")
print(f"    Reaction: {betaine_reductase.reaction}")
print()

# Add TMA transport (cytosol to extracellular)
print("Adding TMA transport reaction...")

if 'cpd00441_e0' not in model.metabolites:
    tma_e = Metabolite('cpd00441_e0', name='Trimethylamine', compartment='e0')
else:
    tma_e = model.metabolites.cpd00441_e0

tma_transport = Reaction('rxn00002_TMA')
tma_transport.name = 'Trimethylamine transport'
tma_transport.subsystem = 'Transport'
tma_transport.lower_bound = -1000  # Reversible
tma_transport.upper_bound = 1000

tma_transport.add_metabolites({
    tma_c: -1.0,
    tma_e: 1.0
})

model.add_reactions([tma_transport])
print(f"  ✓ Added: {tma_transport.id}")
print(f"    Reaction: {tma_transport.reaction}")
print()

# Add or verify TMA exchange reaction
print("Adding/verifying TMA exchange reaction...")

if 'EX_cpd00441_e0' not in model.reactions:
    tma_exchange = Reaction('EX_cpd00441_e0')
    tma_exchange.name = 'Trimethylamine exchange'
    tma_exchange.subsystem = 'Exchange'
    tma_exchange.lower_bound = -1000
    tma_exchange.upper_bound = 1000

    tma_exchange.add_metabolites({
        tma_e: -1.0
    })

    model.add_reactions([tma_exchange])
    print(f"  ✓ Added: {tma_exchange.id}")
    print(f"    Reaction: {tma_exchange.reaction}")
else:
    print(f"  ✓ Already present: EX_cpd00441_e0")

print()

# Verify betaine exchange exists
print("Verifying betaine exchange reaction...")
if 'EX_cpd00540_e0' in model.reactions:
    print(f"  ✓ Present: EX_cpd00540_e0")
    print(f"    Bounds: [{model.reactions.EX_cpd00540_e0.lower_bound}, {model.reactions.EX_cpd00540_e0.upper_bound}]")
else:
    print(f"  ✗ Missing: EX_cpd00540_e0 - This is required for betaine uptake!")

print()
print(f"Gap-filled model stats:")
print(f"  Reactions: {len(model.reactions)} (+{len(model.reactions) - 878})")
print(f"  Metabolites: {len(model.metabolites)} (+{len(model.metabolites) - 906})")
print()

# Save gap-filled model
output_path = f'./models/{model_id}__.RAST.gapfilled.json'
cobra.io.save_json_model(model, output_path)
print(f"✓ Saved gap-filled model: {output_path}")
print()

# Test the gap-filled model
print("="*80)
print("TESTING GAP-FILLED MODEL")
print("="*80)
print()

# Define betaine minimal media
media_betaine_minimal = {
    'EX_cpd00540_e0': 10,     # Betaine (changed to uptake value, not bounds)
    'EX_cpd00001_e0': 1000,   # Water
    'EX_cpd00067_e0': 1000,   # H+
    'EX_cpd00099_e0': 1000,   # Cl-
    'EX_cpd00205_e0': 1000,   # K+
    'EX_cpd00009_e0': 1000,   # Phosphate
    'EX_cpd00048_e0': 1000,   # Sulfate
    'EX_cpd00034_e0': 10,     # Zn2+
    'EX_cpd00058_e0': 10,     # Cu2+
    'EX_cpd00063_e0': 10,     # Ca2+
    'EX_cpd00149_e0': 10,     # Co2+
    'EX_cpd00244_e0': 10,     # Ni2+
    'EX_cpd00030_e0': 10,     # Mn2+
    'EX_cpd10516_e0': 10,     # Fe3+
}

# Filter to only include exchange reactions that exist in the model
test_media = {}
for rxn_id, uptake_rate in media_betaine_minimal.items():
    if rxn_id in model.reactions:
        test_media[rxn_id] = uptake_rate

print("Test 1: Betaine minimal media")
print(f"  Media components: {len(test_media)}")

model.medium = test_media

# Set objective
if 'bio1' in model.reactions:
    model.objective = 'bio1'
    obj_name = 'bio1'
else:
    model.objective = 'ATPM_c0'
    obj_name = 'ATPM_c0'

print(f"  Objective: {obj_name}")

try:
    solution = model.optimize()

    print(f"\nResults:")
    print(f"  Status: {solution.status}")
    print(f"  Growth rate: {solution.objective_value:.6f}")

    if solution.status == 'optimal':
        # Check key fluxes
        if 'EX_cpd00540_e0' in model.reactions:
            betaine_flux = solution.fluxes.get('EX_cpd00540_e0', 0.0)
            print(f"  Betaine uptake: {betaine_flux:.6f}")

        if 'rxn17220_c0' in model.reactions:
            grd_flux = solution.fluxes.get('rxn17220_c0', 0.0)
            print(f"  Betaine reductase flux: {grd_flux:.6f}")

        if 'EX_cpd00441_e0' in model.reactions:
            tma_flux = solution.fluxes.get('EX_cpd00441_e0', 0.0)
            print(f"  TMA production: {tma_flux:.6f}")

        if 'EX_cpd00029_e0' in model.reactions:
            acetate_flux = solution.fluxes.get('EX_cpd00029_e0', 0.0)
            print(f"  Acetate production: {acetate_flux:.6f}")

        print()

        if solution.objective_value > 0.001:
            print("✓ SUCCESS: Model shows growth on betaine media!")
        else:
            print("⚠ WARNING: Model is feasible but shows minimal growth")
            print("   May need additional nutrients or biomass precursors")
    else:
        print(f"✗ FAILED: Optimization status is {solution.status}")

except Exception as e:
    print(f"✗ ERROR: {e}")

print()

# Test 2: Compare to original model (without gap-filling)
print("="*80)
print("COMPARISON: Original vs Gap-filled")
print("="*80)
print()

original_model = cobra.io.load_json_model(model_path)
original_model.medium = test_media
if 'bio1' in original_model.reactions:
    original_model.objective = 'bio1'
else:
    original_model.objective = 'ATPM_c0'

try:
    original_solution = original_model.optimize()
    print(f"Original model:")
    print(f"  Growth rate: {original_solution.objective_value:.6f}")
    print(f"  Has betaine reductase: {'rxn17220_c0' in original_model.reactions}")
    print(f"  Has TMA exchange: {'EX_cpd00441_e0' in original_model.reactions}")
    print()
except:
    print(f"Original model: Failed to optimize")
    print()

print(f"Gap-filled model:")
print(f"  Growth rate: {solution.objective_value:.6f}")
print(f"  Has betaine reductase: {'rxn17220_c0' in model.reactions}")
print(f"  Has TMA exchange: {'EX_cpd00441_e0' in model.reactions}")
print()

growth_improvement = solution.objective_value - original_solution.objective_value
if growth_improvement > 0.001:
    print(f"✓ Growth improvement: +{growth_improvement:.6f} ({100*growth_improvement/max(original_solution.objective_value, 0.001):.1f}% increase)")
elif solution.objective_value > 0.001:
    print(f"✓ Gap-filling enabled betaine metabolism (original had no growth)")
else:
    print(f"⚠ No growth improvement - additional gap-filling may be needed")

print()
print("="*80)
print("GAP-FILLING COMPLETE")
print("="*80)
print()
print(f"Files generated:")
print(f"  1. {output_path} - Gap-filled model")
print()
print(f"Reactions added:")
print(f"  1. rxn17220_c0 - Glycine betaine reductase")
print(f"  2. rxn00002_TMA - TMA transport")
print(f"  3. EX_cpd00441_e0 - TMA exchange (if missing)")
print()
print(f"Next steps:")
print(f"  - If growth is still zero, may need to add:")
print(f"    * Additional electron carriers (ferredoxin)")
print(f"    * Cofactor biosynthesis pathways")
print(f"    * Essential biomass precursors")
print(f"  - Scale gap-filling to all 11 HIGH confidence models")
print(f"  - Validate against experimental flux measurements")
