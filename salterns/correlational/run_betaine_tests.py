#!/usr/bin/env python3
"""
Standalone script to test betaine metabolism in saltern microbiome models
Based on literature review and hypothesis validation framework
"""

import cobra
from modelseedpy import MSBuilder
import json
from tqdm import tqdm
import pandas as pd
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import os

print("="*80)
print("BETAINE METABOLISM MODEL VALIDATION")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# Load taxonomy
print("Loading taxonomy data...")
with open("betaine_reducers_taxonomy.json", "r") as f:
    betaine_reducers_taxonomy = json.load(f)

# Define model IDs
target_model_ids = list(betaine_reducers_taxonomy.keys())
print(f"Found {len(target_model_ids)} models to test\n")

# Load models
print("Loading models...")
models_base = {}
for model_id in tqdm(target_model_ids, desc="Loading models"):
    try:
        model_path = f'./models/{model_id}__.RAST.json'
        if os.path.exists(model_path):
            model = cobra.io.load_json_model(model_path)
            # Add ATPM reaction if not present
            if 'ATPM_c0' not in model.reactions:
                try:
                    from modelseedpy.core.msmodelutl import MSModelUtil
                    MSModelUtil.add_atp_hydrolysis(model, 1.0)
                except:
                    # Fallback: manually add ATPM reaction if needed
                    pass
            models_base[model_id] = model
        else:
            print(f"  Warning: Model file not found: {model_path}")
    except Exception as e:
        print(f"  Error loading {model_id}: {e}")

print(f"\nSuccessfully loaded {len(models_base)} models\n")

if len(models_base) == 0:
    print("ERROR: No models loaded. Check that models directory exists.")
    exit(1)

# Organize models by taxonomy
print("Organizing models by taxonomy...")
model_groups = {
    'Halanaerobiales': [],
    'Peptostreptococcales': [],
    'Tissierellales': [],
    'Synergistales': []
}

for model_id, taxa in betaine_reducers_taxonomy.items():
    if model_id not in models_base:
        continue
    order = taxa.get('Order', '')
    if 'Halanaerobiales' in order:
        model_groups['Halanaerobiales'].append(model_id)
    elif 'Peptostreptococcales' in order:
        model_groups['Peptostreptococcales'].append(model_id)
    elif 'Tissierellales' in order:
        model_groups['Tissierellales'].append(model_id)
    elif 'Synergistales' in order:
        model_groups['Synergistales'].append(model_id)

print("\nModel Organization:")
for group, models in model_groups.items():
    print(f"  {group}: {len(models)} models")

# Define media conditions
print("\nDefining media conditions...")

media_betaine_minimal = {
    'EX_cpd00540_e0': (-10, 1000),
    'EX_cpd00001_e0': (-1000, 1000),
    'EX_cpd00067_e0': (-1000, 1000),
    'EX_cpd00099_e0': (-1000, 1000),
    'EX_cpd00205_e0': (-1000, 1000),
    'EX_cpd00009_e0': (-1000, 1000),
    'EX_cpd00034_e0': (-10, 1000),
    'EX_cpd00058_e0': (-10, 1000),
    'EX_cpd00063_e0': (-10, 1000),
    'EX_cpd00149_e0': (-10, 1000),
    'EX_cpd00244_e0': (-10, 1000),
    'EX_cpd00030_e0': (-10, 1000),
    'EX_cpd10516_e0': (-10, 1000),
    'EX_cpd00048_e0': (-1000, 1000),
}

media_betaine_selenite = media_betaine_minimal.copy()
media_betaine_selenite['EX_cpd00116_e0'] = (-10, 1000)

media_betaine_stickland = media_betaine_selenite.copy()
media_betaine_stickland.update({
    'EX_cpd00035_e0': (-5, 1000),
    'EX_cpd00107_e0': (-5, 1000),
    'EX_cpd00156_e0': (-5, 1000),
    'EX_cpd00119_e0': (-5, 1000),
})

media_amino_acids = {
    'EX_cpd00001_e0': (-1000, 1000),
    'EX_cpd00067_e0': (-1000, 1000),
    'EX_cpd00099_e0': (-1000, 1000),
    'EX_cpd00009_e0': (-1000, 1000),
    'EX_cpd00035_e0': (-10, 1000),
    'EX_cpd00039_e0': (-10, 1000),
    'EX_cpd00041_e0': (-10, 1000),
    'EX_cpd00051_e0': (-10, 1000),
    'EX_cpd00054_e0': (-10, 1000),
    'EX_cpd00060_e0': (-10, 1000),
    'EX_cpd00065_e0': (-10, 1000),
    'EX_cpd00066_e0': (-10, 1000),
    'EX_cpd00069_e0': (-10, 1000),
    'EX_cpd00107_e0': (-10, 1000),
    'EX_cpd00119_e0': (-10, 1000),
    'EX_cpd00129_e0': (-10, 1000),
    'EX_cpd00156_e0': (-10, 1000),
    'EX_cpd00048_e0': (-10, 1000),
}

print("Media conditions defined:")
print(f"  1. Betaine minimal: {len(media_betaine_minimal)} components")
print(f"  2. Betaine + Selenite: {len(media_betaine_selenite)} components")
print(f"  3. Betaine + Amino acids (Stickland): {len(media_betaine_stickland)} components")
print(f"  4. Amino acids only: {len(media_amino_acids)} components")

# Testing function
def test_betaine_metabolism(model, model_id, media_dict, test_name, check_selenite=False):
    """Test betaine metabolism under specified conditions"""
    test_model = deepcopy(model)

    # Filter media and convert tuple format to single values
    test_media = {}
    for rxn_id, bounds in media_dict.items():
        if rxn_id in test_model.reactions:
            if isinstance(bounds, tuple):
                # Use absolute value of lower bound as uptake rate
                test_media[rxn_id] = abs(bounds[0])
            else:
                test_media[rxn_id] = abs(bounds)

    test_model.medium = test_media

    # Set objective
    if 'bio1' in test_model.reactions:
        test_model.objective = 'bio1'
        obj_name = 'bio1'
    else:
        test_model.objective = 'ATPM_c0'
        obj_name = 'ATPM_c0'

    results = {
        'model_id': model_id,
        'test_name': test_name,
        'objective': obj_name,
        'growth_rate': 0.0,
        'betaine_uptake': 0.0,
        'tma_production': 0.0,
        'acetate_production': 0.0,
        'amino_acid_uptake': {},
        'has_betaine_reductase': 'rxs17220_c0' in test_model.reactions,
        'has_tma_transport': 'rxs00002_c0' in test_model.reactions or 'EX_cpd00441_e0' in test_model.reactions,
        'selenite_available': 'EX_cpd00116_e0' in test_media if check_selenite else 'N/A',
        'feasible': False,
        'error': None
    }

    try:
        solution = test_model.optimize()

        if solution.status == 'optimal':
            results['feasible'] = True
            results['growth_rate'] = solution.objective_value

            if 'EX_cpd00540_e0' in test_model.reactions:
                results['betaine_uptake'] = solution.fluxes.get('EX_cpd00540_e0', 0.0)

            tma_exchanges = ['EX_cpd00441_e0', 'rxs00002_c0']
            for ex in tma_exchanges:
                if ex in test_model.reactions:
                    flux = solution.fluxes.get(ex, 0.0)
                    if abs(flux) > 0.001:
                        results['tma_production'] = flux
                        break

            if 'EX_cpd00029_e0' in test_model.reactions:
                results['acetate_production'] = solution.fluxes.get('EX_cpd00029_e0', 0.0)

            aa_exchanges = ['EX_cpd00035_e0', 'EX_cpd00107_e0', 'EX_cpd00156_e0', 'EX_cpd00119_e0']
            for aa_ex in aa_exchanges:
                if aa_ex in test_model.reactions:
                    flux = solution.fluxes.get(aa_ex, 0.0)
                    if abs(flux) > 0.001:
                        results['amino_acid_uptake'][aa_ex] = flux
        else:
            results['error'] = f"Optimization status: {solution.status}"

    except Exception as e:
        results['error'] = str(e)

    return results

# Run tests
print("\n" + "="*80)
print("RUNNING COMPREHENSIVE TESTS")
print("="*80 + "\n")

all_test_results = []
test_configs = [
    ('Betaine minimal', media_betaine_minimal, False),
    ('Betaine + Selenite', media_betaine_selenite, True),
    ('Betaine + Amino acids (Stickland)', media_betaine_stickland, True),
    ('Amino acids only', media_amino_acids, False)
]

for model_id in tqdm(models_base.keys(), desc="Testing models"):
    model = models_base[model_id]
    taxa = betaine_reducers_taxonomy[model_id]
    order = taxa.get('Order', 'Unknown')

    print(f"\n{'='*80}")
    print(f"Model: {model_id[:60]}")
    print(f"Order: {order} | Family: {taxa.get('Family', 'Unknown')}")
    print(f"{'='*80}")

    for test_name, media, check_sel in test_configs:
        print(f"  Testing: {test_name}...")
        result = test_betaine_metabolism(model, model_id, media, test_name, check_sel)
        all_test_results.append(result)

        if result['feasible']:
            print(f"    ✓ Growth: {result['growth_rate']:.4f}")
            print(f"      Betaine uptake: {result['betaine_uptake']:.4f}")
            print(f"      TMA production: {result['tma_production']:.4f}")
            print(f"      Acetate production: {result['acetate_production']:.4f}")
            if result['amino_acid_uptake']:
                print(f"      Amino acid uptake: {len(result['amino_acid_uptake'])} types")
        else:
            print(f"    ✗ No growth - {result['error']}")

print(f"\n{'='*80}")
print(f"Testing complete! Total results: {len(all_test_results)}")
print(f"{'='*80}\n")

# Create results DataFrame
print("Analyzing results...")
results_df = pd.DataFrame(all_test_results)

results_df['Order'] = results_df['model_id'].apply(
    lambda x: betaine_reducers_taxonomy[x].get('Order', 'Unknown')
)
results_df['Family'] = results_df['model_id'].apply(
    lambda x: betaine_reducers_taxonomy[x].get('Family', 'Unknown')
)
results_df['Genus'] = results_df['model_id'].apply(
    lambda x: betaine_reducers_taxonomy[x].get('Genus', 'Unknown')
)
results_df['model_short'] = results_df['model_id'].apply(
    lambda x: x.split('_')[-1].replace('.contigs', '')
)

print(f"\nResults Summary:")
print(f"Total tests: {len(results_df)}")
print(f"Feasible tests: {results_df['feasible'].sum()} ({100*results_df['feasible'].sum()/len(results_df):.1f}%)")
print(f"\nBy taxonomic order:")
print(results_df.groupby('Order')['feasible'].agg(['sum', 'count', lambda x: f"{100*x.sum()/len(x):.1f}%"]))

# Save results
print("\nSaving results...")
results_df.to_csv('betaine_test_results.csv', index=False)
print("  ✓ CSV results saved to: betaine_test_results.csv")

# Generate simple text summary
with open('betaine_test_summary.txt', 'w') as f:
    f.write("="*80 + "\n")
    f.write("BETAINE METABOLISM MODEL VALIDATION - SUMMARY\n")
    f.write("="*80 + "\n")
    f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    f.write(f"Total Models Tested: {len(models_base)}\n")
    f.write(f"Total Tests Run: {len(results_df)}\n")
    f.write(f"Feasible Tests: {results_df['feasible'].sum()} ({100*results_df['feasible'].sum()/len(results_df):.1f}%)\n\n")

    f.write("="*80 + "\n")
    f.write("RESULTS BY TAXONOMIC ORDER\n")
    f.write("="*80 + "\n\n")

    for order in sorted(results_df['Order'].unique()):
        order_df = results_df[results_df['Order'] == order]
        feasible_count = order_df['feasible'].sum()
        total_count = len(order_df)

        f.write(f"{order}\n")
        f.write("-" * 80 + "\n")
        f.write(f"  Models: {len(order_df['model_id'].unique())}\n")
        f.write(f"  Tests: {total_count} ({feasible_count} feasible, {100*feasible_count/total_count:.1f}%)\n")

        feas_df = order_df[order_df['feasible']]
        if len(feas_df) > 0:
            f.write(f"  Average growth rate: {feas_df['growth_rate'].mean():.4f}\n")
            f.write(f"  Average betaine uptake: {feas_df['betaine_uptake'].abs().mean():.4f}\n")
            f.write(f"  Average TMA production: {feas_df['tma_production'].abs().mean():.4f}\n")
        f.write("\n")

print("  ✓ Text summary saved to: betaine_test_summary.txt")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)
print(f"\nGenerated files:")
print(f"  1. betaine_test_results.csv - Full results data")
print(f"  2. betaine_test_summary.txt - Quick summary")
print(f"\nTo view results:")
print(f"  cat betaine_test_summary.txt")
print(f"  open betaine_test_results.csv")
print("="*80)
