#!/usr/bin/env python3
"""
Gapfill Betaine Reductase Pathway in Community Models

This script adds betaine reductase (rxn17220) to all MAG compartments in
community models, enabling betaine utilization as an electron acceptor
in Stickland-type fermentation.

ModelSEED Reaction: rxn17220 (EC 1.21.4.4)
Equation: Betaine + Red-Thioredoxin + Pi + 2H+ → TMA + Acetylphosphate + Ox-Thioredoxin + H2O

Author: Generated for saltern metagenome analysis
"""

import cobra
from cobra import Reaction, Metabolite
import json
from pathlib import Path
from datetime import datetime

# Define betaine MAGs with their sample assignments
BETAINE_MAGS = {
    'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs': 'R1_A_D1',
    'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs': 'R1_A_D1',
    'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs': 'R1_A_D1',
    'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs': 'R1_A_D2',
    'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs': 'R1_A_D2',
    'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs': 'R1_A_D2',
    'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs': 'R1_A_D2',
    'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs': 'R1_B_D1',
    'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs': 'R1_B_D1',
    'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs': 'R1_B_D2',
    'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs': 'R1_C_D1',
    'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs': 'R1_C_D2',
    'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs': 'R1_C_D2',
    'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs': 'R2_A_D1',
    'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs': 'R2_B_D2',
}

# Taxonomic information
TAXONOMY = {
    'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs': 'Clostridia; Tissierellales; Dethiosulfatibacteraceae',
    'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs': 'Halanaerobiia; Halanaerobiales',
    'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs': 'Clostridia; Peptostreptococcales; BM714',
    'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs': 'Clostridia; Peptostreptococcales; Acidaminobacteraceae',
    'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs': 'Halanaerobiia; Halanaerobiales; Halanaerobium saccharolyticum',
    'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs': 'Clostridia; Peptostreptococcales; BM714',
    'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs': 'Clostridia; Peptostreptococcales; BM714',
    'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs': 'Halanaerobiia; Halanaerobiales; Halanaerobium',
    'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs': 'Synergistia; Synergistales',
    'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs': 'Halanaerobiia; Halanaerobiales; Halanaerobium',
    'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs': 'Clostridia; Peptostreptococcales; BM714',
    'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs': 'Clostridia; Tissierellales; Dethiosulfatibacteraceae',
    'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs': 'Halanaerobiia; Halanaerobiales',
    'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs': 'Synergistia; Synergistales',
    'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs': 'Halanaerobiia; Halanaerobiales',
}


def add_betaine_reductase_to_compartment(model, comp):
    """
    Add betaine reductase reaction to a specific compartment.

    Parameters:
    -----------
    model : cobra.Model
        The metabolic model
    comp : str
        Compartment ID (e.g., 'c1', 'c2')

    Returns:
    --------
    bool
        True if reaction was added, False otherwise
    """
    rxn_id = f"rxn17220_{comp}"

    # Skip if already exists
    if rxn_id in model.reactions:
        return False

    try:
        # Get metabolites - using shared extracellular compartment e0
        betaine_ext = model.metabolites.get_by_id('cpd00540_e0')
        tma_ext = model.metabolites.get_by_id('cpd00441_e0')
        thio_red = model.metabolites.get_by_id(f'cpd28060_{comp}')
        thio_ox = model.metabolites.get_by_id(f'cpd27735_{comp}')
        phosphate = model.metabolites.get_by_id(f'cpd00009_{comp}')
        proton = model.metabolites.get_by_id(f'cpd00067_{comp}')
        acetylp = model.metabolites.get_by_id(f'cpd00196_{comp}')
        water = model.metabolites.get_by_id(f'cpd00001_{comp}')

        # Create betaine reductase reaction (EC 1.21.4.4)
        rxn = Reaction(rxn_id)
        rxn.name = f'Betaine reductase [{comp}]'
        rxn.lower_bound = 0
        rxn.upper_bound = 1000
        rxn.subsystem = 'Betaine metabolism'

        # Add metabolites with stoichiometry
        # rxn17220: Betaine + Red-Thioredoxin + Pi + 2H+ → TMA + Acetylphosphate + Ox-Thioredoxin + H2O
        rxn.add_metabolites({
            betaine_ext: -1,      # Glycine betaine [e0]
            thio_red: -1,         # Reduced thioredoxin
            phosphate: -1,        # Phosphate
            proton: -2,           # 2 H+
            tma_ext: 1,           # Trimethylamine [e0]
            acetylp: 1,           # Acetylphosphate
            thio_ox: 1,           # Oxidized thioredoxin
            water: 1,             # H2O
        })

        model.add_reactions([rxn])
        return True

    except KeyError as e:
        return False


def ensure_exchange_reactions(model):
    """
    Ensure exchange reactions exist for betaine and TMA.

    Parameters:
    -----------
    model : cobra.Model
        The metabolic model

    Returns:
    --------
    list
        List of exchange reactions added
    """
    added = []

    # Betaine exchange
    if 'EX_cpd00540_e0' not in model.reactions:
        if 'cpd00540_e0' in model.metabolites:
            betaine_ext = model.metabolites.get_by_id('cpd00540_e0')
            ex_bet = Reaction('EX_cpd00540_e0')
            ex_bet.name = 'Glycine betaine exchange'
            ex_bet.add_metabolites({betaine_ext: -1})
            ex_bet.lower_bound = -1000
            ex_bet.upper_bound = 1000
            model.add_reactions([ex_bet])
            added.append('EX_cpd00540_e0')

    # TMA exchange
    if 'EX_cpd00441_e0' not in model.reactions:
        if 'cpd00441_e0' in model.metabolites:
            tma_ext = model.metabolites.get_by_id('cpd00441_e0')
            ex_tma = Reaction('EX_cpd00441_e0')
            ex_tma.name = 'Trimethylamine exchange'
            ex_tma.add_metabolites({tma_ext: -1})
            ex_tma.lower_bound = 0
            ex_tma.upper_bound = 1000
            model.add_reactions([ex_tma])
            added.append('EX_cpd00441_e0')

    return added


def find_compartments_with_thioredoxin(model):
    """
    Find all compartments that have thioredoxin system.

    Parameters:
    -----------
    model : cobra.Model
        The metabolic model

    Returns:
    --------
    list
        List of compartment IDs with thioredoxin
    """
    compartments = []
    for comp_num in range(0, 100):
        comp = f"c{comp_num}"
        thio_red_id = f"cpd28060_{comp}"
        if thio_red_id in model.metabolites:
            compartments.append(comp)
    return compartments


def test_betaine_utilization(model, compartments=None):
    """
    Test betaine utilization capacity.

    Parameters:
    -----------
    model : cobra.Model
        The metabolic model
    compartments : list, optional
        List of compartments to test

    Returns:
    --------
    dict
        Test results
    """
    results = {
        'betaine_uptake': 0,
        'tma_production': 0,
        'active_compartments': [],
        'feasible': False
    }

    if 'EX_cpd00540_e0' not in model.reactions:
        return results

    # Allow betaine uptake
    model.reactions.get_by_id('EX_cpd00540_e0').lower_bound = -10

    try:
        # Maximize TMA production
        if 'EX_cpd00441_e0' in model.reactions:
            model.objective = 'EX_cpd00441_e0'
            solution = model.optimize()

            if solution.status == 'optimal':
                results['feasible'] = True
                results['tma_production'] = solution.objective_value
                results['betaine_uptake'] = -solution.fluxes.get('EX_cpd00540_e0', 0)

                # Find active betaine reductase reactions
                for rxn in model.reactions:
                    if 'rxn17220' in rxn.id and solution.fluxes[rxn.id] > 0.01:
                        results['active_compartments'].append({
                            'compartment': rxn.id.split('_')[-1],
                            'flux': solution.fluxes[rxn.id]
                        })
    except Exception as e:
        results['error'] = str(e)

    return results


def process_community_model(model_path, output_dir):
    """
    Process a single community model.

    Parameters:
    -----------
    model_path : Path
        Path to the community model
    output_dir : Path
        Output directory for gapfilled model

    Returns:
    --------
    dict
        Processing results
    """
    sample_id = model_path.stem
    print(f"\n{'='*60}")
    print(f"Processing: {sample_id}")

    # Load model
    model = cobra.io.load_json_model(str(model_path))
    initial_rxns = len(model.reactions)

    # Find compartments with thioredoxin
    compartments = find_compartments_with_thioredoxin(model)
    print(f"  Found {len(compartments)} compartments with thioredoxin")

    # Add betaine reductase to each compartment
    added_count = 0
    for comp in compartments:
        if add_betaine_reductase_to_compartment(model, comp):
            added_count += 1

    print(f"  Added {added_count} betaine reductase reactions")

    # Ensure exchange reactions
    exchanges_added = ensure_exchange_reactions(model)
    if exchanges_added:
        print(f"  Added exchange reactions: {exchanges_added}")

    # Test betaine utilization
    test_results = test_betaine_utilization(model, compartments)
    print(f"  Betaine uptake: {test_results['betaine_uptake']:.2f}")
    print(f"  TMA production: {test_results['tma_production']:.2f}")
    print(f"  Active compartments: {len(test_results['active_compartments'])}")

    # Save gapfilled model
    output_path = output_dir / f"{sample_id}_betaine_gapfilled.json"
    cobra.io.save_json_model(model, str(output_path))
    print(f"  Saved to: {output_path}")

    return {
        'sample_id': sample_id,
        'initial_reactions': initial_rxns,
        'final_reactions': len(model.reactions),
        'betaine_reactions_added': added_count,
        'compartments': len(compartments),
        **test_results
    }


def main():
    """Main function to gapfill all community models."""

    base_dir = Path(__file__).parent
    models_dir = base_dir / 'models'
    output_dir = base_dir / 'betaine_models'
    output_dir.mkdir(exist_ok=True)

    print("="*70)
    print("Betaine Reductase Gapfilling for Community Models")
    print("="*70)
    print(f"\nTimestamp: {datetime.now().isoformat()}")
    print(f"ModelSEED Reaction: rxn17220 (EC 1.21.4.4)")
    print("Equation: BET + Red-Thioredoxin + Pi + 2H+ → TMA + AcP + Ox-Thioredoxin + H2O")

    # Get unique samples from betaine MAGs
    samples = set(BETAINE_MAGS.values())
    print(f"\nSamples to process: {sorted(samples)}")

    # Process each community model
    results = []
    for sample in sorted(samples):
        model_path = models_dir / f"{sample}.json"
        if model_path.exists():
            result = process_community_model(model_path, output_dir)
            results.append(result)
        else:
            print(f"\nWarning: Model not found: {model_path}")

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"\nTotal models processed: {len(results)}")
    print(f"Models with feasible betaine utilization: {sum(1 for r in results if r['feasible'])}")

    total_rxns_added = sum(r['betaine_reactions_added'] for r in results)
    print(f"Total betaine reductase reactions added: {total_rxns_added}")

    # MAG-level summary
    print("\n" + "-"*60)
    print("Betaine MAG Summary:")
    for mag_id, sample in sorted(BETAINE_MAGS.items(), key=lambda x: x[1]):
        taxonomy = TAXONOMY.get(mag_id, 'Unknown')
        short_id = mag_id.split('.')[-2] + '.' + mag_id.split('.')[-1].replace('.contigs', '')
        print(f"  {sample}: {short_id}")
        print(f"    Taxonomy: {taxonomy}")

    # Save results
    results_path = output_dir / 'gapfilling_results.json'
    with open(results_path, 'w') as f:
        json.dump({
            'timestamp': datetime.now().isoformat(),
            'betaine_mags': BETAINE_MAGS,
            'taxonomy': TAXONOMY,
            'results': results
        }, f, indent=2)
    print(f"\nResults saved to: {results_path}")

    return results


if __name__ == '__main__':
    main()
