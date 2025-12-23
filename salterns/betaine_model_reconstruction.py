#!/usr/bin/env python3
"""
Betaine Model Reconstruction and Gapfilling Script

This script reconstructs metabolic models for betaine-utilizing organisms
and gapfills the necessary pathways for betaine reduction via Stickland fermentation.

Organisms and their expected pathways:
- Halanaerobium (6 MAGs): Betaine reductase with H2 or amino acids as electron donors
- BM714 Clostridia (4 MAGs): Stickland fermentation with amino acids
- Acidaminobacteraceae (1 MAG): Stickland fermentation
- Synergistales (2 MAGs): Amino acid fermentation + betaine reduction
- Dethiosulfatibacteraceae (2 MAGs): Potentially Stickland-type
"""

import cobra
from cobra import Model, Reaction, Metabolite
import json
import os
from pathlib import Path

# Define betaine-associated MAGs with their taxonomic classifications
BETAINE_MAGS = {
    'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs': {
        'taxonomy': 'Clostridia; Tissierellales; Dethiosulfatibacteraceae; UBA8670',
        'sample': 'R1_A_D1',
        'pathway': 'stickland'
    },
    'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs': {
        'taxonomy': 'Halanaerobiia; Halanaerobiales',
        'sample': 'R1_A_D1',
        'pathway': 'betaine_reductase'
    },
    'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs': {
        'taxonomy': 'Clostridia; Peptostreptococcales; BM714',
        'sample': 'R1_A_D1',
        'pathway': 'stickland'
    },
    'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs': {
        'taxonomy': 'Clostridia; Peptostreptococcales; Acidaminobacteraceae',
        'sample': 'R1_A_D2',
        'pathway': 'stickland'
    },
    'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs': {
        'taxonomy': 'Halanaerobiia; Halanaerobiales; Halanaerobium saccharolyticum',
        'sample': 'R1_A_D2',
        'pathway': 'betaine_reductase'
    },
    'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs': {
        'taxonomy': 'Clostridia; Peptostreptococcales; BM714',
        'sample': 'R1_A_D2',
        'pathway': 'stickland'
    },
    'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs': {
        'taxonomy': 'Clostridia; Peptostreptococcales; BM714',
        'sample': 'R1_A_D2',
        'pathway': 'stickland'
    },
    'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs': {
        'taxonomy': 'Halanaerobiia; Halanaerobiales; Halanaerobium',
        'sample': 'R1_B_D1',
        'pathway': 'betaine_reductase'
    },
    'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs': {
        'taxonomy': 'Synergistia; Synergistales',
        'sample': 'R1_B_D1',
        'pathway': 'stickland'
    },
    'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs': {
        'taxonomy': 'Halanaerobiia; Halanaerobiales; Halanaerobium',
        'sample': 'R1_B_D2',
        'pathway': 'betaine_reductase'
    },
    'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs': {
        'taxonomy': 'Clostridia; Peptostreptococcales; BM714',
        'sample': 'R1_C_D1',
        'pathway': 'stickland'
    },
    'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs': {
        'taxonomy': 'Clostridia; Tissierellales; Dethiosulfatibacteraceae; UBA8670',
        'sample': 'R1_C_D2',
        'pathway': 'stickland'
    },
    'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs': {
        'taxonomy': 'Halanaerobiia; Halanaerobiales',
        'sample': 'R1_C_D2',
        'pathway': 'betaine_reductase'
    },
    'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs': {
        'taxonomy': 'Synergistia; Synergistales',
        'sample': 'R2_A_D1',
        'pathway': 'stickland'
    },
    'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs': {
        'taxonomy': 'Halanaerobiia; Halanaerobiales',
        'sample': 'R2_B_D2',
        'pathway': 'betaine_reductase'
    },
}

# ModelSEED compound definitions for betaine metabolism
METABOLITES = {
    # Cytoplasmic metabolites
    'cpd00540_c0': {'name': 'Glycine betaine', 'formula': 'C5H11NO2', 'charge': 0, 'compartment': 'c0'},
    'cpd00441_c0': {'name': 'Trimethylamine', 'formula': 'C3H9N', 'charge': 0, 'compartment': 'c0'},
    'cpd00196_c0': {'name': 'Acetylphosphate', 'formula': 'C2H3O5P', 'charge': -2, 'compartment': 'c0'},
    'cpd28060_c0': {'name': 'Reduced thioredoxin', 'formula': 'C10H15N4O4R4S2', 'charge': 0, 'compartment': 'c0'},
    'cpd27735_c0': {'name': 'Oxidized thioredoxin', 'formula': 'C10H13N4O4R4S2', 'charge': 0, 'compartment': 'c0'},
    'cpd00009_c0': {'name': 'Phosphate', 'formula': 'HO4P', 'charge': -2, 'compartment': 'c0'},
    'cpd00067_c0': {'name': 'H+', 'formula': 'H', 'charge': 1, 'compartment': 'c0'},
    'cpd00001_c0': {'name': 'H2O', 'formula': 'H2O', 'charge': 0, 'compartment': 'c0'},
    'cpd00005_c0': {'name': 'NADPH', 'formula': 'C21H26N7O17P3', 'charge': -4, 'compartment': 'c0'},
    'cpd00006_c0': {'name': 'NADP+', 'formula': 'C21H25N7O17P3', 'charge': -3, 'compartment': 'c0'},
    'cpd00029_c0': {'name': 'Acetate', 'formula': 'C2H3O2', 'charge': -1, 'compartment': 'c0'},
    'cpd00008_c0': {'name': 'ADP', 'formula': 'C10H12N5O10P2', 'charge': -3, 'compartment': 'c0'},
    'cpd00002_c0': {'name': 'ATP', 'formula': 'C10H12N5O13P3', 'charge': -4, 'compartment': 'c0'},
    # Extracellular metabolites
    'cpd00540_e0': {'name': 'Glycine betaine [e]', 'formula': 'C5H11NO2', 'charge': 0, 'compartment': 'e0'},
    'cpd00441_e0': {'name': 'Trimethylamine [e]', 'formula': 'C3H9N', 'charge': 0, 'compartment': 'e0'},
    'cpd00067_e0': {'name': 'H+ [e]', 'formula': 'H', 'charge': 1, 'compartment': 'e0'},
    'cpd00029_e0': {'name': 'Acetate [e]', 'formula': 'C2H3O2', 'charge': -1, 'compartment': 'e0'},
}

# Reaction definitions for betaine utilization pathways
REACTIONS = {
    # Betaine Reductase (EC 1.21.4.4) - rxn17220
    'rxn17220_c0': {
        'name': 'Betaine reductase',
        'metabolites': {
            'cpd00540_e0': -1,  # Glycine betaine [e] (substrate)
            'cpd28060_c0': -1,  # Reduced thioredoxin (electron donor)
            'cpd00009_c0': -1,  # Phosphate
            'cpd00067_c0': -2,  # 2 H+
            'cpd00441_e0': 1,   # Trimethylamine [e] (product)
            'cpd00196_c0': 1,   # Acetylphosphate
            'cpd27735_c0': 1,   # Oxidized thioredoxin
            'cpd00001_c0': 1,   # H2O
        },
        'lower_bound': 0,
        'upper_bound': 1000,
        'gene_reaction_rule': '',
        'subsystem': 'Betaine metabolism'
    },

    # Thioredoxin reductase (NADPH) - rxn27318
    'rxn27318_c0': {
        'name': 'Thioredoxin reductase (NADPH)',
        'metabolites': {
            'cpd00005_c0': -1,  # NADPH
            'cpd27735_c0': -1,  # Oxidized thioredoxin
            'cpd00006_c0': 1,   # NADP+
            'cpd28060_c0': 1,   # Reduced thioredoxin
        },
        'lower_bound': -1000,
        'upper_bound': 1000,
        'gene_reaction_rule': '',
        'subsystem': 'Thioredoxin system'
    },

    # Acetate kinase (Acetylphosphate -> Acetate + ATP) - rxn00225
    'rxn00225_c0': {
        'name': 'Acetate kinase',
        'metabolites': {
            'cpd00196_c0': -1,  # Acetylphosphate
            'cpd00008_c0': -1,  # ADP
            'cpd00029_c0': 1,   # Acetate
            'cpd00002_c0': 1,   # ATP
        },
        'lower_bound': -1000,
        'upper_bound': 1000,
        'gene_reaction_rule': '',
        'subsystem': 'Acetate metabolism'
    },

    # Betaine transport (proton symport) - rxn12375
    'rxn12375_c0': {
        'name': 'Betaine transport (H+ symport)',
        'metabolites': {
            'cpd00540_e0': -1,  # Glycine betaine [e]
            'cpd00067_e0': -1,  # H+ [e]
            'cpd00540_c0': 1,   # Glycine betaine [c]
            'cpd00067_c0': 1,   # H+ [c]
        },
        'lower_bound': 0,
        'upper_bound': 1000,
        'gene_reaction_rule': '',
        'subsystem': 'Transport'
    },

    # TMA transport (diffusion) - rxn09318
    'rxn09318_c0': {
        'name': 'Trimethylamine transport',
        'metabolites': {
            'cpd00441_c0': -1,  # TMA [c]
            'cpd00441_e0': 1,   # TMA [e]
        },
        'lower_bound': -1000,
        'upper_bound': 1000,
        'gene_reaction_rule': '',
        'subsystem': 'Transport'
    },

    # Acetate export
    'rxn05488_c0': {
        'name': 'Acetate transport',
        'metabolites': {
            'cpd00029_c0': -1,  # Acetate [c]
            'cpd00029_e0': 1,   # Acetate [e]
        },
        'lower_bound': -1000,
        'upper_bound': 1000,
        'gene_reaction_rule': '',
        'subsystem': 'Transport'
    },

    # Exchange reactions
    'EX_cpd00540_e0': {
        'name': 'Glycine betaine exchange',
        'metabolites': {'cpd00540_e0': -1},
        'lower_bound': -1000,
        'upper_bound': 1000,
        'gene_reaction_rule': '',
        'subsystem': 'Exchange'
    },
    'EX_cpd00441_e0': {
        'name': 'Trimethylamine exchange',
        'metabolites': {'cpd00441_e0': -1},
        'lower_bound': 0,
        'upper_bound': 1000,
        'gene_reaction_rule': '',
        'subsystem': 'Exchange'
    },
    'EX_cpd00029_e0': {
        'name': 'Acetate exchange',
        'metabolites': {'cpd00029_e0': -1},
        'lower_bound': 0,
        'upper_bound': 1000,
        'gene_reaction_rule': '',
        'subsystem': 'Exchange'
    },
}


def create_metabolite(met_id, met_data):
    """Create a cobra Metabolite object from definition."""
    met = Metabolite(
        id=met_id,
        name=met_data['name'],
        formula=met_data['formula'],
        charge=met_data.get('charge', 0),
        compartment=met_data['compartment']
    )
    return met


def create_reaction(rxn_id, rxn_data, metabolite_objects):
    """Create a cobra Reaction object from definition."""
    rxn = Reaction(rxn_id)
    rxn.name = rxn_data['name']
    rxn.lower_bound = rxn_data['lower_bound']
    rxn.upper_bound = rxn_data['upper_bound']
    rxn.gene_reaction_rule = rxn_data.get('gene_reaction_rule', '')
    rxn.subsystem = rxn_data.get('subsystem', '')

    # Add metabolites
    mets_to_add = {}
    for met_id, coeff in rxn_data['metabolites'].items():
        if met_id in metabolite_objects:
            mets_to_add[metabolite_objects[met_id]] = coeff
        else:
            # Create new metabolite if not exists
            if met_id in METABOLITES:
                met = create_metabolite(met_id, METABOLITES[met_id])
                metabolite_objects[met_id] = met
                mets_to_add[met] = coeff

    rxn.add_metabolites(mets_to_add)
    return rxn


def add_betaine_pathway_to_model(model, pathway_type='betaine_reductase'):
    """
    Add betaine utilization pathway to an existing model.

    Parameters:
    -----------
    model : cobra.Model
        The metabolic model to modify
    pathway_type : str
        Either 'betaine_reductase' or 'stickland'

    Returns:
    --------
    cobra.Model
        Modified model with betaine pathway
    """
    # Create metabolite objects dictionary
    metabolite_objects = {}

    # Check for existing metabolites in model and add missing ones
    for met_id, met_data in METABOLITES.items():
        if met_id in model.metabolites:
            metabolite_objects[met_id] = model.metabolites.get_by_id(met_id)
        else:
            met = create_metabolite(met_id, met_data)
            metabolite_objects[met_id] = met

    # Core betaine reductase reactions to add
    core_reactions = [
        'rxn17220_c0',   # Betaine reductase
        'rxn27318_c0',   # Thioredoxin reductase
        'rxn00225_c0',   # Acetate kinase
        'rxn09318_c0',   # TMA transport
        'rxn05488_c0',   # Acetate transport
        'EX_cpd00540_e0', # Betaine exchange
        'EX_cpd00441_e0', # TMA exchange
        'EX_cpd00029_e0', # Acetate exchange
    ]

    # Add reactions if not already present
    reactions_added = []
    for rxn_id in core_reactions:
        if rxn_id not in model.reactions:
            rxn = create_reaction(rxn_id, REACTIONS[rxn_id], metabolite_objects)
            model.add_reactions([rxn])
            reactions_added.append(rxn_id)

    return model, reactions_added


def extract_mag_from_community_model(community_model_path, mag_id):
    """
    Extract a single MAG model from a community model.

    Parameters:
    -----------
    community_model_path : str
        Path to the community model JSON file
    mag_id : str
        The MAG identifier to extract

    Returns:
    --------
    cobra.Model or None
        Extracted MAG model or None if extraction fails
    """
    with open(community_model_path, 'r') as f:
        comm_data = json.load(f)

    # Create new model for this MAG
    mag_model = Model(mag_id)

    # Get genes belonging to this MAG
    mag_prefix = mag_id.replace('.contigs', '').replace('_', '.')
    mag_genes = set()
    for gene in comm_data.get('genes', []):
        if mag_prefix in gene['id'] or mag_id.split('.')[0] in gene['id']:
            mag_genes.add(gene['id'])

    # If no genes found, return the full model subset
    # (this is a simplified extraction)
    return None  # Placeholder - full implementation would extract proper MAG


def load_or_create_mag_model(mag_id, mag_info, models_dir):
    """
    Load existing MAG model or create a minimal model for testing.

    Parameters:
    -----------
    mag_id : str
        MAG identifier
    mag_info : dict
        MAG metadata including sample and taxonomy
    models_dir : str
        Directory containing model files

    Returns:
    --------
    cobra.Model
        The MAG model
    """
    sample = mag_info['sample']
    community_model_path = os.path.join(models_dir, f"{sample}.json")

    # Try to load community model and extract MAG
    if os.path.exists(community_model_path):
        print(f"Loading community model from {community_model_path}")
        community_model = cobra.io.load_json_model(community_model_path)

        # Create a copy as the MAG model
        mag_model = community_model.copy()
        mag_model.id = mag_id
        return mag_model
    else:
        # Create minimal model if no community model exists
        print(f"Creating minimal model for {mag_id}")
        mag_model = Model(mag_id)
        mag_model.name = f"{mag_info['taxonomy']} - Betaine utilizer"
        return mag_model


def test_betaine_utilization(model):
    """
    Test if the model can utilize betaine as an electron acceptor.

    Parameters:
    -----------
    model : cobra.Model
        The metabolic model to test

    Returns:
    --------
    dict
        Test results including flux through betaine reductase
    """
    results = {
        'model_id': model.id,
        'has_betaine_reductase': 'rxn17220_c0' in model.reactions,
        'betaine_flux': 0,
        'tma_production': 0,
        'acetate_production': 0,
        'feasible': False
    }

    if 'rxn17220_c0' not in model.reactions:
        return results

    try:
        # Set up for betaine utilization test
        # Allow betaine uptake
        if 'EX_cpd00540_e0' in model.reactions:
            model.reactions.get_by_id('EX_cpd00540_e0').lower_bound = -10

        # Optimize for betaine reductase flux
        model.objective = 'rxn17220_c0'
        solution = model.optimize()

        if solution.status == 'optimal':
            results['feasible'] = True
            results['betaine_flux'] = solution.fluxes.get('rxn17220_c0', 0)
            results['tma_production'] = solution.fluxes.get('EX_cpd00441_e0', 0)
            results['acetate_production'] = solution.fluxes.get('EX_cpd00029_e0', 0)
    except Exception as e:
        results['error'] = str(e)

    return results


def main():
    """Main function to reconstruct and gapfill betaine-utilizing models."""

    # Set up paths
    base_dir = Path(__file__).parent
    models_dir = base_dir / 'models'
    output_dir = base_dir / 'betaine_models'
    output_dir.mkdir(exist_ok=True)

    print("=" * 70)
    print("Betaine Model Reconstruction and Gapfilling")
    print("=" * 70)
    print(f"\nProcessing {len(BETAINE_MAGS)} betaine-utilizing MAGs\n")

    results = []

    for mag_id, mag_info in BETAINE_MAGS.items():
        print(f"\n{'='*60}")
        print(f"Processing: {mag_id}")
        print(f"Taxonomy: {mag_info['taxonomy']}")
        print(f"Pathway type: {mag_info['pathway']}")
        print(f"Sample: {mag_info['sample']}")

        # Load or create model
        model = load_or_create_mag_model(mag_id, mag_info, str(models_dir))

        # Add betaine pathway
        model, reactions_added = add_betaine_pathway_to_model(
            model,
            pathway_type=mag_info['pathway']
        )
        print(f"Added {len(reactions_added)} reactions: {reactions_added}")

        # Test betaine utilization
        test_results = test_betaine_utilization(model)
        test_results['taxonomy'] = mag_info['taxonomy']
        test_results['pathway_type'] = mag_info['pathway']
        results.append(test_results)

        print(f"Betaine flux: {test_results['betaine_flux']:.2f}")
        print(f"TMA production: {test_results['tma_production']:.2f}")
        print(f"Feasible: {test_results['feasible']}")

        # Save model
        output_path = output_dir / f"{mag_id.replace('.contigs', '')}_betaine.json"
        cobra.io.save_json_model(model, str(output_path))
        print(f"Saved to: {output_path}")

    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\nTotal MAGs processed: {len(results)}")
    print(f"MAGs with feasible betaine utilization: {sum(1 for r in results if r['feasible'])}")

    # Group by pathway type
    for pathway in ['betaine_reductase', 'stickland']:
        pathway_results = [r for r in results if r['pathway_type'] == pathway]
        feasible = sum(1 for r in pathway_results if r['feasible'])
        print(f"\n{pathway}: {feasible}/{len(pathway_results)} feasible")

    # Save results summary
    summary_path = output_dir / 'reconstruction_summary.json'
    with open(summary_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {summary_path}")

    return results


if __name__ == '__main__':
    main()
