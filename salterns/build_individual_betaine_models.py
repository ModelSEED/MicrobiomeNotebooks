#!/usr/bin/env python3
"""
Build Individual Betaine-Utilizing MAG Models

Creates metabolic models for each of the 15 betaine-utilizing MAGs with
their literature-supported pathways:

- Halanaerobium (6 MAGs): Betaine reductase (EC 1.21.4.4)
- Clostridia BM714 (4 MAGs): Stickland fermentation with betaine
- Acidaminobacteraceae (1 MAG): Stickland fermentation
- Synergistales (2 MAGs): Amino acid fermentation + betaine reduction
- Dethiosulfatibacteraceae (2 MAGs): Stickland-type fermentation

Key reactions from ModelSEED:
- rxn17220: Betaine reductase
- rxn27318: Thioredoxin reductase
- rxn00225: Acetate kinase
"""

import cobra
from cobra import Model, Reaction, Metabolite
import json
from pathlib import Path
from datetime import datetime

# Define betaine MAGs with taxonomy and pathway type
BETAINE_MAGS = {
    'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs': {
        'sample': 'R1_A_D1',
        'taxonomy': 'Clostridia; Tissierellales; Dethiosulfatibacteraceae; UBA8670',
        'pathway_type': 'stickland',
        'class': 'Dethiosulfatibacteraceae'
    },
    'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs': {
        'sample': 'R1_A_D1',
        'taxonomy': 'Halanaerobiia; Halanaerobiales',
        'pathway_type': 'betaine_reductase',
        'class': 'Halanaerobium'
    },
    'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs': {
        'sample': 'R1_A_D1',
        'taxonomy': 'Clostridia; Peptostreptococcales; BM714',
        'pathway_type': 'stickland',
        'class': 'BM714'
    },
    'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs': {
        'sample': 'R1_A_D2',
        'taxonomy': 'Clostridia; Peptostreptococcales; Acidaminobacteraceae',
        'pathway_type': 'stickland',
        'class': 'Acidaminobacteraceae'
    },
    'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs': {
        'sample': 'R1_A_D2',
        'taxonomy': 'Halanaerobiia; Halanaerobiales; Halanaerobium saccharolyticum',
        'pathway_type': 'betaine_reductase',
        'class': 'Halanaerobium'
    },
    'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs': {
        'sample': 'R1_A_D2',
        'taxonomy': 'Clostridia; Peptostreptococcales; BM714',
        'pathway_type': 'stickland',
        'class': 'BM714'
    },
    'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs': {
        'sample': 'R1_A_D2',
        'taxonomy': 'Clostridia; Peptostreptococcales; BM714',
        'pathway_type': 'stickland',
        'class': 'BM714'
    },
    'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs': {
        'sample': 'R1_B_D1',
        'taxonomy': 'Halanaerobiia; Halanaerobiales; Halanaerobium',
        'pathway_type': 'betaine_reductase',
        'class': 'Halanaerobium'
    },
    'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs': {
        'sample': 'R1_B_D1',
        'taxonomy': 'Synergistia; Synergistales',
        'pathway_type': 'stickland',
        'class': 'Synergistales'
    },
    'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs': {
        'sample': 'R1_B_D2',
        'taxonomy': 'Halanaerobiia; Halanaerobiales; Halanaerobium',
        'pathway_type': 'betaine_reductase',
        'class': 'Halanaerobium'
    },
    'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs': {
        'sample': 'R1_C_D1',
        'taxonomy': 'Clostridia; Peptostreptococcales; BM714',
        'pathway_type': 'stickland',
        'class': 'BM714'
    },
    'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs': {
        'sample': 'R1_C_D2',
        'taxonomy': 'Clostridia; Tissierellales; Dethiosulfatibacteraceae; UBA8670',
        'pathway_type': 'stickland',
        'class': 'Dethiosulfatibacteraceae'
    },
    'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs': {
        'sample': 'R1_C_D2',
        'taxonomy': 'Halanaerobiia; Halanaerobiales',
        'pathway_type': 'betaine_reductase',
        'class': 'Halanaerobium'
    },
    'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs': {
        'sample': 'R2_A_D1',
        'taxonomy': 'Synergistia; Synergistales',
        'pathway_type': 'stickland',
        'class': 'Synergistales'
    },
    'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs': {
        'sample': 'R2_B_D2',
        'taxonomy': 'Halanaerobiia; Halanaerobiales',
        'pathway_type': 'betaine_reductase',
        'class': 'Halanaerobium'
    },
}


def create_metabolite(model, cpd_id, name, formula=None, compartment='c0'):
    """Create and add a metabolite to the model."""
    met_id = f"{cpd_id}_{compartment}"
    if met_id in model.metabolites:
        return model.metabolites.get_by_id(met_id)

    met = Metabolite(
        met_id,
        name=name,
        formula=formula,
        compartment=compartment
    )
    model.add_metabolites([met])
    return met


def add_core_metabolism(model):
    """Add essential core metabolism reactions."""

    # Create compartments
    model.compartments = {'c0': 'cytosol', 'e0': 'extracellular'}

    # Core metabolites needed for betaine pathway
    core_mets = {
        # Cytosolic metabolites
        'cpd00001': ('H2O', 'H2O', 'c0'),
        'cpd00002': ('ATP', 'C10H12N5O13P3', 'c0'),
        'cpd00008': ('ADP', 'C10H12N5O10P2', 'c0'),
        'cpd00009': ('Phosphate', 'HO4P', 'c0'),
        'cpd00067': ('H+', 'H', 'c0'),
        'cpd00005': ('NADPH', 'C21H26N7O17P3', 'c0'),
        'cpd00006': ('NADP+', 'C21H25N7O17P3', 'c0'),
        'cpd00003': ('NAD+', 'C21H26N7O14P2', 'c0'),
        'cpd00004': ('NADH', 'C21H27N7O14P2', 'c0'),
        'cpd00010': ('CoA', 'C21H32N7O16P3S', 'c0'),
        'cpd00022': ('Acetyl-CoA', 'C23H34N7O17P3S', 'c0'),
        'cpd00029': ('Acetate', 'C2H3O2', 'c0'),
        'cpd00196': ('Acetylphosphate', 'C2H3O5P', 'c0'),
        'cpd28060': ('Reduced thioredoxin', None, 'c0'),
        'cpd27735': ('Oxidized thioredoxin', None, 'c0'),
        # Extracellular metabolites
        'cpd00540': ('Glycine betaine', 'C5H11NO2', 'e0'),
        'cpd00441': ('Trimethylamine', 'C3H9N', 'e0'),
        'cpd00029_e': ('Acetate', 'C2H3O2', 'e0'),
        # Amino acids for Stickland (cytosolic)
        'cpd00035': ('L-Alanine', 'C3H7NO2', 'c0'),
        'cpd00033': ('Glycine', 'C2H5NO2', 'c0'),
        'cpd00156': ('L-Valine', 'C5H11NO2', 'c0'),
        'cpd00322': ('L-Isoleucine', 'C6H13NO2', 'c0'),
        'cpd00107': ('L-Leucine', 'C6H13NO2', 'c0'),
        'cpd00013': ('NH3', 'H3N', 'c0'),
        'cpd00011': ('CO2', 'CO2', 'c0'),
    }

    for cpd_id, (name, formula, compartment) in core_mets.items():
        if cpd_id.endswith('_e'):
            actual_id = cpd_id.replace('_e', '')
            create_metabolite(model, actual_id, name, formula, compartment)
        else:
            create_metabolite(model, cpd_id, name, formula, compartment)

    # Add core reactions

    # ATP maintenance
    atpm = Reaction('ATPM')
    atpm.name = 'ATP maintenance'
    atpm.lower_bound = 0
    atpm.upper_bound = 1000
    atpm.add_metabolites({
        model.metabolites.get_by_id('cpd00002_c0'): -1,
        model.metabolites.get_by_id('cpd00001_c0'): -1,
        model.metabolites.get_by_id('cpd00008_c0'): 1,
        model.metabolites.get_by_id('cpd00009_c0'): 1,
        model.metabolites.get_by_id('cpd00067_c0'): 1,
    })

    # Phosphotransacetylase (rxn00173)
    pta = Reaction('rxn00173_c0')
    pta.name = 'Phosphotransacetylase'
    pta.lower_bound = -1000
    pta.upper_bound = 1000
    pta.add_metabolites({
        model.metabolites.get_by_id('cpd00009_c0'): -1,
        model.metabolites.get_by_id('cpd00022_c0'): -1,
        model.metabolites.get_by_id('cpd00010_c0'): 1,
        model.metabolites.get_by_id('cpd00196_c0'): 1,
    })

    # Acetate kinase (rxn00225)
    ack = Reaction('rxn00225_c0')
    ack.name = 'Acetate kinase'
    ack.lower_bound = -1000
    ack.upper_bound = 1000
    ack.add_metabolites({
        model.metabolites.get_by_id('cpd00002_c0'): -1,
        model.metabolites.get_by_id('cpd00029_c0'): -1,
        model.metabolites.get_by_id('cpd00008_c0'): 1,
        model.metabolites.get_by_id('cpd00196_c0'): 1,
    })

    # Thioredoxin reductase (NADPH) - rxn27318
    trxr = Reaction('rxn27318_c0')
    trxr.name = 'Thioredoxin reductase (NADPH)'
    trxr.lower_bound = -1000
    trxr.upper_bound = 1000
    trxr.add_metabolites({
        model.metabolites.get_by_id('cpd00005_c0'): -1,
        model.metabolites.get_by_id('cpd00067_c0'): -1,
        model.metabolites.get_by_id('cpd27735_c0'): -1,
        model.metabolites.get_by_id('cpd00006_c0'): 1,
        model.metabolites.get_by_id('cpd28060_c0'): 1,
    })

    # NADPH regeneration (simplified transhydrogenase)
    nadph_regen = Reaction('NADPH_regen')
    nadph_regen.name = 'NADPH regeneration'
    nadph_regen.lower_bound = 0
    nadph_regen.upper_bound = 1000
    nadph_regen.add_metabolites({
        model.metabolites.get_by_id('cpd00004_c0'): -1,
        model.metabolites.get_by_id('cpd00006_c0'): -1,
        model.metabolites.get_by_id('cpd00003_c0'): 1,
        model.metabolites.get_by_id('cpd00005_c0'): 1,
    })

    # NADH regeneration (ferredoxin-based, simplified)
    nadh_regen = Reaction('NADH_regen')
    nadh_regen.name = 'NADH regeneration (fermentative)'
    nadh_regen.lower_bound = 0
    nadh_regen.upper_bound = 1000
    nadh_regen.add_metabolites({
        model.metabolites.get_by_id('cpd00003_c0'): -1,
        model.metabolites.get_by_id('cpd00067_c0'): -2,
        model.metabolites.get_by_id('cpd00004_c0'): 1,
    })

    # Acetate export
    ace_export = Reaction('ACE_export')
    ace_export.name = 'Acetate export'
    ace_export.lower_bound = 0
    ace_export.upper_bound = 1000
    ace_export.add_metabolites({
        model.metabolites.get_by_id('cpd00029_c0'): -1,
        model.metabolites.get_by_id('cpd00029_e0'): 1,
    })

    model.add_reactions([atpm, pta, ack, trxr, nadph_regen, nadh_regen, ace_export])


def add_betaine_reductase_pathway(model):
    """
    Add betaine reductase pathway for Halanaerobium.

    rxn17220: Betaine + Red-Thioredoxin + Pi + 2H+ → TMA + Acetylphosphate + Ox-Thioredoxin + H2O
    """
    # Betaine transport (e0 -> c0)
    bet_import = Reaction('BET_import')
    bet_import.name = 'Betaine import'
    bet_import.lower_bound = 0
    bet_import.upper_bound = 1000

    # Create cytoplasmic betaine
    bet_c = create_metabolite(model, 'cpd00540', 'Glycine betaine', 'C5H11NO2', 'c0')
    bet_import.add_metabolites({
        model.metabolites.get_by_id('cpd00540_e0'): -1,
        bet_c: 1,
    })

    # Betaine reductase (rxn17220)
    bet_red = Reaction('rxn17220_c0')
    bet_red.name = 'Betaine reductase (EC 1.21.4.4)'
    bet_red.lower_bound = 0
    bet_red.upper_bound = 1000

    # Create cytoplasmic TMA
    tma_c = create_metabolite(model, 'cpd00441', 'Trimethylamine', 'C3H9N', 'c0')

    bet_red.add_metabolites({
        bet_c: -1,
        model.metabolites.get_by_id('cpd28060_c0'): -1,  # Red-thioredoxin
        model.metabolites.get_by_id('cpd00009_c0'): -1,  # Phosphate
        model.metabolites.get_by_id('cpd00067_c0'): -2,  # 2 H+
        tma_c: 1,
        model.metabolites.get_by_id('cpd00196_c0'): 1,   # Acetylphosphate
        model.metabolites.get_by_id('cpd27735_c0'): 1,   # Ox-thioredoxin
        model.metabolites.get_by_id('cpd00001_c0'): 1,   # H2O
    })

    # TMA export
    tma_export = Reaction('TMA_export')
    tma_export.name = 'TMA export'
    tma_export.lower_bound = 0
    tma_export.upper_bound = 1000
    tma_export.add_metabolites({
        tma_c: -1,
        model.metabolites.get_by_id('cpd00441_e0'): 1,
    })

    model.add_reactions([bet_import, bet_red, tma_export])


def add_stickland_pathway(model, include_amino_acid_input=True):
    """
    Add Stickland fermentation pathway with betaine as electron acceptor.

    Stickland: Amino acid (donor) + Betaine (acceptor) →
               Organic acid + CO2 + NH3 + TMA + Acetate
    """
    # First add betaine reductase pathway
    add_betaine_reductase_pathway(model)

    if include_amino_acid_input:
        # Simplified alanine oxidation as Stickland donor half-reaction
        # Alanine → Pyruvate + NH3 + 2[H]
        # Pyruvate → Acetyl-CoA + CO2

        # Create pyruvate
        pyr = create_metabolite(model, 'cpd00020', 'Pyruvate', 'C3H3O3', 'c0')

        # Alanine dehydrogenase
        ala_dh = Reaction('ALA_DH')
        ala_dh.name = 'Alanine dehydrogenase (Stickland donor)'
        ala_dh.lower_bound = 0
        ala_dh.upper_bound = 1000
        ala_dh.add_metabolites({
            model.metabolites.get_by_id('cpd00035_c0'): -1,  # L-Alanine
            model.metabolites.get_by_id('cpd00003_c0'): -1,  # NAD+
            model.metabolites.get_by_id('cpd00001_c0'): -1,  # H2O
            pyr: 1,
            model.metabolites.get_by_id('cpd00013_c0'): 1,   # NH3
            model.metabolites.get_by_id('cpd00004_c0'): 1,   # NADH
            model.metabolites.get_by_id('cpd00067_c0'): 1,   # H+
        })

        # Pyruvate:ferredoxin oxidoreductase (simplified)
        pfor = Reaction('PFOR')
        pfor.name = 'Pyruvate:ferredoxin oxidoreductase'
        pfor.lower_bound = 0
        pfor.upper_bound = 1000
        pfor.add_metabolites({
            pyr: -1,
            model.metabolites.get_by_id('cpd00010_c0'): -1,  # CoA
            model.metabolites.get_by_id('cpd00022_c0'): 1,   # Acetyl-CoA
            model.metabolites.get_by_id('cpd00011_c0'): 1,   # CO2
            model.metabolites.get_by_id('cpd00004_c0'): 1,   # NADH (simplified from Fd)
        })

        # Alanine import
        ala_import = Reaction('ALA_import')
        ala_import.name = 'Alanine import'
        ala_import.lower_bound = 0
        ala_import.upper_bound = 1000
        ala_e = create_metabolite(model, 'cpd00035', 'L-Alanine', 'C3H7NO2', 'e0')
        ala_import.add_metabolites({
            ala_e: -1,
            model.metabolites.get_by_id('cpd00035_c0'): 1,
        })

        # NH3 export
        nh3_export = Reaction('NH3_export')
        nh3_export.name = 'Ammonia export'
        nh3_export.lower_bound = 0
        nh3_export.upper_bound = 1000
        nh3_e = create_metabolite(model, 'cpd00013', 'NH3', 'H3N', 'e0')
        nh3_export.add_metabolites({
            model.metabolites.get_by_id('cpd00013_c0'): -1,
            nh3_e: 1,
        })

        # CO2 export
        co2_export = Reaction('CO2_export')
        co2_export.name = 'CO2 export'
        co2_export.lower_bound = 0
        co2_export.upper_bound = 1000
        co2_e = create_metabolite(model, 'cpd00011', 'CO2', 'CO2', 'e0')
        co2_export.add_metabolites({
            model.metabolites.get_by_id('cpd00011_c0'): -1,
            co2_e: 1,
        })

        model.add_reactions([ala_dh, pfor, ala_import, nh3_export, co2_export])


def add_exchange_reactions(model):
    """Add exchange reactions for extracellular metabolites."""

    exchanges = [
        ('cpd00540_e0', 'Glycine betaine', -1000, 1000),  # Betaine uptake
        ('cpd00441_e0', 'Trimethylamine', 0, 1000),       # TMA export
        ('cpd00029_e0', 'Acetate', 0, 1000),              # Acetate export
    ]

    for met_id, name, lb, ub in exchanges:
        if met_id in model.metabolites:
            ex_rxn = Reaction(f'EX_{met_id}')
            ex_rxn.name = f'{name} exchange'
            ex_rxn.lower_bound = lb
            ex_rxn.upper_bound = ub
            ex_rxn.add_metabolites({model.metabolites.get_by_id(met_id): -1})
            model.add_reactions([ex_rxn])

    # Add Stickland-specific exchanges if metabolites exist
    stickland_exchanges = [
        ('cpd00035_e0', 'L-Alanine', -1000, 0),          # Alanine uptake
        ('cpd00013_e0', 'Ammonia', 0, 1000),             # NH3 export
        ('cpd00011_e0', 'CO2', 0, 1000),                 # CO2 export
    ]

    for met_id, name, lb, ub in stickland_exchanges:
        if met_id in model.metabolites:
            ex_rxn = Reaction(f'EX_{met_id}')
            ex_rxn.name = f'{name} exchange'
            ex_rxn.lower_bound = lb
            ex_rxn.upper_bound = ub
            ex_rxn.add_metabolites({model.metabolites.get_by_id(met_id): -1})
            model.add_reactions([ex_rxn])

    # Add H+ and H2O sinks to balance reactions
    if 'cpd00067_c0' in model.metabolites:
        h_sink = Reaction('H_sink')
        h_sink.name = 'Proton sink'
        h_sink.lower_bound = -1000
        h_sink.upper_bound = 1000
        h_sink.add_metabolites({model.metabolites.get_by_id('cpd00067_c0'): -1})
        model.add_reactions([h_sink])

    if 'cpd00001_c0' in model.metabolites:
        h2o_sink = Reaction('H2O_sink')
        h2o_sink.name = 'Water sink'
        h2o_sink.lower_bound = -1000
        h2o_sink.upper_bound = 1000
        h2o_sink.add_metabolites({model.metabolites.get_by_id('cpd00001_c0'): -1})
        model.add_reactions([h2o_sink])

    # Add phosphate source for betaine reductase
    if 'cpd00009_c0' in model.metabolites:
        pi_source = Reaction('Pi_source')
        pi_source.name = 'Phosphate source'
        pi_source.lower_bound = 0
        pi_source.upper_bound = 1000
        pi_source.add_metabolites({model.metabolites.get_by_id('cpd00009_c0'): 1})
        model.add_reactions([pi_source])

    # Add NADPH source for thioredoxin reductase
    if 'cpd00005_c0' in model.metabolites and 'cpd00006_c0' in model.metabolites:
        nadph_source = Reaction('NADPH_source')
        nadph_source.name = 'NADPH source (simplified)'
        nadph_source.lower_bound = 0
        nadph_source.upper_bound = 1000
        nadph_source.add_metabolites({
            model.metabolites.get_by_id('cpd00006_c0'): -1,
            model.metabolites.get_by_id('cpd00005_c0'): 1,
        })
        model.add_reactions([nadph_source])

    # Add NADP+ sink
    if 'cpd00006_c0' in model.metabolites:
        nadp_sink = Reaction('NADP_sink')
        nadp_sink.name = 'NADP+ sink'
        nadp_sink.lower_bound = 0
        nadp_sink.upper_bound = 1000
        nadp_sink.add_metabolites({model.metabolites.get_by_id('cpd00006_c0'): -1})
        model.add_reactions([nadp_sink])


def build_betaine_model(mag_id, mag_info):
    """Build a complete betaine-utilizing model for a MAG."""

    # Create base model
    model = Model(mag_id)
    model.name = f"Betaine model: {mag_info['taxonomy'].split(';')[-1].strip()}"

    # Add core metabolism
    add_core_metabolism(model)

    # Add pathway based on type
    if mag_info['pathway_type'] == 'betaine_reductase':
        add_betaine_reductase_pathway(model)
    else:  # stickland
        add_stickland_pathway(model)

    # Add exchange reactions
    add_exchange_reactions(model)

    # Set objective to TMA production (main product of betaine reduction)
    if 'EX_cpd00441_e0' in model.reactions:
        model.objective = 'EX_cpd00441_e0'

    return model


def test_betaine_utilization(model):
    """Test if the model can utilize betaine."""
    results = {
        'feasible': False,
        'betaine_uptake': 0,
        'tma_production': 0,
        'acetate_production': 0,
        'objective_value': 0
    }

    try:
        # Set betaine uptake
        if 'EX_cpd00540_e0' in model.reactions:
            model.reactions.get_by_id('EX_cpd00540_e0').lower_bound = -10

        # For Stickland, also provide alanine
        if 'EX_cpd00035_e0' in model.reactions:
            model.reactions.get_by_id('EX_cpd00035_e0').lower_bound = -10

        solution = model.optimize()

        if solution.status == 'optimal':
            results['feasible'] = True
            results['objective_value'] = solution.objective_value

            if 'EX_cpd00540_e0' in model.reactions:
                results['betaine_uptake'] = -solution.fluxes.get('EX_cpd00540_e0', 0)
            if 'EX_cpd00441_e0' in model.reactions:
                results['tma_production'] = solution.fluxes.get('EX_cpd00441_e0', 0)
            if 'EX_cpd00029_e0' in model.reactions:
                results['acetate_production'] = solution.fluxes.get('EX_cpd00029_e0', 0)

    except Exception as e:
        results['error'] = str(e)

    return results


def main():
    """Build all betaine MAG models."""

    base_dir = Path(__file__).parent
    output_dir = base_dir / 'betaine_models' / 'individual'
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Building Individual Betaine-Utilizing MAG Models")
    print("=" * 70)
    print(f"\nTimestamp: {datetime.now().isoformat()}")
    print(f"Output directory: {output_dir}")
    print(f"\nProcessing {len(BETAINE_MAGS)} betaine-utilizing MAGs")

    results = []
    pathway_counts = {'betaine_reductase': 0, 'stickland': 0}

    for mag_id, mag_info in BETAINE_MAGS.items():
        print(f"\n{'='*60}")
        print(f"Processing: {mag_id.split('.')[-2]}")
        print(f"  Taxonomy: {mag_info['taxonomy']}")
        print(f"  Pathway: {mag_info['pathway_type']}")
        print(f"  Sample: {mag_info['sample']}")

        # Build model
        model = build_betaine_model(mag_id, mag_info)
        print(f"  Created model with {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")

        # Test betaine utilization
        test_results = test_betaine_utilization(model)
        print(f"  Feasible: {test_results['feasible']}")
        print(f"  TMA production: {test_results['tma_production']:.2f}")
        print(f"  Acetate production: {test_results['acetate_production']:.2f}")

        # Save model
        short_name = mag_id.split('.')[-2] if '.' in mag_id else mag_id.split('_')[-1]
        output_path = output_dir / f"{short_name}_betaine_model.json"
        cobra.io.save_json_model(model, str(output_path))
        print(f"  Saved to: {output_path.name}")

        pathway_counts[mag_info['pathway_type']] += 1

        results.append({
            'mag_id': mag_id,
            'short_name': short_name,
            'taxonomy': mag_info['taxonomy'],
            'pathway_type': mag_info['pathway_type'],
            'sample': mag_info['sample'],
            'num_reactions': len(model.reactions),
            'num_metabolites': len(model.metabolites),
            **test_results
        })

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\nTotal models built: {len(results)}")
    print(f"  - Betaine reductase pathway: {pathway_counts['betaine_reductase']}")
    print(f"  - Stickland fermentation: {pathway_counts['stickland']}")

    feasible = sum(1 for r in results if r['feasible'])
    print(f"\nModels with feasible betaine utilization: {feasible}/{len(results)}")

    # Save results
    results_path = output_dir / 'model_summary.json'
    with open(results_path, 'w') as f:
        json.dump({
            'timestamp': datetime.now().isoformat(),
            'total_mags': len(BETAINE_MAGS),
            'pathway_counts': pathway_counts,
            'models': results
        }, f, indent=2)
    print(f"\nResults saved to: {results_path}")

    # Print taxonomy breakdown
    print("\n" + "-" * 60)
    print("Taxonomy Breakdown:")
    class_counts = {}
    for r in results:
        cls = BETAINE_MAGS[r['mag_id']]['class']
        class_counts[cls] = class_counts.get(cls, 0) + 1
    for cls, count in sorted(class_counts.items(), key=lambda x: -x[1]):
        feasible_in_class = sum(1 for r in results
                                if BETAINE_MAGS[r['mag_id']]['class'] == cls and r['feasible'])
        print(f"  {cls}: {count} MAGs ({feasible_in_class} feasible)")

    return results


if __name__ == '__main__':
    main()
