import marimo

__generated_with = "0.23.4"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Loading content
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## normalizing abundances
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### volume-weighting DNA measurements
    """)
    return


@app.cell
def _():
    from math import pi
    from pandas import read_csv, Series, DataFrame
    import os
    import cobra
    # from modelseedpy import MSBuilder, MSATPCorrection, MSMedia, MSGapfill
    # import modelseedpy
    import json
    from tqdm import tqdm
    # from glob import glob
    # import cobrakbase
    dna_df = read_csv('../data/SaltPondsDNA.csv').set_index('Sample')
    _sample_dna = dna_df['DNA conc (µg/g soil)'].to_dict()
    _unique_prefixes = {'_'.join(i.split('_')[:-1]) for i in dna_df.index}
    _volume_weighted_sample_dna = Series({prefix: (5 * _sample_dna[prefix + '_D1'] + 10 * _sample_dna[prefix + '_D2']) / 15 for prefix in _unique_prefixes})

    def fluxes_to_emissions(row):
        area_cm2 = pi * (5 / 2) ** 2
        area_m2 = area_cm2 * 0.01 ** 2
        volume_cm3 = 294.52
        ugDNA_gSoil = _volume_weighted_sample_dna.reindex(_row.index).fillna(2)
        DNA_per_biomass = 0.1
        gSoil_cm3 = 0.7  # TODO: LaTeX of the unit conversion
        gBiomass_gSoil = gSoil_cm3 * volume_cm3 * ugDNA_gSoil / DNA_per_biomass / 1000 / 1000  # mmol/hr/g_bio -> umol/(m^2 * d)  via  / area (m^2) * (1000 umol / mmol) * (24 hours / day) * (g biomass / container)
        return _row / area_m2 * 1000.0 * 24 * gBiomass_gSoil
      # defining the dimensions of the sampling system
    def merge_depths_by_DNA(abundances):  # cm^2
        return DataFrame({prefix: (abundances[prefix + '_D1'] * _sample_dna[prefix + '_D1'] + abundances[prefix + '_D2'] * _sample_dna[prefix + '_D2']) / (_sample_dna[prefix + '_D1'] + _sample_dna[prefix + '_D2']) for prefix in _unique_prefixes}, index=abundances.index)  # m^2  # cm^3  # DNA per soil  # 2 is the approximate average from the provided spreadsheet  # estimated from Chris' spreadsheet  # return DataFrame({prefix: mean(abundances[prefix+"_D1"] + abundances[prefix+"_D2"]) for prefix in unique_prefixes},  #                     index=abundances.index)

    return (
        DataFrame,
        Series,
        cobra,
        dna_df,
        fluxes_to_emissions,
        json,
        merge_depths_by_DNA,
        os,
        pi,
        read_csv,
        tqdm,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### normalize the
    """)
    return


@app.cell
def _(Series, display, dna_df, json, merge_depths_by_DNA, read_csv):
    # load data
    total_df = read_csv('../data/Cliff_Sample_Metadata_BGC_NMR.csv').set_index('Sample')
    abundances = read_csv('../data/Cliff_310MAG_relabund.tsv', sep='\t')
    normalized_abundances = abundances.copy()
    # create DNA-normalized abundances for each depth
    #TODO is this necessary?
    for _col, _val in dna_df['DNA conc (µg/g soil)'].items():
        normalized_abundances[_col] = normalized_abundances[_col] / _val
    _unrestored_cols = [_col for _col in normalized_abundances.columns if 'R1_' in _col or 'R2_' in _col]
    unrestored_avg_cols = list(set([_col.split('_D')[0] for _col in _unrestored_cols]))
    # find the 
    normalized_abundances.index.name = 'MAG'  # and "restored" not in col]
    normalized_abundances.to_csv('normalized_MAG_abundances.csv')
    display(normalized_abundances)
    normalized_methane_df = total_df['CH4_umol_m2_d'] / dna_df['DNA conc (µg/g soil)']
    display(normalized_methane_df)
    _unique_prefixes = set(['_'.join(_col.split('_')[:-1]) for _col in normalized_methane_df.index])
    col_merged = {prefix: [prefix + '_D1', prefix + '_D2'] for prefix in _unique_prefixes}
    averaged_normalized_abundances = merge_depths_by_DNA(normalized_abundances)
    unrestored_cols_simplified = [_col for _col in averaged_normalized_abundances.columns if 'R1_' in _col or 'R2_' in _col]
    print(_unrestored_cols, unrestored_avg_cols)
    # averaged_normalized_abundances = DataFrame({prefix: normalized_abundances[cols].mean(axis=1) for prefix, cols in col_merged.items()}, index=normalized_abundances.index)
    averaged_normalized_abundances.to_csv('../data/averaged_normalized_MAG_abundances.csv')
    display(averaged_normalized_abundances)  # and "restored" not in col]
    averaged_unrestored_normalized_abundances = averaged_normalized_abundances[unrestored_avg_cols]
    # averaged_normalized_abundances = normalized_abundances.groupby(prefixes, axis=1).mean()
    averaged_unrestored_normalized_abundances.to_csv('../data/averaged_normalized_abundances_unrestored.csv')
    averaged_normalized_methane = Series({prefix: normalized_methane_df[cols].mean() for prefix, cols in col_merged.items()}, index=col_merged.keys())
    averaged_unrestored_normalized_methane = averaged_normalized_methane[unrestored_avg_cols]
    methane_dic = averaged_unrestored_normalized_methane.to_dict()
    json.dump(methane_dic, open('averaged_normalized_methane_unrestored.json', 'w'), indent=3)
    methane_dic
    return (
        abundances,
        averaged_unrestored_normalized_abundances,
        averaged_unrestored_normalized_methane,
        methane_dic,
        normalized_abundances,
        total_df,
    )


@app.cell
def _():
    media = {'EX_cpd00067_e0': 1000.0,
     'EX_cpd00058_e0': 1000.0,
     'EX_cpd00013_e0': 1000.0,
     'EX_cpd00244_e0': 1000.0,
     'EX_cpd00205_e0': 1000.0,
     'EX_cpd00034_e0': 1000.0,
     'EX_cpd11574_e0': 1000.0,
     'EX_cpd00971_e0': 1000.0,
     'EX_cpd00048_e0': 1000.0,
     'EX_cpd00030_e0': 1000.0,
     'EX_cpd00305_e0': 100.0,
     'EX_cpd00001_e0': 1000.0,
     'EX_cpd10516_e0': 1000.0,
     'EX_cpd00007_e0': 1000.0,
     'EX_cpd00159_e0': 100.0,
     'EX_cpd25960_e0': 1000.0,
     "EX_cpd00009_e0": 100,
     'EX_cpd00063_e0': 1000.0,
     'EX_cpd00149_e0': 1000.0,
     'EX_cpd00254_e0': 1000.0,
     'EX_cpd00099_e0': 1000.0}
    return (media,)


@app.cell
def _():
    mepe = {  #mepn
        'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.20.contigs',
     'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.32.contigs',
     'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.51.contigs',
     'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_metabat.27.contigs',
     'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_metabat.28.contigs',
     'Salt_Pond_MetaG_R1_B_H2O_MG_DASTool_bins_metabat.7.contigs',
     'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_metabat.31.contigs',
     'Salt_Pond_MetaG_R1_C_H2O_MG_DASTool_bins_concoct_out.79.contigs',
     'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.35.contigs',
     'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.49.contigs',
     'Salt_Pond_MetaG_R2_A_D2_MG_DASTool_bins_metabat.14.contigs',
     'Salt_Pond_MetaG_R2_A_D2_MG_DASTool_bins_metabat.16.contigs',
     'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_concoct_out.143.contigs',
     'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_concoct_out.32.contigs',
     'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.10.contigs',
     'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.16.contigs',
     'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.19.contigs',
     'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.38.contigs',
     'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.55.contigs',
     'Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_concoct_out.88.contigs',
     'Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_concoct_out.95.contigs',
     'Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_metabat.17.contigs',
     'Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.18.contigs',
     'Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.45.contigs',
     'Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.47.contigs',
     'Salt_Pond_MetaG_R2_restored_C_black_MG_DASTool_bins_concoct_out.17.contigs',
     'Salt_Pond_MetaG_R2_restored_C_black_MG_DASTool_bins_metabat.8.contigs'}
    mepe_striped = [x.replace(".contigs","") for x in mepe]


    betaine = {'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs',
     'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs',
     'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs',
     'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs',
     'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs',
     'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs',
     'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs',
     'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs',
     'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs',
     'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs',
     'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs',
     'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs',
     'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs',
     'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs',
     'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs'}
    betaine_striped = [x.replace(".contigs","") for x in betaine]
    return betaine, betaine_striped, mepe, mepe_striped


@app.cell
def _(cobra, mepe, os, tqdm):
    models = {}
    for i in tqdm(mepe):
        # filename_model = f'./genome_scale/{i}.json'
        filename_model = f'../models/MePn_consumers/{i}.json'
        if os.path.exists(filename_model):
            models[i] = cobra.io.load_json_model(filename_model)
        else:
            print(f'skip {i}')
    return i, models


@app.cell
def _(media):
    from modelseedpy.biochem import from_local
    # msdb = from_local("../../../shared/code/ModelSEEDDatabase")
    msdb = from_local("../../../ModelSEEDDatabase")

    for cpd in media:
        print(msdb.compounds.get_by_id(cpd.split("_")[1]).name, "\t", cpd.id)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # determine methane emission from methylphosphate as the sole phosphate source on lactate
    """)
    return


@app.cell
def _(display, media, models):
    methane_emissions = {}
    for _ID, _model in models.items():
        _model.medium = media
        _sol = _model.optimize()
        if _sol.objective_value == 0:
            print(f'{_ID} failed')  #, end="\n\n")
            continue
        methane_emissions[_ID] = _sol.fluxes['EX_cpd01024_e0']
    methane_emissions = dict(sorted(methane_emissions.items(), key=lambda item: item[1], reverse=True))  # print(f"{sol.fluxes['EX_cpd01024_e0']} methane flux from {ID}")
    print()
    print('methane emissions in optimum growth')
    display(methane_emissions)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## compare with the abundances of each respective species
    """)
    return


@app.cell
def _(abundances, display):
    display(abundances)
    unrestored = abundances[abundances.index.to_series().str.contains('R1|R2')]
    member_totals = unrestored.sum(axis=1).sort_values(ascending=False)
    most_abundant_members = member_totals[:8]
    display(most_abundant_members)
    _unrestored_cols = [_col for _col in abundances.columns if 'R1_' in _col or 'R2_' in _col]
    display(abundances[_unrestored_cols].sum(axis=1).sort_values(ascending=False)[:8])
    _sample_totals = abundances.sum(axis=0)
    display(_sample_totals)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # notes
    - the Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_metabat.17 organism both generates the most methane emission of all organisms and is the only one of the most abundant members that was defined with a model in the MePn set of organisms
    -
    """)
    return


@app.cell
def _(display, dna_df, total_df):
    def convert_methane(methane):
        from math import pi
        sample_area = pi * (5 / 2) ** 2
        sample_area = sample_area * 0.01 ** 2
        gBiomass_gSoil = 5.9 * 1000
        return methane * sample_area / (1000.0 * 24)
    methane = convert_methane(total_df['CH4_umol_m2_d'])
    dna_df.drop(columns=[_col for _col in dna_df.columns if 'named' in _col.lower()], inplace=True)
    for _sample in dna_df.index:
        dna_df.loc[_sample, 'methane'] = methane[_sample]
        dna_df.loc[_sample, 'CH4 (mmol/hr) / normalized_DNA'] = dna_df.loc[_sample, 'methane'] / dna_df.loc[_sample, 'DNA conc (µg/g soil)']
    display(dna_df)
    dna_df.to_csv('../data/SaltPondsDNA_methane_weighted.csv')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # correlate the normalized abundances of members with the normalized methane flux
    """)
    return


@app.cell
def _(display, mepe_striped, normalized_abundances, read_csv):
    normalized_methane = read_csv('../data/SaltPondsDNA_methane_weighted.csv').set_index('Sample')['CH4 (mmol/hr) / normalized_DNA']
    display(normalized_methane)
    correlations = {}
    for _index, _row in normalized_abundances.iterrows():
        _corr = _row.corr(normalized_methane, method='spearman')
        correlations[_index] = _corr
    correlations = dict(sorted(correlations.items(), key=lambda item: item[1], reverse=True))
    merged_corr = 0
    for _org, _val in correlations.items():
        if _org in mepe_striped:
            merged_corr = merged_corr + _val
    print('max all correlations', max(correlations.values()))
    unrestored_correlations = {_k: _v for _k, _v in correlations.items() if ('R1_' in _k or 'R2_' in _k) and (not 'restored' in _k)}
    print(unrestored_correlations['Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_metabat.17'])
    display(unrestored_correlations)
    import numpy as np
    print(np.mean(list(unrestored_correlations.values())), np.mean(list(unrestored_correlations.values())) ** 2, np.std(list(unrestored_correlations.values())))
    return normalized_methane, np


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## look at the total abundances of MePn consumers and betaine reducers
    """)
    return


@app.cell
def _(betaine_striped, mepe_striped, normalized_abundances):
    merged_normalized_abundances = normalized_abundances.copy()
    merged_normalized_abundances.loc["mepe"] = normalized_abundances.loc[mepe_striped].sum()
    merged_normalized_abundances.loc["betaine_reducers"] = normalized_abundances.loc[betaine_striped].sum()
    merged_normalized_abundances = merged_normalized_abundances.drop(mepe_striped)
    merged_normalized_abundances
    return (merged_normalized_abundances,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # MePn reducers
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # TODO correlate with all the experimental metabolites, including MePn
    """)
    return


@app.cell
def _(merged_normalized_abundances, normalized_methane):
    merged_correlations = {}
    for _index, _row in merged_normalized_abundances.iterrows():
        _corr = _row.corr(normalized_methane, method='spearman')
        merged_correlations[_index] = _corr
    merged_correlations = dict(sorted(merged_correlations.items(), key=lambda item: item[1], reverse=True))  # print(f"{index} correlates {corr} with methane")
    # print("max all correlations", max(merged_correlations.values()))
    # unrestored_merged_correlations = {k:v for k,v in merged_correlations.items() if ("R1_" in k or "R2_" in k) and not "restored" in k}
    # unrestored_merged_correlations
    merged_correlations
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Simulation
    """)
    return


@app.cell
def _(cobra, display, media, models):
    display(cobra.Configuration())
    dontGrow = []
    for _ID, _model in models.items():
        ogMedia = _model.medium
        _model.medium = media
        _sol = _model.optimize()
        if _sol.status == 'infeasible':
            print(_ID, 'is infeasible')
        elif _sol.objective_value == 0:
            print(_ID, "doesn't grow")
            dontGrow.append(_ID)
        else:
            print(_sol)
        _model.medium = ogMedia
    return (dontGrow,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## simulate the MePn-consuming organisms
    """)
    return


@app.cell
def _(
    DataFrame,
    averaged_unrestored_normalized_abundances,
    dontGrow,
    fluxes_to_emissions,
    json,
    media,
    methane_dic,
    models,
    np,
):
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import r2_score
    from optlang.symbolics import Zero
    sols = {}
    _fits, fluxes = ({}, {})
    with open('sols_log.txt', 'w') as _f:
        for multiple in [0.07, 0.33]:
            fluxes[multiple] = {'methane': {}, 'lactate': {}, 'phosphate': {}, 'MePn': {}}
            for _orgID, _row in averaged_unrestored_normalized_abundances.iterrows():
                _ID = _orgID + '.contigs'
                if _ID not in models or _ID in dontGrow:
                    continue
                _model = models[_ID].copy()
                constraint = _model.problem.Constraint(Zero, lb=0, ub=None, name='Phosphate_ratio')
                _model.add_cons_vars(constraint)
                _model.solver.update()
                constraint.set_linear_coefficients(coefficients={_model.reactions.EX_cpd00009_e0.forward_variable: multiple, _model.reactions.EX_cpd00009_e0.reverse_variable: -multiple, _model.reactions.EX_cpd25960_e0.forward_variable: -1, _model.reactions.EX_cpd25960_e0.reverse_variable: 1})
                _model.solver.update()
                _model.medium = media
                _bio = _model.reactions.get_by_id('bio1')
                _model.solver.configuration.timeout = 15
                sols[_ID] = {}
                for _k, _dic in fluxes[multiple].items():
                    fluxes[multiple][_k][_ID] = {}
                _f.write(f'\n\n{_orgID}')
                for _sample, _abundance in _row.items():
                    _bio.bounds = (_abundance, _abundance)
                    _EX = _model.reactions.get_by_id('EX_cpd00159_e0')
                    _EX.bounds = (-1000, 0)
                    _model.objective = _EX
                    _f.write(f'\noptimizing {_sample} at {_abundance}')
                    sols[_ID][_sample] = _model.optimize()
                    if sols[_ID][_sample].status == 'infeasible':
                        _f.write(f'Feasibility Error: {_ID} {_sample}')
                    fluxes[multiple]['methane'][_ID][_sample] = sols[_ID][_sample].fluxes.EX_cpd01024_e0
                    fluxes[multiple]['lactate'][_ID][_sample] = sols[_ID][_sample].fluxes.EX_cpd00159_e0
                    fluxes[multiple]['phosphate'][_ID][_sample] = sols[_ID][_sample].fluxes.EX_cpd00009_e0
                    fluxes[multiple]['MePn'][_ID][_sample] = sols[_ID][_sample].fluxes.EX_cpd25960_e0
                    _bio.bounds = (0, 0)
                    _bio_flux = sols[_ID][_sample].fluxes.bio1
                    if _bio_flux != _abundance:
                        _f.write(f'\nAbundance Error: {_bio_flux} != {_abundance}')
                    _f.write(f'\noptimized at {sols[_ID][_sample].fluxes.EX_cpd00159_e0}')
            _mdl = LinearRegression()
            _xs = np.array(list(methane_dic.values())).reshape(-1, 1)
            _methane_df = DataFrame(fluxes[multiple]['methane'])
            _ys = np.array(fluxes_to_emissions(_methane_df.T.sum()).to_list())
            _mdl.fit(_xs, _ys)
            _y_pred = _mdl.predict(_xs)
            _fits[multiple] = {'eq': f'y = {_mdl.coef_[0]:.2f}x + {_mdl.intercept_:.2f}', 'R2': f'$R^2$ = {r2_score(_ys, _y_pred):.8f}'}
            mepn_df = DataFrame(fluxes[multiple]['MePn'])
            phosphate_df = DataFrame(fluxes[multiple]['phosphate'])
            print('\nMinimum MePn/Phosphate ratio', multiple)
            print('total MePn consumed:', mepn_df.sum(axis=1).to_dict())
            print('total inorganic Phosphate consumed:', phosphate_df.sum(axis=1).to_dict())
            print(_fits[multiple])
    json.dump(_fits, open('fits.json', 'w'), indent=3)
    json.dump(fluxes, open('saltern_fluxes.json', 'w'), indent=3)
    return LinearRegression, r2_score


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### compute the methane emission from substrate addition of MePn
    """)
    return


@app.cell
def _():
    additional_methane_generated = 10 # umol
    hours_equilibrated = 14 * 24 # hours
    volume_added = 10 # mL
    concentration_added = 100 # umol/L
    substrate_added = concentration_added * (volume_added/1000)  # umol
    print("substrate added", substrate_added)
    substrate_flux = substrate_added #  / hours_equilibrated  # not divided by time because all of the substrate was added at once
    additional_methane_flux = additional_methane_generated / hours_equilibrated  # umol/hr
    print(additional_methane_flux, "umol/hr methane flux from addition")
    return (substrate_flux,)


@app.cell
def _(
    DataFrame,
    LinearRegression,
    averaged_unrestored_normalized_abundances,
    display,
    dontGrow,
    media,
    methane_dic,
    models,
    np,
    r2_score,
    substrate_flux,
):
    media['EX_cpd25960_e0'] = substrate_flux / 1000
    display(media)
    print('members', len(averaged_unrestored_normalized_abundances.index))
    fits_addition, fluxes_addition = ({}, {})
    MePn_consumers = {}
    fluxes_addition = {'methane': {}, 'lactate': {}, 'phosphate': {}, 'MePn': {}}
    total_methane_emissions = np.array([0.0] * len(averaged_unrestored_normalized_abundances.columns))
    total_MePn_consumption = np.array([0.0] * len(averaged_unrestored_normalized_abundances.columns))
    for _orgID, _row in averaged_unrestored_normalized_abundances.iterrows():
        _ID = _orgID + '.contigs'
        if _ID not in models or _ID in dontGrow:
            continue
        _model = models[_ID].copy()
        _model.medium = media
        _bio = _model.reactions.get_by_id('bio1')
        _model.solver.configuration.timeout = 15
        for _k, _dic in fluxes_addition.items():
            fluxes_addition[_k][_ID] = {}
        for _sample, _abundance in _row.items():
            _bio.bounds = (_abundance, _abundance)
            _EX = _model.reactions.get_by_id('EX_cpd00159_e0')
            _model.objective = _EX
            _sol = _model.optimize()
            if _sol.status == 'infeasible':
                print(f'Feasibility Error: {_ID} in {_sample}')
                continue
            if abs(_sol.fluxes.EX_cpd25960_e0) > 1e-09:
                MePn_consumers[_sample] = MePn_consumers.get(_sample, set())
                MePn_consumers[_sample].add(_ID)
            if abs(_model.reactions.get_by_id('EX_cpd25960_e0').lower_bound) < abs(_sol.fluxes.EX_cpd25960_e0):
                display(_model.medium)
                raise ValueError(f"MePn flux {_sol.fluxes.EX_cpd25960_e0} is greater than the lower bound                 {_model.reactions.get_by_id('EX_cpd25960_e0').lower_bound} for {_ID} {_sample}")
            fluxes_addition['methane'][_ID][_sample] = _sol.fluxes.EX_cpd01024_e0
            fluxes_addition['lactate'][_ID][_sample] = _sol.fluxes.EX_cpd00159_e0
            fluxes_addition['phosphate'][_ID][_sample] = _sol.fluxes.EX_cpd00009_e0
            fluxes_addition['MePn'][_ID][_sample] = _sol.fluxes.EX_cpd25960_e0
            _bio.bounds = (0, 0)
            _bio_flux = _sol.fluxes.bio1
            if _bio_flux != _abundance:
                print(f'Abundance Error: {_bio_flux} != {_abundance}')
        _xs = np.array(list(methane_dic.values())).reshape(-1, 1)
        _methane_df = DataFrame(fluxes_addition['methane'])
        methane_producing_members = _methane_df.sum()[_methane_df.sum() > 0].index
        if len(methane_producing_members) == 0:
            print(f'Skipping {_ID} because no methane was produced')
            continue
        _ys = np.array(_methane_df.T.sum().to_list())
        _mdl = LinearRegression()
        _mdl.fit(_xs, _ys)
        _y_pred = _mdl.predict(_xs)
        fits_addition = {'eq': f'y = {_mdl.coef_[0]:.2f}x + {_mdl.intercept_:.2f}', 'R2': f'$R^2$ = {r2_score(_ys, _y_pred):.8f}'}
        mepn_addition_df = DataFrame(fluxes_addition['MePn'])
        phosphate_addition_df = DataFrame(fluxes_addition['phosphate'])
        print(f'\nSubstrate addition: {_ID}')
        print(_ys)
        total_methane_emissions = total_methane_emissions + _ys * 1000
        total_MePn_consumption = total_MePn_consumption + np.array(mepn_addition_df.T.sum().to_list()) * 1000
        print('total MePn consumed:', mepn_addition_df.sum(axis=1).to_dict())
        print('total inorganic Phosphate consumed:', phosphate_addition_df.sum(axis=1).to_dict())
        print(fits_addition)
    print('methane producing members', methane_producing_members)
    print('MePn consumers', MePn_consumers)
    print('Total methane emissions from addition (umol/hr):', dict(zip(averaged_unrestored_normalized_abundances.columns, total_methane_emissions)))
    print('Total MePn consumption from addition (umol/hr):', dict(zip(averaged_unrestored_normalized_abundances.columns, total_MePn_consumption)))
    return (MePn_consumers,)


@app.cell
def _(MePn_consumers, display, read_csv):
    display(MePn_consumers)
    mepn_consumers_all_samples = set([_v for x in MePn_consumers.values() for _v in x])
    phylogeny = read_csv('../data/Saltern_phylogeny.csv')
    phylogeny = dict(zip(phylogeny['User Genome'], phylogeny['Classification']))
    mepn_consumers_phylogeny = {mem: phylogeny[mem + '__'] for mem in mepn_consumers_all_samples}
    mepn_consumers_phylogeny
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.55.contigs = Roseovarius
    Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.45.contigs = Roseovarius
    Salt_Pond_MetaG_R2_restored_C_black_MG_DASTool_bins_metabat.8.contigs = Albimonas
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## simulate the betaine-consuming organisms
    """)
    return


@app.cell
def _(betaine, cobra, i, os, tqdm):
    from cobra.core import Reaction

    def _add_atpm(model):
        if 'ATPM_c0' not in _model.reactions:
            atpm = Reaction(f'ATPM_c0', f'ATPM', 'ATPM', 0, 1000)
            atpm.add_metabolites({_model.metabolites.cpd00001_c0: -1, _model.metabolites.cpd00002_c0: -1, _model.metabolites.cpd00008_c0: 1, _model.metabolites.cpd00009_c0: 1, _model.metabolites.cpd00067_c0: 1})
            _model.add_reactions([atpm])
        return _model
    models_1 = {}
    for mdl_name in tqdm(betaine):
        mdl_path = f'models/{mdl_name}__.RAST_gapfilled_betaine.json'
        if os.path.exists(mdl_path):
            print(mdl_name)
            models_1[mdl_name] = _add_atpm(cobra.io.load_json_model(mdl_path))
        else:
            print(f'skip {i}')
    return (models_1,)


@app.cell
def _(models_1):
    list(models_1.values())[0].optimize()
    return


@app.cell
def _(models_1):
    from cobra.io import load_json_model
    from glob import glob
    betaine_reactions = {'rxn17220_c0': 'Betaine reductase (BET + Red-Thioredoxin + Pi → TMA + AcP + Ox-Thioredoxin)', 'rxn07207_c0': 'Glycine reductase (Glycine + Red-Thioredoxin + Pi → AcP + NH3 + Ox-Thioredoxin)', 'rxn27318_c0': 'Thioredoxin reductase (NADPH + Ox-Thioredoxin → NADP + Red-Thioredoxin)'}
    betaine_metabolites = {'cpd00540_c0': 'Glycine betaine (cytoplasmic)', 'cpd00540_e0': 'Glycine betaine (extracellular)', 'cpd00441_c0': 'Trimethylamine (TMA)'}
    betaine_transport = {'EX_cpd00540_e0': 'Betaine exchange'}
    gapfilled_models = sorted(glob('models/*_gapfilled_betaine.json'))
    print(f'{len(gapfilled_models)} gapfilled betaine models\n')
    for short_name, _model in models_1.items():
        print(f'=== {short_name} ({len(_model.reactions)} rxns) ===')
        for rxn_id, desc in {**betaine_reactions, **betaine_transport}.items():
            _present = rxn_id in _model.reactions
            if _present:
                _rxn = _model.reactions.get_by_id(rxn_id)
                print(f'  [YES] {rxn_id}: {desc}  bounds=({_rxn.lower_bound}, {_rxn.upper_bound})')
            else:
                print(f'  [ NO] {rxn_id}: {desc}')
        for met_id, desc in betaine_metabolites.items():
            _present = met_id in _model.metabolites
            print(f"  {('[YES]' if _present else '[ NO]')} {met_id}: {desc}")
        if 'rxn17220_c0' in _model.reactions:
            with _model:
                _model.objective = 'rxn17220_c0'
                _model.medium = {**_model.medium, 'EX_cpd00540_e0': 100}
                _sol = _model.optimize()
                print(f'  >> Betaine reductase max flux: {_sol.objective_value:.4f} ({_sol.status})')
        print()
    return


@app.cell
def _(
    DataFrame,
    LinearRegression,
    averaged_unrestored_normalized_abundances,
    fluxes_to_emissions,
    json,
    media,
    methane_dic,
    models_1,
    np,
    r2_score,
):
    sols_1 = {}
    with open('sols_log_betaine.txt', 'w') as _f:
        fluxes_1 = {'betaine': {}, 'TMA': {}, 'methane': {}, 'lactate': {}, 'phosphate': {}}
        for _orgID, _row in averaged_unrestored_normalized_abundances.iterrows():
            _ID = _orgID + '.contigs'
            if _ID not in models_1:
                continue
            _model = models_1[_ID].copy()
            print(_orgID, _model.slim_optimize(), end=' ')
            media.pop('EX_cpd25960_e0', None)
            _model.medium = {_k: _v for _k, _v in {**media, 'EX_cpd00540_e0': 0.1}.items() if _k in _model.reactions}
            print(_model.medium)
            _bio = _model.reactions.get_by_id('bio1')
            _model.solver.configuration.timeout = 15
            sols_1[_ID] = {}
            for _k, _dic in fluxes_1.items():
                fluxes_1[_k][_ID] = {}
            _f.write(f'\n\n{_orgID}')
            for _sample, _abundance in _row.items():
                _bio.bounds = (_abundance, _abundance)
                _EX = _model.reactions.get_by_id('EX_cpd00159_e0')
                _EX.bounds = (-1000, 0)
                _model.objective = _EX
                _f.write(f'\noptimizing {_sample} at {_abundance}')
                sols_1[_ID][_sample] = _model.optimize()
                if sols_1[_ID][_sample].status == 'infeasible':
                    _f.write(f'Feasibility Error: {_ID} {_sample}')
                fluxes_1['TMA'][_ID][_sample] = sols_1[_ID][_sample].fluxes.EX_cpd00441_e0
                fluxes_1['betaine'][_ID][_sample] = sols_1[_ID][_sample].fluxes.EX_cpd00540_e0
                fluxes_1['lactate'][_ID][_sample] = sols_1[_ID][_sample].fluxes.EX_cpd00159_e0
                fluxes_1['phosphate'][_ID][_sample] = sols_1[_ID][_sample].fluxes.EX_cpd00009_e0
                _bio.bounds = (0, 0)
                _bio_flux = sols_1[_ID][_sample].fluxes.bio1
                if _bio_flux != _abundance:
                    _f.write(f'\nAbundance Error: {_bio_flux} != {_abundance}')
                _f.write(f'\noptimized at {sols_1[_ID][_sample].fluxes.EX_cpd00159_e0}')
        _mdl = LinearRegression()
        _xs = np.array(list(methane_dic.values())).reshape(-1, 1)
        print(fluxes_1)
        _betaine_df = DataFrame(fluxes_1['betaine'])
        TMA_df = DataFrame(fluxes_1['TMA'])
        _ys = np.array(fluxes_to_emissions(_betaine_df.T.sum()).to_list())
        _mdl.fit(_xs, _ys)
        _y_pred = _mdl.predict(_xs)
        _fits = {'eq': f'y = {_mdl.coef_[0]:.2f}x + {_mdl.intercept_:.2f}', 'R2': f'$R^2$ = {r2_score(_ys, _y_pred):.8f}'}
        _methane_df = DataFrame(fluxes_1['methane'])
        print('total TMA flux:', TMA_df.sum(axis=1).to_dict())
        print('total betaine flux:', _betaine_df.sum(axis=1).to_dict())
        print('total methane flux:', _methane_df.sum(axis=1).to_dict())
        print(_fits)
    json.dump(_fits, open('fits_betaine.json', 'w'), indent=3)
    json.dump(fluxes_1, open('saltern_fluxes_betaine.json', 'w'), indent=3)
    return (sols_1,)


@app.cell
def _(DataFrame, np):
    fluxes_2 = {'betaine': {'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs': {'R1_C': np.float64(-0.1), 'R2_B': np.float64(-0.1), 'R2_C': np.float64(-0.1), 'R1_A': np.float64(-0.1), 'R2_A': np.float64(-0.1), 'R1_B': np.float64(-0.1)}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs': {'R1_C': np.float64(-0.004155173944330431), 'R2_B': np.float64(-0.016189024367368715), 'R2_C': np.float64(-0.059260881157459894), 'R1_A': np.float64(-0.003266969474764042), 'R2_A': np.float64(-0.0014186011073906296), 'R1_B': np.float64(-0.004386622400931961)}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(0.0), 'R2_C': np.float64(0.0), 'R1_A': np.float64(0.0), 'R2_A': np.float64(0.0), 'R1_B': np.float64(0.0)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs': {'R1_C': np.float64(-0.1), 'R2_B': np.float64(-0.1), 'R2_C': np.float64(-0.1), 'R1_A': np.float64(-0.1), 'R2_A': np.float64(-0.1), 'R1_B': np.float64(-0.1)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs': {'R1_C': np.float64(-0.04049835052607321), 'R2_B': np.float64(-0.1), 'R2_C': np.float64(-0.1), 'R1_A': np.float64(-0.04309736398630432), 'R2_A': np.float64(-0.02667366807003116), 'R1_B': np.float64(-0.03873150159343811)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs': {'R1_C': np.float64(-0.1), 'R2_B': np.float64(-0.1), 'R2_C': np.float64(-0.1), 'R1_A': np.float64(-0.1), 'R2_A': np.float64(-0.1), 'R1_B': np.float64(-0.1)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(0.0), 'R2_C': np.float64(0.0), 'R1_A': np.float64(0.0), 'R2_A': np.float64(0.0), 'R1_B': np.float64(0.0)}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(0.0), 'R2_C': np.float64(0.0), 'R1_A': np.float64(0.0), 'R2_A': np.float64(0.0), 'R1_B': np.float64(0.0)}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs': {'R1_C': np.float64(-0.010620624142021317), 'R2_B': np.float64(-0.017376308062454465), 'R2_C': np.float64(-0.05538545280165857), 'R1_A': np.float64(-0.02773869657846126), 'R2_A': np.float64(-0.004749024914201657), 'R1_B': np.float64(-0.04266112324083004)}, 'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs': {'R1_C': np.float64(-0.02039362806224093), 'R2_B': np.float64(-0.1), 'R2_C': np.float64(-0.1), 'R1_A': np.float64(-0.01538698570803064), 'R2_A': np.float64(-0.010511101911279047), 'R1_B': np.float64(-0.026983937243229922)}, 'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(0.0), 'R2_C': np.float64(0.0), 'R1_A': np.float64(0.0), 'R2_A': np.float64(0.0), 'R1_B': np.float64(0.0)}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs': {'R1_C': np.float64(-0.1), 'R2_B': np.float64(-0.1), 'R2_C': np.float64(-0.1), 'R1_A': np.float64(-0.1), 'R2_A': np.float64(-0.1), 'R1_B': np.float64(-0.1)}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs': {'R1_C': np.float64(-0.1), 'R2_B': np.float64(-0.1), 'R2_C': np.float64(-0.1), 'R1_A': np.float64(-0.1), 'R2_A': np.float64(-0.09140306704143274), 'R1_B': np.float64(-0.1)}, 'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs': {'R1_C': np.float64(-0.1), 'R2_B': np.float64(-0.1), 'R2_C': np.float64(-0.1), 'R1_A': np.float64(-0.1), 'R2_A': np.float64(-0.1), 'R1_B': np.float64(-0.1)}, 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs': {'R1_C': np.float64(-0.1), 'R2_B': np.float64(-0.1), 'R2_C': np.float64(-0.1), 'R1_A': np.float64(-0.1), 'R2_A': np.float64(-0.1), 'R1_B': np.float64(-0.1)}}, 'TMA': {'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs': {'R1_C': np.float64(0.09999999999999999), 'R2_B': np.float64(0.09999999999999999), 'R2_C': np.float64(0.09999999999999999), 'R1_A': np.float64(0.09999999999999999), 'R2_A': np.float64(0.09999999999999999), 'R1_B': np.float64(0.09999999999999999)}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs': {'R1_C': np.float64(0.00415517394433043), 'R2_B': np.float64(0.016189024367368715), 'R2_C': np.float64(0.05926088115745989), 'R1_A': np.float64(0.0032669694747640434), 'R2_A': np.float64(0.0014186011073906285), 'R1_B': np.float64(0.004386622400931966)}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(0.0), 'R2_C': np.float64(0.0), 'R1_A': np.float64(0.0), 'R2_A': np.float64(0.0), 'R1_B': np.float64(0.0)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs': {'R1_C': np.float64(0.09999999999999998), 'R2_B': np.float64(0.09999999999999998), 'R2_C': np.float64(0.09999999999999998), 'R1_A': np.float64(0.09999999999999998), 'R2_A': np.float64(0.09999999999999998), 'R1_B': np.float64(0.09999999999999998)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs': {'R1_C': np.float64(0.04049835052607321), 'R2_B': np.float64(0.1), 'R2_C': np.float64(0.1), 'R1_A': np.float64(0.04309736398630432), 'R2_A': np.float64(0.026673668070031162), 'R1_B': np.float64(0.03873150159343812)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs': {'R1_C': np.float64(0.09999999999999996), 'R2_B': np.float64(0.1), 'R2_C': np.float64(0.1), 'R1_A': np.float64(0.1), 'R2_A': np.float64(0.1), 'R1_B': np.float64(0.1)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(0.0), 'R2_C': np.float64(0.0), 'R1_A': np.float64(0.0), 'R2_A': np.float64(0.0), 'R1_B': np.float64(0.0)}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(0.0), 'R2_C': np.float64(0.0), 'R1_A': np.float64(0.0), 'R2_A': np.float64(0.0), 'R1_B': np.float64(0.0)}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs': {'R1_C': np.float64(0.010620624142021315), 'R2_B': np.float64(0.017376308062454462), 'R2_C': np.float64(0.05538545280165857), 'R1_A': np.float64(0.02773869657846126), 'R2_A': np.float64(0.004749024914201656), 'R1_B': np.float64(0.04266112324083004)}, 'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs': {'R1_C': np.float64(0.02039362806224093), 'R2_B': np.float64(0.1), 'R2_C': np.float64(0.10000000000000002), 'R1_A': np.float64(0.015386985708030645), 'R2_A': np.float64(0.01051110191127906), 'R1_B': np.float64(0.02698393724322992)}, 'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(0.0), 'R2_C': np.float64(0.0), 'R1_A': np.float64(0.0), 'R2_A': np.float64(0.0), 'R1_B': np.float64(0.0)}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs': {'R1_C': np.float64(0.09999999999999999), 'R2_B': np.float64(0.09999999999999998), 'R2_C': np.float64(0.09999999999999998), 'R1_A': np.float64(0.09999999999999998), 'R2_A': np.float64(0.09999999999999998), 'R1_B': np.float64(0.09999999999999998)}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs': {'R1_C': np.float64(0.09999999999999998), 'R2_B': np.float64(0.09999999999999998), 'R2_C': np.float64(0.09999999999999998), 'R1_A': np.float64(0.09999999999999998), 'R2_A': np.float64(0.09140306704143271), 'R1_B': np.float64(0.09999999999999998)}, 'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs': {'R1_C': np.float64(0.09999999999999998), 'R2_B': np.float64(0.09999999999999998), 'R2_C': np.float64(0.09999999999999998), 'R1_A': np.float64(0.09999999999999998), 'R2_A': np.float64(0.09999999999999998), 'R1_B': np.float64(0.09999999999999998)}, 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs': {'R1_C': np.float64(0.1), 'R2_B': np.float64(0.1), 'R2_C': np.float64(0.1), 'R1_A': np.float64(0.1), 'R2_A': np.float64(0.1), 'R1_B': np.float64(0.1)}}, 'methane': {'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs': {}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs': {}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs': {}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs': {}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs': {}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs': {}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs': {}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs': {}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs': {}, 'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs': {}, 'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs': {}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs': {}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs': {}, 'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs': {}, 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs': {}}, 'lactate': {'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs': {'R1_C': np.float64(-0.011662526176564887), 'R2_B': np.float64(-0.11326810129521071), 'R2_C': np.float64(-0.18607899068481676), 'R1_A': np.float64(-0.013230461708644677), 'R2_A': np.float64(0.0), 'R1_B': np.float64(-0.006491693291684842)}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs': {'R1_C': np.float64(-2.7696619543637135), 'R2_B': np.float64(-10.79091404337583), 'R2_C': np.float64(-39.500779058296395), 'R1_A': np.float64(-2.17762268957782), 'R2_A': np.float64(-0.9455790703820568), 'R1_B': np.float64(-2.923935636579519)}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs': {'R1_C': np.float64(-0.06701737556650532), 'R2_B': np.float64(-0.04721343441407818), 'R2_C': np.float64(-0.09281105592014462), 'R1_A': np.float64(-0.09362834445722439), 'R2_A': np.float64(-0.00048503250127746527), 'R1_B': np.float64(-0.06870555069085246)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs': {'R1_C': np.float64(-0.07664354612157769), 'R2_B': np.float64(-1.1021569232047885), 'R2_C': np.float64(-2.2346700387364167), 'R1_A': np.float64(-0.06606904005511163), 'R2_A': np.float64(-0.00013263035935846384), 'R1_B': np.float64(-0.016273308812607606)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs': {'R1_C': np.float64(-0.04325570720521295), 'R2_B': np.float64(-0.32967928895331156), 'R2_C': np.float64(-0.8269215152445607), 'R1_A': np.float64(-0.03629055405027021), 'R2_A': np.float64(-0.02246082134216034), 'R1_B': np.float64(-0.032614237206514005)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(-0.0026878366793308913), 'R2_C': np.float64(-0.00704304490940093), 'R1_A': np.float64(0.0), 'R2_A': np.float64(0.0), 'R1_B': np.float64(0.0)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs': {'R1_C': np.float64(-0.04450548482458164), 'R2_B': np.float64(-0.014765744622372666), 'R2_C': np.float64(-0.024137307464676277), 'R1_A': np.float64(-0.0507432936665789), 'R2_A': np.float64(-0.0003126084704341648), 'R1_B': np.float64(-0.07414224908928516)}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs': {'R1_C': np.float64(-0.013493001366517144), 'R2_B': np.float64(-0.10700084440023282), 'R2_C': np.float64(-0.20250380099079288), 'R1_A': np.float64(-0.009476200917462), 'R2_A': np.float64(-0.006619184089559963), 'R1_B': np.float64(-0.015174552858336651)}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs': {'R1_C': np.float64(-0.007399275450993545), 'R2_B': np.float64(-0.012105888312789002), 'R2_C': np.float64(-0.03858645365633649), 'R1_A': np.float64(-0.019325253760133988), 'R2_A': np.float64(-0.0033085949558065556), 'R1_B': np.float64(-0.029721549099807183)}, 'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs': {'R1_C': np.float64(-0.014622874189782581), 'R2_B': np.float64(-0.16233984792197698), 'R2_C': np.float64(-0.45171852648335314), 'R1_A': np.float64(-0.01103295379722601), 'R2_A': np.float64(-0.007536791412274483), 'R1_B': np.float64(-0.019348333619132446)}, 'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs': {'R1_C': np.float64(-0.06021723815618566), 'R2_B': np.float64(-0.045501008941393174), 'R2_C': np.float64(-0.09392673467006991), 'R1_A': np.float64(-0.05328186819481229), 'R2_A': np.float64(-0.0013303833392577932), 'R1_B': np.float64(-0.046317913343311756)}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs': {'R1_C': np.float64(-0.13883252834141374), 'R2_B': np.float64(-0.6381136551861526), 'R2_C': np.float64(-1.6647185658341215), 'R1_A': np.float64(-0.3806897478361765), 'R2_A': np.float64(-0.011113124074267088), 'R1_B': np.float64(-0.15840820437459036)}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs': {'R1_C': np.float64(-0.31650730206634603), 'R2_B': np.float64(-0.4974322664919271), 'R2_C': np.float64(-0.5290750736408264), 'R1_A': np.float64(-0.25633787902038185), 'R2_A': np.float64(-0.06997501292276205), 'R1_B': np.float64(-0.11106212972707583)}, 'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs': {'R1_C': np.float64(-0.06666326018475854), 'R2_B': np.float64(-0.290308407312205), 'R2_C': np.float64(-0.5675102409306367), 'R1_A': np.float64(-0.05915067768691815), 'R2_A': np.float64(-0.007923328731048646), 'R1_B': np.float64(-0.047096274061378514)}, 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(-2.8835764421193035), 'R2_C': np.float64(-3.8350698700093413), 'R1_A': np.float64(0.0), 'R2_A': np.float64(-0.19304605304146427), 'R1_B': np.float64(0.0)}}, 'phosphate': {'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs': {'R1_C': np.float64(-0.0005574324466063905), 'R2_B': np.float64(-0.0029014855440902), 'R2_C': np.float64(-0.004581241705157206), 'R1_A': np.float64(-0.0005936049119098031), 'R2_A': np.float64(-0.00019317632137537584), 'R1_B': np.float64(-0.00043814069685490993)}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs': {'R1_C': np.float64(-0.001994896623599438), 'R2_B': np.float64(-0.007772341298478591), 'R2_C': np.float64(-0.02845111499941112), 'R1_A': np.float64(-0.0015684701680200774), 'R2_A': np.float64(-0.0006810695766982843), 'R1_B': np.float64(-0.0021060148946555495)}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs': {'R1_C': np.float64(-0.0010770907374528173), 'R2_B': np.float64(-0.0007588054957519733), 'R2_C': np.float64(-0.0014916419483720388), 'R1_A': np.float64(-0.0015047772570241785), 'R2_A': np.float64(-7.79535172891134e-06), 'R1_B': np.float64(-0.0011042227726043145)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs': {'R1_C': np.float64(-0.0015222971375675142), 'R2_B': np.float64(-0.010598825056064083), 'R2_C': np.float64(-0.02048088159063444), 'R1_A': np.float64(-0.0013995921833846344), 'R2_A': np.float64(-0.0003215426990473287), 'R1_B': np.float64(-0.0008217701194298942)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs': {'R1_C': np.float64(-0.0027849178966335166), 'R2_B': np.float64(-0.021225632668537036), 'R2_C': np.float64(-0.05323941453530227), 'R1_A': np.float64(-0.002336482743742401), 'R2_A': np.float64(-0.0014460876349129779), 'R1_B': np.float64(-0.0020997916517830175)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs': {'R1_C': np.float64(-0.00019460513233706398), 'R2_B': np.float64(-0.0005122467834939206), 'R2_C': np.float64(-0.0005954109050615349), 'R1_A': np.float64(-0.0003697807957552696), 'R2_A': np.float64(-2.5069115308303805e-06), 'R1_B': np.float64(-0.00010764520445513342)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs': {'R1_C': np.float64(-0.0007093745328183139), 'R2_B': np.float64(-0.00023535173775761676), 'R2_C': np.float64(-0.00038472541696232707), 'R1_A': np.float64(-0.00080879919363364), 'R2_A': np.float64(-4.982677720360856e-06), 'R1_B': np.float64(-0.0011817559906855135)}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs': {'R1_C': np.float64(-0.0032220409759517953), 'R2_B': np.float64(-0.025551105773584874), 'R2_C': np.float64(-0.048356590713572054), 'R1_A': np.float64(-0.0022628551515737993), 'R2_A': np.float64(-0.0015806181133913827), 'R1_B': np.float64(-0.0036235845363979794)}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs': {'R1_C': np.float64(-0.00013383400753433488), 'R2_B': np.float64(-0.00021896462138682338), 'R2_C': np.float64(-0.0006979304613766145), 'R1_A': np.float64(-0.0003495445161443979), 'R2_A': np.float64(-5.984403813267218e-05), 'R1_B': np.float64(-0.0005375869640887625)}, 'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs': {'R1_C': np.float64(-0.0004176377939161579), 'R2_B': np.float64(-0.003336019738252453), 'R2_C': np.float64(-0.006270226025070721), 'R1_A': np.float64(-0.00031510757901981546), 'R2_A': np.float64(-0.00021525514736543664), 'R1_B': np.float64(-0.0005525996643186138)}, 'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs': {'R1_C': np.float64(-0.0009595300625384771), 'R2_B': np.float64(-0.0007250346793032699), 'R2_C': np.float64(-0.0014966731844833167), 'R1_A': np.float64(-0.0008490185848200567), 'R2_A': np.float64(-2.1198959763104267e-05), 'R1_B': np.float64(-0.0007380516218908814)}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs': {'R1_C': np.float64(-0.002489727331273735), 'R2_B': np.float64(-0.009072670844306756), 'R2_C': np.float64(-0.022608295913302545), 'R1_A': np.float64(-0.00567857691630589), 'R2_A': np.float64(-0.0008057669758684097), 'R1_B': np.float64(-0.0027478295559068494)}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs': {'R1_C': np.float64(-0.005065370845052625), 'R2_B': np.float64(-0.007265697310873007), 'R2_C': np.float64(-0.007650522623639323), 'R1_A': np.float64(-0.004333617908793496), 'R2_A': np.float64(-0.0019626062189686787), 'R1_B': np.float64(-0.0025668408527540054)}, 'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs': {'R1_C': np.float64(-0.0012247442768691932), 'R2_B': np.float64(-0.0029361981502896077), 'R2_C': np.float64(-0.005021522206546155), 'R1_A': np.float64(-0.0011047259995801456), 'R2_A': np.float64(-0.00028633656112327224), 'R1_B': np.float64(-0.0009121492448756099)}, 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs': {'R1_C': np.float64(-0.0005800686847693273), 'R2_B': np.float64(-0.016909455614019468), 'R2_C': np.float64(-0.021788142206388333), 'R1_A': np.float64(-0.00023733366392565312), 'R2_A': np.float64(-0.0031140323869078643), 'R1_B': np.float64(-0.00034958505215825457)}}}
    _betaine_df = DataFrame(fluxes_2['betaine'])
    _betaine_df.T.sum()
    return


@app.cell
def _(DataFrame, display, np):
    from json import load
    _dic = {'betaine': {'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs': {'R1_C': np.float64(-0.1), 'R2_B': np.float64(-0.1), 'R2_C': np.float64(-0.1), 'R1_A': np.float64(-0.1), 'R2_A': np.float64(-0.1), 'R1_B': np.float64(-0.1)}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs': {'R1_C': np.float64(-0.004155173944330431), 'R2_B': np.float64(-0.016189024367368715), 'R2_C': np.float64(-0.059260881157459894), 'R1_A': np.float64(-0.003266969474764042), 'R2_A': np.float64(-0.0014186011073906296), 'R1_B': np.float64(-0.004386622400931961)}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(0.0), 'R2_C': np.float64(0.0), 'R1_A': np.float64(0.0), 'R2_A': np.float64(0.0), 'R1_B': np.float64(0.0)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs': {'R1_C': np.float64(-0.1), 'R2_B': np.float64(-0.1), 'R2_C': np.float64(-0.1), 'R1_A': np.float64(-0.1), 'R2_A': np.float64(-0.1), 'R1_B': np.float64(-0.1)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs': {'R1_C': np.float64(-0.04049835052607321), 'R2_B': np.float64(-0.1), 'R2_C': np.float64(-0.1), 'R1_A': np.float64(-0.04309736398630432), 'R2_A': np.float64(-0.02667366807003116), 'R1_B': np.float64(-0.03873150159343811)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs': {'R1_C': np.float64(-0.1), 'R2_B': np.float64(-0.1), 'R2_C': np.float64(-0.1), 'R1_A': np.float64(-0.1), 'R2_A': np.float64(-0.1), 'R1_B': np.float64(-0.1)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(0.0), 'R2_C': np.float64(0.0), 'R1_A': np.float64(0.0), 'R2_A': np.float64(0.0), 'R1_B': np.float64(0.0)}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(0.0), 'R2_C': np.float64(0.0), 'R1_A': np.float64(0.0), 'R2_A': np.float64(0.0), 'R1_B': np.float64(0.0)}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs': {'R1_C': np.float64(-0.010620624142021317), 'R2_B': np.float64(-0.017376308062454465), 'R2_C': np.float64(-0.05538545280165857), 'R1_A': np.float64(-0.02773869657846126), 'R2_A': np.float64(-0.004749024914201657), 'R1_B': np.float64(-0.04266112324083004)}, 'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs': {'R1_C': np.float64(-0.02039362806224093), 'R2_B': np.float64(-0.1), 'R2_C': np.float64(-0.1), 'R1_A': np.float64(-0.01538698570803064), 'R2_A': np.float64(-0.010511101911279047), 'R1_B': np.float64(-0.026983937243229922)}, 'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(0.0), 'R2_C': np.float64(0.0), 'R1_A': np.float64(0.0), 'R2_A': np.float64(0.0), 'R1_B': np.float64(0.0)}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs': {'R1_C': np.float64(-0.1), 'R2_B': np.float64(-0.1), 'R2_C': np.float64(-0.1), 'R1_A': np.float64(-0.1), 'R2_A': np.float64(-0.1), 'R1_B': np.float64(-0.1)}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs': {'R1_C': np.float64(-0.1), 'R2_B': np.float64(-0.1), 'R2_C': np.float64(-0.1), 'R1_A': np.float64(-0.1), 'R2_A': np.float64(-0.09140306704143274), 'R1_B': np.float64(-0.1)}, 'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs': {'R1_C': np.float64(-0.1), 'R2_B': np.float64(-0.1), 'R2_C': np.float64(-0.1), 'R1_A': np.float64(-0.1), 'R2_A': np.float64(-0.1), 'R1_B': np.float64(-0.1)}, 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs': {'R1_C': np.float64(-0.1), 'R2_B': np.float64(-0.1), 'R2_C': np.float64(-0.1), 'R1_A': np.float64(-0.1), 'R2_A': np.float64(-0.1), 'R1_B': np.float64(-0.1)}}, 'TMA': {'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs': {'R1_C': np.float64(0.09999999999999999), 'R2_B': np.float64(0.09999999999999999), 'R2_C': np.float64(0.09999999999999999), 'R1_A': np.float64(0.09999999999999999), 'R2_A': np.float64(0.09999999999999999), 'R1_B': np.float64(0.09999999999999999)}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs': {'R1_C': np.float64(0.00415517394433043), 'R2_B': np.float64(0.016189024367368715), 'R2_C': np.float64(0.05926088115745989), 'R1_A': np.float64(0.0032669694747640434), 'R2_A': np.float64(0.0014186011073906285), 'R1_B': np.float64(0.004386622400931966)}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(0.0), 'R2_C': np.float64(0.0), 'R1_A': np.float64(0.0), 'R2_A': np.float64(0.0), 'R1_B': np.float64(0.0)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs': {'R1_C': np.float64(0.09999999999999998), 'R2_B': np.float64(0.09999999999999998), 'R2_C': np.float64(0.09999999999999998), 'R1_A': np.float64(0.09999999999999998), 'R2_A': np.float64(0.09999999999999998), 'R1_B': np.float64(0.09999999999999998)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs': {'R1_C': np.float64(0.04049835052607321), 'R2_B': np.float64(0.1), 'R2_C': np.float64(0.1), 'R1_A': np.float64(0.04309736398630432), 'R2_A': np.float64(0.026673668070031162), 'R1_B': np.float64(0.03873150159343812)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs': {'R1_C': np.float64(0.09999999999999996), 'R2_B': np.float64(0.1), 'R2_C': np.float64(0.1), 'R1_A': np.float64(0.1), 'R2_A': np.float64(0.1), 'R1_B': np.float64(0.1)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(0.0), 'R2_C': np.float64(0.0), 'R1_A': np.float64(0.0), 'R2_A': np.float64(0.0), 'R1_B': np.float64(0.0)}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(0.0), 'R2_C': np.float64(0.0), 'R1_A': np.float64(0.0), 'R2_A': np.float64(0.0), 'R1_B': np.float64(0.0)}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs': {'R1_C': np.float64(0.010620624142021315), 'R2_B': np.float64(0.017376308062454462), 'R2_C': np.float64(0.05538545280165857), 'R1_A': np.float64(0.02773869657846126), 'R2_A': np.float64(0.004749024914201656), 'R1_B': np.float64(0.04266112324083004)}, 'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs': {'R1_C': np.float64(0.02039362806224093), 'R2_B': np.float64(0.1), 'R2_C': np.float64(0.10000000000000002), 'R1_A': np.float64(0.015386985708030645), 'R2_A': np.float64(0.01051110191127906), 'R1_B': np.float64(0.02698393724322992)}, 'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(0.0), 'R2_C': np.float64(0.0), 'R1_A': np.float64(0.0), 'R2_A': np.float64(0.0), 'R1_B': np.float64(0.0)}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs': {'R1_C': np.float64(0.09999999999999999), 'R2_B': np.float64(0.09999999999999998), 'R2_C': np.float64(0.09999999999999998), 'R1_A': np.float64(0.09999999999999998), 'R2_A': np.float64(0.09999999999999998), 'R1_B': np.float64(0.09999999999999998)}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs': {'R1_C': np.float64(0.09999999999999998), 'R2_B': np.float64(0.09999999999999998), 'R2_C': np.float64(0.09999999999999998), 'R1_A': np.float64(0.09999999999999998), 'R2_A': np.float64(0.09140306704143271), 'R1_B': np.float64(0.09999999999999998)}, 'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs': {'R1_C': np.float64(0.09999999999999998), 'R2_B': np.float64(0.09999999999999998), 'R2_C': np.float64(0.09999999999999998), 'R1_A': np.float64(0.09999999999999998), 'R2_A': np.float64(0.09999999999999998), 'R1_B': np.float64(0.09999999999999998)}, 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs': {'R1_C': np.float64(0.1), 'R2_B': np.float64(0.1), 'R2_C': np.float64(0.1), 'R1_A': np.float64(0.1), 'R2_A': np.float64(0.1), 'R1_B': np.float64(0.1)}}, 'methane': {'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs': {}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs': {}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs': {}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs': {}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs': {}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs': {}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs': {}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs': {}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs': {}, 'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs': {}, 'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs': {}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs': {}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs': {}, 'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs': {}, 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs': {}}, 'lactate': {'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs': {'R1_C': np.float64(-0.011662526176564887), 'R2_B': np.float64(-0.11326810129521071), 'R2_C': np.float64(-0.18607899068481676), 'R1_A': np.float64(-0.013230461708644677), 'R2_A': np.float64(0.0), 'R1_B': np.float64(-0.006491693291684842)}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs': {'R1_C': np.float64(-2.7696619543637135), 'R2_B': np.float64(-10.79091404337583), 'R2_C': np.float64(-39.500779058296395), 'R1_A': np.float64(-2.17762268957782), 'R2_A': np.float64(-0.9455790703820568), 'R1_B': np.float64(-2.923935636579519)}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs': {'R1_C': np.float64(-0.06701737556650532), 'R2_B': np.float64(-0.04721343441407818), 'R2_C': np.float64(-0.09281105592014462), 'R1_A': np.float64(-0.09362834445722439), 'R2_A': np.float64(-0.00048503250127746527), 'R1_B': np.float64(-0.06870555069085246)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs': {'R1_C': np.float64(-0.07664354612157769), 'R2_B': np.float64(-1.1021569232047885), 'R2_C': np.float64(-2.2346700387364167), 'R1_A': np.float64(-0.06606904005511163), 'R2_A': np.float64(-0.00013263035935846384), 'R1_B': np.float64(-0.016273308812607606)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs': {'R1_C': np.float64(-0.04325570720521295), 'R2_B': np.float64(-0.32967928895331156), 'R2_C': np.float64(-0.8269215152445607), 'R1_A': np.float64(-0.03629055405027021), 'R2_A': np.float64(-0.02246082134216034), 'R1_B': np.float64(-0.032614237206514005)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(-0.0026878366793308913), 'R2_C': np.float64(-0.00704304490940093), 'R1_A': np.float64(0.0), 'R2_A': np.float64(0.0), 'R1_B': np.float64(0.0)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs': {'R1_C': np.float64(-0.04450548482458164), 'R2_B': np.float64(-0.014765744622372666), 'R2_C': np.float64(-0.024137307464676277), 'R1_A': np.float64(-0.0507432936665789), 'R2_A': np.float64(-0.0003126084704341648), 'R1_B': np.float64(-0.07414224908928516)}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs': {'R1_C': np.float64(-0.013493001366517144), 'R2_B': np.float64(-0.10700084440023282), 'R2_C': np.float64(-0.20250380099079288), 'R1_A': np.float64(-0.009476200917462), 'R2_A': np.float64(-0.006619184089559963), 'R1_B': np.float64(-0.015174552858336651)}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs': {'R1_C': np.float64(-0.007399275450993545), 'R2_B': np.float64(-0.012105888312789002), 'R2_C': np.float64(-0.03858645365633649), 'R1_A': np.float64(-0.019325253760133988), 'R2_A': np.float64(-0.0033085949558065556), 'R1_B': np.float64(-0.029721549099807183)}, 'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs': {'R1_C': np.float64(-0.014622874189782581), 'R2_B': np.float64(-0.16233984792197698), 'R2_C': np.float64(-0.45171852648335314), 'R1_A': np.float64(-0.01103295379722601), 'R2_A': np.float64(-0.007536791412274483), 'R1_B': np.float64(-0.019348333619132446)}, 'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs': {'R1_C': np.float64(-0.06021723815618566), 'R2_B': np.float64(-0.045501008941393174), 'R2_C': np.float64(-0.09392673467006991), 'R1_A': np.float64(-0.05328186819481229), 'R2_A': np.float64(-0.0013303833392577932), 'R1_B': np.float64(-0.046317913343311756)}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs': {'R1_C': np.float64(-0.13883252834141374), 'R2_B': np.float64(-0.6381136551861526), 'R2_C': np.float64(-1.6647185658341215), 'R1_A': np.float64(-0.3806897478361765), 'R2_A': np.float64(-0.011113124074267088), 'R1_B': np.float64(-0.15840820437459036)}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs': {'R1_C': np.float64(-0.31650730206634603), 'R2_B': np.float64(-0.4974322664919271), 'R2_C': np.float64(-0.5290750736408264), 'R1_A': np.float64(-0.25633787902038185), 'R2_A': np.float64(-0.06997501292276205), 'R1_B': np.float64(-0.11106212972707583)}, 'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs': {'R1_C': np.float64(-0.06666326018475854), 'R2_B': np.float64(-0.290308407312205), 'R2_C': np.float64(-0.5675102409306367), 'R1_A': np.float64(-0.05915067768691815), 'R2_A': np.float64(-0.007923328731048646), 'R1_B': np.float64(-0.047096274061378514)}, 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs': {'R1_C': np.float64(0.0), 'R2_B': np.float64(-2.8835764421193035), 'R2_C': np.float64(-3.8350698700093413), 'R1_A': np.float64(0.0), 'R2_A': np.float64(-0.19304605304146427), 'R1_B': np.float64(0.0)}}, 'phosphate': {'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs': {'R1_C': np.float64(-0.0005574324466063905), 'R2_B': np.float64(-0.0029014855440902), 'R2_C': np.float64(-0.004581241705157206), 'R1_A': np.float64(-0.0005936049119098031), 'R2_A': np.float64(-0.00019317632137537584), 'R1_B': np.float64(-0.00043814069685490993)}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs': {'R1_C': np.float64(-0.001994896623599438), 'R2_B': np.float64(-0.007772341298478591), 'R2_C': np.float64(-0.02845111499941112), 'R1_A': np.float64(-0.0015684701680200774), 'R2_A': np.float64(-0.0006810695766982843), 'R1_B': np.float64(-0.0021060148946555495)}, 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs': {'R1_C': np.float64(-0.0010770907374528173), 'R2_B': np.float64(-0.0007588054957519733), 'R2_C': np.float64(-0.0014916419483720388), 'R1_A': np.float64(-0.0015047772570241785), 'R2_A': np.float64(-7.79535172891134e-06), 'R1_B': np.float64(-0.0011042227726043145)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs': {'R1_C': np.float64(-0.0015222971375675142), 'R2_B': np.float64(-0.010598825056064083), 'R2_C': np.float64(-0.02048088159063444), 'R1_A': np.float64(-0.0013995921833846344), 'R2_A': np.float64(-0.0003215426990473287), 'R1_B': np.float64(-0.0008217701194298942)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs': {'R1_C': np.float64(-0.0027849178966335166), 'R2_B': np.float64(-0.021225632668537036), 'R2_C': np.float64(-0.05323941453530227), 'R1_A': np.float64(-0.002336482743742401), 'R2_A': np.float64(-0.0014460876349129779), 'R1_B': np.float64(-0.0020997916517830175)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs': {'R1_C': np.float64(-0.00019460513233706398), 'R2_B': np.float64(-0.0005122467834939206), 'R2_C': np.float64(-0.0005954109050615349), 'R1_A': np.float64(-0.0003697807957552696), 'R2_A': np.float64(-2.5069115308303805e-06), 'R1_B': np.float64(-0.00010764520445513342)}, 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs': {'R1_C': np.float64(-0.0007093745328183139), 'R2_B': np.float64(-0.00023535173775761676), 'R2_C': np.float64(-0.00038472541696232707), 'R1_A': np.float64(-0.00080879919363364), 'R2_A': np.float64(-4.982677720360856e-06), 'R1_B': np.float64(-0.0011817559906855135)}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs': {'R1_C': np.float64(-0.0032220409759517953), 'R2_B': np.float64(-0.025551105773584874), 'R2_C': np.float64(-0.048356590713572054), 'R1_A': np.float64(-0.0022628551515737993), 'R2_A': np.float64(-0.0015806181133913827), 'R1_B': np.float64(-0.0036235845363979794)}, 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs': {'R1_C': np.float64(-0.00013383400753433488), 'R2_B': np.float64(-0.00021896462138682338), 'R2_C': np.float64(-0.0006979304613766145), 'R1_A': np.float64(-0.0003495445161443979), 'R2_A': np.float64(-5.984403813267218e-05), 'R1_B': np.float64(-0.0005375869640887625)}, 'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs': {'R1_C': np.float64(-0.0004176377939161579), 'R2_B': np.float64(-0.003336019738252453), 'R2_C': np.float64(-0.006270226025070721), 'R1_A': np.float64(-0.00031510757901981546), 'R2_A': np.float64(-0.00021525514736543664), 'R1_B': np.float64(-0.0005525996643186138)}, 'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs': {'R1_C': np.float64(-0.0009595300625384771), 'R2_B': np.float64(-0.0007250346793032699), 'R2_C': np.float64(-0.0014966731844833167), 'R1_A': np.float64(-0.0008490185848200567), 'R2_A': np.float64(-2.1198959763104267e-05), 'R1_B': np.float64(-0.0007380516218908814)}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs': {'R1_C': np.float64(-0.002489727331273735), 'R2_B': np.float64(-0.009072670844306756), 'R2_C': np.float64(-0.022608295913302545), 'R1_A': np.float64(-0.00567857691630589), 'R2_A': np.float64(-0.0008057669758684097), 'R1_B': np.float64(-0.0027478295559068494)}, 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs': {'R1_C': np.float64(-0.005065370845052625), 'R2_B': np.float64(-0.007265697310873007), 'R2_C': np.float64(-0.007650522623639323), 'R1_A': np.float64(-0.004333617908793496), 'R2_A': np.float64(-0.0019626062189686787), 'R1_B': np.float64(-0.0025668408527540054)}, 'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs': {'R1_C': np.float64(-0.0012247442768691932), 'R2_B': np.float64(-0.0029361981502896077), 'R2_C': np.float64(-0.005021522206546155), 'R1_A': np.float64(-0.0011047259995801456), 'R2_A': np.float64(-0.00028633656112327224), 'R1_B': np.float64(-0.0009121492448756099)}, 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs': {'R1_C': np.float64(-0.0005800686847693273), 'R2_B': np.float64(-0.016909455614019468), 'R2_C': np.float64(-0.021788142206388333), 'R1_A': np.float64(-0.00023733366392565312), 'R2_A': np.float64(-0.0031140323869078643), 'R1_B': np.float64(-0.00034958505215825457)}}}
    df = DataFrame(_dic['betaine'])
    phylo = load(open('../data/Saltern_phylogeny.json', 'r'))
    _d = dict(sorted(df.sum().to_dict().items(), key=lambda item: item[1]))
    display(_d)
    {phylo.get(_k + '__')['Family'] + ' ' + phylo.get(_k + '__')['Genus']: _v for _k, _v in _d.items()}
    return (load,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## plot the correlational and modeling v experimental figures
    """)
    return


@app.cell
def _(
    DataFrame,
    LinearRegression,
    Pipeline,
    PolynomialFeatures,
    averaged_unrestored_normalized_abundances,
    averaged_unrestored_normalized_methane,
    betaine_striped,
    fluxes_to_emissions,
    json,
    mepe_striped,
    methane_dic,
    np,
    r2_score,
):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogLocator, LogFormatterMathtext, NullFormatter
    from adjustText import adjust_text
    from matplotlib.ticker import FuncFormatter
    from scipy.stats import spearmanr

    def minor_log_formatter(x, pos):
        exp = np.log10(x)
        if abs(exp - round(exp)) < 1e-09:
            return ''
        return f'$10^{{{exp:.1f}}}$'
    _log = False
    files = ['saltern_fluxes.json', 'saltern_fluxes_betaine.json']
    palette = [{'scatter': '#e8671b', 'trendline': '#7f2704', 'label': '#e8671b'}, {'scatter': '#1f77b4', 'trendline': '#0a3d6b', 'label': '#1f77b4'}]
    _fig, _ax = plt.subplots(figsize=(10, 8), dpi=300)
    all_texts = []
    annotations = {}
    all_pt_x = []
    all_pt_y = []
    _label_size = 11
    for _fluxes_file, _color in zip(files, palette):
        fluxes_3 = json.load(open(_fluxes_file))
        betaine_1 = 'betaine' in _fluxes_file
        if betaine_1:
            _fluxes_in_emissions = fluxes_to_emissions(DataFrame(fluxes_3['TMA']).T.sum()).to_dict()
        else:
            _fluxes_in_emissions = fluxes_to_emissions(DataFrame(fluxes_3['0.07']['methane']).T.sum()).to_dict()
        _combined = {_ID: {'fluxes': emissions, 'methane': methane_dic[_ID]} for _ID, emissions in _fluxes_in_emissions.items()}
        _combined_df = DataFrame(_combined).T
        _methane_ary = np.array(list(_combined_df['methane'].values))
        _emission_fluxes = np.array(list(_combined_df['fluxes'].values))
        if betaine_1:
            _emission_fluxes = _emission_fluxes * 0.06
        _xs = np.array(_methane_ary).reshape(-1, 1)
        _ys = np.array(_emission_fluxes)
        xs_flat = _xs.flatten()
        ys_flat = _ys.flatten()
        labels = list(_combined_df.index)
        all_pt_x.extend(xs_flat.tolist())
        all_pt_y.extend(ys_flat.tolist())
        _mdl = LinearRegression()
        _mdl.fit(_xs, _ys)
        _y_pred = _mdl.predict(_xs)
        dataset_label = 'Betaine' if betaine_1 else 'MePn'
        _ax.scatter(xs_flat, ys_flat, color=_color['scatter'], zorder=3, s=60, label=f'{dataset_label} data')
        rng = np.random.default_rng(42)
        for _xi, _yi, lbl in zip(xs_flat, ys_flat, labels):
            angle = rng.uniform(0, 2 * np.pi)
            x_off = _xi * (1 + 0.15 * np.cos(angle))
            y_off = _yi * (1 + 0.15 * np.sin(angle))
            all_texts.append(_ax.text(x_off, y_off, lbl, fontsize=_label_size, ha='center', va='bottom', color=_color['label'], zorder=5))
        _sorted_idx = np.argsort(_methane_ary)
        _ax.plot(_xs[_sorted_idx].flatten(), _y_pred[_sorted_idx], color=_color['trendline'], linestyle='--', linewidth=2, label=f'{dataset_label} trendline')
        anchor_name = 'R2_A' if betaine_1 else 'R1_C'
        if anchor_name in labels:
            idx = labels.index(anchor_name)
            anchor_x = float(xs_flat[idx])
            anchor_y = float(ys_flat[idx])
            anchor_y_pred = float(_mdl.predict([[anchor_x]])[0])
        else:
            idx = len(labels) // 2
            anchor_x = float(xs_flat[idx])
            anchor_y = float(ys_flat[idx])
            anchor_y_pred = float(_mdl.predict([[anchor_x]])[0])
        _eq = f'y = {_mdl.coef_[0]:.2f}x + {_mdl.intercept_:.2f}'
        _r2_text = f'$R^2$ = {r2_score(_ys, _y_pred):.2f}  |  $r_s$ = {spearmanr(_ys, _y_pred)[0]:.2f}'
        annotations[dataset_label] = dict(anchor_x=anchor_x, anchor_y_pred=anchor_y_pred, eq=_eq, r2_text=_r2_text, color=_color, below=betaine_1)
    _ax.set_xscale('log')
    _ax.set_yscale('log')
    _ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0,)))
    _ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0,)))
    _ax.xaxis.set_major_formatter(LogFormatterMathtext(base=10, labelOnlyBase=True))
    _ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10, labelOnlyBase=True))
    log_subs = 10 ** np.arange(0.1, 1.0, 0.1)
    _ax.xaxis.set_minor_locator(LogLocator(base=10, subs=log_subs, numticks=100))
    _ax.yaxis.set_minor_locator(LogLocator(base=10, subs=log_subs, numticks=100))
    _ax.xaxis.set_minor_formatter(FuncFormatter(minor_log_formatter))
    _ax.yaxis.set_minor_formatter(FuncFormatter(minor_log_formatter))
    plt.setp(_ax.get_xticklabels(), rotation=90)
    _ax.grid(True, which='both', linestyle='--', color='gray', alpha=0.3)
    adjust_text(all_texts, x=np.array(all_pt_x), y=np.array(all_pt_y), ax=_ax, expand_text=(2.0, 2.0), expand_points=(2.0, 2.0), force_text=(1.5, 1.5), force_points=(1.5, 1.5), lim=500)
    _fig.canvas.draw()
    renderer = _fig.canvas.get_renderer()
    for txt, _xi, _yi in zip(all_texts, all_pt_x, all_pt_y):
        bb = txt.get_window_extent(renderer=renderer)
        bb_data = bb.transformed(_ax.transData.inverted())
        cx = np.clip(_xi, bb_data.x0, bb_data.x1)
        cy = np.clip(_yi, bb_data.y0, bb_data.y1)
        _ax.plot([cx, _xi], [cy, _yi], color='gray', lw=0.8, zorder=1, clip_on=True)
    for dataset_label, ann in annotations.items():
        anchor_x = ann['anchor_x']
        anchor_y_pred = ann['anchor_y_pred']
        _color = ann['color']
        if ann['below']:
            text_y = anchor_y_pred / 1.5
            va = 'top'
        else:
            text_y = anchor_y_pred * 3.5
            va = 'bottom'
        _ax.text(anchor_x, text_y, ann['eq'] + '\n' + ann['r2_text'], fontsize=_label_size, color=_color['trendline'], ha='center', va=va, bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85, edgecolor=_color['trendline']), zorder=6)
    _axes_size = 18
    _ax.set_xlabel('Experimental measurement $\\left[\\frac{umol~~CH_4}{m^2 \\cdot day \\cdot \\left(\\frac{ug~DNA}{g~soil}\\right)}\\right]$', fontsize=_axes_size)
    _ax.set_ylabel('Metabolic modeling $\\left[\\frac{umol~~CH_4}{m^2 \\cdot day \\cdot \\left(\\frac{ug~DNA}{g~soil}\\right)}\\right]$', fontsize=_axes_size)
    _ax.legend(fontsize=12, loc='lower right')
    plt.tight_layout()
    _fig.savefig('combined_experiment_modeling.png', bbox_inches='tight')
    _fig.show()

    def minor_log_formatter(x, pos):
        exp = np.log10(x)
        if abs(exp - round(exp)) < 1e-09:
            return ''
        return f'$10^{{{exp:.1f}}}$'
    log_subs = 10 ** np.arange(0.1, 1.0, 0.1)
    palette = [{'scatter': '#1f77b4', 'trendline': '#0a3d6b', 'label': '#1f77b4'}, {'scatter': '#e8671b', 'trendline': '#7f2704', 'label': '#e8671b'}]
    _fig, _ax = plt.subplots(figsize=(10, 8), dpi=300)
    all_texts = []
    all_pt_x = []
    all_pt_y = []
    annotations = {}
    rng = np.random.default_rng(42)
    _linear = True
    for _samples, _color in zip(['Betaine', 'MePn'], palette):
        if _samples == 'MePn':
            _sample_totals = averaged_unrestored_normalized_abundances.loc[list(mepe_striped)].sum()
        else:
            _sample_totals = averaged_unrestored_normalized_abundances.loc[list(betaine_striped)].sum()
        _all_xs = np.array(_sample_totals.to_list())
        _xs = _all_xs.reshape(-1, 1)
        _all_ys = averaged_unrestored_normalized_methane.to_list()
        _ys = np.array(_all_ys)
        xs_flat = _all_xs.flatten()
        ys_flat = _ys.flatten()
        labels = list(_sample_totals.index)
        all_pt_x.extend(xs_flat.tolist())
        all_pt_y.extend(ys_flat.tolist())
        if _linear:
            _model = LinearRegression()
            _model.fit(_xs, _ys)
            _y_pred = _model.predict(_xs)
            _coef = _model.coef_[0]
            _intercept = _model.intercept_
        else:
            _poly = PolynomialFeatures(degree=2, include_bias=False)
            _lin = LinearRegression()
            _model = Pipeline([('poly', _poly), ('lin', _lin)])
            _model.fit(_xs, _ys)
            _y_pred = _model.predict(_xs)
            _term_names = _model.named_steps['poly'].get_feature_names_out()
            _coef = _model.named_steps['lin'].coef_
            _intercept = _model.named_steps['lin'].intercept_
        _ax.scatter(xs_flat, ys_flat, color=_color['scatter'], zorder=3, s=60, label=f'{_samples} data')
        for _xi, _yi, lbl in zip(xs_flat, ys_flat, labels):
            angle = rng.uniform(0, 2 * np.pi)
            x_off = _xi * (1 + 0.15 * np.cos(angle))
            y_off = _yi * (1 + 0.15 * np.sin(angle))
            all_texts.append(_ax.text(x_off, y_off, lbl, fontsize=_label_size, ha='center', va='bottom', color=_color['label'], zorder=5))
        _sorted_idx = np.argsort(_all_xs)
        _ax.plot(_all_xs[_sorted_idx], _y_pred[_sorted_idx], color=_color['trendline'], linestyle='--', linewidth=2, label=f'{_samples} trendline')
        if _linear:
            _eq = f'y = {_coef:.2f}x + {_intercept:.2f}'
        else:
            _precision = 4
            _parts = [f'{_intercept:.{_precision}g}']
            for _c, _t in zip(_coef, _term_names):
                _sign = ' + ' if _c >= 0 else ' - '
                _parts.append(f'{_sign}{abs(_c):.{_precision}g}*{_t}')
            _eq = 'y = ' + ''.join(_parts)
        _r2_text = f'$R^2$ = {r2_score(_all_ys, _y_pred):.2f}  |  $r_s$ = {spearmanr(_all_ys, _y_pred)[0]:.2f}'
        mid_idx = len(_sorted_idx) // 2
        anchor_x = float(_all_xs[_sorted_idx][mid_idx])
        anchor_y_pred = float(_y_pred[_sorted_idx][mid_idx])
        annotations[_samples] = dict(anchor_x=anchor_x, anchor_y_pred=anchor_y_pred, eq=_eq, r2_text=_r2_text, color=_color, below=_samples == 'Betaine')
    _ax.set_xscale('log')
    _ax.set_yscale('log')
    _ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0,)))
    _ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0,)))
    _ax.xaxis.set_minor_locator(LogLocator(base=10, subs=log_subs, numticks=100))
    _ax.yaxis.set_minor_locator(LogLocator(base=10, subs=log_subs, numticks=100))
    _ax.xaxis.set_major_formatter(LogFormatterMathtext(base=10, labelOnlyBase=True))
    _ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10, labelOnlyBase=True))
    _ax.xaxis.set_minor_formatter(FuncFormatter(minor_log_formatter))
    _ax.yaxis.set_minor_formatter(FuncFormatter(minor_log_formatter))
    plt.setp(_ax.get_xticklabels(), rotation=90)
    plt.setp(_ax.get_yticklabels(minor=True), fontsize=7, color='gray')
    _ax.grid(True, which='both', linestyle='--', color='gray', alpha=0.3)
    adjust_text(all_texts, x=np.array(all_pt_x), y=np.array(all_pt_y), ax=_ax, expand_text=(2.0, 2.0), expand_points=(2.0, 2.0), force_text=(1.5, 1.5), force_points=(1.5, 1.5), lim=500)
    _fig.canvas.draw()
    renderer = _fig.canvas.get_renderer()
    for txt, _xi, _yi in zip(all_texts, all_pt_x, all_pt_y):
        bb = txt.get_window_extent(renderer=renderer)
        bb_data = bb.transformed(_ax.transData.inverted())
        cx = np.clip(_xi, bb_data.x0, bb_data.x1)
        cy = np.clip(_yi, bb_data.y0, bb_data.y1)
        _ax.plot([cx, _xi], [cy, _yi], color='gray', lw=0.8, zorder=1, clip_on=True)
    box_positions = {'MePn': dict(x=0.72, y=0.3, va='bottom'), 'Betaine': dict(x=0.35, y=0.8, va='top')}
    for _samples, ann in annotations.items():
        _color = ann['color']
        pos = box_positions[_samples]
        _ax.text(pos['x'], pos['y'], ann['eq'] + '\n' + ann['r2_text'], fontsize=_label_size, color=_color['trendline'], ha='center', va=pos['va'], transform=_ax.transAxes, bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85, edgecolor=_color['trendline']), zorder=6)
    _axes_size = 18
    _ax.set_xlabel('Abundances $\\left[\\frac{rel.ab}{\\frac{ug~DNA}{g~soil}}\\right]$', fontsize=_axes_size)
    _ax.set_ylabel('Methane emissions $\\left[\\frac{umol~~CH_4}{m^2 \\cdot day \\cdot \\left(\\frac{ug~DNA}{g~soil}\\right)}\\right]$', fontsize=_axes_size)
    _ax.legend(fontsize=12, loc='lower right')
    plt.tight_layout()
    _fig.savefig('consumers_combined.png', bbox_inches='tight')
    _fig.show()
    return (
        LogFormatterMathtext,
        LogLocator,
        NullFormatter,
        adjust_text,
        fluxes_3,
        plt,
        spearmanr,
    )


@app.cell
def _(fluxes_3, models_1, sols_1):
    for _ID, content in sols_1.items():
        _model = models_1[_ID]
        try:
            print(list(_model.metabolites)[0])
            met = _model.metabolites.get_by_id('cpd01024_c0')
        except KeyError:
            print(f'Metabolite cpd01024_c0 not found in model {_ID}')
            continue
        metReactions = [_rxn.id for _rxn in met.reactions]
        for _sample, _sol in content.items():
            methane_flux = fluxes_3['methane'][_ID][_sample]
            if methane_flux == 0:
                continue
            tma_flux = fluxes_3['TMA'][_ID][_sample]
            betaine_flux = fluxes_3['betaine'][_ID][_sample]
            lactate_flux = fluxes_3['lactate'][_ID][_sample]
            print(_ID, _sample)
            print(methane_flux, tma_flux, betaine_flux, lactate_flux)
            for _rxn in metReactions:
                flux = _sol.fluxes[_rxn]
                stoich = _sol.model.reactions.get_by_id(_rxn).metabolites[met]
                print(f'{_rxn}: stoich={stoich:+.2f}, flux={flux:.4f}')
            non_zero_fluxes = {_k: _f for _k, _f in _sol.fluxes.to_dict().items() if _f != 0}
            break
        break
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ##### Claude updated visualization
    """)
    return


@app.cell
def _(
    DataFrame,
    LinearRegression,
    LogFormatterMathtext,
    LogLocator,
    NullFormatter,
    Series,
    adjust_text,
    json,
    np,
    pi,
    plt,
    r2_score,
    read_csv,
):
    import matplotlib
    matplotlib.use('Agg')
    from sklearn.preprocessing import PolynomialFeatures
    from sklearn.pipeline import make_pipeline, Pipeline

    def fluxes_to_emissions_1(row):
        area_cm2 = pi * (5 / 2) ** 2
        area_m2 = area_cm2 * 0.01 ** 2
        volume_cm3 = 294.52
        ugDNA_gSoil = 2
        DNA_per_biomass = 0.1
        gSoil_cm3 = 1.4
        gBiomass_gSoil = gSoil_cm3 * volume_cm3 * ugDNA_gSoil / DNA_per_biomass / 1000 / 1000
        return _row / area_m2 * 1000.0 * 24 * gBiomass_gSoil
    methane_dic_1 = json.load(open('averaged_normalized_methane_unrestored.json'))
    averaged_unrestored_normalized_abundances_1 = read_csv('averaged_normalized_abundances_unrestored.csv', index_col=0)
    averaged_unrestored_normalized_methane_1 = Series(methane_dic_1)
    betaine_2 = {'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs', 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs', 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs', 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs', 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs', 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs', 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs', 'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs', 'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs', 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs', 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs', 'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs', 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs'}
    betaine_striped_1 = [x.replace('.contigs', '') for x in betaine_2]
    mepe_1 = {'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.20.contigs', 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.32.contigs', 'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.51.contigs', 'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_metabat.27.contigs', 'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_metabat.28.contigs', 'Salt_Pond_MetaG_R1_B_H2O_MG_DASTool_bins_metabat.7.contigs', 'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_metabat.31.contigs', 'Salt_Pond_MetaG_R1_C_H2O_MG_DASTool_bins_concoct_out.79.contigs', 'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.35.contigs', 'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.49.contigs', 'Salt_Pond_MetaG_R2_A_D2_MG_DASTool_bins_metabat.14.contigs', 'Salt_Pond_MetaG_R2_A_D2_MG_DASTool_bins_metabat.16.contigs', 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_concoct_out.143.contigs', 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_concoct_out.32.contigs', 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.10.contigs', 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.16.contigs', 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.19.contigs', 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.38.contigs', 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.55.contigs', 'Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_concoct_out.88.contigs', 'Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_concoct_out.95.contigs', 'Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_metabat.17.contigs', 'Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.18.contigs', 'Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.45.contigs', 'Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.47.contigs', 'Salt_Pond_MetaG_R2_restored_C_black_MG_DASTool_bins_concoct_out.17.contigs', 'Salt_Pond_MetaG_R2_restored_C_black_MG_DASTool_bins_metabat.8.contigs'}
    mepe_striped_1 = [x.replace('.contigs', '') for x in mepe_1]
    _log = False
    _linear = True
    for _fluxes_file in ['saltern_fluxes.json', 'saltern_fluxes_betaine.json']:
        fluxes_4 = json.load(open(_fluxes_file))
        if '0.07' in fluxes_4:
            flux_groups = {'methane': fluxes_4['0.07']['methane']}
        else:
            flux_groups = fluxes_4
        for flux_type, flux_data in flux_groups.items():
            summed = DataFrame(flux_data).T.sum().abs()
            _fluxes_in_emissions = fluxes_to_emissions_1(summed).to_dict()
            _combined = {_ID: {'fluxes': em, 'methane': methane_dic_1[_ID]} for _ID, em in _fluxes_in_emissions.items() if _ID in methane_dic_1}
            _combined_df = DataFrame(_combined).T
            _methane_ary = np.array(list(_combined_df['methane'].values))
            _emission_fluxes = np.array(list(_combined_df['fluxes'].values))
            _fig, _ax = plt.subplots(figsize=(8, 6), dpi=300)
            _xs = _methane_ary.reshape(-1, 1)
            _ys = _emission_fluxes
            if _log:
                _xs = np.log(_xs)
                _ys = np.log(_ys)
                title = 'Methane emissions $\\left[LOG\\left(\\frac{umol~~CH_4}{m^2*day*\\left(\\frac{ug DNA}{g soil}\\right)}\\right)\\right]$'
            else:
                _ax.set_xscale('log')
                _ax.set_yscale('log')
                _ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0,)))
                _ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0,)))
                _ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
                _ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
                fmt = LogFormatterMathtext(base=10, labelOnlyBase=True)
                _ax.xaxis.set_major_formatter(fmt)
                _ax.yaxis.set_major_formatter(fmt)
                _ax.xaxis.set_minor_formatter(NullFormatter())
                _ax.yaxis.set_minor_formatter(NullFormatter())
                plt.setp(_ax.get_xticklabels(), rotation=90)
                title = f'{flux_type.capitalize()} flux correlation'
            _mdl = LinearRegression()
            _mdl.fit(_xs, _ys)
            _ax.scatter(_xs.flatten(), _ys)
            _ax.grid(True, which='both', axis='both', linestyle='--', color='gray', alpha=0.3)
            _label_size = 13
            texts = []
            for _xi, _yi, _label in zip(_xs.flatten(), _ys, list(_combined_df.index)):
                texts.append(_ax.text(_xi, _yi, _label, fontsize=_label_size, ha='center', va='bottom'))
            adjust_text(texts, ax=_ax, arrowprops=dict(arrowstyle='->', color='gray'))
            _y_pred = _mdl.predict(_xs)
            trendline_color = 'green'
            _sorted_idx = np.argsort(_methane_ary)
            _ax.plot(_xs.flatten()[_sorted_idx], _y_pred[_sorted_idx], color=trendline_color, linestyle='--', label='Trendline')
            _eq = f'y = {_mdl.coef_[0]:.2f}x + {_mdl.intercept_:.2f}'
            _r2_text = f'$R^2$ = {r2_score(_ys, _y_pred):.5f}'
            _ax.text(0.3, 0.8, _eq + '\n' + _r2_text, transform=_ax.transAxes, fontsize=_label_size, color=trendline_color, verticalalignment='top')
            _axes_size = 20
            _ax.set_xlabel('Experimental measurement', fontsize=_axes_size)
            _ax.set_ylabel('Metabolic modeling (|flux|)', fontsize=_axes_size)
            _ax.set_title(title, fontsize=_axes_size + 4)
            plt.tight_layout()
            fname = f"experiment_modeling_{_fluxes_file.replace('.json', '')}_{flux_type}.png"
            _fig.savefig(fname)
            plt.close(_fig)
            print(f'Saved {fname}: R2={r2_score(_ys, _y_pred):.5f}')
    for _samples in ['MePn', 'Betaine']:
        _fig, _ax = plt.subplots(figsize=(8, 6), dpi=300)
        if _samples == 'MePn':
            _sample_totals = averaged_unrestored_normalized_abundances_1.loc[list(mepe_striped_1)].sum()
            title = 'MePn consumer correlation'
            export_name = 'MePn_consumers_new.png'
        else:
            _sample_totals = averaged_unrestored_normalized_abundances_1.loc[list(betaine_striped_1)].sum()
            title = 'Betaine consumer correlation'
            export_name = 'Betaine_consumers_new.png'
        _all_xs = np.array(_sample_totals.to_list())
        _xs = _all_xs.reshape(-1, 1)
        _all_ys = averaged_unrestored_normalized_methane_1.to_list()
        _ys = np.array(_all_ys)
        if _log:
            _xs = np.log(_xs)
            _ys = np.log(_ys)
            _ax.set_xlabel('Abundances $\\left[LOG\\left(\\frac{rel.ab}{\\frac{ug DNA}{g soil}}\\right)\\right]$', fontsize=20)
            _ax.set_ylabel('Methane emissions $\\left[LOG\\left(\\frac{umol~~CH_4}{m^2*day*\\left(\\frac{ug DNA}{g soil}\\right)}\\right)\\right]$', fontsize=20)
        else:
            _ax.set_xscale('log')
            _ax.set_yscale('log')
            _ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0,)))
            _ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0,)))
            _ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
            _ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
            fmt = LogFormatterMathtext(base=10, labelOnlyBase=True)
            _ax.xaxis.set_major_formatter(fmt)
            _ax.yaxis.set_major_formatter(fmt)
            _ax.xaxis.set_minor_formatter(NullFormatter())
            _ax.yaxis.set_minor_formatter(NullFormatter())
            plt.setp(_ax.get_xticklabels(), rotation=90)
            _ax.set_xlabel('Abundances $\\left[\\left(\\frac{rel.ab}{\\frac{ug DNA}{g soil}}\\right)\\right]$', fontsize=20)
            _ax.set_ylabel('Methane emissions $\\left[\\left(\\frac{umol~~CH_4}{m^2*day*\\left(\\frac{ug DNA}{g soil}\\right)}\\right)\\right]$', fontsize=20)
        if _linear:
            _model = LinearRegression()
            _model.fit(_xs, _ys)
            _y_pred = _model.predict(_xs)
            _coef = _model.coef_[0]
            _intercept = _model.intercept_
        else:
            _poly = PolynomialFeatures(degree=2, include_bias=False)
            _lin = LinearRegression()
            _model = Pipeline([('poly', _poly), ('lin', _lin)])
            _model.fit(_xs, _ys)
            _y_pred = _model.predict(_xs)
            _term_names = _model.named_steps['poly'].get_feature_names_out()
            _coef = _model.named_steps['lin'].coef_
            _intercept = _model.named_steps['lin'].intercept_
        _ax.scatter(_all_xs, _ys)
        _label_size = 13
        texts = []
        for _xi, _yi, _label in zip(_all_xs, _ys, list(_sample_totals.index)):
            texts.append(_ax.text(_xi, _yi, _label, fontsize=_label_size, ha='center', va='bottom'))
        adjust_text(texts, ax=_ax, arrowprops=dict(arrowstyle='->', color='gray'))
        _sorted_idx = np.argsort(_all_xs)
        sorted_xs = _all_xs[_sorted_idx]
        sorted_ys = _y_pred[_sorted_idx]
        trendline_color = 'green'
        _ax.plot(sorted_xs, sorted_ys, color=trendline_color, linestyle='--', label='Trendline')
        _ax.grid(True, which='both', axis='both', linestyle='--', color='gray', alpha=0.3)
        if _linear:
            _eq = f'y = {_coef:.2f}x + {_intercept:.2f}'
        else:
            _precision = 4
            _parts = [f'{_intercept:.{_precision}g}']
            for _c, _t in zip(_coef, _term_names):
                _sign = ' + ' if _c >= 0 else ' - '
                _parts.append(f'{_sign}{abs(_c):.{_precision}g}*{_t}')
            _eq = 'y = ' + ''.join(_parts)
        _r2_text = f'$R^2$ = {r2_score(_all_ys, _y_pred):.3f}'
        _ax.text(0.2, 0.7, _eq + '\n' + _r2_text, transform=_ax.transAxes, fontsize=_label_size, color=trendline_color, verticalalignment='top')
        _axes_size = 20
        _ax.set_title(title, fontsize=_axes_size + 4)
        plt.tight_layout()
        _fig.savefig(export_name)
        plt.close(_fig)
        print(f'Saved {export_name}')
    print('\nDone!')
    return Pipeline, PolynomialFeatures


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # heatmap co-occurrence assessment
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## IterativeIDs
    """)
    return


@app.cell
def _(display, load):
    from json import dump
    phylogeny_1 = load(open('../data/Saltern_phylogeny.json'))
    _iterativeIDs, iterativeID_levels, iterativeID_phylums = ({}, {}, {})
    for MAG, _taxa in phylogeny_1.items():
        best_level = 'Kingdom'
        for level, taxon in _taxa.items():
            if not isinstance(taxon, str) or level == 'Species' or taxon == '':
                break
            best_level = level
        taxon = _taxa[best_level]
        taxon = taxon + '.1'
        while taxon in _iterativeIDs:
            taxon, _count = taxon.split('.')
            taxon = f'{taxon}.{int(_count) + 1}'
        _iterativeIDs[taxon] = MAG.replace('.contigs__', '')
        iterativeID_levels[taxon] = best_level
        iterativeID_phylums[taxon] = _taxa['Phylum']
    _iterativeIDs = {_v: _k for _k, _v in _iterativeIDs.items()}
    display(list(_iterativeIDs.items())[:5])
    display(list(iterativeID_levels.items())[:5])
    display(list(iterativeID_phylums.items())[:5])
    dump(_iterativeIDs, open('iterativeIDs.json', 'w'))
    dump(iterativeID_levels, open('iterativeID_levels.json', 'w'))
    dump(iterativeID_phylums, open('iterativeID_phylums.json', 'w'))
    print('None taxonomy entries:', [_k for _k, _v in iterativeID_phylums.items() if _v is None])
    betaine_3 = {'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs', 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs', 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs', 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs', 'Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs', 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs', 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs', 'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs', 'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs', 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs', 'Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs', 'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs', 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs'}
    betaine_iterativeIDs = [_iterativeIDs[x.replace('.contigs', '')] for x in betaine_3]
    display(betaine_iterativeIDs)
    mepe_2 = {'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.20.contigs', 'Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.32.contigs', 'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.51.contigs', 'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_metabat.27.contigs', 'Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_metabat.28.contigs', 'Salt_Pond_MetaG_R1_B_H2O_MG_DASTool_bins_metabat.7.contigs', 'Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_metabat.31.contigs', 'Salt_Pond_MetaG_R1_C_H2O_MG_DASTool_bins_concoct_out.79.contigs', 'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.35.contigs', 'Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.49.contigs', 'Salt_Pond_MetaG_R2_A_D2_MG_DASTool_bins_metabat.14.contigs', 'Salt_Pond_MetaG_R2_A_D2_MG_DASTool_bins_metabat.16.contigs', 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_concoct_out.143.contigs', 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_concoct_out.32.contigs', 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.10.contigs', 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.16.contigs', 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.19.contigs', 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.38.contigs', 'Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.55.contigs', 'Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_concoct_out.88.contigs', 'Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_concoct_out.95.contigs', 'Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_metabat.17.contigs', 'Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.18.contigs', 'Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.45.contigs', 'Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.47.contigs', 'Salt_Pond_MetaG_R2_restored_C_black_MG_DASTool_bins_concoct_out.17.contigs', 'Salt_Pond_MetaG_R2_restored_C_black_MG_DASTool_bins_metabat.8.contigs'}
    mepe_iterativeIDs = [_iterativeIDs.get(x.replace('.contigs', ''), x.replace('.contigs', '')) for x in mepe_2]
    display(mepe_iterativeIDs)
    return (
        betaine_iterativeIDs,
        dump,
        iterativeID_phylums,
        mepe_iterativeIDs,
        phylogeny_1,
    )


@app.cell
def _(
    DataFrame,
    Series,
    display,
    dump,
    iterativeID_phylums,
    load,
    np,
    phylogeny_1,
    read_csv,
    spearmanr,
):
    from matplotlib import colors, pyplot, patches
    from matplotlib.patches import Patch
    from numpy import ones, triu, nan, ones_like
    import seaborn as sns
    from statsmodels.stats.multitest import multipletests
    from collections import defaultdict
    from itertools import combinations

    def corr_pvalues(df):
        _n = df.shape[1]
        pvals = DataFrame(ones((_n, _n)), index=df.columns, columns=df.columns)
        corrs = DataFrame(ones((_n, _n)), index=df.columns, columns=df.columns)
        for i in range(_n):
            for _j in range(i + 1, _n):
                _c, p = spearmanr(df.iloc[:, i], df.iloc[:, _j])
                pvals.iloc[i, _j] = p
                corrs.iloc[i, _j] = _c
                pvals.iloc[_j, i] = p
                corrs.iloc[_j, i] = _c
        return (pvals, corrs)

    def merge_depths_by_DNA_1(abundances):
        return DataFrame({prefix: (abundances[prefix + '_D1'] * _sample_dna[prefix + '_D1'] + abundances[prefix + '_D2'] * _sample_dna[prefix + '_D2']) / (_sample_dna[prefix + '_D1'] + _sample_dna[prefix + '_D2']) for prefix in _unique_prefixes}, index=abundances.index)

    def merge_triplicates(abundances):
        return DataFrame({site: (abundances[site + '_A'] + abundances[site + '_B'] + abundances[site + '_C']) / 3 for site in _sites}, index=abundances.index)
    dna_df_1 = read_csv('../data/SaltPondsDNA.csv').set_index('Sample')
    _sample_dna = dna_df_1['DNA conc (µg/g soil)'].to_dict()
    _iterativeIDs = load(open('iterativeIDs.json', 'r'))
    iterativeID_levels_1 = load(open('iterativeID_levels.json', 'r'))
    _iterativeID_taxonomy = {_iterativeIDs.get(_ID, _ID): phylogeny_1.get(_ID, 'Unknown') for _ID in _iterativeIDs}
    _unique_prefixes = {'_'.join(i.split('_')[:-1]) for i in dna_df_1.index}
    _sites = {i.split('_')[0] for i in _unique_prefixes}
    _volume_weighted_sample_dna = Series({prefix: (5 * _sample_dna[prefix + '_D1'] + 10 * _sample_dna[prefix + '_D2']) / 15 for prefix in _unique_prefixes})
    abundances_1 = merge_depths_by_DNA_1(read_csv('../data/Cliff_310MAG_relabund.tsv', sep='\t'))
    abundances_1 = merge_triplicates(abundances_1).T
    _zero_level = 1e-05
    abundances_1 = abundances_1.loc[:, (abundances_1.fillna(0) > _zero_level).sum() >= 3]
    abundances_1.columns = [_iterativeIDs.get(_ID, _ID) for _ID in abundances_1.columns]
    abundances_1.rename(columns={'': 'unbinned'}, inplace=True)
    defined = True
    if defined:
        abundances_1.drop(columns=[_c for _c in abundances_1.columns if _c not in _iterativeIDs.values()], inplace=True)
    display(abundances_1)
    df_1 = abundances_1
    reduced = True
    if reduced:
        reduced_abundances = abundances_1.copy()
        reduced_abundances.columns = ['.'.join(_ID.split('.')[:-1]) for _ID in reduced_abundances.columns]
        reduced_abundances = reduced_abundances.groupby(reduced_abundances.columns, axis=1).sum()
        display(reduced_abundances)
        df_1 = reduced_abundances.rename(columns={'': 'unbinned'})
        iterativeID_phylums_1 = {'.'.join(_ID.split('.')[:-1]): _v for _ID, _v in iterativeID_phylums.items() if '.' in _ID}
    iterativeID_phylums_1 = {_ID: _v for _ID, _v in iterativeID_phylums_1.items() if _ID in df_1.columns}
    pval_matrix, corr_matrix = corr_pvalues(df_1)
    display(corr_matrix)
    _presence = (df_1 > _zero_level).astype(int)
    _cooccurrence = defaultdict(int)
    for _sample in _presence.itertuples(index=False):
        _present = [_col for _col, _val in zip(_presence.columns, _sample) if _val]
        for _pair in combinations(sorted(_present), 2):
            _cooccurrence[_pair] = _cooccurrence[_pair] + 1
    pair_data = []
    for (_a, _b), _count in _cooccurrence.items():
        _rho = corr_matrix.loc[_a, _b]
        if np.isnan(_rho):
            continue
        p_value = pval_matrix.loc[_a, _b]
        pair_data.append((_a, _b, _rho, p_value, _count))
    pvals = np.array([_t[3] for _t in pair_data])
    np.savetxt('pvals.txt', pvals)
    reject, qvals, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
    print(len(pvals), sum(reject))
    level_1 = 'Phylum'
    _orgs = set()
    for (_a, _b, _, _, _), take in zip(pair_data, reject):
        if take:
            _orgs.add(_a)
            _orgs.add(_b)
    dump(list(_orgs), open('significantly_connected_organisms.json', 'w'))
    archaea_phyla = sorted({_phylum for _phylum in iterativeID_phylums_1.values() if _phylum and ('archaeo' in _phylum.lower() or any((_a in _phylum.lower() for _a in ['candidatus thermoplasmatota', 'halobacterota', 'methanobacteriota', 'micrarchaeota'])))})
    bacteria_phyla = sorted({_phylum for _phylum in iterativeID_phylums_1.values() if _phylum and _phylum not in archaea_phyla})
    _n_archaea = len(archaea_phyla)
    _archaea_colors = [pyplot.cm.turbo(i / max(_n_archaea, 1) * 0.15) for i in range(_n_archaea)]
    _n_bacteria = len(bacteria_phyla)
    _bacteria_colors = [pyplot.cm.turbo(0.2 + i / max(_n_bacteria, 1) * 0.8) for i in range(_n_bacteria)]
    taxa_color_map = {}
    for _phylum, _color in zip(archaea_phyla, _archaea_colors):
        taxa_color_map[_phylum] = _color
    for _phylum, _color in zip(bacteria_phyla, _bacteria_colors):
        taxa_color_map[_phylum] = _color
    _iterativeID_color_map = {_ID: taxa_color_map[_phylum] for _ID, _phylum in iterativeID_phylums_1.items() if _phylum}
    DEFAULT_COLOR = 'lightgray'
    _bar_label = 'Phylum'
    row_colors = Series({idx: _iterativeID_color_map.get(idx, DEFAULT_COLOR) for idx in corr_matrix.index}, name=_bar_label)
    col_colors = Series({_col: _iterativeID_color_map.get(_col, DEFAULT_COLOR) for _col in corr_matrix.columns}, name=_bar_label)
    display(row_colors[:5])
    return (
        DEFAULT_COLOR,
        abundances_1,
        archaea_phyla,
        bacteria_phyla,
        col_colors,
        combinations,
        corr_matrix,
        defaultdict,
        defined,
        df_1,
        iterativeID_phylums_1,
        level_1,
        nan,
        ones_like,
        pair_data,
        patches,
        pval_matrix,
        pyplot,
        qvals,
        reduced,
        reject,
        row_colors,
        sns,
        taxa_color_map,
        triu,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## all v all heatmap
    """)
    return


@app.cell
def _(
    DEFAULT_COLOR,
    Series,
    archaea_phyla,
    bacteria_phyla,
    betaine_iterativeIDs,
    col_colors,
    corr_matrix,
    defined,
    df_1,
    level_1,
    load,
    mepe_iterativeIDs,
    nan,
    ones_like,
    patches,
    phylogeny_1,
    pval_matrix,
    reduced,
    row_colors,
    sns,
    taxa_color_map,
    triu,
):
    taxonomies = {}
    for _col in df_1.columns:
        taxonomies[_col] = '|'.join([_v for _k, _v in phylogeny_1.get(_col, {'Unknown': None}).items() if _k != 'Species' and _v is not None])
    taxonomy_series = Series({idx: taxonomies.get(idx, f'Unknown|{idx}') for idx in df_1.columns})
    print(f"Detected depth: {max((len(_t.split('|')) for _t in taxonomy_series))}")
    print(f'df rows: {len(df_1.columns)}')
    print(f'taxonomy_series length: {len(taxonomy_series)}')
    print(f'Sample entries:\n{taxonomy_series.head()}')
    clusterMap = sns.clustermap(corr_matrix, row_colors=row_colors, col_colors=col_colors, cmap='coolwarm_r', center=0, figsize=(90, 100), dendrogram_ratio=(0.1, 0.2))
    clusterMap.figure.subplots_adjust(bottom=0.15, top=0.95)
    clusterMap.ax_row_dendrogram.set_visible(False)
    clusterMap.ax_col_dendrogram.set_visible(False)
    secondary_labels = [tick.get_text() for tick in clusterMap.ax_heatmap.get_yticklabels()]
    clusterMap.ax_heatmap.yaxis.set_ticks_position('left')
    clusterMap.ax_heatmap.yaxis.set_label_position('left')
    clusterMap.ax_heatmap.set_yticklabels(secondary_labels, rotation=0)
    hm_pos = clusterMap.ax_heatmap.get_position()
    fig_w = clusterMap.figure.get_figwidth()
    fig_h = clusterMap.figure.get_figheight()
    strip_w = 0.015
    strip_h = 0.015
    clusterMap.ax_row_colors.set_position([hm_pos.x0 - strip_w, hm_pos.y0, strip_w, hm_pos.height])
    clusterMap.ax_col_colors.set_position([hm_pos.x0, hm_pos.y0 - strip_h, hm_pos.width, strip_h])
    y_pad = strip_w * fig_w * 72 + 15
    clusterMap.ax_heatmap.tick_params(axis='y', pad=y_pad)
    x_pad = strip_h * fig_h * 72 + 15
    clusterMap.ax_heatmap.tick_params(axis='x', pad=x_pad)
    _cbar = clusterMap.ax_cbar
    ticks = _cbar.get_yticks()
    ticks[-1] = corr_matrix.max().max()
    ticks[0] = round(corr_matrix.min().min(), 2)
    _cbar.set_yticklabels(ticks, fontsize=45)
    _cbar.set_yticks(ticks)
    _cbar.set_xlabel('Spearman Correlation', fontsize=45, labelpad=30)
    heatmap_pos = clusterMap.ax_heatmap.get_position()
    clusterMap.ax_cbar.set_position([heatmap_pos.x1 - 0.96, heatmap_pos.y0 + 0.68, 0.06, heatmap_pos.height / 5])
    labelsize = 40
    clusterMap.ax_heatmap.set_yticklabels(clusterMap.ax_heatmap.get_yticklabels(), fontsize=labelsize, rotation=0)
    clusterMap.ax_heatmap.set_xticklabels(clusterMap.ax_heatmap.get_xticklabels(), fontsize=labelsize, rotation=60, ha='right', rotation_mode='anchor')
    color_scheme = {'betaine': {'scatter': '#1f77b4', 'trendline': '#0a3d6b', 'label': '#1f77b4'}, 'mepe': {'scatter': '#e8671b', 'trendline': '#7f2704', 'label': '#e8671b'}}
    iterativeID_levels_2 = load(open('iterativeID_levels.json', 'r'))
    print(betaine_iterativeIDs)
    betaine_mems = betaine_iterativeIDs
    mepe_mems = mepe_iterativeIDs
    if reduced:
        betaine_mems = ['.'.join(x.split('.')[:-1]) for x in betaine_iterativeIDs]
        mepe_mems = ['.'.join(x.split('.')[:-1]) for x in mepe_iterativeIDs]

    def update_labels(labels):
        for _label in labels:
            text = _label.get_text()
            in_mepe = any((x in text for x in mepe_mems))
            in_betaine = any((x in text for x in betaine_mems))
            if in_mepe or in_betaine:
                _orgs.add(text)
                _label.set_fontsize(labelsize * 1.2)
                _label.set_fontweight('bold')
                if in_betaine:
                    _label.set_color(color_scheme['betaine']['label'])
                else:
                    _label.set_color(color_scheme['mepe']['label'])
            if iterativeID_levels_2.get(text) == 'Genus':
                _label.set_fontstyle('italic')
    iterativeID_levels_2 = load(open('iterativeID_levels.json', 'r'))
    ylabels = clusterMap.ax_heatmap.get_yticklabels()
    _orgs = set()
    update_labels(ylabels)
    xlabels = clusterMap.ax_heatmap.get_xticklabels()
    update_labels(xlabels)
    dendrogram_row = clusterMap.dendrogram_row.reordered_ind
    dendrogram_col = clusterMap.dendrogram_col.reordered_ind
    one_triangle = True
    if one_triangle:
        df_reordered = corr_matrix.iloc[dendrogram_row, dendrogram_col]
        _mask = triu(ones_like(df_reordered, dtype=bool), k=1)
        mesh = clusterMap.ax_heatmap.collections[0]
        arr = mesh.get_array().reshape(df_reordered.shape)
        arr[_mask] = nan
        mesh.set_array(arr.ravel())
        clusterMap.ax_cbar.set_position([heatmap_pos.x1 - 0.1, heatmap_pos.y0 + 0.3, 0.06, heatmap_pos.height / 5])
    for _org in _orgs:
        orgIx = corr_matrix.index.get_loc(_org)
        if orgIx not in dendrogram_row:
            continue
        row_pos = dendrogram_row.index(orgIx)
        rect = patches.Rectangle((0, row_pos), len(corr_matrix.columns) if not one_triangle else row_pos + 1, 1, linewidth=6, edgecolor='black', facecolor='none', clip_on=False)
        clusterMap.ax_heatmap.add_patch(rect)
        colIx = corr_matrix.columns.get_loc(_org)
        col_pos = dendrogram_col.index(colIx)
        rect = patches.Rectangle((col_pos, 0) if not one_triangle else (col_pos, len(corr_matrix.index)), 1, len(corr_matrix.index) if not one_triangle else -(len(corr_matrix.index) - col_pos), linewidth=6, edgecolor='black', facecolor='none', clip_on=False)
        clusterMap.ax_heatmap.add_patch(rect)
    pvals_reordered = pval_matrix.iloc[dendrogram_row, dendrogram_col]
    _ax = clusterMap.ax_heatmap
    for i_1 in range(pvals_reordered.shape[0]):
        for _j in range(pvals_reordered.shape[1]):
            if pvals_reordered.iloc[i_1, _j] < 0.05 and (not one_triangle or _mask[i_1, _j] == 0):
                _ax.add_patch(patches.Rectangle((_j, i_1), 1, 1, fill=False, edgecolor='lightgreen', linewidth=8))
    clusterMap.ax_heatmap.set_xlabel('Member ASVs', fontsize=50)
    clusterMap.ax_heatmap.set_ylabel('Member ASVs', fontsize=50)

    def _header_patch(title):
        """Section header: invisible swatch with a bold label."""
        return patches.Patch(color='none', label=f'$\\bf{{{title}}}$')
    fig_width = clusterMap.figure.get_figwidth()
    _handle = []
    if len(archaea_phyla) > 0:
        _archaea_patches = [patches.Patch(color=taxa_color_map.get(p, DEFAULT_COLOR), label=p) for p in archaea_phyla]
        _handle = [_header_patch('Archaea')] + _archaea_patches
    if len(bacteria_phyla) > 0:
        _bacteria_patches = [patches.Patch(color=taxa_color_map.get(p, DEFAULT_COLOR), label=p) for p in bacteria_phyla]
        _handle = _handle + ([_header_patch('Bacteria')] + _bacteria_patches)
    clusterMap.ax_heatmap.legend(handles=_handle, title=f'Taxonomic {level_1}', title_fontsize=8 * (fig_width / 10), loc='lower left', bbox_to_anchor=(0.65, 0.5), fontsize=7 * (fig_width / 10), frameon=True)
    clusterMap.figure.savefig(f"abundance_correlatons{('_one_triangle' if one_triangle else '')}{('_reduced' if reduced else '')}{('_defined' if defined else '')}.png", dpi=300, bbox_inches='tight')
    return betaine_mems, color_scheme, iterativeID_levels_2, mepe_mems


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## co-occurance network
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    #### proposed data processing from Claude
    """)
    return


@app.cell
def _(V, abundances_1, condition, np, plt, spearmanr):
    from skbio.stats.composition import clr
    _abund_pseudo = abundances_1 + 0.5
    _X = clr(_abund_pseudo)
    _R = _X.copy()
    for _c in np.unique(condition):
        _mask = condition == _c
        _R[_mask] = _R[_mask] - _R[_mask].mean(axis=0, keepdims=True)
    _C_marginal = np.corrcoef(_X.T)
    _C_within = np.corrcoef(_R.T)
    _X_means = np.vstack([_X[condition == _c].mean(0) for _c in np.unique(condition)])
    _C_between, _ = spearmanr(_X_means)
    _mask = np.abs(_C_within) > 0.7
    np.fill_diagonal(_mask, False)
    direct_pairs = np.argwhere(np.triu(_mask))
    plt.scatter(_C_marginal[np.triu_indices(V, 1)], _C_within[np.triu_indices(V, 1)], s=2, alpha=0.3)
    plt.plot([-1, 1], [-1, 1], 'k--')
    plt.xlabel('marginal r')
    plt.ylabel('within r')
    from sklearn.covariance import GraphicalLassoCV
    _gl = GraphicalLassoCV().fit(_R)
    _precision = _gl.precision_
    return GraphicalLassoCV, clr


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### co-occurrence network
    """)
    return


@app.cell
def _(
    DEFAULT_COLOR,
    archaea_phyla,
    bacteria_phyla,
    betaine_mems,
    color_scheme,
    defined,
    df_1,
    dump,
    iterativeID_levels_2,
    iterativeID_phylums_1,
    level_1,
    mepe_mems,
    np,
    pair_data,
    patches,
    plt,
    qvals,
    reduced,
    reject,
    taxa_color_map,
):
    import networkx as nx
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    import math
    G = nx.Graph()
    for (_a, _b, _rho, p, _count), q, keep in zip(pair_data, qvals, reject):
        if keep:
            G.add_edge(_a, _b, weight=abs(_rho), rho=_rho, pvalue=p, qvalue=q, cooccurrence=_count)
    print(f'Tests run: {len(pair_data)}')
    print(f'FDR-significant pairs (q < 0.05): {int(reject.sum())}')
    print(f'Edges after FDR correction: {G.number_of_edges()}')
    print(f'Nodes in graph: {G.number_of_nodes()}')
    pos_1 = nx.spring_layout(G, seed=42, iterations=200, k=3)
    edges = G.edges(data=True)
    rho_values = [_d['rho'] for _, _, _d in edges]
    edge_widths = [3 * _d['weight'] for _, _, _d in edges]
    print(sorted(G.degree(), key=lambda x: x[1], reverse=True)[:10])
    norm = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    cmap = cm.RdBu
    edge_colors = [cmap(norm(r)) for r in rho_values]
    width = 40
    height = 30
    _fig, _ax = plt.subplots(figsize=(width, height))
    scale = 5000 * (width / 10) * 2
    compressor = np.sqrt
    relative_abundance = df_1.div(df_1.sum(axis=1), axis=0)
    mean_rel_abund = relative_abundance.mean(axis=0)
    node_sizes = [scale * compressor(mean_rel_abund.get(_n, 0)) for _n in G.nodes()]
    node_colors = [taxa_color_map.get(iterativeID_phylums_1.get(_n, 'Unknown'), DEFAULT_COLOR) for _n in G.nodes()]
    print('unique node colors:', len(set(node_colors)))
    nx.draw_networkx_nodes(G, pos_1, node_size=node_sizes, node_color=node_colors, alpha=0.9, ax=_ax)
    nx.draw_networkx_edges(G, pos_1, width=edge_widths, edge_color=edge_colors, alpha=0.85, ax=_ax)
    significantly_connected_organisms = list(pos_1.keys())
    dump(significantly_connected_organisms, open('significantly_connected_organisms.json', 'w'))
    for _n, (x, y) in pos_1.items():
        text = str(_n)
        abund = mean_rel_abund.get(_n, 0)
        fs = 6 + 14 * compressor(abund) / compressor(mean_rel_abund.max()) * (width / 10)
        fs = max(5, min(fs, 48))
        weight = 'bold'
        style = 'normal'
        _color = 'black'
        in_mepe = any((s in text for s in mepe_mems))
        in_betaine = any((s in text for s in betaine_mems))
        if in_mepe or in_betaine:
            fs = fs * 1.2
            _color = color_scheme['betaine']['label'] if in_betaine else color_scheme['mepe']['label']
        if iterativeID_levels_2.get(text) == 'Genus':
            style = 'italic'
        _ax.text(x, y, text, fontsize=fs, color=_color, fontweight=weight, fontstyle=style, ha='center', va='center')
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    _cbar = plt.colorbar(sm, ax=_ax, shrink=0.6, pad=0.02)
    _cbar.set_label('Spearman ρ', fontsize=10 * (width / 10))
    _cbar.ax.tick_params(labelsize=8 * (width / 10), length=8, width=2)
    legend_entries = [(0.002, '0.2%'), (0.02, '2%'), (0.2, '20%')]
    sizes_pt2 = [scale * compressor(_a) for _a, _ in legend_entries]
    diameters_pt = [2 * np.sqrt(s / np.pi) for s in sizes_pt2]
    breather_pt = 8
    offsets_pt = [0.0]
    for i_2 in range(1, len(legend_entries)):
        offsets_pt.append(offsets_pt[-1] + (diameters_pt[i_2 - 1] + diameters_pt[i_2]) / 2 + breather_pt)
    fig_h_pts = _fig.get_figheight() * 72
    top_y = 0.9
    circle_x = 0.85
    label_x = 0.88
    _ax.text(circle_x, top_y + 0.025, 'Mean rel. abundance', transform=_fig.transFigure, va='bottom', fontweight='bold', fontsize=6 * (width / 10), clip_on=False)
    for (abund, _label), s, off in zip(legend_entries, sizes_pt2, offsets_pt):
        y = top_y - off / fig_h_pts
        _ax.scatter([circle_x], [y], s=s, color='slategray', alpha=0.9, transform=_fig.transFigure, clip_on=False)
        _ax.text(label_x, y, _label, transform=_fig.transFigure, va='center', fontsize=5 * (width / 10), clip_on=False)

    def _header_patch(title):
        """Section header: invisible swatch with a bold label."""
        return patches.Patch(color='none', label=f'$\\bf{{{title}}}$')
    _handle = []
    if len(archaea_phyla) > 0:
        _archaea_patches = [patches.Patch(color=taxa_color_map.get(p, DEFAULT_COLOR), label=p) for p in archaea_phyla]
        _handle = [_header_patch('Archaea')] + _archaea_patches
    if len(bacteria_phyla) > 0:
        _bacteria_patches = [patches.Patch(color=taxa_color_map.get(p, DEFAULT_COLOR), label=p) for p in bacteria_phyla]
        _handle = _handle + ([_header_patch('Bacteria')] + _bacteria_patches)
    _ax.legend(handles=_handle, title=f'Taxonomic {level_1}', title_fontsize=8 * (width / 10), loc='lower left', bbox_to_anchor=(-0.24, 0.1), fontsize=7 * (width / 10), frameon=True)
    _ax.axis('off')
    plt.tight_layout()
    plt.savefig(f"cooccurrence_network_p_value_FDR{('_reduced' if reduced else '')}{('_defined' if defined else '')}.png", dpi=300, bbox_inches='tight')
    plt.show()
    return G, compressor, mean_rel_abund, pos_1


@app.cell
def _(
    DataFrame,
    GraphicalLassoCV,
    Series,
    clr,
    combinations,
    defaultdict,
    display,
    dump,
    iterativeID_phylums_1,
    load,
    np,
    phylogeny_1,
    pyplot,
    read_csv,
    spearmanr,
):
    from scipy.stats import pearsonr
    from itertools import permutations, product

    def merge_depths_by_DNA_2(abundances):
        return DataFrame({prefix: (abundances[prefix + '_D1'] * _sample_dna[prefix + '_D1'] + abundances[prefix + '_D2'] * _sample_dna[prefix + '_D2']) / (_sample_dna[prefix + '_D1'] + _sample_dna[prefix + '_D2']) for prefix in _unique_prefixes}, index=abundances.index)

    def expand_triplicates_long(abundances_depth_merged):
        """Return (12 x V) DataFrame with one row per replicate, plus condition labels."""
        rows = {}
        condition = []
        for site in _sites:
            for rep in ['A', 'B', 'C']:
                sample_id = f'{site}_{rep}'
                rows[sample_id] = abundances_depth_merged[sample_id]
                condition.append(site)
        return (DataFrame(rows).T, np.array(condition))
    dna_df_2 = read_csv('../data/SaltPondsDNA.csv').set_index('Sample')
    _sample_dna = dna_df_2['DNA conc (µg/g soil)'].to_dict()
    _iterativeIDs = load(open('iterativeIDs.json', 'r'))
    iterativeID_levels_3 = load(open('iterativeID_levels.json', 'r'))
    _iterativeID_taxonomy = {_iterativeIDs.get(_ID, _ID): phylogeny_1.get(_ID, 'Unknown') for _ID in _iterativeIDs}
    _unique_prefixes = {'_'.join(i.split('_')[:-1]) for i in dna_df_2.index}
    _sites = {i.split('_')[0] for i in _unique_prefixes}
    abundances_raw = read_csv('../data/Cliff_310MAG_relabund.tsv', sep='\t')
    abundances_depth = merge_depths_by_DNA_2(abundances_raw)
    abundances_2, condition = expand_triplicates_long(abundances_depth)
    _zero_level = 1e-05
    abundances_2 = abundances_2.loc[:, (abundances_2.fillna(0) > _zero_level).sum() >= 3]
    abundances_2.columns = [_iterativeIDs.get(_ID, _ID) for _ID in abundances_2.columns]
    abundances_2.rename(columns={'': 'unbinned'}, inplace=True)
    defined_1 = True
    if defined_1:
        abundances_2.drop(columns=[_c for _c in abundances_2.columns if _c not in _iterativeIDs.values()], inplace=True)
    reduced_1 = True
    if reduced_1:
        abundances_2.columns = ['.'.join(_ID.split('.')[:-1]) for _ID in abundances_2.columns]
        abundances_2 = abundances_2.groupby(abundances_2.columns, axis=1).sum()
        abundances_2.rename(columns={'': 'unbinned'}, inplace=True)
        iterativeID_phylums_2 = {'.'.join(_ID.split('.')[:-1]): _v for _ID, _v in iterativeID_phylums_1.items() if '.' in _ID}
    iterativeID_phylums_2 = {_ID: _v for _ID, _v in iterativeID_phylums_2.items() if _ID in abundances_2.columns}
    display(abundances_2)
    _abund_pseudo = abundances_2.values + 0.5
    _X = clr(_abund_pseudo)
    _taxa = abundances_2.columns.to_numpy()
    V = _X.shape[1]
    _R = _X.copy()
    for _c in np.unique(condition):
        _mask = condition == _c
        _R[_mask] = _R[_mask] - _R[_mask].mean(axis=0, keepdims=True)
    _C_marginal = np.corrcoef(_X.T)
    _C_within = np.corrcoef(_R.T)
    _X_means = np.vstack([_X[condition == _c].mean(0) for _c in np.unique(condition)])
    _C_between, _ = spearmanr(_X_means)
    _C_marginal = DataFrame(_C_marginal, index=_taxa, columns=_taxa)
    _C_within = DataFrame(_C_within, index=_taxa, columns=_taxa)
    _C_between = DataFrame(_C_between, index=_taxa, columns=_taxa)
    display(_C_within)
    EFFECT_THRESHOLD = 0.7
    _presence = (abundances_2 > _zero_level).astype(int)
    _cooccurrence = defaultdict(int)
    for _sample in _presence.itertuples(index=False):
        _present = [_col for _col, _val in zip(_presence.columns, _sample) if _val]
        for _pair in combinations(sorted(_present), 2):
            _cooccurrence[_pair] = _cooccurrence[_pair] + 1
    MIN_COOCCURRENCE = 6

    def within_condition_perm_null(x, y, condition):
        """Return array of |r| under all within-condition permutations of y."""
        cond_levels = np.unique(condition)
        cond_indices = [np.where(condition == _c)[0] for _c in cond_levels]
        null = []
        for perm_combo in product(*[list(permutations(idx)) for idx in cond_indices]):
            y_perm = y.copy()
            for orig_idx, new_idx in zip(cond_indices, perm_combo):
                y_perm[orig_idx] = y[list(new_idx)]
            x_resid = x - np.array([x[condition == _c].mean() for _c in condition])
            y_resid = y_perm - np.array([y_perm[condition == _c].mean() for _c in condition])
            denom = np.sqrt((x_resid ** 2).sum() * (y_resid ** 2).sum())
            null.append(abs((x_resid * y_resid).sum() / denom) if denom > 0 else 0.0)
        return np.array(null)
    pair_data_1 = []
    for (_a, _b), _count in _cooccurrence.items():
        if _count < MIN_COOCCURRENCE:
            continue
        if _a not in _C_within.index or _b not in _C_within.columns:
            continue
        r_within = _C_within.loc[_a, _b]
        if np.isnan(r_within) or abs(r_within) < EFFECT_THRESHOLD:
            continue
        r_marginal = _C_marginal.loc[_a, _b]
        r_between = _C_between.loc[_a, _b]
        null = within_condition_perm_null(_X[:, list(_taxa).index(_a)], _X[:, list(_taxa).index(_b)], condition)
        perm_p = (np.sum(null >= abs(r_within)) + 1) / (len(null) + 1)
        pair_data_1.append((_a, _b, r_within, r_marginal, r_between, _count, perm_p))
    print(f'Pairs surviving |r_within| > {EFFECT_THRESHOLD} and co-occurrence >= {MIN_COOCCURRENCE}: {len(pair_data_1)}')
    interactions_df = DataFrame(pair_data_1, columns=['taxon_a', 'taxon_b', 'r_within', 'r_marginal', 'r_between', 'cooccurrence', 'perm_p'])
    interactions_df.to_csv('direct_interactions.csv', index=False)
    display(interactions_df.head(20))
    iu = np.triu_indices(V, 1)
    _fig, _ax = pyplot.subplots(figsize=(5, 5))
    _ax.scatter(_C_marginal.values[iu], _C_within.values[iu], s=2, alpha=0.3)
    _ax.plot([-1, 1], [-1, 1], 'k--', linewidth=1)
    _ax.axhline(EFFECT_THRESHOLD, color='r', linestyle=':', linewidth=0.8)
    _ax.axhline(-EFFECT_THRESHOLD, color='r', linestyle=':', linewidth=0.8)
    _ax.set_xlabel('marginal $r$ (across 12 samples)')
    _ax.set_ylabel('within-condition $r$ (residuals)')
    _ax.set_title('Coupling vs. condition-response decomposition')
    pyplot.tight_layout()
    pyplot.savefig('marginal_vs_within.png', dpi=200)
    pyplot.show()
    try:
        _gl = GraphicalLassoCV(alphas=10, max_iter=200).fit(_R)
        _precision = DataFrame(_gl.precision_, index=_taxa, columns=_taxa)
        glasso_edges = []
        for i_3, _a in enumerate(_taxa):
            for _j, _b in enumerate(_taxa):
                if _j <= i_3:
                    continue
                if abs(_precision.iloc[i_3, _j]) > 1e-08:
                    p_corr = -_precision.iloc[i_3, _j] / np.sqrt(_precision.iloc[i_3, i_3] * _precision.iloc[_j, _j])
                    glasso_edges.append((_a, _b, p_corr))
        glasso_df = DataFrame(glasso_edges, columns=['taxon_a', 'taxon_b', 'partial_r'])
        glasso_df.to_csv('glasso_direct_edges.csv', index=False)
        print(f'Graphical lasso direct edges: {len(glasso_df)}')
        display(glasso_df.head(20))
    except Exception as e:
        print(f'Graphical lasso failed (likely V > n issue): {e}')
    _orgs = set()
    for _row in pair_data_1:
        _orgs.add(_row[0])
        _orgs.add(_row[1])
    dump(list(_orgs), open('significantly_connected_organisms.json', 'w'))
    level_2 = 'Phylum'
    archaea_phyla_1 = sorted({_phylum for _phylum in iterativeID_phylums_2.values() if _phylum and ('archaeo' in _phylum.lower() or any((_a in _phylum.lower() for _a in ['candidatus thermoplasmatota', 'halobacterota', 'methanobacteriota', 'micrarchaeota'])))})
    bacteria_phyla_1 = sorted({_phylum for _phylum in iterativeID_phylums_2.values() if _phylum and _phylum not in archaea_phyla_1})
    _n_archaea = len(archaea_phyla_1)
    _archaea_colors = [pyplot.cm.turbo(i / max(_n_archaea, 1) * 0.15) for i in range(_n_archaea)]
    _n_bacteria = len(bacteria_phyla_1)
    _bacteria_colors = [pyplot.cm.turbo(0.2 + i / max(_n_bacteria, 1) * 0.8) for i in range(_n_bacteria)]
    taxa_color_map_1 = {}
    for _phylum, _color in zip(archaea_phyla_1, _archaea_colors):
        taxa_color_map_1[_phylum] = _color
    for _phylum, _color in zip(bacteria_phyla_1, _bacteria_colors):
        taxa_color_map_1[_phylum] = _color
    _iterativeID_color_map = {_ID: taxa_color_map_1[_phylum] for _ID, _phylum in iterativeID_phylums_2.items() if _phylum}
    DEFAULT_COLOR_1 = 'lightgray'
    _bar_label = 'Phylum'
    row_colors_1 = Series({idx: _iterativeID_color_map.get(idx, DEFAULT_COLOR_1) for idx in _C_within.index}, name=_bar_label)
    col_colors_1 = Series({_col: _iterativeID_color_map.get(_col, DEFAULT_COLOR_1) for _col in _C_within.columns}, name=_bar_label)
    display(row_colors_1[:5])
    return (
        DEFAULT_COLOR_1,
        V,
        condition,
        defined_1,
        iterativeID_levels_3,
        iterativeID_phylums_2,
        reduced_1,
        taxa_color_map_1,
    )


@app.cell
def within_condition_correlation():
    """
    Within-condition correlation pipeline for 4 conditions x 3 replicates.

    Decomposes co-variation in the (CLR-transformed) abundance matrix into:
      wc_C_marginal  -- Pearson on CLR data across all 12 samples
                        (mixes coupling + condition response)
      wc_C_within    -- Pearson on residuals after subtracting condition means
                        (direct-coupling candidates)
      wc_C_between   -- Spearman on per-condition means (n=4)
                        (shared response to the condition gradient)

    Direct-coupling pairs filtered by:
      (a) effect size: |r_within| >= EFFECT_THRESHOLD
      (b) within-condition permutation p, BH-corrected to q < FDR_ALPHA

    Why permutation p with BH instead of parametric p with BH:
      parametric correlation p-values violate independence assumptions
      (V=140 taxa means each appears in 139 pairs), making BH anti-conservative.
      Permutation p-values are calibrated against the actual data structure,
      so BH on perm-p is much more defensible at this dependence level.

    Outputs to other cells:
      wc_abund_long, wc_condition, wc_taxa
      wc_X, wc_R
      wc_C_marginal, wc_C_within, wc_C_between
      wc_interactions_df  (effect-size + FDR filtered)
      wc_glasso_df        (sparse partial-correlation edges)
    """

    # Cell-local imports & helpers — function scope keeps the namespace clean;
    # only `wc_*` names below the function become cell-level outputs.
    def _wc_run(EFFECT_THRESHOLD=0.7, FDR_ALPHA=0.05, MIN_COOCC=6, ZERO=1e-5):
        import numpy as np
        from pandas import DataFrame, Series, read_csv
        from scipy.stats import spearmanr
        from json import load
        from itertools import combinations, permutations, product
        from collections import defaultdict
        from skbio.stats.composition import clr
        from sklearn.covariance import GraphicalLassoCV
        from statsmodels.stats.multitest import multipletests

        # 1. DNA-weighted depth merge
        dna_df = read_csv("../data/SaltPondsDNA.csv").set_index("Sample")
        sample_dna = dna_df["DNA conc (µg/g soil)"].to_dict()
        abund_raw = read_csv("../data/Cliff_310MAG_relabund.tsv", sep="\t")
        prefixes = {"_".join(s.split("_")[:-1]) for s in dna_df.index
                    if s.endswith(("_D1", "_D2"))}
        sites = sorted({p.split("_")[0] for p in prefixes})
        abund_depth = DataFrame({
            p: (abund_raw[f"{p}_D1"] * sample_dna[f"{p}_D1"]
              + abund_raw[f"{p}_D2"] * sample_dna[f"{p}_D2"])
               / (sample_dna[f"{p}_D1"] + sample_dna[f"{p}_D2"])
            for p in prefixes
        }, index=abund_raw.index)

        # 2. Long form (12 samples) + condition labels
        rows, cond = {}, []
        for site in sites:
            for rep in ("A", "B", "C"):
                sid = f"{site}_{rep}"
                if sid in abund_depth.columns:
                    rows[sid] = abund_depth[sid]
                    cond.append(site)
        abund_long = DataFrame(rows).T
        cond = np.array(cond)

        # 3. Filters: presence + defined IDs + suffix collapse
        abund_long = abund_long.loc[:, (abund_long.fillna(0) > ZERO).sum() >= MIN_COOCC]
        iter_ids = load(open("iterativeIDs.json"))
        abund_long.columns = [iter_ids.get(c, c) for c in abund_long.columns]
        abund_long = abund_long.rename(columns={"": "unbinned"})
        abund_long = abund_long.loc[:, [c for c in abund_long.columns
                                        if c in set(iter_ids.values())]]
        abund_long.columns = [".".join(c.split(".")[:-1]) for c in abund_long.columns]
        abund_long = abund_long.T.groupby(abund_long.columns).sum().T  # avoid axis=1 deprecation
        abund_long = abund_long.rename(columns={"": "unbinned"})

        # 4. CLR + three correlation matrices
        # Pseudocount calibrated to the data scale: the 0.5 default is for
        # count data; here values are relative abundances (0.001-0.5 range),
        # so 0.5 would dwarf real signal. Using half the minimum non-zero
        # abundance places zeros just below the detection limit without
        # distorting non-zero values.
        _vals = abund_long.values
        _nzmin = _vals[_vals > 0].min()
        _pseudocount = _nzmin / 2
        print(f'min non-zero abundance: {_nzmin:.3e}, pseudocount: {_pseudocount:.3e}')
        X = clr(_vals + _pseudocount)
        taxa = abund_long.columns.to_numpy()
        V = X.shape[1]
        C_marginal_arr = np.corrcoef(X.T)
        R = X.copy()
        for c in np.unique(cond):
            m = cond == c
            R[m] -= R[m].mean(axis=0, keepdims=True)
        C_within_arr = np.corrcoef(R.T)
        X_means = np.vstack([X[cond == c].mean(0) for c in np.unique(cond)])
        C_between_arr, _ = spearmanr(X_means)
        C_marginal = DataFrame(C_marginal_arr, index=taxa, columns=taxa)
        C_within   = DataFrame(C_within_arr,   index=taxa, columns=taxa)
        C_between  = DataFrame(C_between_arr,  index=taxa, columns=taxa)

        # 5. Co-occurrence
        presence = (abund_long > ZERO).astype(int)
        cooc = defaultdict(int)
        for sample in presence.itertuples(index=False):
            present = [c for c, v in zip(presence.columns, sample) if v]
            for pair in combinations(sorted(present), 2):
                cooc[pair] += 1

        # 6. Permutation p for every co-occurring pair (no effect-size pre-filter,
        #    so we have a complete set for BH-FDR)
        def perm_null(x, y, condition_arr):
            levels = np.unique(condition_arr)
            idxs   = [np.where(condition_arr == c)[0] for c in levels]
            null = []
            for combo in product(*[list(permutations(idx)) for idx in idxs]):
                yp = y.copy()
                for orig, new in zip(idxs, combo):
                    yp[orig] = y[list(new)]
                xr = x - np.array([x[condition_arr == c].mean() for c in condition_arr])
                yr = yp - np.array([yp[condition_arr == c].mean() for c in condition_arr])
                denom = np.sqrt((xr**2).sum() * (yr**2).sum())
                null.append(abs((xr * yr).sum() / denom) if denom > 0 else 0.0)
            return np.array(null)

        taxa_list = list(taxa)
        candidates = []
        for (a, b), cnt in cooc.items():
            if cnt < MIN_COOCC:                           continue
            if a not in C_within.index or b not in C_within.columns: continue
            rw = C_within.loc[a, b]
            if np.isnan(rw):                              continue
            rm = C_marginal.loc[a, b]
            rb = C_between.loc[a, b]
            null = perm_null(X[:, taxa_list.index(a)],
                             X[:, taxa_list.index(b)], cond)
            p = (np.sum(null >= abs(rw)) + 1) / (len(null) + 1)
            candidates.append((a, b, rw, rm, rb, cnt, p))

        cand_df = DataFrame(candidates,
            columns=["taxon_a","taxon_b","r_within","r_marginal","r_between",
                     "cooccurrence","perm_p"])

        # 7. BH-FDR on permutation p-values
        if len(cand_df) > 0:
            reject, qvals, _, _ = multipletests(cand_df["perm_p"].values,
                                                alpha=FDR_ALPHA, method="fdr_bh")
            cand_df["q_BH"]    = qvals
            cand_df["fdr_pass"] = reject
        else:
            cand_df["q_BH"], cand_df["fdr_pass"] = [], []

        cand_df["effect_pass"] = cand_df["r_within"].abs() >= EFFECT_THRESHOLD
        cand_df["both_pass"]   = cand_df["effect_pass"] & cand_df["fdr_pass"]
        interactions = cand_df[cand_df["both_pass"]] \
            .sort_values(["q_BH", "perm_p"]).reset_index(drop=True)
        interactions.to_csv("direct_interactions.csv", index=False)
        cand_df.to_csv("all_candidate_pairs.csv", index=False)

        # 7.5  maxT permutation: FWER-controlled |r| threshold for the
        #      whole matrix. Independently within-condition-shuffles each
        #      taxon's residuals (so the null preserves the design but
        #      destroys cross-taxon coupling), recomputes the correlation
        #      matrix, takes max |r|, and reports the 95th percentile of
        #      those maxima as a calibrated effect-size threshold.
        rng = np.random.default_rng(42)
        cond_idxs = [np.where(cond == c)[0] for c in np.unique(cond)]
        n_iter = 500
        iu = np.triu_indices(R.shape[1], 1)
        maxT = np.empty(n_iter)
        for it in range(n_iter):
            Rsh = R.copy()
            for col in range(Rsh.shape[1]):
                for idx in cond_idxs:
                    Rsh[idx, col] = Rsh[rng.permutation(idx), col]
            Cperm = np.corrcoef(Rsh.T)
            maxT[it] = np.max(np.abs(Cperm[iu]))
        maxT_95 = float(np.percentile(maxT, 95))
        maxT_99 = float(np.percentile(maxT, 99))
        print(f"maxT (FWER) noise floor: 95th-pct |r| under null = {maxT_95:.3f},"
              f" 99th-pct = {maxT_99:.3f}")
        print(f"  -> {EFFECT_THRESHOLD:.2f} sits at the "
              f"{(maxT < EFFECT_THRESHOLD).mean()*100:.1f} percentile of the null max-|r|")

        # 8. Diagnostic plot
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        iu = np.triu_indices(V, 1)
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(C_marginal_arr[iu], C_within_arr[iu], s=2, alpha=0.3)
        ax.plot([-1, 1], [-1, 1], "k--", lw=1)
        ax.axhline( EFFECT_THRESHOLD, color="r", ls=":", lw=0.8)
        ax.axhline(-EFFECT_THRESHOLD, color="r", ls=":", lw=0.8)
        ax.set_xlabel("marginal r (12 samples, Pearson on CLR)")
        ax.set_ylabel("within-condition r (residuals, Pearson)")
        ax.set_title("Coupling vs. condition-response decomposition")
        ax.set_aspect("equal"); fig.tight_layout()
        fig.savefig("marginal_vs_within.png", dpi=200); plt.close(fig)

        # 9. Sparse direct-interaction network via Graphical Lasso on residuals
        glasso_df = DataFrame(columns=["taxon_a","taxon_b","partial_r"])
        try:
            gl = GraphicalLassoCV(alphas=10, max_iter=200).fit(R)
            P = gl.precision_
            edges = []
            for i, a in enumerate(taxa):
                for j, b in enumerate(taxa):
                    if j <= i: continue
                    if abs(P[i, j]) > 1e-8:
                        pcorr = -P[i, j] / np.sqrt(P[i, i] * P[j, j])
                        edges.append((a, b, pcorr))
            glasso_df = DataFrame(edges, columns=["taxon_a","taxon_b","partial_r"]) \
                .sort_values("partial_r", key=lambda s: s.abs(), ascending=False)
            glasso_df.to_csv("glasso_direct_edges.csv", index=False)
        except Exception as exc:
            print(f"Graphical Lasso failed: {type(exc).__name__}: {exc}")

        # Persist matrices
        C_marginal.to_csv("C_marginal.csv")
        C_within.to_csv("C_within.csv")
        C_between.to_csv("C_between.csv")

        # Console summary
        print(f"Samples: {abund_long.shape[0]}, taxa (post-filter): {V}")
        print(f"Conditions: {dict((c, int((cond == c).sum())) for c in sorted(set(cond)))}")
        n_pairs = len(np.triu_indices(V, 1)[0])
        print(f"Upper-tri pairs: {n_pairs}")
        print(f"  |r_marginal| > {EFFECT_THRESHOLD}: "
              f"{int((np.abs(C_marginal_arr[np.triu_indices(V,1)]) > EFFECT_THRESHOLD).sum())}")
        print(f"  |r_within|   > {EFFECT_THRESHOLD}: "
              f"{int((np.abs(C_within_arr[np.triu_indices(V,1)])   > EFFECT_THRESHOLD).sum())}")
        print(f"Candidates tested (cooc>={MIN_COOCC}): {len(cand_df)}")
        print(f"  pass effect ({EFFECT_THRESHOLD}): {int(cand_df['effect_pass'].sum())}")
        print(f"  pass FDR    (q<{FDR_ALPHA}):    {int(cand_df['fdr_pass'].sum())}")
        print(f"  pass BOTH:                       {len(interactions)}")
        print(f"Graphical Lasso edges: {len(glasso_df)}")

        return {
            "abund_long": abund_long, "condition": cond, "taxa": taxa,
            "X": X, "R": R,
            "C_marginal": C_marginal, "C_within": C_within, "C_between": C_between,
            "interactions": interactions, "all_candidates": cand_df,
            "glasso": glasso_df,
            "maxT_samples": maxT,
            "maxT_95": maxT_95,
            "maxT_99": maxT_99,
        }

    _wc = _wc_run()
    wc_abund_long      = _wc["abund_long"]
    wc_condition       = _wc["condition"]
    wc_taxa            = _wc["taxa"]
    wc_X               = _wc["X"]
    wc_R               = _wc["R"]
    wc_C_marginal      = _wc["C_marginal"]
    wc_C_within        = _wc["C_within"]
    wc_C_between       = _wc["C_between"]
    wc_interactions_df = _wc["interactions"]
    wc_all_candidates  = _wc["all_candidates"]
    wc_glasso_df       = _wc["glasso"]
    wc_maxT_samples    = _wc["maxT_samples"]
    wc_maxT_95         = _wc["maxT_95"]
    wc_maxT_99         = _wc["maxT_99"]

    wc_interactions_df.head(20)

    return wc_C_within, wc_abund_long, wc_all_candidates, wc_glasso_df, wc_taxa


@app.cell
def within_condition_triangle_heatmap(wc_C_within, wc_glasso_df):
    """Triangle heatmap of within-condition correlations.

    After the calibrated CLR pseudocount, no pair passes BH-FDR (q<0.05) at
    this design's permutation-p floor (1/1297 vs 8 180 tests). So we mark
    the effect-size-passing cells (|r_within| >= 0.7) -- the chat's original
    recommendation given the inter-test dependence.
    """
    def _wc_render_heatmap():
        import numpy as np
        from numpy import triu, ones_like, nan
        from json import load
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib import patches
        import seaborn as sns

        LABEL_FONT       = 16
        HIGHLIGHT_FONT   = LABEL_FONT * 1.2
        TITLE_FONT       = 30
        STAR_FONT        = 18
        CBAR_LABEL_FONT  = 18
        CBAR_TICK_FONT   = 14
        STAR_COLOR       = "#90EE90"
        BOX_LW           = 2.0
        EFFECT           = 0.7

        iterativeIDs = load(open("iterativeIDs.json"))
        def _to_short(full_id):
            bare = full_id
            for suf in (".contigs__", ".contigs"):
                if bare.endswith(suf): bare = bare[: -len(suf)]; break
            iid = iterativeIDs.get(bare, bare)
            return ".".join(iid.split(".")[:-1]) if "." in iid else iid

        betaine_full = [
            "Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs",
            "Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs",
            "Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs",
            "Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs",
            "Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs",
            "Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs",
            "Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs",
            "Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs",
            "Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs",
            "Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs",
            "Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs",
            "Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs",
            "Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs",
            "Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs",
            "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs",
        ]
        mepe_full = [
            "Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.20.contigs",
            "Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.32.contigs",
            "Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.51.contigs",
            "Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_metabat.27.contigs",
            "Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_metabat.28.contigs",
            "Salt_Pond_MetaG_R1_B_H2O_MG_DASTool_bins_metabat.7.contigs",
            "Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_metabat.31.contigs",
            "Salt_Pond_MetaG_R1_C_H2O_MG_DASTool_bins_concoct_out.79.contigs",
            "Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.35.contigs",
            "Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.49.contigs",
            "Salt_Pond_MetaG_R2_A_D2_MG_DASTool_bins_metabat.14.contigs",
            "Salt_Pond_MetaG_R2_A_D2_MG_DASTool_bins_metabat.16.contigs",
            "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_concoct_out.143.contigs",
            "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_concoct_out.32.contigs",
            "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.10.contigs",
            "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.16.contigs",
            "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.19.contigs",
            "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.38.contigs",
            "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.55.contigs",
            "Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_concoct_out.88.contigs",
            "Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_concoct_out.95.contigs",
            "Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_metabat.17.contigs",
            "Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.18.contigs",
            "Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.45.contigs",
            "Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.47.contigs",
            "Salt_Pond_MetaG_R2_restored_C_black_MG_DASTool_bins_concoct_out.17.contigs",
            "Salt_Pond_MetaG_R2_restored_C_black_MG_DASTool_bins_metabat.8.contigs",
        ]
        betaine_short = {_to_short(x) for x in betaine_full}
        mepe_short    = {_to_short(x) for x in mepe_full}
        color_betaine = "#1f77b4"
        color_mepe    = "#e8671b"

        cmap_obj = sns.clustermap(
            wc_C_within,
            cmap="coolwarm_r", center=0, vmin=-1, vmax=1,
            figsize=(28, 28),
            dendrogram_ratio=(0.08, 0.08),
            xticklabels=True, yticklabels=True,
        )
        cmap_obj.ax_row_dendrogram.set_visible(False)
        cmap_obj.ax_col_dendrogram.set_visible(False)
        ax = cmap_obj.ax_heatmap
        ax.yaxis.set_ticks_position("left")
        ax.yaxis.set_label_position("left")

        for lbl in ax.get_xticklabels():
            t = lbl.get_text()
            lbl.set_rotation(90)
            if   t in betaine_short:
                lbl.set_color(color_betaine); lbl.set_fontweight("bold")
                lbl.set_fontsize(HIGHLIGHT_FONT)
            elif t in mepe_short:
                lbl.set_color(color_mepe);    lbl.set_fontweight("bold")
                lbl.set_fontsize(HIGHLIGHT_FONT)
            else:
                lbl.set_fontsize(LABEL_FONT)
        for lbl in ax.get_yticklabels():
            t = lbl.get_text()
            lbl.set_rotation(0)
            if   t in betaine_short:
                lbl.set_color(color_betaine); lbl.set_fontweight("bold")
                lbl.set_fontsize(HIGHLIGHT_FONT)
            elif t in mepe_short:
                lbl.set_color(color_mepe);    lbl.set_fontweight("bold")
                lbl.set_fontsize(HIGHLIGHT_FONT)
            else:
                lbl.set_fontsize(LABEL_FONT)

        row_order = cmap_obj.dendrogram_row.reordered_ind
        col_order = cmap_obj.dendrogram_col.reordered_ind
        df_reord = wc_C_within.iloc[row_order, col_order]
        mask = triu(ones_like(df_reord, dtype=bool), k=1)
        mesh = ax.collections[0]
        arr = mesh.get_array().reshape(df_reord.shape)
        arr[mask] = nan
        mesh.set_array(arr.ravel())

        # Mark cells whose pair appears as a direct interaction in the
        # Graphical Lasso network (|partial_r| > 0.05). This is the cross-figure
        # link between the dense Pearson view and the sparse direct-interaction
        # view — the cbar already encodes |r_within| magnitude, so the asterisks
        # now carry strictly new information.
        PARTIAL_FOR_ASTERISKS = 0.05
        direct_pairs = set()
        for _, _row in wc_glasso_df.iterrows():
            if abs(_row["partial_r"]) > PARTIAL_FOR_ASTERISKS:
                direct_pairs.add((_row["taxon_a"], _row["taxon_b"]))
                direct_pairs.add((_row["taxon_b"], _row["taxon_a"]))
        ordered_x = list(df_reord.columns)
        ordered_y = list(df_reord.index)
        nrows, ncols = df_reord.shape
        for i in range(nrows):
            for j in range(ncols):
                if mask[i, j] or i == j:
                    continue
                if (ordered_y[i], ordered_x[j]) in direct_pairs:
                    ax.text(j + 0.5, i + 0.5, "*",
                            ha="center", va="center",
                            color=STAR_COLOR, fontsize=STAR_FONT,
                            fontweight="bold", clip_on=False)

        # Box highlighted rows / cols
        ordered_taxa = list(df_reord.index)
        for i, taxon in enumerate(ordered_taxa):
            if   taxon in betaine_short: box_color = color_betaine
            elif taxon in mepe_short:    box_color = color_mepe
            else: continue
            ax.add_patch(patches.Rectangle((0, i), i + 1, 1,
                fill=False, edgecolor=box_color, linewidth=BOX_LW, clip_on=False))
            ax.add_patch(patches.Rectangle((i, i), 1, nrows - i,
                fill=False, edgecolor=box_color, linewidth=BOX_LW, clip_on=False))

        ax.set_xlabel(""); ax.set_ylabel("")
        ax.set_title("Within-condition Pearson r on CLR residuals  "
                     r"(asterisks: direct interaction in Graphical Lasso, "
                     r"$|r_{\rm partial}| > 0.05$)", fontsize=TITLE_FONT, pad=20)

        hm_pos = ax.get_position()
        cbar_w = hm_pos.width  * 0.025
        cbar_h = hm_pos.height * 0.30
        cbar_x = hm_pos.x0 + hm_pos.width  * 0.78
        cbar_y = hm_pos.y0 + hm_pos.height * 0.62
        cax = cmap_obj.cax
        cax.set_position([cbar_x, cbar_y, cbar_w, cbar_h])
        cax.tick_params(labelsize=CBAR_TICK_FONT)
        cax.set_xlabel(""); cax.set_ylabel("")
        cax.set_title("Pearson r", fontsize=CBAR_LABEL_FONT, pad=10)

        cmap_obj.figure.savefig("within_condition_triangle.png", dpi=600, bbox_inches="tight")
        plt.close(cmap_obj.figure)
        return cmap_obj.figure.get_size_inches().tolist()

    wc_heatmap_size = _wc_render_heatmap()
    print(f"Saved within_condition_triangle.png  ({wc_heatmap_size})")

    return


@app.cell
def within_condition_network(wc_abund_long, wc_all_candidates, wc_taxa):
    """Co-occurrence network from effect-size-passing pairs.

    Single giant component (140 nodes, 1179 edges) is broken up via Louvain
    community detection. Each community gets its own internal spring layout
    and is placed on a generous grid so clusters do not bleed together.
    """
    def _wc_render_network():
        import numpy as np
        import networkx as nx
        from networkx.algorithms import community as nx_comm
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from matplotlib import cm as mcm
        from matplotlib import patches as mpatches
        from json import load

        LABEL_MIN_ABUND = 0.001
        EFFECT = 0.7

        iterativeIDs = load(open("iterativeIDs.json"))
        phylogeny    = load(open("../data/Saltern_phylogeny.json"))
        def _to_short(full_id):
            bare = full_id
            for suf in (".contigs__", ".contigs"):
                if bare.endswith(suf): bare = bare[: -len(suf)]; break
            iid = iterativeIDs.get(bare, bare)
            return ".".join(iid.split(".")[:-1]) if "." in iid else iid

        id_to_phylum = {}
        for full_id, taxonomy in phylogeny.items():
            short = _to_short(full_id)
            ph = taxonomy.get("Phylum", "Unknown") if isinstance(taxonomy, dict) else "Unknown"
            ph = ph or "Unknown"
            if id_to_phylum.get(short, "Unknown") in ("", "Unknown", None):
                id_to_phylum[short] = ph

        present_phyla = sorted({id_to_phylum.get(t, "Unknown") for t in wc_taxa})
        archaea_keywords = ("archae", "halobacterota", "methanobacteriota",
                            "thermoplasmatota", "micrarchaeota", "asgardarchaeota",
                            "nanohaloarchaeota")
        archaea_phyla  = sorted([p for p in present_phyla
                                 if p and any(k in p.lower() for k in archaea_keywords)])
        bacteria_phyla = sorted([p for p in present_phyla
                                 if p and p not in archaea_phyla and p != "Unknown"])

        taxa_color_map = {}
        for i, ph in enumerate(archaea_phyla):
            taxa_color_map[ph] = plt.cm.turbo(0.05 + 0.20 * (i / max(len(archaea_phyla), 1)))
        for i, ph in enumerate(bacteria_phyla):
            taxa_color_map[ph] = plt.cm.turbo(0.30 + 0.65 * (i / max(len(bacteria_phyla), 1)))
        DEFAULT_COLOR = "lightgray"

        betaine_full = ["Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs","Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs","Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs","Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs","Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs","Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs","Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs","Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs","Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs","Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs","Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs","Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs","Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs","Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs","Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs"]
        mepe_full = ["Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.20.contigs","Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.32.contigs","Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.51.contigs","Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_metabat.27.contigs","Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_metabat.28.contigs","Salt_Pond_MetaG_R1_B_H2O_MG_DASTool_bins_metabat.7.contigs","Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_metabat.31.contigs","Salt_Pond_MetaG_R1_C_H2O_MG_DASTool_bins_concoct_out.79.contigs","Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.35.contigs","Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.49.contigs","Salt_Pond_MetaG_R2_A_D2_MG_DASTool_bins_metabat.14.contigs","Salt_Pond_MetaG_R2_A_D2_MG_DASTool_bins_metabat.16.contigs","Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_concoct_out.143.contigs","Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_concoct_out.32.contigs","Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.10.contigs","Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.16.contigs","Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.19.contigs","Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.38.contigs","Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.55.contigs","Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_concoct_out.88.contigs","Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_concoct_out.95.contigs","Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_metabat.17.contigs","Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.18.contigs","Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.45.contigs","Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.47.contigs","Salt_Pond_MetaG_R2_restored_C_black_MG_DASTool_bins_concoct_out.17.contigs","Salt_Pond_MetaG_R2_restored_C_black_MG_DASTool_bins_metabat.8.contigs"]
        betaine_short = {_to_short(x) for x in betaine_full}
        mepe_short    = {_to_short(x) for x in mepe_full}
        color_scheme = {"betaine": "#1f77b4", "mepe": "#e8671b"}

        cand = wc_all_candidates
        eff = cand[cand["effect_pass"]]
        G = nx.Graph()
        for _, row in eff.iterrows():
            G.add_edge(row["taxon_a"], row["taxon_b"], r_within=row["r_within"])
        if G.number_of_edges() == 0:
            print("No edges to plot."); return None
        G.remove_nodes_from(list(nx.isolates(G)))

        # ---- LOUVAIN community detection on a positive-only undirected graph
        Gpos = G.copy()
        for u, v, d in Gpos.edges(data=True):
            d["w"] = max(abs(d["r_within"]), 0.01)
        communities = nx_comm.louvain_communities(Gpos, weight="w", resolution=1.4, seed=42)
        communities = sorted(communities, key=len, reverse=True)
        print(f"Louvain communities: {len(communities)} (sizes: {[len(c) for c in communities[:8]]}...)")

        # ---- Per-community spring layout, very wide grid placement
        n_comm = len(communities)
        cols = int(np.ceil(np.sqrt(n_comm * 1.4)))   # bias toward landscape
        GRID = 18.0
        pos = {}
        for idx, comm in enumerate(communities):
            sub = G.subgraph(comm).copy()
            for u, v, d in sub.edges(data=True):
                d["lw"] = abs(d["r_within"]) ** 0.25
            if len(sub.nodes) <= 1:
                n = next(iter(comm))
                pos[n] = ((idx % cols) * GRID, (idx // cols) * GRID)
                continue
            sub_pos = nx.spring_layout(sub, seed=42, iterations=4000,
                                        k=3.0/max(np.sqrt(len(sub.nodes)),1),
                                        weight="lw")
            xs = np.array([p[0] for p in sub_pos.values()])
            ys = np.array([p[1] for p in sub_pos.values()])
            cx, cy = xs.mean(), ys.mean()
            scale = max((xs - cx).std(), (ys - cy).std(), 1e-6)
            offx = (idx % cols) * GRID
            offy = (idx // cols) * GRID
            for n, (x, y) in sub_pos.items():
                pos[n] = ((x - cx) / scale * 2.5 + offx,
                          (y - cy) / scale * 2.5 + offy)

        rel = wc_abund_long.div(wc_abund_long.sum(axis=1), axis=0).mean(axis=0)
        node_list   = list(G.nodes())
        NODE_BASE   = 800; NODE_SCALE = 28000
        sizes_pt2   = np.array([NODE_BASE + NODE_SCALE * np.sqrt(max(rel.get(n, 0), 0))
                                for n in node_list])
        diam_pt     = 2 * np.sqrt(sizes_pt2 / np.pi)
        node_phyla  = [id_to_phylum.get(n, "Unknown") for n in node_list]
        node_colors = [taxa_color_map.get(p, DEFAULT_COLOR) for p in node_phyla]

        fig, ax = plt.subplots(figsize=(34, 26))
        rhos    = np.array([d["r_within"] for _, _, d in G.edges(data=True)])
        widths  = 0.5 + 2.5 * np.abs(rhos)
        norm    = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
        edge_cmap = mcm.RdBu
        ecolors = [edge_cmap(norm(r)) for r in rhos]
        nx.draw_networkx_edges(G, pos, width=widths, edge_color=ecolors,
                               alpha=0.45, ax=ax)
        nx.draw_networkx_nodes(G, pos, nodelist=node_list,
                               node_size=sizes_pt2, node_color=node_colors,
                               edgecolors="black", linewidths=0.4, ax=ax)
        for n, d_pt in zip(node_list, diam_pt):
            if rel.get(n, 0) < LABEL_MIN_ABUND: continue
            x, y = pos[n]
            if   n in betaine_short: lbl_color, weight = color_scheme["betaine"], "bold"
            elif n in mepe_short:    lbl_color, weight = color_scheme["mepe"],    "bold"
            else:                    lbl_color, weight = "black",                 "normal"
            ax.text(x, y, str(n), fontsize=max(4, 0.45 * d_pt),
                    color=lbl_color, fontweight=weight,
                    ha="center", va="center", clip_on=False)

        ax.set_axis_off()
        ax.set_title("Co-variation pattern: "
                     r"$|r_{\rm within}| \geq 0.7$ "
                     "(effect size, not FWER-calibrated)",
                     fontsize=24, pad=20)

        sm = mcm.ScalarMappable(cmap=edge_cmap, norm=norm); sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.4, pad=0.01, location="right")
        cbar.set_label("within-condition Pearson r", fontsize=20)
        cbar.ax.tick_params(labelsize=16)
        cbar.ax.yaxis.set_label_position("left")
        # Park cbar in lower-right of the figure to free upper-right for the
        # phylum legend.
        cbar.ax.set_position([0.78, 0.08, 0.018, 0.30])

        def header_patch(title):
            return mpatches.Patch(color="none", label=fr"\$\\bf{{{title}}}\$")
        handles = []
        if archaea_phyla:
            handles.append(header_patch("Archaea"))
            handles += [mpatches.Patch(color=taxa_color_map[p], label=p) for p in archaea_phyla]
        if bacteria_phyla:
            handles.append(header_patch("Bacteria"))
            handles += [mpatches.Patch(color=taxa_color_map[p], label=p) for p in bacteria_phyla]
        leg1 = ax.legend(handles=handles, title="Phylum",
                         loc="upper left", bbox_to_anchor=(0.04, 0.84),
                         bbox_transform=fig.transFigure,
                         fontsize=16, title_fontsize=20, frameon=True)
        ax.add_artist(leg1)

        # Horizontal abundance legend: title above, circles in a row, labels below.
        legend_entries = [(0.001, "0.1%"), (0.01, "1%"), (0.05, "5%")]
        sizes_pt2_leg  = [NODE_BASE + NODE_SCALE * np.sqrt(a) for a, _ in legend_entries]
        diameters_pt   = [2 * np.sqrt(s / np.pi) for s in sizes_pt2_leg]
        breather_pt    = 18
        offsets_pt     = [diameters_pt[0] / 2]
        for k in range(1, len(legend_entries)):
            offsets_pt.append(offsets_pt[-1]
                              + (diameters_pt[k - 1] + diameters_pt[k]) / 2 + breather_pt)
        total_pt   = offsets_pt[-1] + diameters_pt[-1] / 2
        fig_w_pts  = fig.get_figwidth()  * 72
        fig_h_pts  = fig.get_figheight() * 72
        x_start    = 0.04           # legend left edge in figure fraction
        y_circles  = 0.91           # vertical center of the circles
        title_y    = y_circles + max(diameters_pt) / fig_h_pts / 2 + 0.01
        label_y    = y_circles - max(diameters_pt) / fig_h_pts / 2 - 0.018

        # Title centered over the circles
        title_x = x_start + (total_pt / fig_w_pts) / 2
        ax.text(title_x, title_y, "Mean rel. abundance",
                transform=fig.transFigure, ha="center", va="bottom",
                fontweight="bold", fontsize=18, clip_on=False)

        for (abund, label), s, off in zip(legend_entries, sizes_pt2_leg, offsets_pt):
            x = x_start + off / fig_w_pts
            ax.scatter([x], [y_circles], s=s, color="slategray", alpha=0.85,
                        edgecolor="black", linewidth=0.4,
                        transform=fig.transFigure, clip_on=False)
            ax.text(x, label_y, label,
                    transform=fig.transFigure, ha="center", va="top",
                    fontsize=16, clip_on=False)

        fig.savefig("within_condition_network.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
        return (G.number_of_nodes(), G.number_of_edges(), len(communities),
                sum(1 for n in node_list if n in betaine_short),
                sum(1 for n in node_list if n in mepe_short),
                sum(1 for n in node_list if rel.get(n, 0) >= LABEL_MIN_ABUND))

    wc_network_stats = _wc_render_network()
    print(f"Saved within_condition_network.png  (nodes, edges, communities, betaine, mepe, labeled = {wc_network_stats})")

    return


@app.cell
def within_condition_glasso_network(wc_abund_long, wc_glasso_df, wc_taxa):
    """Graphical Lasso direct-interaction network.

    Each edge is a non-zero entry of the L1-regularized inverse covariance
    (precision) matrix on within-condition residuals. Edge weight = partial
    correlation: the conditional dependency between two taxa after
    controlling for ALL other taxa simultaneously.

    This is the inferentially defensible network at this design:
    - The Lasso provides sparsity automatically (no effect-size threshold needed)
    - alpha is chosen via cross-validation by GraphicalLassoCV
    - partial correlations remove indirect/transitive coupling that the
      marginal/within-condition matrix doesn't separate

    Edges with positive partial r = direct positive coupling.
    Edges with negative partial r = direct negative coupling (niche partitioning).
    """
    def _wc_render_glasso_network():
        import numpy as np
        import networkx as nx
        from networkx.algorithms import community as nx_comm
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from matplotlib import cm as mcm
        from matplotlib import patches as mpatches
        from json import load

        LABEL_MIN_ABUND = 0.001

        iterativeIDs = load(open("iterativeIDs.json"))
        phylogeny    = load(open("../data/Saltern_phylogeny.json"))
        def _to_short(full_id):
            bare = full_id
            for suf in (".contigs__", ".contigs"):
                if bare.endswith(suf): bare = bare[: -len(suf)]; break
            iid = iterativeIDs.get(bare, bare)
            return ".".join(iid.split(".")[:-1]) if "." in iid else iid

        id_to_phylum = {}
        for full_id, taxonomy in phylogeny.items():
            short = _to_short(full_id)
            ph = taxonomy.get("Phylum", "Unknown") if isinstance(taxonomy, dict) else "Unknown"
            ph = ph or "Unknown"
            if id_to_phylum.get(short, "Unknown") in ("", "Unknown", None):
                id_to_phylum[short] = ph

        present_phyla = sorted({id_to_phylum.get(t, "Unknown") for t in wc_taxa})
        archaea_keywords = ("archae", "halobacterota", "methanobacteriota",
                            "thermoplasmatota", "micrarchaeota", "asgardarchaeota",
                            "nanohaloarchaeota")
        archaea_phyla  = sorted([p for p in present_phyla
                                 if p and any(k in p.lower() for k in archaea_keywords)])
        bacteria_phyla = sorted([p for p in present_phyla
                                 if p and p not in archaea_phyla and p != "Unknown"])

        taxa_color_map = {}
        for i, ph in enumerate(archaea_phyla):
            taxa_color_map[ph] = plt.cm.turbo(0.05 + 0.20 * (i / max(len(archaea_phyla), 1)))
        for i, ph in enumerate(bacteria_phyla):
            taxa_color_map[ph] = plt.cm.turbo(0.30 + 0.65 * (i / max(len(bacteria_phyla), 1)))
        DEFAULT_COLOR = "lightgray"

        betaine_full = ["Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs",
                        "Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs",
                        "Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs",
                        "Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27.contigs",
                        "Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45.contigs",
                        "Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.48.contigs",
                        "Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.50.contigs",
                        "Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_concoct_out.59.contigs",
                        "Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.58.contigs",
                        "Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.73.contigs",
                        "Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_maxbin.047.contigs",
                        "Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_concoct_out.85.contigs",
                        "Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40.contigs",
                        "Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.98.contigs",
                        "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.52.contigs"]
        mepe_full = ["Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.20.contigs",
                     "Salt_Pond_MetaG_R1_B_D1_MG_DASTool_bins_metabat.32.contigs",
                     "Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_concoct_out.51.contigs",
                     "Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_metabat.27.contigs",
                     "Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_metabat.28.contigs",
                     "Salt_Pond_MetaG_R1_B_H2O_MG_DASTool_bins_metabat.7.contigs",
                     "Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_metabat.31.contigs",
                     "Salt_Pond_MetaG_R1_C_H2O_MG_DASTool_bins_concoct_out.79.contigs",
                     "Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.35.contigs",
                     "Salt_Pond_MetaG_R2_A_D1_MG_DASTool_bins_concoct_out.49.contigs",
                     "Salt_Pond_MetaG_R2_A_D2_MG_DASTool_bins_metabat.14.contigs",
                     "Salt_Pond_MetaG_R2_A_D2_MG_DASTool_bins_metabat.16.contigs",
                     "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_concoct_out.143.contigs",
                     "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_concoct_out.32.contigs",
                     "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.10.contigs",
                     "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.16.contigs",
                     "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.19.contigs",
                     "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.38.contigs",
                     "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.55.contigs",
                     "Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_concoct_out.88.contigs",
                     "Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_concoct_out.95.contigs",
                     "Salt_Pond_MetaG_R2_C_D1_MG_DASTool_bins_metabat.17.contigs",
                     "Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.18.contigs",
                     "Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.45.contigs",
                     "Salt_Pond_MetaG_R2_C_D2_MG_DASTool_bins_concoct_out.47.contigs",
                     "Salt_Pond_MetaG_R2_restored_C_black_MG_DASTool_bins_concoct_out.17.contigs",
                     "Salt_Pond_MetaG_R2_restored_C_black_MG_DASTool_bins_metabat.8.contigs"]
        betaine_short = {_to_short(x) for x in betaine_full}
        mepe_short    = {_to_short(x) for x in mepe_full}
        color_scheme = {"betaine": "#1f77b4", "mepe": "#e8671b"}

        # ----- Build graph from Graphical Lasso edges -------------------
        PARTIAL_THRESHOLD = 0.05
        filtered = wc_glasso_df[wc_glasso_df["partial_r"].abs() > PARTIAL_THRESHOLD]
        print(f"|partial_r| > {PARTIAL_THRESHOLD}: kept {len(filtered)} of "
              f"{len(wc_glasso_df)} edges")
        G = nx.Graph()
        for _, row in filtered.iterrows():
            G.add_edge(row["taxon_a"], row["taxon_b"], partial_r=row["partial_r"])
        if G.number_of_edges() == 0:
            print("No Graphical Lasso edges to plot."); return None
        # Drop nodes that ended up with no surviving edges
        G.remove_nodes_from(list(nx.isolates(G)))

        # Louvain communities on positive-weight graph
        Gpos = G.copy()
        for u, v, d in Gpos.edges(data=True):
            d["w"] = max(abs(d["partial_r"]), 0.01)
        communities = nx_comm.louvain_communities(Gpos, weight="w", resolution=1.4, seed=42)
        communities = sorted(communities, key=len, reverse=True)
        print(f"Glasso Louvain communities: {len(communities)} "
              f"(sizes: {[len(c) for c in communities[:8]]}...)")

        # Per-community spring layout
        n_comm = len(communities)
        cols = int(np.ceil(np.sqrt(n_comm * 1.4)))
        GRID = 18.0
        pos = {}
        for idx, comm in enumerate(communities):
            sub = G.subgraph(comm).copy()
            for u, v, d in sub.edges(data=True):
                d["lw"] = abs(d["partial_r"]) ** 0.25
            if len(sub.nodes) <= 1:
                n = next(iter(comm))
                pos[n] = ((idx % cols) * GRID, (idx // cols) * GRID)
                continue
            sub_pos = nx.spring_layout(sub, seed=42, iterations=4000,
                                        k=3.0/max(np.sqrt(len(sub.nodes)),1),
                                        weight="lw")
            xs = np.array([p[0] for p in sub_pos.values()])
            ys = np.array([p[1] for p in sub_pos.values()])
            cx, cy = xs.mean(), ys.mean()
            scale = max((xs - cx).std(), (ys - cy).std(), 1e-6)
            offx = (idx % cols) * GRID
            offy = (idx // cols) * GRID
            for n, (x, y) in sub_pos.items():
                pos[n] = ((x - cx) / scale * 2.5 + offx,
                          (y - cy) / scale * 2.5 + offy)

        rel = wc_abund_long.div(wc_abund_long.sum(axis=1), axis=0).mean(axis=0)
        node_list   = list(G.nodes())
        NODE_BASE   = 800; NODE_SCALE = 28000
        sizes_pt2   = np.array([NODE_BASE + NODE_SCALE * np.sqrt(max(rel.get(n, 0), 0))
                                for n in node_list])
        diam_pt     = 2 * np.sqrt(sizes_pt2 / np.pi)
        node_phyla  = [id_to_phylum.get(n, "Unknown") for n in node_list]
        node_colors = [taxa_color_map.get(p, DEFAULT_COLOR) for p in node_phyla]

        fig, ax = plt.subplots(figsize=(34, 26))
        pcorrs  = np.array([d["partial_r"] for _, _, d in G.edges(data=True)])
        p_max   = max(np.abs(pcorrs).max(), 1e-6)
        widths  = 0.6 + 4.0 * np.abs(pcorrs) / p_max
        norm    = mcolors.TwoSlopeNorm(vmin=-p_max, vcenter=0, vmax=p_max)
        edge_cmap = mcm.RdBu
        ecolors = [edge_cmap(norm(r)) for r in pcorrs]
        nx.draw_networkx_edges(G, pos, width=widths, edge_color=ecolors,
                               alpha=0.7, ax=ax)
        nx.draw_networkx_nodes(G, pos, nodelist=node_list,
                               node_size=sizes_pt2, node_color=node_colors,
                               edgecolors="black", linewidths=0.4, ax=ax)
        for n, d_pt in zip(node_list, diam_pt):
            if rel.get(n, 0) < LABEL_MIN_ABUND: continue
            x, y = pos[n]
            if   n in betaine_short: lbl_color, weight = color_scheme["betaine"], "bold"
            elif n in mepe_short:    lbl_color, weight = color_scheme["mepe"],    "bold"
            else:                    lbl_color, weight = "black",                 "normal"
            ax.text(x, y, str(n), fontsize=max(4, 0.45 * d_pt),
                    color=lbl_color, fontweight=weight,
                    ha="center", va="center", clip_on=False)

        ax.set_axis_off()
        ax.set_title("Direct interactions: Graphical Lasso partial correlations\n"
                     "(within-condition residuals, L1 cross-validated)",
                     fontsize=22, pad=20)

        sm = mcm.ScalarMappable(cmap=edge_cmap, norm=norm); sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.4, pad=0.01, location="right")
        cbar.set_label("partial correlation (controlling for all other taxa)",
                        fontsize=20)
        cbar.ax.tick_params(labelsize=16)
        cbar.ax.yaxis.set_label_position("left")
        # Park cbar in lower-right of the figure to free upper-right for the
        # phylum legend.
        cbar.ax.set_position([0.78, 0.08, 0.018, 0.30])

        def header_patch(title):
            return mpatches.Patch(color="none", label=fr"$\bf{{{title}}}$")
        handles = []
        if archaea_phyla:
            handles.append(header_patch("Archaea"))
            handles += [mpatches.Patch(color=taxa_color_map[p], label=p) for p in archaea_phyla]
        if bacteria_phyla:
            handles.append(header_patch("Bacteria"))
            handles += [mpatches.Patch(color=taxa_color_map[p], label=p) for p in bacteria_phyla]
        leg1 = ax.legend(handles=handles, title="Phylum",
                         loc="upper left", bbox_to_anchor=(0.04, 0.84),
                         bbox_transform=fig.transFigure,
                         fontsize=16, title_fontsize=20, frameon=True)
        ax.add_artist(leg1)

        # Horizontal abundance legend: title above, circles in a row, labels below.
        legend_entries = [(0.001, "0.1%"), (0.01, "1%"), (0.05, "5%")]
        sizes_pt2_leg  = [NODE_BASE + NODE_SCALE * np.sqrt(a) for a, _ in legend_entries]
        diameters_pt   = [2 * np.sqrt(s / np.pi) for s in sizes_pt2_leg]
        breather_pt    = 18
        offsets_pt     = [diameters_pt[0] / 2]
        for k in range(1, len(legend_entries)):
            offsets_pt.append(offsets_pt[-1]
                              + (diameters_pt[k - 1] + diameters_pt[k]) / 2 + breather_pt)
        total_pt   = offsets_pt[-1] + diameters_pt[-1] / 2
        fig_w_pts  = fig.get_figwidth()  * 72
        fig_h_pts  = fig.get_figheight() * 72
        x_start    = 0.04           # legend left edge in figure fraction
        y_circles  = 0.91           # vertical center of the circles
        title_y    = y_circles + max(diameters_pt) / fig_h_pts / 2 + 0.01
        label_y    = y_circles - max(diameters_pt) / fig_h_pts / 2 - 0.018

        # Title centered over the circles
        title_x = x_start + (total_pt / fig_w_pts) / 2
        ax.text(title_x, title_y, "Mean rel. abundance",
                transform=fig.transFigure, ha="center", va="bottom",
                fontweight="bold", fontsize=18, clip_on=False)

        for (abund, label), s, off in zip(legend_entries, sizes_pt2_leg, offsets_pt):
            x = x_start + off / fig_w_pts
            ax.scatter([x], [y_circles], s=s, color="slategray", alpha=0.85,
                        edgecolor="black", linewidth=0.4,
                        transform=fig.transFigure, clip_on=False)
            ax.text(x, label_y, label,
                    transform=fig.transFigure, ha="center", va="top",
                    fontsize=16, clip_on=False)

        fig.savefig("within_condition_glasso_network.png", dpi=600, bbox_inches="tight")
        plt.close(fig)
        return (G.number_of_nodes(), G.number_of_edges(),
                int((pcorrs > 0).sum()), int((pcorrs < 0).sum()),
                len(communities))

    wc_glasso_stats = _wc_render_glasso_network()
    print(f"Saved within_condition_glasso_network.png  "
          f"(nodes, edges, positive_pcorr, negative_pcorr, communities = {wc_glasso_stats})")

    return


@app.cell
def _(
    DEFAULT_COLOR_1,
    G,
    compressor,
    defaultdict,
    defined_1,
    iterativeID_levels_3,
    iterativeID_phylums_2,
    mean_rel_abund,
    pos_1,
    reduced_1,
    taxa_color_map_1,
):
    import plotly.graph_objects as go

    def edge_trace(edges, name, color_fn, width_scale=3):
        _xs, _ys, hovers = ([], [], [])
        for _u, _v, _d in edges:
            x0, y0 = pos_1[_u]
            x1, y1 = pos_1[_v]
            _xs = _xs + [x0, x1, None]
            _ys = _ys + [y0, y1, None]
        return go.Scatter(x=_xs, y=_ys, mode='lines', line=dict(width=2, color=color_fn), hoverinfo='none', name=name)
    pos_edges = [(_u, _v, _d) for _u, _v, _d in G.edges(data=True) if _d['rho'] >= 0]
    neg_edges = [(_u, _v, _d) for _u, _v, _d in G.edges(data=True) if _d['rho'] < 0]
    edge_traces = [edge_trace(pos_edges, 'positive ρ', 'rgba(0,0,200,0.45)'), edge_trace(neg_edges, 'negative ρ', 'rgba(200,0,0,0.45)')]
    nodes_by_phylum = defaultdict(list)
    for _n in G.nodes():
        nodes_by_phylum[iterativeID_phylums_2.get(_n, 'Unknown')].append(_n)

    def rgb_str(rgba):
        r, g, _b, _ = rgba
        return f'rgb({int(r * 255)},{int(g * 255)},{int(_b * 255)})'
    _node_traces = []
    for _phylum, _ns in nodes_by_phylum.items():
        _color = taxa_color_map_1.get(_phylum, DEFAULT_COLOR_1)
        if isinstance(_color, tuple):
            _color = rgb_str(_color)
        _sizes = [40 * compressor(mean_rel_abund.get(_n, 0)) + 5 for _n in _ns]
        _hover = [f"<b>{_n}</b><br>phylum: {iterativeID_phylums_2.get(_n, 'Unknown')}<br>rank: {iterativeID_levels_3.get(_n, '?')}<br>mean rel. abund.: {mean_rel_abund.get(_n, 0):.4f}<br>degree: {G.degree(_n)}" for _n in _ns]
        _node_traces.append(go.Scatter(x=[pos_1[_n][0] for _n in _ns], y=[pos_1[_n][1] for _n in _ns], mode='markers+text', text=[str(_n) for _n in _ns], textposition='middle center', textfont=dict(size=9, color='black'), marker=dict(size=_sizes, color=_color, line=dict(width=1, color='white'), opacity=0.9), name=_phylum, hovertext=_hover, hoverinfo='text'))
    _fig_int = go.Figure(data=edge_traces + _node_traces)
    _fig_int.update_layout(title=f'Co-occurrence network (FDR q<0.05) — {G.number_of_nodes()} nodes, {G.number_of_edges()} edges', showlegend=True, hovermode='closest', xaxis=dict(visible=False), yaxis=dict(visible=False, scaleanchor='x', scaleratio=1), plot_bgcolor='white', width=1400, height=1000, legend=dict(itemclick='toggle', itemdoubleclick='toggleothers'))
    _fig_int.write_html(f"cooccurrence_network_p_value_FDR{('_reduced' if reduced_1 else '')}{('_defined' if defined_1 else '')}.html")
    _fig_int.show()
    return go, nodes_by_phylum


@app.cell
def _(
    DEFAULT_COLOR_1,
    G,
    compressor,
    defaultdict,
    defined_1,
    go,
    json,
    mean_rel_abund,
    nodes_by_phylum,
    pos_1,
    reduced_1,
    taxa_color_map_1,
):
    all_edge_x, all_edge_y = ([], [])
    edge_endpoints = []
    for _u, _v, _d in G.edges(data=True):
        x0, y0 = pos_1[_u]
        x1, y1 = pos_1[_v]
        all_edge_x = all_edge_x + [x0, x1, None]
        all_edge_y = all_edge_y + [y0, y1, None]
        edge_endpoints.append((_u, _v))
    bg_edges = go.Scatter(x=all_edge_x, y=all_edge_y, mode='lines', line=dict(width=1.2, color='rgba(150,150,150,0.35)'), hoverinfo='none', name='all edges')
    hl_edges = go.Scatter(x=[], y=[], mode='lines', line=dict(width=4, color='rgba(255,140,0,0.9)'), hoverinfo='none', name='selected node — edges', showlegend=True)
    _node_traces = []
    for _phylum, _ns in nodes_by_phylum.items():
        _color = taxa_color_map_1.get(_phylum, DEFAULT_COLOR_1)
        if isinstance(_color, tuple):
            _color = f'rgb({int(_color[0] * 255)},{int(_color[1] * 255)},{int(_color[2] * 255)})'
        _sizes = [40 * compressor(mean_rel_abund.get(_n, 0)) + 5 for _n in _ns]
        _hover = [f'<b>{_n}</b><br>phylum: {_phylum}<br>abund: {mean_rel_abund.get(_n, 0):.4f}<br>degree: {G.degree(_n)}' for _n in _ns]
        _node_traces.append(go.Scatter(x=[pos_1[_n][0] for _n in _ns], y=[pos_1[_n][1] for _n in _ns], mode='markers+text', text=[str(_n) for _n in _ns], textposition='middle center', textfont=dict(size=9, color='black'), customdata=[[str(_n)] for _n in _ns], marker=dict(size=_sizes, color=_color, line=dict(width=1, color='white'), opacity=0.9), name=_phylum, hovertext=_hover, hoverinfo='text'))
    _fig_int = go.Figure(data=[bg_edges, hl_edges] + _node_traces)
    _fig_int.update_layout(title=f'Co-occurrence network (FDR q<0.05) — {G.number_of_nodes()} nodes, {G.number_of_edges()} edges', showlegend=True, hovermode='closest', xaxis=dict(visible=False), yaxis=dict(visible=False, scaleanchor='x', scaleratio=1), plot_bgcolor='white', width=1400, height=1000)
    edges_by_node = defaultdict(list)
    for i_4, (_u, _v) in enumerate(edge_endpoints):
        edges_by_node[str(_u)].append(i_4)
        edges_by_node[str(_v)].append(i_4)
    edges_by_node = {_k: _v for _k, _v in edges_by_node.items()}
    edge_segments = [[list(pos_1[_u]), list(pos_1[_v])] for _u, _v in edge_endpoints]
    post_script = f"                                                                                                                                                         \nconst div = document.getElementById('{{plot_id}}');             \nconst edgesByNode = {json.dumps(edges_by_node)};                                                                                                                           \nconst edgeSegments = {json.dumps(edge_segments)};                                                                                                                          \nconst HIGHLIGHT_TRACE_INDEX = 1;   // hl_edges is the second trace                                                                                                         \n                                                                                                                                                                            \nfunction highlight(nodeId) {{                                   \n    const idxs = edgesByNode[nodeId] || [];                                                                                                                                \n    const xs = [], ys = [];                                                                                                                                                \n    for (const i of idxs) {{\n        const [[x0, y0], [x1, y1]] = edgeSegments[i];                                                                                                                      \n        xs.push(x0, x1, null);                                                                                                                                             \n        ys.push(y0, y1, null);\n    }}                                                                                                                                                                     \n    Plotly.restyle(div, {{x: [xs], y: [ys]}}, [HIGHLIGHT_TRACE_INDEX]);\n}}                                                                                                                                                                         \n                                                                \ndiv.on('plotly_click', function(data) {{                                                                                                                                   \n    const pt = data.points && data.points[0];                                                                                                                              \n    if (!pt || !pt.customdata) return;\n    highlight(pt.customdata[0]);                                                                                                                                           \n}});                                                                                                                                                                       \n\n// Click on empty plot area to clear the highlight                                                                                                                         \ndiv.on('plotly_doubleclick', function() {{                                                                                                                                 \n    Plotly.restyle(div, {{x: [[]], y: [[]]}}, [HIGHLIGHT_TRACE_INDEX]);\n}});                                                                                                                                                                       \n"
    _fig_int.write_html(f"cooccurrence_network_p_value_FDR{('_reduced' if reduced_1 else '')}{('_defined' if defined_1 else '')}.html", post_script=post_script, include_plotlyjs='cdn', full_html=True)
    _fig_int.show()
    return


@app.cell
def _(
    DEFAULT_COLOR_1,
    G,
    compressor,
    defined_1,
    iterativeID_phylums_2,
    mean_rel_abund,
    reduced_1,
    taxa_color_map_1,
):
    from pyvis.network import Network
    net = Network(height='900px', width='100%', bgcolor='#ffffff', font_color='black', cdn_resources='in_line')
    net.barnes_hut(spring_length=200, spring_strength=0.005, central_gravity=0.1, gravity=-2000, damping=0.4)
    for _n in G.nodes():
        _phylum = iterativeID_phylums_2.get(_n, 'Unknown')
        _color = taxa_color_map_1.get(_phylum, DEFAULT_COLOR_1)
        if isinstance(_color, tuple):
            _color = f'#{int(_color[0] * 255):02x}{int(_color[1] * 255):02x}{int(_color[2] * 255):02x}'
        net.add_node(_n, label=str(_n), title=f'{_n}<br>phylum: {_phylum}<br>abund: {mean_rel_abund.get(_n, 0):.4f}', size=10 + 80 * compressor(mean_rel_abund.get(_n, 0)), color=_color)
    for _u, _v, _d in G.edges(data=True):
        net.add_edge(_u, _v, value=abs(_d['rho']), color='#1f77b4' if _d['rho'] >= 0 else '#d62728', title=f"ρ={_d['rho']:.2f}, q={_d['qvalue']:.2g}")
    net.show_buttons(filter_=['physics'])
    net.write_html(f"cooccurrence_network_pyvis{('_reduced' if reduced_1 else '')}{('_defined' if defined_1 else '')}.html", notebook=False, open_browser=False)
    return


if __name__ == "__main__":
    app.run()
