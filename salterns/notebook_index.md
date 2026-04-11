# Saltern Notebooks Index

### Root level
| Path | Purpose |
|------|---------|
| `AGORACommunityModeling.ipynb` | Builds probabilistic community metabolic models from MAG annotations with ATP correction |
| `BuildCliffCommModels.ipynb` | Sets up KBase environment, builds initial models from MAG annotations with ANI thresholds |
| `BuildCliffCommModelsNew.ipynb` | Updated community model building workflow |
| `BuildCliffCommModelsNew-Copy1.ipynb` | Copy variant of above |
| `BuildCliffCommModelsNew-Copy2.ipynb` | Copy variant of above |
| `BuildCliffCommModelsNewAF.ipynb` | Another variant of the community model building pipeline |
| `DNA_processing.ipynb` | Correlates DNA concentrations with methane flux, normalizes by DNA density |
| `QA:QC parsing.ipynb` | QA/QC checks on gapfilled model content (e.g., cpd25960 presence) |
| `SMIPPs_AF.ipynb` | Simulates SMIPP phenotypes (uptake/excretion/growth) with probabilistic gapfilling |
| `SMIPPs_CH.ipynb` | Parallel processing variant of SMIPP phenotype simulations |
| `SalternStatistics.ipynb` | Statistical analysis of saltern samples |
| `Visualizations.ipynb` | Dataset visualizations |
| `analysis.ipynb` | General analysis notebook |
| `checkingFeatures.ipynb` | Compares genome annotation features between RAST and GLM4EC methods |
| `nameIDmapping_abundances.ipynb` | Maps metabolomics compounds to ModelSEED IDs, prepares MAG abundance data |

### cliff_mags/
| Path | Purpose |
|------|---------|
| `cliff_mags/ANI.ipynb` | Computes FastANI comparisons between 450 genomes and GTDB references |
| `cliff_mags/ani_graph.ipynb` | Processes and visualizes ANI results |
| `cliff_mags/annotation.ipynb` | Annotates protein sequences using RAST, builds super-model |
| `cliff_mags/core_reconstruction.ipynb` | Analyzes core genome composition at different ANI thresholds (95/85/75%) |
| `cliff_mags/model_analysis.ipynb` | Analyzes metabolic models |
| `cliff_mags/model_analysis-Copy1.ipynb` | Copy variant of above |
| `cliff_mags/motu.ipynb` | Loads core genomic data, analyzes pangenome range at ANI cutoffs |
| `cliff_mags/motu_model.ipynb` | Compares core genome annotations, constructs models with pangenome features |
| `cliff_mags/pangenome_analysis.ipynb` | Pangenome analysis with/without limits across ANI thresholds |
| `cliff_mags/prob_model.ipynb` | Builds probabilistic models from annotations using ANI-based probabilities |
| `cliff_mags/rast.ipynb` | Batch annotates MAG proteins using RAST, stores in MongoDB |
| `cliff_mags/read_function.ipynb` | Reads and processes functional annotation data |
| `cliff_mags/read_functions_clean.ipynb` | Cleaned version of function processing |
| `cliff_mags/reads.ipynb` | Read mapping and processing |
| `cliff_mags/reads_split.ipynb` | Splits paired-end reads, maps to assemblies using CoverM |

### cliff_mags/mac/salterns/
| Path | Purpose |
|------|---------|
| `cliff_mags/mac/salterns/analysis.ipynb` | Loads genomes from KBase, indexes in MongoDB |
| `cliff_mags/mac/salterns/bet_reducers.ipynb` | Betaine reducer analysis |
| `cliff_mags/mac/salterns/build_models.ipynb` | Model building workflow |
| `cliff_mags/mac/salterns/gapfill.ipynb` | Gapfilling workflow |
| `cliff_mags/mac/salterns/index_genomes.ipynb` | Indexes 450+ genomes from KBase workspace |
| `cliff_mags/mac/salterns/mepn.ipynb` | MePn (methylphosphonate) analysis |

### correlational/
| Path | Purpose |
|------|---------|
| `correlational/abundances.ipynb` | Loads taxonomy/MAG abundances, aggregates by family, normalizes |
| `correlational/analysis(1).ipynb` | Flux-methane correlations, abundance-methane regressions, visualization |
| `correlational/bet_reducers_local.ipynb` | Local betaine reducer analysis |
| `correlational/bet_reducers_local copy.ipynb` | Copy of above |
| `correlational/betaine_gapfilling.ipynb` | Gapfilling + betaine reductase addition for 15 betaine MAG models |
| `correlational/test_bet.ipynb` | Test ATP correction and gapfilling on single MAG model (concoct_out.9) |
