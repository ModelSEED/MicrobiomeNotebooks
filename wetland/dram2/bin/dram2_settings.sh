# Shared settings for every DRAM2 launch in this run tree (sourced by run_dram2.sh).
B=/scratch1/afreiburger/wetland_dram2                  # run tree (NFS)
D=/scratch1/fliu/hub_scratch/chenry/DRAM2              # chenry's install (READ-ONLY)
PIPELINE=$D/repo/main.nf                               # DRAM2 v2.0.0-beta17, revision 7754dfaeaf
WORK_ROOT=/scratch/afreiburger/wetland_dram2_work      # local xfs for Nextflow work dirs
MAG_DIR=/scratch1/afreiburger/wetland_coverm_pilot/all_mags

# Annotation databases (all resolved as <launch_dir>/databases/<db> -> $D/databases).
# chenry's KBDL default set is kofam dbcan merops vog; the rest are on disk and wired in beta17.
DRAM2_DBS=${DRAM2_DBS:-"kofam dbcan merops vog camper canthyd methyl sulfur fegenie"}
# pfam: on disk but hard-disabled in beta17 (utils_nfcore_dram_pipeline/main.nf:89-90 errors out)
DISTILL_TOPIC=${DISTILL_TOPIC:-"default,camper"}
DISTILL_ECOSYSTEM=${DISTILL_ECOSYSTEM:-"eng_sys,ag"}

export NXF_VER=24.10.5
export NXF_HOME=$B/nxf_home                            # own home: 24.10.5 jar + nf-schema/nf-validation plugins
export NXF_JAVA_HOME=$D/env/env_nf                     # chenry's JDK 17 (read-only), as in his runs
export NXF_OPTS='-Xms1g -Xmx12g'
export NXF_CONDA_CACHEDIR=$B/conda_cache
export NXF_ANSI_LOG=false
export MAMBA_ROOT_PREFIX=$B/mamba_root
export TMPDIR=$WORK_ROOT/tmp
export NXF_TEMP=$TMPDIR
export PATH=$B/bin:/scratch/afreiburger/metag-hybrid/env/miniforge/condabin:$PATH
# ^ nextflow-classic launcher + micromamba (copies); the user's conda, which the pipeline's
#   startup check calls as `conda config --show channels` (conda-forge, bioconda)
