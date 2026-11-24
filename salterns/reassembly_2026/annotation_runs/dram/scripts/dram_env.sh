# Sourced by every DRAM step on poplar.
export E=/scratch1/afreiburger/salterns_2026/envs/dram15
export PATH=$E/bin:$PATH
export PYTHONNOUSERSITE=1            # ~/.local/lib/python3.10 otherwise leaks an old requests/urllib3 into the env
export PYTHONWARNINGS=ignore::UserWarning
export W=/scratch1/afreiburger/salterns_2026/dram
export DRAMDB=/scratch1/afreiburger/salterns_2026/dram_db
export RAW=$DRAMDB/raw
export MAGS=$W/mags DRAM_OUT=$W/dram_out DRAM_LOGS=$W/logs
export TAX=$W/inputs/gtdb_taxonomy_558.tsv QUAL=$W/inputs/checkm_quality_558.tsv
export SHEET=$W/inputs/custom_distillate_salterns.dram1.tsv
export MAG_LIST=$W/inputs/best_per_species_558.tsv
export TEST_MAG=concoct.bin.101.fasta__R1_sediment   # best new Methanolobus (archaeon; mcr, Pyl mtt/mtb/mtm)
export THREADS=6 P=8 SETUP_THREADS=24                # <= 48 concurrent threads
export MIN_CONTIG=1000
export TMPDIR=$W/tmp
mkdir -p "$TMPDIR" "$DRAM_OUT" "$DRAM_LOGS/mags"
