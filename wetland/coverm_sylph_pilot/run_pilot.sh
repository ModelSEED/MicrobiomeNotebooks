#!/usr/bin/env bash
# CoverM vs sylph pilot on 2 wetland metagenomes (PRJNA869844: SRR21120499, SRR21120498), poplar.
# Reference = the 359 MAGs of PRJNA869844's 9 metagenomes, dereplicated at 95% ANI (dRep + skani, CheckM v1 quality).
# Both tools run on the SAME reference; sylph also against all wetland MAGs if sylph_all/ exists. Timed with /usr/bin/time -v.
set -euo pipefail
cd "$(dirname "$0")"
QC=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/8c703f8abdc13cdc1ef374a2bcf37f14_   # dRep, skani, CheckM2
BIN=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/8c3854a8c5c280d41ee307947faefca0_  # CoverM 0.7.0
MM=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/d501a828539edb92b78db288661dd07f_   # minimap2 2.30
SY=/scratch1/afreiburger/tools/sylph_env                                                        # sylph 1.0.0
T=${T:-32}
SAMPLES="SRR21120499 SRR21120498"
ev() { echo "[$(date '+%F %T')] $*" | tee -a pilot.events; }
tm() { local tag=$1; shift; /usr/bin/time -v -o "time_${tag}.txt" "$@"; \
       ev "TIME $tag wall=$(grep 'Elapsed (wall' time_${tag}.txt | awk '{print $NF}') maxrss_GB=$(awk '/Maximum resident/ {printf "%.1f", $NF/1048576}' time_${tag}.txt) cpu%=$(grep 'Percent of CPU' time_${tag}.txt | awk '{print $NF}')"; }
ev "START (load $(cut -d' ' -f1-3 /proc/loadavg), threads $T)"
for s in $SAMPLES; do [ -f reads/${s}_reads.inter.fastq.gz.done ] || { ev "missing reads for $s"; exit 1; }; done

if [ ! -d derep/dereplicated_genomes ]; then
  tm drep env PATH=$QC/bin:$PATH dRep dereplicate derep -g mags/*.fna --genomeInfo genomeInfo.csv \
      -comp 50 -con 10 -sa 0.95 --S_algorithm skani -p $T
fi
ev "DEREP $(ls derep/dereplicated_genomes | wc -l) of $(ls mags | wc -l) MAGs kept (95% ANI)"

for s in $SAMPLES; do
  [ -s coverm_${s}.tsv ] || tm coverm_${s} env PATH=$BIN/bin:$MM/bin:$PATH coverm genome \
      --interleaved reads/${s}_reads.inter.fastq.gz --genome-fasta-directory derep/dereplicated_genomes -x fna \
      --mapper minimap2-sr -m relative_abundance mean trimmed_mean covered_fraction \
      --bam-file-cache-directory bams_${s} -t $T -o coverm_${s}.tsv
  [ -s coverm_${s}.raw.tsv ] || tm coverm_raw_${s} env PATH=$BIN/bin:$MM/bin:$PATH coverm genome \
      --bam-files bams_${s}/*.bam --separator '~' --min-covered-fraction 0 \
      -m count reads_per_base covered_fraction relative_abundance -t $T -o coverm_${s}.raw.tsv
done

[ -s sylph_derep.syldb ] || tm sylph_db $SY/bin/sylph sketch -g derep/dereplicated_genomes/*.fna -o sylph_derep -t $T
mkdir -p sylph_reads
for s in $SAMPLES; do
  [ -s sylph_reads/${s}_reads.inter.fastq.gz.sylsp ] || tm sylph_sketch_${s} $SY/bin/sylph sketch \
      -r reads/${s}_reads.inter.fastq.gz -d sylph_reads -t $T
done
tm sylph_profile_derep sh -c "$SY/bin/sylph profile sylph_derep.syldb sylph_reads/*.sylsp -t $T -u -o sylph_derep.tsv"
if [ -d all_mags ]; then
  [ -s sylph_all.syldb ] || tm sylph_db_all sh -c "ls all_mags/*.fna > all_mags.list && $SY/bin/sylph sketch --gl all_mags.list -o sylph_all -t $T"
  tm sylph_profile_all sh -c "$SY/bin/sylph profile sylph_all.syldb sylph_reads/*.sylsp -t $T -u -o sylph_all.tsv"
fi
ev "DONE"
