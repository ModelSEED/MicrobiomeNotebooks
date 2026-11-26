#!/usr/bin/env bash
# Step 9 (v2: v1 script, v2 paths; reduced reference = contigs with a hit to any of the 286 v2 profiles): where do the classified reads sit in the sample's own assembly?
# The rollout deleted its BAMs after the depth tables (only R1_A_D1 and SF2_C_D1 keep theirs, in binning_test/), so the
# reads that 08_classify.py assigned to a target / marker / reported class are re-mapped, with the rollout's tools
# (bowtie2 2.5.5 end-to-end --sensitive, samtools 1.25), to the contigs of their own sample that carry at least one
# panel-HMM hit (E <= 1e-5; every gene of every panel family, incl. below-threshold homologs). A read counts as
# assembled if it aligns with >= 95 % identity (1 - NM / aligned read length), as CoverM's cut for the MAG abundances.
# Usage: bash 09_map_assigned.sh <sample>     -> map/<s>.hits.tsv.gz (read, mate, contig, pos, mapq, identity)
set -euo pipefail
s=$1
W=/scratch1/afreiburger/salterns_2026/read_markers_v2
V1=/scratch1/afreiburger/salterns_2026/read_markers
R=/scratch1/afreiburger/salterns_2026/rollout
B=/scratch1/afreiburger/salterns_2026/envs/bt_map/bin
A=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin
export TMPDIR=$W/tmp
T=${T:-8}
cd $W; mkdir -p map/$s
[ -f map/$s.done ] && exit 0
d=map/$s
# reduced reference
cat $V1/hmm/$s.tbl hmm_new/$s.0[0-3].tbl | awk '!/^#/{print $1}' | sed "s/^$s|//; s/_[0-9]*$//" | sort -u > $d/contigs.ids
$A/seqkit grep -f $d/contigs.ids $R/$s/asm/contigs.fa > $d/ref.fa
nice -n 10 $B/bowtie2-build --threads $T --seed 1 $d/ref.fa $d/ref > $d/build.log 2>&1
for m in 1 2; do
  zcat results/per_sample/$s.reads.tsv.gz | awk -F'\t' -v m=$m 'NR>1 && $2==m{print $1}' | sort -u > $d/ids_$m.txt
  $A/seqkit grep -f $d/ids_$m.txt $R/reads/$s/map_$m.fq.gz -o $d/reads_$m.fq.gz
  nice -n 10 $B/bowtie2 -x $d/ref -U $d/reads_$m.fq.gz --end-to-end --sensitive -p $T --no-unal 2> $d/bowtie2_$m.log | \
    $B/samtools view -F 0x904 - | \
    awk -F'\t' -v m=$m 'BEGIN{OFS="\t"} {nm=-1; for(i=12;i<=NF;i++) if($i ~ /^NM:i:/){nm=substr($i,6)}
         q=0; c=$6; while(match(c,/[0-9]+[MIS=X]/)){n=substr(c,RSTART,RLENGTH-1); t=substr(c,RSTART+RLENGTH-1,1); if(t!="S") q+=n; c=substr(c,RSTART+RLENGTH)}
         print $1, m, $3, $4, $5, (q>0 ? 1-nm/q : 0)}' > $d/hits_$m.tsv
done
cat $d/hits_1.tsv $d/hits_2.tsv | gzip > map/$s.hits.tsv.gz
echo -e "$s\t$(wc -l < $d/contigs.ids)\t$(cat $d/ids_1.txt $d/ids_2.txt | wc -l)\t$(cat $d/hits_1.tsv $d/hits_2.tsv | wc -l)" > map/$s.summary
rm -rf $d/ref.* $d/reads_*.fq.gz $d/hits_*.tsv
touch map/$s.done
