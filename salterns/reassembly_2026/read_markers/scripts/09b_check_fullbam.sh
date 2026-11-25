#!/usr/bin/env bash
# Step 9b: check of the reduced-reference re-mapping (09_map_assigned.sh) against the full self-mapping BAMs that the
# binning test kept for R1_A_D1 and SF2_C_D1 (bowtie2 paired-end, all reads vs all contigs >= 1 kb; the contigs are
# byte-identical to rollout/<s>/asm/contigs.fa per README_rollout). For every classified read: contig in the full
# BAM (mate-specific, >= 95 % identity) vs contig in the reduced-reference mapping.
# Usage: bash 09b_check_fullbam.sh <sample>   -> map/<s>.fullbam_check.tsv
set -euo pipefail
s=$1
W=/scratch1/afreiburger/salterns_2026/read_markers
B=/scratch1/afreiburger/salterns_2026/envs/bt_map/bin
BAM=/scratch1/afreiburger/salterns_2026/binning_test/$s/bam/$s.bam
cd $W
d=map/$s; mkdir -p $d
cat $d/ids_1.txt $d/ids_2.txt | sort -u > $d/ids_all.txt
nice -n 10 $B/samtools view -@ 4 -F 0x904 -N $d/ids_all.txt $BAM | \
  awk -F'\t' 'BEGIN{OFS="\t"} {m=(and($2,64)?1:2); nm=-1; for(i=12;i<=NF;i++) if($i ~ /^NM:i:/){nm=substr($i,6)}
       q=0; c=$6; while(match(c,/[0-9]+[MIS=X]/)){n=substr(c,RSTART,RLENGTH-1); t=substr(c,RSTART+RLENGTH-1,1); if(t!="S") q+=n; c=substr(c,RSTART+RLENGTH)}
       print $1, m, $3, $4, $5, (q>0 ? 1-nm/q : 0)}' | gzip > map/$s.fullbam_hits.tsv.gz
echo "$(date -u +%FT%TZ) $s full-BAM lookup done"
