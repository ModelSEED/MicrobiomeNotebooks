#!/usr/bin/env bash
# Threshold test: read-based abundance per contig class. MEGAHIT's `multi=` is k-mer coverage at the k where a contig
# was finished, so short contigs finished at small k carry inflated values next to k141 contigs; reads are the
# unbiased measure, and they also give the right denominator (all of the sample's reads, assembled or not).
# For each re-assembled sample: minimap2 2.28 -x sr (paired, primary alignments only) of all reads against the 1 kb
# re-assembly; per contig, the number of mapped reads (mates counted separately). The KBase assembly is the
# re-assembly's contigs >= 2 kb (compare_reassembly.py: identical), so one mapping serves both floors.
# Output: <run>/<run>.contig_reads.tsv (contig, reads) and <run>/<run>.reads_total (reads in the fastq files).
# Usage (TF): cd ~/tf_runs/wetland_reasm && setsid nohup bash map_reads_tf.sh SRR... > map.log 2>&1 < /dev/null &
set -uo pipefail
B=$HOME/tf_runs/wetland_reasm; MM2=$B/minimap2; THREADS=${THREADS:-8}
for r in "$@"; do
  d=$B/$r; cd "$d" || continue
  until [ -e DONE ] || grep -q "FAIL $r" "$B/events.log" 2>/dev/null; do sleep 300; done
  [ -e DONE ] || continue
  [ -s "$r.contig_reads.tsv" ] && continue
  t0=$(date +%s)
  "$MM2" -x sr --secondary=no -t "$THREADS" "$r.min1000.fa" "${r}_1.fastq.gz" "${r}_2.fastq.gz" 2> "$r.map.log" |
    awk -F'\t' '$13 == "tp:A:P" {c[$6]++} END {for (k in c) print k "\t" c[k]}' > "$r.contig_reads.tsv.tmp"
  # (PAF tag tp:A:P = primary; each mate is its own record. Supplementary records, rare for short reads, also count.)
  n=$(( $(zcat "${r}_1.fastq.gz" | wc -l) / 4 * 2 ))
  echo "$n" > "$r.reads_total"; mv "$r.contig_reads.tsv.tmp" "$r.contig_reads.tsv"
  echo "[$(date '+%F %T')] MAPPED $r $(( $(date +%s) - t0 ))s reads=$n" >> "$B/events.log"
done
