#!/usr/bin/env bash
# Stage DRAM2 inputs as symlinks with DRAM2-safe names and write the name map.
#
# DRAM2 beta17 needs dot-free sample names: with --rename it renames <name>.fna -> <name>.fna but
# re-derives the sample name with dots replaced by '-', and combine_annotations.py recovers each
# gene's input_fasta as commonprefix(<faa stem>, <gene id>); a dotted name breaks that join.
# So  bin.014.fasta__SRR11577120  ->  bin-014-fasta__SRR11577120  (no '-' in the originals,
# so the mapping is 1:1 and reversible).
#
# Chunks: MAGs sorted by FASTA size (desc) and dealt round-robin into N_CHUNKS, so every chunk
# has the same size mix.  Also stages inputs/test3 (3 MAGs of different sizes).
set -euo pipefail
B=/scratch1/afreiburger/wetland_dram2
MAG_DIR=/scratch1/afreiburger/wetland_coverm_pilot/all_mags
N_CHUNKS=${N_CHUNKS:-5}
TEST3="bin.072.fasta__SRR9700688 bin.049.fasta__SRR8556254 bin.020.fasta__SRR11577034"

cd "$B"
[ -e inputs/name_map.tsv ] && { echo "inputs/name_map.tsv exists; remove inputs/ to restage"; exit 1; }
tail -n +2 dram2_hq_multisample.tsv | cut -f1 | while read -r m; do
  f=$MAG_DIR/$m.fna; [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }
  printf '%s\t%s\t%s\n' "$m" "${m//./-}" "$(stat -L -c %s "$f")"
done | sort -t$'\t' -k3,3nr -k1,1 \
  | awk -v n="$N_CHUNKS" 'BEGIN{OFS="\t"; print "mag","dram2_name","chunk","fna_bytes"}
                          {printf "%s\t%s\tchunk_%02d\t%s\n", $1, $2, ((NR-1)%n)+1, $3}' > inputs/name_map.tsv

n=$(tail -n +2 inputs/name_map.tsv | wc -l)
u=$(tail -n +2 inputs/name_map.tsv | cut -f2 | sort -u | wc -l)
[ "$n" = "$u" ] || { echo "sanitized names not unique ($u/$n)"; exit 1; }

tail -n +2 inputs/name_map.tsv | while IFS=$'\t' read -r m s c b; do
  mkdir -p "inputs/$c"; ln -sfn "$MAG_DIR/$m.fna" "inputs/$c/$s.fna"
done
mkdir -p inputs/test3
for m in $TEST3; do ln -sfn "$MAG_DIR/$m.fna" "inputs/test3/${m//./-}.fna"; done

echo "staged $n MAGs"; for d in inputs/*/; do echo "$d $(ls "$d" | wc -l)"; done
