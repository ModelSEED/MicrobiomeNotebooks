#!/usr/bin/env bash
# Step 1b: Prodigal 2.6.3 (-p meta -g 11 -q) on the 520 kept co-assembly MAGs (coasm_qc/mags_kept/*.fa), so the read
# survey can attribute reads to every MAG of the final 585-species set (the 39 per-sample assemblies hold only the 935
# rollout MAGs). Output: genes/coasm520/coasm520.faa/.fna/.gff, gene IDs "coasm|<mag_id>|<contig>_<n>".
set -euo pipefail
W=/scratch1/afreiburger/salterns_2026/read_markers
M=/scratch1/afreiburger/salterns_2026/coasm_qc/mags_kept
A=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin
export TMPDIR=$W/tmp
NP=${NP:-24}
cd $W
d=genes/coasm520; mkdir -p $d/per_mag
ls $M/*.fa | while read f; do b=$(basename $f .fa); [ -f $d/per_mag/$b.ok ] || echo $f; done | \
  xargs -P $NP -I{} bash -c 'f={}; b=$(basename $f .fa); o='$d'/per_mag/$b; nice -n 10 '$A'/prodigal -p meta -g 11 -q -i $f -a $o.faa -d $o.fna -f gff -o $o.gff && touch $o.ok'
n=$(ls $d/per_mag/*.ok | wc -l); [ "$n" = 520 ] || { echo "only $n/520 MAGs done" >&2; exit 1; }
: > $d/coasm520.faa; : > $d/coasm520.fna; : > $d/coasm520.gff
for f in $d/per_mag/*.faa; do b=$(basename $f .faa)
  sed "s/^>/>coasm|$b|/" $f >> $d/coasm520.faa
  sed "s/^>/>coasm|$b|/" $d/per_mag/$b.fna >> $d/coasm520.fna
  grep -v '^#' $d/per_mag/$b.gff | sed "s/^/$b\t/" >> $d/coasm520.gff
done
rm -f $d/per_mag/*.faa $d/per_mag/*.fna $d/per_mag/*.gff
touch $d/DONE
echo "$(date -u +%FT%TZ) coasm520 prodigal done: $(grep -c '^>' $d/coasm520.faa) genes"
