#!/usr/bin/env bash
# Step 1: Prodigal 2.6.3 (-p meta -g 11 -q; the kofam_rollout flags) on ALL contigs (>= 1 kb, asm/contigs.fa)
# of each of the 39 per-sample rollout assemblies. Each assembly is split into 8 pieces so 40 single-CPU jobs run at once.
# Outputs per sample in genes/<s>/: <s>.faa, <s>.fna, <s>.gff (gene IDs "<sample>|<contig>_<n>"), DONE.
# Run detached:  nohup setsid bash scripts/01_prodigal.sh > logs/01_prodigal.nohup 2>&1 < /dev/null &
set -euo pipefail
W=/scratch1/afreiburger/salterns_2026/read_markers
R=/scratch1/afreiburger/salterns_2026/rollout
A=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin
export TMPDIR=$W/tmp
NP=${NP:-40}
cd $W
mkdir -p genes logs
SAMPLES=$(ls $R/reads)
# 1. split
for s in $SAMPLES; do
  d=genes/$s; [ -f $d/DONE ] && continue
  mkdir -p $d
  if [ ! -f $d/SPLIT_OK ]; then
    rm -rf $d/pieces
    nice -n 10 $A/seqkit split2 -p 8 -O $d/pieces $R/$s/asm/contigs.fa 2> $d/split.log
    touch $d/SPLIT_OK
  fi
done
# 2. prodigal on every piece (skip pieces already done)
find genes -path "*/pieces/*.fa" | sort | while read f; do [ -f ${f%.fa}.ok ] || echo $f; done | \
  xargs -P $NP -I{} bash -c 'f={}; b=${f%.fa}; nice -n 10 '$A'/prodigal -p meta -g 11 -q -i $f -a $b.faa -d $b.fna -f gff -o $b.gff && touch $b.ok'
# 3. merge per sample with sample-prefixed IDs; check counts
for s in $SAMPLES; do
  d=genes/$s; [ -f $d/DONE ] && continue
  n_fa=$(ls $d/pieces/*.fa | wc -l); n_ok=$(ls $d/pieces/*.ok | wc -l)
  [ "$n_fa" = "$n_ok" ] || { echo "$s: $n_ok/$n_fa pieces ok" >&2; continue; }
  cat $d/pieces/*.faa | sed "s/^>/>$s|/" > $d/$s.faa
  cat $d/pieces/*.fna | sed "s/^>/>$s|/" > $d/$s.fna
  cat $d/pieces/*.gff | grep -v '^#' > $d/$s.gff
  nc_in=$(grep -c '^>' $R/$s/asm/contigs.fa); ng=$(grep -c '^>' $d/$s.faa)
  echo -e "$s\tcontigs_in\t$nc_in\tgenes\t$ng" > $d/summary.txt
  rm -f $d/pieces/*.faa $d/pieces/*.fna $d/pieces/*.gff $d/pieces/*.fa
  touch $d/DONE
done
ls genes/*/DONE | wc -l > genes/N_DONE   # (counts genes/coasm520/DONE too when 01b ran first)
echo "$(date -u +%FT%TZ) prodigal done: $(cat genes/N_DONE) samples"
