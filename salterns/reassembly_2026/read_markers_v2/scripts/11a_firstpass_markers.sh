#!/usr/bin/env bash
# Step 11a (v2, deliverable 1): extra GTDB-classified taxonomy references for the ribosomal-protein marker families.
# The v1 DB knows the GTDB genus only of rollout / co-assembly MAG members; 25-58 % of the marker reads hit salterns-only
# unbinned clusters (genus_source 'none' in v1 refs/rep_info.tsv.gz). The 1,324 first_pass MAGs (salterns_2026/mags/,
# GTDB R232 classification in meta/mags_three_way.tsv) hold genera that are in no rollout/co-assembly MAG (e.g. the only
# Marinilabilia MAG). Here: Prodigal 2.6.3 (-p meta -g 11 -q, as v1) on every first_pass MAG, then hmmsearch (HMMER 3.4,
# v1 flags) of the 10 marker KOfam profiles (the 8 families; uL4/uL16 are KO pairs) against them.
# Outputs: genes_firstpass/firstpass.faa (IDs "fp|<genome>|<gene>"), hmm/firstpass_markers.{tbl,domtbl}
set -euo pipefail
W=/scratch1/afreiburger/salterns_2026/read_markers_v2
V1=/scratch1/afreiburger/salterns_2026/read_markers
M=/scratch1/afreiburger/salterns_2026/mags
A=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin
P=/scratch1/afreiburger/salterns_2026/kofam_db/profiles
export TMPDIR=$W/tmp
NP=${NP:-24}
cd $W; mkdir -p genes_firstpass/per_mag hmm db
ls $M/*.fna | while read f; do b=$(basename $f .fna); [ -f genes_firstpass/per_mag/$b.ok ] || echo $f; done | \
  xargs -P $NP -I{} bash -c 'f={}; b=$(basename $f .fna); o=genes_firstpass/per_mag/$b; nice -n 10 '$A'/prodigal -p meta -g 11 -q -i $f -a $o.faa -f gff -o /dev/null && touch $o.ok'
n=$(ls genes_firstpass/per_mag/*.ok | wc -l); echo "prodigal ok: $n"
: > genes_firstpass/firstpass.faa
for f in genes_firstpass/per_mag/*.faa; do b=$(basename $f .faa); sed "s/^>/>fp|$b|/" $f >> genes_firstpass/firstpass.faa; done
for k in K02886 K02867 K02982 K02948 K02950 K02996 K02926 K02930 K02878 K02866; do cat $P/$k.hmm; done > db/markers10.hmm
# split the protein file to run 24 single-CPU hmmsearch jobs; -Z fixed as in v1
$A/seqkit split2 -p $NP -O genes_firstpass/split genes_firstpass/firstpass.faa 2> /dev/null
ls genes_firstpass/split/*.faa | xargs -P $NP -I{} bash -c 'f={}; nice -n 10 '$A'/hmmsearch --cpu 1 --noali -E 1e-5 --domE 1e-5 -Z 20000000 --tblout $f.tbl --domtblout $f.domtbl -o /dev/null db/markers10.hmm $f'
cat genes_firstpass/split/*.tbl | grep -v '^#' > hmm/firstpass_markers.tbl
cat genes_firstpass/split/*.domtbl | grep -v '^#' > hmm/firstpass_markers.domtbl
rm -rf genes_firstpass/split genes_firstpass/per_mag/*.faa
echo "$(date) 11a done: $(grep -c '^>' genes_firstpass/firstpass.faa) proteins" > genes_firstpass/DONE
