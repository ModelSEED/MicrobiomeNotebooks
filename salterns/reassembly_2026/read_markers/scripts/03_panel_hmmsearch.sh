#!/usr/bin/env bash
# Step 3: hmmsearch (HMMER 3.4) of the KOfam 2026-10-26 profiles of the panel KOs (meta/ko_panel.tsv: targets, decoys,
# markers) against (i) all proteins of each per-sample assembly, (ii) the 520 co-assembly MAGs, (iii) the external
# references (KEGG GENES + UniProt families). Permissive reporting (-E 1e-5 --domE 1e-5) so that below-threshold homologs
# are visible too; KOfam thresholds are applied afterwards (04_build_db.py). Fixed -Z 20000000 for every input so
# E-values are comparable across inputs. One CPU per job, NP jobs at once, nice 10.
# Usage: bash 03_panel_hmmsearch.sh <input.faa> [<input.faa> ...]   (outputs hmm/<basename>.{tbl,domtbl,done})
set -euo pipefail
W=/scratch1/afreiburger/salterns_2026/read_markers
A=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin
P=/scratch1/afreiburger/salterns_2026/kofam_db/profiles
export TMPDIR=$W/tmp
NP=${NP:-16}
cd $W
mkdir -p hmm db
if [ ! -s db/panel.hmm ]; then
  awk -F'\t' 'NR>1{print $1}' meta/ko_panel.tsv | sort -u | while read k; do cat $P/$k.hmm; done > db/panel.hmm.tmp
  n=$(grep -c '^NAME' db/panel.hmm.tmp); m=$(awk -F'\t' 'NR>1{print $1}' meta/ko_panel.tsv | sort -u | wc -l)
  [ "$n" = "$m" ] || { echo "panel.hmm has $n profiles, panel has $m KOs" >&2; exit 1; }
  mv db/panel.hmm.tmp db/panel.hmm
fi
for f in "$@"; do echo $f; done | while read f; do b=$(basename $f .faa); [ -f hmm/$b.done ] || echo $f; done | \
  xargs -P $NP -I{} bash -c 'f={}; b=$(basename $f .faa); nice -n 10 '$A'/hmmsearch --cpu 1 --noali -E 1e-5 --domE 1e-5 -Z 20000000 \
     --tblout hmm/$b.tbl --domtblout hmm/$b.domtbl -o /dev/null db/panel.hmm $f && touch hmm/$b.done'
echo "$(date -u +%FT%TZ) panel hmmsearch done for $# inputs"
