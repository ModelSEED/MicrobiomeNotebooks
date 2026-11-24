#!/bin/bash
# Validation of the re-pointed scripts before the production launch (kofam_rollout).
#  stage genes : md5 of the 8 first-pass control bins (kofam/meta/controls_chunk01.txt) vs kofam/mags.md5, then
#                call_genes.py (this track's copy) on them into validation/genes (NOT the production genes/ dir).
#  stage search: make_chunks.py (n_chunks=1) on the 8 controls -> validation/chunks/chunk_01; hmmsearch with the
#                68 marker profiles (db/markers.hmm) and -Z = the production total (935 MAGs, chunks/chunk_stats.txt);
#                postprocess.py chunk validation 01. Same commands as the kofam track's test/markers.
#  stage prodparse [chunk=01]: the 68 markers x one production chunk (rollout names) -> validation/prod_chunkXX_markers,
#                to check the binner/unit parsing of rollout names in postprocess.py (chunk 01) and to dry-run the Sec
#                rescue build on rollout names (chunk 21: Bacillota with grd genes).
# The controls are not part of the production run, its genes/ dir or its -Z.
set -euo pipefail
R=/scratch1/afreiburger/salterns_2026/kofam_rollout
S=/scratch1/afreiburger/salterns_2026
ANNOT=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_
export PRODIGAL=$ANNOT/bin/prodigal
PY=$ANNOT/bin/python3
H=$ANNOT/bin/hmmsearch
V=$R/validation
cd $V
case $1 in
genes)
    mkdir -p genes logs
    for d in db meta; do ln -sfn ../$d $d; done
    grep -F -f <(sed 's/$/.fna/' $R/meta/controls_validation.txt) $S/kofam/mags.md5 | sed "s#  #  $S/mags/#" > controls.md5
    [ "$(wc -l < controls.md5)" -eq 8 ]
    md5sum -c controls.md5 > logs/controls_md5_check.txt
    { echo -e "name"; cat $R/meta/controls_validation.txt; } > manifest.tsv
    sed 's#^#'$S'/mags/#; s#$#.fna#' $R/meta/controls_validation.txt \
      | xargs -P 8 -I{} nice -n 10 $PY $R/scripts/call_genes.py {} genes 1000
    [ "$(ls genes/*.done | wc -l)" -eq 8 ] && touch genes/ALL_DONE
    ;;
search)
    mkdir -p chunks hmm results logs
    : > logs/empty_controls.txt
    $PY $R/scripts/make_chunks.py logs/empty_controls.txt manifest.tsv genes chunks 1 > logs/make_chunks.log
    Z=$(awk '$1=="proteins"{print $2}' $R/chunks/chunk_stats.txt); T=14.7
    echo "Z=$Z (production, 935 MAGs)" > logs/Z.txt
    /usr/bin/time -v nice -n 10 $H --cpu 1 --noali -T $T --domT $T -Z $Z --tblout hmm/chunk_01.tblout \
        --domtblout hmm/chunk_01.domtblout -o hmm/chunk_01.out db/markers.hmm chunks/chunk_01.faa > logs/time.log 2>&1
    nice -n 10 $PY $R/scripts/postprocess.py chunk $V 01 > logs/postprocess.log 2>&1
    touch SEARCH_DONE
    ;;
prodparse)
    PC=${2:-01}                       # production chunk to test (01 = name parsing; 21 = Bacillota, Sec dry run)
    P=$V/prod_chunk${PC}_markers
    mkdir -p $P/hmm $P/results
    for d in db meta genes chunks; do ln -sfn ../../$d $P/$d; done
    cd $P
    Z=$(awk '$1=="proteins"{print $2}' chunks/chunk_stats.txt); T=14.7
    /usr/bin/time -v nice -n 10 $H --cpu 1 --noali -T $T --domT $T -Z $Z --tblout hmm/chunk_$PC.tblout \
        --domtblout hmm/chunk_$PC.domtblout -o hmm/chunk_$PC.out db/markers.hmm chunks/chunk_$PC.faa > time.log 2>&1
    nice -n 10 $PY $R/scripts/postprocess.py chunk $P $PC > postprocess.log 2>&1
    nice -n 10 $PY $R/scripts/postprocess.py combine $P >> postprocess.log 2>&1
    touch PRODPARSE_DONE
    ;;
esac
