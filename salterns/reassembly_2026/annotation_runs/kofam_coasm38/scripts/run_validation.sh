#!/bin/bash
# Validation control chunk for kofam_coasm38: 2 rollout MAGs that kofam_rollout already annotated
# (meta/controls_validation.txt: R2_A_H2O__metabat2.bin.034 Natronomonas, phnGHIJKLM 7/7 on one contig;
#  R1_C_D2__maxbin2.Bin.002_sub Bacillota M55B127, grdH by Sec read-through + grdI, MaxBin2).
#  stage genes : md5 of the 2 FASTAs (rollout/mags_kept/<mag>.fa) vs kofam_rollout/mags.md5; call_genes.py (this
#                track's copy) into validation/genes (NOT the production genes/ dir); md5 of the .faa/.gff/.genes.tsv/
#                .done against kofam_rollout/genes/ -> logs/genes_md5_compare.txt (must be identical).
#  stage search: make_chunks.py (1 chunk, one group per MAG) -> validation/chunks/chunk_01; hmmsearch with ALL 27,981
#                profiles (db/kofam_all.hmm) at -Z = this track's production total (the 38 co-assembly MAGs,
#                chunks/chunk_stats.txt), same flags as run_chunk.sh; postprocess.py chunk + combine in validation/.
# kofam_rollout's run_validation.sh searched only the 68 marker profiles (its reference was a marker-only test run);
# here the reference is kofam_rollout's full-DB production output, so the full DB is searched.
# The controls are not part of the production run, its genes/ dir, its chunks or its -Z, and not in results/.
set -euo pipefail
R=/scratch1/afreiburger/salterns_2026/kofam_coasm38
S=/scratch1/afreiburger/salterns_2026
RO=$S/kofam_rollout
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
    grep -F -f <(sed 's/$/.fa/' $R/meta/controls_validation.txt) $RO/mags.md5 | sed "s#  #  $S/rollout/mags_kept/#" > controls.md5
    [ "$(wc -l < controls.md5)" -eq 2 ]
    md5sum -c controls.md5 > logs/controls_md5_check.txt
    { echo -e "name"; cat $R/meta/controls_validation.txt; } > manifest.tsv
    sed 's#^#'$S'/rollout/mags_kept/#; s#$#.fa#' $R/meta/controls_validation.txt \
      | xargs -P 2 -I{} nice -n 10 $PY $R/scripts/call_genes.py {} genes 1000
    [ "$(ls genes/*.done | wc -l)" -eq 2 ] && touch genes/ALL_DONE
    : > logs/genes_md5_compare.txt
    for m in $(cat $R/meta/controls_validation.txt); do
        for suf in faa gff genes.tsv done; do
            a=$(md5sum < genes/$m.$suf | cut -d' ' -f1); b=$(md5sum < $RO/genes/$m.$suf | cut -d' ' -f1)
            echo -e "$m\t$suf\t$a\t$b\t$([ "$a" = "$b" ] && echo identical || echo DIFFERENT)" >> logs/genes_md5_compare.txt
        done
    done
    cat logs/genes_md5_compare.txt
    ;;
search)
    mkdir -p chunks hmm results logs
    : > logs/empty_controls.txt
    $PY $R/scripts/make_chunks.py logs/empty_controls.txt manifest.tsv genes chunks 1 > logs/make_chunks.log
    Z=$(awk '$1=="proteins"{print $2}' $R/chunks/chunk_stats.txt); T=14.7
    : "${Z:?no production chunk_stats.txt}"
    echo -e "Z_new\t$Z\nZ_old\t$(awk '$1=="proteins"{print $2}' $RO/chunks/chunk_stats.txt)" > logs/Z.txt
    p=hmm/chunk_01
    echo "$(date -u +%FT%TZ) start validation chunk host $(hostname) pid $$ Z=$Z T=$T hmm=db/kofam_all.hmm" > logs/search.log
    /usr/bin/time -v nice -n 10 $H --cpu 1 --noali -T $T --domT $T -Z $Z --tblout $p.tblout \
        --domtblout $p.domtblout -o $p.out db/kofam_all.hmm chunks/chunk_01.faa >> logs/search.log 2>&1
    [ "$(tail -1 $p.out)" = "[ok]" ] || { echo "$(date -u +%FT%TZ) hmmsearch FAILED" >> logs/search.log; touch $p.failed; exit 1; }
    touch $p.done
    echo "$(date -u +%FT%TZ) hmmsearch done" >> logs/search.log
    nice -n 10 $PY $R/scripts/postprocess.py chunk $V 01 > logs/postprocess.log 2>&1
    nice -n 10 $PY $R/scripts/postprocess.py combine $V >> logs/postprocess.log 2>&1
    nice -n 10 gzip -f $p.tblout $p.domtblout $p.out
    echo "$(date -u +%FT%TZ) end" >> logs/search.log
    touch SEARCH_DONE
    ;;
esac
