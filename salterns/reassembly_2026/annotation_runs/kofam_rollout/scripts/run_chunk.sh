#!/bin/bash
# One hmmsearch job: all KOfam profiles (2026-10-26 release) x one protein chunk, single CPU, nice 10,
# then the chunk's postprocessing (KO calls, marker/guild tables) and gzip of the raw tables.
# kofam_rollout copy of kofam/scripts/run_chunk.sh: only W changed; -Z is read from this track's chunks/chunk_stats.txt
# (= total proteins of the 935 rollout MAGs; the validation controls are not in it).
# usage: run_chunk.sh <chunk id 01..40> [hmm=db/kofam_all.hmm] [outdir=hmm]
set -uo pipefail
W=/scratch1/afreiburger/salterns_2026/kofam_rollout
ANNOT=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_
cid=$1; HMM=${2:-db/kofam_all.hmm}; OUT=${3:-hmm}
cd $W
Z=$(awk '$1=="proteins"{print $2}' chunks/chunk_stats.txt)   # total proteins of all 935 MAGs -> chunk-independent E-values
T=14.7    # lowest numeric threshold in ko_list 2026-10-26: report every hit that could pass any KO (rule is score > threshold)
p=$OUT/chunk_$cid
[ -e $p.done ] && { echo "chunk $cid already done"; exit 0; }
rm -f $p.failed
echo "$(date -u +%FT%TZ) start chunk $cid host $(hostname) pid $$ Z=$Z T=$T hmm=$HMM" >> logs/chunk_$cid.log
/usr/bin/time -v nice -n 10 $ANNOT/bin/hmmsearch --cpu 1 --noali -T $T --domT $T -Z $Z \
    --tblout $p.tblout --domtblout $p.domtblout -o $p.out $HMM chunks/chunk_$cid.faa >> logs/chunk_$cid.log 2>&1
rc=$?
if [ $rc -ne 0 ] || [ "$(tail -1 $p.out)" != "[ok]" ]; then
    echo "$(date -u +%FT%TZ) hmmsearch FAILED rc=$rc" >> logs/chunk_$cid.log; touch $p.failed; exit 1
fi
touch $p.done
echo "$(date -u +%FT%TZ) hmmsearch done" >> logs/chunk_$cid.log
if [ "$OUT" = hmm ]; then
    nice -n 10 $ANNOT/bin/python3 scripts/postprocess.py chunk $W $cid >> logs/chunk_$cid.log 2>&1 \
        || { echo "postprocess FAILED" >> logs/chunk_$cid.log; touch results/chunks/chunk_$cid.postprocess_failed; }
fi
nice -n 10 gzip -f $p.tblout $p.domtblout $p.out
echo "$(date -u +%FT%TZ) end chunk $cid" >> logs/chunk_$cid.log
