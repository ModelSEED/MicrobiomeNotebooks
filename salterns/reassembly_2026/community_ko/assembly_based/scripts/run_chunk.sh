#!/bin/bash
# One chunk: hmmsearch of all 27,981 KOfam 2026-10-26 profiles x one chunk of unique proteins, single thread
# (--cpu 0 = serial, 1 thread), nice 10; then the chunk's postprocessing and gzip of the tblout.
# Flags = kofam_rollout/scripts/run_chunk.sh (--noali -T 14.7 --domT 14.7, fixed -Z) except: --cpu 0 instead of
# --cpu 1 (same results, one thread instead of worker + reader), no --domtblout (only used there for HMM coverage),
# -o /dev/null (completion is checked on the tblout's '# [ok]' tail line).
# usage: run_chunk.sh <chunk id 0001..>
set -uo pipefail
W=/scratch1/afreiburger/salterns_2026/community_ko
A=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin
cid=$1
cd $W
export TMPDIR=$W/tmp
Z=$(awk '$1=="proteins"{print $2}' chunks/chunk_stats.txt)   # all 12,055,118 genes of the 39 assemblies
T=14.7
p=hmm/chunk_$cid
L=logs/chunks/chunk_$cid.log
if [ ! -e $p.done ]; then
    rm -f $p.failed $p.tblout $p.tblout.gz
    echo "$(date -u +%FT%TZ) start chunk $cid pid $$ Z=$Z T=$T" >> $L
    /usr/bin/time -v nice -n 10 $A/hmmsearch --cpu 0 --noali -T $T --domT $T -Z $Z \
        --tblout $p.tblout -o /dev/null db/kofam_all.hmm chunks/chunk_$cid.faa >> $L 2>&1
    rc=$?
    if [ $rc -ne 0 ] || [ "$(tail -1 $p.tblout)" != "# [ok]" ]; then
        echo "$(date -u +%FT%TZ) hmmsearch FAILED rc=$rc" >> $L; touch $p.failed; exit 1
    fi
    touch $p.done
    echo "$(date -u +%FT%TZ) hmmsearch done" >> $L
fi
if [ ! -e results/chunks/chunk_$cid.postprocessed ]; then
    nice -n 10 $A/python3 -I scripts/02_postprocess_chunk.py $W $cid >> $L 2>&1 \
        || { echo "$(date -u +%FT%TZ) postprocess FAILED" >> $L; touch results/chunks/chunk_$cid.postprocess_failed; exit 1; }
fi
[ -e $p.tblout ] && nice -n 10 gzip -f $p.tblout
echo "$(date -u +%FT%TZ) end chunk $cid" >> $L
