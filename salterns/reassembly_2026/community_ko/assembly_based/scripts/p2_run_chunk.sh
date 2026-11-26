#!/bin/bash
# Two-pass design, pass 2 (confirm): all 27,981 KOfam profiles (db/kofam_all.hmm) x one chunk of pass-1 candidates
# (pass2/chunks/chunk_<id>.faa), exactly the flags and -Z (12,055,118) of the full search (run_chunk.sh), one thread,
# nice 10; then 02_postprocess_chunk.py (unchanged call rule) on the pass2 dir and gzip of the tblout.
# usage: p2_run_chunk.sh <chunk id>
set -uo pipefail
W=/scratch1/afreiburger/salterns_2026/community_ko
A=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin
cid=$1
cd $W
export TMPDIR=$W/tmp
Z=$(awk '$1=="proteins"{print $2}' chunks/chunk_stats.txt)   # all 12,055,118 genes, as in the full search
T=14.7
p=pass2/hmm/chunk_$cid
L=logs/pass2/chunk_$cid.log
if [ ! -e $p.done ]; then
    rm -f $p.failed $p.tblout $p.tblout.gz
    echo "$(date -u +%FT%TZ) start pass2 chunk $cid pid $$ Z=$Z T=$T" >> $L
    /usr/bin/time -v nice -n 10 $A/hmmsearch --cpu 0 --noali -T $T --domT $T -Z $Z \
        --tblout $p.tblout -o /dev/null db/kofam_all.hmm pass2/chunks/chunk_$cid.faa >> $L 2>&1
    rc=$?
    if [ $rc -ne 0 ] || [ "$(tail -1 $p.tblout)" != "# [ok]" ]; then
        echo "$(date -u +%FT%TZ) hmmsearch FAILED rc=$rc" >> $L; touch $p.failed; exit 1
    fi
    touch $p.done
    echo "$(date -u +%FT%TZ) hmmsearch done" >> $L
fi
if [ ! -e pass2/results/chunks/chunk_$cid.postprocessed ]; then
    nice -n 10 $A/python3 -I scripts/02_postprocess_chunk.py $W/pass2 $cid >> $L 2>&1 \
        || { echo "$(date -u +%FT%TZ) postprocess FAILED" >> $L; touch pass2/results/chunks/chunk_$cid.postprocess_failed; exit 1; }
fi
[ -e $p.tblout ] && nice -n 10 gzip -f $p.tblout
echo "$(date -u +%FT%TZ) end pass2 chunk $cid" >> $L
