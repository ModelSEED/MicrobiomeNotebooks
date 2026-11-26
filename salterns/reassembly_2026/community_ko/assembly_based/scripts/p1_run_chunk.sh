#!/bin/bash
# Two-pass design, pass 1 (screen): the 148 profiles of db/pass1_148.hmm (138 panel KOs of ko_panel_methylated.tsv +
# the 10 marker KOs of the 8 ribosomal-protein families; cut from kofam_db/profiles/, byte-identical to their ranges
# in kofam_all.hmm) x one original chunk (chunks/chunk_<id>.faa, unique proteins), same flags and -Z as the full
# search, plus --domtblout (HMM coverage for the five no-threshold panel KOs). One thread (--cpu 0), nice 10.
# usage: p1_run_chunk.sh <chunk id>
set -uo pipefail
W=/scratch1/afreiburger/salterns_2026/community_ko
A=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin
cid=$1
cd $W
export TMPDIR=$W/tmp
Z=$(awk '$1=="proteins"{print $2}' chunks/chunk_stats.txt)
T=14.7
p=pass1/hmm/chunk_$cid
L=logs/pass1/chunk_$cid.log
[ -e $p.done ] && exit 0
rm -f $p.failed $p.tblout $p.domtblout
echo "$(date -u +%FT%TZ) start pass1 chunk $cid pid $$ Z=$Z" >> $L
nice -n 10 $A/hmmsearch --cpu 0 --noali -T $T --domT $T -Z $Z --tblout $p.tblout --domtblout $p.domtblout \
    -o /dev/null db/pass1_148.hmm chunks/chunk_$cid.faa >> $L 2>&1
rc=$?
if [ $rc -ne 0 ] || [ "$(tail -1 $p.tblout)" != "# [ok]" ] || [ "$(tail -1 $p.domtblout)" != "# [ok]" ]; then
    echo "$(date -u +%FT%TZ) FAILED rc=$rc" >> $L; touch $p.failed; exit 1
fi
touch $p.done
echo "$(date -u +%FT%TZ) end pass1 chunk $cid" >> $L
