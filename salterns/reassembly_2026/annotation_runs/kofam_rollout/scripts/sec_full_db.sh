#!/bin/bash
# Re-score the Sec read-through joined proteins of the 935 rollout MAGs against all 27,981 KOfam 2026-10 profiles
# (the bakta_test profile chunks, ../bakta_test/kofam/chunks/chunk_NN.hmm, 39 files whose concatenation is
# byte-identical to kofam/db/kofam_all.hmm; checked 2026-11-19, see kofam_rollout/README.md) with this track's
# settings, and scan a 603-nt window around each TGA for SECIS with Bakta's ncRNA-regions CMs
# (Bakta 1.12 ncRNA-region flags). 8 x hmmsearch --cpu 1 + cmscan --cpu 2, nice 10.
# kofam_rollout copy of kofam/scripts/sec_full_db.sh: W = kofam_rollout/sec_readthrough; -Z = the 935-MAG total
# protein count, read from kofam_rollout/chunks/chunk_stats.txt (the kofam track hard-coded 3911146).
# Inputs in $W (copied from the local build step): all.joined.faa, all.windows.fna.
set -euo pipefail
R=/scratch1/afreiburger/salterns_2026/kofam_rollout
W=$R/sec_readthrough
ANNOT=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_
CH=/scratch1/afreiburger/salterns_2026/bakta_test/kofam/chunks
DB=/scratch1/afreiburger/emilykin/processed/refdata/bakta/db
Z=$(awk '$1=="proteins"{print $2}' $R/chunks/chunk_stats.txt)
: "${Z:?no proteins line in chunk_stats.txt}"
[ "$(ls $CH/chunk_*.hmm | wc -l)" -eq 39 ] || { echo "expected 39 profile chunks in $CH"; exit 1; }
export PATH=$ANNOT/bin:$PATH W Z
cd $W; mkdir -p out
echo "$(date +%FT%T) start Z=$Z" > SEC_START
nice -n 10 cmscan --noali --cut_tc -g --nohmmonly --rfam --cpu 2 -Z 5 --tblout all.windows.cmscan.tbl -o /dev/null \
     $DB/ncRNA-regions all.windows.fna &
ls $CH/chunk_*.hmm | xargs -P 8 -I{} bash -c 'c=$(basename {} .hmm); nice -n 10 hmmsearch --cpu 1 --noali -T 14.7 --domT 14.7 -Z $Z --tblout out/$c.tblout --domtblout out/$c.domtblout -o /dev/null {} all.joined.faa'
wait
cat out/chunk_*.tblout | grep -v '^#' | gzip > all.joined.full.tblout.gz
cat out/chunk_*.domtblout | grep -v '^#' | gzip > all.joined.full.domtblout.gz
echo "$(date +%FT%T) done $(ls out/*.tblout | wc -l) chunks Z=$Z" > SEC_DONE
