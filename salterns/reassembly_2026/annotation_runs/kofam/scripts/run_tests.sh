#!/bin/bash
# Pre-launch tests on chunk 01 (3 single-CPU jobs, nice 10):
#  (a) marker profiles x chunk 01 + postprocess  -> correctness check vs annotation_choice/kofam_control_test
#  (b) 300 random profiles x chunk 01            -> per-chunk timing (layout A: protein chunks x all profiles)
#  (c) 8 of those profiles x ALL proteins        -> layout B (profile chunks x all proteins), same work per profile-residue
set -uo pipefail
W=/scratch1/afreiburger/salterns_2026/kofam
ANNOT=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_
H=$ANNOT/bin/hmmsearch
cd $W
Z=$(awk '$1=="proteins"{print $2}' chunks/chunk_stats.txt); T=14.7
mkdir -p test/markers/hmm test/markers/results test/timing
for d in db meta chunks genes; do ln -sfn ../../$d test/markers/$d; done
$ANNOT/bin/python3 scripts/make_test_hmms.py /scratch1/afreiburger/salterns_2026/kofam_db db > test/timing/subsets.txt
( /usr/bin/time -v nice -n 10 $H --cpu 1 --noali -T $T --domT $T -Z $Z --tblout test/markers/hmm/chunk_01.tblout \
    --domtblout test/markers/hmm/chunk_01.domtblout -o test/markers/hmm/chunk_01.out db/markers.hmm chunks/chunk_01.faa \
    > test/markers/time.log 2>&1 \
  && nice -n 10 $ANNOT/bin/python3 scripts/postprocess.py chunk $W/test/markers 01 > test/markers/postprocess.log 2>&1 ) &
/usr/bin/time -v nice -n 10 $H --cpu 1 --noali -T $T --domT $T -Z $Z --tblout test/timing/r300_chunk01.tblout \
    -o test/timing/r300_chunk01.out db/random300.hmm chunks/chunk_01.faa > test/timing/r300_chunk01.time 2>&1 &
cat chunks/chunk_*.faa > /tmp/afr_kofam_all.faa
/usr/bin/time -v nice -n 10 $H --cpu 1 --noali -T $T --domT $T -Z $Z --tblout test/timing/r8_all.tblout \
    -o test/timing/r8_all.out db/random8.hmm /tmp/afr_kofam_all.faa > test/timing/r8_all.time 2>&1 &
wait
rm -f /tmp/afr_kofam_all.faa
touch test/TESTS_DONE
