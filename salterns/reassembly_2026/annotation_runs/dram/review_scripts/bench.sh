#!/usr/bin/env bash
# Reviewer benchmark: does DRAM's per-MAG KOfam hmmsearch use more than ~1 core? 1,000 KOfam profiles vs the test MAG's 2,279 genes.
source /scratch1/afreiburger/salterns_2026/dram/scripts/dram_env.sh
B=$W/review/bench; mkdir -p $B; cd $B
G=$W/review/annot_concoct.bin.101.fasta__R1_sediment/working_dir/concoct.bin.101.fasta__R1_sediment/tmp/genes.faa
cp $G genes.faa
[ -s sub1000.hmm.h3m ] || { awk '/^HMMER3/{n++} n>1000{exit} {print}' $DRAMDB/kofam_profiles.hmm > sub1000.hmm; hmmpress -f sub1000.hmm > /dev/null; }
echo "profiles: $(grep -c ^NAME sub1000.hmm)"
for c in 1 3 6; do
  /usr/bin/time -f "cpu=$c wall=%e s user+sys=%U+%S s cpu%%=%P" nice -n 10 hmmsearch --domtblout dt_$c.txt --cpu $c -o /dev/null sub1000.hmm genes.faa
done 2>&1
echo BENCH_DONE
