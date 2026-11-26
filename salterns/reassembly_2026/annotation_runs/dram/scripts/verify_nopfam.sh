#!/usr/bin/env bash
# Checks on the non-Pfam raw files before DRAM setup phase A. Exit 1 on any failure.
source /scratch1/afreiburger/salterns_2026/dram/scripts/dram_env.sh
cd "$RAW" || exit 1
die(){ echo "CHECK FAILED: $*"; exit 1; }
for f in kofam_profiles.tar.gz kofam_ko_list.gz dbCAN-HMMdb-V11.txt CAZyDB.08062022.fam-activities.txt \
         CAZyDB.08062022.fam.subfam.ec.txt pepunit.lib viral.1.protein.faa.gz vog.hmm.tar.gz vog.annotations.tsv.gz; do
  [ -s "$f" ] || die "missing $f"; done
[ "$(md5sum kofam_ko_list.gz | cut -d' ' -f1)" = 697af5b3b22a919306a1e60f99cec1a6 ] || die "ko_list is not the 2026-10-26 release"
for f in vog.hmm.tar.gz vog.annotations.tsv.gz; do
  want=$(curl -sS --fail -m 120 "https://fileshare.lisc.univie.ac.at/vog/latest/$f.md5" | grep -o '^[0-9a-f]\{32\}')
  [ -n "$want" ] && [ "$(md5sum $f | cut -d' ' -f1)" = "$want" ] || die "md5 $f"; done
gzip -t kofam_profiles.tar.gz || die "gzip kofam"; gzip -t viral.1.protein.faa.gz || die "gzip viral"
head -1 dbCAN-HMMdb-V11.txt | grep -q '^HMMER3' || die "dbCAN not HMMER3"
head -c 1 pepunit.lib | grep -q '>' || die "pepunit.lib not FASTA"
grep -qi '<html' CAZyDB.08062022.fam-activities.txt CAZyDB.08062022.fam.subfam.ec.txt && die "dbCAN tables are HTML"
echo "kofam profiles: $(tar -tzf kofam_profiles.tar.gz | grep -c '\.hmm$')  VOG hmms: $(tar -tzf vog.hmm.tar.gz | grep -c 'VOG.*\.hmm$')  dbCAN models: $(grep -c '^NAME' dbCAN-HMMdb-V11.txt)"
echo "NOPFAM_RAW_OK"
