#!/bin/bash
# bakta_test orchestrator (detached). Waits on done-markers, then: Bakta batch2 -> Prodigal -> KOfam (hmmsearch, all profiles)
# on Bakta CDS + Prodigal proteins, and eggNOG-mapper on Bakta proteins, plus an eggNOG-DB KO vocabulary check.
# Thread budget <= 48: Bakta batch1 4x8 (+ grdscan 16) -> batch2 2x8 alongside batch1's tail -> KOfam 39x1 + emapper 8 + grep 1.
set -uo pipefail
W=/scratch1/afreiburger/salterns_2026/bakta_test
ANNOT=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_
K=/scratch1/afreiburger/salterns_2026/kofam_db
EGG=/scratch1/afreiburger/emilykin/processed/refdata/eggnog
export PATH=$ANNOT/bin:$PATH
cd $W
log(){ echo "[$(date +%FT%T)] $*"; }

# 1. batch2 (chosen from the grd pre-scan) once grdscan is finished and batch2.txt is written
while ! grep -q GRDSCAN_DONE grdscan/grdscan.log 2>/dev/null || [ ! -s batch2.txt ]; do sleep 60; done
log "bakta batch2: $(tr '\n' ' ' < batch2.txt)"
bash scripts/run_bakta.sh batch2.txt 2 > logs_bakta_batch2.progress 2>&1
while ! grep -q BAKTA_DONE logs_bakta_batch1.progress; do sleep 60; done
log "all Bakta done"
cat batch1.txt batch2.txt > labels.txt

# 2. Prodigal -p meta -g 11 (DRAM's caller/mode), same FASTA as Bakta (all contigs >= 1,000 bp already)
mkdir -p prodigal
prodigal -v > logs/prodigal.version 2>&1
xargs -P 9 -I{} bash -c '/usr/bin/time -v -o logs/prodigal_{}.time nice -n 10 prodigal -p meta -g 11 -q -i fasta/{}.fna -a prodigal/{}.faa -f gff -o prodigal/{}.gff' < labels.txt
log "prodigal done"

# 3. combined protein FASTA: >bakta|<label>|<locus_tag>  and  >prodigal|<label>|<contig>_<n>
mkdir -p kofam/chunks kofam/out
: > kofam/all_proteins.faa
while read l; do
  sed "s/^>\([^ ]*\).*/>bakta|$l|\1/" bakta/$l/$l.faa >> kofam/all_proteins.faa
  sed -e "s/^>\([^ ]*\).*/>prodigal|$l|\1/" -e 's/\*$//' prodigal/$l.faa >> kofam/all_proteins.faa
done < labels.txt
log "proteins: $(grep -c '>' kofam/all_proteins.faa)"

# 4. eggNOG-mapper on Bakta proteins (one job, 8 threads, --dmnd_iterate no as in EmilyKin, --dbmem) - background
mkdir -p eggnog/tmp
: > eggnog/bakta_all.faa
while read l; do sed "s/^>\([^ ]*\).*/>$l|\1/" bakta/$l/$l.faa >> eggnog/bakta_all.faa; done < labels.txt
emapper.py --data_dir $EGG --version > logs/emapper.version 2>&1
( /usr/bin/time -v -o logs/emapper.time nice -n 10 emapper.py --data_dir $EGG -m diamond --itype proteins \
    -i eggnog/bakta_all.faa -o bakta_test --output_dir eggnog --temp_dir eggnog/tmp --cpu 8 \
    --dmnd_iterate no --dbmem --override > logs/emapper.log 2>&1; echo "EMAPPER_EXIT $?" >> logs/emapper.log ) &
EPID=$!
# eggNOG 5.0.2 KO vocabulary (I/O only, 1 thread)
( LC_ALL=C grep -a -o -E 'K(21578|21579|10670|10672|21576|21577|25220|25221|14083|16178|16176|22443|22444|20038|20037|23996|18049|05915|06163)' $EGG/eggnog.db \
    | sort | uniq -c > eggnog/eggnog_db_ko_vocab.txt; echo VOCAB_DONE >> eggnog/eggnog_db_ko_vocab.txt ) &

# 5. KOfam: every profile in the 2026-10 release, 39 chunks x hmmsearch --cpu 1, default reporting thresholds (as DRAM)
ls $K/profiles/*.hmm > kofam/profile_list.txt
split -n l/39 -d -a 2 kofam/profile_list.txt kofam/chunks/list_
for f in kofam/chunks/list_*; do xargs cat < $f > kofam/chunks/$(basename $f | sed 's/list_/chunk_/').hmm; done
log "KOfam start: $(wc -l < kofam/profile_list.txt) profiles"
/usr/bin/time -v -o logs/kofam_total.time bash -c 'ls kofam/chunks/chunk_*.hmm | xargs -P 39 -I{} bash -c "c=\$(basename {} .hmm); /usr/bin/time -v -o logs/kofam_\$c.time nice -n 10 hmmsearch --cpu 1 --domtblout kofam/out/\$c.domtblout -o /dev/null {} kofam/all_proteins.faa"'
cat kofam/out/chunk_*.domtblout | grep -v '^#' | gzip > kofam/all.domtblout.gz
log "KOfam done"
wait $EPID
wait
log "eggnog done: $(tail -1 logs/emapper.log)"
touch PIPELINE_DONE
log PIPELINE_DONE
