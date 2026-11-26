#!/usr/bin/env bash
# Step 3 (v2): HMMER 3.4 hmmsearch, v1 flags (--cpu 1 --noali -E 1e-5 --domE 1e-5 -Z 20000000), of the KOfam 2026-10-26
# profiles of the v2 panel.
#  * salterns inputs (39 per-sample gene sets + coasm520, v1 genes/ via symlink): only the profiles NOT in the v1 panel
#    (db/panel_new.hmm); 04_build_db_v2.py reads v1 hmm/<s>.tbl (62 profiles) + hmm_new/<s>.<chunk>.tbl. With the fixed
#    -Z the per-profile results are independent of which other profiles are searched, so the union equals a fresh run.
#  * ext (external references): the whole v2 panel (db/panel_v2.hmm), via   MODE=ext bash 03_panel_hmmsearch_v2.sh
# Profiles are split into NCH chunks so that NP single-CPU jobs fill the thread budget. nice 10.
set -euo pipefail
W=/scratch1/afreiburger/salterns_2026/read_markers_v2
V1=/scratch1/afreiburger/salterns_2026/read_markers
A=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin
P=/scratch1/afreiburger/salterns_2026/kofam_db/profiles
export TMPDIR=$W/tmp
NP=${NP:-80}; NCH=${NCH:-4}; MODE=${MODE:-salterns}
cd $W; mkdir -p hmm_new db
if [ ! -s db/panel_v2.hmm ]; then
  awk -F'\t' 'NR>1{print $1}' meta/ko_panel_v2.tsv | sort -u | while read k; do cat $P/$k.hmm; done > db/panel_v2.hmm
  awk -F'\t' 'NR>1{print $1}' $V1/meta/ko_panel.tsv | sort -u > db/v1_kos.txt
  awk -F'\t' 'NR>1{print $1}' meta/ko_panel_v2.tsv | sort -u | grep -v -x -f db/v1_kos.txt > db/new_kos.txt
  n=$(wc -l < db/new_kos.txt); per=$(( (n + NCH - 1) / NCH ))
  split -l $per -d db/new_kos.txt db/new_kos.chunk.
  for c in db/new_kos.chunk.*; do while read k; do cat $P/$k.hmm; done < $c > db/panel_new.${c##*.}.hmm; done
  echo "v2 panel $(grep -c '^NAME' db/panel_v2.hmm) profiles; new $(cat db/panel_new.*.hmm | grep -c '^NAME') in $(ls db/panel_new.*.hmm | wc -l) chunks"
fi
if [ "$MODE" = salterns ]; then
  jobs=$(for s in $(ls /scratch1/afreiburger/salterns_2026/rollout/reads) coasm520; do for h in db/panel_new.*.hmm; do c=${h#db/panel_new.}; c=${c%.hmm}; echo "$s $c"; done; done)
  echo "$jobs" | while read s c; do [ -f hmm_new/$s.$c.done ] || echo "$s $c"; done | \
    xargs -P $NP -L 1 bash -c 'nice -n 10 '$A'/hmmsearch --cpu 1 --noali -E 1e-5 --domE 1e-5 -Z 20000000 \
      --tblout hmm_new/$0.$1.tbl --domtblout hmm_new/$0.$1.domtbl -o /dev/null db/panel_new.$1.hmm '$V1'/genes/$0/$0.faa && touch hmm_new/$0.$1.done'
  echo "$(date) salterns hmmsearch done: $(ls hmm_new/*.done | wc -l) jobs" > hmm_new/SALTERNS_DONE
else
  nice -n 10 $A/seqkit split2 -p 16 -O refs/ext_split refs/ext.faa 2> /dev/null
  ls refs/ext_split/*.faa | xargs -P 16 -I{} bash -c 'f={}; nice -n 10 '$A'/hmmsearch --cpu 1 --noali -E 1e-5 --domE 1e-5 -Z 20000000 --tblout $f.tbl --domtblout $f.domtbl -o /dev/null db/panel_v2.hmm $f'
  cat refs/ext_split/*.tbl | grep -v '^#' > hmm_new/ext.tbl; cat refs/ext_split/*.domtbl | grep -v '^#' > hmm_new/ext.domtbl
  rm -rf refs/ext_split; echo "$(date) ext hmmsearch done" > hmm_new/EXT_DONE
fi
