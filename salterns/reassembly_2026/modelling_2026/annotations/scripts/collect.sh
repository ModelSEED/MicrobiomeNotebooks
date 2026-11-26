#!/bin/bash
# Collect annot57 outputs from poplar into modelling_2026/annotations/ and verify them by md5 (run locally).
#   poplar annot57/bakta/<mag_id>/*      -> annotations/bakta/<mag_id>/
#   poplar annot57/eggnog/g{1,2}/*       -> annotations/eggnog/_pooled/g{1,2}/   (input.faa excluded: = prodigal/*.faa minus '*')
#   poplar annot57/prodigal/*            -> annotations/prodigal/
#   poplar annot57/logs/*, driver.log, monitor.log, fasta.md5 -> annotations/logs/
# md5: computed on poplar (find -type f | md5sum), path-mapped, then md5sum -c here.
set -euo pipefail
A=/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/modelling_2026/annotations
W=/scratch1/afreiburger/salterns_2026/annot57
SSH="ssh -o BatchMode=yes"
cd $A
mkdir -p bakta eggnog/_pooled prodigal logs md5
$SSH poplar "cd $W && mkdir -p md5 && find bakta eggnog/g1 eggnog/g2 prodigal logs -type f ! -name input.faa -print0 | xargs -0 md5sum > md5/poplar_md5.txt && md5sum driver.log monitor.log fasta.md5 >> md5/poplar_md5.txt && wc -l md5/poplar_md5.txt"
rsync -a -e "$SSH" poplar:$W/md5/poplar_md5.txt md5/
rsync -a -e "$SSH" poplar:$W/bakta/ bakta/
rsync -a -e "$SSH" --exclude input.faa poplar:$W/eggnog/g1 poplar:$W/eggnog/g2 eggnog/_pooled/
rsync -a -e "$SSH" poplar:$W/prodigal/ prodigal/
rsync -a -e "$SSH" poplar:$W/logs/ poplar:$W/driver.log poplar:$W/monitor.log poplar:$W/fasta.md5 logs/
# map poplar paths to local ones and check
sed -e 's#  eggnog/g#  eggnog/_pooled/g#' -e 's#  driver.log$#  logs/driver.log#' -e 's#  monitor.log$#  logs/monitor.log#' \
    -e 's#  fasta.md5$#  logs/fasta.md5#' md5/poplar_md5.txt > md5/local_paths_md5.txt
md5sum -c --quiet md5/local_paths_md5.txt && echo "MD5_OK $(wc -l < md5/local_paths_md5.txt) files"
