#!/bin/bash
# Selection pre-scan: which new Desulfobacterota / Bacillota bins carry grdH/grdI (and mtgB/mttB)?
# prodigal -p meta -g 11 on each bin (16 parallel, 1 thread each), then hmmsearch with 8 KOfam profiles (2026-10 release).
set -uo pipefail
W=/scratch1/afreiburger/salterns_2026/bakta_test/grdscan
ANNOT=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_
K=/scratch1/afreiburger/salterns_2026/kofam_db
export PATH=$ANNOT/bin:$PATH W
cd $W; mkdir -p faa
ls fasta/*.fna | xargs -P 16 -I{} bash -c 'b=$(basename {} .fna); [ -s faa/$b.faa ] || { nice -n 10 prodigal -p meta -g 11 -q -i {} -a faa/$b.tmp.faa -o /dev/null && sed "s/^>/>${b}|/" faa/$b.tmp.faa > faa/$b.faa && rm faa/$b.tmp.faa; }'
cat faa/*.faa | sed 's/\*$//' > all.faa
echo "PRODIGAL_DONE $(date +%FT%T) $(grep -c '>' all.faa) proteins"
while [ ! -e $K/.download_done ]; do sleep 30; done
cat $K/profiles/{K21578,K21579,K10670,K10672,K21576,K21577,K25220,K14083}.hmm > grd.hmm
nice -n 10 hmmsearch --cpu 16 --domtblout grd.domtblout -o /dev/null grd.hmm all.faa
echo "GRDSCAN_DONE $(date +%FT%T)"
