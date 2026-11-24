#!/bin/bash
# Pool the 38 per-MAG gene tables (kofam_coasm38 copy of kofam_rollout's; only the work dir changed) (gene -> contig -> contig_len -> MAG -> binner ...) into one gzipped table.
cd /scratch1/afreiburger/salterns_2026/kofam_coasm38
{ head -1 genes/$(ls genes | grep -m1 "genes.tsv$"); for f in genes/*.genes.tsv; do tail -n +2 $f; done; } | nice -n 10 gzip > genes_all.tsv.gz.part && mv genes_all.tsv.gz.part genes_all.tsv.gz
