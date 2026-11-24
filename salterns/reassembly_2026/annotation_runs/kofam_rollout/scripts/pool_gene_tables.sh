#!/bin/bash
# Pool the 935 per-MAG gene tables (kofam_rollout copy; only the work dir changed) (gene -> contig -> contig_len -> MAG -> binner ...) into one gzipped table.
cd /scratch1/afreiburger/salterns_2026/kofam_rollout
{ head -1 genes/$(ls genes | grep -m1 "genes.tsv$"); for f in genes/*.genes.tsv; do tail -n +2 $f; done; } | nice -n 10 gzip > genes_all.tsv.gz.part && mv genes_all.tsv.gz.part genes_all.tsv.gz
