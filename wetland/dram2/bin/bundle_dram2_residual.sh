#!/usr/bin/env bash
# Bundle every DRAM2 output that is not already in the repo (per Andrew 2026-11-18: include the residual files of the
# annotation jobs, coding and non-coding). One .tar.gz per chunk x output dir, plus the merged per-MAG tRNA/rRNA
# tables, the final all_mags run (distill/product/tRNA_rRNA), logs and launch configs. Uses pigz -p 32.
set -uo pipefail
D=/scratch1/afreiburger/wetland_dram2
O=/scratch/afreiburger/dram2_residual
mkdir -p "$O"
ev() { echo "[$(date '+%F %T')] $*" >> "$O/bundle.events"; }
pack() {  # pack <tarball name> <dir relative to D>
  local out="$O/$1.tar.gz" src=$2
  [ -s "$out" ] && return 0
  tar -C "$D" -cf - "$src" | pigz -p 32 > "$out.tmp" && mv "$out.tmp" "$out" && ev "packed $1 ($(du -sh "$out" | cut -f1))" || ev "FAIL $1"
}
ev "START"
for c in 01 02 03 04 05; do
  for sub in tRNA_rRNA Prodigal_v2.6.3 HMM_search MMSEQS2 RAW QUAST_v5.2.0 Processed_inputs multiqc pipeline_info HMM_search_kofam_raw; do
    pack "chunk_${c}__${sub}" "runs/chunk_${c}/out/${sub}"
  done
  pack "chunk_${c}__run_meta" "runs/chunk_${c}/command.sh"
done
pack "merged__trnas" "merged/trnas"
pack "merged__rrnas" "merged/rrnas"
pack "all_mags__out" "runs/all_mags/out"
pack "logs" "logs"
pack "launch" "launch"
( cd "$O" && sha256sum *.tar.gz > SHA256SUMS && ls -l *.tar.gz | awk '{print $5"\t"$9}' > SIZES.tsv )
ev "DONE $(ls "$O"/*.tar.gz | wc -l) tarballs, $(du -sh "$O" | cut -f1)"
