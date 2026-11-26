#!/usr/bin/env bash
# Download the raw reads for the tadpole pilot from NCBI SRA (AWS open-data mirror) and convert them to paired FASTQ.
# ENA https was tried first (2026-11-11) but ran at 0.01-1.3 MB/s from poplar; the NCBI S3 mirror gives about 50 MB/s.
# Integrity: vdb-validate (the .sra file's internal MD5/CRC checksums), then the pair count and base count
# of the FASTQ must equal ENA's read_count / base_count (ena_filereport.tsv, fetched from the ENA portal API).
# ENA's fastq_md5 cannot be compared because fasterq-dump writes different header lines from ENA's FASTQ.
# Usage (on poplar, inside /scratch1/afreiburger/salterns_2026/tadpole_pilot/reads):
#   bash ../scripts/download.sh ena_filereport.tsv
# Writes <run>_1.fastq, <run>_2.fastq, <run>.ok (or <run>.FAIL), and download.done when all runs check out.
set -euo pipefail
REPORT=${1:?ena_filereport.tsv}
SRA=/scratch1/afreiburger/metag-hybrid/env/miniforge/envs/sra/bin
fail=0
while IFS=$'\t' read -r run ftp md5 bytes nreads nbases; do
  [ "$run" = run_accession ] && continue
  if [ -e "$run.ok" ]; then echo "$run already done"; continue; fi
  echo "[$(date +%F' '%T)] $run: download .sra"
  nice -n 10 curl -s --retry 5 --retry-delay 20 -C - -o "$run.sra" "https://sra-pub-run-odp.s3.amazonaws.com/sra/$run/$run"
  echo "[$(date +%F' '%T)] $run: vdb-validate"
  if ! nice -n 10 "$SRA/vdb-validate" "$run.sra" > "$run.vdb-validate.log" 2>&1; then
    echo "$run vdb-validate FAILED"; cat "$run.vdb-validate.log"; touch "$run.FAIL"; fail=1; continue
  fi
  tail -3 "$run.vdb-validate.log"
  echo "[$(date +%F' '%T)] $run: fasterq-dump"
  mkdir -p tmp_fqd
  nice -n 10 "$SRA/fasterq-dump" --split-files -e 6 -t tmp_fqd -O . "$run.sra"
  read -r n1 b1 < <(awk 'NR%4==2{n++; b+=length($0)} END{print n, b}' "${run}_1.fastq")
  read -r n2 b2 < <(awk 'NR%4==2{n++; b+=length($0)} END{print n, b}' "${run}_2.fastq")
  echo "$run pairs _1=$n1 _2=$n2 (ENA read_count $nreads); bases $((b1+b2)) (ENA base_count $nbases)"
  if [ "$n1" = "$nreads" ] && [ "$n2" = "$nreads" ] && [ $((b1+b2)) = "$nbases" ]; then
    touch "$run.ok"; rm -f "$run.sra"
  else
    echo "$run count MISMATCH"; touch "$run.FAIL"; fail=1
  fi
done < "$REPORT"
rmdir tmp_fqd 2>/dev/null || true
[ $fail = 0 ] && touch download.done
echo "[$(date +%F' '%T)] download finished, fail=$fail"
