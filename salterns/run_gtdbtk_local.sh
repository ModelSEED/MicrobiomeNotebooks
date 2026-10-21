#!/usr/bin/env bash
# Local GTDB-Tk fallback for the 310 dereplicated salterns MAGs.
#
# Run only if the KBase GTDB-Tk runs fail (queue_mag_qc_jobs.py). Approved by the user on
# 2026-10-20 on that condition. GTDB-Tk 2.8.0 needs reference data r232 and pins numpy~=1.9 /
# pydantic~=1.9, which would downgrade ~/Documents/py_venv (numpy 2.4.6, pydantic 2.11.9), so it
# gets its own uv venv. Its non-Python tools come from Bioconda tarballs (tested on this machine):
#   prodigal 2.6.3, HMMER 3.4, pplacer/guppy 1.1.alpha19 (static), skani 0.4.0
# Input: the 310 derep FASTAs in MAG_contigs/ (index.tsv, population=derep310).
# Output: MAG_QC/gtdbtk_local_r232/ (gtdbtk.bac120.summary.tsv, gtdbtk.ar53.summary.tsv).
# Every step is skipped when already done, so a re-run resumes.
#
# Usage: setsid nohup bash run_gtdbtk_local.sh > MAG_QC/gtdbtk_local_r232.log 2>&1 &
#        CPUS=48 PPLACER_CPUS=16 bash run_gtdbtk_local.sh   # defaults shown
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
VENV="${VENV:-$HOME/Documents/gtdbtk_venv}"
DATA="${DATA:-$HOME/Documents/gtdbtk_data}"
OUT="$HERE/MAG_QC/gtdbtk_local_r232"
# Keep tmp short: multiprocessing puts a Unix socket under it, and AF_UNIX paths are capped at
# 107 bytes — a tmpdir inside $OUT (84 chars) makes Manager() die ("AF_UNIX path too long").
TMP="${TMP_DIR:-/tmp/gtdbtk_salterns}"
CPUS="${CPUS:-48}"
PPLACER_CPUS="${PPLACER_CPUS:-16}"
REF_URL="https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz"
REF_BYTES=60806405195
BIOCONDA="https://conda.anaconda.org/bioconda/linux-64"
PKGS=(prodigal-2.6.3-h516909a_2.tar.bz2 hmmer-3.4-hdbdd923_2.tar.bz2
      pplacer-1.1.alpha19-h9ee0642_2.tar.bz2 skani-0.4.0-hec9b1f2_0.conda)

log() { echo "[$(date '+%F %T')] $*"; }

# ---- 1. dedicated venv with GTDB-Tk -------------------------------------------------------
if [ ! -x "$VENV/bin/gtdbtk" ]; then
  log "creating $VENV with gtdbtk==2.8.0"
  [ -x "$VENV/bin/python" ] || uv venv "$VENV" --python 3.12
  uv pip install --python "$VENV/bin/python" gtdbtk==2.8.0
fi

# ---- 2. external tools from Bioconda into the venv's bin/ ---------------------------------
for f in "${PKGS[@]}"; do
  name="${f%%-*}"
  dir="$VENV/pkgs/$name"
  [ -f "$dir/.done" ] && continue
  log "installing $f"
  rm -rf "$dir" && mkdir -p "$dir"
  curl -sSfL -o "$dir/$f" "$BIOCONDA/$f"
  case "$f" in
    *.tar.bz2) tar -xjf "$dir/$f" -C "$dir" ;;
    *.conda)   mkdir -p "$dir/z" && unzip -q -o "$dir/$f" -d "$dir/z"
               for t in "$dir"/z/pkg-*.tar.zst; do zstd -dq --stdout "$t" | tar -x -C "$dir"; done ;;
  esac
  cp -f "$dir"/bin/* "$VENV/bin/"
  touch "$dir/.done"
done
export PATH="$VENV/bin:$PATH"
for b in prodigal hmmsearch hmmalign pplacer guppy skani; do
  command -v "$b" >/dev/null || { log "missing $b"; exit 1; }
  if ldd "$(command -v "$b")" 2>&1 | grep -q "not found"; then log "$b has unresolved libraries"; exit 1; fi
done

# ---- 3. GTDB r232 reference data (60.8 GB) ----------------------------------------------------
mkdir -p "$DATA"
tarball="$DATA/gtdbtk_r232_data.tar.gz"
if [ ! -f "$DATA/.r232_extracted" ]; then
  if [ "$(stat -c %s "$tarball" 2>/dev/null || echo 0)" -ne "$REF_BYTES" ]; then
    log "downloading r232 reference data (resumable)"
    curl -sSfL -C - -o "$tarball" "$REF_URL"
  fi
  [ "$(stat -c %s "$tarball")" -eq "$REF_BYTES" ] || { log "r232 tarball size mismatch"; exit 1; }
  log "extracting r232 reference data"
  tar -xzf "$tarball" -C "$DATA"
  touch "$DATA/.r232_extracted"
fi
GTDBTK_DATA_PATH="$(find "$DATA" -maxdepth 2 -name metadata -type d -printf '%h\n' -quit)"
export GTDBTK_DATA_PATH
[ -n "$GTDBTK_DATA_PATH" ] || { log "could not locate the extracted r232 data"; exit 1; }
log "GTDBTK_DATA_PATH=$GTDBTK_DATA_PATH"
gtdbtk check_install

# ---- 4. classify the 310 dereplicated MAGs ------------------------------------------------------
mkdir -p "$OUT" "$TMP"  # gtdbtk rejects a --tmpdir that does not already exist
export TMPDIR="$TMP"
batch="$OUT/batchfile.tsv"
awk -F'\t' -v d="$HERE/MAG_contigs" 'NR > 1 && $2 == "derep310" { print d "/" $8 "\t" $4 }' \
    "$HERE/MAG_contigs/index.tsv" > "$batch"
[ "$(wc -l < "$batch")" -eq 310 ] || { log "batchfile does not list 310 MAGs"; exit 1; }
if [ ! -f "$OUT/.classify_done" ]; then
  log "gtdbtk classify_wf on 310 MAGs (cpus=$CPUS, pplacer_cpus=$PPLACER_CPUS)"
  gtdbtk classify_wf --batchfile "$batch" --out_dir "$OUT" --cpus "$CPUS" \
      --pplacer_cpus "$PPLACER_CPUS" --tmpdir "$TMP"
  touch "$OUT/.classify_done"
fi
log "done: $(ls "$OUT"/gtdbtk.*.summary.tsv 2>/dev/null | tr '\n' ' ' || true)"
