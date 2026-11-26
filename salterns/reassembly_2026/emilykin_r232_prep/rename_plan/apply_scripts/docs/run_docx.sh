#!/bin/bash
# Usage: run_docx.sh "<docx name in repo>" <edit script> <tag>
# Backs up the repo DOCX to <stem>.pre_r232.bak.docx (once), unzips the backup into scratch, applies the edit script,
# rezips into scratch and validates (XSD vs original + redlining by author).  Does NOT copy into the repo.
set -euo pipefail
NAME="$1"; SCRIPT="$2"; TAG="$3"
R=$HOME/Documents/EmilyKin
S=/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan/apply_scratch/docs
A=/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan/apply_scripts/docs
SK=/home/freiburger/.claude/skills/synced/49be370f-b8a1-4df3-8125-e9d4302adaab_67b17de6-c8a0-459f-85ef-8dd13a27989a/docx/scripts
PY=$HOME/Documents/py_venv/bin/python
STEM="${NAME%.docx}"
BAK="$R/$STEM.pre_r232.bak.docx"
[ -e "$BAK" ] || cp -p "$R/$NAME" "$BAK"
cmp -s "$BAK" "$R/$NAME" || { echo "repo file differs from backup: refusing to rebuild from a stale backup"; exit 3; }
W=$S/work_$TAG; rm -rf "$W"; mkdir -p "$W/u"
(cd "$W/u" && unzip -q "$BAK")
find "$W/u" -type l -delete
DATE=$(date -u +%Y-%m-%dT%H:%M:%SZ)
PYTHONDONTWRITEBYTECODE=1 nice -n 10 $PY -B "$A/$SCRIPT" "$W/u" "$DATE" | tee "$S/${TAG}_edit.log"
OUT="$W/$TAG.docx"; rm -f "$OUT"
(cd "$W/u" && zip -q -X -r -D "$OUT" '[Content_Types].xml' . )
echo "== XSD (new errors vs original)"; (cd $SK && PYTHONDONTWRITEBYTECODE=1 $PY -B office/validate.py "$OUT" --original "$BAK" 2>&1) | tee "$S/${TAG}_validate.log" | tail -8 || true
echo "== redlining (author)"; (cd $SK && PYTHONDONTWRITEBYTECODE=1 $PY -B office/validate.py "$OUT" --original "$BAK" --author "Claude (GTDB R232 update)" -v 2>&1) | tee "$S/${TAG}_redline.log" | tail -12 || true
