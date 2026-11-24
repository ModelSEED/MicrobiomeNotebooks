#!/usr/bin/env bash
# Reviewer end-to-end test (2026-11-11) of every DRAM 1.5 annotate code path before the real Pfam build finishes:
# a mini Pfam profile DB (3,000 random Pfam 39.0 families < 3 MB, cut from the finished dram_db/pfam.mmsmsa with the same
# msa2profile/createindex settings DRAM uses), a private copy of the description DB with Pfam descriptions filled by
# DRAM's own populate_description_db, and a private config passed with --config_loc. Writes only to $W/review/e2e/.
set -uo pipefail
source /scratch1/afreiburger/salterns_2026/dram/scripts/dram_env.sh
X=$W/review/e2e; mkdir -p "$X/tmp"; export TMPDIR=$X/tmp; cd "$X"
log(){ echo "[$(date '+%F %T')] $*"; }
log "e2e start pid $$"
if [ ! -s pfam_sub.mmspro.dbtype ]; then
  awk '$3 < 3000000 {print $1}' $DRAMDB/pfam.mmsmsa.index | shuf --random-source=<(yes) -n 3000 | sort -n > keys.txt
  nice -n 10 mmseqs createsubdb keys.txt $DRAMDB/pfam.mmsmsa pfam_sub.mmsmsa > mini_pfam.log 2>&1; rc=$?; log "createsubdb rc=$rc"
  nice -n 10 mmseqs msa2profile pfam_sub.mmsmsa pfam_sub.mmspro --match-mode 1 --threads $THREADS >> mini_pfam.log 2>&1; rc=$?; log "msa2profile rc=$rc"
  nice -n 10 mmseqs createindex pfam_sub.mmspro tmp_idx -k 5 -s 7 --threads $THREADS >> mini_pfam.log 2>&1; rc=$?; log "createindex rc=$rc"
fi
cp $DRAMDB/description_db.sqlite description_db.sqlite
python -I - <<'PY'
import json, os
from mag_annotator.database_handler import get_config_loc, DatabaseHandler
X = os.getcwd(); D = os.environ["DRAMDB"]
c = json.load(open(get_config_loc()))
c["search_databases"]["pfam"] = f"{X}/pfam_sub.mmspro"
c["database_descriptions"]["pfam_hmm"] = f"{D}/Pfam-A.hmm.dat.gz"
c["description_db"] = f"{X}/description_db.sqlite"
c["log_path"] = f"{X}/config.log"
c["setup_info"]["pfam"] = {"name": "Pfam db (review mini subset)"}
json.dump(c, open(f"{X}/CONFIG.json", "w"), indent=2)
h = DatabaseHandler(None, f"{X}/CONFIG.json")
h.populate_description_db(select_db=["pfam"])
n = h.session.execute(__import__("sqlalchemy").text("select count(*) from pfam_description")).scalar()
print("pfam descriptions in copy:", n)
PY
log "private config + description DB rc=$?"
n=$TEST_MAG; out=$X/annot_$n; rm -rf "$out"
nice -n 10 /usr/bin/time -v DRAM.py annotate -i "$MAGS/$n.fna" -o "$out" --min_contig_size $MIN_CONTIG \
  --threads $THREADS --gtdb_taxonomy "$TAX" --checkm_quality "$QUAL" --config_loc "$X/CONFIG.json" > annotate.log 2>&1
rc=$?; log "annotate rc=$rc"
rm -rf distill
nice -n 10 DRAM.py distill -i "$out/annotations.tsv" -o distill $( [ -s "$out/trnas.tsv" ] && echo --trna_path "$out/trnas.tsv" ) \
  $( [ -s "$out/rrnas.tsv" ] && echo --rrna_path "$out/rrnas.tsv" ) --custom_distillate "$SHEET" > distill.log 2>&1
rc=$?; log "distill rc=$rc"
python $W/scripts/check_test.py "$out" distill > check_test.log 2>&1
rc=$?; log "check_test rc=$rc (Pfam fraction is low by design: 3,000 of 53,275 families)"
echo E2E_FINISHED
