#!/usr/bin/env bash
# Single detached driver: wait for raw DB files -> DRAM DB setup -> one-MAG test -> 558 MAGs (P x THREADS) -> merge -> distill.
# Stage markers in $W/state/. Re-running resumes: finished stages and MAGs are skipped.
set -uo pipefail
source /scratch1/afreiburger/salterns_2026/dram/scripts/dram_env.sh
S=$W/state; mkdir -p "$S" "$W/env"; cd "$W"
log(){ echo "[$(date '+%F %T')] $*"; }
fail(){ log "FAILED: $*"; echo "$*" > "$S/FAILED"; exit 1; }
rm -f "$S/FAILED"
log "driver pid $$ host $(hostname)"

# 0. environment record
if [ ! -s "$W/env/versions.txt" ]; then
  { DRAM-setup.py version; python -c "import pandas,sqlalchemy,numpy,skbio,setuptools;print('pandas',pandas.__version__,'sqlalchemy',sqlalchemy.__version__,'numpy',numpy.__version__,'scikit-bio',skbio.__version__,'setuptools',setuptools.__version__)";
    prodigal -v 2>&1 | grep -i prodigal; mmseqs version; hmmsearch -h | sed -n 2p; tRNAscan-SE -h 2>&1 | sed -n 2p; barrnap --version 2>&1; } > "$W/env/versions.txt" 2>&1
  /scratch1/afreiburger/metag-hybrid/env/miniforge/bin/conda list -p "$E" --explicit > "$W/env/dram15.explicit.txt" 2>/dev/null
fi

# 1. DRAM databases, phase A (everything but Pfam) as soon as those raw files exist and check out
if [ ! -e "$S/SETUP_A_OK" ]; then
  until [ -e "$RAW/FETCH_OK" ] || [ -d "$RAW/Pfam-A.full.gz.parts" ] || [ -e "$RAW/Pfam-A.full.gz" ]; do
    pgrep -u "$(id -un)" -fx 'bash scripts/fetch_dbs.sh' > /dev/null || nice -n 10 bash scripts/fetch_dbs.sh >> logs/fetch_dbs.retry.log 2>&1
    sleep 60
  done
  bash scripts/verify_nopfam.sh > logs/verify_nopfam.log 2>&1 || fail "phase-A raw file checks (logs/verify_nopfam.log)"
  [ -e "$DRAMDB/database_processing.log" ] && fail "a previous setup attempt left $DRAMDB (not resumable): see README recovery"
  log "setup phase A (no Pfam) start ($SETUP_THREADS threads)"
  nice -n 10 /usr/bin/time -v python scripts/setup_dbs.py "$RAW" "$DRAMDB" "$SETUP_THREADS" nopfam > logs/setup_dbs.nopfam.log 2>&1 \
    || fail "setup phase A (logs/setup_dbs.nopfam.log)"
  touch "$S/SETUP_A_OK"; log "setup phase A OK"
fi

# 2. phase B: Pfam, once fetch_dbs.sh has verified everything (Pfam md5 included)
if [ ! -e "$S/SETUP_OK" ]; then
  tries=0
  until [ -e "$RAW/FETCH_OK" ]; do
    if pgrep -u "$(id -un)" -fx 'bash scripts/fetch_dbs.sh' > /dev/null; then sleep 120; continue; fi
    tries=$((tries+1)); [ $tries -gt 3 ] && fail "fetch_dbs.sh did not finish after 3 attempts (see logs/fetch_dbs*.log)"
    log "running fetch_dbs.sh (attempt $tries)"; nice -n 10 bash scripts/fetch_dbs.sh >> logs/fetch_dbs.retry.log 2>&1
  done
  log "raw files ready (FETCH_OK)"
  if [ ! -e "$S/SETUP_B_STARTED" ]; then
    touch "$S/SETUP_B_STARTED"
    log "setup phase B (Pfam) start ($SETUP_THREADS threads)"
    nice -n 10 /usr/bin/time -v python scripts/setup_dbs.py "$RAW" "$DRAMDB" "$SETUP_THREADS" pfam > logs/setup_dbs.pfam.log 2>&1 \
      || fail "setup phase B (logs/setup_dbs.pfam.log)"
  else
    grep -q "SETUP_FINISHED pfam" logs/setup_dbs.pfam.log || fail "phase B was started before and did not finish: see README recovery"
  fi
  DRAM-setup.py print_config > "$W/env/dram_config.txt" 2>/dev/null || fail "print_config"
  python - <<'PY' > logs/setup_check.log 2>&1 || fail "DRAM config incomplete (logs/setup_check.log)"
import json, os
from mag_annotator.database_handler import get_config_loc
c = json.load(open(get_config_loc()))
flat = {}
def walk(d):
    for k, v in d.items():
        if isinstance(v, dict): walk(v)
        else: flat[k] = v
walk({k: c[k] for k in ("search_databases", "database_descriptions", "dram_sheets")})
flat["description_db"] = c.get("description_db")
need = ["kofam_hmm", "kofam_ko_list", "pfam", "pfam_hmm", "dbcan", "viral", "peptidase", "vogdb", "dbcan_fam_activities",
        "dbcan_subfam_ec", "vog_annotations", "genome_summary_form", "module_step_form", "etc_module_database",
        "function_heatmap_form", "amg_database", "description_db"]
bad = [k for k in need if not (flat.get(k) and os.path.exists(flat[k]))]
for k in need: print(k, flat.get(k))
print("uniref", flat.get("uniref"), "kegg", flat.get("kegg"))
assert not bad, bad
PY
  log "setup OK"; touch "$S/SETUP_OK"
fi

# 3. one-MAG test (annotate + distill with the custom sheet)
if [ ! -e "$S/TEST_OK" ]; then
  bash scripts/run_one_mag.sh "$TEST_MAG" >> logs/run_all.log || fail "test annotate ($DRAM_LOGS/mags/$TEST_MAG.log)"
  rm -rf "$W/test_distill"
  T=$DRAM_OUT/$TEST_MAG
  DRAM.py distill -i "$T/annotations.tsv" -o "$W/test_distill" $( [ -s "$T/trnas.tsv" ] && echo --trna_path "$T/trnas.tsv" ) \
     $( [ -s "$T/rrnas.tsv" ] && echo --rrna_path "$T/rrnas.tsv" ) --custom_distillate "$SHEET" > logs/test_distill.log 2>&1 || fail "test distill (logs/test_distill.log)"
  python scripts/check_test.py "$DRAM_OUT/$TEST_MAG" "$W/test_distill" > logs/test_check.log 2>&1 || fail "test checks (logs/test_check.log)"
  grep -E "Elapsed|Maximum resident" "$DRAM_LOGS/mags/$TEST_MAG.log" >> logs/test_check.log
  log "test OK: $(tail -3 logs/test_check.log | tr '\n' ' ')"; touch "$S/TEST_OK"
fi

# 4. all 558 MAGs, P jobs x THREADS threads, one DRAM process per MAG (failure isolation + resume)
if [ ! -e "$S/ANNOTATE_OK" ]; then
  tail -n +2 "$MAG_LIST" | cut -f1 > "$W/inputs/mag_names.txt"
  for pass in 1 2; do
    log "annotate pass $pass"
    xargs -a "$W/inputs/mag_names.txt" -P "$P" -I{} bash scripts/run_one_mag.sh {} >> logs/run_all.log 2>&1
    nmiss=0; : > "$W/state/missing.txt"
    while read -r n; do [ -s "$DRAM_OUT/$n/annotations.tsv" ] || { echo "$n" >> "$W/state/missing.txt"; nmiss=$((nmiss+1)); }; done < "$W/inputs/mag_names.txt"
    log "pass $pass: missing $nmiss"
    [ $nmiss -eq 0 ] && break
  done
  [ $nmiss -eq 0 ] && touch "$S/ANNOTATE_OK" || log "continuing with $nmiss MAGs missing (state/missing.txt)"
fi

# 5. merge
if [ ! -e "$S/MERGE_OK" ]; then
  rm -rf "$W/merged"
  DRAM.py merge_annotations -i "$DRAM_OUT/*" -o "$W/merged" > logs/merge.log 2>&1 || fail "merge (logs/merge.log)"
  python - <<'PY' > logs/merge_check.log 2>&1 || fail "merge check (logs/merge_check.log)"
import pandas as pd, os
f = os.environ["W"] + "/merged/annotations.tsv"
cols = open(f).readline().rstrip("\n").split("\t")
a = pd.read_csv(f, sep="\t", usecols=[0, cols.index("fasta")], index_col=0, low_memory=False)
n = sum(1 for d in os.listdir(os.environ["DRAM_OUT"]) if os.path.exists(f"{os.environ['DRAM_OUT']}/{d}/annotations.tsv"))
print("genes", len(a), "MAGs in merged", a.fasta.nunique() if "fasta" in a else "?", "MAG dirs", n)
assert "fasta" in a and a.fasta.nunique() == n
PY
  log "merge OK: $(cat logs/merge_check.log)"; touch "$S/MERGE_OK"
fi

# 6. distill (default sheets + Salterns methane), plus a gene-name version of the summary
if [ ! -e "$S/DISTILL_OK" ]; then
  rm -rf "$W/distill"
  nice -n 10 DRAM.py distill -i "$W/merged/annotations.tsv" -o "$W/distill" $( [ -s "$W/merged/trnas.tsv" ] && echo --trna_path "$W/merged/trnas.tsv" ) \
     $( [ -s "$W/merged/rrnas.tsv" ] && echo --rrna_path "$W/merged/rrnas.tsv" ) --custom_distillate "$SHEET" > logs/distill.log 2>&1 || fail "distill (logs/distill.log)"
  touch "$S/DISTILL_OK"
  rm -rf "$W/distill_gene_names"
  nice -n 10 DRAM.py distill -i "$W/merged/annotations.tsv" -o "$W/distill_gene_names" --custom_distillate "$SHEET" \
     --distillate_gene_names > logs/distill_gene_names.log 2>&1 || log "gene-name distill failed (non-fatal)"
fi

nmiss=$(wc -l < "$W/state/missing.txt" 2>/dev/null || echo 0)
if [ -e "$S/ANNOTATE_OK" ]; then touch "$S/DONE"; log "DONE: all MAGs annotated, merged, distilled"
else cp "$W/state/missing.txt" "$S/DONE_INCOMPLETE"; log "DONE_INCOMPLETE: $nmiss MAGs missing (state/DONE_INCOMPLETE)"; fi
