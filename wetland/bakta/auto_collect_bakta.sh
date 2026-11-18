#!/usr/bin/env bash
# After the full Bakta run on the 1,184 DRAM2 representatives finishes (bakta.events "ALL done=<n>/1184"), automatically:
#  1. if any MAG failed: relaunch run_bakta_reps.sh once (it skips DONE MAGs) and wait for that pass;
#  2. bundle every Bakta output on poplar (per MAG: .tsv .gff3 .gbff .embl .json .faa .ffn .fna .txt .log
#     .hypotheticals.*, + logs/) into tarballs of 100 MAGs (pigz -p 32), with SHA256SUMS;
#  3. pull the bundle to bakta/results/residual/ (gitignored) and verify the checksums; pull the per-MAG .tsv tables;
#  4. run merge_bakta.py (coding + non-coding gene table, per-MAG counts, Bakta vs DRAM2) and commit + push the
#     per-MAG tables, the comparison, summary.json and RESIDUAL_MANIFEST.tsv.
# Local, detached. Gives up after 5 days. Log: ~/gtdbtk_kbase/auto_collect_bakta.log
# Usage: setsid nohup bash bakta/auto_collect_bakta.sh > /dev/null 2>&1 < /dev/null &
set -uo pipefail
REPO=$HOME/Documents/MicrobiomeNotebooks
H=$REPO/wetland/bakta
PY="$HOME/Documents/py_venv/bin/python -I"
LOG=$HOME/gtdbtk_kbase/auto_collect_bakta.log
P=/scratch1/afreiburger/wetland_bakta
BUN=/scratch/afreiburger/bakta_bundle
BR=clauded-insights-figure7
N=$(grep -c . "$H/reps.txt")
log() { echo "[$(date '+%F %T')] $*" >> "$LOG"; }
rsh() { timeout 600 ssh -o BatchMode=yes poplar "$@" < /dev/null; }
TRAILER=$'\n\nCo-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>'
deadline=$(( $(date +%s) + 5 * 86400 ))
log "START pid $$ (waiting for 'ALL done=*/$N')"
retried=0; need=1          # number of "ALL done=*/N" lines to wait for (2 after one retry pass)
while :; do
  [ "$(date +%s)" -gt "$deadline" ] && { log "GIVE UP after 5 days"; exit 1; }
  # only the events of the latest launch attempt count (earlier attempts may have failed and been relaunched)
  seg=$(rsh "awk '/LAUNCHER copying/{buf=\"\"} {buf=buf \$0 \"\\n\"} END{printf \"%s\", buf}' $P/bakta.events")
  if echo "$seg" | grep -qE 'LAUNCHER FAIL|FAIL wrong|FAIL DB'; then
    [ "${warned:-0}" = 1 ] || log "WAITING: the latest launch failed: $(echo "$seg" | grep -E 'FAIL' | tail -1)"
    warned=1; sleep 600; continue
  fi
  warned=0
  nall=$(echo "$seg" | grep -cE "ALL done=[0-9]+/$N ")
  running=$(rsh "ps -eo args | grep -c '^bash run_bakta_reps.sh'")
  if [ "${nall:-0}" -ge $need ] && [ "${running:-1}" = 0 ]; then
    last=$(echo "$seg" | grep -E "ALL done=[0-9]+/$N " | tail -1)
    failed=$(echo "$last" | sed -n 's/.*failed=\([0-9]*\).*/\1/p')
    if [ "${failed:-0}" -gt 0 ] && [ $retried = 0 ]; then
      log "$last -> relaunching once for the failed MAGs"
      rsh "cd $P && (setsid nohup bash run_bakta_reps.sh 16 6 >> logs/driver.log 2>&1 < /dev/null &)"
      retried=1; need=$((nall + 1)); sleep 600; continue
    fi
    log "final: $last"
    break
  fi
  sleep 600
done

# 2. bundle on poplar
rsh "mkdir -p $BUN && cd $P/runs && ls -d */ | sed 's#/##' | sort > $BUN/mags.txt && split -l 100 -d -a 2 $BUN/mags.txt $BUN/grp_ && \
     for g in $BUN/grp_*; do b=bakta_runs_\$(basename \$g | sed 's/grp_//'); [ -s $BUN/\$b.tar.gz ] || \
       { tar -cf - -T \$g | pigz -p 32 > $BUN/\$b.tar.gz.tmp && mv $BUN/\$b.tar.gz.tmp $BUN/\$b.tar.gz; }; done && \
     (cd $P && tar -cf - logs bakta.events run_bakta_reps.sh launch_after_db_copy.sh copy_db_parallel.sh reps.txt | pigz -p 8 > $BUN/bakta_logs_scripts.tar.gz) && \
     cd $BUN && sha256sum *.tar.gz > SHA256SUMS" || { log "FAIL bundling"; exit 1; }
log "bundled: $(rsh "ls $BUN/*.tar.gz | wc -l; du -sh $BUN | cut -f1" | tr '\n' ' ')"
# 3. pull + verify; per-MAG tables for the merge
D=$H/results/residual; mkdir -p "$D" "$H/results/runs_tsv"
rsync -a "poplar:$BUN/" "$D/" >> "$LOG" 2>&1 && (cd "$D" && sha256sum -c --quiet SHA256SUMS >> "$LOG" 2>&1) || { log "FAIL pull/verify"; exit 1; }
rsync -a --include='*/' --include='*.tsv' --include='DONE' --exclude='*' "poplar:$P/runs/" "$H/results/runs_tsv/" >> "$LOG" 2>&1 \
  || { log "FAIL pulling the per-MAG tables"; exit 1; }
M=$H/results/RESIDUAL_MANIFEST.tsv
{ printf "tarball\tbytes\tsha256\tmags_or_files\n"
  for f in "$D"/*.tar.gz; do
    printf "%s\t%s\t%s\t%s\n" "$(basename "$f")" "$(stat -c %s "$f")" "$(grep " $(basename "$f")\$" "$D/SHA256SUMS" | cut -d' ' -f1)" \
      "$(tar -tzf "$f" | cut -d/ -f1 | sort -u | wc -l)"
  done; } > "$M"
# 4. merge + commit
(cd "$H" && $PY merge_bakta.py "$H/results/runs_tsv" > "$H/results/merge.log" 2>&1) || { log "FAIL merge_bakta.py (results/merge.log)"; exit 1; }
s=$($PY -c "import json; s=json.load(open('$H/results/summary.json')); t=s['totals']; print(f\"{s['mags']} MAGs, {len(s['missing_or_failed'])} missing/failed; CDS {t['CDS']:,} ({100*s['share_CDS_hypothetical']:.1f}% hypothetical, KO on {t['CDS_with_KO']:,}); tRNA {t['tRNA']:,}, tmRNA {t['tmRNA']:,}, rRNA 5S/16S/23S {t['rRNA_5S']:,}/{t['rRNA_16S']:,}/{t['rRNA_23S']:,}, ncRNA {t['ncRNA']:,}, ncRNA-region {t['ncRNA_region']:,}, CRISPR {t['CRISPR']:,}, sORF {t['sORF']:,}\")")
log "merged: $s"
add=("wetland/bakta/results/bakta_per_mag.tsv" "wetland/bakta/results/bakta_vs_dram2_per_mag.tsv" "wetland/bakta/results/summary.json" "wetland/bakta/results/RESIDUAL_MANIFEST.tsv")
sz=$(stat -c %s "$H/results/bakta_genes_all.tsv.gz"); [ "$sz" -lt 95000000 ] && add+=("wetland/bakta/results/bakta_genes_all.tsv.gz")
for i in 1 2 3 4 5; do
  [ -e "$REPO/.git/index.lock" ] && { sleep 20; continue; }
  git -C "$REPO" add -- "${add[@]}" && git -C "$REPO" commit -q -m "Wetland Bakta: full annotation (coding + non-coding) of the 1,184 DRAM2 representatives"$'\n\n'"$s"$'\n'"Bakta 1.11.4 + DB 6.0-full on poplar; residual outputs (gbff/gff3/json/faa/ffn/...) in bakta/results/residual/ (gitignored, RESIDUAL_MANIFEST.tsv). Collected automatically by bakta/auto_collect_bakta.sh.$TRAILER" -- "${add[@]}" >> "$LOG" 2>&1 && break
  sleep 20
done
for i in 1 2 3 4; do git -C "$REPO" push -q -u origin "$BR" >> "$LOG" 2>&1 && break; sleep $((2 ** i)); done
log "DONE $(git -C "$REPO" log -1 --format=%h)"
