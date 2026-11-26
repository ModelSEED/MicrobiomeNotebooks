#!/bin/bash
# annot57 driver (poplar, detached): eggNOG-mapper 2.1.13 + Bakta 1.12.0 on the 57 modelling species representatives.
# Work dir /scratch1/afreiburger/salterns_2026/annot57 (inputs fasta/<mag_id>.fa and prodigal/<mag_id>.{faa,gff,genes.tsv}
# were copied in and checked beforehand: fasta.md5, protein counts = KOfam per_mag rows).
#
# Thread budget <= 48 at all times:
#   eggNOG-mapper: 2 pooled jobs x --cpu 8 = 16 (Prodigal proteins of the KOfam tracks, groups in scripts/eggnog_group{1,2}.txt)
#   Bakta: 8 workers x --threads 4 = 32 while eggNOG runs; +4 workers (x 4 threads) once both eggNOG jobs have ended.
# Workers claim MAGs with an atomic mkdir in claims/; order = scripts/bakta_order.txt (largest first).
# A Bakta failure writes bakta_FAILED and no worker claims another MAG (stop the tool, report; no tool switching).
# All temp data under $W/tmp (TMPDIR exported; Bakta --tmp-dir, emapper --temp_dir -> DIAMOND --tmpdir; tRNAscan-SE and
# AMRFinderPlus honour TMPDIR). Never /tmp.
set -uo pipefail
W=/scratch1/afreiburger/salterns_2026/annot57
ANNOT=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_
DB=/scratch1/afreiburger/emilykin/processed/refdata/bakta/db
EGG=/scratch1/afreiburger/emilykin/processed/refdata/eggnog
export PATH=$ANNOT/bin:$PATH
export TMPDIR=$W/tmp
cd $W || exit 1
echo $$ > driver.pid
log(){ echo "[$(date -u +%FT%TZ)] $*" >> driver.log; }
log "driver start pid $$ host $(hostname)"

# ---- monitor: /tmp, load, progress every 10 min until ALL_DONE ----
(
  while [ ! -e ALL_DONE ]; do
    {
      echo "== $(date -u +%FT%TZ)"
      uptime
      df -h /tmp | tail -1
      echo "bakta claimed=$(ls claims | wc -l) ok=$(cat logs/bakta_worker_*.progress 2>/dev/null | grep -c '^OK') fail=$(cat logs/bakta_worker_*.progress 2>/dev/null | grep -c '^FAIL')"
      echo "emapper: $(for g in 1 2; do printf 'g%s=%s ' $g "$(grep -h EMAPPER_EXIT logs/emapper_g$g.log 2>/dev/null || echo running)"; done)"
      echo "my tmp entries: $(ls -A tmp | wc -l); my /tmp files: $(find /tmp -maxdepth 2 -user afreiburger 2>/dev/null | wc -l)"
    } >> monitor.log
    sleep 600
  done
) &
MON=$!

# ---- versions ----
bakta --version > logs/bakta.version 2>&1
emapper.py --data_dir $EGG --version > logs/emapper.version 2>&1
cp $DB/version.json logs/bakta_db_version.json

# ---- eggNOG-mapper: 2 pooled jobs on the KOfam-track Prodigal proteins (trailing '*' stripped) ----
EPIDS=()
for g in 1 2; do
  mkdir -p eggnog/g$g tmp/emapper_g$g
  : > eggnog/g$g/input.faa
  while read m; do sed 's/\*$//' prodigal/$m.faa >> eggnog/g$g/input.faa; done < scripts/eggnog_group$g.txt
  log "emapper g$g: $(wc -l < scripts/eggnog_group$g.txt) MAGs, $(grep -c '>' eggnog/g$g/input.faa) proteins"
  ( /usr/bin/time -v -o logs/emapper_g$g.time nice -n 10 emapper.py --data_dir $EGG -m diamond --itype proteins \
      -i eggnog/g$g/input.faa -o annot57_g$g --output_dir eggnog/g$g --temp_dir $W/tmp/emapper_g$g --cpu 8 \
      --dmnd_iterate no --dbmem --override < /dev/null > logs/emapper_g$g.log 2>&1
    echo "EMAPPER_EXIT $?" >> logs/emapper_g$g.log ) &
  EPIDS+=($!)
done

# ---- Bakta workers ----
bakta_worker(){
  local wid=$1
  while read m; do
    [ -e bakta_FAILED ] && { echo "STOP worker $wid: bakta_FAILED exists"; break; }
    mkdir claims/$m 2>/dev/null || continue
    if [ -s bakta/$m/$m.json ]; then echo "SKIP $m"; continue; fi
    rm -rf bakta/$m tmp/bakta_$m; mkdir -p tmp/bakta_$m
    echo "START $m $(date -u +%FT%TZ)"
    if /usr/bin/time -v -o logs/bakta_$m.time nice -n 10 bakta --db $DB --threads 4 --skip-plot --keep-contig-headers \
         --tmp-dir $W/tmp/bakta_$m --output bakta/$m --prefix $m fasta/$m.fa < /dev/null > logs/bakta_$m.log 2>&1 \
       && [ -s bakta/$m/$m.json ]; then
      echo "OK $m $(date -u +%FT%TZ)"
    else
      echo "FAIL $m $(date -u +%FT%TZ)"; echo "$m $(date -u +%FT%TZ)" >> bakta_FAILED
    fi
    rm -rf tmp/bakta_$m
  done < scripts/bakta_order.txt
  echo "WORKER_END $wid $(date -u +%FT%TZ)"
}
BPIDS=()
for i in 1 2 3 4 5 6 7 8; do
  bakta_worker $i > logs/bakta_worker_$i.progress 2>&1 &
  BPIDS+=($!)
  sleep 20
done
log "bakta: 8 workers started (pids ${BPIDS[*]}); emapper pids ${EPIDS[*]}"

# ---- after eggNOG: hand its 16 threads to 4 more Bakta workers ----
for p in "${EPIDS[@]}"; do wait $p; done
log "emapper done: $(grep -h EMAPPER_EXIT logs/emapper_g1.log logs/emapper_g2.log | tr '\n' ' ')"
if grep -q 'EMAPPER_EXIT 0' logs/emapper_g1.log && grep -q 'EMAPPER_EXIT 0' logs/emapper_g2.log; then touch EGGNOG_DONE; else touch EGGNOG_FAILED; fi
for i in 9 10 11 12; do
  bakta_worker $i > logs/bakta_worker_$i.progress 2>&1 &
  BPIDS+=($!)
done
log "bakta: 4 more workers started (total pids ${BPIDS[*]})"
for p in "${BPIDS[@]}"; do wait $p; done
nok=$(cat logs/bakta_worker_*.progress | grep -c '^OK')
log "bakta done: OK=$nok of $(wc -l < scripts/bakta_order.txt)"
[ "$nok" -eq "$(wc -l < scripts/bakta_order.txt)" ] && touch BAKTA_DONE
touch ALL_DONE
log "ALL_DONE"
{ echo "== final $(date -u +%FT%TZ)"; uptime; df -h /tmp | tail -1; } >> monitor.log
