#!/usr/bin/env bash
# Push the downloaded assemblies (kbase_mags.py fetch -> ~/gtdbtk_kbase/assemblies, md5-checked, .part until complete)
# to poplar as they arrive, every 5 min; once the fetch has exited and all 781 are synced (md5 re-checked on poplar), drop
# inputs/ALL_INPUTS_SYNCED so the Bakta driver can finish. Local, detached.
# Usage: setsid nohup bash bakta/sync_assemblies.sh > ~/gtdbtk_kbase/asm_sync.log 2>&1 < /dev/null &
set -uo pipefail
L=$HOME/gtdbtk_kbase/assemblies; R=poplar:/scratch1/afreiburger/wetland_bakta_asm/inputs
while :; do
  rsync -a --include='*.fna' --exclude='*' "$L/" "$R/" || echo "[$(date '+%F %T')] rsync failed; retry"
  n=$(ls "$L"/*.fna 2>/dev/null | wc -l)
  echo "[$(date '+%F %T')] synced $n"
  if ! pgrep -f "kbase_mags.py fetch --list asm_all.tsv" > /dev/null && [ -s "$HOME/gtdbtk_kbase/asm_all_manifest.tsv" ]; then
    rsync -a --include='*.fna' --exclude='*' "$L/" "$R/"
    awk -F'\t' 'NR>1{n=split($4,p,"/"); print $5"  "p[n]}' "$HOME/gtdbtk_kbase/asm_all_manifest.tsv" > /tmp/asm_md5.txt
    if [ "$(wc -l < /tmp/asm_md5.txt)" -eq 781 ] && ssh -o BatchMode=yes poplar 'cd /scratch1/afreiburger/wetland_bakta_asm/inputs && md5sum -c --quiet' < /tmp/asm_md5.txt; then
      ssh -o BatchMode=yes poplar 'touch /scratch1/afreiburger/wetland_bakta_asm/inputs/ALL_INPUTS_SYNCED' < /dev/null
      echo "[$(date '+%F %T')] ALL 781 synced and md5-verified on poplar"; exit 0
    fi
    echo "[$(date '+%F %T')] fetch ended but not all 781 verified; check asm_fetch.log"; exit 1
  fi
  sleep 300
done
