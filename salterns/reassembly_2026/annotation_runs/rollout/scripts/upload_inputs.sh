#!/usr/bin/env bash
# Run LOCALLY (frankensteinsmonster): upload the 2026 comparison inputs to poplar.
#  - comparison/mags_old_vs_new.tsv -> rollout/inputs/mags_old_vs_new.tsv
#  - for every per-sample unit except the two pilot samples (their copies come from the binning test), all 2026 DAS Tool
#    bins work/dastool/<unit>/attempt1-*/out/<unit>_DASTool_bins/*.fa -> rollout/inputs/old/<unit>/, md5-checked on
#    poplar against a local md5 list, then .ok (the gate of the checkm2:old / gunc:old tasks).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
R26=$(cd "$HERE/../../.." && pwd)                       # salterns/reassembly_2026
RR=/scratch1/afreiburger/salterns_2026/rollout
SSH="ssh -o BatchMode=yes"
STAGE=$(mktemp -d)
trap 'rm -rf "$STAGE"' EXIT
scp -o BatchMode=yes -q "$R26/comparison/mags_old_vs_new.tsv" poplar:$RR/inputs/mags_old_vs_new.tsv
n=0
for u in $(awk -F'\t' 'NR>1{print $1}' "$R26/manifest.tsv"); do
  case $u in R1_A_D1|SF2_C_D1) continue;; esac
  d=$(ls -d "$R26"/work/dastool/$u/attempt1-*/out/${u}_DASTool_bins)
  [ "$(ls -d "$R26"/work/dastool/$u/attempt* | wc -l)" = 1 ] || { echo "$u: more than one DAS Tool attempt"; exit 1; }
  mkdir -p "$STAGE/$u"; cp "$d"/*.fa "$STAGE/$u/"
  (cd "$STAGE/$u" && md5sum *.fa > ../$u.md5 && mv ../$u.md5 ./.md5)
  echo "$d" > "$STAGE/$u/.source"
  n=$((n+1))
done
echo "staged $n units"; du -sh "$STAGE"
rsync -a -e "$SSH" "$STAGE/" poplar:$RR/inputs/old/
$SSH poplar "cd $RR/inputs/old && for u in \$(ls); do [ -f \$u/.md5 ] || continue; (cd \$u && md5sum -c --quiet .md5 && n=\$(ls *.fa | wc -l) && [ \"\$n\" = \"\$(wc -l < .md5)\" ] && touch .ok && echo \"\$u \$n bins ok\") || echo \"\$u FAILED\"; done" < /dev/null
