#!/usr/bin/env bash
# READ-ONLY inventory of the EmilyKin GTDB-Tk R220 run on poplar (nothing is written on poplar).
# Needs the ControlMaster session from the user's Duo login (memory: poplar-annotation-envs):
#     tmux new -s poplar-auth "ssh -fN -o ClearAllForwardings=yes poplar"     # user types password + Duo once
# Saves to ../poplar_inventory/: config.yaml, reprocess_291 env.sh, mag/ listing, gtdbtk.log heads (GTDB-Tk version +
# the 'not using the reference data intended for this release' warning), the server-native R220 summaries (the repo copy
# is a best-of-both hybrid), poplar's MAG md5s (compared with ../local_catalogue_md5.txt) and assign_mag_iterativeids.py.
# Never runs du/find on /scratch1 (hangs on poplar).
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
OUT="$HERE/../poplar_inventory"; mkdir -p "$OUT"
S=(ssh -o BatchMode=yes -o ClearAllForwardings=yes -o ConnectTimeout=20 poplar)
"${S[@]}" true || { echo "poplar not reachable (log in once via the tmux command above)"; exit 2; }
P=/scratch1/afreiburger/emilykin
PROC=$P/processed
"${S[@]}" "cat $PROC/config.yaml"                         > "$OUT/config.yaml"            2>&1
"${S[@]}" "cat $P/EmilyKin/reprocess_291/env.sh"          > "$OUT/reprocess_291_env.sh"   2>&1
"${S[@]}" "ls -la $PROC/mag/; for d in $PROC/mag/gtdbtk*; do echo == \$d; ls -la \$d \$d/classify 2>&1; done" > "$OUT/mag_dirs.txt" 2>&1
for d in gtdbtk_291 gtdbtk; do
  "${S[@]}" "head -4 $PROC/mag/$d/gtdbtk.log; grep -m5 -i 'reference data\|data version\|r232\|r220' $PROC/mag/$d/gtdbtk.log" > "$OUT/${d}.log_head.txt" 2>&1
  for dom in bac120 ar53; do
    "${S[@]}" "cat $PROC/mag/$d/classify/gtdbtk.$dom.summary.tsv" > "$OUT/${d}.$dom.summary.tsv" 2>/dev/null || rm -f "$OUT/${d}.$dom.summary.tsv"
  done
done
"${S[@]}" "source $P/EmilyKin/reprocess_291/env.sh; \$ENV_G/bin/gtdbtk --version; cat \$GTDBTK_DATA_PATH/metadata/metadata.txt" > "$OUT/gtdbtk_env_version.txt" 2>&1
"${S[@]}" "cd $P/EmilyKin/dereplicated_genomes && md5sum *.fa"                      > "$OUT/poplar_catalogue_md5.txt" 2>&1
"${S[@]}" "cd $P/EmilyKin && git log -1 --format='%h %ad %s' -- dereplicated_genomes && git status --short dereplicated_genomes | head -20" > "$OUT/poplar_repo_state.txt" 2>&1
"${S[@]}" "cat $PROC/mag/abundance/assign_mag_iterativeids.py"                    > "$OUT/assign_mag_iterativeids.py" 2>&1
# compare poplar's catalogue with the local best-of-both one (committed in EmilyKin at 6541627)
n_diff=$(join <(awk '{print $2, $1}' "$OUT/poplar_catalogue_md5.txt" | sort) <(awk '{print $2, $1}' "$HERE/../local_catalogue_md5.txt" | sort) | awk '$2 != $3' | wc -l)
n_pop=$(wc -l < "$OUT/poplar_catalogue_md5.txt")
echo "poplar MAG files: $n_pop; differing from the local best-of-both catalogue: $n_diff" | tee "$OUT/catalogue_compare.txt"
grep -h "GTDB-Tk v\|reference data" "$OUT"/*.log_head.txt | sort | uniq -c
