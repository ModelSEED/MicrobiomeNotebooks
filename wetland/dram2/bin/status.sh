#!/usr/bin/env bash
# Progress report for the wetland DRAM2 run.   bin/status.sh
B=/scratch1/afreiburger/wetland_dram2
cd "$B" || exit 1
echo "== driver (logs/driver.log)"; tail -n 4 logs/driver.log 2>/dev/null
pgrep -u afreiburger -f 'run_all_chunks.sh' >/dev/null && echo "driver: RUNNING" || echo "driver: not running"
echo "== chunks"
for d in inputs/chunk_*/ inputs/test3/; do
  c=$(basename "$d"); n=$(ls "$d" | wc -l); st="pending"
  [ -e "runs/$c/FAILED" ] && st="FAILED($(cat runs/$c/FAILED))"
  [ -e "runs/$c/DONE" ] && st="DONE $(cat runs/$c/DONE)"
  t=$(ls -t runs/$c/out/pipeline_info/execution_trace_*.txt 2>/dev/null | head -1)
  if [ -n "$t" ] && [ ! -e "runs/$c/DONE" ]; then
    st="running"
    for p in CALL_GENES TRNA_SCAN HMM_SEARCH_KOFAM HMM_SEARCH_VOG MMSEQS_SEARCH_MEROPS KOFAM_HMM_FORMATTER COMBINE_ANNOTATIONS DISTILL PRODUCT_HEATMAP; do
      k=$(awk -F'\t' -v p=":$p" 'index($4,p) && ($5=="COMPLETED"||$5=="CACHED")' "$t" | wc -l)
      st="$st $p=$k"
    done
  fi
  printf '%-10s %4s MAGs  %s\n' "$c" "$n" "$st"
done
echo "== nextflow processes"; pgrep -u afreiburger -fa 'wetland_dram2/launch/nfcore_offline' | grep -o 'runs/[a-z0-9_]*/out' | sort -u
echo "== CPU used by the DRAM2 Nextflow process trees (5 s sample, cores)"
roots=$(pgrep -u afreiburger -f 'wetland_dram2/launch/nfcore_offline' | tr '\n' ' ')
if [ -n "$roots" ]; then
  pids=$(ps -e -o pid=,ppid= | awk -v r="$roots" 'BEGIN{n=split(r,a," "); for(i=1;i<=n;i++) keep[a[i]]=1}
         {par[$1]=$2; ids[NR]=$1} END{ch=1; while(ch){ch=0; for(i in ids){p=ids[i]; if(!(p in keep) && (par[p] in keep)){keep[p]=1; ch=1}}} for(p in keep) print p}')
  snap(){ for p in $pids; do awk '{print $14+$15}' /proc/$p/stat 2>/dev/null; done | awk '{s+=$1} END{print s+0}'; }
  a=$(snap); sleep 5; b=$(snap)
  echo "$(( (b-a) / (5*$(getconf CLK_TCK)) )) cores (approx; tasks starting/ending in the window are missed)"
else
  echo "no DRAM2 nextflow running"
fi
uptime
