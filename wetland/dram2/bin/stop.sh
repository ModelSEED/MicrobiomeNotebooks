#!/usr/bin/env bash
# Gracefully stop the wetland DRAM2 driver and the running chunk (resumable: rerun the driver).
# Matches by command-line prefix, so it never hits the shell that runs this script.
u=$(id -un)
drv=$(pgrep -u "$u" -f '^(/usr/bin/|/bin/)?bash (/scratch1/afreiburger/wetland_dram2/)?bin/run_all_chunks\.sh')
jvm=$(pgrep -u "$u" -f '^/scratch1/fliu/hub_scratch/chenry/DRAM2/env/env_nf/bin/java .*wetland_dram2/launch/nfcore_offline')
[ -n "$drv" ] && kill -TERM $drv && echo "driver stopped: $drv"
[ -n "$jvm" ] && kill -TERM $jvm && echo "nextflow (SIGTERM, tasks are killed, cache kept): $jvm"
[ -z "$drv$jvm" ] && echo "nothing running"
