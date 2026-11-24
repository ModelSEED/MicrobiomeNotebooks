#!/usr/bin/env bash
# Quick status: bash /scratch1/afreiburger/salterns_2026/dram/scripts/status.sh
W=/scratch1/afreiburger/salterns_2026/dram
L=$W/logs
echo "markers: raw=$(ls /scratch1/afreiburger/salterns_2026/dram_db/raw/FETCH_OK 2>/dev/null | wc -l) $(ls $W/state 2>/dev/null | tr '\n' ' ')"
[ -f $W/state/FAILED ] && echo "FAILED: $(cat $W/state/FAILED)"
[ -f $W/driver.pid ] && { p=$(cat $W/driver.pid); ps -p $p > /dev/null && echo "driver $p running" || echo "driver $p not running"; }
nd=0; nf=0
[ -f $L/run_all.log ] && { nd=$(grep -c ' done ' $L/run_all.log); nf=$(grep -c FAILED $L/run_all.log); }
echo "MAGs done: $nd  FAILED lines: $nf  of 558; DRAM annotate processes: $(pgrep -u $(id -un) -f 'DRAM.py annotate' | wc -l)"
tail -3 $L/driver.log 2>/dev/null
[ -f $L/setup_dbs.pfam.log ] && grep -E '^20' $L/setup_dbs.pfam.log | tail -2
uptime; free -g | sed -n 2p
