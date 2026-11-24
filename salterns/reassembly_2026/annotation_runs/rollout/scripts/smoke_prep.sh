#!/usr/bin/env bash
# Prepare the smoke-test root: the 2026 DAS Tool bins used as SMK_A's "old" set (a copy of the 6 SF2_C_D1 bins that the
# binning test already holds); SMK_B deliberately gets none (tests the 'never uploaded' path). Run once on poplar.
set -euo pipefail
S=/scratch1/afreiburger/salterns_2026/rollout/smoke
BT=/scratch1/afreiburger/salterns_2026/binning_test
mkdir -p $S/inputs/old $S/logs $S/scripts
rm -rf $S/inputs/old/SMK_A; cp -r $BT/inputs/old/SF2_C_D1 $S/inputs/old/SMK_A; touch $S/inputs/old/SMK_A/.ok
cp /scratch1/afreiburger/salterns_2026/rollout/scripts/{driver.py,rolib.py,btlib.py,status.sh} $S/scripts/
ls -la $S/inputs/old/SMK_A | head
