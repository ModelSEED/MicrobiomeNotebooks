#!/usr/bin/env python3
"""Write the attempt-1 task scripts a driver config would generate for one sample, without running anything
(migration check: compare the same sample's scripts on poplar and frankensteinsmonster).
usage: python -I gen_task_scripts.py DRIVER_DIR CONFIG SAMPLE OUTDIR
DRIVER_DIR holds driver.py; use a dry-run copy root (dryrun_migration.sh) as the config's root, never a live root.
"""
import importlib.util
import os
import sys

ddir, cfgp, sample, out = sys.argv[1:5]
spec = importlib.util.spec_from_file_location("rodriver", os.path.join(ddir, "driver.py"))
drv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(drv)
cfg = drv.load_cfg(cfgp)
d = drv.Driver(cfg)
os.makedirs(out, exist_ok=True)
n = 0
for name, t in sorted(d.tasks.items()):
    if t.get("sample") != sample:
        continue
    b, _ = drv.body(d, t, 1)
    with open(os.path.join(out, f"{drv.safe(name)}.1.sh"), "w") as fh:
        fh.write(drv.header(cfg, t, 1) + b + drv.FOOTER)
    n += 1
print(f"{n} task scripts for {sample} -> {out}")
