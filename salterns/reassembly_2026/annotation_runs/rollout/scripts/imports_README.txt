imports/  (R2A migration, 2026-11-16; full record: README.md "R2A migrated to frankensteinsmonster" in the
MicrobiomeNotebooks repo, salterns/reassembly_2026/annotation_runs/rollout/)

R2A_A_D1 R2A_A_D2 R2A_B_D1 R2A_B_D2 R2A_C_D1 R2A_C_D2 R2A_B_H2O are assembled, binned and QC'd on frankensteinsmonster
(~/Documents/salterns_rollout_local/rollout) by this same driver.py with the same commands (paths and threads differ;
SPAdes -m 120). R2A_A_H2O and R2A_C_H2O (> 9 Gbp, 200g heap) run here as ordinary samples.

When a sample finishes there, its push:<s> task copies <s>/, qc/new/<s>/ and mags_kept/<s>__*.fa into
imports/<s>.part/ with MANIFEST.md5 and PUSHED.json, runs md5sum -c here, renames the directory to imports/<s>/ and
writes imports/<s>/READY last. The driver here then runs import:<s> (held as 'import_wait' until READY exists): md5sum -c
again, copy into <s>/ and qc/new/<s>/, verify-import (one CheckM2/GUNC/CheckM1 row per bin), and re-apply the KEEP
rule, which must give the same kept.txt and byte-identical kept FASTAs as the remote host; then keep:<s> runs here as for
every sample, and GTDB-Tk / report / DONE wait for all of it.

If the remote host gives up on a sample, write imports/<s>/FAILED (one line why): the driver marks import:<s> failed,
and the run ends FAILED for that sample instead of waiting forever. To run such a sample here instead, see README.md.
Do not edit files under imports/<s>/ by hand: the import verifies them against MANIFEST.md5.
