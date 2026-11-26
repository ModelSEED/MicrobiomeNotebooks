# KBase RASTtk route: cancelled 2026-11-20 ~19:50 UTC (Andrew's decision)

**Why:** KBase's njs and bigmem queues were backed up, with other jobs waiting 15–28 h to start. All 57 genomes would have taken about 5–6 days. Andrew chose instead to annotate the 57 modelling genomes through the public RAST annotation service (ModelSEEDpy `RastClient`; RAST stages annotate_proteins_kmer_v2, annotate_proteins_kmer_v1 and annotate_proteins_similarity) on the Prodigal genes the KOfam calls use, and to cancel the KBase route.

**Done to cancel:**
- The watcher (PID 2020663, scope `salterns-m26-rast-watch-20261120t192259z`) was stopped with SIGTERM while sleeping between checks (wchan `hrtimer_nanosleep`).
- EE2 job `6b009653b3a2d40340bd94aa` (AssemblyUtil `save_assemblies_from_fastas`, the 2 pilot genomes) went from `queued` to `terminated` through `execution_engine2.cancel_job`.

**State left behind:**
- No `m26_*` object was ever saved in ws 275645.
- The 2 pilot FASTAs remain in the KBase Blobstore. They are unreferenced and harmless.
- `ledger.json` still shows those 2 genomes as `uploading`. It is not used any more; do not resume the watcher.
