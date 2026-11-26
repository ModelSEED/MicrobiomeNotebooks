# Brief: GROW site multi-sample binning + co-assembly pilots (salterns-style, parallel effort)

The user approved this on 2026-11-18. It is a **parallel** effort to the per-sample GROW run (`grow/assembly_binning_2026/`, watcher PID in `ps`, cwd there).
- It reuses that run's trimmed and merged reads and its per-sample MEGAHIT assemblies.
- It re-does only binning (with site multi-sample coverage) and adds co-assembly pilots.
- All outputs are kept **separate**: from the per-sample run's MAGs, and from GROWdb's 3,825/2,093 MAGs. Never merge or dereplicate them together.

## Target and inputs
- **Output narrative: ws 277632** (`freiburgermsu:narrative_1795046239204`, "GROWdb metagenomes: site multi-sample binning and co-assembly (2026)"). Everything this effort makes is saved there.
- **Inputs, by ref from ws 276048:** the per-sample trimmed/merged reads library and the `<sample>.MEGAHIT.Assembly` of each sample. Take them from `grow/assembly_binning_2026/ledger.json`, read only.
  - Never write that ledger, take its lock, or touch its watcher.
- **Groups:** `site_groups.tsv` (this folder). `site_group` is the grouping; empty means unpaired. `coassembly_pilot` marks the 8 co-assembly pilot groups.
  - **25 groups, 67 samples:**
    - 20 NASQAN/USGS station time series. One station per group, NOT river-wide; the Mississippi/Missouri stations stay separate.
    - 5 WHONDRS same-site pairs: Red Butte Creek, Columbia 300 area, Columbia River 1–2, White Clay Creek, Logan River.
  - **The 66 unpaired samples are out of scope:** no re-binning, no co-assembly, no coverage borrowed from other samples. Their per-sample results stand as they are (user decision).

## Workflow: the salterns pipeline (`salterns/reassembly_2026/pipeline.py`), adapted
Same apps, versions and parameters as salterns and the GROW per-sample run, unless noted.
1. **sitecov binning, for each of the 67 grouped samples:**
   - Run MetaBAT2 (beta, bigmem), MaxBin2 (release, extreme) and CONCOCT (release, extreme) on the sample's existing MEGAHIT assembly, with `reads_list` = the trimmed/merged libraries of **every sample in its site group** (salterns `coverage_group`). Groups have 2–5 samples.
   - Then local DAS Tool 1.1.7, the KBase extract beta, and CheckM beta with the reduced tree. Collect bins ≥50% complete and <10% contaminated.
   - Names: `<sample>.sitecov.{MetaBAT2,MaxBin2,CONCOCT,DASTool}.BinnedContigs` and `<sample>.sitecov.DASTool_bins.AssemblySet`.
   - Extracted bins get `<bid>__<sample>.sitecov` (or another name that cannot collide).
2. **Co-assembly pilots, only for the 8 `coassembly_pilot` groups:**
   - Groups: the 5 WHONDRS pairs, plus USGS 03609750 Tennessee (2 samples, 7 weeks apart), 03378500 Wabash (3, within 2011) and 04264331 St. Lawrence (4, 2014–2016). They are 45–71 Gbp each.
   - Steps:
     - Merge the group's per-sample libraries with `kb_ReadsUtilities` merge, validating that read and base counts equal the sums.
     - Co-assemble with MEGAHIT/run_megahit (release 1.2.9) **meta-large**, using the salterns co-assembly settings and its retry ladder (min-count 2 with k 27..127, then defaults).
     - Bin the co-assembly with group coverage, using the same chain as step 1.
   - Names: `<group>.MEGAHIT.coassembly.Assembly`, `<group>.coassembly.*`.
   - The NASQAN pilots test co-assembly across time; strain drift may fragment them. The WHONDRS same-day pairs are the closest thing GROW has to salterns' replicate cores.
3. **GTDB-Tk:** local 2.8.0 / R232, as its own batch with separate output directories for the sitecov MAGs and the co-assembly MAGs.
   - Never run it while another GTDB-Tk scope is active (GROW per-sample or salterns). Use your own scope name.

## Scheduling: this effort queues behind the per-sample GROW run
This is the same arrangement the per-sample run had behind salterns.
- **While the per-sample run still has assemble or binner work pending or live:** cap this driver at **extreme 3, bigmem 6, njs 6**.
- **After that:** use the totals (12/20/15) minus the per-sample run's live jobs.
- Count live jobs from both ledgers.
- **Local:** at most 1 DAS Tool worker while the per-sample run's DAS Tool is active, 2 after. Local steps need ≥100 GB free.
- **Start conditions:** a group starts once every member's MEGAHIT assembly and reads are done in the per-sample ledger.
  - Run one small NASQAN pair through sitecov first as a pilot, then fan out.
  - The co-assembly pilots can go once their members' reads are done; they need no per-sample assemblies.
- **Transient errors:** EE2 "Unexpected Job Error" / `monitor_jobrunner_logs` kills are transient (user decision of 2026-11-12). They don't count toward retries and get a 72-h budget.

## Deliverables
- **Code and records:** a new driver in this folder, with its own ledger, lock, watcher and systemd scope (`grow-sitecov-watch-<utc stamp>`), plus `README.md` and `manifest.tsv`. Copy the needed parts of the per-sample driver; don't modify that driver.
- **Comparison table** at the end, on the 67 grouped samples: MAGs ≥50/<10 per sample for per-sample vs sitecov. For the 8 groups: co-assembly MAGs vs the union of their members' sitecov MAGs. Also report completeness, contamination and N50 distributions.
- **Memory:** a section in `grow-mucc-data-inventory.md`.
- **Operating rules:** the same as `grow/assembly_binning_2026/BRIEF.md`. No `pkill -f` with self-matching patterns; restart by PID while sleeping; new scope name each start; keep everything in the repo (/tmp is lost on reboot); never print the token; no commits unless asked.

## Update 2026-11-23 (user decision): hybrid local/KBase mode
KBase queue waits reached 35–48 h. The user approved moving the stalled work (remaining assemblies, binners, extract, CheckM) to this machine.
- **How:** same tools, versions and parameters; a new micromamba env, made like ~/Documents/dastool_env.
- **No double runs:** the KBase jobs taken over are cancelled.
- **Results:** every result is pushed back by direct Workspace save into the SAME narrative (per-sample → 276048, sitecov → 277632), with provenance saying the step ran locally. CheckM results go on the bins as metadata.
- **Local limits:** ≤48 threads and ≤160 GB RAM, ≥100 GB free disk, never alongside GTDB-Tk. The per-sample run keeps priority.
- **TF as a second node (also 2026-11-23):** `ssh TF` gives 24 cores, 91 GB RAM, 722 GB disk, but only ~5 MB/s to and from it.
  - **TF's work:** steps with small inputs and heavy compute: CheckM; MaxBin2/CONCOCT/MetaBAT2 from coverage tables computed here; MEGAHIT for samples under ~15 Gbp.
  - **Stays on FM:** read mapping, large MEGAHIT runs and every KBase push. No KBase token on TF.
  - **Limits on TF:** ≤20 threads, ≤60 GB RAM.
- **Correction, same day: poplar REPLACES TF as the compute node.** TF is not used.
  - **Method stays the same:** our tools and versions (MEGAHIT 1.2.9, MetaBAT2+MaxBin2+CONCOCT, DAS Tool 1.1.7, CheckM 1.1.2), installed in new envs under /scratch/afreiburger using metag-hybrid's mamba. metag-hybrid's own MetaBAT2-only/CheckM2 path is NOT used, so the 133 samples stay consistent.
  - **Poplar limits:** ≤64 threads (shared, load ~500); keep ≥200 GB free on /scratch (622 GB free); pull inputs straight from KBase.
  - **Token:** passed by stdin into the process environment, never written to a file. FM stays the coordinator.
- **FINAL compute plan (user, 2026-11-23; supersedes the poplar redirect):** FM + TF by strength, poplar as overflow only.
  - **FM:** trim, mapping and coverage, all MEGAHIT, MetaBAT2, DAS Tool, extract, and every KBase push. Cap ≈48 threads.
  - **TF:** CheckM, and MaxBin2/CONCOCT from coverage tables; no reads moved there. Cap ≈20 threads at nice ≥10, with no KBase token.
  - **poplar:** overflow only, mainly co-assemblies or OOM retries. Cap ≈32 threads, keep ≥200 GB free, pinned envs, metag-hybrid left untouched.
- **2026-11-24 (user, option b):** sitecov gets a reserved local share: 12 of 48 FM threads (~40 GB), plus a fair share of TF for CheckM, MaxBin2 and CONCOCT. Per-sample keeps 36 threads (~120 GB) and first claim on idle capacity.
- **2026-11-25 (user):** providencecreek MEGAHIT on FM went over memory (161 GB) and was stopped.
  - **Poplar:** the overflow runner is built now. providencecreek is resumed there from its checkpoints (or rerun if a resume is unsafe). Large libraries and the co-assemblies go to poplar.
  - **FM:** keeps the full pipeline but assembles only small samples, under a conservative size threshold, with a memory watchdog that requeues a job to poplar.
  - **Pilot:** one small FM assembly is started now as the end-to-end pilot.
