# Brief: assemble and bin the GROWdb raw metagenomes (separate Claude session)

The user asked (2026-11-09): "automatically also start/queue the same assembling/binning pipeline that is in progress for
the salterns project in a separate claude session". You are that session. Work in `grow/assembly_binning_2026/`.

## Inputs (already in place, owned by another session; read only)
- **The reads.** KBase narrative **ws 275929** (`freiburgermsu:narrative_1794152940209`) is being filled with 171 interleaved
  PairedEndLibraries `<sample>.<SRR>.reads`, one per SRA run, from 133 river samples. They are the metagenomes behind the 3,825
  GROWdb MAGs (ws 145226). `grow/raw_metagenomes_sra/README.md` has the provenance.
- **Run list:** `grow/raw_metagenomes_sra/runs.tsv` gives sample, run, spots, bases and source. 46 samples have 2–3 runs.
- **Import status:** `grow/raw_metagenomes_sra/kbase_import_ledger.json` is written atomically by the import watcher
  (`kbase_import.py watch`, PID in `ps`). A run is imported once `runs[SRR].state == "done"`; the object ref is in
  `runs[SRR].ref`. As of 2026-11-09 18:45 UTC, 57 of 171 were done and the rest are expected over the next 1–2 days.
  - **Never** write that ledger, and never stop that watcher.
  - Once all imports finish, it files a request to add ws 275929 to the KBase Org `grow`.
- **Local FASTQ copies** of every run (ENA, MD5-verified) are in `/mnt/cold/Documents/GROW_raw_metagenomes/<sample>/`. Use them if a
  local step needs reads.

## What to build: the salterns pipeline, adapted
The template is `salterns/reassembly_2026/pipeline.py` with its `README.md`, plus the memory note `salterns-reassembly-2026`.
That pipeline is running live in another watcher.
- **Do not modify it**, its ledger, or its watcher (PID from `pgrep -af "[p]ipeline.py watch"`).
- Copy what you need into a new driver here, with its own ledger, its own lock and its own watcher.

Keep the same apps, versions, settings and safeguards:
- **Steps:** import is already done, then Trimmomatic, then metaSPAdes per sample (falling back to MEGAHIT on OOM, disk errors or repeated failure), then MetaBAT2, MaxBin2 and CONCOCT.
- **Refinement and QC:** local DAS Tool 1.1.7 (`~/Documents/dastool_env`; KBase's 1.1.2 has a scoring bug), then extract bins, then KBase CheckM (beta, reduced tree), then keep bins ≥50% complete and <10% contaminated.
- **Taxonomy:** local GTDB-Tk 2.8.0 / R232, using the same environment, scope and mem_guard pattern as salterns.
- **Safeguards:** write-ahead ledger, recovery of lost submits, pilots, output validation and retry budgets.

## GROW-specific decisions
The user confirmed decisions 2, 3 and 5 on 2026-11-09: no co-assembly, per-sample coverage, and GROW stays behind salterns in the KBase queue.
1. **Unit = sample (133).** For multi-run samples, merge the trimmed runs with `kb_ReadsUtilities` merge (validate that read and base counts equal the sums), then assemble the merged library. Start a sample as soon as all of its runs are `done` in the import ledger; don't wait for all 171.
2. **No co-assembly.** The samples are independent rivers, not replicate cores; the salterns site co-assembly does not apply.
3. **Binning coverage: per-sample** (each sample's own library), as GROWdb itself did.
   - Multi-sample coverage across 133 libraries would be far heavier.
   - Treat it as a possible later improvement and note it in the README.
4. **Outputs go to a separate narrative: ws 276048** (`freiburgermsu:narrative_1794250260895`, "GROWdb metagenomes: assembly and binning (2026)").
   - The user decided this on 2026-11-09 to keep the reads narrative from growing too large.
   - Read the libraries from ws 275929 by ref, and write every output (trimmed and merged reads, assemblies, BinnedContigs, extracted bins, CheckM reports) to ws 276048. KBase jobs here run with ws 276048 as their workspace.
   - Name objects as salterns does, with the sample as the unit.
5. **Share KBase capacity with salterns; salterns has priority.**
   - Salterns is still in metaSPAdes and binning, and the KBase extreme pool is about 12 slots, the bottleneck there.
   - Cap this driver at **extreme 3, bigmem 6, njs 6** while the salterns ledger still has assemble or binning work pending or live.
   - Raise the caps only once salterns has cleared those stages.
   - Count live jobs from both ledgers when deciding.
6. **Local resources.**
   - DAS Tool workers: at most 2 here.
   - **Never run GTDB-Tk while the salterns GTDB-Tk scope (`gtdbtk-salterns-reassembly`) is active.** Use your own scope name, e.g. `gtdbtk-grow-assembly`.
   - Check free disk before local steps.
7. **Pilot before fan-out.** Take one small single-run sample and the largest sample (by total bases, about 50 Gbp or more after merging) through trim → assemble → bin. Fan out once both work.
   - The large samples are 2–5× the salterns ones, so expect metaSPAdes OOM and rely on the MEGAHIT fallback.

## Operating rules (learned the hard way; see memory)
- **Processes:**
  - Never use `pkill -f` or `pgrep -f` with a pattern that also appears in your own command line; it kills your own shell. Use the bracket trick or PIDs.
  - Start the watcher detached in its own systemd user scope, with a new name each start (`grow-assembly-watch-<utc stamp>`) and `OOMPolicy=continue`, as the salterns README shows. Restart it only by PID, while it sleeps.
- **After a reboot** (the machine rebooted 2026-11-09 11:55 CST): /tmp scratchpads are lost, so keep everything in the repo. Check that the watchers are alive and restart them per README.
- **Code:** KBUtilLib is imported from `~/Documents/KBUtilLib/src` and is not pip-installed. The token is in `~/.kbase/token`; never print it. Only use `~/Documents/py_venv/bin/python`.
- **Git:** don't commit or push unless the user asks.

## Deliverables
- `grow/assembly_binning_2026/`: the driver, `README.md` (decisions, pools, how to run and restart), `manifest.tsv`, and the ledger.
- The pilot, then a running watcher that keeps advancing automatically as imports land, through to GTDB-Tk.
- Add a section to the memory file `grow-mucc-data-inventory.md`, and a pointer line in `MEMORY.md` if you create a new memory.
- Ask the user only about decisions that genuinely block you. Otherwise proceed with the defaults above and report what you did.

## Update 2026-11-10 (user decision): MEGAHIT for every sample
Assemble every sample with MEGAHIT (release 1.2.9, meta-sensitive, bigmem); metaSPAdes is not used.
- **Changes:** cancel the eastriver metaSPAdes job and re-assemble deadrun with MEGAHIT. The metaSPAdes-based deadrun outputs are superseded.
- **Reasons:**
  - The river samples are high-diversity with low per-genome coverage.
  - Memory: 13 samples exceed 50 Gbp.
  - MEGAHIT runs on bigmem, which frees the extreme pool for MaxBin2 and CONCOCT.
  - One assembler is used throughout, and GROWdb itself used mostly MEGAHIT.

## Update 2026-11-12 (user decisions)
1. **KBase errors as transient:** an EE2 "Unexpected Job Error", including kills by `monitor_jobrunner_logs.py`, is treated as a transient KBase outage. It does not count toward retries.
   - The 3 tasks that gave up after the outages of Nov 11 and 12 are reset.
2. **Pilot gate relaxed:** assembly and binning fan out now to every sample under 50 Gbp. Only the 13 samples of 50 Gbp or more wait for the eastriver pilot.

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
- **2026-11-25 (user):**
  - **Poplar:** runs multiple assemblies concurrently, under dynamic limits: ≥30% idle, ≤96 threads, ≥400 GB MemAvailable, ≥200 GB /scratch free.
  - **FM:** the assembly threshold is raised to ~16 GB gz, one large FM assembly at a time, with the watchdog cap at 120 GB.
  - **Routing:** each assembly goes to the first available route of FM, poplar or KBase.
- **2026-11-25 (user): more work to poplar.** Up to 5 concurrent MEGAHIT jobs (≤160 threads) at nice 10, sharing CPU equally with the user's EmilyKin work.
  - **New start rule:** start a job when other users' normal-priority load leaves room and load < ~500, replacing the 30%-idle rule.
  - **Disk:** staging and outputs on /scratch1 (8.1 TB free), keeping the /scratch ≥200 GB floor.
  - **Routing:** queued KBase assemblies are taken over onto poplar, largest first.
