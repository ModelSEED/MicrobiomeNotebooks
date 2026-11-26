# hybrid/: GROW local/KBase hybrid mode (user-approved 2026-11-23)

KBase queue waits reached 35–48 h, even for small extract and CheckM jobs. The user approved moving the stalled
work of both GROW drivers off KBase:
- the per-sample run, `../assembly_binning_2026`, ws 276048;
- the site multi-sample / co-assembly run, `../multisample_binning_2026`, ws 277632.

The steps run on our own machines, and every result is saved back into the SAME narrative, so each narrative holds all of
its results. The decision is recorded in both `BRIEF.md` files (Update 2026-11-23).

## What runs where
| Step | Before | Now |
|---|---|---|
| MEGAHIT (per-sample meta-sensitive; co-assembly meta-large) | KBase MEGAHIT app | first route that can take it (user-confirmed 2026-11-25): **FM** for per-sample libraries ≤ 16 GB gz (one at a time), **poplar** (concurrent jobs) for any assembly, else a **KBase** job as the fallback queue; never two routes at once |
| MetaBAT2, MaxBin2, CONCOCT (mapping + binning) | KBase apps | first route that can take it (user-confirmed 2026-11-25): **FM** (within the run's share), then **TF** for MaxBin2 / CONCOCT (mapping + binning there), else a **KBase** job as the fallback queue (taken over while queued); never two routes at once |
| DAS Tool 1.1.7 | already local (FM) | unchanged |
| extract + CheckM | KBase MetagenomeUtils + kb_Msuite | **FM**, with the CheckM part on **TF** when FM is full |
| GTDB-Tk 2.8.0 / R232 | already local (FM) | unchanged |
| trim, merge | KBase | all done in KBase before the switch, so nothing to move |

- **Machines:**
  - **FM:** this machine. It coordinates, keeps the ledgers and makes every KBase push.
  - **TF:** `ssh TF`. It only computes, and never gets the KBase token. It runs CheckM, and since 2026-11-25 also MaxBin2 and CONCOCT.
  - **poplar:** overflow for MEGAHIT (user-confirmed 2026-11-25): large per-sample libraries, the sitecov co-assemblies, and steps stopped by FM's watchdog. Up to 5 concurrent jobs of 32 threads (≤ 160 of ours, user-confirmed 2026-11-25, third decision), `nice 10` (the same priority as the user's own EmilyKin work there), under `/scratch/afreiburger/grow_hybrid` (metag-hybrid untouched). See "poplar overflow" below.
- **KBase jobs:**
  - **Queued** jobs are cancelled when a local slot is free for that unit, and the step starts locally in the same cycle. The ledger records the cancel (`status: terminated, takeover: true`), and it is not counted against the unit.
  - **Running** jobs are left to finish and are adopted as before.
  - No new KBase jobs are submitted for these stages, except MEGAHIT and the three binners as the fallback queue: an assembly that neither FM nor poplar can take in this cycle gets a KBase job (user-confirmed 2026-11-25). While that job is still queued, FM or poplar takes it over (cancel, same cycle) as soon as one has room. A pinned task (`route` in the ledger: the FM pilot, or an overflow that resumes on poplar) never goes to KBase.

## Binners on TF (user-confirmed 2026-11-25)
- **Why:** on 2026-11-25 no binner ran anywhere for hours. FM's per-sample share was taken by the weststlouiscreek MEGAHIT plus DAS Tool / extract, and sitecov's 12 / 40 by one CONCOCT. FM was also loaded past its 64 cores by other users' jobs, and the binners had no other route.
- **Routing (`glue.try_run`):** each binner step goes to the first route that can take it:
  - **FM:** within the run's share (16 threads; 12 for sitecov).
  - **TF:** MaxBin2 and CONCOCT.
  - **KBase:** a job as the fallback queue (MetaBAT2 too). It is taken over (cancel confirmed) as soon as FM or TF has room while it's still queued.
- **What runs on TF:**
  - **MaxBin2:** the whole step runs on TF: bowtie2 2.2.3 `-U` per library into MaxBin's getsam through a FIFO, then `run_MaxBin.pl -abund_list` (`hw.maxbin_tf`).
  - **CONCOCT:** de-interleave, bowtie2 → samtools sort/index, `cut_up_fasta`, `concoct_coverage_table`, then `concoct` all run on TF (`hw.concoct_tf`). `merge_cutup_clustering` and `extract_fasta_bins` run on FM afterwards, as before.
  - **FM's part:** it only fetches the reads from KBase (md5-checked), ships them to TF (md5-checked there), polls and saves. It counts as 1 thread / 2 GB and skips FM's idle-core check.
  - **TF limits:** each TF binner uses 10 threads (MaxBin2 8 GB, CONCOCT 20 GB) within TF's half-and-half shares and its 20-thread / 50 GB cap. It runs only when TF's own load and memory allow, at `nice 10`.
- **Same builds:**
  - **MaxBin2:** TF's tree (MaxBin 2.2.4, bowtie2 2.2.3, HMMER 3.1b1, FragGeneScan 1.30) is md5-identical to FM's.
  - **CONCOCT:** TF's `env_concoct` has the same 65-package set as FM's.
  - **Mapping:** TF's `env_bt2` (bowtie2 2.3.0, samtools 1.3.1) was created 2026-11-25 from FM's explicit package list, and its package set is identical.
- **Pilots:** `maxbin@TF` and `concoct@TF` run one at a time until one is done. A failed one holds only the TF route of that stage (FM and KBase still take it) until `pilot_cleared` is set on it. The first MaxBin2 TF pilots started at 4:07 PM CST: columbia_0074 (per-sample) and NASQAN2015_150.sitecov.

## Same tools, versions and parameters
`~/Documents/metag_env/VERSIONS.md` has the full table. The env is a user-level micromamba, made like `~/Documents/dastool_env`, with the app source at the exact commits KBase ran in `src/`. The same pinned envs are on TF in `~/metag_env`. `hw.py` replays each app's own command sequence, read from the app source and the job logs:
- **MEGAHIT 1.2.9:** the official static build kb_megahit uses, `--min-contig-len 1000`, with the driver's presets and ladders.
- **MetaBAT2 (metabat app bd317234):** the app repo's own static binaries: MetaBAT 2.12.1, bbmap 37.90 (`fast`, `-Xmx30g`), jgi_summarize (`--minContigLength 2500 --minContigDepth 1`) and the app's depth-table merge. samtools 1.8, as shipped.
- **MaxBin2 (kb_maxbin 3574f0b2):** MaxBin 2.2.4 built from source with its own auxiliary versions (bowtie2 2.2.3, HMMER 3.1b1, FragGeneScan 1.30). Settings: `-prob_threshold 0.8 -markerset 40 -min_contig_length 1000`.
- **CONCOCT (kb_concoct 936c9bfe):** CONCOCT 1.1.0, with bowtie2 2.3.0 and samtools 1.3.1 (the Debian-stretch versions of the app image). Settings: cut-up 10 kb, `-l 2500 -k 4 -c 400 -i 500 --total_percentage_pca 90`, merge/extract scripts.
- **CheckM 1.1.2 (kb_Msuite 7495e02e):** `lineage_wf --reduced_tree`, with HMMER 3.1b2, Prodigal 2.6.3, pplacer 1.1.alpha19 and data 2015_01_16.
- **Inputs:** the stored trimmed/merged PairedEndLibrary (interleaved .fastq.gz, md5-checked) and the stored Assembly FASTA, both downloaded from ws 276048. This is exactly what the KBase apps received. Nothing is re-trimmed.

**Recorded differences** (none changes what the tools compute):
- **Threads and memory:** MEGAHIT uses 24–32 threads (63 on KBase) and `-m` as bytes (0.9 on KBase).
- **No SAM on disk:** mapper output is streamed into `samtools view | sort` instead of being written to a SAM file first, so the records are the same.
- **MEGAHIT input:** it gets the interleaved file as `--12` (KBase passes the same pairs as `-1/-2`).
- **MaxBin2:** it runs `bowtie2 -U`, and its own `getsam` reads that output through a FIFO. Then `run_MaxBin.pl -abund_list` does what `-reads_list` does internally.
- **CheckM:** Python 3.8 (3.6 on KBase).
- **CONCOCT:** headers are cut to the contig id instead of kb_concoct's rename-and-revert, so the ids are unchanged.
- **Already random on KBase:** MetaBAT2 2.12.1 is not seeded, and kb_concoct draws a random bowtie2-build seed, so neither is bit-reproducible on KBase either.

**Validation (2026-11-24):** local CheckM was compared with KBase's kb_Msuite on the 62 DAS Tool bins of eastriver_0100. The completeness, contamination, marker lineage and 50/10 call were **identical for all 62 bins**, on FM and on TF (`validation/checkm_eastriver_{FM,TF}/comparison.json`).

## FM assembly size limit and memory watchdog (user-confirmed 2026-11-25; raised the same day)
- **Why:** the first local MEGAHIT, providencecreek_0007 (23.5 GB gz, 46.6 Gbp), grew to 161 GB RSS in its k29 iterate step after 19.7 h, with `-m` at 99 GB, because `-m` caps only the SdBG build. That was over FM's 160 GB limit, with 32 GB of swap in use and systemd-oomd active, so it was stopped by hand at 3:48 AM CST Nov 25.
- **Size rule:** per-sample MEGAHIT runs on FM only when the library is ≤ 16 GB gz (`glue.FM_MAX_LIB_GB`; 10 GB in the first decision, raised to 16 the same morning). providencecreek used ≥ 7 GB of RAM per GB gz, so a 16 GB library may need ~110 GB. Larger libraries and co-assemblies go to poplar (or KBase).
- **One at a time:** at most one MEGAHIT on FM. The gate is asked for the expected peak, 7 GB per GB gz + 10 GB (at most 120 GB; `glue.FM_GB_PER_GZ_GB`), so binners can still use the rest of the per-sample share. If the assembly then grows into that room, the watchdog's MemAvailable floor stops it, and it resumes on poplar.
- **Watchdog (`hw.run_watched`):**
  - An FM MEGAHIT runs in its own session, and its total RSS is checked every 10 s.
  - Above 120 GB RSS (the per-sample share; 100 GB in the first decision), or with FM MemAvailable below 20 GB, it is stopped (SIGTERM, then SIGKILL).
  - The attempt is recorded as `terminated` + `overflow`. That does not count against the unit and does not move the variant ladder.
  - The task is routed to poplar, which resumes it from the kept checkpoints.
- **Order:** while an FM assembly waits for room, the same run starts no new binners, so the assembly gets the next free slot. Running steps are not preempted. The hold applies only when yielding our own FM steps would make that room (MemAvailable plus our steps' memory ≥ request + 30 GB, and the totals fit); memory held by other users' jobs does not hold the binners.

## poplar overflow (user-confirmed 2026-11-25)
- **Layout:**
  - `/scratch/afreiburger/grow_hybrid/`, with `env/megahit/` (the same MEGAHIT 1.2.9 static build as FM; every file's md5 checked identical), `bin/poplar_runner.py` (standard library only) and `runs/<attempt>_<unit>/`.
  - **Since the third decision of 2026-11-25:** each job's run dir (job.json, STATE, logs, staged reads, finished outputs) is `/scratch1/afreiburger/grow_hybrid_stage/runs/<attempt>_<unit>/` on /scratch1 (NFS, ~8 TB free). MEGAHIT's work dir (`-o` and `--tmp-dir`, or the copied checkpoints of a resume) is `/scratch/afreiburger/grow_hybrid/runs/<attempt>_<unit>/` on local /scratch when the gate finds room there. Otherwise it is `/scratch1/afreiburger/grow_hybrid_stage/work/<attempt>_<unit>/` (slower, but /scratch had only ~350 GB free, 97% used). On success the runner moves final.contigs.fa, the MEGAHIT log and options.json into the run dir's `out/` and removes the work dir at once.
  - The two jobs started before that (providencecreek_0007, redbuttecreek_0086) keep the older layout, everything under `/scratch/afreiburger/grow_hybrid/runs/`.
  - metag-hybrid and its envs are not touched.
- **Runner:** `poplar_runner.py <rundir>` reads `job.json`, writes `STATE` and `state.json`, and detaches (double fork), so the ssh call returns at once.
  - **Token:** the KBase token comes only on the ssh stdin, and only when the job downloads inputs. It stays in the runner's memory: never in a file, never in the command's environment (tested 2026-11-25).
  - **Inputs:** they are pulled from KBase on poplar, md5-checked.
  - **Disk:** if the work dir's file system falls below 150 GB free while running (/scratch; 500 GB for a work dir on /scratch1), the job is stopped.
- **FM side (`hw.run_megahit_poplar`):**
  - It polls `STATE` every 5 min. A dropped MFA ControlMaster only delays polling, for up to 72 h.
  - It copies `final.contigs.fa` back and checks its md5, then removes the poplar run dir (and any work-dir remains) at once. It then saves the Assembly from FM (provenance: "MEGAHIT on poplar"). After a resume, the FM checkpoints it came from are removed once the Assembly is saved.
- **Resume (`megahit_resume_plan`):**
  - **Checks first:** a stopped run is continued with `megahit --continue` only after these checks: options.json and checkpoints.txt are present; the checkpoints are a complete prefix that maps onto MEGAHIT 1.2.9's step order (0 library file, 1 build library, 2 first graph, 3 assemble k_min, then local / iterate / build_graph / assemble per k); the recorded paths are the old run's.
  - **Copy:** the files are copied with rsync, and the file lists and sizes are compared afterwards. The read file and reads.lib.bin are also md5-compared.
  - **Path rewrite:** only the recorded paths in options.json, reads.lib and reads.lib.lib_info are rewritten.
  - **Partial edges:** if the killed step was an iterate, the partial `k<next>` edges dir is removed, because build_graph would otherwise read it. iterate recreates it.
  - **Settings kept:** `--continue` keeps the recorded options (threads 32, `-m` 99 GB, meta-sensitive).
  - **Fallback:** if any check fails, the step fails. Setting `no_resume: true` on the task makes it rerun from the KBase library on poplar (`-m` 270 GB).
- **Gate (`gate.can_poplar`; concurrent jobs, user-confirmed 2026-11-25, rules as of the third decision that day):** another poplar MEGAHIT starts whenever ALL of these hold:
  - Our poplar threads stay ≤ 160 (5 jobs of 32; 96 in the second decision).
  - **CPU:** the runnable threads of OTHER users at normal priority (nice ≤ 0; median of 5 one-second `ps` samples) leave at least the job's 32 threads of the 336 cores free. The user's own niced work (EmilyKin Snakemake, ~144 cores at nice 10) doesn't count: our jobs run at the same nice 10 and share CPU with it fairly. This replaced "CPU idle ≥ 30%", which that niced load kept tripping.
  - **Load:** the 1-min load average stays below 1.5 × cores (504).
  - poplar MemAvailable minus the new job's 300 GB budget stays ≥ 400 GB.
  - **Disk:** the work dir goes on local /scratch if /scratch keeps ≥ 200 GB free after the new job's projected work size AND the still-unused part of the running /scratch jobs' projections. Each work dir's current size is measured with `du` in the same check. Otherwise the work dir goes on /scratch1, which must keep ≥ 1000 GB free. The projection is 7 GB per GB gz (reads now sit on /scratch1) for a fresh run, and 250 GB for a resume. providencecreek's dir was 170 GB at k29 for 23.5 GB gz, reads included.
  - Per-sample has priority on poplar. The per-sample driver offers its assemblies largest library first (FM can't take those), so queued KBase jobs are taken over largest first. The sitecov co-assemblies come after the per-sample queue.
  - The FM side counts as 1 thread / 4 GB and is exempt from the FM shares.
- **Pilots:** `assemble` (FM) and `assemble@poplar` are separate stage pilots. On poplar, only one *fresh* run (pulling its reads from KBase) at a time until one has finished; a resume from checkpoints does not count, so the first fresh run starts beside providencecreek.

## No duplicated work (user-confirmed 2026-11-25)
No task ever runs on KBase and on FM / TF / poplar at the same time, in either direction. `glue.py` enforces this for both drivers:
1. **One live attempt per task, across all routes.**
   - **Local start after a KBase job:** before a local attempt starts for a task that has a KBase job, the job is cancelled, and EE2 must then report it `terminated` (or `error`).
   - **If EE2 shows it `running` or `completed` after the cancel:** nothing starts locally. KBase keeps it, and its output is adopted through the normal completion check.
   - **If the cancel isn't confirmed yet:** nothing starts this cycle, and it is tried again next cycle.
2. **No KBase (re)submit while any attempt of the task is live,** local or KBase (a guard in both drivers' submit path). A local result is pushed by its worker before `result.json` is written, so "pending push" is a live local attempt.
3. **Every cycle, right after the poll (`glue.dedup`):**
   - **Several live attempts on one task:** the one furthest along is kept. Running beats queued; on a tie the local one wins, then the older one. The others are cancelled (EE2, confirmed) or killed (FM processes holding the workdir, and the run on TF / poplar). They are recorded as `terminated` + `dedup`, not counted against retries, and logged.
   - **EE2 jobs live in the workspace (276048 / 277632) that the ledger doesn't know:** if the job's params match a task with nothing live, it is adopted. Otherwise it is cancelled. Only the driver's own apps are touched; another app's job is logged and left alone.
   - **Jobs the ledger records as finished but EE2 still shows live:** cancelled if the task is done or has another live attempt, else re-adopted.
   - **Read-only check:** the coordinating session's `dup_audit.py` (exit 1 on any finding) reports the same conditions.
4. **Before every push (`hw.pre_push`, for the Assembly, the BinnedContigs and extract's bins + AssemblySet):**
   - If another version of the output name was saved after the worker started, the local result is not pushed (`superseded`). The first one is kept, so the narrative never gets two competing outputs for one task.
   - **Assemblies and binners:** the driver then validates and adopts the existing object.
   - **extract / CheckM:** the task is held (`dedup_hold`) for a person to look at.
   - Workers that were already running when this went in (weststlouiscreek, the two poplar MEGAHITs, one sitecov CONCOCT) started before this check existed. Their KBase jobs had been cancelled before they started.

## Pushing to KBase (no KBase jobs)
- `kbwrite.py` has its own write allowlist: Blobstore upload in DataFileUtil's request shape (md5-verified), `AbstractHandle.persist_handle`, and `Workspace.save_objects`. It writes only to ws 276048 and 277632, and the read-only `kbrpc.py` is unchanged.
- The objects are built exactly as the KBase code at the used commits builds them, and saved under the drivers' object names:
  - **Assembly:** as AssemblyUtil 59855bcb, with contig md5/gc and the md5 of the sorted contig md5s.
  - **BinnedContigs:** as MetagenomeUtils 49489d11 `file_to_binned_contigs`, including each app's summary file.
  - **Bin Assemblies plus AssemblySet:** as extract.
- **Provenance:** each object's provenance names the KBase app it stands in for, plus "run LOCALLY on host … with …" and the tool versions.
- **CheckM in the narrative:**
  - Each bin Assembly carries `CheckM completeness`, `CheckM contamination` and `CheckM marker lineage` as object metadata.
  - A `KBaseReport.Report` named `<unit>.CheckM_local.report` has the same `CheckM_summary_table.tsv.zip` and `full_output.zip` file links as kb_Msuite's report, so the drivers' `collect_checkm_select.py` reads it unchanged.
  - For units whose extract had already run in KBase, the CheckM values are in that report only, since existing objects' metadata cannot be changed.

## Files
- `hw.py`: the worker, one step per process. The driver writes `<workdir>/spec.json` and the worker writes `result.json`. It never touches a ledger.
- `glue.py`: the driver side (launch, takeover, poll, stage pilots), shared by both drivers.
- `gate.py` and `slots.json`: the shared resource gate.
  - **FM:** at most 48 threads and 160 GB for our steps (DAS Tool workers included), MemAvailable ≥ request + 30 GB, enough idle cores by load average, ≥ 100 GB free disk plus the step's need, and never while a GTDB-Tk scope is active.
  - **TF:** at most 20 threads and 50 GB, nice 10, load and free disk checked over ssh.
  - **Shares (user-confirmed 2026-11-24; replaces the strict per-sample priority of 2026-11-23):**
    - FM: sitecov has 12 threads / 40 GB reserved, per-sample the other 36 threads / 120 GB. A step that keeps its run inside its own share may always start (if the totals allow).
    - A run may go beyond its share (borrow) only while the other run reported no ready-but-waiting steps in its last cycle (< 2 h). Borrowing stops when the owner has work again; running steps are never preempted, so a borrowed share returns as those steps finish.
    - TF: half each (10 threads / 25 GB), on the same rule, using the counts of waiting TF-eligible steps.
    - Binner, extract and CheckM steps are capped at the run's share in threads (12 for sitecov); threads only change speed. A co-assembly (150 GB) exceeds the sitecov share and runs only on borrowed capacity or, later, poplar.
    - Per-sample keeps priority for KBase submissions and poplar overflow.
- `cache/reads/`: downloaded libraries. They are removed once neither driver needs them (`needs_*.json`).
- `validation/`: comparison runs against KBase outputs.
- **Work dirs:** `<driver>/work/hybrid/<stage>/<unit>/attempt*/`, holding `spec.json`, `worker.log`, the tool logs and `result.json`. Large intermediates are deleted on success.

**Stage pilots:** each local stage runs one step at a time until one has finished `done`, then fans out.

**Times:** the logs are in UTC; reports to the user are in US Central.
