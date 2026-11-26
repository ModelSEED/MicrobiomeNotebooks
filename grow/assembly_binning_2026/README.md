# assembly_binning_2026: assembling and binning the 133 GROWdb river metagenomes

This directory rebuilds MAGs for the 133 GROWdb samples, the source of the 3,825 MAGs in ws 145226, from the raw SRA reads (171 runs, 3.72 Tbp, all 2x151). It uses the same pipeline as `salterns/reassembly_2026` (driver v2.2): the same apps, tags, settings and safeguards. Where GROW differs, the choice is recorded below. The brief this was built from is `BRIEF.md`.

- **Reads (input, read only):** ws **275929** `freiburgermsu:narrative_1794152940209`, one interleaved PairedEndLibrary `<sample>.<SRR>.reads` per run.
  - They are imported by `../raw_metagenomes_sra/kbase_import.py`, whose ledger this driver only reads.
- **Outputs:** ws **276048** `freiburgermsu:narrative_1794250260895`, "GROWdb metagenomes: assembly and binning (2026)".
  - It holds the trimmed and merged reads, assemblies, BinnedContigs, extracted bins and CheckM reports.
  - KBase jobs run with ws 276048 as their workspace and read the libraries from 275929 by absolute ref.
  - The user chose this on 2026-11-09 to keep the reads narrative from growing too large.
- **Manifest:** `manifest.tsv` has one row per run (171): sample, run, source, BioProject, BioSample, spots, bases, runs per sample, sample bases and MAGs per sample. `make_manifest.py` builds it from `../raw_metagenomes_sra/runs.tsv`.
  - 96 samples have 1 run, 36 have 2 and 1 has 3, so 37 samples are multi-run. BRIEF.md says 46; `runs.tsv` gives 37.

## Decisions (2026-11-09)

Items 2, 3, 5 and 7 were **confirmed by the user** on 2026-11-09 (relayed by the import session grow-a2): no co-assembly, per-sample coverage, GROW behind salterns at reduced caps, and no GTDB-Tk alongside salterns'.

1. **The unit is the sample (133).**
   - Trimmomatic runs per run.
   - Multi-run samples have their trimmed runs merged with `kb_ReadsUtilities` merge. The merge is accepted only if `read_count` and `total_bases` equal the sums of the inputs; otherwise the sample is not assembled.
   - A sample starts as soon as all of its runs are `done` in the import ledger.
   - A run the import watcher gave up on (`failed`) is not mirrored as failed. Its sample, and therefore GTDB-Tk, waits for a re-import. The watcher logs such runs.
2. **No co-assembly.** The samples are independent rivers, not replicate cores.
3. **Binning coverage is per-sample:** each sample's own library, trimmed and merged, as GROWdb did.
   - *A possible later improvement:* multi-sample coverage, for example within a river or basin. Across 133 libraries it would be far heavier (bowtie2 of up to 133 × 3.7 Tbp).
4. **Names follow salterns, with the sample as the unit:**
   - `<sample>.<SRR>.trim_paired` and `<sample>.runsN.trim_paired.merged.reads`;
   - `<sample>.MEGAHIT.Assembly` (the superseded deadrun pilot was `<sample>.metaSPAdes.Assembly`);
   - `<sample>.{MetaBAT2,MaxBin2,CONCOCT,DASTool}.BinnedContigs` and `<sample>.DASTool_bins.AssemblySet`;
   - bins are named `<bid>__<sample>`.
5. **KBase capacity is shared with salterns, and salterns has priority.** Each cycle the driver reads `../../salterns/reassembly_2026/ledger.json` (read only, never locked).
   - **Reduced caps:** while any salterns import, trim, merge, assemble, coassemble or binner task is not finished, this driver is capped at **extreme 3, bigmem 6, njs 6**.
   - **After salterns clears those stages:** the caps become the salterns totals (12/20/15) minus the live salterns jobs, and for njs also minus the import watcher's live imports. They never drop below 3/6/6.
   - Cap changes are logged.
6. **Submission order is downstream first,** so that finishing a sample (MetaBAT2, the MEGAHIT fallback) is never starved by the 168 pending trims in the same bigmem pool.
7. **Local resources:**
   - At most **2** DAS Tool workers.
   - GTDB-Tk runs in its own scope **`gtdbtk-grow-assembly`** and never starts while the salterns GTDB-Tk task is not finished or its scope `gtdbtk-salterns-reassembly` is active.
   - Local steps (DAS Tool, CheckM collection, GTDB-Tk) start only with **≥ 100 GB** free.
8. **Sample pilots:**
   - The two pilots are `deadrun_2019_sw_WHONDRS-S19S_0006` (1 run, 13.0 Gbp) and `eastriver_2019_sw_WHONDRS-S19S_0100` (2 runs, 67.7 Gbp, the largest sample).
   - **Assembly and the three binners** fan out to the other samples only after both pilots are assembled and binned (all three binners finished, at least one with output). This follows the brief.
   - **Trim and merge** fan out once both pilots have passed that stage. *This differs slightly from the brief:* the steps are cheap, they change nothing downstream, and the other samples' libraries are then ready when the assembly gate opens. Bigmem stays capped at 6, and downstream work goes first (item 6).
   - A pilot task that ends `failed` or `empty` holds its gate, and the watcher logs it. Inspect it before loosening anything.
   - On top of this, DAS Tool, extract and CheckM each keep salterns' single-unit pilot.
9. **MEGAHIT for every sample (user decision, 2026-11-10; relayed by grow-a2 and recorded in BRIEF.md).**
   - **Assembler:** every sample is assembled with MEGAHIT/run_megahit (release, 1.2.9), meta-sensitive, on bigmem. metaSPAdes is not used at all.
   - **Floor:** 10 Mbp, as before.
   - **Retries:** these follow the salterns co-assembly ladder. After an intrinsic MEGAHIT failure (including an out-of-memory error) or an assembly under 10 Mbp, the next try uses min-count 2 with k 27..127, then the MEGAHIT defaults. Other failures get 3 counted tries.
   - **Reasons:**
     - The river samples are high-diversity with low per-genome coverage.
     - Memory: 13 samples exceed 50 Gbp.
     - MEGAHIT runs on bigmem, which frees the extreme pool for MaxBin2 and CONCOCT.
     - One assembler is used throughout.
     - GROWdb itself used mostly MEGAHIT.
   - **Effect on the caps:** MEGAHIT now shares the 6 bigmem slots with trim and MetaBAT2. Downstream-first ordering (item 6) gives it the next free slot, but bigmem, not extreme, is now this driver's bottleneck while salterns is busy.
   - **Pilot gate:** the gate for assembly and binning applies again to the MEGAHIT pilots. The other 131 samples are assembled once both pilots' MEGAHIT assemblies are binned.
   - **The metaSPAdes pilot outputs are superseded, not deleted:**
     - Superseded objects: the deadrun assembly 276048/28/1, MetaBAT2 56/1, MaxBin2 52/1, CONCOCT 46/1, DAS Tool 58/1, extract 91/1 and the CheckM report 101/1. The CheckM job had already finished (11 MAGs ≥ 50/< 10).
     - In the ledger these tasks are under `superseded`, the refs under `superseded_refs`, and the event under `history`.
     - The driver never adopts a superseded object version. The new binner, DAS Tool and extract outputs reuse the object names, so they become version 2 of the same objects.
     - The local files are in `work/superseded_metaspades/`.
     - The ledger before the switch is `ledger.pre_megahit_switch_20261110.json`.
     - eastriver's metaSPAdes job 6af330f5eacecf21bd99a735 (2nd attempt; the 1st had errored) was cancelled at 2026-11-10 20:49Z.
     - **Kept for a MEGAHIT-vs-metaSPAdes comparison (user, 2026-11-10).** Do not delete or overwrite these, and never adopt them:
       - The import session copied the v1 objects to names the driver never uses: `<sample>.metaSPAdes_pilot.{CONCOCT,MaxBin2,MetaBAT2,DASTool}.BinnedContigs` and `….metaSPAdes_pilot.DASTool_bins.AssemblySet` (276048/138–142), plus the 16 bins as `metaSPAdes_pilot.<bin name>`.
       - The assembly 276048/28 (`….metaSPAdes.Assembly`) and the CheckM report 101 keep their own unique names.
       - The ref map is in `work/superseded_metaspades/kbase_preserved_copies.json`.
       - Keep `work/superseded_metaspades/` (the CheckM selection: 11 MAGs ≥ 50/< 10).
       - An optional local GTDB-Tk run of those 11 MAGs, as a separate set, has not been done.

10. **Transient KBase errors and a relaxed pilot gate (user decisions, 2026-11-12; relayed by grow-a2 and recorded in BRIEF.md).**
    - **Transient EE2 errors:**
      - An EE2 job that ends with "Unexpected Job Error" (including "killed by ./monitor_jobrunner_logs.py") is a KBase-side outage, not a data problem. Such attempts are marked `transient` and are not counted toward MAX_TRIES or the hard cap of 6.
      - A run of them gives up 72 h after the first one ended. The budget runs from the end, because an EE2 job can sit in the queue for hours.
      - Evidence: two mass failures, 2026-11-11T21:54Z and 2026-11-12T10:55Z, hit many unrelated trim, merge and MEGAHIT jobs at once.
      - Past attempts were reclassified as transient (32 attempts).
      - Reset to retry: the 3 tasks that had given up, trim:SRR35424631 and merge of oksrukuyikcreek_0050 and orco_0018.
      - "No output generated" still counts as a real failure.
    - **Pilot gate:** deadrun proved the whole MEGAHIT chain (13 MAGs). Assembly and binning therefore fan out to every sample under 50 Gbp (summed over its runs). Only the 13 samples of 50 Gbp or more still wait for the eastriver pilot to be assembled and binned.
    - **Trim priority:** with fewer than 10 trims left to run, trim is submitted before the other stages. Otherwise ~120 queued MEGAHIT jobs would hold bigmem and starve the samples those trims block.
    - **Caps:** salterns finished on 2026-11-11 (its watcher was stopped), so the caps are at their post-salterns levels: extreme 12, bigmem 20, njs 15. These are the shared totals minus live salterns and import jobs, of which there are none.
    - **Ledger backup:** the ledger before this change is `ledger.pre_transient_gate_20261112.json`.

## Pipeline (`pipeline.py`)

| Stage | App (tag), client group | Notes |
|---|---|---|
| import | read-only mirror of `../raw_metagenomes_sra/kbase_import_ledger.json` | Never submitted here. |
| trim | kb_trimmomatic/run_trimmomatic (release), bigmem | Per run. Trimmomatic 0.39 PE: TruSeq3-PE-2 2:30:10, LEADING:3, TRAILING:3, SLIDINGWINDOW:4:15, MINLEN:50. If trimming fails, the raw import is used. |
| merge | kb_ReadsUtilities/KButil_Merge_MultipleReadsLibs_to_OneLibrary (release), njs | Multi-run samples only. The read and base sums are validated. |
| assemble | MEGAHIT/run_megahit (release, 1.2.9), bigmem | meta-sensitive for every sample (since 2026-11-10; item 9). Retry ladder: min-count 2 with k 27..127, then defaults. Output `<sample>.MEGAHIT.Assembly`. The floor is 10 Mbp. |
| metabat, maxbin, concoct | metabat/run_metabat (beta, bigmem), kb_maxbin/run_maxbin2 (release, extreme), kb_concoct/run_kb_concoct (release, extreme) | `reads_list` is the sample's own library. |
| dastool | local DAS Tool 1.1.7 (`local/dastool_worker.py`, `~/Documents/dastool_env`) | diamond, score 0.5, duplicate penalty 0.6, megabin penalty 0.5, 8 threads. Saves `<sample>.DASTool.BinnedContigs` to ws 276048. |
| extract | MetagenomeUtils/extract_bins_as_assemblies (beta), njs | |
| checkm | kb_Msuite/run_checkM_lineage_wf (beta, CheckM 1.1.2), njs | Reduced tree. `local/collect_checkm_select.py` keeps bins ≥ 50 % complete and < 10 % contaminated. |
| gtdbtk | local GTDB-Tk 2.8.0 / R232 (`local/run_gtdbtk_reassembly.sh`, `~/gtdbtk_kbase/env.sh`) | One batch, once the CheckM of **every one of the 133 samples** is finished and collected. Runs in a systemd user scope (MemoryMax 170G, no swap) under `local/mem_guard.sh`, after MemAvailable ≥ 150 GB. Logs go to `~/gtdbtk_kbase/grow_assembly_2026/`; summaries are copied to `gtdbtk/`. |

The safeguards are the same as salterns:
- a write-ahead ledger (`ledger.json`) with the exact parameters, tag and service version of every attempt;
- lost submits reconciled with EE2 before any resubmission;
- existing outputs validated before they are adopted;
- 3 counted tries per task (cancels and transport errors excluded) and a hard cap of 6;
- 24-hour transient budgets for local workers;
- a GTDB-Tk relaunch after an OOM kill.

`local/` holds copies of the salterns helpers, changed only to write to ws 276048. `dastool_worker.py` now passes `--ws-id 276048` to `dastool_to_binnedcontigs.py`.

## Hybrid local/KBase mode (user-approved 2026-11-23)
The KBase queue waits reached 35–48 h, so the stalled work no longer goes to KBase. It runs on our own machines with the same tools, versions and parameters (`../hybrid/README.md`), and every result is saved back into **ws 276048** under the usual object names, with provenance saying it ran locally.
- **Stages run locally:** MEGAHIT, the three binners, extract and CheckM (`HYBRID` in `pipeline.py`). DAS Tool and GTDB-Tk were already local.
- **Takeover:** a queued KBase job is cancelled when a local slot is free for its unit (`takeover: true`, not counted). Running jobs finish in KBase.
- **No duplicates (user-confirmed 2026-11-25):** one live attempt per task across KBase / FM / TF / poplar, enforced every cycle (`../hybrid/README.md`).
- **Where:**
  - **FM** coordinates and makes every KBase push.
  - **TF** runs CheckM when FM is full.
  - **MEGAHIT routing (user-confirmed 2026-11-25):** each remaining assembly goes to the first route that can take it: FM (libraries ≤ 16 GB gz, one at a time, memory watchdog at 120 GB RSS / 20 GB MemAvailable), poplar (up to 5 concurrent 32-thread jobs, ≤ 160 threads of ours, since the third decision of 2026-11-25; largest libraries first), else a KBase MEGAHIT job as the fallback queue (taken over by FM or poplar while still queued). Never two routes at once for a sample.
- **Limits:** one shared gate (`../hybrid/gate.py`) covers both drivers: ≤ 48 threads and ≤ 160 GB on FM, ≥ 100 GB free, no step while GTDB-Tk runs. **Shares (user-confirmed 2026-11-24):** sitecov has a reserved 12 threads / 40 GB of the FM budget, and the per-sample run the other 36 threads / 120 GB. Either run may borrow the other's idle share only while the other reports no ready-but-waiting steps, and stops borrowing as soon as it does; running steps are never preempted. TF is split half and half (10 threads / 25 GB each) on the same rule. Per-sample keeps priority for KBase submissions and poplar overflow.
- **Pilots:** each local stage first runs one step until it is `done`.
- **CheckM:** local CheckM matched KBase's kb_Msuite exactly on 62 bins, on FM and on TF. The values are on each bin Assembly as metadata and in `<unit>.CheckM_local.report`.

## Running

```bash
cd grow/assembly_binning_2026
~/Documents/py_venv/bin/python pipeline.py status
~/Documents/py_venv/bin/python pipeline.py plan --through gtdbtk      # what would start now; submits nothing
U=grow-assembly-watch-$(date -u +%Y%m%dt%H%M%Sz)                      # a new scope name for every start
setsid nohup systemd-run --user --scope --unit=$U -p OOMPolicy=continue -p ManagedOOMPreference=avoid \
    ~/Documents/py_venv/bin/python -u pipeline.py watch --through gtdbtk --interval 1800 >> watch.log 2>&1 < /dev/null &
```

- **Finding the watcher:** the salterns watcher has the same command line. To find ours, run `pgrep -af "[p]ipeline.py watch"` and keep the PID whose `readlink /proc/<pid>/cwd` is this directory. Never use `pkill -f`/`pgrep -f` with a pattern that also appears in your own command line.
- **Starting it:** start it from a shell that has the user bus (`XDG_RUNTIME_DIR`, `DBUS_SESSION_BUS_ADDRESS`), with a new scope name every start.
- **Restarting it:** stop the python PID while it sleeps (`watch.log` written 1–29 min ago). Never stop the scope itself: that would also kill running DAS Tool workers and the GTDB-Tk memory guard. Live workers are adopted from the ledger.
- **One writer:** only one `advance` or `watch` may run at a time (`ledger.json.lock`).
- **After a reboot** (this machine reboots unexpectedly): check the watcher, and restart it as above if it is dead.
- **Never write** the import ledger, and never stop the import watcher or the salterns watcher.
- **Keep `work/` until the end:**
  - `work/dastool/*/attempt*/<sample>.save_objects.json` is read until extract has finished;
  - `work/checkm/*/fasta` and the batchfiles are read until GTDB-Tk has finished.

## History
- 2026-11-09 18:55Z: the first `advance` submitted the 3 pilot trims (deadrun SRR35181236; eastriver SRR35180665 and SRR35180666) to ws 276048. The watcher was started in scope `grow-assembly-watch-20261109t185523z` (PID 219790). At that point 59 of 171 runs were imported, and 35 samples were ready.
- 2026-11-10 20:49Z: switched to MEGAHIT for every sample (item 9). The watcher was stopped by PID while it slept, eastriver's metaSPAdes was cancelled, the deadrun metaSPAdes chain was superseded, and the watcher was restarted in scope `grow-assembly-watch-20261110t205038z` (PID 3817250). At that point 130 of 171 runs were imported, 29 trims and 9 merges were done, and both pilot libraries were ready for MEGAHIT.
- 2026-11-12 19:32Z: transient EE2 errors and the relaxed pilot gate (item 10). The watcher was restarted in scope `grow-assembly-watch-20261112t193255z` (PID 1995190). Its first cycle submitted 18 MEGAHIT assemblies (all samples < 50 Gbp), the reset trim and 2 merges. With eastriver's MEGAHIT (attempt 3, 6af59b16…) that makes 20 live bigmem jobs.
- 2026-11-24, hybrid mode (user decision 2026-11-23):
  - Both watchers were stopped by PID at 9:19 PM CST on Nov 23 while the local mode was built, then restarted at 7:59 AM CST on Nov 24. This run's watcher was restarted three more times by 8:28 AM CST (`grow-assembly-watch-20261124t142833z`) to load fixes: a short TMPDIR for CheckM, DAS Tool counted in the gate, smaller budgets and TF routing.
  - Its first hybrid cycle took over the queued KBase MEGAHIT of providencecreek_0007, which now runs locally.
  - CheckM pilot: NASQAN2015_227, report 276048/2795/1, 25 bins (17 at ≥ 50/< 10). collect_checkm_select.py read it unchanged.
  - CheckM fan-out: 3 queued KBase CheckM jobs were taken over, 1 on FM and 2 on TF.
- 2026-11-24 9:52 AM CST, reserved sitecov share (user decision 2026-11-24, relayed by the coordinating session): the gate's strict per-sample priority was replaced by shares (FM 36/120 per-sample, 12/40 sitecov; TF half each). Both watchers were restarted by PID to load it, this one in scope `grow-assembly-watch-20261124t155156z`. Backups: `../hybrid/gate.pre_share_20261124.py` and `../hybrid/glue.pre_share_20261124.py`.
- 2026-11-25 6:26–6:31 AM CST, poplar overflow + FM size limit (user decision 2026-11-25, relayed by the coordinating session; details in `../hybrid/README.md`):
  - **FM:** per-sample MEGAHIT runs on FM only for libraries ≤ 10 GB gz, under a memory watchdog (100 GB RSS / 20 GB MemAvailable; 120 GB since the second decision that morning). The first FM sample is weststlouiscreek_0027 (8.2 GB gz). It waits for the two running per-sample binners, and no new binners start before it.
  - **poplar:** 16 larger libraries go there, 1 at a time (32 threads).
  - **providencecreek_0007:** its attempt 2 was re-marked `terminated` + `overflow` (backup `ledger.pre_poplar_20261125.json`). Attempt 3 resumes it on poplar from the FM checkpoints (`megahit --continue`; 5 checkpoints done, next the k29 iterate).
  - **Watchers:** restarted by PID, the last time in scope `grow-assembly-watch-20261125t123125z`. Backups: `../hybrid/{hw,glue,gate}.pre_poplar_20261125.py`.
- 2026-11-25 9:12 AM CST, concurrent poplar + FM limit 16 GB + KBase fallback (user decision 2026-11-25, second one that morning, relayed by the coordinating session; details in `../hybrid/README.md`):
  - **poplar:** several MEGAHIT jobs at once (≤ 96 threads of ours; CPU idle ≥ 30%; MemAvailable − 300 ≥ 400 GB; `/scratch` ≥ 200 GB free after the projected peaks). The run dir is deleted as soon as the contigs are back and md5-verified.
  - **FM:** threshold raised from 10 to 16 GB gz, one assembly at a time, watchdog at 120 GB RSS (the per-sample share) with the 20 GB MemAvailable floor.
  - **KBase:** the fallback queue for assemblies no local route can take (`pipeline.py` falls through to the submit when `glue.try_run` returns False).
  - **First cycle (9:12 AM CST):** weststlouiscreek_0027 (8.2 GB gz, pinned to FM as the end-to-end pilot: `route: FM` in the ledger, backup `ledger.pre_concurrent_20261125.json`) started on FM, 24 threads, 67 GB requested. redbuttecreek_0086 (11.8 GB gz) started as the first fresh poplar run (the poplar pilot) beside the providencecreek resume. The other 14 got KBase MEGAHIT jobs (queued). They are taken over as FM or poplar get room.
  - **Watcher:** restarted by PID, scope `grow-assembly-watch-20261125t151247z`, and again at 9:17 AM CST (`grow-assembly-watch-20261125t151713z`; once more as `grow-assembly-watch-20261125t151745z` to count KBase-queued samples as waiting for poplar, keeping per-sample's poplar priority) so that binners are held for an FM assembly only when yielding our own steps would make room. Backups: `../hybrid/{hw,glue,gate}.pre_concurrent_20261125.py`, `pipeline.pre_kbase_fallback_20261125.py`.
- 2026-11-25 3:36 PM CST, more of the assemblies to poplar (user decision 2026-11-25, third that day, relayed by the coordinating session; details in `../hybrid/README.md`):
  - **Concurrency:** poplar now runs up to 5 MEGAHIT jobs (≤ 160 of our threads), at nice 10 like the user's own EmilyKin work there.
  - **CPU rule:** "CPU idle ≥ 30%" became "other users' normal-priority runnable threads leave ≥ 32 cores free, and 1-min load < 1.5 × cores".
  - **Disk:** reads, control files and outputs are staged on /scratch1 (`/scratch1/afreiburger/grow_hybrid_stage`). MEGAHIT's work dir stays on local /scratch only while the 200 GB floor allows, otherwise it goes on /scratch1.
  - **Order:** assemblies are offered largest library first, so queued KBase jobs move to poplar (cancelled while queued) largest first. The sitecov co-assemblies come after that.
  - **Takeover candidates:** by then KBase had started all but 2 of the 14 fallback jobs, which can no longer move. sangamon, southforkpolouse and redbuttecreek_0029 finished there. Still queued: yakimariver (16.9 GB gz) and whiteclaycreek2 (14.5). They move to poplar once the fresh-run pilot, redbuttecreek_0086, has finished, as long as KBase hasn't started them by then.
  - **Watcher:** restarted by PID, scope `grow-assembly-watch-20261125t213601z`. Backups: `../hybrid/*.pre_scratch1_20261125.py`, `pipeline.pre_scratch1_20261125.py`.
- 2026-11-25 3:44 PM CST, no duplicated work (user decision 2026-11-25, relayed by the coordinating session; details in `../hybrid/README.md`, "No duplicated work"): one live attempt per task across KBase / FM / TF / poplar. A local start needs EE2's confirmed cancel. KBase (re)submits are blocked while anything is live. A dedup sweep runs every cycle (duplicates, EE2 orphans, jobs the ledger records as finished), and a pre-push check keeps the first of two competing outputs. A tested run on copies of both ledgers found nothing to fix, as did the coordinating session's audit at ~4 PM CST. The watcher was restarted by PID, scope `grow-assembly-watch-20261125t214430z`. Backups: `../hybrid/{glue,hw}.pre_dedup_20261125.py`, `pipeline.pre_dedup_20261125.py`.
- 2026-11-25 4:04–4:07 PM CST, binners on every route (user decision 2026-11-25, relayed by the coordinating session; details in `../hybrid/README.md`, "Binners on TF"): no binner had run for hours. FM's share was taken by the FM MEGAHIT plus DAS Tool / extract, FM was overloaded by other users, and the binners had no other route. Each binner step now goes to FM within the share, then TF (MaxBin2 / CONCOCT, mapping included), else a KBase job (fallback queue, taken over while queued; never two routes). The first cycle submitted KBase binner jobs up to the pool caps (extreme 12 / bigmem 20, shared with its 12 KBase MEGAHITs). The first MaxBin2 TF pilot started at 4:07 PM CST. The watcher was restarted by PID, scope `grow-assembly-watch-20261125t220737z`. Backups: `../hybrid/{glue,hw}.pre_tfbin_20261125.py`.
- 2026-11-25 3:48 AM CST: the providencecreek_0007 local MEGAHIT (pilot, 32 threads, `-m` 99 GB) grew to 161 GB in its k-iteration step after 19.7 h, with 32 GB swap in use and FM MemAvailable at 3 GB. That is over the 160 GB limit, and systemd-oomd would have hit other users' jobs first (our scope is OOM-avoid), so it was stopped by hand (SIGTERM, rc 241, recorded as an intrinsic failure). The failed pilot holds the local assemble stage (no automatic retry). Its MEGAHIT checkpoints are kept in the attempt's `tmp/megahit_out` for a possible `--continue` on a machine with more memory (poplar).
  - Backups before the change: `pipeline.pre_hybrid_20261124.py` and `ledger.pre_hybrid_20261124.json`.
