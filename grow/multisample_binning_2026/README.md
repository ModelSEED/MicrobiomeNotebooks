# multisample_binning_2026: GROW site multi-sample binning ("sitecov") and co-assembly pilots

This is a parallel effort to the per-sample GROW run (`../assembly_binning_2026`), approved by the user on 2026-11-18. The brief is in `BRIEF.md`.
- It re-bins the 67 GROW samples that share a site, using **site multi-sample coverage**.
- It also co-assembles **8 pilot site groups**.
- It reuses the per-sample run's trimmed/merged reads libraries and its MEGAHIT assemblies (ws 276048, by ref), taken from that run's `ledger.json` and `manifest.tsv`. Those files are only read: never written or locked here.

Everything this effort makes goes to **ws 277632** `freiburgermsu:narrative_1795046239204`, "GROWdb metagenomes: site multi-sample binning and co-assembly (2026)". It is kept separate from the per-sample MAGs and from GROWdb's 3,825/2,093 MAGs; never merge or dereplicate them together.

## Scope
- **`manifest.tsv`:** one row per grouped sample, 67 samples in 25 site groups. `make_manifest.py` builds it from `site_groups.tsv`, adding each sample's SRA runs.
  - 20 NASQAN/USGS station time series, one station per group (not river-wide).
  - 5 WHONDRS same-site pairs.
- **Out of scope (user decision):** the 66 unpaired samples get nothing new.
- **Units:**
  - **`<sample>.sitecov`** (67): the sample's per-sample MEGAHIT assembly, binned with `reads_list` = the libraries of every sample in its site group (2–5 libraries).
  - **`<group>.coassembly`** (8 `coassembly_pilot` groups): the 5 WHONDRS pairs and USGS 03609750 Tennessee, 03378500 Wabash and 04264331 St. Lawrence.
    - The member libraries are merged (read and base sums validated), then co-assembled with MEGAHIT meta-large.
    - The co-assembly is binned with the members' **individual** libraries, as in salterns, not the merged library.
- **Names:**
  - sitecov: `<sample>.sitecov.{MetaBAT2,MaxBin2,CONCOCT,DASTool}.BinnedContigs` and `<sample>.sitecov.DASTool_bins.AssemblySet`; bins are `<bid>__<sample>.sitecov`.
  - co-assembly: `<group>.n<N>.trim_paired.merged.reads` and `<group>.MEGAHIT.coassembly.Assembly`, then `<group>.coassembly.*` as above; bins are `<bid>__<group>.coassembly`.

## Pipeline (`pipeline.py`)

It uses the same apps, versions and parameters as salterns and the per-sample run, copied from the per-sample driver (`../assembly_binning_2026/pipeline.py`, unchanged).

| Stage | App (tag), pool | Notes |
|---|---|---|
| merge | kb_ReadsUtilities merge (release), njs | Co-assembly groups only. Sums are validated. |
| coassemble | MEGAHIT/run_megahit (release, 1.2.9), bigmem | meta-large. The salterns ladder applies after an intrinsic failure or an assembly < 50 Mbp: min-count 2 with k 27..127, then defaults. |
| metabat, maxbin, concoct | metabat (beta, bigmem), kb_maxbin (release, extreme), kb_concoct (release, extreme) | `reads_list` is the site group's member libraries. |
| dastool | local DAS Tool 1.1.7 | Saves `<unit>.DASTool.BinnedContigs` to ws 277632. |
| extract | MetagenomeUtils extract (beta), njs | |
| checkm | kb_Msuite CheckM (beta, reduced tree), njs | Bins ≥ 50 % complete and < 10 % contaminated are collected. |
| gtdbtk | local GTDB-Tk 2.8.0 / R232 | Two batches, `sitecov` (all 67 units) and `coassembly` (all 8), one at a time (whichever is ready first) in scope `gtdbtk-grow-sitecov`. Output goes to `~/gtdbtk_kbase/grow_sitecov_2026/{sitecov,coassembly}/`, with summaries in `gtdbtk/{sitecov,coassembly}/`. |

**Scheduling (the brief: this effort queues behind the per-sample run):**
- **Caps:** while the per-sample ledger has any `assemble` or binner task unfinished, the caps are **extreme 3, bigmem 6, njs 6**. After that they are 12/20/15 minus the per-sample run's live jobs, never below 3/6/6.
  - The per-sample driver does not count these jobs, so the combined load can exceed 12/20/15. The brief accepts that.
- **Local DAS Tool:** at most 1 worker while the per-sample run has a DAS Tool worker live, 2 otherwise. Local steps need ≥ 100 GB free.
- **Pilot:** the Mobile River pair (USGS 02470500, the smallest NASQAN pair: 2 samples, 24.6 Gbp) goes first. Every other binner (sitecov and co-assembly) waits until both Mobile units have all three binners finished, at least one with output.
  - The co-assembly merges and MEGAHIT runs need only reads, so they go at once.
  - DAS Tool, extract and CheckM keep the single-unit stage pilots.
- **Order:** downstream stages are submitted first.
- **Transient errors:** an EE2 "Unexpected Job Error" (including `monitor_jobrunner_logs.py` kills) is transient. It is not counted and gets a 72-h budget, as decided by the user on 2026-11-12.
- **GTDB-Tk:** never starts while the per-sample GTDB-Tk (`gtdbtk:ALL`) is unfinished, or while scope `gtdbtk-grow-assembly` or `gtdbtk-salterns-reassembly` is active. The per-sample run has priority.

**Deviation from the brief, with the reason.** The brief says a group starts once every member's assembly and reads are done. Instead, each `<sample>.sitecov` unit starts once **its own** assembly is done and **all its group members' libraries** are settled.
- NASQAN2011_239 (USGS 07374000 Mississippi) failed per-sample assembly (< 10 Mbp on all three MEGAHIT tries). Under the literal rule its group would never start.
- Under this rule, NASQAN2011_239 gets no sitecov unit (`empty`), but its reads still feed the coverage of its two group mates.
- A unit with fewer than 2 member libraries would only repeat the per-sample binning, so it is skipped.

## Comparison (`compare.py`, deliverable at the end)
`compare.py` reads the local CheckM collections of both runs and writes `comparison/`. It can be re-run at any time; missing sides show as NA.
- **`per_sample_vs_sitecov.tsv`:** MAGs ≥ 50/< 10 per sample, per-sample vs sitecov, with median completeness, contamination and N50.
- **`coassembly_vs_sitecov_union.tsv`:** for the 8 groups, co-assembly MAGs vs the union of the members' sitecov MAGs.
- **`distributions.tsv`:** completeness, contamination and N50 quartiles per MAG set.
- At the start (2026-11-18 6:04 PM CST), the per-sample side had 30 of 67 samples collected, 279 MAGs.

## Hybrid local/KBase mode (user-approved 2026-11-23)
The KBase queue waits reached 35–48 h, so the stalled work no longer goes to KBase. It runs on our own machines with the same tools, versions and parameters (`../hybrid/README.md`), and every result is saved back into **ws 277632** under the usual object names, with provenance saying it ran locally.
- **Stages run locally:** the co-assemblies, the three binners, extract and CheckM (`HYBRID` in `pipeline.py`). DAS Tool and GTDB-Tk were already local.
- **Takeover:** a queued KBase job is cancelled when a local slot is free for its unit (`takeover: true`, not counted). Running jobs finish in KBase.
- **No duplicates (user-confirmed 2026-11-25):** one live attempt per task across KBase / FM / TF / poplar, enforced every cycle (`../hybrid/README.md`).
- **Where:**
  - **FM** coordinates and makes every KBase push.
  - **TF** runs CheckM when FM is full.
  - **poplar** is overflow for co-assemblies only.
- **Limits:** one shared gate (`../hybrid/gate.py`) covers both drivers: ≤ 48 threads and ≤ 160 GB on FM, ≥ 100 GB free, no step while GTDB-Tk runs. **Shares (user-confirmed 2026-11-24):** sitecov has a reserved 12 threads / 40 GB of the FM budget, and the per-sample run the other 36 threads / 120 GB. Either run may borrow the other's idle share only while the other reports no ready-but-waiting steps, and stops borrowing as soon as it does; running steps are never preempted. TF is split half and half (10 threads / 25 GB each) on the same rule. Per-sample keeps priority for KBase submissions and poplar overflow. Binner steps of this run use 12 threads to fit the share. The co-assemblies run on poplar (concurrent jobs), or as KBase MEGAHIT jobs when poplar cannot take them; a queued KBase job is taken over by poplar (user-confirmed 2026-11-25).
- **Pilots:** each local stage first runs one step until it is `done`.
- **CheckM:** local CheckM matched KBase's kb_Msuite exactly on 62 bins, on FM and on TF. The values are on each bin Assembly as metadata and in `<unit>.CheckM_local.report`.

## Running

```bash
cd grow/multisample_binning_2026
~/Documents/py_venv/bin/python pipeline.py status
~/Documents/py_venv/bin/python pipeline.py plan --through gtdbtk      # submits nothing
U=grow-sitecov-watch-$(date -u +%Y%m%dt%H%M%Sz)                       # a new scope name for every start
setsid nohup systemd-run --user --scope --unit=$U -p OOMPolicy=continue -p ManagedOOMPreference=avoid \
    ~/Documents/py_venv/bin/python -u pipeline.py watch --through gtdbtk --interval 1800 >> watch.log 2>&1 < /dev/null &
```

- **Finding the watcher:** the per-sample watcher (and the salterns one, if it runs) has the same command line. Find ours with `pgrep -af "[p]ipeline.py watch"`, keeping the PID whose `readlink /proc/<pid>/cwd` is this directory. Never use `pkill -f`/`pgrep -f` with a pattern that also appears in your own command.
- **Restarting:** stop the python PID while it sleeps (`watch.log` written 1–29 min ago). Never stop the scope, which would kill DAS Tool workers and the GTDB-Tk memory guard. Start again with a new scope name, from a shell that has the user bus.
- **One writer:** `ledger.json.lock` allows only one writer.
- **Keep `work/`:** it is needed until GTDB-Tk and `compare.py` have finished.

## History
- 2026-11-18 6:04 PM CST (2026-11-19 00:04Z): the first `advance` submitted the Mobile River pilot binners (MetaBAT2 ×2, CONCOCT ×2, MaxBin2 ×1, the extreme cap of 3) and 6 of the 8 co-assembly merges (njs cap 6). NASQAN2011_239.sitecov went `empty` (no per-sample assembly). The watcher was started in scope `grow-sitecov-watch-20261119t000430z` (PID 4189947).
- 2026-11-24, hybrid mode (user decision 2026-11-23): the watcher was stopped by PID at 9:19 PM CST on Nov 23 and restarted at 8:07 AM CST on Nov 24 (and again at 8:33 AM CST with the TF-routing code, scope `grow-sitecov-watch-20261124t143209z`). Its local steps wait behind the per-sample run (gate priority), so its queued KBase jobs keep their place in the KBase queue until then. Backups: `pipeline.pre_hybrid_20261124.py` and `ledger.pre_hybrid_20261124.json`. Later restarted at 8:44 AM CST (`grow-sitecov-watch-20261124t144430z`).
- 2026-11-24 9:52 AM CST, reserved share (user decision 2026-11-24): this run now has 12 threads / 40 GB of FM reserved (and half of TF), instead of waiting behind the per-sample run. The watcher was restarted by PID, scope `grow-sitecov-watch-20261124t155158z`. The share is applied as per-sample steps finish (none preempted).
- 2026-11-25 6:31 AM CST, poplar overflow (user decision 2026-11-25): the 2 co-assemblies now run only on poplar (MEGAHIT 1.2.9, 32 threads, inputs pulled from KBase there; see `../hybrid/README.md`). Per-sample has priority on poplar. The watcher was restarted by PID, scope `grow-sitecov-watch-20261125t123127z`.
- 2026-11-25 9:12 AM CST, concurrent poplar + KBase fallback (user decision 2026-11-25, second that morning): poplar runs several MEGAHIT jobs at once, and an assembly no local route can take gets a KBase job (fallback queue), taken over while still queued. Both co-assemblies (RedButteCreek, WhiteClayCreek) got KBase MEGAHIT jobs in the first cycle, since per-sample has priority on poplar. The watcher was restarted by PID, scope `grow-sitecov-watch-20261125t151249z`. Backup: `pipeline.pre_kbase_fallback_20261125.py`.
- 2026-11-25 3:36 PM CST, more poplar (user decision 2026-11-25, third): poplar runs up to 5 jobs and stages on /scratch1 (`../hybrid/README.md`). Both co-assemblies (RedButteCreek 24.4 GB gz, WhiteClayCreek 25.8 GB gz merged) are still queued on KBase. They move to poplar after the per-sample queue (per-sample has priority there), unless KBase starts them first. The watcher was restarted by PID, scope `grow-sitecov-watch-20261125t213603z`. Backup: `../hybrid/*.pre_scratch1_20261125.py`.
- 2026-11-25 3:44 PM CST, no duplicated work (user decision 2026-11-25, relayed by the coordinating session; details in `../hybrid/README.md`, "No duplicated work"): one live attempt per task across KBase / FM / TF / poplar. A local start needs EE2's confirmed cancel. KBase (re)submits are blocked while anything is live. A dedup sweep runs every cycle (duplicates, EE2 orphans, jobs the ledger records as finished), and a pre-push check keeps the first of two competing outputs. A tested run on copies of both ledgers found nothing to fix, as did the coordinating session's audit at ~4 PM CST. The watcher was restarted by PID, scope `grow-sitecov-watch-20261125t214432z`. Backups: `../hybrid/{glue,hw}.pre_dedup_20261125.py`, `pipeline.pre_dedup_20261125.py`.
- 2026-11-25 4:04–4:07 PM CST, binners on every route (user decision 2026-11-25, relayed by the coordinating session; details in `../hybrid/README.md`, "Binners on TF"): no binner had run for hours. FM's share was taken by the FM MEGAHIT plus DAS Tool / extract, FM was overloaded by other users, and the binners had no other route. Each binner step now goes to FM within the share, then TF (MaxBin2 / CONCOCT, mapping included), else a KBase job (fallback queue, taken over while queued; never two routes). The first cycle submitted KBase binner jobs up to the pool caps (extreme 3 / bigmem 6 while per-sample is busy). The first MaxBin2 TF pilot started at 4:07 PM CST. The watcher was restarted by PID, scope `grow-sitecov-watch-20261125t220739z`. Backups: `../hybrid/{glue,hw}.pre_tfbin_20261125.py`.
