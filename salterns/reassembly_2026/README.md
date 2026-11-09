# reassembly_2026: re-assembling and re-binning the 39 salterns metagenomes

The current saltern MAGs (310 dereplicated plus 140) were binned in about 2018, one sample at a time, from a Tringe-lab SPAdes reassembly that is not public. `provenance/PROVENANCE.md` traces where they came from.

This directory rebuilds them from the raw JGI reads in a fresh KBase narrative, using the newest KBase app versions. It adds multi-sample coverage, site co-assemblies of the sediment cores, and a local DAS Tool 1.1.7. The resulting MAGs are classified with a local GTDB-Tk 2.8.0 against GTDB R232.

- **Narrative:** ws 275645, `freiburgermsu:narrative_1793847843029`, "Salterns metagenomes: reassembly and binning (2026)".
- **Inputs:** `manifest.tsv`, holding 39 HiSeq 2500 2x151 runs (JGI proposal 1890; 344 Gbp). They were imported from the AWS ODP mirror of SRA.

## Why re-process (decisions of 2026-11-04/05)

- **Binning was per-sample.** The original binning used single-sample coverage and DAS Tool, and KBase's DAS Tool 1.1.2 has a silent scoring bug that was fixed in 1.1.3. Multi-sample coverage within each site is the main expected gain. See `provenance/KBASE_VS_UPSTREAM.md`.
- **Assembly is the ceiling in the cores.** In the SF2 and R2A sediment cores only 14–15 % of reads assemble (`provenance/assembly_vs_bins.tsv`). Co-assembling the six replicate cores of each site is the one in-KBase lever on that ceiling.
- **Where each step runs:**
  - Everything runs in KBase with the latest app versions, except DAS Tool, which runs locally at 1.1.7 (`~/Documents/dastool_env`).
  - GTDB-Tk runs locally, in the established wetland environment (`~/gtdbtk_kbase/env.sh`).
  - DRAM2 annotation of the new MAGs is run separately (KBUtilLib DRAM2Utils on h100).
- **What happens next:** the prFBA augmentation, and every analysis that uses the MAGs, are to be regenerated on the new MAGs. Existing MAG augmentations are not reused.

## Pipeline (`pipeline.py`)

The unit is one assembly: 39 per-sample assemblies plus 4 site co-assemblies `<site>_sediment`, which are binned and refined separately.

| Stage | App (tag), client group | Notes |
|---|---|---|
| import | kb_uploadmethods/import_sra_as_reads_from_web (beta), njs | Gives an interleaved PairedEndLibrary. Both mates carry identical deflines, so pairing is positional. |
| trim | kb_trimmomatic/run_trimmomatic (release), bigmem | Trimmomatic 0.39 PE with ILLUMINACLIP TruSeq3-PE-2 2:30:10, LEADING:3, TRAILING:3, SLIDINGWINDOW:4:15 and MINLEN:50. The output is `<sample>.<srr>.trim_paired`. If trimming fails, the raw import is used. |
| assemble | kb_SPAdes/run_metaSPAdes (beta, 3.15.3), extreme | k 33–127. The run falls back to MEGAHIT 1.2.9 meta-sensitive (bigmem) after an out-of-memory or out-of-disk error, or after 2 failures. |
| merge | kb_ReadsUtilities/KButil_Merge_MultipleReadsLibs_to_OneLibrary (release), njs | Merges the six trimmed cores of each site. The read and base totals must equal the sums of the inputs. |
| coassemble | MEGAHIT/run_megahit (release, 1.2.9), bigmem | Runs meta-large. If MEGAHIT itself fails, the retry uses min-count 2 with k 27..127, and then the MEGAHIT defaults. |
| metabat, maxbin, concoct | metabat/run_metabat (beta), kb_maxbin/run_maxbin2 (release), kb_concoct/run_kb_concoct (release) | Coverage comes from every library of the unit's site (`coverage_group`): R1 9, R2 12 (including the 3 R2_restored samples), R2A 9 and SF2 9. |
| dastool | local DAS Tool 1.1.7 (`local/dastool_worker.py`) | Settings: diamond, score 0.5, duplicate penalty 0.6, megabin penalty 0.5. It saves `<unit>.DASTool.BinnedContigs` to the narrative. |
| extract | MetagenomeUtils/extract_bins_as_assemblies (beta), njs | Turns each refined bin into an Assembly, plus an AssemblySet when there is more than one bin. |
| checkm | kb_Msuite/run_checkM_lineage_wf (beta, CheckM 1.1.2), njs | Runs with the reduced tree. `local/collect_checkm_select.py` then keeps bins ≥ 50 % complete and < 10 % contaminated and writes their FASTA. |
| gtdbtk | local GTDB-Tk 2.8.0 / R232 (`local/run_gtdbtk_reassembly.sh`) | Same settings as the wetland runs: 48 CPUs, 8 pplacer threads, `--min_perc_aa 10`. It runs in a systemd user scope (MemoryMax 170G, no swap) watched by `local/mem_guard.sh`, starts once every CheckM unit has finished and MemAvailable is at least 150 GB, and copies its summaries to `gtdbtk/`. |

**How the driver runs jobs:**
- **Retries and limits:**
  - A failed attempt is retried up to 3 times; cancels do not count.
  - Live KBase jobs are capped per client group (extreme 12, bigmem 20, njs 15), and local DAS Tool workers at 3.
  - The first unit of dastool, extract and checkm runs alone as a pilot before the rest of the stage fans out.
- **Ledger:** every attempt is recorded in `ledger.json` with its exact parameters, tag and service version. The v1 run (imports only) is kept in `ledger.v1_backup.json`.

```bash
cd salterns/reassembly_2026
~/Documents/py_venv/bin/python pipeline.py status
~/Documents/py_venv/bin/python pipeline.py plan --through gtdbtk      # what would start now; submits nothing
U=salterns-reassembly-watch-$(date -u +%Y%m%dt%H%M%Sz)                # a new scope name for every start
setsid nohup systemd-run --user --scope --unit=$U -p OOMPolicy=continue -p ManagedOOMPreference=avoid \
    ~/Documents/py_venv/bin/python -u pipeline.py watch --through gtdbtk --interval 1800 >> watch.log 2>&1 < /dev/null &
```

**Running the watcher:**
- **Its own scope:** the watcher runs in its own systemd user scope with OOMPolicy=continue, so an out-of-memory kill of a DAS Tool child does not take the watcher with it.
  - Start it from a shell that has the user bus (`XDG_RUNTIME_DIR`, `DBUS_SESSION_BUS_ADDRESS`), which GTDB-Tk's `systemd-run`/`systemctl --user` need.
  - Give every start a new scope name: DAS Tool workers outlive a restart inside the old scope, and the new watcher adopts them from the ledger.
- **One writer:** only one `advance` or `watch` can run at a time; `ledger.json.lock` enforces this.
- **Restarting:**
  - Stop the watcher by its PID while it sleeps (`watch.log` written 1–29 min ago). Do not use `pkill -f` with a pattern that also appears in your own command line.
  - Never stop the watcher's scope (`systemctl --user stop`): that would also kill running DAS Tool workers and the GTDB-Tk memory guard. A guard that dies anyway is restarted by the next cycle, and it adopts a freeze left behind.
  - A restart is safe: the state is in `ledger.json`, and the intent to submit is saved before every submit.
- **Assembly floors:** a per-sample assembly must reach 10 Mbp and a site co-assembly 50 Mbp.
  - The site floor was lowered from 100 Mbp because single SF2/R2A cores assemble to only 19–85 Mbp.
  - A smaller output is rejected and retried with the next assembler or setting; it is never adopted.
- **Do not delete `work/` while the run is going:**
  - `work/dastool/<unit>/attempt*/<unit>.save_objects.json` is read until extract has finished.
  - `work/checkm/<unit>/fasta` and the batchfiles are read until GTDB-Tk has finished.
- **GTDB-Tk:** it starts automatically after the last CheckM unit is collected. It uses one batch of all bins ≥ 50 % complete and < 10 % contaminated, and waits for MemAvailable ≥ 150 GB. The run and its `mem_guard` logs go to `~/gtdbtk_kbase/salterns_reassembly_2026/`.

## Files

- `manifest.tsv`: the 39 runs (sample, site, coverage group, layer, SRR, IMG, GOLD, BioSample, size and ODP URL).
- `pipeline.py`, `ledger.json`, `watch.log`: the driver, its state and its log.
- `local/`: local steps and helpers.
  - `kbrpc.py`: read-only KBase RPC.
  - `fetch_assembly.py`: Blobstore FASTA download, verified against the Workspace object.
  - `binnedcontigs_to_contig2bin.py`, `dastool_to_binnedcontigs.py`, `dastool_worker.py`: the local DAS Tool step.
  - `collect_checkm_select.py`: the CheckM collection.
  - `run_gtdbtk_reassembly.sh`, `mem_guard.sh`: the GTDB-Tk run.
  - the exact-parameter records of the new stages.
- `provenance/`: where the old MAGs came from, the KBase-vs-upstream tool versions, read and assembly statistics.
- `work/` (git-ignored): large local intermediates such as assembly FASTA, DAS Tool bins and CheckM tables.
- `gtdbtk/`: the GTDB-Tk R232 summaries of the new MAGs, once the run finishes.
