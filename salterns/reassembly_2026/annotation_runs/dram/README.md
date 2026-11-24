# Track `dram`: DRAM 1.5.0 on the 558 best-per-species re-assembly MAGs (poplar)

*Started 2026-11-11 by a Claude subagent. DRAM2 on h100 isn't reachable, so this runs DRAM 1.5.0 on poplar instead, as `annotation_choice/README.md` §6.2 allows ("if the wetland project stays on DRAM 1.5, use DRAM 1.5 here too, with `--min_contig_size 1000` and the `.dram1.tsv` sheet").*

**Status (2026-11-13): DONE.** `state/DONE` at 02:15 CST; 558/558 MAGs annotated, merged (1,645,407 genes) and distilled. The results were copied to `results/` and md5-checked. They were validated against the KOfam track: identical gene calls and marker KOs, with only the Sec-rescue *grd* differences. The manuscript summary is in [`RESULTS.md`](RESULTS.md). The launch-time status below is kept for the record.

**Status at the end of this session (2026-11-11 12:20 CST):** the environment is built, the 558 FASTAs are on poplar (md5-checked), all raw database files are downloaded and verified, every database except Pfam is built, and the Pfam build is running. The detached driver on poplar (PID 4007531) finishes the Pfam build, then runs a one-MAG test, annotates all 558 MAGs, merges the results and runs distill. Check progress with `status.sh` (below). The **Results** section at the end is what was known when this session ended.

## Where things are

| what | where |
|---|---|
| work dir (poplar) | `/scratch1/afreiburger/salterns_2026/dram/` (`$W`) |
| DRAM env | `/scratch1/afreiburger/salterns_2026/envs/dram15` |
| DRAM databases | `/scratch1/afreiburger/salterns_2026/dram_db/` |
| raw DB downloads (kept; md5 and Last-Modified in `MANIFEST.tsv`) | `/scratch1/afreiburger/salterns_2026/dram_db/raw/` |
| MAG FASTAs (558, `<name>.fna`; md5 list `$W/mags558.md5`) | `$W/mags/` |
| inputs: MAG list, taxonomy, quality, custom sheet | `$W/inputs/` (copies of the files in this folder) |
| per-MAG DRAM output | `$W/dram_out/<MAG>/` (`annotations.tsv`, `genes.faa/fna/gff`, `trnas.tsv`, `rrnas.tsv`, `genbank/`, `scaffolds.fna`) |
| merged / distilled | `$W/merged/`, `$W/distill/` (`metabolism_summary.xlsx`, `product.html`, `product.tsv`, `genome_stats.tsv`), `$W/distill_gene_names/` |
| one-MAG test distill | `$W/test_distill/` |
| logs | `$W/logs/`: `driver.log`, `fetch_dbs.log`, `verify_nopfam.log`, `setup_dbs.nopfam.log`, `setup_dbs.pfam.log`, `setup_check.log`, `test_check.log`, `run_all.log` (one line per MAG), `mags/<MAG>.log` (DRAM log + `/usr/bin/time -v`), `merge*.log`, `distill*.log` |
| stage markers | `$W/state/`: `SETUP_A_OK`, `SETUP_B_STARTED`, `SETUP_OK`, `TEST_OK`, `ANNOTATE_OK`, `MERGE_OK`, `DISTILL_OK`, **`DONE`**; on trouble `FAILED` (one line: what failed and which log) or `DONE_INCOMPLETE` (list of MAGs without output) |
| env record | `$W/env/versions.txt`, `$W/env/dram15.explicit.txt`, `$W/env/dram_config.txt` (`DRAM-setup.py print_config`, written after setup) |
| scripts (same as `scripts/` here) | `$W/scripts/` |

Local files in this folder:
- `best_per_species_558.tsv`: the 558 MAGs with unit, species cluster, completeness, contamination, score, `tie_in_cluster`, GTDB R232 genus/species, size, N50 and local FASTA path.
- `gtdb_taxonomy_558.tsv` (`user_genome`, `classification`; identical to `reassembly_2026/gtdbtk/` R232 for all 558) and `checkm_quality_558.tsv` (`Bin Id`, `Completeness`, `Contamination`; KBase CheckM, from `comparison/mags_old_vs_new.tsv`). These are DRAM's `--gtdb_taxonomy` / `--checkm_quality` inputs.
- `scripts/`: everything run on poplar.

## How to check status

```bash
ssh -o BatchMode=yes poplar 'bash /scratch1/afreiburger/salterns_2026/dram/scripts/status.sh'
# details
ssh -o BatchMode=yes poplar 'cd /scratch1/afreiburger/salterns_2026/dram; cat logs/driver.log; ls state; tail logs/run_all.log'
```

`status.sh` prints the stage markers, whether the driver PID is alive, MAGs done/failed, running DRAM processes, `uptime` and `free -g`.

**Processes started this session (poplar):**
- `fetch_dbs.sh`: PID 3693645. **Finished:** `dram_db/raw/FETCH_OK` at 12:15; log `logs/fetch_dbs.log`.
- `driver.sh`: PID **4007531** (also in `$W/driver.pid`), started 2026-11-11 11:04 CST with `setsid -f bash scripts/driver.sh > logs/driver.log`. It is the only thing that needs to keep running. The finish marker is `$W/state/DONE`.

Both are detached (`setsid`, stdin from `/dev/null`), so a dropped ssh or tmux session doesn't stop them.

## What the driver does (`scripts/driver.sh`)

Each stage is skipped if its marker exists, so re-running the driver resumes where it stopped.

0. Writes the tool versions and an explicit conda list to `$W/env/`.
1. **Setup phase A** (everything except Pfam) as soon as the non-Pfam raw files are present. First `verify_nopfam.sh` checks them:
   - the KOfam `ko_list` md5 is the 2026-10-26 release;
   - the VOG md5s match the server's `.md5` files;
   - gzip integrity;
   - the dbCAN file is HMMER3 and the CAZy tables are not HTML;
   - `pepunit.lib` is FASTA.

   Then it runs `setup_dbs.py ... nopfam`.
2. **Setup phase B** (Pfam) after `FETCH_OK`. `fetch_dbs.sh` also checks the Pfam md5 against EBI's `md5_checksums`. Then `setup_dbs.py ... pfam`, `DRAM-setup.py print_config`, and a check that every configured path exists (`logs/setup_check.log`).
3. **Test** on one MAG, `concoct.bin.101.fasta__R1_sediment`: the best new *Methanolobus* (98.4/2.6), chosen because it should show *mcrA* and the Pyl methylamine genes.
   - Runs `DRAM.py annotate`, then `DRAM.py distill` with the custom sheet, into `$W/test_distill`.
   - `check_test.py` requires:
     - the KO, Pfam, CAZy, MEROPS and taxonomy/quality columns;
     - ≥20 % of genes with a KO and ≥20 % with Pfam;
     - *mcrA* present;
     - no HTML in the description columns;
     - GTDB taxonomy and completeness > 0 (DRAM silently substitutes the bin name, 0 and 100 on a name mismatch);
     - a 65-row "Salterns methane" sheet in `metabolism_summary.xlsx`;
     - `product.html`/`product.tsv`/`genome_stats.tsv`.

   The driver stops with `state/FAILED` if any check fails.
4. **All 558 MAGs:** one `DRAM.py annotate` process per MAG, 8 at a time × 6 threads (≤48 threads), `nice -n 10`. The `xargs` pass runs twice, so failures get one retry; finished MAGs are skipped.
   - Why per MAG rather than 8 batches of ~70: DRAM 1 annotates the FASTAs of one call one after another and writes nothing until the end. One bad MAG in a 70-MAG batch would lose a day of work, and a restart couldn't resume. The CPU use is the same.
5. **Merge:** `DRAM.py merge_annotations -i "$W/dram_out/*" -o $W/merged`. Checks that the number of distinct `fasta` values equals the number of per-MAG outputs.
6. **Distill:** `DRAM.py distill -i merged/annotations.tsv -o distill --trna_path ... --rrna_path ... --custom_distillate inputs/custom_distillate_salterns.dram1.tsv`. Then a second, non-fatal run with `--distillate_gene_names` into `distill_gene_names/`, which lists gene IDs instead of counts and is useful for checking marker hits.
7. Writes `state/DONE`, or `state/DONE_INCOMPLETE` listing the MAGs without output.

## Exact commands

**Environment** (`$W/make_env.sh`, then the fix):
```bash
export CONDA_PKGS_DIRS=/scratch1/afreiburger/salterns_2026/envs/pkgs
nice -n 10 /scratch1/afreiburger/metag-hybrid/env/miniforge/bin/mamba create -y -p /scratch1/afreiburger/salterns_2026/envs/dram15 \
  -c conda-forge -c bioconda "dram=1.5.0" "python=3.10" "pandas>=1.5,<2" "sqlalchemy<2"
# fix 1: setuptools 84 has no pkg_resources, which DRAM 1.5.0's database_handler imports
mamba install -y -p .../envs/dram15 -c conda-forge -c bioconda "setuptools<81" "pandas>=1.5,<2" "sqlalchemy<2" "dram=1.5.0"
# fix 2: always run with PYTHONNOUSERSITE=1 (scripts/dram_env.sh); otherwise ~/.local/lib/python3.10 on poplar injects another requests/urllib3
```

**Versions** (`$W/env/versions.txt`):

| component | version |
|---|---|
| DRAM | 1.5.0 (bioconda) |
| Python | 3.10.22 |
| pandas | 1.5.3 |
| SQLAlchemy | 1.4.54 |
| numpy | 1.26.4 |
| scikit-bio | 0.5.8 |
| setuptools | 80.10.2 |
| Prodigal | 2.6.3 |
| MMseqs2 | 18.8cc5c |
| HMMER | 3.4 |
| tRNAscan-SE | 2.0.13 |
| barrnap | 1.10.5 |

**Per-MAG annotate** (`scripts/run_one_mag.sh`):
```bash
nice -n 10 /usr/bin/time -v DRAM.py annotate -i mags/<MAG>.fna -o dram_out/<MAG> --min_contig_size 1000 --threads 6 \
  --gtdb_taxonomy inputs/gtdb_taxonomy_558.tsv --checkm_quality inputs/checkm_quality_558.tsv
```
Everything else is DRAM's default:
- prodigal `meta`, translation table 11;
- bit score 60, RBH bit score 350;
- KOfam with KOfam's own thresholds (`sig_scores`: `-` thresholds are skipped, so *mpnS*, *phnY\** and *phnZ1* can't be called; see `annotation_choice/README.md` §2a);
- Pfam (MMseqs2 profiles), dbCAN, MEROPS;
- no KEGG Genes, no UniRef90;
- **no VOGDB search**: `--use_vogdb` is opt-in in DRAM 1.5 and is not passed. VOGDB is still built, so it can be added later.

`--low_mem_mode` is not passed; without KEGG Genes and UniRef it changes nothing.

## Databases: one release of each, recorded

All files were downloaded 2026-11-11 to `dram_db/raw/`; see `MANIFEST.tsv` there for URL, Last-Modified, bytes and md5.

| DB | release / file | source used | note |
|---|---|---|---|
| **KOfam** | `ko_list.gz` **2026-10-26** (md5 `697af5b3b22a919306a1e60f99cec1a6`, identical to `annotation_choice/sources/kegg/kofam_ko_list_2026-10-26.gz`); `profiles.tar.gz` Last-Modified **2026-10-27** (1,567,491,980 B, md5 `bc68ad2a52ed48812b5496d15c215a7d`, 27,981 profiles) | `https://www.genome.jp/ftp/db/kofam/` | **Use this pair for every KOfam call in the paper.** It is byte-identical to the KOfam track's `/scratch1/afreiburger/salterns_2026/kofam_db/{ko_list.gz,profiles.tar.gz}` (both md5s checked 2026-11-11), so DRAM and the 1,324-MAG KOfam run share one release. The `annotation_choice/` control test used the 2026-05-28 build; thresholds moved between them (e.g. *phnJ* 246.6 domain → 194.6 full). |
| Pfam | **39.0** (2026-09, UniProtKB 2026_02), `Pfam-A.full.gz` 40.5 GB, md5 `1ca35677412407da16a340bdd293e773`; `Pfam-A.hmm.dat.gz` md5 `13c095351722642b3c51b56d7ecf9a0a` | `ftp.ebi.ac.uk/pub/databases/Pfam/current_release/` | DRAM builds an MMseqs2 profile DB from the full alignments (`convertmsa`, `msa2profile --match-mode 1`, `createindex -k 5 -s 7`) |
| dbCAN | **V11**, `dbCAN-HMMdb-V11.txt` + `CAZyDB.08062022.fam-activities.txt` + `CAZyDB.08062022.fam.subfam.ec.txt` (DRAM 1.5 defaults) | `https://pro.unl.edu/dbCAN2/download/Databases/V11/` | DRAM's `bcb.unl.edu` URLs now redirect to the dbCAN2 home page (HTML with HTTP 200), so DRAM would have stored HTML as these files |
| MEROPS | `pepunit.lib`, Last-Modified 2023-02-22 (MEROPS 12.x current_release) | `ftp.ebi.ac.uk/pub/databases/merops/current_release/` | |
| RefSeq viral | `viral.1.protein.faa.gz`, 2026-09-04 (DRAM downloads only file 1) | NCBI | not used by `annotate` for MAGs |
| VOGDB | release **237** (`vog.hmm.tar.gz` 2026-09-09) + `vog.annotations.tsv.gz`, md5-checked | `fileshare.lisc.univie.ac.at/vog/latest/` (DRAM's `csb.univie.ac.at` redirects there) | built but not searched (no `--use_vogdb`) |
| DRAM sheets | genome summary, module step, function heatmap, AMG, ETC forms | downloaded by DRAM from GitHub `WrightonLabCSU/DRAM` `master` at setup time | paths in `env/dram_config.txt` |
| UniRef90, KEGG Genes | not used | | |

**How the setup differs from a plain `DRAM-setup.py prepare_databases --skip_uniref`.** Each change was needed for the setup to finish or give correct files:
1. **Raw files are downloaded first** by `fetch_dbs.sh`, using 8 parallel HTTP range requests per file, resumable. Each file is checked for size and, where the source publishes one, md5. The files are then passed to DRAM. Single-stream speeds were 0.4 MB/s from genome.jp and 4.6 MB/s from EBI; Pfam alone is 40.5 GB.
2. **`setup_dbs.py` calls `mag_annotator.database_processing.prepare_databases()` directly**, because:
   - `prepare_databases()` honours only arguments ending in `_loc`, so `dbcan_fam_activities`/`dbcan_subfam_ec` are always downloaded from the dead bcb.unl.edu URLs. `download_dbcan_fam_activities` and `download_dbcan_subfam_ec` are replaced by copies of the verified pro.unl.edu files, and `download_vog_annotations` by the md5-checked file.
   - **VOG release 237 puts the HMMs under `hmm/`**, while DRAM's `process_vogdb` globs only the top level. It would have merged nothing and failed `hmmpress` after the Pfam build. `process_vogdb` is replaced by the same steps with a recursive glob, asserting more than 1,000 models.
   - **Two phases** through DRAM's own `select_db`:
     - phase A builds everything except Pfam while Pfam downloads, so any failure there shows up early;
     - phase B adds Pfam to the same config and description database. `populate_description_db` clears and refills only the selected tables.
3. `--threads 24` for MMseqs2 during setup.

## The 558 MAGs

- One per `species_cluster` among the 1,324 `new|` rows of `comparison/mags_old_vs_new.tsv`, maximising completeness − 5 × contamination.
- 511 bacteria and 47 archaea, 1,622 Mb.
- **Ties:** 4 clusters have an exact score tie (`tie_in_cluster` = 2).
  - The list uses the existing `best_of_cluster` flag in `annotation_choice/short_contig_fraction_per_mag.tsv`, so these are the same 558 MAGs the annotation write-up counts. Every flagged MAG has its cluster's maximum score.
  - An N50 tie-break would pick a different member in 3 of the 4 tied clusters:
    - `maxbin2.Bin.037.fasta__R2_sediment` instead of `metabat2.bin.198.fasta__R1_sediment`;
    - `metabat2.bin.265.fasta__R2A_sediment` instead of `metabat2.bin.004.fasta__R2A_B_D2`;
    - `metabat2.bin.153.fasta__R1_sediment` instead of `metabat2.bin.081.fasta__R2_sediment`.

## Custom distillate (`annotation_choice/custom_distillate_salterns.dram1.tsv`): validated, unchanged

Checked against DRAM 1.5.0's `summarize_genomes.py` in this env:
- `summarize_genomes()` concatenates the sheet onto `genome_summary_form` and drops `potential_amg`.
- `fill_genome_summary_frame()` matches `gene_id` exactly (or as `<id>.`) against IDs parsed from `ko_id`, `kegg_id`, `kegg_hit` EC, `pfam_hits`, `cazy_best_hit` and `peptidase_family`.
- `write_summarized_genomes_to_xlsx()` writes one Excel tab per `sheet`, sorted by header, subheader, module and gene_id.

The sheet has exactly the required columns (`gene_id, gene_description, module, sheet, header, subheader, potential_amg`), 65 unique KO IDs, and no blanks. All 65 KOs exist in the 2026-10-26 `ko_list` with a real threshold, K28979 (*tdm*) included.

A synthetic run through DRAM's own functions produced a "Salterns methane" tab with the right counts (*phnJ* 2, *grdH* 1, *mttB* 1 in one genome; *mcrA* 1, *mtbB* 1 in the other). **No change was needed.** The real test-MAG distill repeats this check.

The sheet affects **only `metabolism_summary.xlsx`** (the new "Salterns methane" tab). `product.html`/`product.tsv` come from DRAM's fixed module, ETC and function-heatmap forms. So the K14083 (*mttB*) over-call in Rhodobacteraceae will still light "trimethylamine => dimethylamine" in the product heatmap. Apply the *mcrA*/Pyl rule from `annotation_choice/README.md` §6.3 when reading it.

## Expected time

- Downloads: ~1.5 h, of which Pfam is most.
- Setup:
  - phase A: under 1 h (KOfam/VOG `hmmpress`, MEROPS/viral MMseqs2 DBs, description DB);
  - phase B: several hours for Pfam 39 (`convertmsa` on 40.5 GB gz is single-threaded, then `msa2profile` and `createindex`).
- Annotation: KOfam dominates at about 0.6 CPU-s per protein. For ~1.66 M genes that is roughly 280 CPU-h plus Pfam, dbCAN, MEROPS and tRNAscan. At 48 threads on a host at load ~420/336 cores, expect **~1.5–3 days**.
- Merge and distill: 1–3 h.

## If something fails

- **`state/FAILED` exists:** read the log it names.
- **Phase A failed:** `rm -rf /scratch1/afreiburger/salterns_2026/dram_db/*` except `raw/`, using `find ... -mindepth 1 -maxdepth 1 ! -name raw -exec rm -rf {} +`. Then fix the cause and re-run the driver.
- **Phase B (Pfam) failed:** remove `dram_db/database_files`, `dram_db/pfam.*`, `dram_db/tmp` and `$W/state/SETUP_B_STARTED`, then re-run the driver. Phase A is kept.
- **Re-run the driver:**
  ```bash
  cd /scratch1/afreiburger/salterns_2026/dram && setsid -f bash scripts/driver.sh >> logs/driver.log 2>&1 < /dev/null; pgrep -fx 'bash scripts/driver.sh' > driver.pid
  ```
  Finished stages and MAGs are skipped.
- **A MAG that keeps failing:** read `logs/mags/<MAG>.log`. The driver continues without it and lists it in `state/DONE_INCOMPLETE`.
- **Stop everything:** `kill $(cat $W/driver.pid)`, then `pkill -u afreiburger -f 'DRAM.py annotate'`.

## Results so far (2026-11-11 12:20 CST, end of this session)

| step | state | when / numbers |
|---|---|---|
| env `dram15` | done | DRAM 1.5.0 works (`DRAM-setup.py version` → 1.5.0) after the setuptools<81 pin and `PYTHONNOUSERSITE=1` |
| 558 FASTAs on poplar | done | `md5sum -c` OK for all 558 |
| custom sheet vs DRAM 1.5 distill code | done | valid as-is (synthetic run through DRAM's functions) |
| raw DB downloads | **done**, `FETCH_OK` 12:15 | Pfam md5 OK vs EBI; Pfam 39.0 (53,275 families); KOfam 27,981 profiles; VOG 49,128; dbCAN 699 models |
| setup phase A (no Pfam) | **done**, `SETUP_A_OK` 12:17 | 70 min, mostly NFS small-file I/O for the KOfam and VOG merges. Description DB: dbCAN 660 families with activities/EC (e.g. `GH51 endoglucanase (EC 3.2.1.4)…`), MEROPS 1,227,939, RefSeq viral 722,738, VOG 49,128 |
| setup phase B (Pfam) | **running**, started 12:17 | At 13:11, `mmseqs convertmsa` had written 200 GB of `pfam.mmsmsa.0` (~57 MB/s) and indexed ~37k of 53,275 families, so it should finish around 13:40. Expect `convertmsa` ~1.5–2 h, then `msa2profile` + `createindex` (24 threads), so ~3–6 h in total. `dram_db/pfam.mmsmsa.*` will be a few hundred GB; it isn't needed after setup and can be deleted once the run is finished |
| one-MAG test, 558 annotations, merge, distill | queued in the driver | runs unattended. Watch for `state/TEST_OK`, then `logs/run_all.log` (expect ~1.5–3 days), then `state/DONE` |

**Next steps for whoever picks this up:**
1. Confirm `state/TEST_OK` and read `logs/test_check.log`: fraction of genes with a KO, *mcrA*/*mtbB*/*mtmB*/*mttB*/*pylB* calls, and the "Salterns methane" row counts for the *Methanolobus* test MAG. Then read the per-MAG time and peak RSS at the bottom of that file.
2. When `state/DONE` appears, copy back `merged/annotations.tsv` (gzip it), `merged/trnas.tsv`, `merged/rrnas.tsv`, `distill/`, `distill_gene_names/metabolism_summary.xlsx`, `env/` and `logs/` to `annotation_runs/dram/results/` here:
   ```bash
   rsync -a -e "ssh -o BatchMode=yes" poplar:/scratch1/afreiburger/salterns_2026/dram/{distill,env} results/
   ```
3. **Watch item:** DRAM's description DB is SQLite on NFS (`dram_db/description_db.sqlite`). The one-MAG test runs alone, so the first time 8 `DRAM.py annotate` processes read it together is the main run. Once `logs/run_all.log` shows the first MAGs finishing, run `grep -l "database is locked\|disk I/O error" logs/mags/*.log`. If anything shows up, lower `P` in `scripts/dram_env.sh` and re-run the driver. The second `xargs` pass retries failed MAGs, so nothing is lost.
4. When reading `product.html`, apply the guild caveats from `annotation_choice/README.md` §3/§6.3: K14083 *mttB* in non-methanogens, Sec-truncated *grdH*, Pyl fragments.

## Review (2026-11-11, second agent, from 13:19 CST)

*A reviewer checked this track before the long compute. This section records what was checked, the results, and anything changed. Nothing above it was edited.*

### Your tmux message: "poplar connection ended" didn't mean the connection ended

- The `poplar-conn` tmux pane runs `ssh -MN -o ControlPersist=yes ... poplar; echo 'poplar connection ended'; exec bash`.
- With `ControlPersist=yes`, ssh moves the master connection into the background once you've logged in. The foreground command then returns, so the `echo` prints "poplar connection ended" every time, including when the login worked.
- The master is still up: `ssh -O check poplar` reports `Master running (pid=2604712)`, started 10:27 together with the tmux session. Every BatchMode ssh call by the agents goes through it.
- The poplar jobs don't depend on it anyway. The driver (PID 4007531) has PPID 1, its own session and no terminal, and its child `setup_dbs.py pfam` is running at nice 10. A dropped tmux or ssh session doesn't stop either.

### Checks that passed (nothing changed)

| check | how | result |
|---|---|---|
| scripts on poplar = local `scripts/` | `md5sum scripts/*` on both hosts | all 8 identical. `driver.sh` was last written 11:04:36, 7 s before the driver started, so the running bash isn't reading an edited file. Only `status.sh` changed later (13:14), and it isn't a running file |
| the 558 MAGs | recomputed from `comparison/mags_old_vs_new.tsv` (1,324 `new\|` rows, max completeness − 5 × contamination per `species_cluster`) | 558 clusters with one MAG each. Every listed MAG has its cluster's top score. The 4 tie clusters are as described. The list is exactly the `best_of_cluster` set in `short_contig_fraction_per_mag.tsv`. Completeness, contamination and unit match the source |
| FASTAs on poplar | md5 of all 558 local `work/checkm/<unit>/fasta/<name>.fna` vs `$W/mags/*.fna` | identical (digest of the sorted md5 list `f9c10643…` on both hosts) |
| DRAM inputs | md5 of `inputs/*` vs the files here | identical |
| taxonomy / quality joins | DRAM 1.5 matches `get_fasta_name()` (basename minus `.fna`) against column 1 of each table; the fallbacks are the bin name, or completeness 0 and contamination 100 | all 558 names are in both tables and unique, with no special characters. All 558 classifications are byte-identical to `reassembly_2026/gtdbtk/gtdbtk.{bac120,ar53}.summary.tsv` |
| KOfam release | md5 | `dram_db/kofam_ko_list.tsv` = KOfam track's `kofam_db/ko_list` (`5dca1ddf…`). `raw/` `ko_list.gz` and `profiles.tar.gz` md5s = KOfam track's `md5sums.txt`. `ko_list.gz` = `annotation_choice/sources/kegg/kofam_ko_list_2026-10-26.gz` |
| KOfam profiles vs `ko_list` | `grep ^NAME dram_db/kofam_profiles.hmm` vs `ko_list` column 1 | 27,981 profiles, no duplicates, **every profile has a `ko_list` row**, so `sig_scores`' `score_db.loc[target_id]` can't raise `KeyError` mid-run. 523 `ko_list` rows have no profile (harmless). `score_type` values are only `full` 17,950, `domain` 8,684 and `-` 1,870, so no `ValueError` |
| custom sheet KOs | 65 `gene_id`s vs the 2026-10-26 `ko_list` | all present, none with a `-` threshold (K28979 *tdm*: 1,048.07 full) |
| KO-calling rule | read DRAM 1.5.0 as installed (`utils.py` `sig_scores`, `annotate_bins.py` `kofam_hmmscan_formater`, `run_hmmscan`) | `hmmsearch --domtblout --cpu T` with default reporting; a hit passes if score **>** threshold, using the full or domain score as `score_type` says; `-` KOs are skipped; best KO per gene = lowest full E-value. This is the same rule as `annotation_choice/kofam_control_test/run_kofam_controls.py` (strict `>`) and the KOfam track's `postprocess.py` (strict `>`, best KO by E-value) |
| shared config not rewritten by `annotate` | grep of `annotate_bins.py`/`database_handler.py` | `annotate` never calls `write_config`; `filter_db_locs` only edits the in-memory copy, so 8 parallel annotates can't corrupt `CONFIG` |
| phase B leaves phase A intact | read `prepare_databases` / `populate_description_db` | `clear_config=False` keeps the phase-A paths. The description DB refill is limited to `select_db`, so only `pfam_description` is cleared and refilled |
| merge with failed MAGs | read `merge_annotations_cmd` | directories without `annotations.tsv` are skipped with a warning, not an error. The driver's merge check counts only directories that have one, so a twice-failed MAG doesn't block merge/distill |
| SQLite description DB on NFS (v4.2) | `review/sqlite_conc.py`: 8 processes × 18,024 DRAM-style `get_descriptions` lookups (peptidase, viral, dbCAN) at once, through DRAM's own `DatabaseHandler` | no lock or I/O errors; 22–61 s per process (about 2–3 ms per ID). The watch item above is now low risk for 8 jobs |
| driver detached | `ps` | PID 4007531: PPID 1, own session (SID = PID), no tty. Child chain `time` → `setup_dbs.py pfam` → `mmseqs` at nice 10 |
| Pfam build | `wc -l dram_db/pfam.mmsmsa.index` | `convertmsa` finished at about 13:24 with all **53,275** Pfam 39.0 families (244 GB `pfam.mmsmsa`). `msa2profile` (24 threads) was running at 13:35. The whole phase B finished at 14:07 (see the status section below) |
| thread use | `ps` | this track: 1 thread (`convertmsa`), then 24 (`msa2profile`), plus 6 for the review smoke test. The ~40 single-thread `hmmsearch` chunks belong to the KOfam track, not this one |

### Measured: DRAM's per-MAG steps use about 1.1 cores whatever `--threads` is

- **Smoke annotate** of the test MAG (2,279 genes at ≥1,000 bp; `review/annotate.log`), all databases except Pfam:

  | step | wall time |
  |---|---|
  | prodigal + mmseqs gene DB | 30 s |
  | **KOfam `hmmsearch`** | **24 min** (13:29:44 → 13:53:59) |
  | MEROPS forward + reverse best hits + descriptions | 5 min |
  | dbCAN | 14 s |

  Totals: 30 min wall, 1,944 CPU-s (107 % CPU), **peak RSS 2.5 GB**. It then stopped at the merge step for a reason specific to running without Pfam (see the next section).
- **`hmmsearch` thread scaling** (`review/bench.sh`): 1,000 KOfam profiles vs the same 2,279 proteins, nice 10, at load ~550.

  | `--cpu` | wall | CPU | %CPU | domtblout lines |
  |---|---|---|---|---|
  | 1 | 79 s | 60 s | 75 % | 5,488 |
  | 3 | 64 s | 67 s | 105 % | 5,488 |
  | 6 | 48 s | 67 s | 139 % | 5,488 |

  The CPU time is constant and the output identical. Each profile is a short job against a one-MAG target set, on an oversubscribed host, so extra threads give little.
- **Estimate.** At ~1,022 genes/Mb the average MAG (2.9 Mb, ~3,000 genes) needs about 30 min of KOfam and 6 min of MEROPS. Pfam (MMSeqs2, which does use its threads) is still unmeasured, so call it ~45–60 min per MAG. At P = 8 that is 558 × ~50 min / 8 ≈ **2–3 days**, in line with the estimate above. Only ~10 of the 48 allocated threads will be busy during the KOfam steps.
- **Option, not applied:** `P=16 THREADS=3` keeps the same 48 nominal threads and should be about 1.5× faster in the KOfam phase: (16 / 64 s) ÷ (8 / 48 s), from single runs on a host with a swinging load.
  - The cost: peak RSS is 24 GB per job in the Pfam step (measured in the driver's test, below), so 16 jobs could reach ~380 GB if they peak together; and 16 SQLite readers instead of 8.
  - Why it wasn't applied: P is read once, when the driver starts, so changing it means killing and relaunching the driver. That isn't justified for speed alone.
  - **To apply it now that annotation is running.** In-flight MAGs are lost and redone, up to ~45 min each.
    1. Edit `scripts/dram_env.sh` to `THREADS=3 P=16`. Write a new file and `mv` it into place.
    2. `kill 4007531 3831773` (driver and its `xargs`), then `pkill -u afreiburger -f 'DRAM.py annotate'`, then make sure `pgrep -af 'DRAM.py annotate'` is empty. Otherwise a relaunched run could annotate the same MAG twice into one directory.
    3. Relaunch the driver:
       ```bash
       cd /scratch1/afreiburger/salterns_2026/dram && nohup setsid bash scripts/driver.sh >> logs/driver.log 2>&1 < /dev/null & echo $! > driver.pid
       ```
       Finished MAGs are skipped, and `run_one_mag.sh` deletes and redoes partial ones.
    4. Watch `free -g` while 16 jobs are in the Pfam step.

### Comparability with the KOfam track (1,324 MAGs)

- **Same:** release, threshold rule (strict `>`, full or domain score per `score_type`, `-` never called), and best KO per gene by lowest full E-value.
- **E-values:** they differ in scale (DRAM: Z = genes in one MAG; KOfam track: `-Z 3911146`), but within one gene and one search Z is the same for every profile, so the E-value ranking of KOs is the same.
- **One real difference, E = 0 ties:** when several KOs reach E = 0 for one gene, DRAM takes the first domtblout row (`kofam_hmmscan_formater`: `best_hit["target_id"].iloc[0]`), which is profile order in `kofam_profiles.hmm`. The KOfam track breaks ties by the higher full score. This only affects very strong hits to several KOs, so a disagreement there isn't a threshold bug.
- **Contig minimum:** both use ≥1,000 bp contigs, but DRAM runs its own prodigal (`meta`, table 11) per MAG, so gene IDs differ (`<MAG>_<contig>_<n>`).

### Minor notes

- `barrnap --version` prints 1.10.5, while the env's explicit list has the `barrnap-1.10.6` package. The binary reports the older string.
- The env has **altair 6.2.2**, while DRAM 1.5.0 was written for altair 4. `product.html` is checked by the smoke distill below.
- `status.sh` counts any `DRAM.py annotate` of this user, so it showed the review smoke test while that ran.
- The test distill (driver line 78) and `merge_annotations` (line 102) run without `nice`. Both are short and single-threaded, so this isn't worth relaunching the driver.

### Distill works with altair 6.2.2 (synthetic table)

- `review/synth_distill.sh` runs `DRAM.py distill --custom_distillate inputs/custom_distillate_salterns.dram1.tsv` on a 2-genome, 14-gene table: an archaeon gA with *mcrABG*, *mtbB*, *mtmB* and *mttB*, and a bacterium gB with *phnGHIJ*, *grdI* and *mttB*.
- Result: exit code 0, and `genome_stats.tsv`, `metabolism_summary.xlsx`, `product.html` (vega-lite 6.4.1 embedded) and `product.tsv` were all written.
- `metabolism_summary.xlsx` has the 6 default tabs plus **"Salterns methane"**, 65 rows. The counts are exactly the KOs put in: gA *mcrA/B/G*, *mtbB*, *mtmB*, *mttB* = 1; gB *phnG/H/I/J*, *grdI*, *mttB* = 1.
- `product.tsv` lights "trimethylamine => dimethylamine" from K14083 alone, which is the *mttB* caveat above.
- DRAM 1.5.0 has no `methyl` database for `annotate` (distill only parses a `methyl_id` column), so "methyl if installed" from `annotation_choice/README.md` §6.2 can't apply here.

### Smoke test 1 (no Pfam) stopped at merge, for a reason specific to the test

- **What ran:** `review/smoke.sh`, i.e. `DRAM.py annotate` on `concoct.bin.101.fasta__R1_sediment` with exactly the flags of `run_one_mag.sh`, against the live config, where Pfam is still `null`. Output went to `review/`.
- **What happened:** prodigal, KOfam, MEROPS and dbCAN all ran. DRAM then failed in `generate_annotated_fasta` (`annotate_bins.py` line 485) with `AttributeError: 'Series' object has no attribute 'pfam_hits'`.
- **Why it doesn't affect the real run:**
  - DRAM grades a gene **D** when it has a Pfam, CAZy *or* MEROPS hit. For a D gene it then writes `annotation.pfam_hits` into the FASTA header, and with Pfam unconfigured that column doesn't exist.
  - With Pfam configured, every MAG has Pfam hits, so the column always exists.
  - It does mean **DRAM 1.5.0 can't run with Pfam switched off**. Nobody should try that as a shortcut.
- **A bug in my script, not the track's:** `smoke.sh` printed `rc=0` for each step because `$?` was read after a `$(date)` substitution. The real exit code was 1 (`time -v` says "Exit status: 1"). The track's `run_one_mag.sh` captures `rc=$?` correctly, and `e2e.sh` below does too.
- Because of that failure, tRNAscan, barrnap, the GFF/GenBank writers, the taxonomy and quality join, `annotations.tsv`, distill and `check_test.py` hadn't been exercised on real output. Smoke test 2 covers them.

### Smoke test 2: end to end through every annotate step, check_test passes except for Pfam coverage (expected)

- **Setup (`review/e2e.sh`, 14:01–14:33):**
  - A mini Pfam profile DB: 3,000 random Pfam 39.0 families under 3 MB, cut from the finished `dram_db/pfam.mmsmsa` with `mmseqs createsubdb`, then DRAM's own `msa2profile --match-mode 1` and `createindex -k 5 -s 7`.
  - A private copy of `description_db.sqlite`, with Pfam descriptions filled by DRAM's own `populate_description_db(select_db=['pfam'])`.
  - A private config, passed with `--config_loc`.
  - Then `DRAM.py annotate` with the `run_one_mag.sh` flags on the test MAG, `DRAM.py distill --custom_distillate`, and `scripts/check_test.py`.
  - The live config, description DB and work directories weren't touched.
- **Annotate:** exit 0 in 29:53 wall, 110 % CPU, **peak RSS 3.0 GB**.
  - Steps: prodigal 14 s, KOfam 25.5 min, MEROPS 5.5 min, Pfam (mini) 35 s, dbCAN 24 s, merge plus FASTA/GFF/GenBank 1 min 15 s, tRNAscan-SE and barrnap within the final minute.
  - Wrote `annotations.tsv`, `genes.faa/fna/gff`, `scaffolds.fna`, `genbank/`, `trnas.tsv` and `rrnas.tsv`.
- **Distill:** exit 0. Sheets: MISC, carbon utilization, Transporters, Energy, Organic Nitrogen, carbon utilization (Woodcroft), **Salterns methane**, rRNA, tRNA. `product.html`, `product.tsv` and `genome_stats.tsv` were written.
- **`check_test.py` output:**
  - **Columns:** all 25 expected, including `ko_id`, `kegg_hit`, `pfam_hits`, `cazy_*`, `peptidase_*`, `bin_taxonomy`, `bin_completeness` and `bin_contamination`.
  - **Taxonomy/quality join:** `d__Archaea;…;g__Methanolobus;s__Methanolobus sp020055655`, completeness 98.37, contamination 2.61.
  - **KOs:** 1,093 of 2,279 genes (0.480). *mcrA/B/G*, *mtbB*, *mtmB*, *mttB*, *mtbA* (×2), *mtbC*, *mttC* and *pylB/C/D* are all present, plus *pylS* in the sheet. *mtmC* is absent.
  - **"Salterns methane" (65 rows), non-zero rows:** *phnD*, *phnP*, *betT*, *opuA*, *opuBD* ×2, *opuC*, *proV/W*, *proX* ×3, *mcrA/B/G*, *mtaA* ×3, *mtbA* ×2, *mtbB/C*, *mtmB*, *mttB/C* and *pylB/C/D/S*.
  - **No HTML** in the Pfam, CAZy or MEROPS columns, e.g. `GT81 NDP-Glc: glucosyl-3-phosphoglycerate synthase (EC 2.4.1.266)…` and `MER0030061 - family M38 …`.
  - **The only failure:** `Pfam fraction < 0.2` (52 genes), expected with 3,000 of 53,275 families. With the full Pfam a methanogen has well over 20 % of genes with a Pfam hit.
- **KOfam calls are reproducible:** applying DRAM's `kofam_hmmscan_formater` to smoke test 1's KOfam output gives the same 1,093 genes with a KO.

### The real run passed its own test and is annotating (status at 14:47 CST)

| stage | result |
|---|---|
| setup phase B (Pfam) | **done 14:07:28**, 1 h 50 min instead of the 3–6 h estimated. Peak RSS 254 GB, mostly the mapped 244 GB MSA file. `pfam.mmspro` 264 MB plus a 13 GB `.idx`. 53,275 Pfam descriptions in the live description DB. `logs/setup_check.log` shows all 17 configured paths exist; `uniref` and `kegg` are `None`, as intended |
| driver's one-MAG test | **`TEST_OK` 14:44:50, `TEST_CHECK OK`**. The same numbers as smoke test 2, plus the full Pfam: **1,789 of 2,279 genes with a Pfam hit (78 %)**. 36 min wall per MAG (full Pfam search 6 min), **peak RSS 24 GB** (Pfam MMSeqs2 search, which maps the 13 GB index) |
| 558-MAG annotation | pass 1 started 14:44:50. `xargs -P 8` (PID 3831773, child of the driver) runs 8 `DRAM.py annotate` processes, each at nice 10 and in its KOfam `hmmsearch --cpu 6` at 84–205 % CPU. Status: 1 done (the test MAG), 0 failed, no "database is locked", "disk I/O error" or Traceback in `logs/mags/*.log` |
| memory | 8 jobs × up to 24 GB peak ≈ 190 GB at worst. `free -g`: 852 GB free, 2.4 TB available |

**Expected finish:** at ~36 min for 2,279 genes, the average MAG (~3,000 genes) takes about 45 min, so 557 MAGs / 8 ≈ **50–60 h**. Then merge and distill take 1–3 h. Expect `state/DONE` on about **13–14 Nov 2026**, depending on load.

**Only process to keep alive:** driver **PID 4007531** (unchanged). Check it with `bash scripts/status.sh` as above.

**Small things, none worth relaunching for:**
- Each DRAM Python process has 128 LWPs. These are idle OpenBLAS/numpy pools created at import, not work threads. CPU per job stays around 1–2 cores.
- With `--cpu 6`, `hmmsearch` has 7 threads (6 workers and the master), but it uses about 1.2 cores (see the benchmark).

### Reviewer verdict

**Sound. Nothing on poplar was changed.**
- No script, input, database, environment or job of this track was edited, killed or relaunched.
- The scripts, the 558-MAG list, the FASTAs, the taxonomy and quality inputs, the KOfam release (shared with the KOfam track) and the threshold rule all check out.
- Every DRAM code path has now run on real output, and the run is past the gating test.

**Corrections to the text above:**
- the Pfam build took 1 h 50 min, not 3–6 h;
- annotation needs ~36–45 min per MAG at ~1.1–1.5 cores each;
- the SQLite-on-NFS watch item is low risk (8 concurrent readers tested);
- DRAM 1.5.0 can't run with Pfam unconfigured (rank-D FASTA header bug), which is irrelevant here because Pfam is configured.

**Review files on poplar:** `$W/review/`, 1.17 GB in total, of which 1.12 GB is `e2e/` (mini Pfam profile DB and description-DB copy). The scripts are `smoke.sh`, `e2e.sh`, `bench.sh`, `sqlite_conc.py` and `synth_distill.sh`, with their logs and outputs. Nothing in the pipeline reads them, so they can be deleted at any time: `rm -rf /scratch1/afreiburger/salterns_2026/dram/review`. Local copies of the scripts are in `review_scripts/` here (not `scripts/`, which stays identical to the copy on poplar).
