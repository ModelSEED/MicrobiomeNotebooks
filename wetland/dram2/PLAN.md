# Wetland DRAM2 run on poplar (afreiburger) — plan and record

Goal: DRAM2 (rename → call → annotate → distill → product) for the 1,184 wetland MAGs in
`dram2_hq_multisample.tsv` (one representative per species cluster; CheckM2 >90 % complete,
<5 % contaminated; cluster seen in ≥2 samples), using chenry's native DRAM2 install read-only.

## What is reused from chenry's install (read-only, nothing written there)

| Item | Path | How it is used |
|---|---|---|
| Pipeline | `/scratch1/fliu/hub_scratch/chenry/DRAM2/repo/main.nf` | DRAM2 v2.0.0-beta17, revision 7754dfaeaf |
| Databases | `.../DRAM2/databases` | `launch/<run>/databases` is a symlink to it (config resolves `${launchDir}/databases/<db>`) |
| JDK 17 | `.../DRAM2/env/env_nf` | `NXF_JAVA_HOME` |
| annotate env | `.../DRAM2/nxf-scratch/conda_cache/env-fd46c7768a790ad1-f25d581e073c803fd2df039b7aadc6c8` | hmmer 3.3.2, mmseqs2 13.45111 (built at the /scratch1 path, used by chenry's KBDL runs) |
| multiqc env | `.../DRAM2/nxf-scratch/conda_cache/env-152dac747f07927d-eec4e05336245ea50db6412788f8c4da` | multiqc 1.25.1 |

## What is our own (under `/scratch1/afreiburger/wetland_dram2`)

- `nxf_home/` — own NXF_HOME. chenry's `nextflow-24.10.5-one.jar` copies are mode 600 (unreadable
  to afreiburger), so the identical official release jar was fetched from
  `https://www.nextflow.io/releases/v24.10.5/nextflow-24.10.5-one.jar` (34,364,570 bytes, same size
  as chenry's; sha256 bac779dd8e0eed9f78e2475e4327ae0fc095594c6ad379375133922192a50974). Plugins
  nf-schema-2.0.1 and nf-validation-1.1.4 copied from chenry's `.nextflow_home/plugins`.
- `bin/nextflow-classic`, `bin/micromamba` — copies of chenry's launcher script and micromamba 2.8.1
  (micromamba only activates envs; no env is ever created by Nextflow).
- `conda_cache/dram2_beta17_{call,collect_rna,rename,distill,product}` — built with the user's conda
  (`/scratch/afreiburger/metag-hybrid/env/miniforge/condabin/conda`, 26.3.2) from the exact beta17
  `environment.yml` files (`envs_yml/`; explicit package lists in `envs_yml/*.explicit.txt`).
  Reason: chenry's call / collect_rna / rename envs (`nxf_conda/env-8393b9af…`, `env-21ab46b9…`,
  `env-7a8564ac…`) were built on popjup under `/storage/chenry/DRAM2`, which does not exist on
  poplar (`quast.py` shebang and `tRNAscan-SE.conf` `bin_dir` point there), and distill / product
  envs were never built.
- `launch/wetland_dram2.config` — the `-c` config: local executor capped at 56 CPU slots / 400 GB
  (at 64 slots the measured real use was 59-67 cores: 30 hmmsearch = 50 cores + ~9 cores of
  Nextflow's per-task trace-polling bash; chunk_01 was restarted with -resume at 56 slots);
  label caps small 2 CPU/8 GB, medium 4/48, big 8/64, huge 2/16 (MMseqs2; at 4 CPUs the 4-slot
  searches starved behind the 2-slot hmmsearch stream: 1 of 948 had run), longer wall limits;
  every module's conda env pinned to an explicit prefix; QUAST confined to 4 CPUs with `taskset`
  (QUAST would otherwise start cpu_count()//4 = 84 workers); per-MAG MMseqs2 query indexes and
  the raw HMM_SEARCH/PARSE_HMM files not published (every database's raw file has the same name,
  `<mag>_hmmsearch.out`, so the published copies overwrote each other) EXCEPT KOfam's raw
  `--domtblout`, published to `runs/<chunk>/out/HMM_search_kofam_raw/` (~45 MB/MAG, ~50 GB total;
  the only place with per-domain scores, needed to apply KOfam's thresholds properly — delete if
  not wanted); BLAS thread pools set to 1.
- `bin/janitor.sh` — runs beside each chunk; for tasks that finished OK it truncates hmmsearch's
  full-alignment stdout (.command.out/.command.log, up to ~400 MB each for KOfam) and deletes the
  unused ~0.9 GB per-MAG query index (`<mag>.mmsdb.idx`) once all 4 MMseqs searches of that MAG
  are done. 3-MAG test work dir: 4.67 GB → 0.3 GB. If /scratch free space drops below 120 GB it
  writes `runs/<chunk>/LOWDISK` and stops the chunk (resumable); the driver then stops.
- The user's conda (`/scratch/afreiburger/metag-hybrid/env/miniforge/condabin`) is on PATH because the
  pipeline's startup check runs `conda config --show channels` (chenry had env_nf/bin on PATH).
- `launch/nfcore_offline/nfcore_custom.config` — empty; passed as `--custom_config_base` so the
  launch does not fetch nf-core institutional configs from GitHub.
- `inputs/name_map.tsv` and `inputs/<chunk>/<name>.fna` symlinks to
  `/scratch1/afreiburger/wetland_coverm_pilot/all_mags/<mag>.fna`. DRAM2 needs dot-free names
  (`combine_annotations.py` recovers each gene's `input_fasta` as the common prefix of the .faa
  name and the renamed gene id), so `bin.014.fasta__SRR11577120` → `bin-014-fasta__SRR11577120`.
- Work dirs on local xfs: `/scratch/afreiburger/wetland_dram2_work/<run>`.

## Pipeline options

`--rename --call --annotate --input_fasta inputs/<chunk> --fasta_fmt '*.fna' --threads 2`
`--use_kofam --use_dbcan --use_merops --use_vog` (chenry's set) plus, if the 3-MAG test passes,
`--use_camper --use_canthyd --use_methyl --use_sulfur --use_fegenie`;
`--distill_topic default,camper --distill_ecosystem eng_sys,ag --product`.
Prodigal mode meta, table 11, `min_contig_len` 2500 (DRAM default). A 150-MAG sample of the
list has no contigs <2,500 bp, so this filter drops nothing.
Not used: pfam (on disk, but beta17 hard-disables it: `subworkflows/local/utils_nfcore_dram_pipeline/main.nf:89-90`
errors out with "PFAM database is currently disabled ... due to a bug"); kegg and uniref (not on disk); viral (`use_viral` is an undeclared TODO param in
beta17); tcdb, card, metals, antismash, dram_db (on disk but not wired into beta17 annotate).

## Steps

1. 3-MAG test (`runs/test3`): bin.072.fasta__SRR9700688 (0.55 Mb), bin.049.fasta__SRR8556254
   (3.7 Mb), bin.020.fasta__SRR11577034 (8.3 Mb). Check annotations, distill, product; time it.
   DONE (see Results).
2. Full run: 5 chunks of 237/237/237/237/236 MAGs (sorted by size, dealt round-robin so every
   chunk has the same size mix), one Nextflow launch per chunk, run sequentially by
   `bin/run_all_chunks.sh` (56 CPU slots). Each chunk: `runs/<chunk>/DONE` after verification
   (Nextflow exit 0; every input MAG present in `RAW/raw-annotations.tsv`; per MAG annotation rows ==
   genes in its Prodigal .faa and >= 1 KOfam hit; no "Error loading" in RAW/combine_annotations.log;
   distillate.xlsx and product.tsv present; warns if raw KOfam files != MAG count). Restart = rerun the driver (DONE chunks skipped,
   others resumed with `-resume`).
3. Merge: `bin/merge_and_distill.sh` → `merged/raw-annotations.all.tsv` (with a `mag` column) and
   one more DRAM2 launch (`runs/all_mags`) that runs distill + product over all MAGs.

## Results of the 3-MAG test (2026-11-15, load ~450-550 on 336 threads)

- Wall 51 min (13:20 → 14:12), Nextflow exit 0; 626 / 3,429 / 7,409 genes; all 3 MAGs in
  `RAW/raw-annotations.tsv` (11,464 rows, 47 columns); `DISTILL/distillate.xlsx`,
  `DISTILL/summarized_genomes.tsv`, `PRODUCT/product.{tsv,html}`, `tRNA_rRNA/`, `QUAST_v5.2.0/` present.
- Gene ids are `<mag>_<contig>_<n>` (e.g. `bin-072-fasta__SRR9700688_k141_1025319_1`).
- Hits per MAG (626 / 3,429 / 7,409 genes): kofam 441 / 2,784 / 4,963; vogdb 153 / 718 / 1,362;
  dbcan 27 / 65 / 348; merops 21 / 112 / 452; cant_hyd 19 / 285 / 393; camper 1 / 70 / 106;
  fegenie 7 / 19 / 57; sulfur 0 / 13 / 24; methyl 0 / 38 / 22.
- Cost per MAG (slot-seconds = wall × requested CPUs, what the 64-CPU cap packs):
  2,318 / 4,891 / 8,882  →  slot_s ≈ 1,700 + 0.97 × genes. KOfam hmmsearch is ~65 % of it
  (6 min / 31 min / 46 min wall at ~1-1.9 cores); VOG ~20 %.
- Projection for the 1,184 MAGs (≈3.85 M genes at 0.913 genes/kb): ≈1,520-1,600 slot-h (MMseqs
  label now 4 CPUs) → ≈5.5 h per 237-MAG chunk at 56 slots + ~1 h tail per chunk → ≈32 h for all
  5 chunks (range ~1.2-1.8 days with load), + ~0.5 h for the final merge/distill.
- Run log: chunk_01 started 14:19 at 64 slots; stopped 14:50 (SIGTERM) and resumed 15:41 at 56 slots
  (1,034 tasks cached); stopped again 16:01 and resumed 16:02 (1,483 cached) with MMseqs at 2 CPUs,
  raw KOfam publishing and the per-MAG verification. Exact CPU of the Nextflow tree at 56 slots
  (utime+stime+cutime+cstime over the JVM and descendants, 312 s window): 42.0 cores. Validated only
  at host load ~550-600; on a quiet host hmmsearch can use >2 cores per 2-slot task. Measured at 56 slots: 42-44 cores over 20 s windows, 54-68 over 5 s windows
  (whole Nextflow process tree; short-lived wrapper children not counted). To change the slot count,
  edit `executor.cpus` in launch/wetland_dram2.config; it applies from the next chunk launch.
- The all-MAG distill/product path (`bin/merge_and_distill.sh`, `--annotations` + `--rrnas/--trnas`,
  no `--call`) was tested on the test3 output (`runs/test3_alldistill`): product.tsv and
  summarized_genomes.tsv are identical to the chunk-level ones. It needs `--validate_params false`
  (beta17's schema types `--rrnas/--trnas` as files but the code globs a directory).

## Caveat: KOfam calls are NOT threshold-filtered in DRAM2 beta17

`bin/kofam_hmm_formatter.py` keeps the best KOfam hit per gene with E ≤ 1e-5 and never applies
KOfam's per-KO adaptive thresholds (`kofam_ko_list.tsv` threshold/score_type), unlike KofamScan
and DRAM1. In the test 67-81 % of genes got a KO; only ~52 % of those hits reach the KO threshold
(1,450 / 2,784 for bin.049; 257 / 441; 2,458 / 4,963). `kofam_score_rank` is the hit's
full-sequence score. `bin/merge_dram2.py` therefore adds `kofam_ko_threshold`,
`kofam_ko_score_type`, `kofam_passes_ko_threshold` (score ≥ threshold; for "domain"-type KOs an
upper-bound approximation since the domain score is not kept). DRAM2's own distill/product use the
unfiltered `kofam_id`, so module completeness from DRAM2 is optimistic.

## Exact command per chunk (from `runs/<chunk>/command.sh`)

```
cd /scratch1/afreiburger/wetland_dram2/launch/<chunk>     # databases -> chenry/DRAM2/databases
source /scratch1/afreiburger/wetland_dram2/bin/dram2_settings.sh
nextflow-classic -log logs/<chunk>/nextflow.log run /scratch1/fliu/hub_scratch/chenry/DRAM2/repo/main.nf \
  -profile conda -c /scratch1/afreiburger/wetland_dram2/launch/wetland_dram2.config \
  -work-dir /scratch/afreiburger/wetland_dram2_work/<chunk> -resume -ansi-log false \
  --rename --call --annotate --input_fasta inputs/<chunk> --fasta_fmt '*.fna' --outdir runs/<chunk>/out \
  --threads 2 --use_kofam --use_dbcan --use_merops --use_vog --use_camper --use_canthyd --use_methyl \
  --use_sulfur --use_fegenie --distill_topic default,camper --distill_ecosystem eng_sys,ag --product \
  --custom_config_base /scratch1/afreiburger/wetland_dram2/launch/nfcore_offline
```

## Operating it

- Progress: `bin/status.sh` (per-chunk DONE/FAILED/running with task counts, CPU in use);
  `tail logs/driver.log`; per chunk `logs/<chunk>/stdout.log`, `logs/<chunk>/nextflow.log`,
  `runs/<chunk>/janitor.log`, trace `runs/<chunk>/out/pipeline_info/execution_trace_*.txt`.
- Restart after any stop: `cd /scratch1/afreiburger/wetland_dram2 && setsid nohup bin/run_all_chunks.sh > logs/driver.log 2>&1 < /dev/null &`
  (DONE chunks skipped; others `-resume`; single-instance lock `logs/driver.lock`).
- After `runs/ALL_CHUNKS_DONE` appears:
  `setsid nohup bin/merge_and_distill.sh > logs/merge_and_distill.log 2>&1 < /dev/null &`
  → `merged/raw-annotations.all.tsv`, `runs/all_mags/out/DISTILL/{distillate.xlsx,summarized_genomes.tsv}`,
  `runs/all_mags/out/PRODUCT/product.{tsv,html}`, `runs/all_mags/DONE`.
- Stop (graceful, resumable): `bin/stop.sh`.

## Deviations from the original brief (all deliberate)

- Nextflow 24.10.5 jar downloaded from nextflow.io (chenry's copies are mode 600); same size, sha256 above.
- Five module envs built in our own tree with the user's conda (call, collect_rna, rename: chenry's
  copies point at /storage/chenry/DRAM2, which does not exist on poplar; distill, product: never built).
- `--rename` used and sample names changed from dots to dashes (DRAM2 needs both to keep each gene's
  MAG identity in combine_annotations.py).
- Pfam dropped (hard-disabled in beta17).
- Executor 56 slots (not 64) and MMseqs label 2 CPUs / 16 GB (not 16 / 120) — see above.
- `bin/janitor.sh` truncates hmmsearch stdout logs and deletes unused MMseqs query indexes in the
  work dir of finished tasks; the driver deletes a chunk's whole work dir after its DONE marker.
- chunk_01 was stopped and resumed twice (≈1 h of its in-flight hmmsearch work redone).
