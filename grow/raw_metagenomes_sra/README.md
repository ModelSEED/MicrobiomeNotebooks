# GROWdb raw metagenomes (the samples behind the 3,825 MAGs of KBase ws 145226)

`runs.tsv` lists 171 SRA runs from 133 samples: every metagenome the GROWdb MAGs were binned from.
- **NASQAN and Yukon:** 65 NASQAN samples and 2 Yukon samples, in BioProject PRJNA946291 (GROWdb, Colorado State).
- **WHONDRS 2019:** 66 surface-water samples, one JGI BioProject per sample (PRJNA7103xx–7104xx and PRJNA13166xx–13173xx).
  - SRA sample names don't contain "WHONDRS"; these runs were found through the BioProject titles.
  - The 11 WHONDRS runs that CSU re-deposited in PRJNA946291 are byte-identical duplicates (same spots and bases) and are left out.

The columns are sample → BioProject → BioSample → run, with spots and bases, ENA `fastq_ftp`/`fastq_md5`, the AWS ODP `.sra` URL and size, and the number of MAGs from each sample.
`GROW_collection_samples.tsv` is the KBase **GROW collection** sample table: 133 samples with IGSN, coordinates, collection date and program (ROMEO = NASQAN, WHONDRS, PARTNERS). The collection names WHONDRS samples `S19S_NNNN`.

## Local copy: `download_sra.py`
- **Where:** `/mnt/cold/Documents/GROW_raw_metagenomes/<sample>/<run>_{1,2}.fastq.gz`, downloaded from ENA. That is 342 files, 1.74 TB.
- **Checks:** each file must match ENA's size and MD5, and a `.md5` file is written next to it. The run can be restarted and resumes partial files. Progress goes to `download.log` there.
- **Speed:** ENA caps this machine at about 23 MB/s in total.
- **Start command:** `python download_sra.py runs.tsv /mnt/cold/Documents/GROW_raw_metagenomes 16`

## KBase copy: `kbase_import.py`
- **Narrative:** ws **275929** `freiburgermsu:narrative_1794152940209`, "GROWdb raw metagenomes: 133 MAG source samples (SRA)". Created 2026-11-08.
- **App:** `kb_uploadmethods/import_sra_as_reads_from_web` (beta), the same app and parameters as `salterns/reassembly_2026`. KBase pulls each run from the AWS ODP mirror itself, so nothing is uploaded from this machine.
- **Objects:** one interleaved PairedEndLibrary per run, `<sample>.<SRR>.reads`. Multi-run samples are not merged.
- **Acceptance:** a library is accepted only if `read_count == 2 × spots` and `total_bases == bases` from SRA.
- **Pilot:** the largest run (SRR35181274, 49 Gbp, 14.8 GB `.sra`) runs alone first. The other runs fan out once it validates.
- **Job cap:** at most 10 live njs jobs; the salterns watcher caps its own njs jobs at 15.
- **Expected size:** about 2.7 TB in KBase (0.72 bytes/base × 3.72 Tbp, measured on the salterns imports). For comparison, the wetland narrative 261700 already holds 7.4 TB of reads.

```bash
~/Documents/py_venv/bin/python kbase_import.py status
# the watcher (started 2026-11-08) logs to kbase_import.log and exits when every run is done or failed; restart it with:
setsid nohup ~/Documents/py_venv/bin/python -u kbase_import.py watch --interval 900 >> kbase_import.log 2>&1 < /dev/null &
```
Only one advance or watch may run at a time (`kbase_import_ledger.json.lock`).

**Last step, attaching to the GROW organization (automatic):**
- **When:** once all 171 libraries are `done`, the watcher sets the narrative description and calls `POST groups/group/grow/resource/workspace/275929`.
- **Approval:** you are a Member of the KBase Org `grow`, not an admin, so this files a request that an Org admin must approve. The admins are owner mikayla_borton, plus growdb, janakakbase, allenbh, egdow, elishawc, sjurgensen and zcrockett. The result (request id) is recorded under `org_add` in the ledger.
- **If runs fail:** if any run ends `failed`, the request is NOT filed and the watcher logs why.
- **Not the Collections service:** the KBase *Collections* service (`/collections/GROW`) is different. It is curated by KBase staff (service admins only) and has no data product for reads.
- **Private narrative:** the Org also holds a private narrative you cannot read, "GROWdb US River Systems - Transcriptomics" (ws 109710).

