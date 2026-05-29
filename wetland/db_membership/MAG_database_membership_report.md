# Wetland MAG database membership: SPIRE & MGnify

**Question:** Are the MAGs derivable from the samples in `Sample_metadata_wetland_only.csv`
present in the [SPIRE](https://spire.embl.de) or [MGnify](https://www.ebi.ac.uk/metagenomics) databases?

**Date:** 2026-05-26
**Inputs:** `Sample_metadata_wetland_only.csv` (1,400 sample rows; 1,393 SRA runs + 7 non-INSDC; 774 BioProjects)
**Pipeline:** [`check_mag_databases.py`](../check_mag_databases.py) (cached/resumable)

---

## TL;DR

| | SPIRE | MGnify |
|---|---|---|
| **Runs with MAG(s)** | **612 / 1,400** | **0 / 1,400** |
| Unique biosamples covered | 578 | — |
| **MAGs linked (de-duplicated)** | **10,243** | 0 |
| Matched by | `SAMN` biosample (via ENA crosswalk) | per-run `/runs/{SRR}` — all 404; 0/774 BioProjects present |

- **SPIRE** holds ~10,243 MAGs derivable from 578 of these samples.
- **MGnify** contains **none** of these samples (verified two independent ways).
- **Direct SRR matching does not work for SPIRE here** — all hits required the run→biosample crosswalk.

---

## Why the CSV needs a crosswalk

The CSV is keyed by **SRA run** (`SRR…`) and **BioProject** (`PRJNA…`), and contains **no MAG identifiers**.
Neither database is queried by MAG name from this side — both are indexed by sample/run/study accession:

- **SPIRE** indexes MAGs by INSDC **biosample** (`derived_from_sample`); sequencing runs for one
  biosample are pooled before assembly. So a run must first be mapped to its biosample.
- **MGnify** normalizes run accessions, so `SRR…` can be queried directly.

A direct intersection of the 1,393 run accessions against SPIRE's 443 run-keyed entries returned
**0 matches** — confirming SPIRE stores these MAGs under biosample IDs, not runs.

---

## Method

1. **Index SPIRE.** Downloaded `spire_v1_genome_metadata.tsv.gz` (1,158,553 MAG rows across 73,675
   source samples); built a map `derived_from_sample → [genome_id, taxonomy, completeness, contamination]`.
2. **ENA crosswalk.** For each BioProject, queried the ENA Portal `filereport` API to map
   `run_accession → sample_accession` (biosample) + secondary + study accession.
   → 774 BioProjects resolved to 2,618 run→biosample mappings.
3. **SPIRE match.** For each run, looked up its biosample (then secondary, then raw run) in the SPIRE index.
4. **MGnify — BioProject level.** Queried `studies?accession={PRJNA}` (the filter matches the BioProject form).
5. **MGnify — run level.** Queried `/runs/{SRR}` (200 = present, 404 = absent), then `/runs/{SRR}/analyses`
   for present runs to record analysis types (assembly analyses are what yield MAGs).

All API endpoints were validated against known-present controls before the full run.

---

## Results

### SPIRE coverage by habitat

| Habitat | Runs with SPIRE MAGs | Total runs |
|---|---|---|
| Peatlands | 448 | 884 |
| Wetlands | 126 | 423 |
| Permafrost regions | 38 | 93 |
| **Total** | **612** | **1,400** |

- Unique biosamples with MAGs in SPIRE: **578**
- Unique SPIRE MAGs linked (de-duplicated): **10,243**
- All 612 matches were via the `SAMN` biosample key (0 via raw run, 0 via secondary).

### MGnify

- Per-run: **1,393 / 1,393 returned HTTP 404** (genuine "not found"), 0 present.
- Per-BioProject: **0 / 774** present.
- Conclusion: no MGnify MAGs derive from this wetland sample set.

### Non-INSDC rows

7 rows (`DZ1B`, `DZ2M`, `H03B`, `QM01A_25-30`, `Z2-45-50`, `LH01`, `LH02`; BioProject `"This study"`)
are newly-sequenced samples without INSDC accessions — unmatchable by accession in either database,
correctly reported as misses.

---

## Caveats

- A SPIRE match means SPIRE assembled MAG(s) from that **biosample**. Because runs are pooled per
  biosample, the MAGs are not attributable to an individual sequencing run.
- "In MGnify" here means the run/study is registered and analyzed by MGnify. MGnify also maintains
  biome-specific **genome catalogues** (`/api/v1/genome-catalogues`) assembled at the catalogue level;
  none of these samples are present in MGnify at the run or study level, so they contribute no MGnify MAGs.

---

## Artifacts

| File | Description |
|---|---|
| [`mag_database_membership.csv`](mag_database_membership.csv) | One row per run: biosample, `spire_mag_count`, `spire_mag_ids`, `spire_match_key`, `mgnify_study_present`, `mgnify_study`, `mgnify_run_present`, `mgnify_run_analysis_types` |
| [`../check_mag_databases.py`](../check_mag_databases.py) | Reproducible, cached/resumable pipeline |
| `cache/spire_v1_genome_metadata.tsv.gz` | SPIRE MAG table (92 MB) |
| `cache/ena_bioproject.json` | ENA run→biosample crosswalk cache |
| `cache/mgnify_studies.json`, `cache/mgnify_runs.json` | MGnify lookup caches |

### Data sources

- SPIRE genome metadata: `https://black.embl.de/~fullam/spire/metadata/spire_v1_genome_metadata.tsv.gz`
- SPIRE: Schmidt et al., *SPIRE: a Searchable, Planetary-scale mIcrobiome REsource*, NAR 2024 ([PMC10767986](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10767986/))
- ENA Portal API: `https://www.ebi.ac.uk/ena/portal/api/`
- MGnify API: `https://www.ebi.ac.uk/metagenomics/api/v1/`
