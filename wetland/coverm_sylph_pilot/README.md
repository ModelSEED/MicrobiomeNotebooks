# CoverM vs sylph pilot: MAG abundance for the wetland metagenomes (2026-11-13)

**Question.** Should per-sample MAG abundance (read mapping with CoverM, or k-mer profiling with sylph) be added to the
wetland MAG pipeline (KBase narratives 261700 reads → 272312 MAGs), and with which tool?

**Short answer.** Yes, add abundance, and use **sylph as the default** for the full sample × MAG matrix:
- **Same answer:** sylph gives the same relative abundances as CoverM (Spearman ρ = 0.99, Bray–Curtis 0.05).
- **Far cheaper:** about 1/65 of CoverM's CPU and 1/18 of its memory (about 1/14 when profiling against all 12,859 MAGs), with no BAM files.
- **No dereplication needed:** it can profile against all 12,859 MAGs at once.

Keep **CoverM for targeted work that needs alignments**: strain/SNV analysis (inStrain), coverage-based MAG QC, and
read-recruitment plots, on a chosen subset of samples.

## Setup
| | |
|---|---|
| Samples | SRR21120499 (10.7 Gbp) and SRR21120498 (11.4 Gbp): PRJNA869844, "Wetlands"; the KBase read libraries in ws 261700 (interleaved FASTQ.gz, md5 recorded) |
| Reference for the head-to-head | all 359 MAGs of PRJNA869844's 9 metagenomes, dereplicated at 95% ANI with dRep (skani; CheckM v1 ≥50% / ≤10%): **188 representatives** (327 passed the CheckM v1 filter) |
| CoverM | 0.7.0, `coverm genome`, minimap2-sr 2.30, defaults (min covered fraction 10%), 32 threads; second pass on the same BAM with `--min-covered-fraction 0` for counts |
| sylph | 1.0.0 (Bioconda), `sylph sketch` (reads as single-end), `sylph profile -u` (estimate unknown), 32 threads |
| Host | poplar (336 cores), heavily loaded: 1-minute load average 440–690 during the run (~560 on average, about 1.7× its 336 cores), so the wall times below are pessimistic |

Scripts: this folder (originals in `~/gtdbtk_kbase/coverm_pilot/`); run directory
`poplar:/scratch1/afreiburger/wetland_coverm_pilot`; outputs and per-MAG comparison tables in
`~/gtdbtk_kbase/coverm_pilot/results/`; `summary.json` and `compare_<sample>.tsv` are copied to `results/` here.

## Results
| Per sample | SRR21120499 | SRR21120498 |
|---|---|---|
| Reads explained: CoverM, 188 MAGs (after its default 10% covered-fraction filter; all reads mapped to the 188 MAGs in brackets) | 67.9% (68.5%) | 66.0% (66.7%) |
| Reads explained: sylph sequence abundance, 188 MAGs | 75.0% | 70.3% |
| Reads explained: sylph, **all 12,859 MAGs** (no dereplication) | 81.6% | 72.4% |
| MAGs detected: CoverM / sylph / both | 109 / 103 / 98 | 110 / 99 / 96 |
| Relative abundance where both detect: Spearman / Pearson | 0.989 / 0.999 | 0.988 / 1.000 |
| CoverM trimmed-mean depth vs sylph coverage: Spearman | 0.995 | 0.995 |
| Bray–Curtis, CoverM vs sylph profile | 0.051 | 0.047 |
| Detected by only one tool | CoverM-only 11 (each ≤0.22%, ≤4.8×); sylph-only 5 (each ≤0.015%, ANI 96–99%) | CoverM-only 14 (≤0.12%); sylph-only 3 (≤0.012%) |
| sylph detections at <99% ANI (a related strain, not the MAG itself) | 18 (0.7% of abundance) | 13 (0.6%) |
| sylph vs all MAGs: own-sample MAGs found / MAGs from other projects | 48 of 48 / 2 (0.01% of reads) | 49 of 49 / 4 (0.02%) |

| Cost per sample (~11 Gbp) | CoverM | sylph |
|---|---|---|
| Wall time (32 threads, loaded host) | 17–18 min | 1.3–1.5 min to sketch reads; profiling both samples in one call took 7 s (33 s vs all 12,859 MAGs, ~26 s of it loading the database) |
| CPU | 14–16 core-min per Gbp | 0.2 core-min per Gbp (**~65× less**) |
| Peak memory | 25 GB (with only a 188-MAG reference) | 1.4 GB to sketch reads; 0.5 GB to profile against the 188 MAGs (1.8 GB against all 12,859) |
| Reference preparation | dereplication required (4.5 min here); minimap2 index rebuilt every run | 12,859-MAG database built once in 2.7 min (1.5 GB RAM) |
| Intermediate files | BAM per sample | a small sketch per sample; reads can be deleted right after |

## What it means for the whole data set (702 metagenomes with MAGs, 12.2 Tbp; ~6.7 TB gzipped in KBase, estimated from the two pilot files)
- **sylph:** about 45 core-hours of sketching plus minutes of profiling. Downloading the reads is the bottleneck: about
  23 h at the ~80 MB/s measured from KBase to poplar. Each sample can be downloaded, sketched and deleted, so disk use
  stays at a few files. The output is the sample × MAG abundance matrix (taxonomic and sequence abundance, plus ANI) for
  all 12,859 MAGs.
- **CoverM:** about 3,000 core-hours even with a small reference. A useful all-MAG reference (the census's ≈6.1k
  species-level representatives, ~19 Gbp) would need roughly a 60 GB minimap2 index per job (rough estimate scaled from
  the pilot reference's index) and BAMs to manage. That is feasible on
  poplar but costs ~2 orders of magnitude more than sylph, for the same abundance answer.

## Value to the pipeline
1. **Prevalence and abundance of every MAG across all 702 metagenomes.** This is the basis for relating MAGs and guilds
   (e.g. methane-cycle MAGs) to habitat and biogeography. In this pilot, MAGs from other projects contributed only 0.01–0.02%
   of reads, so wetland MAGs are largely study-specific; that is a finding the full matrix would quantify.
2. **How much of each community the MAG collection captures** (here 66–82% of reads). This is the context any functional
   census needs.
3. **Representative selection (handoff Phase 5).** Rank cluster members by prevalence/abundance as well as quality, and
   flag clusters whose MAGs recruit reads at <99% ANI (strain heterogeneity).
4. **Not needed for it:** BAMs. Use CoverM only where alignments are the point (inStrain/SNVs, coverage-based QC), on
   selected samples.

## Caveats
- Only two samples, from one project. Both tools ran against MAGs from that project, and sylph also against all MAGs.
- sylph treated the interleaved pairs as single-end reads. That is acceptable for profiling, but `-1/-2` would be cleaner.
- The timings come from a host at ~1.7× oversubscription and are pessimistic for both tools. The CPU ratio is what to rely on.
- CoverM's defaults explain the disagreements. It applies no read-identity or alignment-length filter, so reads from
  related strains are counted, which gives its extra low-abundance detections. Its 10% minimum covered fraction gives
  zero coverage for the few divergent (95–99% ANI) genomes that sylph reports.
