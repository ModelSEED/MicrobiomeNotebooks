# Are the 2026 metaSPAdes assemblies more fragmented than the assemblies behind the original MAGs?

Run 2026-11-11. KBase was only read; nothing was saved to KBase.

## Short answer

**Yes, clearly, compared with the assembly the original MAGs came from.** Compared with the public JGI assemblies, they hold more sequence, in somewhat shorter contigs.

Each sample has three assemblies (`provenance/PROVENANCE.md`):

| | Assembly | Public? | What this directory measured |
|---|---|---|---|
| (ii) | JGI IMG-team "SPAdes version 2" reassembly, 2018 (Ga0209xxx; `provenance/assembly_vs_bins.tsv`) | yes, 38/39 | exact lengths of every contig ≥ 500 bp, from the CoverM BAM headers (33 samples) |
| (iii) | Tringe-lab Tadpole + SPAdes reassembly, about 2018. **All 450 original MAGs come from this one.** | no | contigs ≥ 2.5 kb, reconstructed from the MAG contig names (39 samples, see Methods) |
| new | 2026 KBase metaSPAdes 3.15.3, scaffolds ≥ 1 kb | ours | exact lengths (39 samples) |

The task's premise needs one correction. `assembly_vs_bins.tsv` describes (ii), which is not the assembly the MAGs were binned from. Measured on contigs ≥ 2.5 kb (what the binners see), the three assemblies rank like this:

| Contigs ≥ 2.5 kb | (ii) JGI 2018 | new 2026 | (iii) MAG source |
|---|---|---|---|
| Total Mbp, 33 samples with (ii) | 2,069 | 3,208 | ≈ 5,447 |
| Total Mbp, all 39 samples | – | 3,691 | ≈ 6,300 |
| Median per-sample ratio, new / other | 1.62 (new > (ii) in 33/33) | 1 | **0.57** (new < (iii) in 39/39; range 0.31–0.74) |
| Median N50 of these contigs, 33 samples | 6.1 kb | 4.8 kb | 8.4 kb |
| Mbp in contigs ≥ 10 kb, 33 samples | 706 | 722 | ≈ 2,390 |
| Mbp in contigs ≥ 20 kb, 33 samples | 372 | 298 | ≈ 1,513 |

- **Against the MAG-source assembly (iii):** the 2026 assemblies hold 57 % as much sequence in contigs ≥ 2.5 kb (median per sample) and about 30 % as much in contigs ≥ 10 kb. Their N50 on contigs ≥ 2.5 kb is 0.64× (median).
  - The shortfall is in the long contigs, and it explains why the 2026 MAGs have half the N50 of the original MAGs.
  - It is worst in the SF2 cores (median 0.42, range 0.31–0.60) and smallest in R2_restored (0.65–0.74).
- **Against JGI (ii): a trade-off.**
  - The 2026 assemblies hold more sequence at every threshold up to 10 kb: 1.9× the Mbp in contigs ≥ 1 kb and 1.6× in contigs ≥ 2.5 kb.
  - (ii) has the longer contigs: N50 of contigs ≥ 2.5 kb 6.1 vs 4.8 kb, and 372 vs 298 Mbp in contigs ≥ 20 kb.
  - Neither comes close to (iii) (8.4 kb; ≈ 1,513 Mbp ≥ 20 kb). Switching to the public JGI assemblies would trade completeness for somewhat longer contigs and would not recover the original contiguity.
  - (iii) itself holds 2.7× (median) the ≥ 2.5 kb Mbp of (ii). This matches PROVENANCE's "2–8× more contigs ≥ 2.5 kb".

## Why: the most likely cause is read preparation, not the assembler version

**What is the same in both runs (measured).**
- **Final k:** in SPAdes output, contigs that are adjacent in the graph overlap by exactly the final k. Counting shared contig-end k-mers gives a final **k = 127 for (iii)** in all 7 samples tested, far above background (e.g. R1_A_D1: 1,251 shared ends at k = 127 against about 350 at any other k). The 2026 assemblies, as a positive control, also give 127 (`final_k_overlaps.tsv`).
- **Input reads:** JGI holds one raw FASTQ per library, the same file as in SRA (PROVENANCE). Zhou et al.'s "~20 Gb per sample, 700 Gb total" is not borne out by SRA (344 Gbp), so I assume the same reads. This is not proven.

**What differs (measured).** I paired original and 2026 MAGs of the same species (skani ANI ≥ 95 %) from the same sample's own assembly: 190 pairs in 35 samples.
- **k-mer coverage:** SPAdes stores the k = 127 k-mer coverage in every contig name. In (iii) the same genome has a median **2.9× higher** k-mer coverage than in 2026 (IQR 2.5–3.1; per-sample medians 1.5–3.3; `kmer_cov_pairs.tsv`). With the same reads and the same k, only the reads fed to SPAdes can do this. At k = 127 a 151-bp read yields 25 k-mers (about 22 after the 2026 trimming), whereas a merged 245-bp fragment yields 119.
- **The gap grows with coverage** (`pair_contiguity.tsv`; N50 is over each MAG's contigs ≥ 2.5 kb):

| Quartile of 2026 k-mer coverage | 2026 coverage | (iii) coverage | Ratio | (iii) N50 | 2026 N50 | 2026 / (iii) N50 |
|---|---|---|---|---|---|---|
| Q1 | 0.92 | 2.04 | 2.27 | 8.3 kb | 4.6 kb | 0.56 |
| Q2 | 1.28 | 3.57 | 2.75 | 23.3 kb | 7.9 kb | 0.34 |
| Q3 | 1.63 | 5.04 | 3.06 | 41.7 kb | 9.4 kb | 0.23 |
| Q4 | 3.01 | 9.59 | 3.21 | 55.2 kb | 8.1 kb | 0.18 |

- **Reading of the table:**
  - In (iii), contig N50 keeps rising with a genome's coverage, to 55 kb. In 2026 it plateaus at about 8–9 kb however well covered the genome is.
  - So the 2026 fragmentation is not simple coverage starvation. Well-covered genomes break at repeats or strain variants that the original assembly got past.
- **The plateau is in the assembly, not the binning.** A binning-free check (`contiguity_by_coverage.py`, pooled over 39 samples) puts every 2026 contig ≥ 2.5 kb into a class of its own k-mer coverage; (iii) MAG contigs are classed at their coverage divided by 2.9.

| k-mer coverage (2026 scale) | 2026, all assembly contigs: N50 / % bp in ≥ 20 kb | 2026 MAG contigs | (iii) MAG contigs |
|---|---|---|---|
| < 1 | 3.5 kb / 0.1 % | 4.1 kb / 0.2 % | 7.7 kb / 12 % |
| 1–1.5 | 5.7 kb / 6 % | 6.9 kb / 8 % | 19.6 kb / 49 % |
| 1.5–2.5 | 7.7 kb / 19 % | 9.2 kb / 22 % | 35.6 kb / 70 % |
| 2.5–5 | 6.9 kb / 20 % | 7.7 kb / 23 % | 44.9 kb / 71 % |
| ≥ 5 | 7.1 kb / 20 % | 8.0 kb / 22 % | 27.4 kb / 60 % |

  The whole 2026 assembly stops at about 7–8 kb N50 above coverage 1.5, so the binners did not leave long contigs behind.

**What the original read preparation was (documented).**
- Zhou et al. 2022 (ISME J, PMC8692437), Methods: "2 × 150 paired-end ... with ~400 bp inserts with reads overlapped and merged", and "Each metagenome was re-assembled (de novo) using Tadpole error correction in BBTools and the SPAdes assembler to generate longer contigs".
- **Merge rate on these reads:** BBMerge 39.13 on 300,000 R1_A_D1 pairs merges 38.9 % with defaults (mean insert 240 bp) and ≥ 45.4 % with Tadpole extension (`rem extend2`; `merge_rate_R1_A_D1.txt`). That second figure is a lower bound, because extension needs the k-mer depth of the full read set.
  - The BBMerge paper reports 77.5 % correctly merged real-world pairs in its extension mode, and a near doubling of SPAdes NA50 from merged reads (Bushnell et al. 2017, PLoS ONE 12:e0185056).
  - Two readings fit the coverage ratio rising from 2.3 to 3.2 across the quartiles:
    - Tadpole extension merges more pairs in well-covered genomes.
    - 2026 metaSPAdes splits abundant, strain-rich genomes into partial-coverage paths that (iii) collapsed. Tadpole correction of minor variants, or a non-meta SPAdes mode (see below), would both collapse them.

**Pipeline differences, all three assemblies**

| Step | (ii) JGI "SPAdes v2", 2018 | (iii) Tringe lab, about 2018 (MAG source) | 2026 KBase |
|---|---|---|---|
| Reads | raw/unfiltered (README) | BBTools-filtered (Zhou: "filtering and trimming"), **pairs merged** | SRA raw → Trimmomatic 0.39 TruSeq3-PE-2, SLIDINGWINDOW 4:15, MINLEN 50 (kept 65–96 % of bases, median 92 %) |
| Error correction | bfc r181, k = 21 | **Tadpole** (BBTools) | BayesHammer (SPAdes; changed 0.2 % of bases in R1_A_D1) |
| Assembler | SPAdes 3.10.0-dev `--meta --only-assembler -k 21,33,55,77,99,127` | "the SPAdes assembler" (Zhou): **version and mode (`--meta` or default multicell) not reported**; final k = 127 (measured here) | metaSPAdes 3.15.3 `--meta -k 33,55,77,99,127` |
| Output | `assembled.fna` (scaffolds; IMG) | NODE contigs or scaffolds (unknown) | `scaffolds.fasta` filtered ≥ 1,000 bp (5 N per 100 kb) |
| Insert size seen by SPAdes | – | "~400 bp" nominal (Zhou) | 300 ± 67 bp (282–309 across samples), 4–45 % of pairs on long edges |

- **Sources for (ii):** the README text is inferred from a sibling README of the same 2018-03-21 IMG batch ([portal.nersc.gov/cfs/m342/BioIT/1006528/README.txt](https://portal.nersc.gov/cfs/m342/BioIT/1006528/README.txt): "raw/unfiltered reads were read corrected using bfc (r181)", "SPAdes version: 3.10.0-dev", "--meta --only-assembler -k 21,33,55,77,99,127"), as recorded in PROVENANCE.
- **The later JGI workflow does not describe (ii).** That workflow is rqcfilter → bbcms → metaSPAdes 3.13.0 assembly-only, k = 33–127 (Clum et al. 2021, *mSystems* 6:e00804-20, [PMC8269246](https://pmc.ncbi.nlm.nih.gov/articles/PMC8269246/)), and it postdates these 2018 assemblies. It is not evidence for how (ii) was made, and nothing local says (ii) used bbcms.
- **SPAdes version is unlikely to be the cause.** The core metaSPAdes algorithm changed little between 3.10 and 3.15 (`provenance/KBASE_VS_UPSTREAM.md`), and the (ii) assemblies are the least contiguous of the three.
  - (ii) shares the 2026 run's k-mer list (plus k = 21) and, like it, used unmerged pairs. It differs by bfc correction and by running on raw, untrimmed reads (no quality trimming at all).
  - So its contiguity says little about whether quality trimming hurt the 2026 run.
- **Quality trimming is a minor factor at most.** The ratio to (iii) correlates only weakly with the % of bases Trimmomatic kept (Spearman ρ 0.36), and the two samples that lost a third of their bases (SF2_B_H2O, SF2_C_D2) sit at 0.52 and 0.60, close to the median.

**Most likely causes, in order:**
1. Unmerged pairs.
2. BayesHammer instead of Tadpole correction (Tadpole also collapses low-frequency variants, which simplifies the graph of abundant, strain-rich genomes).
3. The SPAdes mode. Zhou says only "SPAdes". If (iii) ran default multicell SPAdes rather than `--meta`, it would merge strain variants of abundant genomes into longer consensus contigs. These are chimeric at the strain level, but contiguous.
4. Quality trimming shortening reads at k = 127.

None of these is proven by a controlled rerun. The pilot below is that test.

## What would fix it, and what it costs

| Option | Availability | Expected effect | Cost |
|---|---|---|---|
| A. Bin the public JGI assemblies (ii) | 38/39 on the JGI portal (`jgi_portal.tsv` `sp2_assembled_fna`; none for R1_C_H2O; login, possibly tape restore); 38 ≥ 500 bp copies listed under `/scratch/chenry/...` in `reads/metagenome_assembly.txt` (presence not verified) | **A trade-off, not a fix.** In 33/33 samples (ii) has less sequence in contigs ≥ 2.5 kb (median 0.62× of 2026), so fewer and less complete bins. Its contigs are somewhat longer (N50 6.1 vs 4.8 kb), but still far from (iii) (8.4 kb). | Not worth doing. The 2015 MEGAHIT 1.0.3 assemblies (i) (38/39) were not measured here; whether JGI merged reads for them is unknown. |
| B. Get assembly (iii) from the Tringe lab (Zhou / Tringe; the "JGI Salterns Modeling" Drive) | Not public: not in IMG, NCBI or KBase | Best: the original contiguity, re-binned with multi-sample coverage | An email. Re-binning in KBase and locally takes about 3–4 days, using the existing pipeline from the binner stage on. |
| C. Re-assemble with the original style of read preparation, outside KBase. Steps: BBDuk → `bbmerge-auto.sh rem extend2 ecct` → Tadpole correction → metaSPAdes, given the merged and unmerged reads, k = 127 included | Reads in SRA. No BBMerge or Tadpole app was found in the KBase catalog review (`provenance/KBASE_VS_UPSTREAM.md`; BBTools there sits only inside RQCFilter and jgi_mg_assembly). SPAdes and BBTools need a local micromamba env (BBTools is a tarball, works with the system Java 21). | Expected to recover much of the gap. Unproven, so pilot first. | **Pilot (R1_A_D1 + SF2_C_D1): about 1 day.** Full run: per sample about 6–12 h at ≤ 24 threads (estimate: the KBase run took about 3 h of BayesHammer plus 2.5 h of assembly at 64 threads); peak RAM 62–154 GB (median 101) in the 2026 KBase logs. On this shared machine that allows 1–2 at a time, about 1–2 weeks for 39 samples, plus about 170 GB of SRA downloads and about 100 GB of scratch per run. Then re-binning, CheckM and GTDB-Tk (about 3–4 days). |

**Pilot decision rule.**
- **Targets:** R1_A_D1's (iii) values here are 268.9 Mbp in contigs ≥ 2.5 kb and 131.4 Mbp in contigs ≥ 10 kb; the 2026 run has 162.8 and 32.6.
- **Arm 1 (read preparation):** BBTools preparation, then metaSPAdes 3.15.3 `--meta`. This uses the same SPAdes version and mode as KBase, so it isolates read preparation.
- **Arm 2 (SPAdes mode), only if arm 1 falls short:** the same reads in default multicell SPAdes mode. A strain-collapsed result may need GUNC or CheckM2 scrutiny.
- **Go / stop:** if an arm reaches about ≥ 85 % of the (iii) values, roll it out. If neither does, the remaining gap is in something unrecorded (SPAdes version or settings), and B is the only route.

## Methods and checks

- **2026 lengths** (`fetch_2026_lengths.py`):
  - **Source:** Workspace.get_objects2 with `included` contig lengths, for the ledger's 43 assembly refs.
  - **Check:** contig counts and bp equal the Assembly metadata. In all 39 metaSPAdes assemblies the NODE index is exactly the length rank (0 violations; max NODE = number of contigs).
  - **Bases < 1 kb:** contigs below 1 kb were never saved (`min_contig_length` 1000). So `n_contigs_ge500`, `asm_Mbp_ge500`, `median_len`, `N50_ge500`, `pct_bp_ge1000` and `pct_bp_ge2500` (denominators ≥ 500 bp) are left blank for "new" in `assembly_vs_jgi.tsv`. `new_pct_of_ge1000_in_ge2500` is the comparable shape measure.
- **(ii) lengths** (`jgi_v2_lengths.py`): BAM @SQ headers. Contig counts and Mbp ≥ 500 match `assembly_vs_bins.tsv` exactly, e.g. R1_A_D1 has 209,409 contigs and 294.1 Mbp. Six samples have no BAM: R1_C_D1, R1_C_H2O, R2_A_D2, R2_A_H2O, R2_C_D2, SF2_A_D1.
- **(iii) reconstruction** (`reconstruct_mag_source.py`):
  - **Principle:** SPAdes numbers NODE_1, NODE_2, … by decreasing length, so every MAG contig `NODE_n_length_L` is an exact point (rank n, length L) on the assembly's rank–length curve. All 450 MAGs give 459–9,444 points per sample, with 0 rank-order violations. Every contig between two sampled ranks has a length between theirs, which gives hard bounds; the estimate interpolates linearly.
  - **Validation on 2026,** where the truth is known, rebuilt from the 2026 MAG contigs subsampled to the same number of points as the originals: estimate / true Mbp ≥ 2.5 kb has a median of 1.000 (range 0.963–1.003); lower bound / true ranges 0.954–0.999; estimated N50 / true has a median of 1.000 (0.960–1.001) (`reconstruction_validation_2026.tsv`). The head (contigs longer than the first MAG contig) is counted at that contig's length, so (iii) values lean low.
- **Run statistics** (`spades_run_stats.py`): Trimmomatic retention from the reads objects; SPAdes version, command, insert size and % of pairs on long edges from the EE2 job logs (gzipped in `spades_logs/`).
- **Coverage and k checks:** `kmer_coverage_check.py` (final k from contig-end overlaps; same-genome k-mer coverage), `pair_contiguity.py`, and `contiguity_by_coverage.py` (binning-free). The overlap test was uninformative for 2026 R2A_B_D1 (too few adjacent MAG contigs); it gave k = 127 in the other 13 runs. `merge_rate_check.sh` holds the BBMerge commands; BBTools ran from a tarball in the session scratchpad.
- **Co-assemblies (context):** the four 2026 MEGAHIT sediment co-assemblies hold 0.71–1.31× the summed (iii) ≥ 2.5 kb Mbp of their six cores: R1 0.71, R2 0.84, R2A 1.29, SF2 1.31 (`coassembly_2026.tsv`). They are not directly comparable, because a co-assembly collapses shared genomes.

## Files

- **Results:**
  - `assembly_vs_jgi.tsv`: per sample. It has the `assembly_vs_bins.tsv` statistics plus N50 ≥ 500 / ≥ 1 kb / ≥ 2.5 kb and Mbp ≥ 5, 10 and 20 kb for (ii) and new, and the (iii) reconstruction. Also the ratios, binned Mbp, trimming and SPAdes statistics, and JGI read recruitment.
  - `assembly_vs_jgi_summary.tsv`: by sample type and site.
  - `mag_source_reconstruction.tsv`, `reconstruction_validation_2026.tsv`, `kmer_cov_pairs.tsv`, `pair_contiguity.tsv`, `contiguity_by_coverage.tsv`, `final_k_overlaps.tsv`, `spades_run_stats.tsv`, `coassembly_2026.tsv`, `merge_rate_R1_A_D1.txt`.
- **Data:** `lengths_2026.tsv.gz` and `lengths_jgi_v2.tsv.gz`, one row per contig.
- **Rerun order:**
  1. `fetch_2026_lengths.py`, `jgi_v2_lengths.py`, `spades_run_stats.py`
  2. `reconstruct_mag_source.py`
  3. `compare_assemblies.py`
  4. `kmer_coverage_check.py`
  5. `pair_contiguity.py`
  6. `contiguity_by_coverage.py`

  All run with `~/Documents/py_venv/bin/python -I`.
