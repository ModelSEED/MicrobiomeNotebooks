# Tadpole pilot: results

Evaluated 2026-11-12 on poplar. Nothing was committed, and nothing was written to KBase. The setup, the exact commands and the decision rule are in `README.md`.

## Review corrections (2026-11-12)

These corrections follow an adversarial review of this file. The text below already includes them, and the earlier version is kept at `eval/before_review_2026-11-12/RESULTS.md`. No data or script changed; the review re-read the existing tables.

- **Contiguity.**
  - Arm C's median NGA50 on R1_A_D1's genomes is **6.5–9.2 kb**. That is 0.60–0.85× the MAGs' median N50 of 10.8 kb and 2.4–3.4× KBase's 2.7 kb.
  - The earlier single figure (9.2 kb, 0.85×) is the top of that range, and "nearly piece for piece" is dropped.
- **Chimera signal.**
  - 25 of Arm C's 34 cross-MAG contigs involve metabat.6, a *Fodinibius* bin with 12 % contamination. They are therefore more likely mis-binning or strain splits in the original binning than assembly chimeras.
  - The "3–4×" rate is dropped. Per contig the figures are 0.34 % vs 0.04–0.07 %.
- **The 30 lost MAGs.**
  - Fragmentation is at most a partial explanation. 7 of the 30 come from R2_restored_DShore, the *best*-recovered sample (2026 / (iii) = 0.741, rank 39 of 39 counting from the worst).
  - Across thirds of the samples ranked by that ratio, the lost fraction falls only from 9.1 % to 7.1 % to 5.5 % (Denominator: all 450 original MAGs. With the matched denominator, the 275 medium/high-quality dereplicated MAGs the 30 lost ones come from, it is 15.7 %, 13.6 %, 8.2 %: about 1.9x from the worst to the best third. Added 2026-11-12 after verification.)
  - "Their profile is the one Arm C helps most" is dropped. The first rollout batch is now framed as a discriminating test.
- **Rollout gain.**
  - The gain is a single-sample extrapolation from R1_A_D1, which ranks 29 of 39 (above the median). So ≈ 5,650–5,870 Mbp at ≥ 2.5 kb is optimistic.
  - At ≥ 10 kb it is ≈ 2,650–2,740 Mbp (3.25–3.36×), not 2,650–2,790 (3.3–3.4×).
- **Factor attribution.**
  - The "merging" effect includes Tadpole extension inside BBMerge (34 % of R1_A_D1 reads fully extended).
  - Tadpole correction was only tested on unmerged reads under `--meta`. Arm B = KBase shows that it has no *net* effect there.
  - The effects depend on the order of the changes, and there is no Arm D.
  - The mode adds +53 Mbp vs merging's +34 Mbp at ≥ 2.5 kb, and about the same at ≥ 10 kb (+52 vs +49).
  - "(iii) used default-mode SPAdes" is consistent with the data but not shown. Arm C joins more aggressively than (iii): it has 0.82× (iii)'s Mbp in the 2.5–5 kb band and 1.05× at ≥ 20 kb. An older SPAdes `--meta` was not tested.
- **SF2_C_D1 is not a low-depth sample.**
  - It has 9.09 Gbp raw against R1_A_D1's 9.03 Gbp. Its low extension rate and k-mer counts point to high diversity and low per-genome coverage.
  - "The gap is in the assembler, not the input" is dropped.
- **Resources.**
  - "≤ 300–400 GB with two lanes" holds only if the 200 GB Java heap is lowered. The heap is fixed with `-Xms` = `-Xmx`, so it is reserved in full.
  - The Arm A vs Arm C runtime (10.6 vs 9.2 h) is not a like-for-like comparison, because the host load differed.
- **SF2_C_D1 Arm C** was started by the evaluator, outside the original pilot design. The run is described under "Pending", with its PIDs (driver 3634245, chain 3634387, `spades.py` 3634657), the done marker `SF2_C_D1/.done/C1_spades_default` and the evaluation commands.
  - At 18:41 CST on 2026-11-12, all three processes were alive (1 h 08 min elapsed) and the marker did not exist yet.
  - `.done/ALL` already exists from the Arm A/B run, so wait for the `C1_spades_default` marker, not for `ALL`.

## Short answer

**For R1_A_D1, yes: the (iii) assembly can be regenerated at about the right size and with most of its contiguity.** The recipe is:

- merged and Tadpole-corrected reads;
- SPAdes in **default (multicell) mode**, not `--meta`.

That is Arm C. On R1_A_D1 it holds 0.93× (iii)'s Mbp in contigs ≥ 2.5 kb and 1.02× in contigs ≥ 10 kb, so it passes the 85 % go rule.

For the same original genomes, its median NGA50 is 6.5–9.2 kb. That is 0.60–0.85× the MAGs' median N50 of 10.8 kb and 2.4–3.4× KBase's (section 2). This is one sample; it ranks 29 of 39 by 2026 / (iii), above the median.

**Which factor matters (R1_A_D1, Mbp in contigs ≥ 2.5 kb / ≥ 10 kb):**

| Change | From → to | Effect |
|---|---|---|
| Tadpole correction instead of BayesHammer, with BBDuk instead of Trimmomatic (tested only on unmerged reads, under `--meta`) | 2026 KBase → Arm B | **no net effect**: 162.8 / 32.6 → 164.1 / 32.5 |
| Pair merging, including Tadpole extension inside BBMerge (`bbmerge-auto rem extend2 ecct`; 34 % of reads fully extended) | Arm B → Arm A | +34 Mbp / **2.5×** (32.5 → 81.5); +49 Mbp at ≥ 10 kb |
| SPAdes default mode instead of `--meta` | Arm A → Arm C | +53 Mbp / **1.6×** (81.5 → 133.8); +52 Mbp at ≥ 10 kb |

Each effect is measured on top of the change before it, so the attribution depends on the order. There is no Arm D (default mode on unmerged reads).

- **Merging is what makes the reads look like (iii)'s.** For the same genome, SPAdes k-mer coverage in Arm A and Arm C equals (iii)'s (ratio 0.96–1.00). Without merging, in Arm B and in the 2026 KBase run, it is 0.34–0.41. So (iii) was assembled from merged reads. The merging step here also extends reads with Tadpole, and the pilot does not separate the two.
- **Merging alone does not close the gap.**
  - At ≥ 2.5 kb the mode adds more than merging (+53 vs +34 Mbp). At ≥ 10 kb the two add about the same (+52 vs +49 Mbp).
  - That (iii) was run in default mode is consistent with this, but not shown. Arm C joins more aggressively than (iii): it holds 0.82× (iii)'s Mbp in the 2.5–5 kb band but 1.05× at ≥ 20 kb.
  - An older SPAdes `--meta`, such as the 3.11–3.12 likely used in 2018, was not tested.
- **It is not a bit-for-bit regeneration.** Only 9.7 % of (iii)'s MAG contigs reappear end to end in Arm C (at ≥ 99.9 % identity). Some detail still differs: the SPAdes version (likely 3.11–3.12 in 2018), the BBTools version or parameters, or the k list.

**SF2_C_D1, the worst-recovered sample, is not resolved.**

- Arm A reached only 0.38× / 0.21× of (iii).
- Arm C was not part of the original design for this sample. The evaluator launched it on 2026-11-12 at 17:33 CST, and it is still running (see "Pending").
- So the go decision holds for R1-type samples. For the six SF2-type samples (2026 / (iii) ≤ 0.45), it waits on this run.

## 1. Assembly statistics

`scaffolds.fasta` (what KBase saved), lengths with N's. `contigs.fasta` differs by < 1 % in every arm (`eval/pilot_assemblies.tsv`). The (iii) values are the reconstruction from MAG contig names, which leans slightly low.

**R1_A_D1** (SRR5209943, 9.03 Gbp; 2026 / (iii) = 0.61, rank 29 of 39 from the worst, a little above the 39-sample median of 0.57)

| Assembly | Mbp ≥ 1 kb | ≥ 2.5 kb | ≥ 5 kb | ≥ 10 kb | ≥ 20 kb | N50 of ≥ 2.5 kb | Largest | Contigs ≥ 1 kb | k-mer cov.* |
|---|---|---|---|---|---|---|---|---|---|
| (iii) target | – | **268.9** | 189.0 | **131.4** | 87.2 | 9,579 | – | – | – |
| 2026 KBase (metaSPAdes, BayesHammer) | 354.1 | 162.8 | 78.6 | 32.6 | 12.4 | 4,843 | 177,593 | 164,365 | – |
| B: Tadpole, metaSPAdes | 361.9 | 164.1 | 78.0 | 32.5 | 12.7 | 4,791 | 177,366 | 169,613 | 1.27 |
| A: merge + Tadpole, metaSPAdes | 364.5 | 197.7 | 126.4 | 81.5 | 51.1 | 7,294 | 401,108 | 145,847 | 4.42 |
| **C: merge + Tadpole, default SPAdes** | **393.9** | **250.5** | **184.6** | **133.8** | **91.9** | **11,499** | **458,402** | 129,813 | 3.56 |
| C / (iii) | – | **0.93** | 0.98 | **1.02** | 1.05 | 1.20 | | | |
| A / (iii) | – | 0.74 | 0.67 | 0.62 | 0.59 | 0.76 | | | |

\* Length-weighted median SPAdes k-mer coverage (k = 127) of contigs ≥ 2.5 kb.

**SF2_C_D1** (SRR5210468, 9.09 Gbp; 2026 / (iii) = 0.31, the worst of the 39)

| Assembly | Mbp ≥ 1 kb | ≥ 2.5 kb | ≥ 5 kb | ≥ 10 kb | ≥ 20 kb | N50 of ≥ 2.5 kb | Largest | Contigs ≥ 1 kb | k-mer cov. |
|---|---|---|---|---|---|---|---|---|---|
| (iii) target | – | **59.3** | 31.1 | **18.2** | 9.0 | 5,344 | – | – | – |
| 2026 KBase | 76.9 | 18.3 | 6.7 | 2.0 | 0.5 | 4,029 | 57,307 | 46,758 | – |
| B: Tadpole, metaSPAdes | 79.8 | 18.7 | 6.6 | 2.1 | 0.6 | 3,903 | 36,850 | 48,646 | 0.91 |
| A: merge + Tadpole, metaSPAdes | 80.4 | 22.5 | 8.8 | 3.8 | 1.9 | 4,148 | 95,434 | 46,763 | 2.58 |
| C: default SPAdes | *running* | | | | | | | | |
| A / (iii) | – | 0.38 | 0.28 | 0.21 | 0.22 | 0.78 | | | |

**No loss of low-coverage sequence is visible at the 1 kb floor in default mode.** Arm C holds *more* sequence in contigs ≥ 1 kb than Arm A (393.9 vs 364.5 Mbp), in fewer contigs. Not tested: sequence below 1 kb (KBase never saved it), and genomes below binning depth (the genome check uses only binned genomes).

### Read preparation (full data; `eval/pilot_prep.tsv`)

| Step | R1_A_D1 | SF2_C_D1 |
|---|---|---|
| BBDuk adapters (rqcfilter2 defaults): bases kept | 95.57 % | 85.86 % (12.9 % of reads adapter-trimmed; Trimmomatic kept 84.2 %) |
| BBDuk artifacts and phiX: reads removed | 914 (0.00 %) | 376,504 (0.72 %) |
| BBMerge `rem k=62 extend2=50 ecct`: pairs joined | **62.5 %** (18.0 M of 28.8 M) | **47.9 %** (12.4 M of 25.9 M) |
| ambiguous / no solution | 37.5 % / 0.01 % | 52.1 % / 0.01 % |
| reads fully extended by Tadpole inside BBMerge | 34.2 % | **2.8 %** |
| mean insert (mode) | 267.7 bp (277) | 248.0 bp (272) |
| non-singleton 62-mers (BBMerge table) | 357 M | 87 M |
| Tadpole, Arm A unmerged pairs: reads with errors found; bases corrected | 10.1 %; 2.91 M | 1.7 %; 0.53 M |
| Tadpole, Arm A merged reads | 10.1 %; 1.41 M | 3.1 %; 0.25 M |
| Tadpole, Arm B pairs | 11.5 %; 8.49 M | 2.1 %; 1.20 M |
| quality trim after correction (`qtrim=r trimq=10`): bases kept, A / B | 98.5 % / 99.0 % | 99.2 % / 99.4 % |
| SPAdes mean read length, A / B | 202.0 / 148.2 | 179.8 / 148.9 |
| SPAdes insert size (all arms) | 301–302 ± 67 bp | 295–296 ± 65 bp |

**Merging lifted the rate from the subsample estimate.** On the full data, extension merges 62.5 % of R1_A_D1 pairs, against the 45 % lower bound from 300,000 pairs.

**SF2_C_D1 is not a low-depth sample. It looks like a high-diversity one.**
- It has as much raw data as R1_A_D1 (9.09 vs 9.03 Gbp).
- Tadpole's extension needs per-genome k-mer depth, and it fully extended only 2.8 % of reads. The BBMerge tables hold a quarter of R1_A_D1's 62-mers. Both point to many genomes, each at low coverage.
- Arm A's k-mer coverage for SF2_C_D1's own genomes equals (iii)'s (section 2). So (iii) most likely had no more reads for this sample than SRA does.
- Whether the remaining gap is the assembler mode is what the pending Arm C run tests.

### Runtimes and peak memory (`eval/pilot_resources.tsv`)

Measured at nice 10 on poplar with a host load of about 320–560. BBTools steps used 24 threads (12 for B1 and B2), and each SPAdes run used 12.

| Step | R1_A_D1 wall | SF2_C_D1 wall | Peak memory, R1 / SF2 |
|---|---|---|---|
| separator fix (awk) | 0.10 h | 0.10 h | – |
| BBDuk, two passes | 0.22 h | 0.18 h | < 1 GB used (16 GB heap) |
| BBMerge with extension | 0.59 h | 0.43 h | 68 / 53 GB used (200 GB heap) |
| Tadpole, two passes (Arm A) | 1.10 h | 0.96 h | 107 / 106 GB used |
| quality trim | 0.02 h | 0.02 h | – |
| **Preparation, Arm A / C path** | **2.0 h** | **1.7 h** | |
| Tadpole (Arm B), 12 threads | 0.65 h | 0.72 h | 147 / 144 GB used |
| metaSPAdes, Arm A | 10.6 h | 14.3 h† | 26 / 53 GB RSS |
| metaSPAdes, Arm B | 10.7 h | 14.8 h | 36 / 55 GB RSS |
| **default SPAdes, Arm C** | **9.2 h** | running | **37** / – GB RSS |

- "Used" is the largest Java heap figure in each step log. The heaps were fixed at `-Xms` = `-Xmx` = 200 GB, so that is what each sample reserved.
- † Relaunching the driver for SF2_C_D1's Arm C truncated `logs/armA_chain.out`, which held Arm A's start and done lines (see "Gotchas"). The start comes from the driver's "Arm A chain" line (12:46:44), and the end is the mtime of `.done/A4_metaspades` (03:04:52).
- Arm C took 9.2 h and Arm A 10.6 h on the same reads. This is not a like-for-like comparison: the runs did not overlap, and the host load differed.
- SF2_C_D1's assemblies took about 35–40 % longer than R1_A_D1's, although SF2_C_D1 assembles to a quarter of the size.

**Disk.**

| Item | R1_A_D1 | SF2_C_D1 |
|---|---|---|
| raw FASTQ (fasterq-dump, uncompressed) | 20.7 GB | 20.8 GB |
| `qc/` | 11 GB | 11 GB |
| `armA/` reads | 13 GB | 12 GB |
| `armA/` SPAdes inputs (`reads_{1,2}`, `ecc_merged`) | 5.4 GB | 4.9 GB |
| SPAdes output directory at the end | 9.4–22 GB | 4.8–17 GB |

- `qc/` holds `trim_*` and `clean_*`. The `trim_*` files are deletable after the second BBDuk pass.
- SPAdes k-mer counting fills `tmp/` temporarily: Arm C on SF2_C_D1 had 77 GB in `tmp/` 13 min after its start.
- The pilot directory now holds 77 GB (R1_A_D1), 54 GB (SF2_C_D1, plus Arm C) and 39 GB of reads.

## 2. Genome-level check: the same genomes in each assembly

**Method** (`scripts/genome_check.py`, `scripts/make_mag_list.py`).

- **Query:** the contigs of every original MAG binned from these samples' own (iii) assemblies.
  - R1_A_D1: 11 derep310 + 12 nonderep140 MAGs = 23 MAGs, 73.1 Mbp, 8,497 contigs.
  - SF2_C_D1: 1 + 4 = 5 MAGs, 17.7 Mbp, 1,955 contigs.
  - Every MAG contig is a (iii) contig ≥ 2.5 kb, so the MAG itself is the (iii) reference for contiguity.
- **Target:** each assembly cut to sequences ≥ 1 kb (KBase's floor). The 2026 KBase assemblies were fetched with `local/fetch_assembly.py`.
  - Their md5 matched the Blobstore node, the handle and the Assembly object.
  - Their per-contig md5s matched too: R1_A_D1 `275645/238/1`, md5 `df2841aa…`; SF2_C_D1 `275645/283/1`, md5 `14da8978…`.
- **Mapping:** `minimap2 -x asm10 -c --secondary=no`. Identity is matches / alignment block length; blocks < 500 bp are ignored.
- **Metrics:**
  - **Covered:** the fraction of MAG bp covered at ≥ 99 % (and ≥ 95 %) identity.
  - **Binnable:** the same at ≥ 99 %, on target contigs ≥ 2.5 kb (the MetaBAT2/CONCOCT floor).
  - **Pieces:** target contigs carrying ≥ 1 kb of the MAG at ≥ 99 %.
  - **NGA50:** the pieces sorted by the MAG bp they carry; the value at 50 % of MAG size. Compare it with the MAG's own N50.
  - **Identical:** MAG contigs reproduced end to end, both ends within 5 bp, at ≥ 99.9 % identity.
  - **k-mer coverage ratio:** target / (iii), from the SPAdes contig names; length-weighted median over alignments ≥ 5 kb.
  - **Cross-MAG joins:** target contigs carrying ≥ 1 kb from two *different* original MAGs. DAS Tool bins are disjoint, so these are candidate chimeras.

**R1_A_D1, 23 original MAGs** (median MAG N50 10.8 kb; 8,497 (iii) contigs)

| Assembly | Covered ≥ 99 % | Covered ≥ 95 % | **Binnable** | **Pieces** | **Median NGA50** | Identical contigs | k-mer cov. ratio | Cross-MAG joins |
|---|---|---|---|---|---|---|---|---|
| 2026 KBase | 78.5 % | 87.1 % | 51.8 % | 20,897 | 2,712 | 1.6 % | 0.39 | 8 |
| B | 78.3 % | 87.3 % | 51.1 % | 21,160 | 2,623 | 2.3 % | 0.34 | 9 |
| A | 81.9 % | 90.1 % | 60.6 % | 18,463 | 3,036 | 3.5 % | 0.99 | 13 |
| **C** | **88.6 %** | **95.8 %** | **83.2 %** | **9,940** | **9,185** | 9.7 % | 0.96 | 34 |
| (iii) MAGs themselves | 100 % | 100 % | 100 % | 8,497 | 10,805 (N50) | 100 % | 1 | – |

- The derep310 subset alone (11 MAGs, 3,900 contigs, median N50 12.8 kb) gives the same picture.
  - KBase: 59.0 % binnable, 9,505 pieces, NGA50 2.9 kb.
  - Arm C: 86.3 % binnable, 4,731 pieces, NGA50 9.2 kb.
- Measuring contiguity at ≥ 95 % identity, which tolerates strain-consensus differences, changes little: Arm C's median NGA50 is 9.6 kb against the MAGs' 10.8 kb.

**Per genome: merging helps only well-covered genomes, while default mode helps at every coverage.** The MAGs are split at a (iii) k-mer coverage of 2.2, which gives groups of 10 and 13 (the median is 2.6):

| R1_A_D1 MAGs | Median MAG N50 | Binnable: KBase / B / A / **C** | Median NGA50: KBase / B / A / **C** |
|---|---|---|---|
| (iii) k-mer cov. < 2.2 (n = 10) | 6.4 kb | 45 / 45 / 47 / **76 %** | 2.0 / 2.1 / 2.3 / **4.4 kb** |
| (iii) k-mer cov. ≥ 2.2 (n = 13) | 19.6 kb | 56 / 55 / 70 / **88 %** | 3.2 / 3.0 / 4.3 / **18.4 kb** |

**Examples: for well-covered genomes Arm C comes close to (iii)'s piece count and NGA50.**
- Each cell is pieces / NGA50, against the MAG's own contigs / N50.
- A piece can carry several MAG contigs (the per-MAG `joins` column), so equal counts do not mean contig-for-contig equivalence.

| MAG ((iii) k-mer cov.) | (iii) MAG | 2026 KBase | Arm A | **Arm C** |
|---|---|---|---|---|
| metabat.45 (6.1) | 160 / 41.7 kb | 724 / 3.9 kb | 381 / 11.1 kb | **171 / 41.7 kb** |
| metabat.47 (13.5) | 141 / 29.6 kb | 776 / 3.3 kb | 224 / 24.9 kb | **154 / 29.6 kb** |
| metabat.51 (7.5) | 64 / 75.2 kb | 229 / 18.9 kb | 117 / 35.5 kb | **71 / 71.4 kb** |
| concoct_out.3 (3.9) | 140 / 35.0 kb | 751 / 4.1 kb | 626 / 5.2 kb | **150 / 31.4 kb** |
| concoct_out.36 (1.4) | 585 / 5.9 kb | 1,010 / 1.7 kb | 1,110 / 1.7 kb | 682 / 4.0 kb |

**SF2_C_D1, 5 original MAGs** (median MAG N50 14.4 kb; 1,955 contigs). Arm C pending.

| Assembly | Covered ≥ 99 % | ≥ 95 % | Binnable | Pieces | Median NGA50 | k-mer cov. ratio |
|---|---|---|---|---|---|---|
| 2026 KBase | 57.9 % | 75.1 % | 26.5 % | 4,859 | 1,273 | 0.41 |
| B | 57.0 % | 75.7 % | 25.0 % | 4,879 | 1,139 | 0.35 |
| A | 62.3 % | 82.0 % | 34.8 % | 4,745 | 1,464 | **1.00** |

- Even Arm A covers only 62 % of SF2_C_D1's original genomes at ≥ 99 % identity, and only 35 % on binnable contigs.
- Yet its k-mer coverage on matched contigs equals (iii)'s, so (iii) did not have more reads for these genomes. The pending Arm C run tests whether the assembler mode closes the gap.

**Answer to "does any arm reproduce (iii)'s contiguity for the same genomes?"**

- **Arm C comes closest on R1_A_D1, but does not reproduce (iii).**
  - It needs 1.17× as many pieces as (iii) had contigs. KBase needs 2.46× and Arm A 2.17×.
  - Its median NGA50 is 6.5–9.2 kb, 0.60–0.85× the MAG N50 of 10.8 kb. The upper value is `genome_check.py`'s 9,185 bp; the lower is the review's estimate.
    - KBase reaches 2.7 kb (0.25×) and Arm A 3.0 kb (0.28×). So Arm C is 2.4–3.4× KBase.
  - It puts 83 % of the original genomes on binnable contigs, against 52 % for KBase.
- **It comes closest for the best-covered genomes.** At (iii) k-mer coverage ≥ 6 (metabat.44, .45, .47, .51), Arm C's NGA50 is 0.95–1.09× the MAG N50. At coverage 3.5–5 it is 0.62–0.94×.
- **For low-coverage genomes Arm C doubles KBase's contiguity** but stays below (iii) (NGA50 4.4 vs 6.4 kb).
- **Arms A and B do not reproduce (iii).**

**Cross-MAG contigs.** Arm C has more contigs that carry ≥ 1 kb from two different original MAGs:

- Arm C: 34 contigs. That is 0.34 % of the contigs carrying MAG sequence and 0.8 % of their bp.
- Other assemblies: 8–13 contigs, 0.04–0.07 % of contigs and about 0.1 % of bp.

**These are mostly not assembly chimeras.**
- 25 of the 34 involve metabat.6, a *Fodinibius* bin with 12 % contamination.
- So most of them more likely reflect mis-binning or strain splits in the original binning.
- No chimera-rate multiplier is claimed.

A chimerism check (GUNC or CheckM2) on rollout MAGs is still prudent.

### The 30 "lost" original MAGs

**None of the 30 comes from R1_A_D1 or SF2_C_D1.** These are the medium- or high-quality original MAGs with no 2026 MAG above 80 % ANI (`comparison/README.md`), reproduced by `make_mag_list.py`.

- **By source sample:**
  - R2_restored_DShore (sediment): 7;
  - SF2_A_H2O, R1_B_H2O and R1_C_H2O: 3 each;
  - 11 other samples: 1–2 each (14 MAGs).
  - 13 of the 30 come from water samples, 7 from R2_restored_DShore and 10 from other sediments. Both pilot samples are sediments.
- **The pilot cannot measure their recovery directly.** Mapped against both pilot samples, their genomes are essentially absent:
  - R1_A_D1: 2.0–2.2 % of their bp at ≥ 99 % identity, in every arm and in KBase.
  - SF2_C_D1: 0.3–0.4 %.
- **One partial exception:** the R1_A_D2 UBA9300 MAG. 83 % of it aligns at ≥ 95 % identity in R1_A_D1 Arm C, but only 12 % lies on binnable contigs at ≥ 99 %. That is a related strain at low coverage.
- **Next most present:** two DSBN01, a Natrinema and a CAZUKV01 MAG (R1 sediment and water), each about 40–45 % present at ≥ 95 % in Arm C.

**Fragmentation is at most a partial explanation.**

- **The source samples argue against it.**
  - 7 of the 30 come from R2_restored_DShore, which is the *best*-recovered sample (2026 / (iii) = 0.741 at ≥ 2.5 kb, rank 39 of 39 counting from the worst).
  - Ranking the samples by that ratio, the lost fraction of original MAGs falls only weakly from the worst third to the best: 9.1 %, 7.1 %, 5.5 % (Denominator: all 450 original MAGs. With the matched denominator, the 275 medium/high-quality dereplicated MAGs the 30 lost ones come from, it is 15.7 %, 13.6 %, 8.2 %: about 1.9x from the worst to the best third. Added 2026-11-12 after verification.)
- The lost-30 MAGs have a median (iii) k-mer coverage of 2.1 (IQR 1.6–3.0) and a median N50 of 8.6 kb. That resembles R1_A_D1's low-coverage half (above), where Arm C raised the binnable fraction from 45 % to 76 %.
  - But that alone does not show they were lost to fragmentation.
- The pilot cannot say how many a rollout would bring back. The first rollout batch is designed as the discriminating test (section 3).

## 3. Conclusions and rollout

**Can the Tringe-lab assemblies be regenerated?**

- **Functionally yes, for R1_A_D1.**
  - Arm C reaches 93–105 % of (iii) at every length threshold. It recovers the same genomes at 0.60–0.85× (iii)'s contiguity (median NGA50 6.5–9.2 vs 10.8 kb).
  - This is one sample, which ranks 29 of 39 by 2026 / (iii), above the median. It is not shown to be typical.
- **Exactly, no:** under 10 % of (iii) contigs are reproduced end to end.
- **The factors:**
  - **The SPAdes mode and the merging both matter.**
    - At ≥ 10 kb they add about the same (+52 vs +49 Mbp). At ≥ 2.5 kb the mode adds more (+53 vs +34 Mbp).
    - Each was measured on top of the previous change, so the split depends on the order. There is no Arm D (default mode without merging).
    - The "merging" step includes Tadpole extension inside BBMerge (34 % of R1_A_D1 reads fully extended).
  - **Swapping BayesHammer for Tadpole and Trimmomatic for BBDuk has no net effect on unmerged reads under `--meta`:** Arm B equals the KBase run. Tadpole correction was not tested on its own with merging or in default mode.
  - The k-mer coverage test shows that (iii) used merged reads.
  - **That (iii) used default-mode SPAdes is consistent with the data, not shown.**
    - Arm C joins more aggressively than (iii): it holds 0.82× (iii)'s Mbp in the 2.5–5 kb band and 1.05× at ≥ 20 kb.
    - An older SPAdes `--meta` (3.11–3.12) was not tested.
    - The k = 127 overlap found earlier fits either mode.
- **For the worst SF2 sample**, merging alone closes almost none of the gap (0.31 → 0.38). Whether default mode closes the rest is the open question (pending run).

**Decision rule** (README: ≥ 85 % of (iii) at both ≥ 2.5 kb and ≥ 10 kb):

- R1_A_D1 **passes with Arm C** (0.93 / 1.02) and fails with Arm A (0.74 / 0.62).
- SF2_C_D1 fails with Arm A (0.38 / 0.21); Arm C is pending.
- **Recommendation:** if SF2_C_D1 Arm C also reaches about ≥ 0.85, roll Arm C out to all 39 samples.
- If it improves SF2_C_D1 substantially but stays below 0.85, roll out anyway. The 33 samples outside the SF2-type group carry 94 % of (iii)'s ≥ 2.5 kb sequence (5,897 of 6,300 Mbp). Ask the Tringe lab for (iii) only for the SF2-type samples.

**What a 39-sample rollout of Arm C would cost on poplar.** Arm A and Arm B assemblies are not needed.

| Item | Per sample | 39 samples |
|---|---|---|
| Download (NCBI ODP mirror) + vdb-validate + fasterq-dump | about 0.3–0.5 h | about 15 h, can run ahead |
| Preparation: BBDuk ×2, BBMerge, Tadpole ×2, trim (24 threads) | 1.7–2.0 h | about 75 h of 24-thread lanes |
| Default SPAdes, 12 threads | 9.2 h (R1-type); SF2-type probably 12–15 h (metaSPAdes took 14.3 h) | about 390–590 h of 12-thread lanes (9.2–15 h each) |
| **Wall time** | about 11–17 h | **about 8–12 days with two 12-thread SPAdes lanes plus preparation at 8 threads (≤ 32 threads)**; about 6–8 days with three lanes (≤ 48 threads, as in the pilot) |
| **Memory** | BBTools heap use ≤ 107 GB on the Arm C path. The 'used' figure is noisy, though: Arm B's single Tadpole pass on the same reads reported 147 GB. `-Xmx150g` should cover a 9 Gbp sample; keep 200 GB for the largest samples (up to 12.4 Gbp), or lower it only after one more full-size run. SPAdes ≤ 37–55 GB | ≤ 300–400 GB with two lanes **only if the Java heap is lowered** from the current 200 GB. The driver sets `-Xms` = `-Xmx`, so each BBTools step reserves its full heap. Poplar showed 2.4 TB available. |
| **Disk, in flight** | ≈ 20.7 GB raw FASTQ + 11 GB BBDuk + 13 GB merged/corrected reads + up to ≈ 80 GB SPAdes `tmp/` ≈ 125 GB peak | ≈ 250–375 GB with 2–3 lanes |
| **Disk, kept** | ≈ 5 GB SPAdes inputs + ≈ 2–3 GB of FASTA and graphs after deleting `K*/` and `tmp/`. The raw FASTQ can go after BBDuk (re-downloadable). | ≈ 300 GB; /scratch1 has 12 TB free (83 % used, shared) |
| Then: re-binning, CheckM, GTDB-Tk (existing pipeline from the binner stage on) | | about 3–4 days (assembly_vs_jgi README estimate) |

**What it would gain.**

| Measure | Now | Projected | Basis |
|---|---|---|---|
| Mbp in contigs ≥ 2.5 kb, 39 samples | 3,691 | **≈ 5,650–5,870 (+53–59 %), optimistic** | A single-sample extrapolation from R1_A_D1, which ranks 29 of 39 by 2026 / (iii), above the median. Upper: every sample at R1_A_D1's 0.93× of (iii) (6,300 Mbp). Lower: the six SF2-type samples unchanged (159 Mbp now). |
| Mbp in contigs ≥ 10 kb | 815 | **≈ 2,650–2,740 (3.25–3.36×)** | the same scenarios at 1.0× of (iii) (2,740 Mbp; SF2-type samples 26 Mbp now) |
| Original genomes on binnable contigs (R1_A_D1) | 52 % | 83 % | genome check |
| Median NGA50 of the same genomes (R1_A_D1) | 2.7 kb | 6.5–9.2 kb | genome check (9.2 kb) and the review's lower estimate; (iii) MAG N50 10.8 kb |
| 30 lost MAGs | 0 recovered | not measurable here | fragmentation is at most a partial explanation (section 2). The first batch tests it. |

**Before or at the start of the rollout:**

1. Wait for SF2_C_D1 Arm C (see "Pending").
2. **Make the first batch a discriminating test of the lost MAGs.** Assemble two samples with Arm C and bin them:
   - **R2_restored_DShore** (7 of the 30). It is already the best-recovered sample by bulk sequence (0.741).
     - If Arm C brings its 7 back, the loss was in the assembly despite good bulk recovery.
     - If it does not, fragmentation does not explain them, and absence, strain heterogeneity or binning become the candidates.
   - **One water sample** (SF2_A_H2O, 3 of the 30). Water samples are untested by this pilot.
3. **Driver changes for the rollout.** Leave the running copy alone; work on a copy of `run_sample.sh`.
   - Arm C only.
   - `XMX_BIG=150g` for samples up to about 9 Gbp, and 200 GB for the larger ones, until a full-size run confirms the lower value.
   - `fix_plus` as `awk … > tmp && mv`, checking both mates (README review point 3).
   - Delete `trim_*`, `merged`/`unmerged`/`ecc_u_*` and SPAdes `K*/`/`tmp/` when each step finishes.
   - A per-launch name for `armA_chain.out` (see "Gotchas").
4. Check every rolled-out MAG for chimeras (GUNC or CheckM2) before it replaces a 2026 MAG. Arm C has more cross-MAG contigs (34 vs 8–13), although 25 of the 34 involve one contaminated original bin (metabat.6).
5. Optional, for the methods text: default-mode SPAdes on Arm B's unmerged reads (Arm D, about 10 h). It would separate the mode effect from merging, since the factor effects above depend on the order of the changes. The rollout does not depend on it.

## Pending: SF2_C_D1 Arm C

| Item | Value |
|---|---|
| Command | `run_sample.sh SF2_C_D1 SRR5210468 armC`. Every finished step is skipped, so only `C1_spades_default` runs: the same `spades.py --only-assembler -k 33,55,77,99,127 -t 12 -m 400` on Arm A's reads, as for R1_A_D1. |
| Started by | the evaluator, outside the original pilot design |
| Started | 2026-11-12 17:33:52 CST |
| Status at 18:41 CST | all three PIDs alive (1 h 08 min elapsed); `.done/C1_spades_default` not yet present |
| Done marker | `SF2_C_D1/.done/C1_spades_default`. `.done/ALL` already exists from the Arm A/B run, so do not wait on it. |
| Driver PID | 3634245 (`SF2_C_D1/driver.pid`) |
| Chain PID | 3634387 |
| SPAdes PID | 3634657 |
| Log | `SF2_C_D1/driver.armC.log` |
| Previous driver log | `driver.armAB.log`, `driver.armAB.pid` (copies) |
| Expected end | about 06:00–09:00 CST on 2026-11-13, if like metaSPAdes on this sample (14.3 h) |

**When `SF2_C_D1/.done/C1_spades_default` exists:**

```
cd /scratch1/afreiburger/salterns_2026/tadpole_pilot
PY=/scratch1/afreiburger/salterns_2026/envs/spades3153/bin/python
nice -n 10 $PY -I scripts/eval_assemblies.py --root $PWD --targets scripts/targets.tsv --out-prefix $PWD/eval/pilot
nice -n 10 $PY -I scripts/genome_check.py --root $PWD --mags scripts/genome_check_mags.tsv \
    --mag-dir eval/genome_check/mags --kbase-dir eval/genome_check/kbase2026 --out eval/genome_check \
    --minimap2 /scratch1/afreiburger/metag-hybrid/env/miniforge/envs/assembly/bin/minimap2 --threads 8
# locally:
rsync -a -e "ssh -o BatchMode=yes" --exclude='genome_check/asm' --exclude='genome_check/paf' \
    --exclude='genome_check/mags' --exclude='genome_check/kbase2026' \
    poplar:/scratch1/afreiburger/salterns_2026/tadpole_pilot/eval/ eval/
```

- `genome_check.py` picks up `armC` by itself once the `.done` marker exists. It reuses the cached PAFs and maps only the new assembly.
- Then fill in the SF2_C_D1 rows above and apply the decision rule.

## Gotchas found during the evaluation

- **Relaunching `run_sample.sh` truncates `logs/armA_chain.out`.** The script redirects the Arm A/C chain into it with `>`. SF2_C_D1's Arm A start and done lines were lost that way. `eval_assemblies.py` now falls back to the driver's "Arm A chain" line and to the `.done/A4_metaspades` mtime.
- **The relaunch also overwrites `driver.pid` and `STATUS`.** This is harmless: `STATUS` reads `C1_spades_default` and returns to `ALL_DONE` at the end.
- **poplar's NFS file mtimes run about 3–4 min behind its `date`.** Wall times here come from the driver's own timestamps, except the one fallback above.

## Files

Local, under `annotation_runs/tadpole/`; mirrored in poplar's `tadpole_pilot/scripts/` and `eval/`.

- **`scripts/eval_assemblies.py`**, extended here with `<prefix>_resources.tsv`: wall time per step from the driver logs, and peak memory per step (largest Java heap used, or SPAdes peak RSS).
- **`scripts/make_mag_list.py`** writes `scripts/genome_check_mags.tsv`: 58 MAGs, the 28 from the pilot samples plus the lost 30. Run it locally with `~/Documents/py_venv/bin/python -I`.
- **`scripts/genome_check.py`**: the genome-level check, stdlib only.
- **`eval/pilot_assemblies.tsv`, `eval/pilot_prep.tsv`, `eval/pilot_resources.tsv`, `eval/logs/<sample>/`**: copies of the step logs and driver logs.
- **`eval/genome_check/genome_check_mags_by_target.tsv`** (per MAG × assembly), **`genome_check_summary.tsv`** (per sample × assembly × population) and `run.log`.
- **poplar only** (large): `eval/genome_check/asm/` (≥ 1 kb targets and query FASTA), `paf/`, `mags/` (the 58 MAG FASTAs) and `kbase2026/` (the two KBase assemblies with their fetch logs).
