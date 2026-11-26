# Where the current saltern MAGs came from (researched 2026-11-05)

This summarises the research workflow; its full output, with a source for every claim, is in `research_results_full.json`. The tables here were copied from that run's scratchpad.

## Timeline

- **2014-07-11:** R1, R2 and R2A cores (D1 0–5 cm, D2 5–15 cm) and water were sampled.
- **2014-07-17:** SF2 was sampled.
- **2015-03-10:** the three R2_restored samples (C_black, D_shore, H2O) were collected, in a separate later campaign.
- **2015-03-02:** JGI proposal 1890 (PI Tringe) was approved.
- **Dec 2015 and Jan 2016:** the 39 libraries were sequenced on two HiSeq 2500 runs, 10057 (18 libraries) and 10092 (21 libraries). Each library is 2×151 bp.
- **2017-01 to 2021-05:** the reads were released to SRA in waves. They total 344.10 Gbp, 6.1–12.4 Gbp per run (`reads_summary.tsv`).

## Three assemblies per sample

| | Assembly | Made | Use |
|---|---|---|---|
| (i) | JGI standard: MEGAHIT v1.0.3 | 2015-12 to 2016-01 | IMG KO profiles. This is the assembly the manuscript's Metagenomics paragraph describes. None exists for R1_C_H2O. |
| (ii) | JGI IMG-team "SPAdes version 2" reassembly | 2018-03/04 | IMG 3300026094–3300026251, the assemblies behind the Henry-lab CoverM BAMs (`reads/coverm`). |
| (iii) | Tringe lab per-sample reassembly: Tadpole error correction plus SPAdes, version not reported | about 2018 | **The source of all 450 MAG contigs.** It is not public, and none of 129,258 MAG contigs matches (ii) by NODE index and length. |

## Binning

- **Method:** per-sample binning of contigs ≥ 2.5 kb with CONCOCT, MaxBin and MetaBAT, then DAS Tool, then CheckM (completeness > 70 %, contamination < 25 %), then dRep (gANI > 96.5 %, 60 % overlap).
- **Result:** dRep kept 310 MAGs and set 140 aside.
- **Sources of the 310:** 209 from the 24 cores, 64 from the 12 water samples and 37 from the 3 R2_restored samples.
- **Date:** bounded by BBTools "accessed 16 April 2018" in Genes 2021 and the KBase upload on 2020-11-18.

## Read recruitment

Reads were mapped to assembly (ii) (contigs ≥ 500 bp) and competitively to the 450 MAGs (`assembly_vs_bins.tsv`). Each cell gives % of reads to the assembly, then % to the MAGs.

| Site | Cores | Water |
|---|---|---|
| R1 | 53.9 / 52.3 | 73.1 / 39.8 |
| R2 | 41.4 / 38.9 | 85.9 / 56.9 |
| R2A | **15.1** / 9.3 | 40.7 / 23.3 |
| SF2 | **14.2** / 12.0 | 38.0 / 8.7 |

The three R2_restored samples come to 47.1 / 37.2.

**Interpretation, with caveats:**
- **Sediment cores:** in the cores, where methane was measured, the MAGs already recruit nearly every read that assembles. The low binned fraction in the SF2 and R2A cores is therefore an assembly ceiling: about 85 % of their reads never assemble into contigs of 500 bp or more. A better binner alone cannot lift it. Lifting it would take deeper sequencing, long reads, or per-site co-assembly of the replicate cores.
- **Water samples:** assembled reads remain unbinned in the water samples, a gap of 17–33 percentage points, so better (multi-sample) binning can add genomes there.
- **Caveats:**
  - This measures assembly (ii), not the more contiguous MAG-source assembly (iii).
  - MAG recruitment is competitive mapping against MAGs from all 39 samples.
  - 6 of the 39 samples have no BAM.

## Discrepancies with the manuscript

- **Assembly:** the Methods say "assembled with MEGAHIT", which is true only of assembly (i). The MAGs come from (iii).
- **Number of samples:** the MAGs come from 39 metagenomes, not 24.
- **Sampling date:** the R2_restored samples date from 2015-03-10, not July 2014.
- **Sequencing yield:** Zhou et al. (2022) report "~20 Gb per sample, 700 Gb in total", but SRA holds 344 Gbp.
- **Unbinned fractions:** the Table 3 caption ("68 %") and Discussion P96 ("47/58 %") disagree.
- **Citations:** CONCOCT is cited as Lin & Liao 2016 (MyCC) in Zhou et al. and in Genes 2021.
- **Data availability:** the 450 MAGs are not public (KBase ws 77580, 163264 and 155805 are private).
