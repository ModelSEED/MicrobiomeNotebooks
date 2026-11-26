# Salterns MAGs: original study vs 2026 re-assembly

Both sets: same 39 metagenomes, KBase CheckM (lineage_wf) and GTDB-Tk 2.8.0 / GTDB R232. The original 310 MAGs are dereplicated; the 1,324 new MAGs are not (39 per-sample assemblies + 4 per-site sediment co-assemblies), so species clusters (skani, ANI >= 95 %, aligned fraction >= 50 %) are the fair comparison.

## 1. Numbers and quality

|                        | original (310, dereplicated)   | re-assembly (all 1,324)   | re-assembly, best MAG per species cluster   |
|:-----------------------|:-------------------------------|:--------------------------|:--------------------------------------------|
| MAGs                   | 310                            | 1324                      | 558                                         |
| species clusters       | 296                            | 558                       | 558                                         |
| high quality (>=90/<5) | 130                            | 323                       | 177                                         |
| medium quality         | 145                            | 1001                      | 381                                         |
| below MQ (<50 or >=10) | 35                             | 0                         | 0                                           |
| median completeness %  | 92.0                           | 83.2                      | 84.5                                        |
| median contamination % | 3.48                           | 3.51                      | 2.88                                        |
| median contigs         | 241                            | 475                       | 392                                         |
| median N50 (kb)        | 19.1                           | 6.4                       | 8.4                                         |
| median size (Mb)       | 3.27                           | 2.78                      | 2.9                                         |
| total size (Gb)        | 1.02                           | 3.72                      | 1.62                                        |
| Bacteria / Archaea     | 280 / 30                       | 1208 / 116                | 511 / 47                                    |
| with a GTDB species    | 176                            | 802                       | 226                                         |
| no GTDB genus          | 20                             | 74                        | 52                                          |

Species clusters across both sets: 604; in both: 250, original only: 46, re-assembly only: 308.

For the 250 species clusters found by both, the best re-assembly MAG vs the best original MAG: better quality score (completeness - 5 x contamination) in 105, worse in 144; median change in completeness -1.9 points; median N50 ratio 0.50x.

Original species clusters NOT recovered by the re-assembly (58 MAGs): Aquiluna (83/20.4), Aquiluna (85/24.2), Aquiluna (84/7.8), Aquiluna (83/3.4), Pontimonas (87/2.4), S36-B12 (89/4.3), S36-B12 (80/3.3), Kariarchaeum (87/4.2), Halanaerobium (88/11.8), Halanaerobium (80/20.2), Halanaerobium (76/5.5), Christiangramia (87/2.9), Draconibacterium (88/19.7), Fodinibius (85/11.2), Fodinibius (93/7.7), MAG-120531 (86/12.5), MAG-120531 (92/20.3), MAG-120531 (82/3.1), Psychroflexus (99/0.6), Rhodohalobacter (87/20.8), UBA3478 (87/19.2), UBA3478 (77/10.2), UBA9300 (78/3.9), VMCP01 (93/3.9), Yeosuana (83/9.0), Sulfurovum (77/3.5), Sulfurovum (96/15.1), Humimonas (93/8.3), Hikarchaeaceae (82/8.5), Haladaptataceae (80/4.3), Haloarcula (87/23.0), Halonotius (85/16.0), Halonotius (95/12.6), Natrinema (82/4.2), CAZUKV01 (77/2.9), D2472 (84/8.6), DASQZJ01 (100/5.2), DATFLG01 (85/2.8), Guyparkeria (88/1.4), Halofilum (92/5.0), JBBDSH01 (84/4.9), Marinobacterium_E (93/3.8), Planktomarina (84/3.1), Roseovarius (80/7.8), Roseovarius (91/9.4), Roseovarius (87/5.2), Roseovarius (81/19.9), Roseovarius (83/21.4), Spiribacter (87/7.8), Thiohalophilus (79/7.7), VMDJ01 (90/2.7), Vreelandella (92/9.3), DSBN01 (87/14.4), DSBN01 (83/6.0), DSBN01 (85/12.0), DSBN01 (83/3.5), DSBN01 (93/11.4), Chordibacter (82/10.6)

Of these 58 original MAGs, 21 were below MIMAG medium quality in the 2026 re-QC. Of the 37 medium- or high-quality ones, 4 have a re-assembly MAG at >= 95 % ANI that covers < 50 % of the genome (fragment recovered), and 30 have no re-assembly MAG above 80 % ANI. By sample type: water 17, sediment 12, R2_restored sediment 8.

**Why the re-assembly MAGs are more fragmented.** Every original MAG contig is >= 2,500 bp. In the re-assembly, MetaBAT2 and CONCOCT also used 2,500 bp, but MaxBin2 used 1,000 bp (pipeline.py), and DAS Tool kept far more MaxBin2 bins than in 2020. MaxBin2 bins have the lowest N50. MetaBAT2/CONCOCT bins are also below the original N50 (19.1 kb), so the metaSPAdes assemblies themselves are more fragmented than the 2020 ones; the cause is not established (read processing differs: Trimmomatic here vs JGI's pipeline in 2020).

| binner   |   original MAGs |   re-assembly MAGs |   re-assembly median N50 (kb) |   re-assembly median contigs |
|:---------|----------------:|-------------------:|------------------------------:|-----------------------------:|
| concoct  |             155 |                307 |                           7.5 |                        405   |
| maxbin   |              14 |                533 |                           3.9 |                        828   |
| metabat  |             141 |                484 |                           8.5 |                        355.5 |

## 2. Taxa of interest (GTDB R232 on both sides)

| taxon                                      |   original MAGs |   original species clusters |   re-assembly MAGs |   re-assembly species clusters |   new species clusters not in original | original HQ / re-assembly HQ MAGs   |
|:-------------------------------------------|----------------:|----------------------------:|-------------------:|-------------------------------:|---------------------------------------:|:------------------------------------|
| Archaea                                    |              30 |                          30 |                116 |                             47 |                                     24 | 11 / 16                             |
| Halobacteriota (haloarchaea + methanogens) |              24 |                          24 |                 80 |                             30 |                                     12 | 10 / 13                             |
| Methanogens (Methanosarcinales)            |               1 |                           1 |                  5 |                              2 |                                      1 | 1 / 3                               |
| Methanolobus                               |               1 |                           1 |                  4 |                              1 |                                      0 | 1 / 3                               |
| Methanohalophilus                          |               0 |                           0 |                  1 |                              1 |                                      1 | 0 / 0                               |
| Rhodobacteraceae                           |              25 |                          24 |                122 |                             38 |                                     19 | 14 / 23                             |
| Roseovarius                                |              12 |                          11 |                 41 |                             12 |                                      5 | 4 / 2                               |
| Marivita                                   |               3 |                           3 |                 17 |                              3 |                                      0 | 2 / 9                               |
| Albimonas                                  |               1 |                           1 |                  5 |                              1 |                                      0 | 1 / 2                               |
| Halomonadaceae                             |               2 |                           2 |                  3 |                              1 |                                      0 | 0 / 0                               |
| Halomonas                                  |               1 |                           1 |                  3 |                              1 |                                      0 | 0 / 0                               |
| Halanaerobiaceae                           |               3 |                           2 |                  3 |                              2 |                                      2 | 0 / 2                               |
| Halanaerobium                              |               3 |                           2 |                  1 |                              1 |                                      1 | 0 / 1                               |
| Desulfobacterota (sulfate reducers)        |              28 |                          27 |                118 |                             50 |                                     23 | 14 / 31                             |

**Methanogens:** original Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_concoct_out.21.conti (R1_A_D2, 98/2.6, Methanolobus sp020055655); re-assembly metabat2.bin.187_sub.fasta__R1_sediment (R1_sediment, 55/0.7, Methanohalophilus portucalensis); re-assembly concoct.bin.097.fasta__R1_A_D2 (R1_A_D2, 93/2.9, Methanolobus sp020055655); re-assembly maxbin2.Bin.045.fasta__R1_C_D1 (R1_C_D1, 61/2.6, Methanolobus sp020055655); re-assembly concoct.bin.101.fasta__R1_sediment (R1_sediment, 98/2.6, Methanolobus sp020055655); re-assembly concoct.bin.079.fasta__R2_sediment (R2_sediment, 98/3.9, Methanolobus sp020055655)

**Albimonas:** original Salt_Pond_MetaG_R2_restored_C_black_MG_DASTool_bins_metabat. (R2_restored_C_black, 99/2.4, Albimonas sp035268645); re-assembly concoct.bin.061.fasta__R1_C_D1 (R1_C_D1, 67/2.5, Albimonas sp035268645); re-assembly concoct.bin.157_sub.fasta__R1_sediment (R1_sediment, 98/3.7, Albimonas sp035268645); re-assembly concoct.bin.018.fasta__R2_B_D2 (R2_B_D2, 52/1.9, Albimonas sp035268645); re-assembly metabat2.bin.062.fasta__R2_restored_C_black (R2_restored_C_black, 78/5.9, Albimonas sp035268645); re-assembly concoct.bin.058_sub.fasta__R2_sediment (R2_sediment, 95/2.4, Albimonas sp035268645)

**Halomonas:** original Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.26.contigs__ (R1_A_D1, 97/5.2, Halomonas fontilapidosi); re-assembly maxbin2.Bin.045.fasta__R1_A_D2 (R1_A_D2, 79/3.0, Halomonas fontilapidosi); re-assembly metabat2.bin.336.fasta__R1_sediment (R1_sediment, 84/2.8, Halomonas fontilapidosi); re-assembly metabat2.bin.443.fasta__R2_sediment (R2_sediment, 75/2.5, Halomonas fontilapidosi)

## 3. Phyla

| phylum              |   original MAGs |   re-assembly MAGs |   original species |   re-assembly species |
|:--------------------|----------------:|-------------------:|-------------------:|----------------------:|
| Pseudomonadota      |              92 |                406 |                 91 |                   159 |
| Bacteroidota        |              68 |                244 |                 65 |                    95 |
| Desulfobacterota    |              28 |                118 |                 27 |                    50 |
| Halobacteriota      |              24 |                 80 |                 24 |                    30 |
| Patescibacteriota   |               4 |                 52 |                  4 |                    26 |
| Actinomycetota      |              19 |                 52 |                 16 |                    23 |
| Bacillota           |              15 |                 88 |                 13 |                    22 |
| Planctomycetota     |               4 |                 36 |                  4 |                    21 |
| Gemmatimonadota     |               5 |                 22 |                  5 |                    14 |
| Acidobacteriota     |               5 |                 17 |                  5 |                    12 |
| Chloroflexota       |               3 |                 14 |                  3 |                    11 |
| Spirochaetota       |              12 |                 23 |                  9 |                    10 |
| Fidelibacterota     |               4 |                 13 |                  4 |                     9 |
| Verrucomicrobiota   |               3 |                 12 |                  3 |                     9 |
| Deinococcota        |               4 |                 22 |                  4 |                     8 |
| Thermoproteota      |               2 |                  9 |                  2 |                     6 |
| Myxococcota_A       |               0 |                  8 |                  0 |                     5 |
| Nanobdellota        |               0 |                  7 |                  0 |                     4 |
| Krumholzibacteriota |               3 |                 14 |                  3 |                     4 |
| Bacillota_I         |               0 |                  3 |                  0 |                     3 |

## 4. By site (original = dereplicated, so counts are after cross-sample dereplication)

| site            |   original MAGs |   re-assembly MAGs |
|:----------------|----------------:|-------------------:|
| R1              |              84 |                267 |
| R1 co-assembly  |               0 |                155 |
| R2              |              76 |                253 |
| R2 co-assembly  |               0 |                199 |
| R2A             |              67 |                116 |
| R2A co-assembly |               0 |                138 |
| R2_restored     |              37 |                 55 |
| SF2             |              46 |                 85 |
| SF2 co-assembly |               0 |                 56 |

Files: comparison/mags_old_vs_new.tsv (per MAG: set, sample, CheckM, R232 lineage, contigs, N50, species cluster); comparison/skani_edges.tsv.
