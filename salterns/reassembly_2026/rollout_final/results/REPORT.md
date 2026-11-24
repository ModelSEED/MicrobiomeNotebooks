# Salterns Arm C rollout: summary

Generated 2026-11-19T12:26:57Z by `rolib.py report` from `/scratch1/afreiburger/salterns_2026/rollout`.

KEEP rule: CheckM2 completeness >= 50 AND CheckM2 contamination < 10 AND GUNC pass (`pass.GUNC == True`; a missing GUNC row does not pass). HQ = kept with CheckM2 >= 90 / < 5. CheckM1 `--reduced_tree` values are informational only.

Per-sample assemblies only. The 2026 site co-assemblies (R1/R2/R2A/SF2 sediment) are NOT redone; their 2026 MAGs remain as they are.

- Samples: 39; complete: 39; failed: 0; running: 0.
- NEW DAS Tool bins: 1,394; kept: 935 (HQ 473, MQ 462); dropped: 459.
- Kept MAGs with a GTDB-Tk (R232) classification: 935.
- Assembly mode: default 39.

- Assembled, binned and QC'd on frankensteinsmonster:~/Documents/salterns_rollout_local/rollout with the same driver and recipe, then imported here (md5-verified; the KEEP rule re-applied here gave the same kept set): R2_B_D1, R2_C_D1, R2_C_D2, R2_restored_C_black, R2A_A_D1, R2A_A_D2, R2A_B_D1, R2A_B_D2, R2A_C_D1, R2A_C_D2, R2_A_H2O, R2_B_H2O, R2_C_H2O, R2_restored_H2O, R2A_B_H2O.

- 2026 DAS Tool bins of the same units scored here with CheckM2 + GUNC: 1,074 bins, 603 pass the same KEEP rule (HQ 145).

## Per sample

| sample | status | mode | asm Mbp >=2.5 kb | asm N50 >=1 kb | NEW bins | kept | HQ | kept median N50 | 2026 KBase MAGs (HQ) | 2026 bins CheckM2-kept (HQ) | 2026 kept median N50 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| R2_restored_D_shore | complete | default | 187.6 | 3,987 | 38 | 26 | 12 | 26,367.5 | 19 (0) | 17 (5) | 6,365 |
| SF2_A_H2O | complete | default | 120.1 | 2,737 | 30 | 11 | 6 | 28,240 | 13 (0) | 4 (0) | 3,132.5 |
| R1_A_D1 | complete | default | 250.5 | 4,303 | 71 | 54 | 26 | 16,456.5 | 45 (0) | 37 (10) | 4,911 |
| R1_A_D2 | complete | default | 237.7 | 4,876 | 64 | 49 | 32 | 20,914 | 45 (0) | 37 (10) | 6,033 |
| R1_B_D1 | complete | default | 245.0 | 3,399 | 64 | 51 | 34 | 25,861 | 39 (0) | 31 (9) | 5,526 |
| R1_B_D2 | complete | default | 245.3 | 4,865 | 67 | 53 | 31 | 21,711 | 46 (0) | 40 (12) | 5,733.0 |
| R1_C_D1 | complete | default | 183.6 | 3,993 | 46 | 32 | 21 | 19,146.0 | 28 (0) | 22 (5) | 4,384.5 |
| R1_C_D2 | complete | default | 182.0 | 3,900 | 47 | 35 | 20 | 22,844 | 29 (0) | 20 (5) | 5,252.5 |
| SF2_A_D1 | complete | default | 49.6 | 2,250 | 10 | 4 | 3 | 20,529.5 | 3 (0) | 3 (1) | 5,293 |
| SF2_A_D2 | complete | default | 95.8 | 2,566 | 23 | 12 | 8 | 41,348.0 | 12 (0) | 9 (4) | 6,404 |
| SF2_B_D1 | complete | default | 69.4 | 2,067 | 16 | 7 | 5 | 13,493 | 7 (0) | 7 (0) | 5,176 |
| SF2_B_D2 | complete | default | 67.1 | 2,116 | 17 | 11 | 5 | 10,675 | 7 (0) | 6 (2) | 8,298.5 |
| SF2_C_D1 | complete | default | 43.9 | 2,164 | 11 | 6 | 2 | 11,331.5 | 4 (0) | 4 (0) | 4,851.5 |
| SF2_C_D2 | complete | default | 107.5 | 3,577 | 27 | 22 | 12 | 11,238.0 | 16 (0) | 16 (5) | 4,777.5 |
| R1_A_H2O | complete | default | 60.2 | 1,796 | 18 | 6 | 1 | 42,742.5 | 9 (0) | 5 (1) | 26,207 |
| R1_B_H2O | complete | default | 94.7 | 5,506 | 23 | 14 | 4 | 24,893.0 | 13 (0) | 7 (1) | 6,129 |
| R1_C_H2O | complete | default | 102.0 | 3,258 | 33 | 13 | 6 | 17,013 | 13 (0) | 8 (1) | 6,135.5 |
| SF2_B_H2O | complete | default | 156.4 | 3,575 | 33 | 16 | 3 | 22,935.5 | 13 (0) | 7 (0) | 4,352 |
| SF2_C_H2O | complete | default | 200.6 | 4,066 | 36 | 16 | 5 | 28,423.5 | 10 (0) | 6 (0) | 4,138.5 |
| R2_A_D1 | complete | default | 290.1 | 3,828 | 75 | 57 | 28 | 16,095 | 51 (0) | 43 (8) | 5,631 |
| R2_A_D2 | complete | default | 139.1 | 3,024 | 31 | 17 | 6 | 25,417 | 16 (0) | 12 (0) | 5,508.0 |
| R2_B_D1 | complete | default | 164.1 | 2,823 | 48 | 31 | 14 | 11,345 | 26 (0) | 16 (3) | 3,519.0 |
| R2_B_D2 | complete | default | 355.3 | 3,531 | 90 | 71 | 35 | 15,678 | 54 (0) | 43 (15) | 5,738 |
| R2_C_D1 | complete | default | 240.8 | 3,915 | 60 | 48 | 23 | 13,569.0 | 39 (0) | 31 (7) | 4,690 |
| R2_C_D2 | complete | default | 151.8 | 3,202 | 42 | 30 | 17 | 11,106.0 | 25 (0) | 22 (5) | 5,954.5 |
| R2_restored_C_black | complete | default | 226.4 | 3,803 | 53 | 41 | 23 | 18,640 | 32 (0) | 27 (6) | 5,887 |
| R2A_A_D1 | complete | default | 98.8 | 2,173 | 25 | 15 | 7 | 8,555 | 15 (0) | 13 (2) | 4,340 |
| R2A_A_D2 | complete | default | 77.5 | 1,868 | 20 | 18 | 9 | 8,010.0 | 9 (0) | 9 (3) | 6,702 |
| R2A_B_D1 | complete | default | 29.0 | 1,633 | 10 | 5 | 1 | 5,285 | 2 (0) | 1 (0) | 3,965 |
| R2A_B_D2 | complete | default | 67.3 | 2,091 | 19 | 15 | 3 | 5,617 | 12 (0) | 11 (2) | 5,047 |
| R2A_C_D1 | complete | default | 94.8 | 2,373 | 26 | 19 | 8 | 9,379 | 16 (0) | 13 (3) | 6,614 |
| R2A_C_D2 | complete | default | 211.7 | 2,632 | 52 | 34 | 19 | 24,902.5 | 31 (0) | 28 (12) | 7,795.5 |
| R2_A_H2O | complete | default | 156.7 | 4,946 | 27 | 21 | 10 | 25,748 | 14 (0) | 9 (2) | 6,904 |
| R2_B_H2O | complete | default | 130.0 | 4,922 | 25 | 16 | 7 | 25,094.5 | 15 (0) | 9 (2) | 6,815 |
| R2_C_H2O | complete | default | 145.9 | 4,773 | 25 | 18 | 6 | 28,888.0 | 13 (0) | 6 (2) | 8,446.5 |
| R2_restored_H2O | complete | default | 75.5 | 5,796 | 7 | 4 | 0 | 10,834.0 | 4 (0) | 5 (0) | 4,510 |
| R2A_A_H2O | complete | default | 130.2 | 3,958 | 34 | 13 | 7 | 29,315 | 10 (0) | 9 (1) | 3,327 |
| R2A_B_H2O | complete | default | 121.5 | 4,897 | 28 | 11 | 7 | 71,771 | 10 (0) | 4 (0) | 2,444.0 |
| R2A_C_H2O | complete | default | 85.1 | 2,811 | 23 | 13 | 7 | 27,633 | 11 (0) | 6 (1) | 2,507.5 |

Columns: '2026 KBase MAGs' = the KBase-CheckM-selected 2026 MAGs of the same unit (`comparison/mags_old_vs_new.tsv`, KBase CheckM 1.1.2). '2026 bins CheckM2-kept' = all 2026 DAS Tool bins of the unit re-scored here with CheckM2 + GUNC under the same KEEP rule (CheckM2 on both sides).

Files: `samples_summary.tsv`, `bins_all.tsv`, `mags_kept.tsv`, `bins_dropped.tsv`, `old2026_bins_checkm2_gunc.tsv`; kept FASTAs in `../mags_kept/`; GTDB-Tk output in `../gtdbtk/`.
