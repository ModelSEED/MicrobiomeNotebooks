# Salterns binning test: decision NO_GO

Generated 2026-11-13T14:04:09Z. SF2 assembly deadline 2026-11-14T06:18:00Z.

- only R1_A_D1 passes: NO_GO; the SF2-type decision needs Andrew

## R1_A_D1: complete — PASS

| criterion | metric | NEW | OLD | threshold | pass |
|---|---|---|---|---|---|
| a | MQ+ MAGs (CheckM2 completeness >= 50, contamination < 10) | 62 | 47 | NEW >= 1.2 x OLD = 56.40 | True |
| b | HQ MAGs (CheckM2 >= 90, < 5) | 27 | 11 | NEW >= 1.2 x OLD = 13.20 or NEW >= OLD + 5 = 16 | True |
| c | median N50 of MQ+ MAGs (bp) | 15,416 | 4,656 | NEW >= 1.5 x OLD = 6984.0 | True |
| e | originals (derep310 of R1_A_D1, n=11) recovered by an MQ+ MAG at skani ANI >= 95 and max AF >= 50 | 10 | 8 | NEW >= OLD | True |
| d | median CheckM2 contamination of MQ+ (%) | 3.34 | 3.48 | NEW <= OLD + 1.0 = 4.48 | True |
| d | GUNC pass fraction of MQ+ | 0.87 | 0.79 | NEW >= OLD - 0.05 = 0.737 | True |

CheckM1 (informational): MQ+ NEW 60 vs OLD 45 (a True); HQ NEW 29 vs OLD 10 (b True).
KBase CheckM link (OLD bins): `{"status": "ok", "bins_matched": 60, "bins_in_kbase_table": 60, "kbase_mq": 45, "kbase_hq": 10, "checkm1_here_mq": 45, "checkm2_here_mq": 47, "median_abs_diff_completeness_kbase_vs_checkm1": 0.0026473709239098753, "median_abs_diff_contamination_kbase_vs_checkm1": 0.002351103003277011}`

| set | bins | MQ+ | HQ | Mbp in MQ+ | median completeness | median contamination | median N50 | median contigs | GUNC pass | originals recovered |
|---|---|---|---|---|---|---|---|---|---|---|
| NEW | 71 | 62 | 27 | 173.67 | 93.62 | 3.34 | 15,416 | 269.00 | 0.87 | 10 ['Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.44.contigs__', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.5.contigs__', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.54.contigs__', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs__', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.15.contigs__', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.18.contigs__', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.26.contigs__', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.45.contigs__', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs__', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs__'] |
| OLD | 60 | 47 | 11 | 121.36 | 81.55 | 3.48 | 4,656 | 514 | 0.79 | 8 ['Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.44.contigs__', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.5.contigs__', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.54.contigs__', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.9.contigs__', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.26.contigs__', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.45.contigs__', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.47.contigs__', 'Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs__'] |

## SF2_C_D1: complete — FAIL

| criterion | metric | NEW | OLD | threshold | pass |
|---|---|---|---|---|---|
| a | MQ+ MAGs (CheckM2 completeness >= 50, contamination < 10) | 7 | 4 | NEW >= 1.2 x OLD = 4.80 | True |
| b | HQ MAGs (CheckM2 >= 90, < 5) | 2 | 0 | NEW >= 1.2 x OLD = 0.00 or NEW >= OLD + 5 = 5 | True |
| c | median N50 of MQ+ MAGs (bp) | 9,133 | 4,852 | NEW >= 1.5 x OLD = 7277.2 | True |
| e | originals (derep310 of SF2_C_D1, n=1) recovered by an MQ+ MAG at skani ANI >= 95 and max AF >= 50 | 0 | 0 | NEW >= OLD | True |
| d | median CheckM2 contamination of MQ+ (%) | 4.67 | 1.70 | NEW <= OLD + 1.0 = 2.70 | False |
| d | GUNC pass fraction of MQ+ | 0.86 | 1.00 | NEW >= OLD - 0.05 = 0.950 | False |

CheckM1 (informational): MQ+ NEW 7 vs OLD 4 (a True); HQ NEW 2 vs OLD 0 (b True).
KBase CheckM link (OLD bins): `{"status": "ok", "bins_matched": 6, "bins_in_kbase_table": 6, "kbase_mq": 4, "kbase_hq": 0, "checkm1_here_mq": 4, "checkm2_here_mq": 4, "median_abs_diff_completeness_kbase_vs_checkm1": 0.002646820441174924, "median_abs_diff_contamination_kbase_vs_checkm1": 0.0025811243363423886}`

| set | bins | MQ+ | HQ | Mbp in MQ+ | median completeness | median contamination | median N50 | median contigs | GUNC pass | originals recovered |
|---|---|---|---|---|---|---|---|---|---|---|
| NEW | 11 | 7 | 2 | 23.83 | 88.21 | 4.67 | 9,133 | 425 | 0.86 | 0 [] |
| OLD | 6 | 4 | 0 | 11.32 | 64.03 | 1.70 | 4,852 | 583.00 | 1.00 | 0 [] |

## Rule (pre-registered)

```
(fixed before any result was seen; implemented literally):
Per sample S in {R1_A_D1, SF2_C_D1}, compare NEW (Arm C assembly binned like 2026) vs OLD (the 2026 pipeline's DAS Tool
output for the same unit: all its DAS Tool bins), both scored here with the same tools: CheckM2 (primary), CheckM1
lineage_wf --reduced_tree (secondary), GUNC.
 a) MQ+ MAGs (CheckM2 completeness >= 50, contamination < 10): NEW >= 1.2 x OLD
 b) HQ MAGs (>= 90, < 5): NEW >= 1.2 x OLD or NEW >= OLD + 5
 c) median N50 of MQ+ MAGs: NEW >= 1.5 x OLD
 d) median CheckM2 contamination of MQ+ MAGs: NEW <= OLD + 1.0 point; and GUNC pass fraction of MQ+ MAGs:
    NEW >= OLD - 0.05
 e) original species recovered (skani ANI >= 95 and max AF >= 50 between an MQ+ MAG and an original MAG of sample S;
    originals = derep310 MAGs whose sample is S): NEW >= OLD
 S passes if a, b, c, d and e all pass. CheckM1-based a/b are reported (informational).
 GO for the full salterns rollout only if BOTH samples pass. If only R1_A_D1 passes: NO_GO (the SF2-type decision
 needs Andrew). If SF2_C_D1's Arm C assembly fails or never finishes within 24 h: decide on R1_A_D1 alone and mark
 the decision R1_ONLY_PASS (if R1_A_D1 passes) rather than GO.
```

## Conventions

```
(not part of the rule text, stated in DECISION.json):
 - A median or fraction over an empty MQ+ set is undefined. Undefined OLD with defined NEW -> the criterion passes;
   undefined NEW -> it fails.
 - GUNC pass = pass.GUNC == True in GUNC.<db>.maxCSS_level.tsv; a genome missing from GUNC's output counts as not
   passing.
 - A genome missing from the CheckM2 report is not MQ+.
 - A sample whose processing failed (a required step failed after retries) is INCOMPLETE; the overall decision is then
   INCOMPLETE (needs Andrew), except that a failed/timed-out SF2 assembly gives R1_ONLY_PASS / NO_GO as the rule says.
 - (review, 2026-11-13) Every threshold is evaluated in exact rational arithmetic on the values as printed (CheckM2
   completeness/contamination as decimal strings, GUNC pass fractions as ratios, N50 medians), so a value exactly at a
   threshold gets the inclusive inequality the rule states (floats made e.g. 1.36 <= 0.36 + 1.0 and 0.15 >= 0.2 - 0.05
   false). Reported numbers are floats of the exact values.
 - (review, 2026-11-13) GUNC 1.1.1 names a genome basename.split(file_suffix)[0], i.e. it cuts the file name at its
   first '.fa'; GUNC rows are joined to genomes on that key (it changes only originals named '*.fa_assembly'). Two
   genomes of one set with the same key would get no GUNC value (not passing) and are listed in qc_status.
```
