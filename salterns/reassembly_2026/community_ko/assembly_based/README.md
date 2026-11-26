# Assembly-based community KO profile of the 39 salterns metagenomes (panel + marker KOs, two-pass KOfam)

*Track `community_ko/assembly_based`, run on poplar on 2026-11-25, 14:56–20:18 CST (all times US Central, CST; poplar logs are in UTC and were converted). Poplar work dir: `/scratch1/afreiburger/salterns_2026/community_ko/`; this folder holds the scripts, the small results and the checks. Nothing was committed and nothing was written to KBase. No other track's files were modified (read_markers, kofam_rollout, kofam, kofam_coasm38, rollout, binning_test were only read; the read-based survey `read_markers_v2/` was not touched), and no other user's or session's process was touched.*

## Status

- **Full search (stopped).** 2026-11-25 15:16:25 CST: all 27,981 KOfam profiles × all 9,764,756 unique proteins in 521 chunks (`run_all.sh`, 90 single-thread jobs). poplar's load rose to 750–1,150 on 336 cores soon after the launch, so each chunk took ≈ 3.4–4 h (projected ≈ 20 more hours). At 18:53 CST, on the coordinator's request, the run was stopped (my driver, its `xargs` queue, the 90 chunk shells and their `hmmsearch` processes, killed by PID from the driver's process tree; `logs/kill_tree.txt`). 4 chunks (0027, 0046, 0049, 0063) had finished all profiles; they and the 86 partial tblouts are kept in `full_run_stopped/` on poplar and serve as the equivalence check.
- **Two-pass design (used).** Pass 1 (screen, 18:55:23–19:00:22 CST, 5 min): the 148 profiles of the 138 panel KOs + 10 marker KOs on all 521 chunks. Candidate selection (19:00–19:01:35 CST): 105,559 representatives (35.5 M residues; 142,487 genes). Pass 2 (confirm, 19:01:43–19:51:45 CST, 50 min): all 27,981 profiles on the candidates, 180 chunks, 90 single-thread jobs, 0 failed. Aggregation 19:51:45–19:57:27 CST, re-run 20:04:42–20:10:22 CST with the v1 reading of the marker rescue for the GE (see Genome equivalents). Checks (`04`, `07`, `09`) 19:53–20:03 CST; local tables 20:12 CST.
- **DONE** (`ALL_DONE` on poplar, 20:10:22 CST). The equivalence check against the 4 finished full-run chunks and the gene-level check against the kofam_rollout MAG track both show 0 differences (Checks).
- **Scope change:** `results/ko_by_sample.tsv.gz` covers the **138 panel KOs + 10 marker KOs** (148 KOs × 39 samples), **not every KEGG ortholog** as first planned; the full-DB calls exist only for the 105,559 pass-1 candidates. A complete all-KO profile needs the full search of all 9.76 M unique proteins (≈ 1,900 thread-hours at poplar's load of 2026-11-25). `scripts/run_all.sh` can be re-launched for that; first move the 4 finished chunks' files from `full_run_stopped/hmm/` and `full_run_stopped/results_chunks/` back to `hmm/` and `results/chunks/` so they are skipped, and note that `03_aggregate.py` now reads the pass-2 calls (`P2DIR`).

## Methods

### Inputs

| input | where (poplar unless noted) |
|---|---|
| Gene calls | `read_markers/genes/<s>/<s>.faa`: Prodigal 2.6.3 `-p meta -g 11 -q` on all contigs ≥ 1 kb of the 39 per-sample rollout assemblies (`read_markers/scripts/01_prodigal.sh`), 12,055,118 proteins, IDs `<sample>|<contig>_<n>` |
| KOfam | `kofam/db/kofam_all.hmm` (KOfam 2026-10-26, 27,981 profiles, md5 `6fcb36e4…`; linked as `db/kofam_all.hmm`, not rebuilt), `db/hmm_info.tsv` (copy of `kofam_rollout/db/hmm_info.tsv`, md5 `88c4b270…`: threshold and score type per KO from `ko_list` 2026-10-26) |
| Contig depth | own-library mean depth per contig, `rollout/<s>/depth/<s>.d1000.txt` (`jgi_summarize_bam_contig_depths`, column `totalAvgDepth`; every file has one BAM, i.e. 5 columns); `binning_test/<s>/depth/` for R1_A_D1 and SF2_C_D1 |
| Contig status | `read_markers/meta/contig_status.tsv.gz`: kept MAG (935 rollout MAGs) / other DAS Tool bin / unbinned, per contig ≥ 1 kb |
| Domain of binned contigs | kept MAG: GTDB R232 domain (`rollout/results/mags_kept.tsv`); other DAS Tool bin: CheckM1 marker-set lineage (`rollout/<s>/keep/bins.tsv`) mapped by name (archaeal taxa → Archaea, `root` → unknown, any other → Bacteria; agrees with GTDB for every kept bin, see Checks) |
| Read-based GE depth | `read_markers/results/genome_equivalents.tsv` (`GE_depth_mean`, 8 ribosomal-protein families) |
| Panel | `../ko_panel_methylated.tsv` (138 KOs, incl. the 30 FigS4 methanogenesis-pathway KOs) |
| MAG KO calls | `annotation_runs/{kofam_rollout,kofam,kofam_coasm38}/results/all_mag_ko_counts.tsv`, `all_no_threshold_hits.tsv.gz` (local) |

### Steps (scripts in `scripts/`, identical copies on poplar; md5 in `poplar_copies/scripts_md5_*.txt`)

| step | script | what | where / threads |
|---|---|---|---|
| 0 | `00_smoke_compare.py` | smoke test before the full launch (see Checks) | poplar, 1 thread |
| 1 | `01_make_chunks.py` | global de-duplication: one representative per distinct protein sequence (md5 of the Prodigal protein incl. `*`), first occurrence in sorted sample order; `meta/genes/<s>.tsv.gz` (every gene: contig, partial, aa length, representative); 521 chunks of ≈ 4.5 M residues. 12,055,118 proteins → 9,764,756 unique (19.0 % duplicates), 2.90 → 2.35 G residues | poplar, 1 thread, 3.4 min |
| (2) | `run_all.sh` → `run_chunk.sh` → `02_postprocess_chunk.py` | the stopped full search (all profiles × all chunks); 4 chunks finished | poplar, ≤ 90 threads |
| P1 | `run_pass1.sh` → `p1_run_chunk.sh` | pass 1: HMMER 3.4 `hmmsearch --cpu 0 --noali -T 14.7 --domT 14.7 -Z 12055118 --tblout … --domtblout … -o /dev/null db/pass1_148.hmm chunk.faa` on all 521 chunks. `db/pass1_148.hmm` = the 138 panel KOs + 10 marker KOs, concatenated from `kofam_db/profiles/<KO>.hmm`; md5 `d92b324a…`, identical to the same 148 records cut from `kofam_all.hmm` (`logs/pass1_hmm_check.txt`) | poplar, `xargs -P 90`, 1 thread per job, 5 min |
| sel | `08_select_candidates.py` | candidates for pass 2 (any pass-1 hit with score ≥ 0.5 × the KO's threshold — includes every passing hit and the marker-rescue range — or a no-threshold panel KO hit with full E ≤ 1e-30); 180 pass-2 chunks | poplar, 1 thread, 1 min |
| P2 | `run_pass2.sh` → `p2_run_chunk.sh` → `02_postprocess_chunk.py` | pass 2: the full-search command (`db/kofam_all.hmm`, same flags, same `-Z 12055118`, tblout only) on the candidate chunks; per-chunk KOfam call rule (below) | poplar, `xargs -P 90`, 1 thread per job |
| 3 | `03_aggregate.py` | representative calls → every gene; × contig depth, status, domain; per-sample profile of the 148 KOs; GE | poplar, 1 thread |
| 4 | `04_validate_vs_rollout.py` | gene-level check against the kofam_rollout MAG track | poplar, 1 thread |
| 9 | `09_two_pass_check.py` | equivalence of the two-pass design with the full search on the 4 finished full-run chunks | poplar, 1 thread |
| 7 | `07_nothreshold_panel.py` | supplement for the five panel KOs without a KOfam threshold | poplar, 1 thread |
| 5 | `05_panel_tables.py` | `panel_status_split.tsv`, `panel_mag_attribution*.tsv`, `nothreshold_panel_*.tsv` | local |
| 6 | `06_validate_vs_reads.py` | comparison with the v1 read survey | local |

Every poplar job ran under `nice -n 10`, detached (`nohup setsid`), with `TMPDIR` on /scratch1 (`community_ko/tmp`). `hmmsearch --cpu 0` runs serially in one thread (no reader thread), so the track used ≤ 90 busy threads plus idle shells, under the 96-thread cap.

**Flags vs the kofam_rollout track** (`annotation_runs/kofam_rollout/scripts/run_chunk.sh`): identical `--noali -T 14.7 --domT 14.7`; `-Z` = all genes of this set (12,055,118) instead of 2,732,576; `--cpu 0` instead of `--cpu 1` (same results, see Checks); no `--domtblout` in the full-DB search (it served only the HMM-coverage column there; the call uses the tblout's full score and best-single-domain score); `-o /dev/null`.

**Why the two-pass design is exact for the 148 KOs.** A gene's score and E-value against a profile do not depend on which other profiles are searched (same flags and `-Z`). A gene can only get a panel or marker KO as its best KO if it passes that KO (score > threshold), and only be rescued as a marker if its marker score is ≥ 0.5 × threshold; both are visible in pass 1, so every such gene is a candidate, and pass 2 gives it the same best KO, passing set and top-scoring profile as the full search. A non-candidate has no panel or marker KO in either design. For any other KO the two-pass counts are partial (only candidates were searched against it), so the KO outputs hold only the 148 KOs. Checked on the 4 finished full-run chunks (Checks).

### KOfam call rule (= `kofam_rollout/scripts/postprocess.py`, DRAM 1.5 `sig_scores`)

- A (gene, KO) hit passes when score **>** the KO's `ko_list` threshold; score = full-sequence score for `score_type == full`, best single-domain score (tblout column 9) for `domain`.
- KOs with threshold `-` (1,347 of 27,981) are never called.
- Best KO per gene = the passing hit with the lowest full-sequence E-value; ties → higher full score → KO id. Only the best KO is counted (one KO per gene); all passing KOs are kept in the per-gene table.
- Identical proteins get the representative's call (exact: hmmsearch scores do not depend on the target set, E-values only through the fixed `-Z`).

### Genome equivalents (GE)

- **Assembly-internal GE depth** = mean over the 8 ribosomal-protein marker families of the read_markers survey (uL2 K02886, uL11 K02867, uS3 K02982, uS9 K02996, uS11 K02948, uS12 K02950, uL4 K02926+K02930, uL16 K02878+K02866) of each family's summed assembled depth (Σ over the family's genes of their contig's own-library mean depth).
- Marker genes = best-KO calls of the 10 KOs **plus the v1 rescue**: a gene with no passing KO (over all 27,981 profiles) whose top-scoring profile among the 62 profiles of the v1 read_markers panel (`meta/v1_panel_kos.txt`) is one of the 10 marker KOs and whose score (by that KO's score type) is ≥ 0.5 × its threshold (KOfam under-calls archaeal ribosomal proteins; `read_markers/scripts/04_build_db.py`, `marker_relaxed_0.5thr`). The rescue is used only for the GE denominator; `ko_by_sample` counts strict best-KO calls.
- Why the v1 reading: the first aggregation took the top-scoring profile over **all** 27,981 profiles. That loses archaeal uS9, uS11 and uL2 genes whose top profile is the eukaryote-named archaeal homolog KO (K02960 S16e, K02955 S14e, K02938 L8e; scores below those KOs' own thresholds), e.g. a uS9 and a uS11 gene on archaeal other-bin contigs at depth ≈ 465 in R1_A_H2O; the marker CV of the R1 waters was then 30–39 % (`poplar_copies/ge_assembly.run1_allprofile_rescue.tsv`). With the v1 reading it is 15–20 %. Both readings are in `results/marker_rescue_check.tsv` and `results/ge_assembly.tsv` (`*_allprofile`).
- **Read-based GE depth** = `GE_depth_mean` of `read_markers/results/genome_equivalents.tsv`.

### Per-sample KO profile (`results/ko_by_sample.tsv.gz`: the 138 panel KOs + 10 marker KOs × 39 samples, zeros filled; **not every KEGG KO**, see the two-pass design)

| column | definition |
|---|---|
| `gene`, `role`, `kofam_threshold`, `definition` | panel gene name (marker family for markers), panel / marker, `ko_list` threshold (`-` = never called), KO definition |
| `n_genes` | genes (Prodigal proteins on contigs ≥ 1 kb) whose best KO (over all 27,981 profiles) is the KO |
| `sum_depth` | Σ own-library mean depth of their contigs (≈ summed coverage of the assembled copies) |
| `copies_per_GE_assembly` | `sum_depth` / assembly-internal GE depth (≈ copies per assembled genome) |
| `copies_per_GE_readGE` | `sum_depth` / read-based GE depth (≈ assembled copies per genome in the community; comparable to the read survey's copies per GE) |
| `gene_fraction` | `n_genes` / all Prodigal genes of the sample (IMG-like unweighted profile) |

`results/ko_status_domain_by_sample.tsv.gz`: the same `n_genes` / `sum_depth` split by contig status (`kept_mag`, `other_bin`, `unbinned`) and, for binned contigs, domain (`kept_mag_Archaea`, `kept_mag_Bacteria`, `other_bin_Archaea`, `other_bin_Bacteria`, `other_bin_unknown`).

### Panel tables (local)

- `results/panel_status_split.tsv`: per panel KO × sample: `n_genes`, `sum_depth`, `copies_per_GE_readGE`, depth and depth-weighted share on kept-MAG / other-bin / unbinned contigs, the binned depth by domain, archaeal / bacterial shares of the binned depth, and archaeal / bacterial binned copies per read-based GE. Unbinned contigs carry no domain.
- `results/panel_mag_attribution.tsv`: per panel KO × final-585 species that carries it (best-KO count > 0 in the representative MAG, or in any MAG of the species; members from `comparison_rollout/mags_three_way.tsv`, sets rollout + coassembly = 1,455 MAGs). Calls: kofam_rollout (935 rollout MAGs), first-pass kofam (482 co-assembly MAGs; names `concoct.bin.002.fasta__R1_sediment` → `R1_sediment__concoct.bin.002`), kofam_coasm38 (the other 38). GTDB R232 taxonomy of the representative; `mean_read_fraction_12_unrestored_core_samples` = simple mean of the 12 samples R1/R2 × A–C × D1/D2 in `abundance_final/final585_relabund_24sediment.tsv` (fraction of QC reads). The `mean_rel_abund_12cores` column of `guilds_final585.csv` is a different quantity (all 12 cores of the four sites, D1/D2 merged DNA-weighted) and is carried along for reference only. `results/panel_mag_attribution_summary.tsv`: per KO, number of carrier species and their summed read fraction.
- **No-threshold panel KOs** (K16967 dmoA, K23996 mpnH/phnY\*, K05915 phnZ, K18049 mpnS, K12905 hepD): never called by the KOfam rule, so they are 0 in the main tables. Supplement (`07_nothreshold_panel.py`, `05_panel_tables.py`): the v1 mpnS/hepD rule extended to the five — a gene with no passing KO gets the KO if that hit has full E ≤ 1e-30 and HMM coverage ≥ 0.5 (coverage from a targeted hmmsearch of the five profiles with `--domtblout`), reported with and without v1's extra condition that the KO is the gene's top-scoring profile → `results/nothreshold_panel_by_sample.tsv`; the MAG side from the MAG tracks' `all_no_threshold_hits.tsv.gz` with the same E / coverage / no-passing-KO rule (top-scoring condition not checkable there) → `results/nothreshold_panel_mag_attribution.tsv`.

## Re-run (exact commands)

```bash
W=/scratch1/afreiburger/salterns_2026/community_ko; PYP=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin/python3
# poplar: chunks (once), pass 1 (+ candidate selection), pass 2 (+ 03_aggregate -> ALL_DONE); each driver is restart-safe
ssh -o BatchMode=yes poplar "cd $W && nice -n 10 $PYP -I scripts/01_make_chunks.py $W /scratch1/afreiburger/salterns_2026/read_markers/genes 4500000"
ssh -o BatchMode=yes poplar "cd $W && nohup setsid bash scripts/run_pass1.sh > logs/run_pass1.nohup 2>&1 < /dev/null &"   # -> pass1/SELECT_DONE
ssh -o BatchMode=yes poplar "cd $W && nohup setsid bash scripts/run_pass2.sh > logs/run_pass2.nohup 2>&1 < /dev/null &"   # -> ALL_DONE
# poplar: checks and supplement (1 thread each)
ssh -o BatchMode=yes poplar "cd $W && for s in 04_validate_vs_rollout 09_two_pass_check 07_nothreshold_panel; do nice -n 10 $PYP -I scripts/\$s.py $W; done"
# local (from this folder)
rsync -a -e "ssh -o BatchMode=yes" poplar:$W/results/{ko_by_sample.tsv.gz,ko_status_domain_by_sample.tsv.gz,ge_assembly.tsv,marker_rescue_check.tsv,checkm1_lineage_domain_map.tsv,checkm1_lineage_domain_other_bins.tsv,aggregate_qc.tsv,aggregate_qc.txt,validation_vs_kofam_rollout.txt,validation_vs_kofam_rollout_diffs.tsv,nothreshold_panel_by_sample_long.tsv,two_pass_equivalence_check.txt,two_pass_equivalence_diffs.tsv} results/
~/Documents/py_venv/bin/python scripts/05_panel_tables.py
~/Documents/py_venv/bin/python scripts/06_validate_vs_reads.py
```

## Results

**Scale.** 12,055,118 genes (9,764,756 distinct proteins). Pass 1: 1,285,809 hits of the 148 profiles. Candidates: 105,559 distinct proteins = 142,487 genes (reasons, overlapping: passes a panel KO 26,555, passes a marker KO 12,681, panel KO at 0.5–1 × threshold 69,334, marker KO at 0.5–1 × threshold 1,686, no-threshold panel KO at E ≤ 1e-30 167). Of the 142,487 searched genes, 100,177 have a passing KO (any KO) and **56,324 genes have a panel or marker KO as best KO** (0.47 % of all genes; 21,398 of them are marker genes), plus 1,581 marker-rescue genes (v1 reading). 131 of the 148 KOs are found in ≥ 1 sample; never found: AOX, cntB, ddd (K28072), ddhA, ddhC, fwdH, mdh, mtsA, mxaF, mxaI, torA, torZ, and the five no-threshold KOs (dmoA, mpnH/phnY\*, phnZ, mpnS, hepD; see the supplement).

**Genome equivalents** (`results/ge_assembly.tsv`). Assembly-internal GE depth / read-based GE depth (≈ the assembled fraction of the community): 0.38–0.62 in the 12 unrestored cores, 0.11–0.44 in the SF2 / R2A sediment cores, 0.39–0.81 in the waters, 0.51 in R2_restored C_black / D_shore. Coefficient of variation of the 8 marker-family depths per sample: **3.2–21.0 %, median 6.9 %** (unrestored cores 3.2–21.0 %, median 6.1 %; other sediment cores 4.4–14.8 %; waters 3.9–20.9 %, median 14.1 %), vs 3.8–13.3 % (median 6.8 %) for the read-based GE. The highest CVs: R2_A_D2 21.0 % (uL11 at 0.67 × the family mean), R2A_A_H2O 20.9 %, R1_A_H2O 19.9 %. Use `copies_per_GE_readGE` for comparisons with the read survey and across samples; `copies_per_GE_assembly` is the copies per assembled genome.

**Preview, `results/ko_by_sample.tsv.gz`** (3 of 39 samples):

| sample | KO | gene | n_genes | sum_depth | copies_per_GE_assembly | copies_per_GE_readGE | gene_fraction |
|---|---|---|---|---|---|---|---|
| R1_A_D1 | K00399 | mcrA | 1 | 9.48 | 0.0063 | 0.0039 | 2.2e-06 |
| R1_A_D1 | K06163 | phnJ | 17 | 142.8 | 0.095 | 0.059 | 3.7e-05 |
| R1_A_D1 | K21578 | grdI | 9 | 152.2 | 0.102 | 0.063 | 1.9e-05 |
| R1_A_D1 | K16176 | mtmB | 7 | 60.4 | 0.040 | 0.025 | 1.5e-05 |
| R1_A_D1 | K00200 | fwdA | 6 | 29.7 | 0.020 | 0.012 | 1.3e-05 |
| R2_B_D1 | K06163 | phnJ | 13 | 79.7 | 0.116 | 0.044 | 3.5e-05 |
| R2_B_D1 | K21578 | grdI | 5 | 27.1 | 0.040 | 0.015 | 1.3e-05 |
| SF2_A_H2O | K06163 | phnJ | 1 | 4.5 | 0.0038 | 0.0017 | 3.2e-06 |

**Panel, 12 unrestored-core samples pooled** (`results/panel_status_split.tsv`; copies per read-based GE = Σ depth / Σ read GE depth; shares are depth-weighted):

| gene | KO | genes | assembled copies per GE (%) | kept MAG / other bin / unbinned | archaeal share of the binned depth | archaeal binned copies per GE (%) |
|---|---|---|---|---|---|---|
| mcrA | K00399 | 4 | 0.10 | 0.72 / 0.13 / 0.15 | 0.84 | 0.07 |
| mtbB | K16178 | 10 | 0.41 | 0.76 / 0.20 / 0.04 | 0.21 | 0.08 |
| mtmB | K16176 | 34 | 0.79 | 0.35 / 0.09 / 0.55 | 0.79 | 0.28 |
| mttB | K14083 | 1,189 | 39.1 | 0.55 / 0.11 / 0.33 | 0.00 | 0.08 |
| phnJ | K06163 | 204 | 7.0 | 0.49 / 0.16 / 0.34 | 0.08 | 0.38 |
| grdI | K21578 | 78 | 3.3 | 0.66 / 0.15 / 0.19 | 0 | 0 |
| cntA | K22443 | 77 | 2.6 | 0.53 / 0.17 / 0.30 | 0 | 0 |
| fwdA | K00200 | 49 | 0.91 | 0.25 / 0.20 / 0.55 | 0.18 | 0.07 |
| ftr | K00672 | 47 | 0.91 | 0.35 / 0.17 / 0.49 | 0.26 | 0.12 |
| mtrA | K00577 | 10 | 0.27 | 0.57 / 0.15 / 0.28 | 1.00 | 0.20 |
| cdhA | K00192 | 9 | 0.17 | 0.42 / 0.08 / 0.50 | 0.84 | 0.07 |
| mch | K01499 | 101 | 3.8 | 0.47 / 0.27 / 0.26 | 0.68 | 1.9 |
| mer | K00320 | 189 | 6.2 | 0.50 / 0.21 / 0.30 | 0.65 | 2.9 |
| ackA | K00925 | 1,077 | 36.0 | 0.57 / 0.14 / 0.29 | 0.01 | 0.21 |

- The FigS4 methanogenesis-pathway KOs are mostly not methanogen genes here: fwdA and ftr are 74–82 % bacterial in the binned part (their MAG carriers are mainly Planctomycetota, Pseudomonadota and Bacillota species), and **archaeal is not methanogen**: among the species representatives the archaeal carriers of mch and mer are Halobacteria (26 of 27 and 34 of 35 archaeal species; the other one is *Methanolobus*, class Methanosarcinia), and mtrA is also carried by 4 Nitrososphaeria species and 1 E2 species. For methanogen abundance use mcrA, or the species of `panel_mag_attribution.tsv`. The unbinned part (15–55 % of the depth) has no domain.
- `results/panel_mag_attribution.tsv`: all 1,455 MAGs of the 585 species have KOfam calls (935 kofam_rollout, 482 first-pass kofam, 38 kofam_coasm38). 121 of the 138 panel KOs are carried by ≥ 1 species representative. Examples: mcrA 1 species (*Methanolobus*); phnJ 35 species (rep) / 39 (any MAG), led by *Roseovarius*, *Marivita*; grdI 15 / 17.

**Supplement: no-threshold panel KOs** (`results/nothreshold_panel_by_sample.tsv`, `nothreshold_panel_mag_attribution.tsv`). Of 2,287 hits of the five profiles, 168 have E ≤ 1e-30, 98 of those genes pass no KO, and 88 distinct proteins have HMM coverage ≥ 0.5: dmoA (K16967) 47, phnZ (K05915) 40, mpnS (K18049) 1; only 1 (the mpnS gene) is also its gene's top-scoring profile (v1's extra condition). Expanded to genes: dmoA 48, phnZ 46, mpnS 1 gene over the 39 samples. MAG side: only dmoA, in 10 MAGs of 5 Pseudomonadota species (*JACKII01*, *CAJQQK01*, *Roseovarius*, two unnamed genera). No hepD or mpnH/phnY\* gene passes the rule.

## Checks

1. **Smoke test before the launch** (`validation/smoke/smoke_compare.txt`): the first 1,000 proteins of kofam_rollout chunk 01 with this track's command (`--cpu 0`, tblout only) at kofam_rollout's `-Z`: all 22,433 tblout rows have identical scores and E-values (1,075 differ only in the `inc` column, the number of domains passing the inclusion threshold, which depends on the chunk's composition through domZ and is not used); the 464 genes with a passing KO have identical best KO, score and passing set. Fixed cost of one full-DB run: 2.5 min (`s1.log`).
2. **Pass-1 HMM** = the same 148 records cut from `kofam_all.hmm` (md5 `d92b324a95cd33d073f1f5e81ed07019` both; `poplar_copies/pass1_hmm_check.txt`).
3. **Two-pass vs full search** (`results/two_pass_equivalence_check.txt`, `09_two_pass_check.py`), on the 4 chunks that finished all profiles in the full run (0027, 0046, 0049, 0063; 87,728 distinct proteins, 730 of them candidates, 282 with a panel/marker best KO): tblout rows of the 148 profiles identical in all 18 fields (8,095 / 8,095); panel/marker label identical for 87,728 / 87,728; marker rescue (both readings) identical for all; complete call rows of the candidates identical in all 12 columns (730 / 730); 0 proteins that pass a panel/marker KO or are rescued in the full run were missed by pass 1. **0 differences.**
4. **Gene level vs the kofam_rollout MAG track** (`results/validation_vs_kofam_rollout.txt`): all 2,732,576 rollout-MAG genes exist here under the same contig and gene number, with the identical protein, all on contigs of the same kept MAG, and no kept-MAG gene here is missing there. Panel/marker label identical for 2,732,576 / 2,732,576 genes (19,281 with a panel/marker KO on both sides); for the 48,392 of them searched in pass 2, best KO, score and the full list of passing KOs with scores are identical for all. 0 rollout genes passing a panel/marker KO were missed by pass 1.
5. **Domain of other bins**: the CheckM1-lineage → domain rule agrees with GTDB R232 for 935 / 935 kept bins; the 459 other bins → 407 Bacteria, 52 Archaea, 0 unknown (`results/checkm1_lineage_domain_*.tsv`).
- Timing note: checks 4 (`04_validate_vs_rollout.py`) and the supplement (`07_nothreshold_panel.py`) ran on the `genes_ko/` files of the first aggregation (19:51–19:57 CST). The re-run (20:04–20:10 CST) changed only the GE denominator and renamed the two rescue columns of `genes_ko/`; every column these two scripts read (`searched_full_db`, `ko`, `score`, `all_passing`, `status`, `mag_id`, `rep`, `depth`, `domain`) is the same in both runs.
6. **Code-path test** (before the production numbers existed): a 380-protein test chunk of R1_A_D1 (120 marker + 60 target genes of v1 + 200 random) through postprocessing and the aggregation: all 181 v1-labelled genes got the same KO or marker-family rescue as v1 (poplar `test_agg/`).

## Validation against the v1 read survey (`results/validation_vs_read_survey*.tsv`)

Assembly side = this track's summed contig depth / read-based GE depth (markers: family depth incl. the v1 rescue); read side = `read_markers/results/per_sample_KO.tsv` copies per GE. Assembled fraction = assembly / read. Labels not compared: grdB, pmoA, mmoX (not panel KOs, so only partially searched here) and mpnS (no KOfam threshold).

- **Against v1's own assembly complement** (`read_markers/results/assembly_depth_complement.tsv`): identical for every target label (pooled ratio 1.000; same genes, same depth); markers 0.9994–1.0003 (≤ 5 genes per family differ over 39 samples: v1 called the best KO among its 62 panel profiles with score ≥ threshold, this track among all 27,981 with score > threshold).
- **Assembled fraction of the read-based copies per GE**, pooled over samples (all 39 / 12 unrestored cores / 12 SF2+R2A sediment cores / 13 waters):

| gene | all 39 | unrestored cores | SF2 + R2A cores | waters | per-sample range (samples with ≥ 20 reads) | Spearman ρ (39 samples) |
|---|---|---|---|---|---|---|
| 8 markers | 0.42–0.50 | 0.47–0.56 | 0.19–0.24 | 0.55–0.69 | 0.10–1.00 | (≈ 1 copy per GE everywhere) |
| mcrA | 0.33 | 0.34 | 0 | 0 | 0–1.01 | 0.57 |
| mcrB / mcrG | 0.24 / 0.29 | 0.31 / 0.37 | 0 | 0 | 0–0.72 | 0.46 / 0.55 |
| phnJ | 0.47 | 0.54 | 0.08 | 0.51 | 0–0.69 | 0.93 |
| phnG / H / I / K / M | 0.38–0.54 | 0.44–0.62 | 0.05–0.09 | 0.40–0.62 | 0–0.90 | 0.80–0.93 |
| phnL | 0.18 | 0.24 | 0.01 | 0.24 | 0–0.41 | 0.79 |
| grdI | 0.47 | 0.49 | 0.13 | 0.21 | 0–0.64 | 0.87 |
| cntA | 0.33 | 0.48 | 0.28 | 0.27 | 0–0.63 | 0.74 |
| cutC | 0.11 | 0.33 | 0.01 | 0 | 0–0.62 | 0.42 |
| mtbB / mtmB | 0.13 / 0.16 | 0.48 / 0.51 | 0.11 / 0.11 | 0 / 0.01 | 0–1.04 | 0.67 / 0.75 |
| mttB | 0.33 | 0.53 | 0.19 | 0.43 | 0–0.95 | 0.66 |
| mtaB / mtaA / mtbA | 0.01–0.03 | 0–0.13 | ≤ 0.02 | 0 | 0–0.37 | 0.18–0.46 |

- Agreement: where a gene is assembled, the assembly-based and read-based abundances track each other across samples (ρ 0.74–0.93 for the C–P lyase genes, grdI, cntA, mtmB), and the assembled fraction of the guild genes in the unrestored cores (0.33–0.62) is in the range of the universal markers (0.47–0.56). The SF2 / R2A sediment cores assemble only 11–44 % of their genomes, and their guild genes are even less assembled (≤ 0.28).
- mtgB (1.6 in the unrestored cores) and grdH (0.70; up to 2.0 per sample) exceed 1 because the read survey under-counts them (82 % of mtgB reads are called mttB, 47 % of grdH reads grdB; read_markers validation); compare the merged MttB superfamily and GrdB family there. mtaB / mtaA / mtbA are hardly assembled at KOfam thresholds, and the read side also counts KOfam-negative ("twilight") homologs of these genes.

## Caveats

- **Assembly limits.** 39–80 % of the sediment marker reads do not assemble into contigs ≥ 1 kb (read_markers); here the assembled GE is 11–62 % of the read GE in the sediment cores. Assembly-based copies per read GE are therefore lower bounds; a KO's assembled fraction differs by habitat (above: phnJ 0.54 in the unrestored cores vs 0.08 in the SF2/R2A cores) and by gene (mtmB, mcrA, cutC assemble worse than the markers in several sediments). Normalising by the assembly-internal GE (`copies_per_GE_assembly`) assumes the gene is assembled as well as the markers, which this check shows is not generally true.
- **Depth** is the own-library mean contig depth of `jgi_summarize_bam_contig_depths` (reads ≥ 97 % identity), so strain variants below that identity are not counted, and a gene inherits its contig's mean depth.
- **Only the best KO per gene is counted** (DRAM rule); a gene that passes a panel KO but has a better non-panel KO is counted for that other KO (0 such genes in the 4 checked chunks).
- **KOfam rule inheritance**: `-T 14.7` reporting cut (3 domain-type KOs with thresholds < 15.7 could lose a hit inside a 1-bit window; none is a panel KO), no Sec read-through rescue (grdH/grdB in Bacillota can be split by Prodigal at the selenocysteine codon and under-called; the MAG track rescued 34 such calls), KOs with threshold `-` never called (five panel KOs: supplement).
- **Partial genes** on short contigs are counted like complete ones if they pass the threshold; short fragments that do not pass are lost.
- **Domain split** covers binned contigs only; other-bin domains come from CheckM1 marker-set lineages (checked against GTDB on the kept bins, not on the other bins themselves). Archaeal ≠ methanogen.
- **E-values** use `-Z 12,055,118` (MAG tracks: 2,732,576 and 3,911,146), so E-values differ from the MAG tracks by a constant factor; scores, passing sets and best KOs do not (check 4).
- `ko_by_sample` covers only the 148 panel + marker KOs (two-pass design).

## Problems

- The full all-KO search was too slow at poplar's load (750–1,150 on 336 cores from 16:00 CST; ≈ 3.4–4 h per chunk instead of the ≈ 2.2 h measured at launch) and was stopped after 3.6 h on the coordinator's request (≈ 325 thread-hours spent, 4 of 521 chunks finished); replaced by the two-pass design.
- The first aggregation used the all-profile reading of the marker rescue, which under-counts archaeal uS9/uS11/uL2 (see Genome equivalents); the aggregation was re-run with v1's reading.
- One status monitor (my own, on this workstation) briefly used `pgrep -u afreiburger -fc kofam_all.hmm` inside an ssh command whose own command line contained that string, against the brief's rule on self-matching `pgrep -f` patterns. It only counted processes, killed nothing, and was stopped within seconds and replaced by `pgrep -x hmmsearch`. Processes were stopped only by PID from my driver's process tree.

## Files

| here (`community_ko/assembly_based/`) | content |
|---|---|
| `scripts/` | every script that ran (identical copies on poplar): `00_smoke_compare.py`, `01_make_chunks.py`, `02_postprocess_chunk.py`, `run_all.sh` + `run_chunk.sh` (stopped full search), `run_pass1.sh` + `p1_run_chunk.sh`, `08_select_candidates.py`, `run_pass2.sh` + `p2_run_chunk.sh`, `03_aggregate.py`, `04_validate_vs_rollout.py`, `07_nothreshold_panel.py`, `09_two_pass_check.py`, `05_panel_tables.py`, `06_validate_vs_reads.py`, `v1_panel_kos.txt` |
| `results/ko_by_sample.tsv.gz` | 148 KOs × 39 samples: n_genes, sum_depth, copies_per_GE_assembly, copies_per_GE_readGE, gene_fraction |
| `results/ko_status_domain_by_sample.tsv.gz` | the same by contig status × binned domain |
| `results/ge_assembly.tsv`, `marker_rescue_check.tsv` | GE per sample (8 family depths, mean, median, CV; read GE; both rescue readings) |
| `results/panel_status_split.tsv` | 138 panel KOs × 39 samples: kept MAG / other bin / unbinned, binned domain split |
| `results/panel_mag_attribution.tsv`, `_summary.tsv`, `_check.txt` | carrier species of each panel KO (final 585), genus, mean read fraction in the 12 unrestored-core samples |
| `results/nothreshold_panel_*.tsv` | supplement for the five no-threshold panel KOs |
| `results/validation_vs_read_survey.tsv`, `_summary.tsv` | check against the v1 read survey |
| `results/validation_vs_kofam_rollout.txt`, `_diffs.tsv` (empty) | gene-level check against the MAG track |
| `results/two_pass_equivalence_check.txt`, `_diffs.tsv` (empty) | two-pass vs full search |
| `results/aggregate_qc.txt`, `.tsv`, `checkm1_lineage_domain_*.tsv` | counts per sample, domain map |
| `validation/smoke/` | smoke test |
| `poplar_copies/` | logs (`run_all.log`, `run_pass1.log`, `run_pass2.log`, `03_aggregate*.log`, `04/07/08/09` logs, `kill_tree.txt`), chunk statistics, `pass2/candidates.tsv.gz` + `chunk_stats.txt`, script md5s at each launch, the first aggregation's GE table |

On poplar only (`/scratch1/afreiburger/salterns_2026/community_ko/`): `genes_ko/<s>.tsv.gz` (**per gene**, all 12,055,118 genes: sample, contig, representative, `searched_full_db`, best KO, score, all passing KOs, marker rescue (both readings), top-scoring KO, contig depth, status, MAG, domain; 231 MB), `meta/genes/`, `chunks/` (521 chunks, 2.7 GB), `pass1/` (tblout + domtblout + `hits.tsv.gz`), `pass2/` (candidate chunks, tblouts, calls), `full_run_stopped/` (4 finished + 86 partial full-DB chunks, 4.3 GB), `nothr/`, `smoke/`, `test_agg/` (code-path test), `logs/`.
