# Read-based survey of methane-relevant genes in the 39 salterns metagenomes

*Track `read_markers`, set up and run on poplar on 2026-11-24 (all times CST). Nothing was committed, nothing was written to KBase, and no other user's or session's process was touched. Poplar work dir: `/scratch1/afreiburger/salterns_2026/read_markers/` (reads, DIAMOND hit tables, gene calls and the reference DB stay there; this folder holds the scripts, the small metadata, the validation tables and the results).*

**Status: DONE** (2026-11-24 17:46 CST). Results below.

## Why

In the sediments, 44–83 % of the reads do not assemble into contigs ≥ 1 kb, and the final 585 species MAGs recruit only 11–49 % of the reads. MAG-based guild abundances therefore cannot say how much of a gene lies outside the MAGs. This track counts the genes directly in the reads, then asks where the counted reads sit: in a kept MAG, in another bin, on an unbinned contig, or nowhere in the assembly.

## Target genes (KEGG KOs; final list in `meta/ko_panel.tsv`)

| group | genes (KO) | note |
|---|---|---|
| methanogenesis | mcrA K00399, mcrB K00401, mcrG K00402 | K00399 also holds the alkyl-CoM reductases (Acr) of ANME / alkane oxidisers; the read taxonomy shows which |
| methylotrophic methyltransferases | mtaB K04480, mtaA K14080, mtbB K16178, mtmB K16176, mtbA K14082 | |
| MttB family (reported separately, never called methanogenesis) | mttB K14083, mtgB K25220, `mttB_family_other` | K14083 includes the pyrrolysine-free MttB superfamily that is widespread in bacteria; `mtgB` (glycine betaine MT) is a characterised Pyl-free member; `mttB_family_other` = UniProt Pfam PF06253 members that pass neither KO. Reported as one merged class (see validation) |
| C–P lyase | phnG K06166, phnH K06165, phnI K06164, phnJ K06163, phnK K05781, phnL K05780, phnM K06162 | phnJ is the guild marker |
| betaine reductase | grdH K21579, grdI K21578 | grdI is the read-level proxy (grdH reads are not separable from grdB, see validation) |
| glycine reductase (decoy class, reported) | grdB K10672 | paralog of grdH |
| TMA production | cutC K20038, cutD K20037, cntA K22443, cntB K22444 | |
| MePn synthesis | mpnS K18049 | KOfam has no threshold; see the mpnS/HEPD rule; not validated (too few references) |
| methanotrophy | pmoA (K10944, methanotroph references), mmoX K16157 | K10944 lumps pmoA with amoA; split by reference taxonomy, `amoA` and `cuMMO_other` reported as decoy classes and merged into `CuMMO_all` |

**Decoy classes** (close paralogs; a read whose best hit is a decoy is not counted for the target): grdB/grdE/grdF/grdG vs grdH/grdI; PflB, HpdB, BssA, GrcA and other glycyl-radical enzymes (`GRE_other`) vs CutC; PflA/PflX vs CutD; Cmo, VanA, GbcA and other Rieske oxygenases (`Rieske_other`) vs CntA; VanB, GbcB vs CntB; PrmA, TmoA, alkene MO, BmoX and other soluble di-iron monooxygenases (`SDIMO_other`) vs MmoX; amoA and other CuMMOs vs pmoA; HEPD vs MpnS; HemE (uroporphyrinogen decarboxylase family) vs MtaA/MtbA; PhnC, DdpD, DdpF, ABC.CD.A, ABC-2.A ATPases vs PhnK/PhnL; eukaryotic ribosomal proteins (`euk_marker`) vs the markers; and `<gene>_like` = salterns proteins that hit a target profile but score < 0.5 × its KOfam threshold (complete genes only; e.g. 37,066 ABC ATPases that are not PhnK).

**Universal single-copy markers (normalisation).** Eight universal ribosomal-protein families present once in every bacterial and archaeal genome: uL2 (K02886), uL11 (K02867), uS3 (K02982), uS9 (K02996), uS11 (K02948), uS12 (K02950), uL4 (K02926 bacterial + K02930 archaeal) and uL16 (K02878 bacterial + K02866 archaeal).
- Choice: in the 473 HQ rollout MAGs these KOs (or KO pairs) are present in ≥ 0.82 of archaea and ≥ 0.94 of bacteria and are multi-copy in ≤ 1.4 %. KEGG GENES assigns archaeal uL4 and uL16 to the eukaryote-named KOs K02930 / K02866 (checked on *Halobacterium*, *Methanosarcina*, *Nitrosopumilus*), so those two are merged pairs.
- Not used: L3, L22, L24, S10, S19 of the Hug et al. 16-protein set, which KOfam under-calls in the salterns archaea (0–10 % of HQ archaeal MAGs), and S17 (archaea carry both uS17 and eS17; *Methanosarcina* has two K02961 genes).
- GE depth = mean of the eight families' read depths (median also reported, and every family's depth is in `results/genome_equivalents.tsv` so a deviating marker is visible).

## What ran (in order; scripts in `scripts/`, identical copies on poplar)

| step | script | what | where / threads |
|---|---|---|---|
| 1 | `01_prodigal.sh` | Prodigal 2.6.3 `-p meta -g 11 -q` (the kofam_rollout flags) on **all contigs ≥ 1 kb** of the 39 per-sample rollout assemblies (`rollout/<s>/asm/contigs.fa`, each split in 8 pieces): **12,055,118 genes** | poplar, 40 × 1 CPU, nice 10 |
| 1b | `01b_prodigal_coasm.sh` | the same on the **520 kept co-assembly MAGs** (`coasm_qc/mags_kept/*.fa`), so reads can be attributed to every MAG of the final 585-species set: 1,585,384 genes | poplar, 16 × 1 CPU |
| 2 | `02_fetch_refs.py` | external references (see "Reference database") | local (internet) |
| 3 | `03_panel_hmmsearch.sh` | HMMER 3.4 `hmmsearch --cpu 1 --noali -E 1e-5 --domE 1e-5 -Z 20000000` of the 62 KOfam 2026-10-26 profiles of the panel (`db/panel.hmm`, cut from `kofam_db/profiles/`) against the 39 assemblies, the co-assembly MAGs and the external set | poplar, 24 × 1 CPU |
| 4 | `04_build_db.py` | `ext` (external FASTA + taxonomy), `status` (`meta/contig_status.tsv.gz`: every contig ≥ 1 kb of every sample → kept MAG / other DAS Tool bin / unbinned, from `rollout/<s>/keep/bins.tsv`), `label` (KOfam calls → class labels) | poplar, 1 CPU |
| 5 | `05_cluster_db.sh` | K10944 resolution, per-label MMseqs2 clustering (99 % id for targets/markers/reported classes, 95 % for decoys; 90 % for the hold-out units), DIAMOND DB | poplar, 16 threads |
| 5b | `05b_rep_info.py` | per-representative taxonomy and MAG membership (`refs/rep_info.tsv.gz`), gene lengths (`refs/label_lengths.tsv`, overrides in `meta/length_overrides.tsv`) | poplar |
| 7–8 | `07_simulate.py`, `run_validation.sh`, `08_classify.py validate`, `08a_exclude_reps.py`, `08b_validation_report.py`, `08c_pick_cutoffs.py` | hold-out validation and cut-off choice (`cutoffs.tsv`) | poplar, 16 threads |
| 6 | `06_diamond_reads.sh` | DIAMOND 2.2.8 `blastx --mid-sensitive -k 5 -e 1e-3 -b 4 -c 1` of **every QC read** (`rollout/reads/<s>/map_1.fq.gz`, `map_2.fq.gz`; each mate a separate query; 2,084,152,888 reads) against the survey DB | poplar, 8 jobs × 8 threads (6 + 1 + 1 drivers) |
| 8 | `08_classify.py sample` (via `run_post.sh`) | best-hit classification with decoys at the chosen cut-offs → `results/per_sample/<s>.counts.tsv`, `<s>.reads.tsv.gz` | poplar, 1 CPU per sample |
| 9 | `09_map_assigned.sh` (via `run_post.sh`), `09b_check_fullbam.sh` | re-mapping of the classified reads to their own sample's contigs (bowtie2 2.5.5 / samtools 1.25, as the rollout) and the check against the two surviving full BAMs | poplar, 4 threads per sample |
| 10d | `10d_methanol_taxonomy.py` (2026-11-25) | genus of the *mtaB* and *mtaA* reads (not covered by 10b) → `results/methanol_taxonomy.tsv`. In the 12 unrestored-core samples, 21 % of the 2,559 *mtaB* reads are archaeal (*Methanohalophilus* 10 %, *Methanococcoides* 6 %); the rest are bacterial (Bacteroidota *Draconibacterium* 27 %), with 43 % of all *mtaB* reads below 70 % identity | poplar, 1 CPU, nice 10, about 20 min |
| 10 | `10_summarize.py`, `10b_split_and_taxonomy.py`, `10c_key_tables.py` | the result tables | poplar, 1 CPU |

Thread budget: every poplar job ran at nice 10, detached in tmux, with `TMPDIR` and DIAMOND `-t` under `read_markers/tmp/` on /scratch1 (never /tmp). This track's thread total stayed ≤ 64, apart from three one-thread checks of ≤ 2 min (see "Problems") (load 180–540 on 336 cores during the run, with other jobs of the same account, which were not touched). `run_post.sh` starts a post-processing job only while 8 × (DIAMOND slots still needed) + 5 × (running post jobs) + 5 ≤ 64.

**Read files.** All 39 samples' `map_1/map_2.fq.gz` exist under `rollout/reads/<s>/` (18 are symlinks into `binning_test/reads/`); the qtrim `result_reads` in each `read_counts.json` equals the sample's total in `abundance_final/qc_read_totals.tsv` (asserted in `10_summarize.py`). Mean read length 147.3–149.0 nt.

## Reference database (`db/survey.faa`, 137,012 representatives; md5 in `poplar_copies/db_md5.txt`)

Sources, all labelled with the KOfam 2026-10-26 thresholds and the DRAM best-KO rule (lowest full E-value among passing KOs):

1. **KEGG GENES** members of every panel KO (KEGG REST `link/genes/<KO>` + `get/<genes>/aaseq`, fetched 2026-11-24; organism taxonomy from KEGG `br08601`). Prokaryote genes only, except for the marker KOs, whose eukaryote genes (≤ 1,000 per KO) become the `euk_marker` decoy class. Large KOs were randomly subsampled (`random.Random(1)`): 3,000 genes for pflB, pflA, hemE and the four ABC ATPase KOs; for the markers every archaeal gene + 3,000 bacterial genes. The KEGG KO is the label of a KEGG gene (KEGG orthology, not a KOfam call). 84,414 sequences.
   - Deviation from the brief: UniProt no longer exposes KO cross-references (`xref:ko-K00399` returns 0 hits in UniProt 2026_04, and the ID-mapping service maps only ~55 % of KEGG genes to UniProtKB), so KEGG GENES is the KO-defined source.
2. **UniProtKB (2026_04) Pfam families** around selected targets, Bacteria + Archaea, for breadth and taxonomy: PF02249/PF02745 (McrA), PF02241/PF02783 (McrB), PF02240 (McrG), PF06007 (PhnJ), PF02461 (CuMMO subunit A), PF02332 (SDIMO), PF07355 (GrdB family), PF02901 (glycyl-radical enzymes), PF00848 (Rieske oxygenase α), PF12176 (MtaB), PF06253 (MttB), PF09505 (MtbB), PF05369 (MtmB), plus the name queries for MpnS and HEPD. 46,753 sequences fetched, 31,819 labelled (14,934 passed nothing and fell under no family rule). Each member is labelled by its passing KOfam KO; if none passes, family rules apply (McrA/B/G, PhnJ, MtaB/MtbB/MtmB families: the family's KO if it is the top-scoring profile with E ≤ 1e-10; SDIMO / GRE / Rieske / MttB families: `SDIMO_other` / `GRE_other` / `Rieske_other` / `mttB_family_other` unless the top profile is the target KO; PF02461 by taxonomy into pmoA / amoA / cuMMO_other).
3. **Salterns proteins** of the 39 assemblies (all contigs ≥ 1 kb) and the 520 co-assembly MAGs that pass a panel KO. Markers: if nothing passes, the top-scoring marker profile at ≥ 0.5 × its threshold also counts (KOfam under-calls archaeal ribosomal proteins). Decoy-KO proteins are kept only if they also hit a target profile (E ≤ 1e-5). Proteins that hit a target profile but pass nothing are listed in `refs/unassigned_homologs.tsv`; complete ones scoring < 0.5 × the target threshold become `<gene>_like` decoys; the 0.5–1 × "twilight" ones are left out of the DB (and used in the validation). MpnS/HEPD (no threshold): top-scoring panel profile K18049 or K12905, full E ≤ 1e-30, HMM coverage ≥ 0.5. 145,441 labelled proteins (163,835 target homologs passed nothing).
4. K10944 proteins of the salterns set (1) are resolved to pmoA / amoA / cuMMO_other by their best DIAMOND hit among the labelled external CuMMO references.

Each label is clustered separately (MMseqs2 easy-cluster, coverage ≥ 0.8): 99 % identity for targets, markers and reported classes, 95 % for decoys. Every representative keeps its member list (`refs/clu/clusters99.tsv`): the source of each member, and for salterns members the sample, contig status, MAG and species cluster (`comparison_rollout/mags_three_way.tsv`).

**Classification-time exclusions** (`refs/excluded_reps.tsv`): 284 marker representatives longer than their family's 99th length percentile (fusions or mis-predicted genes: in the validation, background reads of non-ribosomal genes hit e.g. a 263-aa "uS11" and a 344-aa KEGG "uS9" at 50–70 % identity), and 69 MpnS/HEPD name-query references that fail the K18049/K12905 rule (`gene:hepD` also returns heptaprenyl diphosphate synthase genes; two "MpnS" entries have no K18049 hit). Hits to them are ignored.

## Classifier and cut-offs (`cutoffs.tsv`)

Per read (each mate separately), from DIAMOND's top 5 references (E ≤ 1e-3 reported):
1. keep hits with **E ≤ 1e-5** and **aligned length ≥ 30 aa** (90 of the read's ~148 nt), ignoring the excluded representatives;
2. the best hit (highest bitscore) gives the label;
3. if a hit of a *different* label scores as high or higher (margin 0 bits) and either label is a counted class, the read is **ambiguous** and not counted;
4. otherwise the read counts for the label if the best hit's identity is ≥ the label's minimum: **none for most labels**, **50 % for the eight markers**, **60 % for cutC, cntA, cntB and mpnS**.

Why these values (from the hold-out grid, `validation/grid_*.tsv`, `grid_summary.tsv`):
- E ≤ 1e-5 vs 1e-3, and aligned length 25 vs 30 aa, changed macro recall/precision by < 0.5 %; 30 aa / 1e-5 is the more conservative pair.
- A bitscore margin of 2–5 bits lowered recall (0.890 → 0.883) without raising precision.
- A global identity cut lowered recall without raising precision for most labels (e.g. phnJ precision 0.962 at 0 %, 0.966 at 60 %), so it is applied only where it removed specific false positives:
  - markers at 50 %: within the background genes, reads of non-ribosomal proteins hit marker references at 40–74 % (after the length-outlier exclusion, 2 such reads remain out of 138,035 background reads);
  - cutC / cntA / cntB at 60 %: reads of the KOfam-unassigned ("twilight", 0.5–1 × threshold) glycyl-radical enzymes and Rieske oxygenases that are not in the DB were assigned to cutC / cntA at < 60 % (cutC: 248 → 61 reads; cntA 60 → 14); recall costs 1.5 % (cutC) and 6 % (cntA);
  - mpnS at 60 %: no hold-out data (too few references), so a conservative cut.
- DIAMOND sensitivity: `--mid-sensitive` was used for production. On the hold-out reads it matched `--sensitive` within 2 % recall for every label, and beat the default mode for the divergent ones (mtmB 0.901 vs 0.821, phnH 0.943 vs 0.903, cutD 0.864 vs 0.832). Cost on real reads (8 M-read benchmark, 16 threads): default 2:02 min, mid-sensitive 2:52 min; 14 % more reads with any hit.

## Validation (`validation/`)

**Design** (`07_simulate.py`).
- Hold-out unit: a 90 %-identity cluster of one label (MMseqs2, coverage 0.8). Per label, 20 % of its 90 % clusters (≥ 1 when ≥ 3 clusters, ≤ 60) were held out (`random.Random(1)`): 3,148 clusters.
- Validation DB: the production representatives minus every representative in a held-out cluster (131,677 of 137,012), so no held-out gene has a ≥ 90 %-identical relative of its own label left. The median identity of a held-out gene to its nearest remaining same-label reference is 70–90 % (column `median_divergence_identity`).
- Reads: 148 nt, uniform start, random strand, 0.3 % substitutions, ~5× coverage, from one member per held-out cluster: 862 salterns genes (CDS ± 150 nt of their own contig, so reads straddle gene ends as real ones do) and 1,524 KEGG genes (CDS from KEGG `ntseq`). A read is "true" for the label if it overlaps the CDS by ≥ 90 nt.
- Negatives: 4,000 random salterns genes with no panel-HMM hit (± 150 nt; 138,035 reads), and 142 KOfam-"twilight" target homologs (0.5–1 × threshold, not in the DB; 6,216 reads, truth unknown, reported separately as `uh_assigned`).
- 227,672 simulated reads, run against the validation DB with the production DIAMOND flags in three modes.

**Per-label results at the chosen cut-offs, `--mid-sensitive`** (`validation/validation_per_label_mid-sensitive.tsv`; recall by divergence in `validation_recall_by_divergence_mid-sensitive.tsv`; confusion in `validation_confusion_mid-sensitive.tsv`). Precision = TP / (TP + reads of other held-out labels or background genes assigned to the label) in the simulated mix.

| label | held-out genes | reads | recall | precision | note |
|---|---|---|---|---|---|
| mcrA | 9 | 499 | 0.998 | 1.000 | recall 1.0 down to 70–80 % identity to the nearest reference |
| mcrB / mcrG | 19 / 17 | 833 / 427 | 1.000 / 0.995 | 1.000 / 1.000 | |
| mtaB / mtaA | 15 / 28 | 689 / 950 | 0.977 / 0.883 | 1.000 / 0.926 | mtaA loses 5 % to HemE, 4 % to mtaA_like |
| mtbB / mtmB / mtbA | 10 / 19 / 7 | 380 / 756 / 267 | 0.984 / 0.901 / 0.689 | 0.906 / 0.987 / 0.944 | mtbA: 27 % of its reads go to mtaA or HemE |
| mttB / mtgB | 24 / 27 | 1,205 / 1,248 | 0.945 / 0.106 | 0.499 / 0.956 | 82 % of the mtgB reads are called mttB → **reported only as the merged MttB superfamily** |
| phnG / phnH / phnI / phnJ | 57 / 59 / 50 / 23 | 858 / 1,160 / 1,818 / 685 | 0.954 / 0.943 / 0.973 / 0.997 | 0.975 / 0.987 / 0.981 / 0.962 | |
| phnK / phnL / phnM | 58 / 57 / 55 | 1,515 / 1,389 / 2,114 | 0.986 / 0.970 / 0.989 | 0.986 / 0.992 / 0.991 | the 10 background "false positives" are flank reads of neighbouring phn genes |
| grdI | 15 | 667 | 0.992 | 0.934 | read-level proxy for betaine reductase |
| grdH | 12 | 460 | 0.439 | 0.701 | 47 % of the grdH reads are called grdB, 9 % of grdB reads grdH → **grdH is not used as a read-level proxy**; reported with grdB as the merged GrdB family |
| cutC / cutD | 21 / 25 | 1,750 / 770 | 0.983 / 0.864 | 0.993 / 0.992 | cutD loses 10 % to PflA |
| cntA / cntB | 23 / 15 | 830 / 473 | 0.912 / 0.898 | 0.977 / 0.977 | |
| pmoA | 6 | 148 | 0.601 | 1.000 | 24 % of pmoA reads go to `cuMMO_other`: the taxonomy-based pmoA/amoA split is soft → **reported as the merged CuMMO class with read taxonomy** |
| mmoX | 1 | 52 | 1.000 | 1.000 | one held-out gene only |
| markers uL2, uL11, uS3, uS9, uS11, uS12, uL4, uL16 | 46–55 each | 673–1,268 | 0.911–0.996 (mean 0.963) | 0.956–0.997 | background "false positives" are almost all flank reads of neighbouring ribosomal-protein genes (operons); 2 reads inside non-ribosomal CDSs |

- **Background genes**: at the chosen cut-offs only the markers and phnK/phnL received background reads (62 reads in total), and 60 of those overlap the background gene's own CDS by < 90 nt, i.e. they are reads of the neighbouring gene in a ribosomal-protein or phn operon. Inside non-target CDSs: 2 reads of 138,035 (≈ 14 per million).
- **Twilight homologs** (not in the DB, truth unknown): reads assigned mostly to mtaB, mtbB, mtbA, phnH, phnI at ≥ 70 % identity, i.e. these genes are close to KEGG references of those KOs but fall below the KOfam threshold. They would be counted by the read survey; this is a definitional difference to KOfam, not a classifier error, but it means read counts of mtaB, mtbA, mtbB, phnH and phnI can include KOfam-negative homologs.
- **Not validated**: mpnS (15 references, no hold-out), hepD, amoA/cuMMO_other (1–2 held-out genes).

## Normalisation and units (`10_summarize.py`)

Per sample *s* (reads = single mates; both mates counted, as in `qc_read_totals.tsv`):
- `L_nt` = 3 × the median aa length of the label's complete references (KEGG genes; `refs/label_lengths.tsv`, overrides for fragment-dominated classes in `meta/length_overrides.tsv`);
- `r_s` = mean QC read length; `L_eff = L_nt + r_s − 2 × 90 + 1` (read start positions that give ≥ 30 aligned aa);
- `depth = reads × r_s / L_eff` (summed coverage of all copies of the gene family);
- **GE depth** = mean `depth` of the eight marker families (median reported too, and each family's depth);
- **copies per GE** = `depth / GE depth` (≈ the fraction of genomes carrying the gene, if single-copy);
- **RPKM** = `reads / (L_nt / 1e3) / (total QC reads / 1e6)`.

**MAG side**, same unit: for each of the 585 species, `genomes per GE = (CoverM reads / MAG size in bp) / (GE depth / r_s)` (both sides in reads per bp; CoverM: ≥ 95 % identity, ≥ 75 % aligned, `abundance_final/final585_read_count.tsv`); a guild's value is the sum over its species (`correlational/final585/guilds_final585.csv`). The manuscript's unit, the read fraction (`final585_relabund_24sediment.tsv`), is given alongside. A read fraction underweights small genomes: *Methanolobus* (2.32 Mb) gets a smaller read share than a 3.5 Mb genome at the same cell abundance, so read-based copies per GE must be compared with MAG genomes per GE, not with the read fraction.

Read-level guild proxies: methanogen = mcrA; methylamine user = mean(mtbB, mtmB); MePn consumer = phnJ; betaine reducer = grdI (validated; grdH is not separable from grdB at read level); TMA producer = cutC + cntA.

## Deviations from the brief

1. **UniProt `xref:ko-` no longer works** (0 hits in UniProt 2026_04), so the KO-defined references come from KEGG GENES (REST), with UniProt Pfam families added for breadth and labelled by KOfam.
2. **No full re-mapping.** The rollout deleted its BAMs after the depth tables (only R1_A_D1 and SF2_C_D1 keep theirs, in `binning_test/`). Instead (a) the classified reads themselves are re-mapped with bowtie2 to the contigs of their own sample that carry any panel-HMM hit (the "reduced reference"; checked against the two full BAMs), and (b) the assembly-based complement uses the existing own-library depth tables (`depth/<s>.d1000.txt`).
3. **Markers**: 8 universal ribosomal-protein families instead of a larger set, chosen for single copy in both domains in the salterns MAGs; uL4 and uL16 are merged bacterial + archaeal KO pairs.
4. **K10944 and K14083 split / merged**: pmoA vs amoA by reference taxonomy, then reported as one CuMMO class with read taxonomy; mttB, mtgB and the PF06253 'other' members reported as one MttB superfamily (validation shows they cannot be separated at read level). grdH is reported, but grdI is the betaine-reductase proxy.
5. **Read survey on all reads**: no subsampling was needed (DIAMOND mid-sensitive on all 2,084,152,888 reads took 4 h 22 min at 64 threads, 13:19–17:41 CST; 15–38 min per mate).

## Files

| here | content |
|---|---|
| `scripts/` | every script that ran (identical copies on poplar) |
| `meta/ko_panel.tsv` | the KO panel: KO, gene / class, role (target / decoy / decoy_reported / marker), group, KEGG cap |
| `meta/length_overrides.tsv` | gene lengths for fragment-dominated reported classes |
| `cutoffs.tsv` | the classification cut-offs |
| `validation/` | hold-out grids (3 DIAMOND modes), per-label metrics, recall by divergence, confusion, twilight-homolog assignments, held-out genes |
| `refs_small/` | `rep_info.tsv.gz` (per DB representative: label, members, MAGs, species clusters, genus), `label_lengths.tsv`, `excluded_reps.tsv`, `clusters99.tsv.gz`, `labeled_meta.tsv.gz`, `unassigned_homologs.tsv.gz`, `labels.tsv`, `k10944_resolution.tsv`, `fetch_log.tsv` |
| `results/per_sample_KO.tsv` | **per sample × label**: reads, ambiguous / low-identity best hits, L_nt, L_eff, depth, RPKM, copies per GE (mean- and median-normalised); plus three merged classes |
| `results/genome_equivalents.tsv` | per sample: GE depth (mean, median, CV) and each marker family's depth |
| `results/per_guild_read_vs_MAG.tsv` | per sample × guild: MAG read fraction, MAG genomes per GE, read-based proxy copies per GE, fold |
| `results/key_unrestored_cores.tsv`, `key_unrestored_cores_summary.tsv` | the 12 samples of the six unrestored cores: read vs MAG per guild, placement of the proxy reads, MAG-attributable fraction |
| `results/read_taxonomy.tsv`, `mcrA_taxonomy_by_sample.tsv`, `mcrA_taxonomy_unrestored_pooled.tsv`, `CP_lyase_taxonomy_unrestored_pooled_*.tsv` | genus of the classified reads (best hit), identity bin, placement |
| `results/binned_unbinned_split.tsv`, `split_per_gene_by_habitat.tsv` | per sample × label: reads on kept-MAG contigs / other bins / unbinned contigs / unassembled |
| `results/mag_attributable_reads.tsv` | reads whose best-hit cluster contains a MAG of a guild species at ≥ 95 % identity |
| `results/assembly_depth_complement.tsv`, `assembly_vs_reads.tsv` | KOfam-called assembled genes × own-library contig depth, by contig status, vs the read survey |
| `results/fullbam_check.tsv` | reduced-reference re-mapping vs the two full BAMs |
| `results/per_sample/*.counts.tsv` | per-sample classification counts (the per-read tables `*.reads.tsv.gz` stay on poplar) |
| `poplar_copies/` | logs (`logs/`), DIAMOND `/usr/bin/time` reports (`*.time`), re-mapping summaries, `db_md5.txt` (survey DB, panel HMM, labelled FASTA), `scripts_md5.txt` |

On poplar only (large): `genes/` (Prodigal calls, 20 GB), `hmm/`, `refs/` (labelled FASTA), `db/survey.{faa,dmnd}`, `dmnd/` (DIAMOND hit tables), `map/`, `results/per_sample/*.reads.tsv.gz`, `val/` (simulated reads).

## Results

**Scale.** 2,084,152,888 reads of 39 samples were classified: 2,177,966 marker reads, 1,073,932 target reads, 37,009 reads of reported decoy classes. The eight marker families agree well within each sample: the coefficient of variation of their depths is 3.8–13.3 % (median 6.8 %; `results/genome_equivalents.tsv`), so the GE normalisation is not driven by one marker.

**Placement check.** For the two samples with a surviving full self-mapping BAM, the placement of every classified read (kept MAG / other bin / unbinned / unassembled) from the reduced-reference re-mapping agrees with the full BAM for 98.95 % (R1_A_D1, 104,514 reads) and 99.62 % (SF2_C_D1, 83,671 reads) of reads; unassembled 34.5 vs 34.1 % and 85.7 vs 85.8 % (`results/fullbam_check.tsv`). The assembly-depth complement (KOfam-called assembled genes × own-library contig depth, `results/assembly_vs_reads.tsv`) gives assembled fractions of 0.33–0.56 of the read-based copies per GE in the 12 unrestored cores, close to the re-mapping's 0.35–0.63 (it is lower because jgi depths count only ≥ 97 %-identity reads and miss below-threshold partial genes).

### Read-based vs MAG-based guild abundance in the six unrestored cores (12 samples: R1_A/B/C and R2_A/B/C, D1 and D2)

Units: copies (read side) or genomes (MAG side) per genome equivalent; `results/key_unrestored_cores.tsv`, `key_unrestored_cores_summary.tsv`, all samples in `per_guild_read_vs_MAG.tsv`.

| guild (read proxy) | MAG species | MAG read fraction, median (%) | MAG genomes per GE, median (%) | read copies per GE, median (%) | fold read / MAG, median (range) | fold, pooled over 12 | proxy reads on kept-MAG / other-bin / unbinned contigs / unassembled (%) | proxy reads attributable to the guild's MAGs (≥ 95 %) |
|---|---|---|---|---|---|---|---|---|
| methanogen (mcrA) | 1 | 0.08 | 0.13 | 0.21 | 1.8 (0.9–14.6) | 1.74 | 27 / 6 / 10 / 57 | 49 % |
| methylamine user (mean mtbB, mtmB) | 1 | 0.08 | 0.13 | 1.18 | 8.6 (4.4–84.6) | 7.36 | 19 / 5 / 15 / 60 | 25 % |
| MePn consumer (phnJ) | 39 | 5.78 | 6.23 | 12.23 | 2.1 (1.7–2.4) | 2.04 | 29 / 9 / 20 / 42 | 43 % |
| betaine reducer (grdI) | 18 | 2.42 | 3.60 | 6.28 | 1.7 (1.3–2.4) | 1.67 | 42 / 10 / 12 / 37 | 72 % |
| TMA producer (cutC + cntA) | 0 | 0 | 0 | 6.96 | — (no MAG species) | — | 34 / 7 / 16 / 43 | — |
| one-sided TMA candidates (cutC + cntA) | 36 | 4.89 | 5.40 | 6.96 | 1.4 (1.2–3.2) | 1.38 | 34 / 7 / 16 / 43 | 38 % |

- "kept-MAG contigs" = contigs of a kept MAG of the sample's own assembly; "attributable" = best-hit cluster contains a protein of a MAG of the guild's species (any sample, rollout or co-assembly) at ≥ 95 % identity.
- **Methylamine users**: the read proxy is not methanogen-specific. 32 % of the mtbB and 54 % of the mtmB reads are archaeal (*Methanolobus*, *Methanohalophilus*, *Methanococcoides*); the rest are bacterial Pyl-methyltransferases (Bacillota JAAYHV01, *Formimonas*, *Acetohalobium*, Desulfobacterota). Methanosarcinales also carry several mtmB/mtbB copies per genome. So the 7–85× fold is not a methanogen fold; for methanogens use mcrA.
- **TMA**: no MAG species is a two-sided TMA producer, but cutC + cntA reads amount to 1.7–9.5 % of genomes. cntA reads are 94 % bacterial (*Roseovarius* 46 %, *Marivita* 21 %), cutC reads 99.9 % bacterial (Bacillota BM714 39 %). Note the cntA KO also holds YeaW-type choline/betaine monooxygenases.
- Other genes per sample × KO (RPKM, copies per GE): `results/per_sample_KO.tsv`. pmoA: 56 reads in all 39 samples (mostly *Methylotuvimicrobium*-like), mmoX 2, amoA 2, mpnS 52 in the 12 cores (best hits to dinoflagellate MpnS; unvalidated); all negligible.

### Is total mcrA several-fold the *Methanolobus* MAG's share?

**Not overall.** Pooled over the 12 unrestored-core samples, mcrA copies per GE are **1.74 ×** the *Methanolobus* MAG's genomes per GE (median 1.8 ×). Per sample:

| sample | mcrA reads | mcrA copies per GE (%) | *Methanolobus* MAG read fraction (%) | *Methanolobus* genomes per GE (%) | fold | mcrA reads of genus *Methanolobus* | MAG-attributable (≥ 95 %) | on kept-MAG contigs | unassembled |
|---|---|---|---|---|---|---|---|---|---|
| R1_A_D1 | 154 | 0.58 | 0.260 | 0.395 | 1.5 | 69 % | 59 % | 68 % | 24 % |
| R1_A_D2 | 216 | 0.84 | 0.267 | 0.416 | 2.0 | 53 % | 44 % | 55 % | 45 % |
| R1_B_D1 | 49 | 0.19 | 0.072 | 0.115 | 1.6 | 73 % | 49 % | 0 % | 57 % |
| R1_B_D2 | 72 | 0.28 | 0.066 | 0.103 | 2.7 | 62 % | 53 % | 0 % | 36 % |
| R1_C_D1 | 59 | 0.31 | 0.131 | 0.212 | 1.4 | 73 % | 64 % | 0 % | 24 % |
| R1_C_D2 | 44 | 0.22 | 0.120 | 0.187 | 1.2 | 91 % | 82 % | 0 % | 100 % |
| R2_A_D1 | 66 | 0.22 | 0.035 | 0.052 | 4.2 | 26 % | 20 % | 0 % | 100 % |
| R2_A_D2 | 31 | 0.11 | 0.006 | 0.007 | 14.6 | 0 % | 0 % | 0 % | 100 % |
| R2_B_D1 | 36 | 0.18 | 0.117 | 0.198 | 0.9 | 86 % | 75 % | 0 % | 100 % |
| R2_B_D2 | 42 | 0.12 | 0.034 | 0.054 | 2.2 | 29 % | 19 % | 0 % | 100 % |
| R2_C_D1 | 45 | 0.20 | 0.091 | 0.154 | 1.3 | 58 % | 56 % | 0 % | 82 % |
| R2_C_D2 | 16 | 0.11 | 0.017 | 0.028 | 3.7 | 44 % | 44 % | 0 % | 100 % |

- The excess is several-fold only where *Methanolobus* is nearly absent (R2_A_D1 4.2 ×, R2_A_D2 14.6 ×, R2_C_D2 3.7 ×); there the totals are 16–66 reads, so the per-sample folds are noisy (Poisson ± 12–25 %).
- Across these 12 samples, mcrA ranges 0.11–0.84 % of genomes; the *Methanolobus* MAG 0.007–0.42 %.
- Only the R1_A cores (where the *Methanolobus* MAG was assembled) place mcrA reads on kept-MAG contigs; elsewhere the *Methanolobus* mcrA reads are unassembled or on unbinned contigs, but still match the MAG at ≥ 95 % identity.

### mcrA taxonomy (`results/mcrA_taxonomy_unrestored_pooled.tsv`, `mcrA_taxonomy_by_sample.tsv`)

Genus of the best-hit reference (GTDB R232 for salterns MAG members, KEGG/UniProt organism otherwise), 12 unrestored-core samples pooled, 830 mcrA reads:

| genus | reads | share | best hit ≥ 95 % id | unassembled |
|---|---|---|---|---|
| *Methanolobus* | 479 | 57.7 % | 464 | 149 |
| *Methanohalophilus* | 143 | 17.2 % | 134 | 130 |
| *Methanogenium* | 54 | 6.5 % | 40 | 54 |
| *Methanococcoides* | 49 | 5.9 % | 35 | 49 |
| uncultured "methanogenic archaeon" (UniProt) | 24 | 2.9 % | 17 | 14 |
| *Ca.* Methanospirare (Methanophagales, ANME-1-related) | 13 | 1.6 % | 0 | 13 |
| Methanomassiliicoccales (*Methanomethylophilus*, *Methanomassiliicoccus*, *Methanogranum*, …) | 25 | 3.0 % | 2 | 25 |
| others (*Methanohalobium*, *Methanomethylovorans*, *Methanoplanus*, *Methanoculleus*, *Methanofollis*, *Methanothrix*, …) | 43 | 5.2 % | | |

- **Yes, methanogens other than *Methanolobus* are present in the unbinned reads**: 42 % of the mcrA reads (351), mostly other methylotrophic / halophilic Methanosarcinaceae (*Methanohalophilus*, *Methanococcoides*) and *Methanogenium*. 93 % of them (326) are unassembled; 69 % hit a reference at ≥ 95 % identity, i.e. they are close relatives of cultured genera, not novel lineages.
- In the other 12 sediments (153 reads) *Methanolobus* is 40 %, ANME-1-related *Ca.* Methanospirare 18 %, *Methanococcoides* 12 %; in the 15 waters (63 reads, all unassembled) *Methanolobus* 54 %, *Methanohalophilus* 32 %.
- C–P lyase (phnJ, 19,445 reads, 12 cores): 87 % bacterial, led by *Roseovarius* (37 %) and *Marivita* (14 %); 13 % archaeal (Halobacteria: *Haloarcula* 3.5 %, *Natrinema*, *Halonotius*, *Natronomonas*) (`CP_lyase_taxonomy_unrestored_pooled_phnJ.tsv`, `_phnGHIJKLM.tsv`).

### Binned / unbinned / unassembled split per gene (`results/split_per_gene_by_habitat.tsv`, per sample in `binned_unbinned_split.tsv`)

Share of the classified reads on kept-MAG contigs / other DAS Tool bins / unbinned contigs ≥ 1 kb / unassembled (%):

| gene | 12 unrestored cores | all 24 sediments | 15 waters |
|---|---|---|---|
| markers (8 families) | 34 / 9 / 15 / 43 | 25 / 7 / 13 / 56 | 36 / 22 / 16 / 27 |
| mcrA | 27 / 6 / 10 / 57 | 23 / 10 / 13 / 54 | 0 / 0 / 0 / 100 |
| phnJ | 29 / 9 / 20 / 42 | 26 / 8 / 18 / 48 | 26 / 12 / 30 / 32 |
| grdI | 42 / 10 / 12 / 37 | 41 / 9 / 12 / 38 | 12 / 10 / 6 / 72 |
| cutC | 40 / 4 / 13 / 44 | 18 / 3 / 9 / 70 | 2 / 2 / 4 / 92 |
| cntA | 30 / 9 / 19 / 42 | 26 / 5 / 14 / 55 | 13 / 7 / 19 / 62 |
| mtbB | 32 / 9 / 8 / 51 | 5 / 6 / 17 / 72 | 0 / 0 / 5 / 95 |
| mtmB | 12 / 3 / 19 / 65 | 3 / 2 / 12 / 82 | 1 / 0 / 2 / 97 |
| MttB family (mttB) | 31 / 6 / 16 / 47 | 18 / 4 / 13 / 65 | 31 / 5 / 7 / 57 |

So in the unrestored cores, 37–57 % of each guild gene's reads are not in any contig ≥ 1 kb, about the same as for the universal markers (43 %); methylotrophic methyltransferases (mtmB, mtbB) and mcrA are somewhat less assembled than the community average, grdI somewhat better.

## Problems and caveats

- **Read-level confusions** (validation): grdH vs grdB (grdH not used; grdI is), mttB vs mtgB (reported only merged), pmoA vs other CuMMOs (merged), mtbA vs mtaA/HemE (mtbA recall 0.69). Twilight homologs (below the KOfam threshold) of mtaB, mtbB, mtbA, phnH and phnI are counted by the read survey at ≥ 70 % identity; the read side is therefore slightly broader than KOfam for those genes.
- **phnM and phnL** read counts are 2–5 × phnJ (also in the MAGs: KOfam K06162 is found in 126 rollout MAGs vs K06163 in 78); use phnJ for the C–P lyase guild.
- **MAG-side units**: CoverM counts reads at ≥ 95 % identity to the 585 representatives only; strain variants below 95 % fall to "unbinned" on the MAG side but are counted on the read side, which inflates the fold a little. Genome sizes are MAG assembly sizes (consistent with CoverM's reads per MAG).
- **Small counts** for mcrA (16–216 reads per core sample): per-sample folds carry ± 7–25 % Poisson error.
- **mpnS** is not validated (15 references; 90 % of the 52 reads best-hit dinoflagellate MpnS sequences).
- **Thread budget**: three short checks (a 300k-line classification smoke test, a 300-contig bowtie2 test, and the summary scripts on the first 17 samples, each one thread for ≤ 2 min) ran while the 8 DIAMOND jobs held the 64-thread budget, so the track briefly used 65 threads.
- **A local process kill**: early on, a local `pkill -f "02_fetch_refs.py meta/ko_panel"` matched its own command line and killed that shell together with the first reference-fetch process (both mine, on this workstation; nothing on poplar and nothing of anyone else). The fetch was restarted, with cached files reused. Afterwards only PIDs were used.
- The KEGG fetch ran as two instances on the marker KOs (to save time); the KEGG gene lists are therefore re-derived from the fetched FASTA files and br08601 (`04_build_db.py ext`), not from `kegg_genes.tsv`.
