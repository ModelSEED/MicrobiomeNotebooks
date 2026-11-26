# Read survey v2: the methylated-compound / methanogenesis KO panel (138 KOs) in the 39 salterns metagenomes

*Track `read_markers_v2`, started on poplar on 2026-11-25 (all times CST). v1 (`../read_markers/`, poplar `read_markers/`) is only read, never written. Nothing is committed and nothing is written to KBase. Poplar work dir: `/scratch1/afreiburger/salterns_2026/read_markers_v2/`; this folder holds the scripts, the small metadata and the results.*

**Status: DONE** (2026-11-26 2:30 AM CST). The genus profile (section 1) was delivered at 3:40 PM CST on 2026-11-25. The v2 DB, validation, the full re-run of all 2,084,152,888 reads, the summaries under the output contract and the v1 → v2 regression check (section 6) are complete. Total wall time: 2026-11-25 3:10 PM – 2026-11-26 2:25 AM CST.

Tools: DIAMOND 2.2.8, HMMER 3.4, MMseqs2 18.8cc5c, Prodigal 2.6.3, seqkit 2.13.0, bowtie2 2.5.5, samtools 1.25, Python 3.11.15 / pandas 3.0.3 / numpy 2.4.5 (poplar env `emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_`); the workstation steps used `~/Documents/py_venv`. Thread budget: this track never ran more than 96 threads on poplar (all jobs at nice 10, detached in tmux, `TMPDIR` on /scratch1). Only my own PIDs were ever killed (a leftover `find`, a superseded fetch, its shell).

## 1. Read-based genus profile from the ribosomal-protein markers (`results/marker_genus_profile.tsv`)

Computed from the **v1 per-read tables** (`read_markers/results/per_sample/<s>.reads.tsv.gz`, the 2,177,966 reads that v1 classified to the 8 marker families) with `scripts/11a_firstpass_markers.sh`, `11b_kegg_org_assembly.py`, `11b2_map_org_gtdb.py` and `11c_marker_genus_profile.py`. It needed no new DIAMOND run. 11a ran 3:10–3:28 PM CST (24 single-CPU jobs, nice 10) and 11c 3:28–3:35 PM CST (1 CPU + an 8-thread DIAMOND blastp).

**Columns** (one row per sample × genus, all 39 samples, 56,625 rows):

| column | meaning |
|---|---|
| `sample` | |
| `genus_harmonized` | GTDB R232 genus; GTDB suffixes kept (`Clostridium_F`); Halomonas and Vreelandella kept separate; `unassigned` = marker reads without a genus (see below) |
| `genus_source_name` | the name(s) as they came from the references, `;`-separated, most reads first. For salterns MAGs it is the GTDB name; for KEGG genes it is the KEGG organism genus (NCBI naming). For example, `Halomonas;Billgrantia` in the Halomonas row means that some KEGG organisms named *Billgrantia* are GTDB R232 *Halomonas*. |
| `marker_reads` | marker reads of the sample whose best hit carries this genus |
| `share_of_marker_reads` | `marker_reads` / all marker reads of the sample. The denominator includes the `unassigned` row, so the shares sum to 1. |
| `copies_per_GE` | **depth-based genomes per genome equivalent**, as v1 `10_summarize.py`: for each marker family *f*, depth(g, f) = reads(g, f) × r_s / L_eff(f), with L_eff = 3 × median aa length + r_s − 180 + 1 and r_s = the sample's mean read length; copies_per_GE(g) = mean over the 8 families of depth(g, f) / the sample's GE depth (v1 `results/genome_equivalents.tsv`, `GE_depth_mean`). It sums to 1 over all rows of a sample, `unassigned` included. It differs from the share only through the families' different lengths (≤ 10 % relative). |
| extra columns (appended) | `share_of_assigned_marker_reads` (the denominator excludes `unassigned`); `reads_best_hit_ge95` (reads whose best hit is ≥ 95 % identical: genus-level confident); `reads_tier_salterns_MAG` / `reads_tier_kegg` / `reads_tier_transferred` (how the genus was obtained, below); `n_marker_families_with_reads` |

**Genus of a marker read** = genus of its best-hit v1 DB representative (99 % cluster), resolved in this order:
1. **salterns_MAG**: majority GTDB R232 genus of the cluster's binned salterns members (rollout kept MAGs + co-assembly MAGs), the v1 `10b` rule. If those MAGs have no GTDB genus (novel genus; 438 marker clusters), the read counts as `unassigned`.
2. **kegg**: a cluster with KEGG GENES members but no binned salterns member gets the majority **GTDB R232** genus of the members' genomes. Each KEGG organism was mapped to its genome assembly (KEGG REST `get/gn:`, `meta/kegg_org_assembly.tsv`) and then to the full R232 taxonomy (`bac120_taxonomy_r232.tsv.gz` + `ar53_taxonomy_r232.tsv.gz` from data.gtdb.ecogenomic.org, MD5 checked; accession matched without version, GCA↔GCF). This mapped 10,332 of the 11,829 KEGG organisms (`meta/kegg_org_gtdb.tsv`). The remaining 1,497 keep their KEGG organism genus.
3. **transferred**: a salterns-only cluster without any binned member (v1 `genus_source = none`; 7,842 marker clusters, which hold 25–58 % of the marker reads in v1) gets the genus of its best DIAMOND blastp hit (`--sensitive`) among same-family marker proteins that carry a genus. Those proteins are the binned rollout/co-assembly members, the KEGG members, and the marker proteins of the **1,324 first_pass MAGs** (Prodigal + the 10 marker KOfam profiles, top marker profile ≥ 0.5 × threshold, GTDB R232 classification from `mags_three_way.tsv`; 9,160 proteins). The genus is taken **only at ≥ 95 % identity**, because ribosomal proteins of neighbouring genera are often 85–95 % identical. Result: 2,024 clusters via rollout/co-assembly MAG proteins, 703 via first_pass MAGs, 415 via KEGG, and 4,700 stay `unassigned`.
v1 did not run this transfer for the marker families (`05b` transferred only for target classes), so without it the waters would be 50–60 % unassigned.

**Coverage.** `unassigned` is 4–44 % of the marker reads per sample (median 26 %), highest in the SF2 waters (41–44 %) and R2_restored_D_shore (40 %).

**Name notes.**
- **Halomonas / Vreelandella**: these are separate GTDB R232 genera, and both occur among the salterns MAGs. KEGG organisms named *Halomonas* map to GTDB *Vreelandella* (21 of 38), *Halomonas* (14) and others. `genus_source_name` shows the original names.
- **Saligentibacter** occurs in no reference (GTDB, KEGG or UniProt); it is presumably a spelling of *Salegentibacter*, which is a GTDB genus of rollout MAGs here (species cluster `rollout|R1_A_D1__concoct.bin.034`).
- **Marinilabilia** exists only as one first_pass MAG (`first_pass|metabat2.bin.296.fasta__R2_sediment`). It is not in the rollout or co-assembly sets, and KEGG has no *Marinilabilia* genome among the marker references. Its reads therefore come entirely from the first_pass transfer (73 % of them at ≥ 95 % best-hit identity).
- **Idiomarina** has no salterns MAG. Its reads come from KEGG references only.
- **Clostridium**: only KEGG references, split over many GTDB `Clostridium_*` genera. Only 16 % of these reads have a best hit ≥ 95 %, so in these samples the reads are mostly other Clostridia that land on *Clostridium* references. Treat the rows as "Clostridium-like".

**The six unrestored cores** (`results/marker_genus_profile_unrestored_focus.tsv`; copies per GE in %, mean of D1 and D2):

| genus | R1_A | R1_B | R1_C | R2_A | R2_B | R2_C | reads (12 samples) | best hit ≥ 95 % |
|---|---|---|---|---|---|---|---|---|
| Halomonas + Vreelandella | 0.51 | 0.22 | 0.32 | 0.06 | 0.15 | 0.17 | 1,545 | 91 % |
| Salegentibacter | 0.67 | 0.25 | 0.34 | 0.04 | 0.17 | 0.13 | 1,812 | 98 % |
| Idiomarina | 0.011 | 0.002 | 0.004 | 0.000 | 0.005 | 0.005 | 33 | 73 % |
| Marinilabilia | 0.03 | 0.05 | 0.04 | 0.02 | 0.07 | 0.16 | 366 | 73 % |
| Clostridium (all GTDB suffixes) | 0.016 | 0.024 | 0.018 | 0.021 | 0.035 | 0.040 | 173 | 16 % |

Other files: `results/marker_genus_profile_name_map.tsv` (harmonized ↔ source names, reads, tiers), `refs_small/marker_rep_genus.tsv.gz` (per marker representative: tier, genus, transfer hit and identity).

## 2. v2 panel and reference DB

### Panel (`meta/ko_panel_v2.tsv`, `scripts/00_make_panel_v2.py`)

286 KOs, of which 140 are counted targets and 10 are marker KOs:
- the 62 v1 KOs, unchanged (incl. `kegg_cap`). Four of them change role from v1 decoy to v2 counted target: K00499 CMO, K00479 gbcA, K21832 gbcB, K12905 hepD.
- the 108 + 30 KOs of `community_ko/ko_panel_methylated.tsv` (the FigS4 methanogenesis-pathway KOs added 2026-11-25 included). Every one of the 138 is a counted label; 110 of them were not in v1.
- 114 new decoy KOs: close paralogs of the new targets, chosen from KEGG `list/ko` (2026-11-25) and grouped by paralog family. For example, Mo-bisPGD enzymes vs torA/torZ/dmsA/ddhA/fwdB; MDR / Fe-ADH alcohol dehydrogenases vs frmA/fdhA/mdo/mdh/gbsB; ALDHs vs betB; DAAO-family oxidases (incl. PIPOX K00306, a second sarcosine oxidase) vs soxB/dmgo/soxA (monomeric); GcvT/YgfZ vs dmdA/tdm/dmgdh/sardh/soxA; glnA vs gmaS; gltB/gltD/gltS vs mgsABC/fwdC; metH/mtsB vs the corrinoid proteins and bhmt; UROD-type methyltransferases (mtgA, mtqA, mtvA, mtoA) vs mtaA/mtbA/mtsA; acsE vs mtrH; FMO/BVMO vs tmm; OYE vs tmd; AMP-binding ligases vs acs/dmdB; ACADs vs dmdC; crotonases vs dmdD; buk/tdcD vs ackA; ptb vs pta/eutD; cooS vs cdhA; acsB vs cdhC; 2-hydroxyacid dehydrogenases vs fdh; purN/fmt vs purU; mmuM vs bhmt; exaA/qheDH vs mxaF/xoxF; M24 peptidases vs dddP; the AcdA family vs dddX; CoA-transferase family III vs dddD; gsmt/bsmB vs gsdmt/sdmt; other phosphonate dioxygenases vs phnY/phnZ; fae-hps vs fae; LLM monooxygenases vs dmoA/ssuD/sfnG/mer; 4Fe-4S ferredoxins vs fwdH/fwdF/fwdG; pstS/phnS vs phnD. The full list with the gene each guards is in the panel's `note` column.

**Labels (output contract of 2026-11-25).** Each label is the panel `gene` string, unique per KO. A repeated symbol becomes `<gene>_<KO>` (`pta_K13788`). v1 labels are unchanged. Exception: the panel's explicit placeholders (`ddd (K28067)` … `ddd (K28073)`, "gene symbol to confirm from KEGG", and `mpnH/phnY*`) get the KEGG symbols confirmed on 2026-11-25 (table below). Labels containing spaces or slashes (`soxA (monomeric)`, `ssuD/msuD`, `pta (K15024)`) keep them; file names use `label_safe`.

**KEGG symbols confirmed (KEGG REST `get/ko:`, 2026-11-25; `meta/ko_symbols.tsv`):**

| KO | KEGG symbol | KEGG name | COG / note |
|---|---|---|---|
| K28067 | dddP | dimethylpropiothetin dethiomethylase [EC:4.4.1.3] | COG0006 (M24 peptidase family) |
| K28068 | dddW | dimethylpropiothetin dethiomethylase [EC:4.4.1.3] | COG1917 (cupin) |
| K28069 | dddY | dimethylpropiothetin dethiomethylase [EC:4.4.1.3] | |
| K28072 | dddX | DMSP lyase | COG1042 (acyl-CoA synthetase, NDP-forming) |
| K28073 | dddD | DMSP lyase | COG1804 (CaiB/BaiF CoA-transferase family III) |
| K23996 | phnY | methylphosphonate hydroxylase [EC:1.14.11.71] | the panel's `mpnH/phnY*`; label `phnY` |
| K28979 | tdm | trimethylamine-oxide aldolase [EC:4.1.2.32] | COG0404 (GcvT) + COG3665 |
| K24071 | gsmt-sdmt | glycine/sarcosine/dimethylglycine N-methyltransferase [EC:2.1.1.156 2.1.1.157] | label `gsdmt` (panel symbol) |
| K18897 | sdmt | sarcosine/dimethylglycine N-methyltransferase [EC:2.1.1.157] | |

**KOfam coverage.** All 286 KOs have a KOfam 2026-10-26 profile. Five counted KOs have no KOfam threshold: mpnS K18049 and hepD K12905 keep the v1 rule (top-scoring panel profile, full E ≤ 1e-30, HMM coverage ≥ 0.5); dmoA K16967, phnY K23996 and phnZ K05915 use the same rule plus a minimum full score derived from the KEGG references (`meta/synthetic_thresholds.tsv`, see below).

### External references (`scripts/02_fetch_refs_v2.py`, 2026-11-25 3:16–~4:35 PM CST)

- **v1 KOs:** v1's raw KEGG FASTA and link files and v1's UniProt sets are copied unchanged, so every v1 reference is byte-identical. K12905 has no KEGG gene in v1 either.
- **New KOs:** KEGG GENES via REST, as v1. Targets keep **every archaeal gene + 1,500 random bacterial genes**, the v1 marker rule, so that archaeal homologs of the shared methanogenesis/C1 genes are all present. New decoys get **1,000 random prokaryote genes**. Seeds are per KO (`random.Random(int(KO[1:]))`); v1 used one shared `Random(1)` whose draws depend on the panel order. *Deviation:* a first attempt with v1's cap of 3,000 would have fetched 265k genes (~3 h at KEGG's rate) and roughly doubled the DB, so the caps were lowered before any capped KO was fetched.
- **Eukaryote genes** of every counted target KO (≤ 300 random, per-KO seed): 5,213 genes, class `<gene>_euk` (decoy_reported). v1 had this only for the markers (`euk_marker`). The 15 waters are Dunaliella-rich, so algal homologs of generic genes (glyA, folD, fhs, betB, frmA, acs, …) would otherwise be counted as prokaryotic copies.
- **UniProtKB reviewed (Swiss-Prot) members** of the Pfam families of the generic targets (Bacteria + Archaea; folD PF00763/PF02882, purU PF00551, fhs PF01268, glyA PF00464, sdaA PF03313/PF03315, fdh PF00389/PF02826, amtB PF00909, fghA PF00756, frmA/fdhA/mdo PF08240/PF00107, ssuD PF00296, ackA PF00871, pta/eutD PF01515, acs PF00501). Members that pass no panel KO and whose top profile is not the family's target become `<family>_other` decoys, as v1's `SDIMO_other` / `GRE_other` / `Rieske_other`.
- To keep each host under KEGG's request rate, the fetch ran as two instances with the same per-KO seeds: poplar fetched the targets, eukaryote sets and UniProt, and this workstation fetched the new decoys, which were then copied to poplar.

The fetch ran 3:16–4:38 PM CST. Result: 228,214 KEGG prokaryote genes (84,414 of them v1 copies) + 5,213 eukaryote genes of 20 target KOs, and 17 reviewed UniProt sets (12–806 entries each). Every requested gene was received (`fetch_log.tsv`). `04_build_db_v2.py ext` (4:40 PM CST) wrote 285,164 external proteins.

### HMM search (`scripts/03_panel_hmmsearch_v2.sh`)

HMMER 3.4 with the v1 flags (`--cpu 1 --noali -E 1e-5 --domE 1e-5 -Z 20000000`), nice 10.
- **Salterns** (3:30–5:14 PM CST): the 224 profiles not in v1, split into 4 chunks, against the 39 per-sample gene sets and the 520 co-assembly MAGs (v1 `genes/`, 13.6 M proteins). 160 single-CPU jobs, 80 at a time; the four co-assembly chunks took up to 104 min (poplar load 600–1,200 on 336 cores from other users). Labelling uses v1's 62-profile tables together with these. With the fixed `-Z`, each profile's hits do not depend on the other profiles searched, so the union equals a fresh 286-profile search. All 320 output tables end with `[ok]`.
- **External set** (4:40–4:57 PM CST): all 286 profiles, 16 single-CPU jobs.

### KOs without a KOfam threshold (`scripts/04b_synthetic_thresholds.py`, `meta/synthetic_thresholds.tsv`)

For dmoA K16967, phnY K23996 and phnZ K05915, the threshold (full score) = 1 + the highest score that any KEGG gene of another panel KO reaches on that profile, i.e. KOfam-like separation from the nearest paralogs. This threshold is used exactly like a KOfam threshold in the best-KO rule.

| KO | references of the KO itself (full score) | closest other panel genes | threshold |
|---|---|---|---|
| dmoA K16967 | KEGG links K16967 to no genome gene (only GenBank ADU77278). The one characterised DmoA is Swiss-Prot **E9JFX9** (*Hyphomicrobium sulfonivorans*) in the reviewed LLM set: 1,169.6 | LadA K20938 ≤ 850.7 | 851.7 |
| phnY K23996 | 2 KEGG genes (*Gimesia maris*, *Planctomycetes* Pan161): 624.9–625.0 | phnY_AEP K21195 ≤ 500.7 | 501.7 |
| phnZ K05915 | the same 2 genomes: 478.1–478.6 | tmpB K24290 ≤ 147.2 | 148.2 |

### Labelling and clustering (`04_build_db_v2.py label`, 5:26–5:37 PM CST; `05_cluster_db_v2.sh`, 5:43–6:09 PM CST)

- 613,088 labelled proteins: 347,6xx salterns, 233,427 KEGG, 31,733 UniProt (v1: ~261k). 204k of them are in `_like` / `_other` decoy classes. The largest are phnK_like 37k and phnL_like 22k (both as in v1), sdmt_like 19k (the K18897 profile reaches generic class-I SAM methyltransferases), dmdD_like 10k (crotonases), frmA_like 9k (MDR) and dmdC_like 8k (acyl-CoA dehydrogenases).
- K10944 resolution: 1 salterns protein (→ pmoA). **Eukaryote resolution** (new): of the 10,449 salterns target proteins that are not on a kept or co-assembly MAG, 95 best-hit a eukaryote reference and became `<gene>_euk` (amtB 35, betA 28, glyA 17, sardh 4, folD 3, frmA 2, fghA 2, CMO 2, MTO 1, acs 1).
- Clustering as in v1 (MMseqs2 easy-cluster, coverage 0.8; 99 % for targets / markers / reported classes, 95 % for decoys, 90 % for the hold-out units). It ran 8 labels at a time × 4 threads and was concatenated in label order. 431 labels, **338,105 representatives** (v1: 137,012).
- **v1 ↔ v2 identity check** (`05c_check_markers_v1_v2.py`, `results/check_labels_v1_v2.tsv`): the 8 marker families and `euk_marker` have identical members and identical 99 % clusters in v1 and v2, and so does every v1 counted target class (mcrA/B/G, mtaA/B, mtbA/B, mtmB, mttB, mtgB, mttB_family_other, phnG–M, grdH/I, grdB, cutC/D, cntA/B, mpnS, hepD, pmoA, mmoX, amoA, cuMMO_other). Differences are only in decoy classes: `Rieske_other` 4,953 → 279 members and `cntB_like` 3,954 → 68, `cutD_like` 3,335 → 314, `phnM_like` 9,490 → 8,368, because their members now pass a new panel KO or are below-threshold homologs of a new target. gbcA/gbcB have the same members but are clustered at 99 % instead of 95 % as counted targets. CMO was the v1 decoy `cmo`. Any v1 → v2 change in a v1 label's counts can therefore only come from the new competitors in the DB (new targets and decoys), not from changed v1 references.

*Why a threshold and not the v1 mpnS/hepD rule:* the first labelling pass applied the v1 rule (top panel profile, E ≤ 1e-30, coverage ≥ 0.5, plus the floor) and gave dmoA zero proteins. E9JFX9 also passes the KOfam threshold of LadA (795.8 ≥ the LadA threshold), and the best-KO rule takes a passing KO before any no-threshold rule. As a thresholded KO it is called dmoA (E = 0 vs 6e-236). mpnS and hepD keep the v1 rule unchanged.

## 3. Validation (`scripts/07_simulate.py`, `run_validation.sh`, `08_classify.py validate`, `08b`, `08c`-style grid summary, `08d_identity_cuts.py`; 6:10–6:48 PM CST on poplar, 32 threads; final tables recomputed on this workstation)

**Design = v1** (see `../read_markers/README.md`). Hold-out unit = a 90 % cluster of one label; 20 % of each label's 90 % clusters (≥ 1 when ≥ 3 clusters, ≤ 60) are held out (`random.Random(1)`), and the validation DB is the production DB without every representative of a held-out cluster. 148-nt reads with 0.3 % substitutions, ~5× coverage. v2 numbers: **14,205 held-out genes simulated** (4,865 salterns genes with ±150 nt of contig, 9,340 KEGG genes from KEGG `ntseq`), 133,414 background reads from 4,000 random salterns genes without a hit to any of the 286 profiles, and 1,035 KOfam-twilight homologs (0.5–1 × threshold, not in the DB): **793,002 simulated reads**.

**Top-k test.** The validation reads were run with the production flags at `-k 5` and at `-k 25`. Mean target recall was 0.9192 vs 0.9196 and precision 0.9548 vs 0.9544; only gsdmt and grdH moved by more than 1 point, and both are in confused pairs. **Production keeps v1's `-k 5`**, so no DIAMOND flag deviates from v1. On real reads, 0.58 % of reads with hits had all 5 hits tied at the top bitscore (v1 DIAMOND table of R1_A_D1).

**Global cut-offs** (`validation/grid_summary.tsv`, EMAX 1e-5, LMIN 30; the grid was reduced to MARGIN 0/2 × global PMIN 0/50/60/70, because v1 showed EMAX and LMIN to be insensitive). As in v1, margin 0 and no global identity cut give the best mean F1 over the targets (0.940); margin 2 or global identity cuts cost more recall than they gain in precision. The v1 cut-offs are kept.

**Per-label identity cuts** (v1 rule: only where an identity cut removes identified false positives at a few % recall; `validation/identity_cuts.tsv`). Added in v2 (`cutoffs.tsv`):

| label | PMIN | effect in the hold-out test |
|---|---|---|
| mtmC | 70 % | mttC reads called mtmC: precision 0.762 → 0.982, recall 1.000 → 0.991 |
| mtbA (v1 label) | 60 % | mtbA_like / mtaA_like reads: precision 0.898 → 0.988, recall 0.688 → 0.664 |
| amtB | 60 % | background / amtB_like / amtB_euk reads at 50–60 %: precision 0.909 → 0.963, recall 0.961 → 0.924 |
| glyA | 60 % | precision 0.950 → 0.977, recall 0.981 → 0.963 |
| sdaA | 50 % | background false positives 37 → 22, recall 0.973 → 0.951 |
| dmgo_euk | 60 % | background false positives 17 → 4, recall unchanged |

These cuts did not help, so none was set: cntB vs dmmB, gsdmt vs sdmt / gsmt, mttC vs mtbC, eutD vs pta, fhs (false positives from fhs_like / fhs_euk), mgsC.

**Background.** At the final cut-offs, 90 of the 133,414 background reads are assigned to a counted label (v1: 62 of 138,035). Before the four v2 cuts it was 206, and 180 of those lay inside a background CDS, at 35–59 % identity, mostly to amtB, glyA, sdaA, mgsC and fhs (`validation/background_fp_final.tsv`). For the generic genes, read counts can therefore include a few reads of distant homologs without a panel-HMM hit (order 10⁻⁴ of all non-target reads).

**Merged reporting classes** (`meta/merged_v2.json`; individual labels are still reported; merged rows in `per_sample_KO.tsv`; `validation/validation_merged_classes_v2.tsv`):

| merged class | recall | precision | why |
|---|---|---|---|
| MttB_superfamily(mttB+mtgB+other) | 0.951 | 0.957 | v1 merge: 90 % of mtgB reads are called mttB |
| GrdB_family(grdH+grdB) | 0.913 | 0.862 | v1 merge: 39 % of grdH reads are called grdB; grdI stays the betaine-reductase proxy |
| CuMMO_all(pmoA+amoA+other) | 0.858 | 1.000 | v1 merge |
| Methylamine_corrinoid_proteins(mtbC+mtmC+mttC) | 0.924 | 0.978 | mttC precision 0.836 (80 mtbC reads); mtmC needed a 70 % cut |
| GB_methylation(gsdmt+sdmt) | 0.873 | 0.952 | gsdmt recall 0.628 (20 % to the decoy gsmt K18896, 12 % to sdmt), sdmt precision 0.881 |
| DMGDH_SARDH(dmgdh+sardh) | 0.982 | 0.951 | sardh: 1 held-out gene, 58 % of its reads called dmgdh, 35 % ambiguous |
| Pta_family(pta+pta_K13788+pta (K15024)+eutD) | 0.976 | 0.980 | eutD recall 0.721 (21 % to pta) |
| MtaA_MtsA(mtaA+mtsA) | 0.905 | 0.964 | mtsA recall 0.597 (32 % to mtaA; 2 held-out genes) |
| CntB_DmmB_reductases(cntB+dmmB) | 0.968 | 0.987 | cntB precision 0.884 (52 dmmB reads); 3 % of dmmB reads are called cntB |

**Flagged labels** (precision < 0.90 or recall < 0.80 with ≥ 50 reads; 16 of 128 validated targets): cntB (precision 0.884, from dmmB); eutD (recall 0.721, to pta); **fwdC (recall 0.611) and fwdF (0.788)**, lost to divergence (median identity of a held-out gene to its nearest remaining reference 77–78 %), not to other labels, so their precision is 0.99; grdH (0.485, to grdB; v1 same); gsdmt (0.628 / 0.860); msmB (recall 0.642; 5 held-out genes at a median 47 % identity to the remaining references); mtbA (recall 0.664, 13 % to hemE; v1 0.689); mtbB (precision 0.853, from mtbB_like and flank reads; v1 0.906); mtgB / mttB (merged, as v1); mtsA (recall 0.597, to mtaA); mttC (precision 0.836); pmoA (merged, as v1); sardh (merged); sdmt (precision 0.881). **The three 4Fe-4S ferredoxin KOs are not precision-flagged:** fwdH 0.983 / 0.983, fwdG 0.864 / 0.987, fwdF 0.788 / 0.993 (recall / precision; the `fwdF_like` 5.7k and `fwdH_like` 4.2k decoy classes absorb the generic ferredoxins). In real reads, though, generic ferredoxins are far more abundant than in the simulated mix, so treat them with caution (see the archaeal share in the results).

**Not validated** (no held-out gene because the label has < 3 clusters at 90 % identity, or < 50 simulated reads): AOX, CMO, dddW (3 genes, 43 reads: 0.860 / 0.949), dddX, ddhB (1 gene, 0.568), ddhC (1 gene, 0.962), dmoA (1 reference), hepD, mdh, mpnS (1 gene, 47 reads, recall 0.04), mttB_family_other, phnY (2 KEGG references), phnZ (1 gene, 0.478).

**Twilight homologs** (KOfam 0.5–1 × threshold, not in the DB): as in v1, reads of such genes are counted for some labels at moderate identity: pta 958, fhs 992, folD 509, mtrH 475, fwdD 336, mttC 353, mtbB 288, fwdE 283 of the simulated twilight reads. The read survey is therefore somewhat broader than KOfam for these genes.

**Per-label results at the final cut-offs** (`validation/validation_summary_v2.tsv`; all labels incl. decoy_reported and markers; confusion in `validation_confusion_final.tsv`, recall by divergence in `validation_recall_by_divergence_final.tsv`). Markers: recall 0.925–0.985, precision ≥ 0.962. "(v1)" = label also in v1.

| label | KO | held-out genes | reads | recall | precision | note |
|---|---|---|---|---|---|---|
| ackA | K00925 | 40 | 1608 | 0.963 | 0.972 |  |
| acs | K01895 | 53 | 3184 | 0.970 | 0.941 |  |
| amtB | K03320 | 34 | 1598 | 0.924 | 0.964 |  |
| betA | K00108 | 58 | 3272 | 0.967 | 0.929 |  |
| betB | K00130 | 60 | 2959 | 0.904 | 0.968 |  |
| bhmt | K00544 | 19 | 679 | 0.985 | 0.996 |  |
| cdhA | K00192 | 15 | 1145 | 0.988 | 0.998 |  |
| cdhB | K00195 | 17 | 284 | 0.986 | 0.993 |  |
| cdhC | K00193 | 16 | 782 | 0.899 | 0.923 |  |
| cdhD | K00194 | 39 | 1502 | 0.860 | 0.956 |  |
| cdhE | K00197 | 42 | 1905 | 0.984 | 0.991 |  |
| cntA (v1) | K22443 | 23 | 858 | 0.977 | 0.982 |  |
| cntB (v1) | K22444 | 15 | 473 | 0.867 | 0.884 | precision<0.90; FP from: dmmB 52; vanB 2 |
| cutC (v1) | K20038 | 25 | 2043 | 0.968 | 0.998 |  |
| cutD (v1) | K20037 | 25 | 764 | 0.811 | 0.998 |  |
| dddD | K28073 | 10 | 837 | 0.952 | 0.999 |  |
| dddL | K16953 | 26 | 559 | 0.945 | 0.993 |  |
| dddP | K28067 | 13 | 582 | 1.000 | 0.968 |  |
| dddY | K28069 | 4 | 165 | 0.867 | 1.000 |  |
| ddhA | K16964 | 1 | 92 | 1.000 | 0.989 |  |
| dgcA | K21833 | 14 | 974 | 0.977 | 0.995 |  |
| dmdA | K17486 | 23 | 883 | 0.969 | 0.998 |  |
| dmdB | K20034 | 60 | 3308 | 0.984 | 0.985 |  |
| dmdC | K20035 | 46 | 2731 | 0.974 | 1.000 |  |
| dmdD | K20036 | 59 | 1533 | 0.989 | 0.999 |  |
| dmgdh | K00315 | 56 | 4598 | 0.989 | 0.941 |  |
| dmgo | K00309 | 50 | 4140 | 0.982 | 0.979 |  |
| dmmA | K22342 | 49 | 1659 | 0.995 | 0.997 |  |
| dmmB | K22343 | 51 | 1598 | 0.954 | 0.973 |  |
| dmmC | K22344 | 53 | 947 | 0.872 | 0.999 |  |
| dmsA | K07306 | 60 | 4755 | 0.918 | 0.968 |  |
| eutD | K04020 | 11 | 355 | 0.721 | 0.992 | recall<0.80; lost to: pta 21% |
| fae | K10713 | 59 | 1034 | 0.996 | 1.000 |  |
| fdh | K00122 | 33 | 1264 | 0.998 | 0.994 |  |
| fdhA | K00148 | 60 | 2318 | 0.997 | 0.997 |  |
| fghA | K01070 | 57 | 1587 | 0.994 | 0.983 |  |
| fhs | K01938 | 58 | 3224 | 0.992 | 0.956 |  |
| folD | K01491 | 52 | 1526 | 0.973 | 0.946 |  |
| frmA | K00121 | 60 | 2201 | 0.981 | 0.953 |  |
| ftr | K00672 | 59 | 1783 | 0.985 | 0.991 |  |
| fwdA | K00200 | 53 | 2919 | 0.982 | 0.998 |  |
| fwdB | K00201 | 55 | 2284 | 0.927 | 0.993 |  |
| fwdC | K00202 | 60 | 1637 | 0.611 | 0.990 | recall<0.80 |
| fwdD | K00203 | 32 | 429 | 0.986 | 0.991 |  |
| fwdE | K11261 | 50 | 1178 | 0.938 | 0.997 |  |
| fwdF | K00205 | 30 | 1078 | 0.788 | 0.993 | recall<0.80 |
| fwdG | K11260 | 19 | 177 | 0.864 | 0.987 |  |
| fwdH | K00204 | 8 | 120 | 0.983 | 0.983 |  |
| gbcA | K00479 | 31 | 1280 | 0.993 | 0.977 |  |
| gbcB | K21832 | 60 | 2302 | 0.954 | 0.996 |  |
| gbsB | K11440 | 6 | 239 | 0.908 | 0.904 |  |
| gfa | K03396 | 41 | 819 | 0.948 | 0.990 |  |
| glyA | K00600 | 38 | 1645 | 0.963 | 0.977 |  |
| gmaS | K01949 | 6 | 260 | 0.896 | 0.979 |  |
| grdH (v1) | K21579 | 10 | 425 | 0.485 | 0.907 | recall<0.80; lost to: grdB 39%; ambiguous 7%; grdH_like 5% |
| grdI (v1) | K21578 | 13 | 577 | 0.948 | 0.986 |  |
| gsdmt | K24071 | 13 | 702 | 0.628 | 0.860 | precision<0.90 recall<0.80; lost to: gsmt 20%; sdmt 12%; ambiguous 5%; FP from: sdmt 44; gsmt 23; none 5 |
| mauA | K15228 | 29 | 541 | 0.952 | 1.000 |  |
| mauB | K15229 | 32 | 1246 | 0.954 | 1.000 |  |
| mch | K01499 | 55 | 1823 | 0.988 | 0.997 |  |
| mcrA (v1) | K00399 | 14 | 783 | 0.997 | 1.000 |  |
| mcrB (v1) | K00401 | 19 | 828 | 0.998 | 1.000 |  |
| mcrG (v1) | K00402 | 16 | 405 | 0.993 | 1.000 |  |
| mddA | K21310 | 52 | 1362 | 0.974 | 0.999 |  |
| mdo | K17067 | 7 | 298 | 0.970 | 0.997 |  |
| mer | K00320 | 54 | 1839 | 0.956 | 0.973 |  |
| mgdA | K22084 | 53 | 2200 | 0.946 | 0.978 |  |
| mgdB | K22085 | 60 | 574 | 0.960 | 0.944 |  |
| mgdC | K22086 | 60 | 5794 | 0.951 | 0.967 |  |
| mgdD | K22087 | 59 | 1104 | 0.874 | 0.959 |  |
| mgsA | K22081 | 57 | 1705 | 0.963 | 0.990 |  |
| mgsB | K22082 | 55 | 1304 | 0.902 | 0.966 |  |
| mgsC | K22083 | 47 | 2332 | 0.958 | 0.963 |  |
| mmoX (v1) | K16157 | 1 | 52 | 1.000 | 1.000 |  |
| msmA | K16968 | 5 | 205 | 0.985 | 0.976 |  |
| msmB | K16969 | 5 | 81 | 0.642 | 1.000 | recall<0.80 |
| mtaA (v1) | K14080 | 28 | 942 | 0.903 | 0.937 |  |
| mtaB (v1) | K04480 | 21 | 965 | 0.996 | 1.000 |  |
| mtaC | K14081 | 14 | 364 | 0.984 | 0.997 |  |
| mtbA (v1) | K14082 | 7 | 256 | 0.664 | 0.988 | recall<0.80; lost to: hemE 13%; mtbA_like 5% |
| mtbB (v1) | K16178 | 9 | 357 | 0.958 | 0.853 | precision<0.90; FP from: mtbB_like 32; none 27 |
| mtbC | K16179 | 37 | 768 | 0.835 | 0.921 |  |
| mtgB (v1) | K25220 | 32 | 1529 | 0.033 | 0.360 | precision<0.90 recall<0.80; lost to: mttB 90%; ambiguous 5%; FP from: mttB 86; mtgB_like 3 |
| mtmB (v1) | K16176 | 18 | 740 | 0.932 | 0.986 |  |
| mtmC | K16177 | 5 | 109 | 0.991 | 0.982 |  |
| MTO | K17285 | 48 | 2263 | 0.974 | 0.990 |  |
| mtrA | K00577 | 40 | 905 | 0.901 | 0.995 |  |
| mtrB | K00578 | 16 | 163 | 0.932 | 1.000 |  |
| mtrC | K00579 | 18 | 487 | 0.977 | 0.998 |  |
| mtrD | K00580 | 14 | 347 | 0.986 | 1.000 |  |
| mtrE | K00581 | 13 | 385 | 0.982 | 1.000 |  |
| mtrF | K00582 | 13 | 92 | 0.989 | 1.000 |  |
| mtrG | K00583 | 12 | 89 | 0.921 | 0.988 |  |
| mtrH | K00584 | 20 | 624 | 0.885 | 0.974 |  |
| mtsA | K16954 | 2 | 72 | 0.597 | 0.977 | recall<0.80; lost to: mtaA 32% |
| mttB (v1) | K14083 | 25 | 1256 | 0.904 | 0.432 | precision<0.90; FP from: mtgB 1377; mtgB_like 60; none 29 |
| mttC | K14084 | 30 | 641 | 0.819 | 0.836 | precision<0.90; FP from: mtbC 80; mttC_like 17; mtmC_like 3 |
| mxaF | K14028 | 11 | 679 | 1.000 | 0.987 |  |
| mxaI | K14029 | 16 | 149 | 1.000 | 1.000 |  |
| phnD | K02044 | 53 | 1747 | 0.843 | 0.984 |  |
| phnG (v1) | K06166 | 60 | 919 | 0.970 | 0.989 |  |
| phnH (v1) | K06165 | 59 | 1152 | 0.900 | 0.981 |  |
| phnI (v1) | K06164 | 50 | 1847 | 0.980 | 0.996 |  |
| phnJ (v1) | K06163 | 20 | 589 | 0.995 | 0.923 |  |
| phnK (v1) | K05781 | 59 | 1604 | 0.968 | 0.992 |  |
| phnL (v1) | K05780 | 60 | 1413 | 0.987 | 0.985 |  |
| phnM (v1) | K06162 | 49 | 1895 | 0.978 | 0.992 |  |
| pmoA (v1) | K10944 | 2 | 53 | 0.604 | 0.711 | precision<0.90 recall<0.80; lost to: ambiguous 11%; amoA 11%; cuMMO_other 11%; FP from: cuMMO_other 13 |
| pta | K00625 | 60 | 2119 | 0.952 | 0.919 |  |
| pta (K15024) | K15024 | 59 | 1230 | 0.986 | 0.986 |  |
| pta_K13788 | K13788 | 54 | 3621 | 0.986 | 0.988 |  |
| purU | K01433 | 55 | 1630 | 0.979 | 0.997 |  |
| sardh | K00314 | 1 | 81 | 0.000 | 0.000 | precision<0.90 recall<0.80; lost to: dmgdh 58%; ambiguous 35%; dmgo_like 5%; FP from: none 1 |
| sdaA | K01752 | 43 | 1566 | 0.952 | 0.985 |  |
| sdmt | K18897 | 39 | 1147 | 0.908 | 0.881 | precision<0.90; FP from: gsdmt 87; none 25; gsmt 17 |
| sfnG | K17228 | 59 | 2117 | 0.997 | 1.000 |  |
| soxA | K00302 | 59 | 5649 | 0.950 | 0.967 |  |
| soxA (monomeric) | K00301 | 59 | 2236 | 0.936 | 0.998 |  |
| soxB | K00303 | 60 | 2474 | 0.983 | 0.969 |  |
| soxD | K00304 | 60 | 610 | 0.975 | 0.949 |  |
| soxG | K00305 | 59 | 1149 | 0.946 | 0.978 |  |
| ssuD/msuD | K04091 | 60 | 2202 | 0.954 | 0.998 |  |
| tdm | K28979 | 33 | 2548 | 0.957 | 0.998 |  |
| tmd | K00317 | 43 | 3017 | 0.980 | 0.997 |  |
| tmm | K18277 | 47 | 2138 | 0.999 | 1.000 |  |
| torA | K07811 | 19 | 1572 | 0.949 | 0.989 |  |
| torZ | K07812 | 51 | 4124 | 0.967 | 0.970 |  |
| xoxF | K23995 | 59 | 3591 | 0.983 | 0.998 |  |

## 4. Full read survey

**DIAMOND** (`scripts/06_diamond_reads.sh`): DIAMOND 2.2.8 `blastx --mid-sensitive -k 5 -e 1e-3 -b 4 -c 1`, the v1 flags unchanged, of every QC read (`rollout/reads/<s>/map_1.fq.gz`, `map_2.fq.gz`; each mate a separate query; 2,084,152,888 reads) against `db/survey.dmnd` (338,105 representatives; md5 in `poplar_copies/db_md5.txt`). 6 samples at a time × 16 threads, largest samples first, nice 10, `TMPDIR` and `-t` on /scratch1. **Ran 2026-11-25 6:58 PM – 2026-11-26 1:54 AM CST (6 h 56 min).** Per mate: up to 56:54 min; maximum RSS 40.5 GB per job; 19 GB of hit tables. Post-processing ran 1:36–2:02 AM CST: 15,183,990 counted reads (v1: 3.29 M), re-mapped to 486,803 panel-carrying contigs. Summaries ran 2:02–2:10 AM CST, and the lost-read attribution (`12 attribute`, 5 samples) until 2:16 AM CST. At 2:19 AM CST, 10 / 10c / 12 were re-run once, after adding explicit zero rows for AOX (no reference in the DB); no other number changed. `12b_lost_read_evalues.py` ran at 2:24 AM CST.
- *Benchmark* (`scripts/00_bench_v1_v2.sh`, 6:10–6:27 PM CST): the first 8 M reads of R1_A_D1 `map_1`, run against the v1 and v2 DBs at the same time, 16 threads each. v1 DB 12:07 min (310 % CPU), v2 DB 16:30 min (403 % CPU): 1.36× wall time, 1.8× CPU time. Reads with any hit: 3.2 % (v1) vs 8.2 % (v2).
- Poplar's load was 600–1,200 on 336 cores from other users during the run, so each 16-thread job got about 4 cores (the first mate, SF2_C_H2O mate 1, took 48 min at 415 % CPU).

**Post-processing** (`scripts/run_post.sh`): per sample, once both mates are through DIAMOND, `08_classify.py sample` at the final cut-offs, then `09_map_assigned.sh` (bowtie2 2.5.5 / samtools 1.25 re-mapping of the classified reads to the sample's contigs with a hit to any of the 286 profiles; ≥ 95 % identity). Thread budget: a post job (5 threads) starts only while 16 × (DIAMOND slots still needed, ≤ 6) + 5 × (running post jobs) + 5 ≤ 96.

**Summaries** (`scripts/run_summaries.sh`): `10_summarize.py`, `10b_split_and_taxonomy.py`, `10c_key_tables.py`, `10d_methanol_taxonomy.py`, `12_regression_v1_v2.py`, 1 CPU each. Before the run, the v2 versions were dry-tested on this workstation with v1's per-sample tables of R1_A_D1 and R2_A_D1 as input. They reproduce v1's `per_sample_KO.tsv` copies per GE and `binned_unbinned_split.tsv` exactly (maximum difference 0), and all contract files keep the v1 columns.

### Output contract (coordinator, 2026-11-25)

`results/` holds the v1 file names with the v1 columns, computed from v2. New columns are appended at the end; nothing is renamed or dropped:
- `per_sample_KO.tsv`: v1 columns + `ko` (KEGG KO of target / marker labels, `K02926,K02930` for uL4; empty for decoy_reported and merged rows), `reads_archaeal`, `copies_per_GE_archaeal` (reads whose best-hit representative is archaeal; copies_per_GE × that share), `reads_methanogen_lineage`, `copies_per_GE_methanogen_lineage`, `reads_bacterial`, `reads_eukaryotic`, `reads_domain_unknown`. `label` = the gene symbol (v1 labels unchanged; v2 labels as in section 2), plus the v1 and v2 merged classes.
- `read_taxonomy.tsv`: v1 columns + `order`, `archaeal_lineage`; now covers every counted non-marker label.
- `key_unrestored_cores.tsv`, `mcrA_taxonomy_unrestored_pooled.tsv`, `mcrA_taxonomy_by_sample.tsv`, `binned_unbinned_split.tsv`, `methanol_taxonomy.tsv`, `CP_lyase_taxonomy_unrestored_pooled_phnJ.tsv`, `split_per_gene_by_habitat.tsv`: v1 columns unchanged (`split_per_gene_by_habitat.tsv` has rows for every counted label).
- New: `methanogenesis_pathway_archaeal_share.tsv` (FigS4 and methyltransferase labels per sample: archaeal and methanogen-lineage share and copies per GE), `regression_v1_v2_*.tsv`, `marker_genus_profile*.tsv`.
- **Taxonomy names change in v2**: KEGG references now carry the GTDB R232 genus of their genome (v1 used the KEGG organism name), so genus columns use GTDB R232 throughout wherever the genome is in R232.
- Methanogen lineage = best-hit representative in Archaea whose order, family or genus name starts with *Methano* (or *Methermicoccus*). ANME lineages (order Methanophagales; families Methanoperedenaceae, Methanocomedenaceae, Methanogasteraceae; genera Methanoperedens, Methanocomedens, Methanogaster, Methanophaga, Methanospirare) are excluded and counted only as archaeal. The domain of salterns-only clusters comes from their best blastp hit (`transferred`, any identity, as v1 for the genus), which is reliable for domain only above ~50–60 % identity.

## 5. Results (all 39 samples; `results/`)

- **AOX (K17066) has no prokaryotic reference**: KEGG lists only eukaryote genes for K17066 (yeast methanol oxidases), which form `AOX_euk` (1,008 reads in all 39 samples), and no salterns protein passes the K17066 profile (131 below-threshold homologs form `AOX_like`). Its rows in `per_sample_KO.tsv` are explicit zeros (group tagged "no prokaryotic reference in the DB"), so every panel KO has a row.
- `results/per_sample_KO.tsv`: every counted label (140 targets + `mttB_family_other`, 8 markers, the `_euk` / v1 decoy_reported classes) and the 9 merged classes × 39 samples. Columns: reads, ambiguous / low-identity best hits, L_nt, L_eff, depth, RPKM, copies per GE (mean- and median-normalised), and the appended `ko`, archaeal and methanogen-lineage columns. `results/summary_per_label_v2.tsv` gives the per-label medians: 12 unrestored cores, 24 sediments, 15 waters, and the archaeal / methanogen-lineage / eukaryotic fractions.
- **Methanogenesis-pathway genes** (FigS4 KOs; per sample in `results/methanogenesis_pathway_archaeal_share.tsv`). In the 12 unrestored-core samples (median copies per GE), the methanogen-specific genes are at mtrB–G 0.16–0.20 %, mtrA 0.66 %, mtrH 0.99 %, cdhA–C 0.22–0.30 %, fwdF 0.70 %, fwdG 0.50 % and fwdH 0.01 % (59 reads in all 39 samples). 83–100 % of the mtrB–G, cdhA–C, fwdG and fwdH reads best-hit methanogen lineages, against 67 % for fwdF and 44–46 % for mtrA/mtrH. Their methanogen-lineage copies per GE are 0.0014–0.0048 (fwdH 0.0001), in line with mcrA (median 0.0021). The shared genes are dominated by non-methanogens: ackA, pta, pta_K13788, pta (K15024) and eutD (≤ 0.1 % methanogen reads); acs (32 % archaeal, 92 % of them Halobacteriales; 1.8 % methanogen); mch (77 % archaeal, 97 % Halobacteriales; 2.1 % methanogen); mer (64 % archaeal, 98 % Halobacteriales; 1.3 % methanogen); fwdE (15 % archaeal; 2 % methanogen); fwdA–D (11–39 % archaeal, nearly all methanogen); ftr (17 % archaeal, 84 % of them Methanosarcinales); cdhD/cdhE (3 % / 11 % archaeal; the bacterial Wood–Ljungdahl acsD/acsC are in the same KOs). **For these, use `copies_per_GE_methanogen_lineage` (or `_archaeal`), not `copies_per_GE`.**
- The most abundant new targets in the 12 cores (copies per GE, median) are amtB 1.44, glyA 1.13, acs 1.08, phnD 1.08, folD 0.89, sdaA 0.80, ackA 0.73, fhs 0.60, mgsC 0.59, dmdC 0.57, dmgdh 0.54 and betA 0.49. Values > 1 are expected for multi-copy generic genes.
- The guild key tables (`key_unrestored_cores*.tsv`, `per_guild_read_vs_MAG.tsv`) are essentially unchanged from v1 (medians over the 12 cores, v1 → v2): methanogen fold read/MAG 1.81 → 1.81; methylamine users 8.64 → 8.50; MePn consumers 2.05 → 2.05; betaine reducers 1.69 → 1.69; one-sided TMA candidates 1.35 → 1.36. The mcrA taxonomy is unchanged (section 6).
- **Like-for-like with v1 cut-offs** for a label that got a v2 identity cut (mtbA; also the new labels mtmC, amtB, glyA, sdaA, dmgo_euk): `reads + low_identity_best` in `per_sample_KO.tsv` (v1 had no identity cut for these, so its `low_identity_best` was 0).
- `methanol_taxonomy.tsv` covers 38 samples, as in v1: R2_restored_H2O has no mtaA / mtaB read in either survey.
- `read_taxonomy.tsv` now covers 163 labels: 1,230,189 rows, 138 MB (above GitHub's 100 MB file limit if this folder is ever pushed). For convenience, `read_taxonomy_v1labels.tsv` (same columns, v1 labels only, 12 MB) sits alongside it; the full file is still the contract file. No `read_taxonomy_v2.tsv` / `binned_unbinned_split_v2.tsv` aliases were written: the contract names were used from the start.

## 6. Regression check v1 → v2 (`scripts/12_regression_v1_v2.py`; `results/regression_*.tsv`)

- **Inputs identical where they should be**: the 8 marker families, `euk_marker` and every v1 counted class have identical references and clusters (section 2). The v2 summary code reproduces v1 exactly on v1 inputs (section 4).
- **GE depth** (mean of the 8 marker depths) v2 / v1 − 1: −0.12 to −0.42 % per sample (median −0.20 %).
- **Shared labels**, per sample × label with ≥ 20 v1 reads (1,160 pairs, 39 labels; v1 has 42 labels incl. 3 merged classes): **median |relative deviation| of copies per GE = 0.36 %; maximum 96 %** (cntB, R1_C_D1). Excluding the four labels with explained large shifts (mtbA, cntB, cutD, mtmB in three waters): the median is 0.30 % and, also excluding mpnS / phnL / mtaA, the maximum is 9.1 % (phnH, SF2_A_D1).
- **Key genes** (pooled reads v1 → v2; median / max |Δ| over samples with ≥ 20 reads; per sample in `regression_key_genes_copies_per_GE_wide.tsv` and `regression_v1_v2_per_sample.tsv`):

| gene | reads v1 → v2 | median \|Δ\| | max \|Δ\| | note |
|---|---|---|---|---|
| mcrA / mcrB / mcrG | 1,046 → 1,045 / 751 → 749 / 504 → 503 | 0.2 % | 0.4 % / 3.4 % / 0.2 % | |
| phnJ | 30,962 → 30,889 | 0.2 % | 2.8 % | |
| grdI / grdH | 16,733 → 16,724 / 2,298 → 2,294 | 0.2 % | 0.4 % | |
| cutC / cntA | 18,402 → 18,308 / 31,794 → 32,229 | 0.2 % / 1.2 % | 2.4 % / 4.9 % | cntA gains reads that v1 sent to the old `Rieske_other` / `cmo` decoys |
| mtaB / mtbB / mtgB | 8,429 → 8,340 / 17,803 → 17,779 / 18,800 → 18,641 | 0.6 % / 0.2 % / 0.8 % | 4.0 % / 2.2 % / 4.0 % | |
| mttB / MttB superfamily | 587,034 → 582,016 / 606,912 → 601,701 | 0.7 % | 1.3 % / 1.1 % | ≤ 1 % of mttB reads go to new labels (dmgdh, acsE, metH, betA_euk, dmdC_like) |
| mtaA | 14,835 → 14,254 | 3.9 % | 12.9 % | new UROD-family competitors (mtsA target; mtgA / mtqA / mtvA / mtoA decoys) |
| mtmB | 21,603 → 21,109 | 1.7 % | 78 % | the large shifts are three waters with 23–44 v1 reads: low-identity hits that now fail E ≤ 1e-5 or best-hit `soxA (monomeric)_like` |
| **mtbA** | 14,479 → 9,018 | 40 % | 64 % | **the new 60 % identity cut** (section 3): 37.7 % of v1's mtbA reads were < 60 % identical to their best hit. Like-for-like with v1 rules: `reads + low_identity_best` = 14,189 (0.98 × v1; e.g. R1_A_D1 774 → 772, R2_A_D1 919 → 910) |
| **cutD** | 5,031 → 3,422 | 35 % | 49 % | v1 counted ferredoxin-domain reads: they now best-hit the ferredoxin decoys (fwdF_like, fwdH_like, korD, mvhD, rnfB) or tmd_like |
| **cntB** | 488 → 149 | 83 % | 96 % | most v1 cntB reads now best-hit the new target dmmB (DMA monooxygenase reductase) |
| the 8 markers | e.g. uL2 435,071 → 434,358 | 0.05–0.3 % | 0.3–2.4 % | |

- **Why reads move** (`regression_lost_reads_<sample>.tsv`, 5 samples): most v1-counted reads that v2 does not count for the same label have **no passing hit** in v2. They are hits to unchanged v1 references near the E ≤ 1e-5 cut-off, and DIAMOND E-values scale with DB size (the v2 DB has ~2.5× the letters). Checked for R1_A_D1 (`12b_lost_read_evalues.py`, `regression_lost_read_evalues_R1_A_D1.tsv`; markers, phnL, phnM, mttB, cutD, mtmB): the 416 lost reads had a v1 best-hit E-value of median 5.4e-6 (90th percentile 8.2e-6). For the 344 whose v1 representative is still in their v2 top 5, the same read–reference hit has a median E of 1.8e-5 (median ratio v2/v1 2.97), and 311 of them now fall between 1e-5 and 1e-4. The other 72 lost their v1 reference from the top 5 to other hits (mostly cutD / mttB reads that now best-hit decoys). This explains the small losses of the markers, phnL (−5.7 %) and phnM (−1.3 %) and part of mtmB. The rest go to new counted labels or decoys, as listed in the table.
- **Pooled mcrA genus shares, 12 unrestored cores** (`regression_mcrA_genus_pooled.tsv`; v1 names mapped to GTDB R232): 830 reads in both. *Methanolobus* 57.7 % → 57.7 %, *Methanohalophilus* 17.2 → 17.0 %, *Methanogenium* 6.5 → 6.9 %, *Methanococcoides* 5.9 → 6.1 %, uncultured "methanogenic archaeon" 2.9 → 2.9 %, *Methanomethylophilus* 1.9 → 1.9 %, *Ca.* Methanospirare 1.6 → 1.6 %. The only genus that changes is *Methanomassiliicoccus* (0 → 6 reads, 0.7 %); in v1 those reads had another Methanomassiliicoccales name.

## 7. Deviations from v1 and from the brief

1. **KEGG caps** for new KOs: targets keep every archaeal gene + 1,500 bacterial genes, decoys 1,000 genes (v1: 3,000 for large KOs), with per-KO seeds. The v1 KOs' files were copied unchanged.
2. **Eukaryote classes** `<gene>_euk` (KEGG eukaryote genes of the target KOs, plus a eukaryote resolution of unbinned salterns proteins). v1 had only `euk_marker`.
3. **No-threshold KOs** dmoA / phnY / phnZ: thresholds derived from the KEGG / Swiss-Prot references (section 2). mpnS / hepD keep the v1 rule.
4. **Taxonomy**: KEGG references carry the GTDB R232 genus of their genome; marker clusters also get the v1-style transfer (05b), with first_pass MAG proteins added to the transfer references.
5. **New identity cuts** (mtmC, mtbA, amtB, glyA, sdaA, dmgo_euk). mtbA is a v1 label, so its counts drop by 38 % (section 6).
6. **Validation**: the grid was reduced to MARGIN 0/2 × PMIN 0/50/60/70 at EMAX 1e-5 / LMIN 30, and only `--mid-sensitive` was run (production mode; v1 compared three modes). A `-k 25` run was added as a check.
7. **Thread budget 96** (v1: 64): DIAMOND ran as 6 × 16 threads (v1: 8 × 8).
8. Not done: no full re-mapping (as v1); no KBase writes; no commits.

## 8. Caveats

- **Generic genes** (amtB, glyA, folD, fhs, sdaA, fdh, frmA, fghA, purU, ackA, pta, acs …): read counts include genes of any taxon. KOfam-twilight homologs are partly counted (e.g. fhs, folD, pta in the twilight test), and a few reads of distant homologs without any panel-HMM hit are counted at ~10⁻⁴ of non-target reads (section 3). Copies per GE > 1 are normal for multi-copy genes.
- **Archaeal / methanogen shares** come from the best-hit reference. For salterns-only clusters the domain is transferred from the nearest classified protein at any identity, so it is reliable for domain but not for lineage below ~60 % identity. `reads_domain_unknown` counts reads whose best-hit cluster has no domain (decoy-only salterns clusters are not counted anyway).
- **Flagged / merged / not-validated labels**: section 3. Use the merged classes for mttB / mtgB, grdH / grdB, pmoA / CuMMO, the methylamine corrinoid proteins, gsdmt / sdmt, dmgdh / sardh, the Pta family, mtaA / mtsA and cntB / dmmB. mpnS, phnY, phnZ, dmoA, dddX, AOX, CMO, mdh and hepD could not be validated; their (small) counts are unvalidated.
- **The three ferredoxin KOs** (fwdH, fwdF, fwdG) passed the precision bar in the simulation (≥ 0.98), and their reads are 90–100 % archaeal in the cores. Their real-data specificity rests on the `_like` decoys and could not be tested beyond the simulation.
- **Larger DB, larger E-values**: about 0.1–0.7 % of v1-counted reads per label fall below E ≤ 1e-5 in v2, purely from the DB size (section 6). The effect is nearly uniform across the marker families, so copies per GE are almost unaffected (GE depth −0.2 %).
- **Runtime**: poplar was oversubscribed by other users throughout (load 450–1,500 on 336 cores). Each 16-thread job got about 4–8 cores.

## 9. Files

| here | content |
|---|---|
| `scripts/` | every script that ran (identical copies on poplar; md5 in `poplar_copies/scripts_md5_*.txt`) |
| `meta/` | `ko_panel_v2.tsv`, `ko_symbols.tsv`, `synthetic_thresholds.tsv`, `merged_v2.json`, `kegg_org_assembly.tsv`, `kegg_org_gtdb.tsv`, fetch subsets |
| `cutoffs.tsv` | the classification cut-offs (v1 + 6 per-label identity cuts) |
| `validation/` | `validation_summary_v2.tsv` (per label: recall, precision, flags, confusions), `validation_merged_classes_v2.tsv`, per-label / confusion / divergence / twilight tables (`*_final.tsv` = final cut-offs, `*_mid-sensitive_k5/k25.tsv` = v1 cut-offs), `grid_*.tsv`, `grid_summary.tsv`, `identity_cuts.tsv`, `background_fp_final.tsv`, `holdout_genes.tsv` |
| `refs_small/` | `rep_info.tsv.gz`, `labeled_meta.tsv.gz`, `clusters99.tsv.gz`, `label_lengths.tsv`, `labels.tsv`, `excluded_reps.tsv`, `euk_resolution.tsv`, `k10944_resolution.tsv`, `marker_rep_genus.tsv.gz`, fetch logs |
| `results/` | contract files (section 4) + `marker_genus_profile*.tsv`, `methanogenesis_pathway_archaeal_share.tsv`, `summary_per_label_v2.tsv`, `regression_*.tsv`, `check_labels_v1_v2.tsv`, `genome_equivalents.tsv`, `per_guild_read_vs_MAG.tsv`, `assembly_*`, `mag_attributable_reads.tsv`, `per_sample/*.counts.tsv` |
| `poplar_copies/` | logs, DIAMOND `/usr/bin/time` reports, benchmark reports, re-mapping summaries, DB md5 |

On poplar only (large): `refs/` (labelled FASTA, clusters), `db/survey.{faa,dmnd}`, `hmm_new/`, `dmnd/` (19 GB), `map/`, `results/per_sample/*.reads.tsv.gz`, `val/`, `refs_raw/`, `genes_firstpass/`.
