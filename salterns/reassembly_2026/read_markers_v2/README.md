# Read survey v2: the methylated-compound / methanogenesis KO panel (138 KOs) in the 39 salterns metagenomes

*Track `read_markers_v2`, started on poplar on 2026-11-25 (all times CST). v1 (`../read_markers/`, poplar `read_markers/`) is only read, never written. Nothing is committed and nothing is written to KBase. Poplar work dir: `/scratch1/afreiburger/salterns_2026/read_markers_v2/`; this folder holds the scripts, the small metadata and the results.*

**Status (2026-11-25 3:40 PM CST):** deliverable 1 (`results/marker_genus_profile.tsv`) is DONE. The v2 reference DB is being built (KEGG fetch and HMM search running); validation, the full DIAMOND re-run and the summaries follow.

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

## 2. v2 panel and reference DB (in progress)

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

### Salterns HMM search (`scripts/03_panel_hmmsearch_v2.sh`, 3:30–~3:55 PM CST)

HMMER 3.4 with the v1 flags (`--cpu 1 --noali -E 1e-5 --domE 1e-5 -Z 20000000`). The 224 profiles not in v1, split into 4 chunks, were searched against the 39 per-sample gene sets and the 520 co-assembly MAGs (v1 `genes/`, 13.6 M proteins): 160 single-CPU jobs, 80 at a time, nice 10. Labelling uses v1's 62-profile tables together with these. With the fixed `-Z`, each profile's hits do not depend on the other profiles searched, so the union equals a fresh 286-profile search.
