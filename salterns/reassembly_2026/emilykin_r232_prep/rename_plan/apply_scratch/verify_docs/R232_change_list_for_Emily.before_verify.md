# GTDB R232 update: change list for Emily

Prepared 2026-11-12 for Emily Kin (from Andrew Freiburger's GTDB R232 update of the EmilyKin MAG catalogue).

**What happened.** All 291 MAGs were re-classified with GTDB-Tk v2.8.0 against GTDB release R232, the release this
GTDB-Tk version expects. The earlier taxonomy was GTDB-Tk v2.7.2 with GTDB r220, which ran with a version-mismatch
warning. The MAG iterativeIDs were then re-derived with the unchanged naming rule: 51 of 291 MAG IDs change, and no ASV
ID changes. The old-to-new pairs are listed by bin in `mag_iterativeID_r220_to_r232.csv` (repo root).

**How to use this list.** Section 1 is keyed by page and line of your newest draft, `paper3/papert3_071626.pdf`
(Word export of 2026-10-21; there is no DOCX of it on disk). The same edits were made as tracked changes, with Word
comments, in the newest editable files: `papert3_070626 (2)_correction.docx` and `SI-paper3 (3)_correction.docx`
(author "Claude (GTDB R232 update)"). Items marked **science** change a scientific statement rather than a name.

## 1. Main text (`paper3/papert3_071626.pdf`)

| # | page | line(s) | now | change to | type |
|---|---|---|---|---|---|
| 1 | 7 | 229–230 | "generated using GTDB-Tk v2.7.2 against Genome Taxonomy Database release r220 (50, 51)" | "generated using GTDB-Tk v2.8.0 (classify_wf) against Genome Taxonomy Database release R232 (50, 51)". Consider citing the R232 release; the plan proposes Parks et al. 2025, *Nucleic Acids Res.*, doi:10.1093/nar/gkaf1040 (please check this reference before adding it) | version |
| 2 | 8–9 | 307–310 | "spanned 14 bacterial phyla: … Patescibacteria (21) … Myxococcota (3) …" | "spanned **15** bacterial phyla: … **Patescibacteriota** (21) … **Myxococcota (2), Myxococcota_A (1)** …". R232 renames Patescibacteria (the same 21 MAGs) and places one former Myxococcota MAG (coasm_bin.258) in Myxococcota_A. Your choice: if suffixed GTDB names are merged (as Figure 1 does for genera), keep 14 and "Myxococcota (3)" | number |
| 3 | 9 | 321–322 | "with four streamlined Patescibacteria (0.67–1.17 Mb)" | "Patescibacteriota". The closed genomes still span five phyla | name |
| 4 | 10 | 345 | "For each genome: GTDB taxonomy (lowest two named ranks)" | "GTDB R232 taxonomy" | version |
| 5 | 10–12 | Table 1 | 28 of the 54 lineage strings | the R232 lineages in section 1.1 below (by bin) | name |
| 6 | 13 | Figure 1 (image) | R220 genus bands | replace the image with the re-rendered `figures/percell_composition.png` (same script, data and colours). Legend changes: SCMT01 → *Lysobacter* (GTDB *Lysobacter_I*), *Azospira* → *Proximibacter*, *Giesbergeria* → *Simplicispira*, "UBA5532 (Patescibacteria)" → "UBA5532 (Patescibacteriota)". JJ008 is now CAN_5_bin.112 only (coasm_bin.492 is CADEDV01 and is pooled in Others). "Unclassified" falls from 1.89 % to 0.57 % (11 MAGs gained a genus) and becomes the last band. The *Competibacter*/*Accumulibacter* values do not change. The legend text (l350–358) needs no change | figure |
| 7 | 13 | 359–360 | "48 were substantially divergent … (ANI <85% …; a further 77 could not be confidently placed …" | **28** … a further **66** … (Table S6 rebuilt on R232 with the same skani settings; 23 of the 28 have no skani hit at all) | number |
| 8 | 15 | 389–391 | "three of these were closely related to existing references (≥97% ANI), whereas the near-complete coasm_bin.476 (99.98% completeness, 0.20% contamination) aligned poorly to any reference and likely represents a novel Accumulibacter lineage." | "all six belong to GTDB R232 species, and five are closely related to their species representative (≥97% ANI), including the near-complete coasm_bin.476 (99.98% completeness, 0.20% contamination; *Accumulibacter* sp032229525, 98.68% ANI)." See section 3 | **science** |
| 9 | 15 | 408–412 | "The IC genome (coasm_bin.250) falls near Ca. A. delftensis … The IA genome (coasm_bin.185) sits near proposed Ca. A. regalis …" | optional: the R232 genome-level species now agree with the ppk1 clades. coasm_bin.250 is *Accumulibacter delftensis* (98.60 % ANI), coasm_bin.185 is *A. regalis* (97.58 %) and coasm_bin.347 is *A. phosphatis* (98.14 %). You could add "(GTDB R232: *A. delftensis*)" / "(GTDB R232: *A. regalis*)" | **science** (optional) |
| 10 | 15 | 416–418 | "The ID genome (coasm_bin.476) is the strongest novelty point: … and it aligned poorly to any reference (nearest only to the sparsely characterized BPBW2884/BPBW2899)" | e.g. "The ID genome (coasm_bin.476) belongs to a recognized ppk1-defined Type I subclade that remains poorly characterized compared with IC and IIA: its ppk1 is nearest only to the sparsely characterized BPBW2884/BPBW2899, although the genome falls in the GTDB species *Accumulibacter* sp032229525 (98.68% ANI); …" The genome-level "aligned poorly to any reference" no longer holds; the ppk1-level statement does | **science** |
| 11 | 15 | 422 | "The nearest GTDB references are largely metagenome-derived drafts" | please re-check against the R232 representatives: coasm_bin.185 GCA_000585075.1 (*A. regalis*), coasm_bin.250 GCF_012939955.1 (*A. delftensis*), coasm_bin.347 GCF_000024165.1 (*A. phosphatis*), coasm_bin.476 GCF_032229525.1 (sp032229525, added to GTDB in R226), CAN_1_bin.98 GCF_024235695.1 (sp014879755), CAN_4_bin.64 GCF_029255145.1 (sp003538495, 96.05 %) | author check |
| 12 | 28 | 711–712 | "including a novel near-complete Clade ID lineage MAG nearest only to sparsely characterized references." | e.g. "including a near-complete Clade ID MAG (GTDB R232 *Accumulibacter* sp032229525) whose ppk1 is nearest only to sparsely characterized references." | **science** |
| 13 | 1 | 17–18 (abstract) | "including a near-complete, sparsely characterized Clade ID lineage" | no change needed (true at the ppk1 level); optional "… Clade ID lineage (GTDB *Accumulibacter* sp032229525)" | optional |
| 14 | 24 | Figure 4 | x-axis labelled by MAG bin | no change. (The older ID-labelled Figure 4 in the working DOCX would need 17 relabels; see the comment there) | none |

Unchanged on purpose (checked against the R232 data): Figure 2, Figure 3 and Figure 5; Tables 2–3; "eight
*Competibacter* MAGs" and "six *Accumulibacter* MAGs" (p9 l329–332; coasm_bin.365 was never in the GTDB
*Competibacter* genus: r220 JARSSO01, R232 *Contendibacter*); the guild counts (four canonical PAO, five canonical GAO,
41 strong putative PAO-like, 118 other); the nosZ counts; the 5-class denitrifier counts; all network results (the
MAG link of the network analyses was frozen by bin, so renaming MAGs changes no network node, edge or module).

### 1.1 Table 1 lineages that change (pp. 10–12; keyed by bin)

| bin | what changes | R232 lineage |
|---|---|---|
| CAN_1_bin.174 | p__ Patescibacteria → Patescibacteriota; c__ JAEDAM01 → Gracilibacteria | `d__Bacteria;p__Patescibacteriota;c__Gracilibacteria;o__BD1-5;f__UBA2023;g__UBA5532;s__UBA5532 sp002426885` |
| CAN_1_bin.195 | p__ Patescibacteria → Patescibacteriota; c__ ABY1 → Patescibacteriia; o__ BM507 → Patescibacteriales; f__ UBA917 → Patescibacteriaceae; g__ UBA919 → Patescibacterium | `d__Bacteria;p__Patescibacteriota;c__Patescibacteriia;o__Patescibacteriales;f__Patescibacteriaceae;g__Patescibacterium;s__` |
| CAN_1_bin.2 | o__ SJA-28 → Tepidaquicellales; f__ B-1AR → Tepidaquicellaceae; g__ (blank) → JAFDZO01 | `d__Bacteria;p__Bacteroidota;c__Ignavibacteria;o__Tepidaquicellales;f__Tepidaquicellaceae;g__JAFDZO01;s__` |
| CAN_1_bin.48 | p__ Patescibacteria → Patescibacteriota; c__ Paceibacteria → Minisyncoccia; g__ (blank) → CAYHLH01 | `d__Bacteria;p__Patescibacteriota;c__Minisyncoccia;o__UBA6257;f__GWB1-50-10;g__CAYHLH01;s__` |
| CAN_3_bin.221 | o__ SJA-28 → Tepidaquicellales | `d__Bacteria;p__Bacteroidota;c__Ignavibacteria;o__Tepidaquicellales;f__OLB5;g__OLB5;s__OLB5 sp002344125` |
| CAN_3_bin.227 | c__ Kapaibacteriia → Ignavibacteria | `d__Bacteria;p__Bacteroidota;c__Ignavibacteria;o__Kapaibacteriales;f__Kapaibacteriaceae;g__Kapaibacterium;s__` |
| CAN_4_bin.208 | g__ Azospira_A → Proximibacter | `d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Rhodocyclaceae;g__Proximibacter;s__` |
| CAN_5_bin.177 | g__ SCMT01 → Lysobacter_I; s__ SCMT01 sp019510125 → Lysobacter_I sp019510125 | `d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Lysobacter_I;s__Lysobacter_I sp019510125` |
| CAN_5_bin.2 | g__ (blank) → JBDHTQ01 | `d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__UBA2400;f__UBA2400;g__JBDHTQ01;s__` |
| CAN_5_bin.70 | s__ Azonexus sp017302355 → (blank) | `d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Rhodocyclaceae;g__Azonexus;s__` |
| CAN_5_bin.90 | p__ Patescibacteria → Patescibacteriota; c__ JAEDAM01 → Gracilibacteria | `d__Bacteria;p__Patescibacteriota;c__Gracilibacteria;o__BD1-5;f__UBA2023;g__UBA5532;s__UBA5532 sp002426885` |
| coasm_bin.10 | s__ (blank) → Leptovillus sp964405285 | `d__Bacteria;p__Chloroflexota;c__Anaerolineae;o__Promineifilales;f__Promineifilaceae;g__Leptovillus;s__Leptovillus sp964405285` |
| coasm_bin.171 | f__ Burkholderiaceae → Burkholderiaceae_C | `d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae_C;g__Inhella;s__Inhella sp017302655` |
| coasm_bin.185 | s__ Accumulibacter sp000585075 → Accumulibacter regalis | `d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Rhodocyclaceae;g__Accumulibacter;s__Accumulibacter regalis` |
| coasm_bin.2 | f__ Burkholderiaceae → Burkholderiaceae_C | `d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae_C;g__Ottowia;s__` |
| coasm_bin.250 | s__ Accumulibacter phosphatis_G → Accumulibacter delftensis | `d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Rhodocyclaceae;g__Accumulibacter;s__Accumulibacter delftensis` |
| coasm_bin.321 | s__ (blank) → Dokdonella sp964390215 | `d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Xanthomonadales;f__Rhodanobacteraceae;g__Dokdonella;s__Dokdonella sp964390215` |
| coasm_bin.337 | o__ SJA-28 → Tepidaquicellales; f__ B-1AR → Tepidaquicellaceae | `d__Bacteria;p__Bacteroidota;c__Ignavibacteria;o__Tepidaquicellales;f__Tepidaquicellaceae;g__UBA9002;s__` |
| coasm_bin.347 | s__ Accumulibacter sp017302455 → Accumulibacter phosphatis | `d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Rhodocyclaceae;g__Accumulibacter;s__Accumulibacter phosphatis` |
| coasm_bin.365 | g__ JARSSO01 → Contendibacter | `d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Competibacterales;f__Competibacteraceae;g__Contendibacter;s__` |
| coasm_bin.40 | g__ Flavipsychrobacter → JAFIGV01 | `d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Chitinophagales;f__Chitinophagaceae;g__JAFIGV01;s__` |
| coasm_bin.489 | s__ (blank) → UBA2287 sp964479665 | `d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Tenuifilaceae;g__UBA2287;s__UBA2287 sp964479665` |
| coasm_bin.521 | c__ Kapaibacteriia → Ignavibacteria; s__ (blank) → PGYR01 sp041595255 | `d__Bacteria;p__Bacteroidota;c__Ignavibacteria;o__Kapaibacteriales;f__UBA2268;g__PGYR01;s__PGYR01 sp041595255` |
| coasm_bin.53 | s__ (blank) → OLB9 sp964501555 | `d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Chitinophagales;f__Saprospiraceae;g__OLB9;s__OLB9 sp964501555` |
| coasm_bin.563 | f__ Burkholderiaceae → Burkholderiaceae_C | `d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae_C;g__Hydrogenophaga;s__` |
| coasm_bin.578 | c__ Verrucomicrobiae → Verrucomicrobiia | `d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiia;o__Chthoniobacterales;f__Terrimicrobiaceae;g__QYQL01;s__QYQL01 sp014879575` |
| coasm_bin.579 | o__ Flavobacteriales → Flavobacteriales_F; s__ (blank) → PHOS-HE28 sp035441065 | `d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales_F;f__PHOS-HE28;g__PHOS-HE28;s__PHOS-HE28 sp035441065` |
| coasm_bin.88 | o__ Rhodospirillales_A → Rhodospirillales; f__ UXAT02 → Defluviicoccaceae | `d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Rhodospirillales;f__Defluviicoccaceae;g__Defluviicoccus;s__` |

The other 26 rows of Table 1 do not change.

## 2. Supplementary Information (`SI-paper3 (3)_correction.docx`; tracked changes there)

Your newest SI is not on disk, so these items are keyed by SI section and the start of the paragraph.

| where | change | type |
|---|---|---|
| SI-3, ppk1 paragraph ("All six Accumulibacter MAGs carried a recoverable, full-length ppk1 …") | sentence added after "… coasm_bin.476 to clade ID (100/100).": "At genome level, GTDB R232 assigns coasm_bin.185, coasm_bin.250 and coasm_bin.347 to *Accumulibacter regalis*, *A. delftensis* and *A. phosphatis*, matching these clades, and coasm_bin.476 to *Accumulibacter* sp032229525 (98.68% ANI)." Not changed, but flagged in a comment: "CAN_1_bin.98 was a divergent, candidate-novel Type I lineage" refers to its ppk1 only; its genome is the known GTDB species *Accumulibacter* sp014879755 (99.48 % ANI) in both releases | **science** |
| SI-6 ("A quality audit of the 188 overwrites …") | "*Giesbergeria* CAN_5_bin.40" → "*Simplicispira* CAN_5_bin.40"; "closest GTDB r220 reference genomes" + ", labelled with their GTDB R232 names," | name |
| SI-7 iterativeID paragraph ("… On the best-of-both 291 catalogue this yields …") | "155 ASV-inherited, 62 GTDB-genus, and 74 GTDB-lowest" → "on the catalogue classified with GTDB R232 … **138 / 70 / 83**" | number |
| Figure S4 legend | "61 GTDB r220 species-representative reference genomes" + "(labelled with their GTDB R232 names; GCA_017302355.1 is not in R232)"; "rooted on the **Patescibacteriota** clade"; CAN_4_bin.208 is *Proximibacter* in R232, and its reference *Azospira_A restricta* (GCF_016858125.1) is *Proximibacter restrictus* in R232 (the reference was chosen on r220, so "r220 *Azospira_A*" stays for the selection rule); "best skani ANI to a GTDB **R232** representative (Table S6)". The tree itself was not rebuilt (relabel only) | name |
| Figure S4 image | replace with the re-rendered `fig1_tree/fig1_tree.png`. It has the same tree, 8 MAG tip labels change (e.g. Ignavibacteria.21 → Tepidaquicellales.1_m, Azospira.3 → Proximibacter.2_m, Ca_Competibacter.32 → Contendibacter.1_m), 9 reference labels show R232 names, and the red circles follow the rebuilt Table S6: 10 MAGs change ANI class, e.g. CAN_5_bin.70 ≥95 → 85–95 because its r220 reference left GTDB | figure |
| Figure S5 image | replace with the re-rendered `clade_classify/out/nosz_top5_tree.png`. The only change is that rank 5 (CAN_3_bin.221) is now labelled Tepidaquicellales.1_m; the legend is unchanged | figure |
| Tables S1, S3–S7 (`MAG_target_gene_filled_UPDATED (1).xlsx`) | updated by bin, with a backup `….pre_r232.bak.xlsx`. MAG iteration ID: 82 cells per full sheet (51 R232 renames plus 31 IDs that still carried an older map). GTDB rank and lineage columns changed (130 lineages). Table S6 was rebuilt on R232: ≥95 % **121** / 85–95 % **76** / <85 % **28** / unresolved **66** (r220: 96 / 70 / 48 / 77), reference database "GTDB R232 (GTDB-Tk v2.8.0, classify_wf)". Table S7 classifications are unchanged. The per-cell, s16 and 5-class columns were not touched. `cladeexcel.xlsx` (S4 source) and `abundance and coverage (1).xlsx` were updated the same way | table |

## 3. Scientific changes to be aware of

1. **coasm_bin.476 is no longer a candidate novel species.** R232 places it in *Accumulibacter* sp032229525 (GTDB-Tk
   98.68 % ANI, AF 0.84 to GCF_032229525.1, a species representative added to GTDB in R226; skani AF_query 72.9 %). Its
   clade ID ppk1 is still nearest only to sparsely characterized references, so the ppk1-level novelty stands.
2. **All six *Accumulibacter* MAGs now have a GTDB species, and five are ≥97 % ANI** (coasm_bin.185 97.58,
   coasm_bin.250 98.60, coasm_bin.347 98.14, CAN_1_bin.98 99.48, coasm_bin.476 98.68; CAN_4_bin.64 96.05). The R232
   species of the three named-clade MAGs match their ppk1 clades (IA *regalis*, IC *delftensis*, IIA *phosphatis*).
3. **Species-level resolution rises from 97 to 122 MAGs** (26 gained, 1 lost: CAN_5_bin.70, whose r220 reference
   genome left GTDB). Genus-level: 273 → 283 (172 → 183 genera). RED-flagged novelty: 49 → 21.
4. **Table S6 (skani ANI novelty):** ≥95 % 96 → 121, 85–95 % 70 → 76, <85 % 48 → 28, unresolved 77 → 66.
5. **One Table S6 call is an artefact of the skani settings.** coasm_bin.535 is "low_AF_unresolved", because
   `-n 3` keeps only its three highest-ANI hits, all tiny alignments. GTDB-Tk assigns it to *Novosphingobium*
   sp017302615 (98.4 % ANI, AF 0.90).
6. **Novel-bacteria module (if you use it):** coasm_bin.260 (Ca_Sarcinithrix.1) now has a GTDB species (JAAEKA01
   sp049294015, 98.33 %), so six of the seven "species-novel" MAGs remain novel. The novel-bacteria bubble figure now
   uses the R232 set. The metabolic-potential figure and Tables N2–N3 still use the r220 set. Which set to keep is open.
7. **15 phyla** if Myxococcota_A is counted separately (item 2 above).
8. **Figure S4 red circles:** 10 MAGs change ANI class (from the rebuilt Table S6).

## 4. Name changes you may meet elsewhere

- Old MAG IDs that are also ASV names (for example Azospira.3, Ca_Competibacter.32, Ignavibacteria.5/.21,
  Ca_Falkowbacteria.7) remain valid **as ASV names**. Network figures and the SPIEC-EASI write-up label nodes by ASV,
  so they keep these names. The SPIEC-EASI document now says which MAG stands behind such a node, e.g. "Azospira.3
  [MAG CAN_4_bin.208, Proximibacter.2_m under GTDB R232]".
- coasm_bin.365 is *Contendibacter* in R232 (MiDAS spells it *Contendobacter*). Its GAO call is unchanged (family
  Competibacteraceae). The Table S7 rule now accepts both spellings, so no classification changed.
- Full map: `mag_iterativeID_r220_to_r232.csv` (bin, old_id, new_id, case_old, case_new).

## 5. Open choices

1. 15 vs 14 phyla (item 2).
2. Whether to cite the R232 release paper (item 1).
3. Figure S4: relabel only (done; the tree keeps its 61 r220 references) or rebuild on R232 references.
4. Novel-bacteria module: keep the r220 set of seven or move to the R232 set.
5. SI ppk1 text: "candidate-novel" wording for CAN_1_bin.98 (pre-existing; flagged, not changed).
