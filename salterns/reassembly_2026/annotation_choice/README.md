# Annotating the 2026 re-assembly MAGs: DRAM or Bakta?

*Written 2026-11-11 for the saltern methane manuscript. Scope: the 1,324 re-assembly MAGs (≥50 % complete, <10 % contaminated) in 558 species clusters (skani ANI ≥95 %, AF ≥50 %), from `../comparison/`. Nothing was run on KBase or poplar. poplar needs a Duo login that wasn't available.*

## Short answer

1. **The paper needs threshold-based KO calls (KOfam), not either tool's defaults.**
   - The guilds (MePn consumers, betaine reducers, TMA routes, methanogens) are gene-presence rules.
   - DRAM's KO layer is KOfam with KEGG's per-KO score thresholds, which is what those rules need. Its other layers (dbCAN, MEROPS, Pfam) add little to a methane and MePn paper.
   - **Neither DRAM distillate covers the paper's markers out of the box.**
     - DRAM 1.5 has no C-P lyase, betaine reductase, CutC or CntA rows.
     - DRAM2 puts *grdH*, *cutC* and *cntA* only in a sheet that loads with `--sum_custom`.
   - I wrote a custom distillate sheet for either version.
2. **Bakta adds little for this paper.**
   - Bakta emits a KEGG KO cross-reference only when a protein is ≥90 % identical, over ≥80 % of its length, to a UniRef90 seed whose cluster carries a KOfam label. Those labels exist for 24.3 M of 135.3 M clusters.
   - These MAGs are mostly novel: 374 of the 558 best-per-species MAGs are in GTDB placeholder genera, and only 226 have a GTDB species. Expect sparse KOs.
   - The `light` database gets no KO, EC, COG or GO from the UniRef layer; it has only UniRef50 cluster information for CDSs.
   - Bakta is built for bacteria only. 116 MAGs (47 species) are archaea, including the *Methanolobus* the paper depends on.
   - Bakta's real strength is INSDC-ready GFF3/GenBank/EMBL output. Neither GenBank nor ENA requires annotation for MAG submission, and NCBI will run PGAP on request.
3. **I tested DRAM's KO layer locally on 38 original MAGs, using RAST as truth, and on 8 new ones** (section 3). The 38 are every RAST guild member in `guilds_r232.csv` plus 4 other taxa. This found problems that matter more than the choice of tool:
   1. **MePn (*phnJ*) mostly holds: KOfam calls it in 21 of the 22 RAST MePn consumers.**
      - 18 have all of *phnGHIJKLM*.
      - 2 lack only *phnH* (Marivita.2, Haloquadratum.4).
      - Phaeodactylibacter.1 has a lone *phnJ*.
   2. **The exception is Albimonas.1, and the problem is not an assembly artefact.**
      - Its *phnI* is split in two, and its *phnJ* (371 aa, 37 % HMM coverage) scores 134. The threshold is 247 in the 2026-05 KOfam release used here; the 2026-10 release uses a 195 full-sequence threshold.
      - The identical broken *phnGI(I)J* fragment recurs in 2 of the 5 re-assembly *Albimonas* bins, from independent assemblies.
      - The best new bin (95.2/2.4) has *phnG, I, K, M, N* and a just-below-threshold *phnL*, but **no *phnJ* and no *phnH* gene at all**.
      - **No *Albimonas* bin has a KOfam-callable *phnJ*.** Whether *Albimonas* is a MePn consumer is a manuscript-level question.
   3. **The selenocysteine *grdH* breaks the betaine-reducer rule.** Gene callers stop at the Sec codon. The shortened GrdH then scores just under its threshold (602–606 vs 619) and is labelled as glycine-reductase *grdB* instead.
      - Of the 9 RAST betaine reducers, KOfam *grdH* finds 3.
      - *grdI* (K21578, no Sec) finds 7, as does *grdH* or *grdI*.
      - The other 2 (Halanaerobium.3 at 76 % complete; BM714.4 at 25 % contamination) have only fragments of both genes.
   4. **KO K14083 (*mttB*) is called in all 16 Rhodobacteraceae controls and in 23 of 37 non-methanogens.** These are the RAST "trimethylamine methyltransferase family" proteins, which lack pyrrolysine.
      - DRAM 1's product heatmap would light "trimethylamine → dimethylamine" for *Roseovarius*, *Albimonas* and *Marivita*.
      - In 5 anaerobes the same family is called *mtgB* (K25220, betaine demethylase).
   5. **KOfam finds TMA-route genes that RAST missed.**
      - A CntA/YeaW (with no CntB) in 6 *Roseovarius*. RAST calls it "Rieske 2Fe-2S domain protein".
      - A CutC in Desulfotignum.1 (1,031 vs 990, so near the threshold). RAST calls it "pyruvate formate-lyase", and its activating enzyme scores below the *cutD* threshold.
      - These are candidates to verify, not guild members.
   6. **DRAM's default 2,500 bp contig minimum removes data.** It drops 12 % of all MAG bp: 30 % of MaxBin2 bp, and >25 % of the bp in 299 MaxBin2 MAGs.
      - In two MaxBin2 *Roseovarius* bins it removes 37 % and 47 % of genes, along with their *proVWX*/*opu* transporter calls.
4. **Recommendation (section 6).**
   - Run KOfam on all 1,324 MAGs, with contigs ≥1,000 bp. Locally this is about 25 h on 24 threads.
   - Run DRAM2 (already on h100), with the custom distillate, on the 558 best-per-species MAGs for the metabolism figure.
   - Revise two guild rules: betaine reducer = *grdH* or *grdI*; methylamine route = *mttB/mtbB/mtmB* only in an *mcrA* genome or with an in-frame UAG.
   - Recheck *Albimonas* as a MePn consumer before the paper relies on it.
   - Skip Bakta for now. Run it (or request PGAP) only on whatever set you deposit with annotation.
   - The 8-MAG test on poplar (section 7) should measure Bakta's KO yield on novel taxa, which can't be estimated without its database.

## 1. What the guilds use now

`correlational/r232/rerun_correlational_r232.py` (`GENES`, `build_guilds`) matches regular expressions against the **RAST function strings** of the 310 original MAGs. The results are in `guilds_r232.csv`, columns `ev_*`:

| guild | rule (RAST string) | KO equivalent |
|---|---|---|
| MePn consumer | `c-p lyase\|ec 4.7.1.1` (*phnJ*) | K06163 |
| (context) `ev_phn_accessory` | PhnG/H/I/L synthase, PhnM diphosphatase, PhnN | K06166/5/4, K05780, K06162, K05774 |
| betaine reducer | `betaine reductase component b` (matches both GrdH β and GrdI α) | K21579 / K21578 |
| TMA producer | `carnitine monooxygenase` (*cntA*) or genes "clustered with choline TMA-lyase" | K22443 (= *yeaW*; *yeaX* = K22444); K20038/K20037 for CutC/D themselves |
| methanogen | `methyl.coenzyme m reductase` (*mcrA*) | K00399 |
| `ev_mtbB_mtmB` | DMA/MMA methyltransferases and corrinoid proteins | K16178/K16176 (+ K16179/K16177) |

The new MAGs have no RAST genomes yet; they are Assemblies in ws 275645. Whatever annotation is chosen replaces RAST for the guild calls on the new set, unless RAST is also run for model building (section 4).

## 2. Marker coverage by tool

The full table is **`marker_ko_table.tsv`**: 68 markers, one row per KO. Columns:
- KEGG name and EC;
- KOfam threshold and score type in the 2026-10-26 and 2026-05-28 releases;
- DRAM 1.5.0 and DRAM2 distillate rows;
- the KO→ModelSEED and EC→ModelSEED reactions;
- the `guilds_r232.csv` column each marker feeds.

KO numbers were checked against `rest.kegg.jp/list/ko` on 2026-11-11. Corrections to the task's list:
- Betaine reductase protein B is **GrdH (β, Sec) = K21579** and **GrdI (α) = K21578**; GrdB (K10672) is the *glycine* reductase.
- *yeaW/yeaX* share K22443/K22444 with *cntA/cntB* (KEGG links `eco:b1802/b1803` to them).
- *mtgB* (K25220) is the betaine:corrinoid methyltransferase of the non-pyrrolysine MttB family, a direct betaine route that does not go through TMA.

### 2a. Does the tool call the gene at all?

**DRAM 1.5.0 and DRAM2 (KOfam layer).**
- A gene gets a KO when its HMM score beats that KO's KOfam threshold. DRAM 1.5's `sig_scores` skips KOs whose threshold is `-` (`sources/dram1_v1.5.0/mag_annotator_utils.py`).
- 1,870 of 28,504 KOs have no threshold, including three of ours:
  - **mpnS** K18049 (MePn synthesis);
  - the oxidative-MePn **phnY\*** K23996;
  - **phnZ1** K05915.
  - DRAM can't call these. In KofamScan they'd need a manual cutoff.
- DRAM2 `hmm_parser.py` keeps only `score_type` full or domain. If h100's preformatted `kofam_ko_list.tsv` stores the dash as empty, DRAM2 instead fills threshold 0 and calls every reported hit (open item).
- Every other marker has a threshold, so DRAM calls *phnC–P*, *phnX/W/Y/Z*, the whole *grd* complex, *mtgB*, the transporters, *cutC/D*, *cntA/B*, *torA/Z*, *dmsA*, *tmm*, *mcrABG*, *mttBC/mtbBC/mtmBC/mtbA*, *mtaABC* and *pylBCDS*.
- **Thresholds move between KOfam releases**, so the KOfam release must be recorded in Methods. From 2026-05 to 2026-10:
  - *phnJ* 246.6 (domain) → 194.6 (full);
  - *grdI* 713.6 → 505.9;
  - *grdB* 412.6 → 278.0;
  - *mttB* 160.7 → 218.9.

**Bakta 1.12.2 / DB v6.0.** Read from `bakta/psc.py`, `pscc.py` and `features/annotation.py`, copied in `sources/bakta_v1.12.2/`.
- `db_xrefs` emitted on a CDS:
  - `SO`, `UniParc`, `RefSeq` (`WP_`), `UniRef` (100/90/50);
  - `EC`, `COG` (ID and category), `GO`, `KEGG:Knnnnn`;
  - `Pfam`/`NCBIFam` (hypotheticals), `IS`, `VFDB`, `amrfinder`, `UserProtein`/`UserHMM`.
- **`KEGG`, `EC`, `COG` and `GO` come only from a valid PSC**, i.e. a UniRef90-seed hit with ≥90 % identity and ≥80 % query and subject coverage, or from an exact UPS/IPS match.
- A 50–90 % identity hit becomes a PSCC hit: a UniRef50 product name and `UniRef:UniRef50_*`, with no KO, EC or gene symbol.
- KO cross-references are pre-assigned to 24,267,514 of 135,274,518 UniRef90 clusters (Bakta README, DB v6.0).
- So Bakta's KO coverage falls exactly where these MAGs are novel. Its product names may still identify markers ("methyl-coenzyme M reductase subunit alpha"), but that is string matching again, the RAST problem.
- Bakta's `--hmms` option takes user HMMs, applies the trusted cutoff (TC) in each HMM header, and adds a 1e-6 E-value cap. This is a way to add the KOfam markers inside Bakta. KOfam profiles have no TC line, so the `ko_list` thresholds would have to be written in.

**eggNOG-mapper 2.1.13 / eggNOG DB 5.0.2** (on poplar).
- Orthology-based `KEGG_ko`, EC, KEGG module/pathway, CAZy, PFAMs and BiGG.
- Its KO vocabulary is fixed by the eggNOG 5.0 build, so recent KOs such as *mtgB* K25220 and *phnY\** K23996 can't appear. Whether *grdH/grdI* and *cntA/B* are covered is part of the test.
- In EmilyKin, eggNOG left one MAG nearly empty (3 of 42 target KOs vs 26 with KOfam). KOfam added 198 MAG×KO calls over eggNOG across 291 MAGs (`~/Documents/EmilyKin/MANUSCRIPT_methods_results_discussion.md` §2.6).

### 2b. Does the distillate show it?

**DRAM 1.5.0** (`metabolism_summary.xlsx` and `product.html`):
- Shown: *phnCDE*, *proVWX* and *opuA/BD/C* as transporters; *torA* (Energy > C1; heatmap "TMAO reductase"); *mcrABG* and the *mtt/mtb/mtm/mta* genes (Energy > C1-methane; heatmap "Key functional gene", "trimethylamine => dimethylamine", …).
- **Not shown:** any C-P lyase gene (*phnG–N*), any *grd*, *cutC/D*, *cntA/B*, *tmm*, *dmsA*, *mtgB* or *pyl*.

**DRAM2 v2.0.0-beta42**, the four topic sheets (`--sum_topics default`; pass it explicitly, because I did not confirm what an empty `sum_topics` resolves to):
- *phnG–M* including *phnJ*;
- *grdI*, *grdC/D*, *grdB/E*;
- *mtgB/mtgA*, the transporters, *mcr*, *mtt/mtb/mtm*, *mta*.

The **methylotrophy sheet** (`methylotrophy_distillate.tsv`) is shipped but loaded only with `--sum_custom`. It is the only sheet with:
- *grdH*, *grdA*, *cutC/D*, *cntA/B* ("YeaW/YeaX … QA cleavage to TMA"), *tmm*, *dmd-tmd* and *pylBCDS*;
- 48 rows that need DRAM2's separate `methyl` database (METHYL_0001 …, e.g. MtgC, RamQ, MtqA).

No DRAM2 sheet has *phnN/O/P*, *phnX/W/Y/Z*, *mpnS*, *torZ* or *tdm*.

**Deliverable:** `custom_distillate_salterns.dram1.tsv` (`DRAM.py distill --custom_distillate`) and `custom_distillate_salterns.dram2.tsv` (`--sum_custom`). They hold the 65 markers that have a KOfam threshold, as one "Salterns methane" sheet grouped like the table.

## 3. Local test of DRAM's KO layer (`kofam_control_test/`)

**Method.**
- DRAM's KOfam step was re-implemented: `prodigal -p meta`, contigs ≥2,500 bp, `hmmsearch --domtblout` against the 67 marker profiles, DRAM's `sig_scores` rule, and the best KO per gene by E-value.
- Profiles are the KOfam 2026-05-28 build in `~/Documents/EmilyKin/kofam_new/`, with its matching `ko_list`. K28979 isn't in that build.
- Controls:
  - **38 original 2020 MAGs**, matched gene-by-gene to their RAST CDSs. These are every RAST guild member in `guilds_r232.csv` (22 MePn consumers, 9 betaine reducers, 1 methanogen, 4 TMA-route candidates), plus Halomonas.1, Desulfobacter.2 and the archaea Haloquadratum.4 and Natrinema.1.
  - **8 re-assembly bins**, each at ≥2,500 and ≥1,000 bp:
    - two MaxBin2 *Roseovarius*;
    - the best new *Methanolobus*;
    - all 5 bins in Albimonas.1's species cluster.
- 186,714 proteins in total. The 15-MAG first pass (57,250 proteins) took 61 CPU-s.

**Results:**
- `guild_retention.tsv`: per guild member, RAST rule vs KOfam.
- `genes_vs_rast.tsv`: every marker hit or RAST marker CDS.
- `control_summary.tsv`: RAST vs KOfam per MAG and marker.
- `min_contig_effect.tsv`: KOfam calls in the new bins at ≥2,500 vs ≥1,000 bp.

| finding | evidence |
|---|---|
| **MePn: 21 of 22 retained** | KOfam calls *phnJ* in 21 of the 22 RAST MePn consumers: 18 with *phnGHIJKLM*, Marivita.2 and Haloquadratum.4 without *phnH*, and Phaeodactylibacter.1 with a lone *phnJ*. KOfam also finds the archaeal PhnG that RAST missed (Haloquadratum.4, Natrinema.1). It calls no *phnJ* that RAST lacked. |
| **Albimonas has no callable *phnJ* in any assembly** | **Albimonas.1:** the locus is on a 5.4 kb contig; *phnI* is split (342 + 163 aa); *phnJ* (371 aa, 37 % HMM coverage) scores 134 vs 247 (2026-05 domain threshold; the 2026-10 threshold is 195, full sequence). **Re-assembly:** the same fragment, with identical gene lengths, recurs in `concoct.bin.157_sub` (R1 sediment co-assembly, MEGAHIT) and `metabat2.bin.062` (R2_restored_C_black, metaSPAdes). The best bin, `concoct.bin.058_sub` (95.2/2.4), has *phnG*, a single 494-aa *phnI*, *phnK*, *phnL* (226 vs 228), *phnN* and *phnM* on one contig, but no gene in the *phnH/J* position. RAST labelled Albimonas.1's fragment C-P lyase. |
| **grdH (Sec) fails in most betaine reducers** | Of the 9 RAST betaine reducers, KOfam *grdH* finds 3, *grdI* 7, and *grdH* or *grdI* 7. The shortened GrdH (332–350 aa, cut at the Sec codon) scores 602–606 vs 619 and is called *grdB* (Halanaerobium.1, M55B127.1). DATFJG01.1 passes at 625, with one full-length 437-aa copy at 770. *grdI* scores 796–840 vs 714. Halanaerobium.3 (76 % complete) and BM714.4 (25 % contaminated) have only *grdH/grdI* fragments; their *grdA/C/D* are present. |
| **K14083 (*mttB*) over-calls outside methanogens** | Called in **all 16 Rhodobacteraceae** controls and in 23 of the 37 non-methanogens, e.g. Roseovarius.8 6 genes (scores 423–560) and J153.1 7. RAST calls these "Trimethylamine methyltransferase family protein", the non-pyrrolysine MttB family. They still pass the newer 218.9 threshold. In UBA8670.2, M55B127.1/.2 and Desulfotignum.1/.2 the same family scores as *mtgB* (K25220, e.g. 556–601 vs 536), the betaine demethylase. |
| **Pyl genes are found as fragments** | Methanolobus.1: *mtbB* (375 aa) passes; both *mtmB* fragments pass, which doubles the copy count; the *mttB* N-terminal fragment passes at 229 (2026-05 threshold 161; 2026-10 threshold 219, marginal). The C-terminal pieces don't. *mcrA*, *mtbA*×2, *mttC/mtbC/mtmC* and *pylBCD* all pass. Pyl genes (pyrrolysine, UAG codon) reach KOfam only as fragments. |
| **TMA routes RAST missed** | **CntA/YeaW** (K22443, no CntB) in Roseovarius.2, 5, 7, 8, 10 and 11, e.g. 452 vs 305; RAST calls it "Rieske 2Fe-2S domain protein". J153.1's lone CntA is confirmed (396). **CutC** (K20038) in Desulfotignum.1 at 1,031 vs 990, a gene RAST calls "Pyruvate formate-lyase"; its activating enzymes score 329–338, below the 396 *cutD* threshold. Glycyl-radical enzymes are easily confused, so treat this as a candidate. |
| **2,500 bp minimum** | **R1_A_D2** (new MaxBin2 *Roseovarius*): 2,512 vs 3,972 genes; *proV* 2→0, *proX* 2→0, *proW* 1→0, *opuBD* 2→1, *opuC* 1→0, *mttB* 5→4. **R1_B_D2:** 2,058 vs 3,872 genes; *proW* 1→0, *proX* 2→1, *mttB* 2→0. The *phn* loci were on long contigs in both. The new *Methanolobus* (CONCOCT) is unaffected. |

What this means:
- **The weakest part of any automated guild call is Sec/Pyl genes and broken loci**, and Bakta and RAST call genes the same way, so they face the same problem.
  - Bakta's ≥80 % coverage requirement for a PSC hit may push the shortened *grdH* and the Pyl pieces down to product-only PSCC hits. The poplar test should check this.
- **A KO-based distillate will mis-assign *mttB***. Report *mttB/mtbB/mtmB* as a methylamine-methanogenesis capacity only in genomes with *mcrA* or with the pyrrolysine machinery (*pylBCDS*), or after checking for an in-frame UAG.

## 4. What happens downstream

**Metabolic models (ModelSEED/KBase).** The existing models were built from RAST, and the ModelSEED template is keyed to SSO terms. For models on the new MAGs that are comparable with the 2020 ones, **run RAST and add DRAM/KOfam/Bakta terms as extra annotation events**. Merge Metabolic Annotations followed by Build Metabolic Model can work from KO/EC events alone, but the resulting reaction set would not be directly comparable with the 2020 models.
- KBase `MergeMetabolicAnnotations/import_annotations` (release ec971d11) takes a 2-column TSV of `gene_id → term` for EC, KO, RO, META, MSRXN, GO or SSO. The gene IDs must be the KBase Genome's own.
  - Then use `merge_metabolic_annotations`, then Build Metabolic Model.
  - So run DRAM or KOfam **on the RAST proteins**: `DRAM.py annotate_genes`, DRAM2 `--input_genes` (KBUtilLib `DRAM2Utils` already does this), or `bakta_proteins`.
  - Otherwise map the calls by coordinates. Bakta's `--regions` can also force RAST's CDS coordinates.
- **KO events would not add the reactions the paper hinges on.**
  - `cb_annotation_ontology_api/data/kegg_95_0_ko_seed.tsv` (KEGG 95, 2020) maps no *grd*, *cut*, *cnt*, *mtgB* or transporter KO to a reaction.
  - RAST's betaine-reductase SSOs map to none either, which is why `build_betaine_models.py` adds rxn17220 by hand.
  - **EC 1.21.4.4 does map** (rxn07209, rxn11678, rxn17220, rxn33055). Bakta emits EC only on PSC/IPS hits. DRAM carries the KO's EC in `kegg_hit` (DRAM1) or `kofam_EC` (DRAM2).
  - For *phnJ*, KO K06163 → rxn40081, while RAST's SSO → rxn26459 (the SAM-dependent version). Importing KOs next to RAST would create a duplicate C-P lyase reaction.

**Submission (GenBank/ENA).**
- Annotation is optional for MAGs: GenBank can annotate with PGAP on request ([NCBI metagenome guide](https://www.ncbi.nlm.nih.gov/genbank/metagenome/)), and ENA takes MAG assemblies through Webin-CLI ([ENA MAG docs](https://ena-docs.readthedocs.io/en/latest/submit/assembly/metagenome/mag.html)).
- Bakta's `--compliant` GFF3/GBFF/EMBL is checked against table2asn and Webin-CLI, but it's for bacteria.
- DRAM 1 writes per-MAG `genbank/` and `genes.gff` that are not INSDC-validated. DRAM2 ≥beta38 no longer writes a GFF from called genes (changelog afd40b1); from FASTA it keeps prodigal GFFs.
- No KBase app runs Bakta (release catalog scanned 2026-11-11). `kb_DRAM` 935498b6 runs DRAM 1.4.0 with 2022-09 databases and a 2,500 bp minimum.

**Metabolism figure.** DRAM's `product.html` and `metabolism_summary` are the main thing Bakta can't give. They are only useful here with the custom sheet and the *mttB* caveat from section 3.

**Consistency with wetland.**
- `wetland/dram_all_mags_runbook.md` picks **DRAM 1.5.0**: KOfam, Pfam, dbCAN, MEROPS, no UniRef, 2,500 bp minimum, for parity with KBase's DRAM 1.4.0 and narrative cell C's DRAM-1 columns.
- It argues *against* DRAM2.
- The salterns memory note, though, says **you run DRAM2 on h100**:
  - `/scratch1/fliu/hub_scratch/chenry/DRAM2`, pinned 2026-06-18;
  - Nextflow 24.10.5, whereas current DRAM2 needs ≥25.04, so h100 is an older beta;
  - default DBs kofam, dbcan, merops, vog; no UniRef, KEGG or Pfam, and `methyl` not listed.
- The KO calls are comparable either way: both use KOfam with the same rule. What differs is the summary sheets and column names.
- Choose one for both projects if a cross-project figure is planned. Either way, record the KOfam release.

## 5. Cost

| | per MAG | 558 best-per-species (1,622 Mb, ~1.66 M genes) | all 1,324 (3,724 Mb, ~3.8 M genes) | database / where |
|---|---|---|---|---|
| **KOfam only** (DRAM's KO layer; hmmsearch, local) | ≈0.6 CPU-s/protein, so ≈0.5 CPU-h for a 2.9 Mb MAG | ≈270 CPU-h ≈ **11 h on 24 threads** | ≈630 CPU-h ≈ **26 h on 24 threads** | KOfam 1.5 GB tar; a 2026-05-28 copy is already at `~/Documents/EmilyKin/kofam_new/` |
| **DRAM on KBase** (measured, wetland) | 21 min + 4 min/Mb at 30 threads, serial | ≈300 node-h | ≈710 node-h | KBase queue (unreliable, see memory notes) |
| **DRAM 1.5 local** | KOfam above plus Pfam/dbCAN/MEROPS; <50 GB RSS/process | ~1–2 days on 24 threads (estimate) | ~2–4 days | ~30 GB DBs; setup ~128 GB RAM; not installed here or on poplar |
| **DRAM2 on h100** | unmeasured | the user's existing plan | 2.3× | installed (kofam, dbcan, merops, vog) |
| **Bakta 1.12 full** (measured, poplar) | 134 MAGs (median 3.5 Mb) in 2 h 50 min at 10 jobs × 6 threads, ≈12.7 slot-min/MAG on a loaded host | ≈100 slot-h, ≈10 h at 10 jobs | ≈220 slot-h, ≈22 h at 10 jobs | full DB 30 GB download / 84 GB unpacked; on poplar (EmilyKin env `544db00d…`); the `light` DB gets no KO/EC from UniRef |
| **eggNOG-mapper 2.1.13** | unmeasured; needed `--dmnd_iterate no` and the DB in /dev/shm | | | on poplar |

**How the numbers were derived.**
- The KOfam figure scales the measured 61 CPU-s for 67 profiles × 57,250 proteins by total profile length (16.6 M vs 29.5 k positions, ×~560). Treat it as ±2×.
- Genes: 1,022 per Mb in the controls.
- Bakta: the EmilyKin log `reprocess_291/bestofboth/bakta_revert.log`, scaled by Mb. Novel MAGs make Bakta slower, because fewer proteins are exact matches and more go to DIAMOND.

## 6. Recommendation

**Annotation:**
1. **KO calls for the guilds: KOfam on all 1,324 MAGs, contigs ≥1,000 bp.**
   - Run it here: Prodigal plus HMMER from `~/Documents/gtdbtk_venv`, `nice` and ≤24 threads, about one day. Or run DRAM2's kofam module.
   - Use 1,324 rather than 558 because the C-P lyase is an accessory, mobile operon that can differ within a 95 % ANI species. Copying the representative's call to the whole cluster would hide that, and guild abundance is computed per MAG.
   - Use the RAST proteins instead of Prodigal if KBase models will be rebuilt, so that there is one gene-ID space.
2. **DRAM2 on h100 for the 558 best-per-species MAGs**, with `custom_distillate_salterns.dram2.tsv`. This gives the metabolism figure. Settings:
   - `--min_contig_len 1000 --anno_dbs kofam,dbcan,merops`, plus `methyl` if installed;
   - `--summarize --sum_topics default --sum_custom custom_distillate_salterns.dram2.tsv`.
   - If the wetland project stays on DRAM 1.5, use DRAM 1.5 here too, with `--min_contig_size 1000` and the `.dram1.tsv` sheet.
3. **Guild rules for the new MAGs:**
   - betaine reducer = *grdH* **or** *grdI* (K21579 | K21578), with *grdA/C/D* as support;
   - MePn consumer = *phnJ* (K06163). Report the *phnGHIJKLM* completeness next to it and flag loci split across contigs. Resolve *Albimonas* first: tblastn reference PhnH/PhnJ against its locus in Albimonas.1 and `concoct.bin.058_sub` to tell a frameshifted pseudogene from a divergent gene;
   - methanogen = *mcrA*;
   - methylamine use = *mttB/mtbB/mtmB* in an *mcrA* genome. Bacterial K14083/K25220 is "MttB family (quaternary-amine demethylase?)", not TMA methanogenesis;
   - TMA producer = *cutC* + *cutD*, or *cntA* + *cntB*. The Desulfotignum.1 CutC and the *Roseovarius* CntA-without-CntB are candidates to check, not members.
4. **Bakta: not now.** Run it, or request PGAP, only for the MAG set you deposit with annotation.
   - It is the better tool for INSDC files, 16S extraction and product names.
   - It is not designed for the 116 archaeal MAGs, and its KO yield on novel taxa is the open question for the poplar test.
   - If you do run it, use the full database, not `light`.
5. **Before choosing between Bakta + eggNOG and DRAM2,** run the 8-MAG test below once poplar is reachable.

## 7. The 8-MAG test on poplar

**Inputs.** Run each tool on the same contigs, at ≥1,000 bp:
1. Bakta 1.12 full DB, default and `--meta`;
2. eggNOG-mapper `--dmnd_iterate no`;
3. full KOfam, run with `hmmsearch` on Bakta's `.faa` (the EmilyKin method);
4. DRAM2 on h100 or DRAM 1.5.

For the controls, use the original MAG and its re-assembly counterpart (`../comparison/mags_old_vs_new.tsv`, `species_cluster`), so that RAST is ground truth:

| taxon | original (RAST truth) | why |
|---|---|---|
| *Roseovarius* HQ | Roseovarius.8 (99.1/3.1) | C-P lyase reference; 6 MttB-family genes |
| *Albimonas* | Albimonas.1 (99.0/2.4) + new `concoct.bin.058_sub.fasta__R2_sediment` (95.2/2.4) | *phn* locus without a callable *phnJ* in any assembly (KOfam ≠ RAST) |
| *Methanolobus* | Methanolobus.1 (98.4/2.6) + new `concoct.bin.101.fasta__R1_sediment` | Pyl *mttB/mtbB/mtmB*, *mcr*; archaeon |
| *Halomonas* | Halomonas.1 (97.0/5.2) | betaine/osmolyte transporters, *gbcA/B* |
| *Marivita* | Marivita.2 (99.1/0.7) | second Rhodobacteraceae MePn consumer |
| Halanaerobiaceae | Halanaerobium.1 (88.3/11.8) | Sec *grdH* fails KOfam threshold |
| Desulfobacterota | Desulfobacter.2 (99.2/2.7) | MttB-family false positive, *dsr* context |
| haloarchaeon | Haloquadratum.4 (83.5/3.7) | archaeal C-P lyase, PhnG missed by RAST |

Add M55B127.1 (Sec *grdH* plus *mtgB*) and the two MaxBin2 *Roseovarius* bins if cheap.

**Measure, per tool and per MAG:**
1. CDS count and N50 of called genes. Also: do gene coordinates match RAST and Prodigal, and can they be mapped to KBase gene IDs?
2. **Fraction of CDS with any KO**, split by tool:
   - Bakta `KEGG:` cross-reference;
   - eggNOG `KEGG_ko`;
   - KOfam.

   Report the Bakta fraction separately for named-genus and placeholder-genus MAGs. Also report the fractions with EC, with a non-hypothetical product, and Bakta's split into IPS/UPS, valid PSC (≥90 %) and PSCC-only.
3. **The 68 markers from `marker_ko_table.tsv` compared with RAST truth:** found, missed and extra calls per tool. Specifically:
   - the shortened Sec *grdH* (does Bakta give a KO, a PSCC product, or a hypothetical?);
   - the Pyl pieces of *mttB/mtbB/mtmB*;
   - the split *phnI/phnJ* in Albimonas.1;
   - *mtgB* versus *mttB* naming in bacteria;
   - *mpnS/phnY\*/phnZ1*, which have no threshold;
   - h100's `kofam_ko_list.tsv` row for K18049 (`-` or blank).
4. Whether eggNOG's 5.0.2 vocabulary has K21578/K21579, K22443/K22444, K25220 and K23996.
5. For archaea (Methanolobus, Haloquadratum): Bakta's tRNA/rRNA counts versus `tRNAscan-SE -A`/barrnap `--kingdom arc`, and whether `--compliant` output passes `table2asn` and Webin-CLI.
6. Wall time, CPU time, peak RSS and output size per MAG. Use these to rescale section 5.
7. Using DRAM2's output: does `--sum_custom custom_distillate_salterns.dram2.tsv` run, and does the K14083 caveat show up in the heatmap as expected?

## Open items

- poplar access (Duo) for section 7. Also check the EmilyKin Bakta DB version (`refdata/bakta/db/version.json`) and the eggNOG DB version there.
- On h100: the DRAM2 revision, whether the `methyl` and `pfam` DBs are installed, and the K18049/K23996/K05915 rows of `kofam_ko_list.tsv`.
- Which DRAM version for both projects (wetland runbook = DRAM 1.5; salterns note = DRAM2).
- Whether RAST is to be run on the new MAGs for models (KBase queue) before or after the KO layer.
- The new MAGs' 2,500 bp question also depends on the separate MaxBin2 re-binning decision. If MaxBin2 is re-run at 2,500 bp, DRAM's default minimum stops mattering.
- **Albimonas as a MePn consumer** (manuscript-level). No assembly gives a KOfam-callable *phnJ*. Settle it with a tblastn/frame check of the locus.
- The control test used the 2026-05-28 KOfam build. Re-score with the current release before using any number in Methods, because thresholds moved (section 2a).

## Files

| file | what |
|---|---|
| `marker_ko_table.tsv` | 68 markers: KO, EC, KOfam thresholds (2 releases), DRAM1/DRAM2 distillate rows, KO/EC→ModelSEED, guild column |
| `custom_distillate_salterns.dram1.tsv`, `.dram2.tsv` | custom distillate sheets (65 markers with a threshold) |
| `short_contig_fraction_per_mag.tsv` | per MAG (1,324): contigs and bp below 2,500 bp, binner, species cluster, best-of-cluster flag |
| `kofam_control_test/` | `run_kofam_controls.py`, `summarize_controls.py`, `guild_retention.tsv`, `genes_vs_rast.tsv`, `control_summary.tsv`, `marker_calls_by_mag.tsv`, `min_contig_effect.tsv`, `summary.txt` |
| `build_marker_table.py`, `markers_input.tsv` | rebuild the table and sheets from `sources/` |
| `sources/` | cached evidence:<br>• KEGG KO list 2026-11-11 and KOfam `ko_list` 2026-10-26 and 2026-05-28;<br>• DRAM 1.5.0 sheets and `utils.py` (master, last commit 2025-05-15);<br>• DRAM2 v2.0.0-beta42 sheets, `distill.py`, `hmm_parser.py`, `nextflow.config` and CHANGELOG (dev, 2026-11-09);<br>• Bakta v1.12.2 `psc.py`, `pscc.py`, `constants.py`, `annotation.py` and README (main; DB v6.0) |

**Other sources:**
- KEGG REST: <https://rest.kegg.jp/list/ko>, <https://rest.kegg.jp/link/ko/eco:b1802+eco:b1803>.
- KOfam: <https://www.genome.jp/ftp/db/kofam/>.
- DRAM: <https://github.com/BortonWrightonLabs/DRAM> (DRAM1 = `master`, DRAM2 = `dev`).
- Bakta: <https://github.com/oschwengers/bakta>; Schwengers et al. 2021, *Microb. Genom.* 7:000685.
- eggNOG-mapper v2: Cantalapiedra et al. 2021, *MBE* 38:5825.
- The KBase release app catalog (NarrativeMethodStore `list_methods`, read 2026-11-11).
- `~/Documents/cb_annotation_ontology_api/data/` (`kegg_95_0_ko_seed.tsv`, `EC_translation.tsv`, `SSO_reactions.json`).
- `~/Documents/ModelSEEDDatabase/Biochemistry/reaction_0*.tsv`.
- `~/Documents/EmilyKin/` (`reprocess_291/`, `MANUSCRIPT_methods_results_discussion.md`).
- `~/Documents/KBUtilLib/src/kbutillib/dram2_utils.py`.
- `wetland/dram_all_mags_runbook.md`.
