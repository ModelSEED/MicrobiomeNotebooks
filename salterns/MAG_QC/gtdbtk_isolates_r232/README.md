# GTDB R232 reclassification of the two *Halomonas* isolates (I1, I4)

Done 2026-10-26 in response to review comment 3 ("re-classify the two isolates against R232; GTDB has split *Halomonas*").
The draft (P47) states that the isolates were classified with GTDB-Tk v1.7.0 in KBase narrative 139158.

## Bottom line

| Isolate | Name in draft (P67, Fig. 5/6) | NCBI deposit (BioProject PRJNA992430) | GTDB R232 (GTDB-Tk 2.8.0) | C-P lyase (*phnG–M*) |
|---|---|---|---|---|
| **I1** (`data/I1_SPAdes.Assembly.fa`) | *Halomonas* sp. D (Fig. 6); "*Halomonas* sp. (isolate)" (Fig. 5) | GCA_032334135.1 / GCF_032334135.1, "*Halomonas* sp. I1", strain I1 | g__*Halomonas*; s__*Halomonas* sp032334135 | **present**, full *phnCDEFGHIJKLMNP* cluster |
| **I4** (`data/I4_SPAdes.Assembly.fa`) | *H. saccharevitans* | GCA_032333995.1 / GCF_032333995.1, "*Halomonas saccharevitans*", strain I4 | g__*Halomonas*; s__*Halomonas saccharevitans* | absent |

- **Both isolates stay in *Halomonas*.** The R232 split does not affect either isolate's name at genus level, and the draft's "*Halomonas*" wording for the isolates (P47, P48, P67, P93, P98, P122, P123) can stay.
- **"*Halomonas* sp. D" = I1**, the C-P lyase carrier. This is confirmed by gene content (below) and agrees with the Figure 6 legend, the Figure 5 heat map, and the GTDB-Tk result for I4.
- **The changes are elsewhere:**
  - Figure 5 references: 11 of its 23 *Halomonas*-labelled reference genomes are in other genera in R232, which affects the P67 sentence "Six other *Halomonas* genomes ... while 17 did not".
  - The methods version string (P47).
  - The naming of I1. The draft calls it "sp. D", but NCBI calls it "sp. I1".

## 1. GTDB-Tk 2.8.0 / R232 run

- **How it was run:** `run_gtdbtk_isolates.sh` (repo root) ran `gtdbtk classify_wf` against reference data r232 (`~/Documents/gtdbtk_data/release232`). The log is `gtdbtk.log`; the result is `gtdbtk.bac120.summary.tsv`.
- **Classification method:** for both genomes it was `ani_screen` ("classification based on ANI only"). Both were assigned by skani ANI to a species representative, without pplacer placement.
- **I1:**
  - Closest reference: GCF_032334135.1, ANI 100.0, AF 0.999.
  - **This is a self-match, not an independent identification.** GCF_032334135.1 is I1's own deposited assembly. It has the same BioProject PRJNA992430 and strain I1, and its total length of 3,937,892 bp equals that of `data/I1_SPAdes.Assembly.fa`.
  - **What it means:** R232 holds the I1 genome as the representative of its own species cluster, s__*Halomonas* sp032334135, in g__*Halomonas*. No species with a validly published name is within the species radius.
  - The only other reference GTDB-Tk reports is GCF_029081005.1, s__*Halomonas elongata*_B, at 92.7% ANI (AF 0.757), below the 95% radius. I1 is therefore a distinct *Halomonas* species with no published name.
- **I4:**
  - Closest reference: GCF_900116405.1, the *H. saccharevitans* type strain CGMCC 1.6493 (NCBI: "assembly from type material"), at ANI 98.49 and AF 0.912.
  - I4's own deposit, GCF_032333995.1, is also in GTDB. The GTDB API taxon history gives s__*Halomonas saccharevitans* in both R226 and R232. It is a cluster member, not the representative, so it is absent from the reps-only taxonomy in the local GTDB-Tk data.
  - The next closest species are *H. pelophila* (93.2% ANI) and *Halomonas* sp043496965 (92.5% ANI).
- **I1 history:** the GTDB API taxon history for GCF_032334135.1 is s__*Halomonas* sp032334135 in both R226 and R232.

## 2. Which isolate is "*Halomonas* sp. D"? (gene content)

- **Method:**
  - Genes were predicted with Prodigal V2.6.3 (`~/Documents/gtdbtk_venv/bin/prodigal`, single mode). This gave 3,639 proteins for I1 and 3,105 for I4.
  - The proteins were searched with `hmmsearch --cut_tc` (HMMER in the same venv) against two HMM sets:
    - Pfam PhnG (PF06754.13), PhnH (PF05845.13), PhnI (PF05861.13) and PhnJ (PF06007.12).
    - TIGRFAMs 15.0 phosphonate models: PhnC TIGR02315, PhnD TIGR03431, PhnE TIGR01097, TIGR01098, PhnF TIGR02325, PhnG TIGR03293, PhnH TIGR03292, PhnK TIGR02323, PhnL TIGR02324, PhnM TIGR02318, PhnN TIGR02322, PhnP TIGR03307, and TIGR03223.
  - Both HMM sets were taken from the local copy at `~/Documents/emilykin_tools/gtdbtk_r220_markers/` (`Pfam-A.hmm.gz`, `tigrfam_full.hmm`).
- **I1:** one contiguous, same-strand cluster on contig `NODE_6_length_145494` (bp 68,112–78,719, prodigal genes 6_73 to 6_85) in the order *phnC, phnD, phnE, phnF, phnG, phnH, phnI, phnJ, phnK, phnL, phnM, phnN, phnP*.
  - PhnJ is gene 6_80, a complete ORF of 277 aa (834 nt including the stop codon). It scores 481.7 bits against PF06007 (E = 9.7e-146).
  - I1 also has a second PhnCDE-family transporter locus (`NODE_9`, genes 24–27).
- **I4:** no hits above trusted cutoffs to PhnF, G, H, I, J, K, L, M, N or P. The C-P lyase is absent.
- **Conclusion:** I1 = *Halomonas* sp. D (C-P lyase carrier, CH₄ from MPn in Figure 6). I4 = *H. saccharevitans* (no C-P lyase, no CH₄).
- **Caveat for anyone comparing with Figure 5:** I4 does have PhnC-family (TIGR02315 ×2) and PhnE-family (TIGR01097 ×4) ABC-transporter genes above trusted cutoffs, in two loci on `NODE_1` (genes 537–539 and 701–703). It also has three TIGR01098 binding-protein hits: two at the start of those `NODE_1` loci (genes 536 and 700) and the third on `NODE_3` (gene 258). It has no TIGR03431 (PhnD) hit. A later check with TIGR01422 (PhnX, phosphonoacetaldehyde hydrolase) and TIGR00686 (PhnA) also finds *phnX* (`NODE_1` gene 788) and *phnA* (`NODE_7` gene 30) above trusted cutoffs in I4; neither is part of the C-P lyase, so this does not change the conclusion. Figure 5, which is RAST-based, scores this isolate as lacking *phnC/D/E*. These TIGRFAM families also include phosphite and other phosphonate-like transporters, so the hits are not evidence of MPn transport. They are an annotation-scope difference, not an error in the figure. No change is proposed.
- **Repo state:** no isolate RAST/KEGG annotation was found in the repo; the KBase narrative is not accessible (section 4).

## 3. NCBI records (NCBI Datasets v2 API, fetched 2026-10-26)

| | GCF_032334135.1 (I1) | GCF_032333995.1 (I4) | GCF_900116405.1 (closest ref. of I4) |
|---|---|---|---|
| Organism | *Halomonas* sp. I1 | *Halomonas saccharevitans* | *Halomonas saccharevitans* |
| Strain | I1 | I4 | CGMCC 1.6493 (type material) |
| BioProject / BioSample | PRJNA992430 / SAMN37482663 | PRJNA992430 / SAMN37482664 | PRJEB17565 / SAMN04487956 |
| Submitter, release | LBNL, 2023-10-03 | LBNL, 2023-10-03 (same BioProject) | DOE JGI, 2016-11-02 |
| Length / contigs | 3,937,892 bp / 40 scaffolds, 53 contigs (scaffold level) | 3,406,698 bp / 10 scaffolds, 18 contigs (scaffold level) | 3,681,289 bp / 94 |
| BioSample attributes | soil; "USA: San Francisco"; collected 2021-10-07; pure culture | soil; "USA: San Francisco"; collected 2021-10-07; pure culture | geo_loc_name "missing" |

- **The BioProject holds exactly these two genomes.** The PRJNA992430 dataset report lists only I1 and I4 (GCA and GCF of each).
- **Lengths match the repo assemblies.** Both deposits match the repo assemblies in total length: 3,937,892 and 3,406,698 bp. The repo FASTA files have 40 and 10 sequences, the same as NCBI's 40 and 10 scaffolds. NCBI's 53 and 18 contigs are those scaffolds split at their N runs: the repo FASTAs have 13 and 8 runs of at least 10 N, and 40 + 13 = 53, 10 + 8 = 18.
- **NCBI's CheckM for I1** is 99.89% complete and 0.75% contaminated.
- **The draft wording at P47 needs a small correction.** It says the sequences are "available on NCBI GenBank via accession PRJNA992430", but that is a BioProject. The assembly accessions are GCA_032334135.1 and GCA_032333995.1.
- **Naming:** the draft calls I1 "*Halomonas* sp. D", but NCBI's organism name is "*Halomonas* sp. I1". Figure 5 labels it "*Halomonas* sp. (isolate)". The authors should use one name, or state the equivalence once.

## 4. How GTDB R232 splits *Halomonas*

- **Data used:** the local R232 taxonomy (`~/Documents/gtdbtk_data/release232/taxonomy/gtdb_taxonomy.tsv`). It holds **species representatives only**: its 199,923 lines equal the 189,801 bacterial plus 10,122 archaeal representative lines. The counts below are therefore species clusters, not genomes or validly named species.
- **f__Halomonadaceae has 342 species clusters in 25 genera.** The largest are:
  - g__*Halomonas*: 117;
  - g__*Vreelandella*: 101;
  - *Salinicola*: 18;
  - *Modicisalibacter*: 17;
  - *Marinospirillum*: 13;
  - *Oceanospirillum*: 12;
  - *Kushneria* and *Chromohalobacter*: 11 each;
  - *Onishia*, *Franzmannia* and *Cobetia*: 6 each;
  - *Litchfieldella*: 3, *Litchfieldella*_A: 1;
  - *Larsenimonas* and *Bisbaumannia*: 3 each;
  - smaller genera.
- **There is no "*Halomonas*_A" or other suffixed *Halomonas* genus in R232** (`grep -c g__Halomonas_` returns 0, and likewise 0 for `g__Vreelandella_`). The split is into separately named genera.
- **The type species stays in *Halomonas*:** *H. elongata* (GCF_000196875.2) is in g__*Halomonas*.
- **Examples that moved:** *H. venusta* and *H. aquamarina* are now *Vreelandella venusta* and *V. aquamarina*.
- **NCBI has adopted many of these transfers.** It now names, for example, *Vreelandella boliviensis*, *V. stevensii*, *V. jeotgali*, *V. subterranea*, *Litchfieldella anticariensis*/*xinjiangensis* and *Modicisalibacter luteus*.
- **Literature (unverified details):** the transfers follow the taxogenomic revision of the family by de la Haba et al. (2023, *Front Microbiol*; citation details unverified).
- **Both isolates are in g__*Halomonas* (confirmed above):** I1 = s__*Halomonas* sp032334135 and I4 = s__*Halomonas saccharevitans*.

## 5. Figure 5 reference genomes (50) mapped to R232

- **Why the mapping is by name.** KBase narrative/workspace 139158 is not readable by this account:
  - `list_objects` and `get_workspace_info` fail with "User freiburgermsu may not read workspace 139158";
  - `get_permissions_mass` returns only `{'freiburgermsu': 'n'}`, with no public `'*': 'r'` entry.

  So the SpeciesTree object, with its genome refs, could not be read.
- **Inference: the narrative does not appear to be public.** That conflicts with P47's "publicly accessible with narrative ID 139158".
- **How the labels were mapped:**
  1. The 50 reference labels were read off the Figure 5 image in the draft (`word/media/image2.png`): 52 tips = 2 isolates + 50 references.
  2. Each label was resolved to RefSeq assemblies with the NCBI Datasets API, matching on strain where the label gives one.
  3. Each assembly was mapped with the GTDB API taxon history (R232).
  4. KBase's public RefSeq workspace (ReferenceDataManager, ws 19217) was used to check ambiguous names.
- **The 23 *Halomonas*-labelled references match P67's numbers.** There are 23 *Halomonas*-labelled references besides the two isolates, which matches P67's "six ... while 17".

*Halomonas*-labelled references ("MPn genes" = multiple *phn* transport and C-P lyase genes in Figure 5, as read from the image):

| Figure 5 label | Assembly (NCBI) | Current NCBI name | GTDB R232 | MPn genes |
|---|---|---|---|---|
| *H. daqiaonensis* | GCF_900109725.1 | *H. daqiaonensis* | *Halomonas daqiaonensis* | no |
| *H. salina* | GCF_000758175.1 (the KBase "*Halomonas salina*" object) | *H. halophila* | *Halomonas halophila* | no (*phnP* only) |
| *H. saccharevitans* | GCF_900116405.1 (type) | *H. saccharevitans* | *Halomonas saccharevitans* | no |
| *H. shengliensis* | GCF_900104135.1 | *H. shengliensis* | *Halomonas shengliensis* | no |
| *Halomonas* sp. PBN3 | GCF_000475415.1 | *Halomonas* sp. PBN3 | *Halomonas shengliensis* | **yes** |
| *H. halodenitrificans* DSM 735 | GCF_000620045.1 | same | *Halomonas halodenitrificans* | no (*phnC* only) |
| *H. korlensis* | GCF_900116705.1 | same | *Halomonas korlensis* | no (*phnC* only) |
| *H. smyrnensis* AAD6 | GCF_000265245.1 | same | *Halomonas smyrnensis* | no |
| *H. halocynthiae* DSM 14573 | GCF_000425725.1 | same | *Halomonas halocynthiae* | no |
| *H. huangheensis* | GCF_001431725.1 / GCF_000470745.1 (same strain) | same | *Halomonas huangheensis* | **yes** |
| *Halomonas* sp. S2151 | GCF_000967645.1 | same | *Halomonas litopenaei* | **yes** |
| *Halomonas* sp. KM-1 | GCF_000246875.1 | same | *Halomonas aerodenitrificans* | no (*phnC* only) |
| *H. stevensii* S18214 | GCF_000275725.1 | *Vreelandella stevensii* | ***Vreelandella** stevensii* | no |
| *H. boliviensis* LC1 | GCF_000236035.1 / GCF_002265845.1 | *Vreelandella boliviensis* | ***Vreelandella** boliviensis* | no |
| *H. campaniensis* | GCF_002211105.1 (SK03), GCF_000696485.1 (LS21) **or** GCF_014193375.1 (type, 5AG) | *H. campaniensis* | SK03: ***Vreelandella** campaniensis*_A; LS21: ***Vreelandella** alkaliphila*; 5AG: *Halomonas campaniensis* | no |
| *Halomonas* sp. GFAJ-1 | GCF_000236625.2 | same | ***Vreelandella*** sp002966495 | no |
| *H. subterranea* | GCF_900111305.1 | *Vreelandella subterranea* | ***Vreelandella** subterranea* | no |
| *Halomonas* sp. HL-93 | GCF_900086985.1 | same | ***Vreelandella*** sp900086985 | no (*phnC* only) |
| *H. jeotgali* Hwa | GCF_000334215.1 | *Vreelandella jeotgali* | ***Vreelandella** jeotgali* | **yes** |
| *H. xinjiangensis* | GCF_000759345.1 | *Litchfieldella xinjiangensis* | ***Litchfieldella**_A xinjiangensis* | no |
| *H. anticariensis* FP35 | GCF_000409775.1 / GCF_000428505.1 | *Litchfieldella anticariensis* | ***Litchfieldella** anticariensis* | no |
| *H. zincidurans* B6 | GCF_000731955.1 | *H. zincidurans* | ***Modicisalibacter** zincidurans* | **yes** |
| *H. lutea* DSM 23508 | GCF_000378505.1 | *Modicisalibacter luteus* | ***Modicisalibacter** luteus* | **yes** |

**Summary: 12 of the 23 stay in *Halomonas*, and 11 move.** The 11 are 7 to *Vreelandella*, 2 to *Litchfieldella*/*Litchfieldella*_A and 2 to *Modicisalibacter*. These counts hold for two of the three candidate *H. campaniensis* genomes (see the inference below).

| Figure 5 group | Still *Halomonas* in R232 | Moved to another genus |
|---|---|---|
| Six with MPn genes | 3: PBN3, *huangheensis*, S2151 | 3: *jeotgali* → *Vreelandella*; *zincidurans* and *lutea* → *Modicisalibacter* |
| 17 without | 9 | 8: 6 → *Vreelandella*, 2 → *Litchfieldella* |

- **Inference: the *H. campaniensis* genome.** KBase's public RefSeq set (ws 19217) holds three genomes named "*Halomonas campaniensis*": SK03 (GCF_002211105.1), LS21 (GCF_000696485.1, a complete genome) and the type strain 5AG (GCF_014193375.1). SK03 and LS21 are both *Vreelandella* in R232 (*V. campaniensis*_A and *V. alkaliphila*); only 5AG is *Halomonas*. Which one the tree used therefore rests on tree topology. In Figure 5, *H. campaniensis* sits inside the clade with *H. stevensii*, *H. boliviensis*, GFAJ-1, *H. subterranea* and HL-93, with *H. jeotgali* as sister; its sister tip is GFAJ-1 (*Vreelandella* sp002966495). All of these are *Vreelandella* in R232. That placement fits SK03 or LS21 (both *Vreelandella*), not the *Halomonas* type-strain genome. The topology supports the genus only: it cannot tell SK03 from LS21, so the tip's species is not identified.
  - If the KBase tree object can be opened, check this genome ref. If it turns out to be 5AG, the counts become 13 stay and 10 move (17 without: 10 *Halomonas*, 5 *Vreelandella*, 2 *Litchfieldella*).
  - The six-with-MPn-genes breakdown does not depend on this.
- **Other references.** All are Halomonadaceae or other Pseudomonadales genera that keep their genus, apart from GTDB-suffixed genera and renamed species:
  - *Chromohalobacter salexigens* DSM 3043 (GCF_000055785.1) is *C. israelensis* in R232 and NCBI.
  - "*C. japonicus*" resolves in NCBI to *C. beijerinckii*. All 7 of the 8 candidate genomes that are in R232 are g__*Chromohalobacter*.
  - Gammaproteobacteria bacterium MFB021 (GCA_000737985.1) is *Salinicola* sp000737985 in R232.
  - *Marinobacter hydrocarbonoclasticus* ATCC 49840 (GCF_000284615.1) is *M. nauticus*.
  - *Marinobacter nanhaiticus* is *Marinobacter*_A.
  - *Marinobacterium* sp. AK27 (GCF_000705555.1) is *Marinobacterium*_D *lacunae*, and *Marinobacterium litorale* is *Marinobacterium*_D.
  - The rest (*Aidingimonas*, *Halotalea*, *Carnimonas*, *Zymobacter*, *Kushneria*, *Marinobacter*, *Nitrincola*, *Oceanospirillum*, *Marinospirillum*) map to their genus in R232.
- **P67 needs rewording.** "Six other *Halomonas* genomes also contained multiple MPn transport and degradation genes, while 17 did not" stays true of the labels in Figure 5. Under R232, though, only three of the six (PBN3 = *H. shengliensis*, *H. huangheensis*, S2151 = *H. litopenaei*) are *Halomonas*; the other three are *Vreelandella* or *Modicisalibacter*.
- **Figure 5 itself:** a caption note covers the change. Relabelling the tree requires the KBase SpeciesTree object, which is not accessible.

## 6. MAG-level *Halomonas* statements (task item 5)

- **The R232 MAG changes** (`correlational/iterativeIDs_r232_vs_2020.csv`, `correlational/r232/guilds_r232.csv`):
  - *Halomonas*.1 (`R1_A_D1 ... metabat.26`) stays *Halomonas* (s__*Halomonas fontilapidosi*). It is 97.0% complete with 5.21% contamination, and its mean relative abundance over the 12 cores is 8.2 × 10⁻⁴.
  - 2020 *Halomonas*.2 (`R1_B_H2O ... metabat.5`) is now ***Vreelandella*.1**, with no species assigned. It is 92.27% complete with 9.27% contamination, and its mean relative abundance is 8.3 × 10⁻⁵.
- **The original draft makes no MAG-level statement about *Halomonas*.** Every *Halomonas* mention concerns the isolates or the genus as a whole.
- **Figure S6 (P144) is unaffected.** Its metagenomic panel comes from IMG/M contig mapping to the IMG reference database, not from the MAGs.
- **The R232 review docx already reports the MAG change.** `make_r232_review_docx.py` adds a Results paragraph saying one of the two *Halomonas* MAGs is now *Vreelandella*. It also uses "*Halomonas*/*Vreelandella*" in the new co-occurrence text and the Figure 8 caption.
- **Context for the "case study" framing (P67, P93, P98):**
  - Neither Halomonadaceae MAG carries *phnJ* or any accessory *phn* gene (the `ev_phnJ` and `ev_phn_accessory` evidence columns in `guilds_r232.csv` are empty).
  - Neither isolate is represented among the 450 MAGs. `skani dist` of I1 and I4 against all 450 `MAG_contigs/*.fna` returned a single pair above skani's screening threshold: I4 vs *Halomonas*.1, at 85.65% ANI with AF of about 25%. That is not the same species. I1 matched nothing.
  - So the C-P lyase *Halomonas* is supported by the isolate experiment, not by the metagenome-assembled genomes (inference from these two checks).

## 7. Manuscript edits (applied)

These edits were applied as tracked changes by `clauded_insights/make_r232_review_docx.py` to `Metabolic Manuscript - R232 review (tracked changes).docx`; they are listed, with rationales, in `clauded_insights/r232_review_edits_2026-10-26.json` (origin "halomonas", plus the orchestrator additions noted below).
- **P47:**
  - The original clause "classified taxonomically with GTDB-Tk v1.7.0 (68)" is kept, because its citations are hyperlinks and cannot be placed inside a tracked deletion. A new sentence is inserted after the KBase sentence: the isolate assemblies were later reclassified with GTDB-Tk v2.8.0 against GTDB release R232, outside KBase.
  - A comment notes that reference 68 is the GTDB-Tk v1 paper (Chaumeil et al. 2020), so GTDB-Tk v2 needs its own citation. The orchestrator anchored this comment at P47 instead of P245, because P245's text is all hyperlinks.
  - The clause "available on NCBI GenBank via accession PRJNA992430" is replaced to give the BioProject and the assembly accessions with strain names.
  - A comment notes that the narrative is not public.
  - These replace the P47 comment that was previously pending in `make_r232_review_docx.py`; the builder now points to the edits JSON for P47.
- **P67:** the isolates' R232 placement and the R232 status of the six MPn-positive references are added, plus comments on the naming and on the MAG evidence.
- **P122:** a caption note says that 11 of the 23 *Halomonas*-labelled references are in other genera in R232. The orchestrator also added a caption comment that these counts assume a *Vreelandella* genome for the *H. campaniensis* tip and gives the 5AG alternative (10 of 23). That comment should name the tip as SK03 or LS21 (both *Vreelandella* in R232), not SK03 alone (section 5).
- **P93 and P173:** comments only.
- **No change:** P48, P98, P123.

## Reproduce

Intermediates (prodigal proteins/GFF, HMM tables, NCBI/GTDB API JSON, skani table) were kept in session scratch only. To regenerate:

```bash
B=~/Documents/gtdbtk_venv/bin
$B/prodigal -i data/I1_SPAdes.Assembly.fa -a I1.faa -q -o /dev/null      # same for I4
# PhnG/H/I/J models extracted from Pfam-A.hmm.gz; TIGR models from tigrfam_full.hmm (see section 2)
$B/hmmsearch --cut_tc --tblout I1_phn.tbl pfam_phn.hmm I1.faa
$B/skani dist -q data/I1_SPAdes.Assembly.fa data/I4_SPAdes.Assembly.fa -r MAG_contigs/*.fna
curl "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/bioproject/PRJNA992430/dataset_report"
curl "https://gtdb-api.ecogenomic.org/genome/GCF_032333995.1/taxon-history"
```
