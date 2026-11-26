# bakta_test: Bakta vs Prodigal gene calls for the KO layer (9 re-assembly MAGs, poplar)

*Track of the 2026-11-11 annotation runs. Question: should the KOfam calls for the 1,324 re-assembly MAGs use Bakta's CDS
(to get a full-length selenocysteine GrdH) instead of Prodigal? Test design: `../../annotation_choice/README.md` section 7.*

**Status (2026-11-11 12:00 CST): RUNNING on poplar.** All 9 Bakta runs finished OK (last at 11:56; 30–45 min each at load
~500 on 336 cores). Prodigal is done (53,219 Bakta + Prodigal proteins). KOfam (27,981 profiles in 39 chunks) started at
11:59 and is writing output. eggNOG DIAMOND and the eggNOG-DB vocabulary grep are running. Expected end ~13:00–14:00
(`PIPELINE_DONE`). *(Reviewer, 12:08: KOfam was at ~7,700 of 27,981 profiles after 8.3 min, so the KOfam step should
end ~12:30–12:45; `PIPELINE_DONE` then waits for eggNOG-mapper. See "Review" at the end of this file.)* After that, run `scripts/analyze_bakta_test.py` by hand (command below). The Sec-GrdH question is already answered (Results 1); the remaining outputs fill the KO-yield and
naming tables. Note on the poplar login: the ControlMaster (`ssh -MN -o ControlPersist=yes`, pid 2604712 on
frankensteinsmonster) forked to the background after the Duo login; "poplar connection ended" in tmux was the foreground
wrapper exiting. BatchMode ssh has worked throughout.

## Where things are

| what | path |
|---|---|
| poplar work dir | `/scratch1/afreiburger/salterns_2026/bakta_test/` |
| KOfam DB (this track downloaded it; shared with the kofam track) | `/scratch1/afreiburger/salterns_2026/kofam_db/` (`ko_list`, `profiles/K*.hmm`, `md5sums.txt`, `README`) |
| local copy of scripts, MAG table, results | this folder (`scripts/`, `mags.tsv`, `grdscan_ranking.tsv`, `results/`) |

## Versions

| tool / DB | version |
|---|---|
| Bakta | 1.12.0 (EmilyKin env `544db00dd5c254ecfa7c6335967a046f_`; `annotation_choice/sources/` holds 1.12.2 source) |
| Bakta DB | v6.0 **full**, 2025-02-24 (`/scratch1/afreiburger/emilykin/processed/refdata/bakta/db`, `version.json`) |
| Pyrodigal (inside Bakta) | 3.7.1 |
| Prodigal | 2.6.3 (same env; `-p meta -g 11`, as DRAM and the kofam track) |
| HMMER | 3.4 (Aug 2023) |
| tRNAscan-SE / Infernal / AMRFinderPlus | 2.0.12 / 1.1.5 / 4.2.7 |
| eggNOG-mapper | 2.1.13, eggNOG DB 5.0.2 (files dated 2020-11/2021-03), DIAMOND 2.2.0 |
| **KOfam** | **genome.jp release: `ko_list.gz` Last-Modified 2026-10-26 05:26 GMT, `profiles.tar.gz` 2026-10-27 01:41 GMT** (**27,981 HMM profiles** searched; the `profiles/` directory has 27,983 entries because it also holds `prokaryote.hal` and `eukaryote.hal`; 28,504 ko_list rows, 523 of them without a profile, e.g. the rRNA KOs K01977–9; all 68 marker KOs have a profile and a ko_list row). md5 `ko_list.gz` 697af5b3b22a919306a1e60f99cec1a6, `profiles.tar.gz` bc68ad2a52ed48812b5496d15c215a7d. The unzipped `ko_list` is byte-identical (md5 5dca1ddf…) to `annotation_choice/sources/kegg/kofam_ko_list_2026-10-26.gz`. |

KO rule (identical to `annotation_choice/kofam_control_test/run_kofam_controls.py`, i.e. DRAM's `sig_scores`): `hmmsearch --domtblout`
with default reporting thresholds; a hit is significant when the score selected by `score_type` (full-sequence or per-domain)
exceeds the `ko_list` threshold; KOs with threshold `-` are never called; each gene gets the significant KO with the lowest
full-sequence E-value. Thresholds that matter here (2026-10): grdH K21579 624.03, grdI K21578 505.90, grdB K10672 277.97,
mttB K14083 218.93, mtgB K25220 536.00, phnJ K06163 194.63 (all full-sequence).

## The MAGs (`mags.tsv`)

| label | bin (`reassembly_2026/work/checkm/<sample>/fasta/<bin>.fna`) | why |
|---|---|---|
| Roseovarius_mb121 | metabat2.bin.121_sub.fasta__R1_sediment (96.1/2.8) | best-quality new *Roseovarius* |
| Albimonas_cc058 | concoct.bin.058_sub.fasta__R2_sediment (95.2/2.4) | Albimonas.1 cluster; *phn* locus without a callable *phnJ* |
| Methanolobus_cc101 | concoct.bin.101.fasta__R1_sediment (98.4/2.6) | best HQ *Methanolobus*; Pyl *mttB/mtbB/mtmB*; archaeon |
| Halomonas_mb336 | metabat2.bin.336.fasta__R1_sediment (83.5/2.8) | best new *Halomonas* (Halomonas.1 cluster) |
| Marivita_cc101 | concoct.bin.101_sub.fasta__R2_restored_C_black (92.8/1.7) | best new *Marivita* (Marivita.2 cluster) |
| Haloquadratum_mx640 | maxbin2.Bin.640_sub.fasta__R2_sediment (91.5/3.3) | haloarchaeon, Haloquadratum.4 cluster (archaeal C-P lyase) |
| Halanaerobium_mb206 | metabat2.bin.206.fasta__R1_sediment (97.2/2.4) | Halanaerobiaceae with a Sec *grdH*: Prodigal's 350-aa GrdH ends at TGA, frame reads on 88 codons, K21579 623.6 vs 624.03 |
| Desulfo_M55B175_mb289 | metabat2.bin.289.fasta__R2_sediment (99.4/2.6) | Desulfobacterota: **no Desulfobacterota bin carries grdH/grdI** (pre-scan below), so the best one with *grd* genes: *grdA* + a Sec *grdB*-like gene + 4 MttB-family genes (one passes *mtgB*, 539.1 vs 536) |
| M55B127_mb089 | metabat2.bin.089.fasta__R2_sediment (97.2/1.8) | extra, per section 7 ("add M55B127.1"): Sec *grdH* 599.3 + 6 *mtgB* |

All contigs in these bins are ≥1,000 bp (min 1,001–2,614 bp), so no length filter was applied; Bakta and Prodigal see identical input.

### Selection pre-scan for *grd* genes (`grdscan_ranking.tsv`)

`kofam_control_test` had no new Desulfobacterota or Halanaerobiaceae bins, so 184 new bins (all 118 Desulfobacterota + 66
Bacillota in Halanaerobiaceae, M55B127, Acidaminobacteraceae, BM714, Dethiosulfatibacteraceae) were scanned on poplar:
Prodigal `-p meta -g 11` (543,724 proteins) and `hmmsearch` with 8 KOfam 2026-10 profiles (K21578, K21579, K10670, K10672,
K21576, K21577, K25220, K14083). Script `scripts/grdscan.sh`, parser `scripts/parse_grdscan.py`. Bins were ranked by raw
score (not pass/fail), and the stop codon of the best *grdH* gene was checked in the contig.

- **No Desulfobacterota bin has grdH or grdI above threshold** (113 with any hit; best grdI 251 vs 506). 4 pass *mtgB*; 34
  pass *mttB* (e.g. *Desulfotignum* metabat2.bin.079: 17 MttB-family genes, up to 610).
- **The Sec-truncated GrdH is the normal case in the Bacillota.** Of the 41 bins whose best *grdH* hit scores ≥500, 35
  have a best *grdH* gene of 323–353 aa that ends at **TGA** with the frame continuing for 43–88 codons (Prodigal calls the
  C-terminal ~65–79 aa as a separate gene). Their scores, 538–635, straddle the 624.03 threshold: only 2 pass (the two
  BM714 bins, 635.3); Halanaerobium metabat2.bin.206 fails by 0.4. 32 of the 35 have a passing *grdI*.
- A **full-length 431–437-aa grdH (TAA stop)** in 5 bins (UBA8670, DATFJG01) scores 770–784 and passes: no Sec codon in
  these copies.

## Commands (all on poplar, inside the work dir)

```bash
ANNOT=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_
export PATH=$ANNOT/bin:$PATH
# KOfam download (16 parallel byte ranges; genome.jp gives ~0.4 MB/s per connection)
bash /scratch1/afreiburger/salterns_2026/kofam_db/download.sh
# Bakta, per MAG, default options (scripts/run_bakta.sh; 4 then 2 jobs in parallel)
bakta --db /scratch1/afreiburger/emilykin/processed/refdata/bakta/db --threads 8 --skip-plot \
      --tmp-dir tmp/bakta_<label> --output bakta/<label> --prefix <label> fasta/<label>.fna
# Prodigal (DRAM / kofam-track settings)
prodigal -p meta -g 11 -q -i fasta/<label>.fna -a prodigal/<label>.faa -f gff -o prodigal/<label>.gff
# KOfam: all 27,981 profiles in 39 chunk files, one hmmsearch --cpu 1 per chunk, against Bakta+Prodigal proteins
hmmsearch --cpu 1 --domtblout kofam/out/chunk_NN.domtblout -o /dev/null kofam/chunks/chunk_NN.hmm kofam/all_proteins.faa
# eggNOG-mapper on Bakta proteins (one job)
emapper.py --data_dir /scratch1/afreiburger/emilykin/processed/refdata/eggnog -m diamond --itype proteins \
    -i eggnog/bakta_all.faa -o bakta_test --output_dir eggnog --temp_dir eggnog/tmp --cpu 8 --dmnd_iterate no --dbmem --override
```

`scripts/pipeline.sh` chains everything after Bakta batch 1 (waits on done-markers, so it survives this agent ending).
Every step is wrapped in `/usr/bin/time -v` (`logs/*.time`). Thread budget ≤48: Bakta 4×8 + pre-scan 16, then Bakta 2×8
next to batch 1's tail, then KOfam 39×1 + eggNOG 8 + one `grep` of the eggNOG DB for KO vocabulary.

## Processes (poplar)

| job | PID | log / done marker |
|---|---|---|
| KOfam download | 3495089 (finished 10:49) | `kofam_db/download.log`, `kofam_db/.download_done` |
| grd pre-scan | 3537061 (finished 10:50) | `grdscan/grdscan.log` (`GRDSCAN_DONE`) |
| Bakta batch 1 (6 MAGs) | 3501139 (finished 11:55, all OK) | `logs_bakta_batch1.progress` (`BAKTA_DONE`), `logs/bakta_<label>.log`, `bakta/<label>/<label>.log` |
| orchestrator (batch 2 Bakta → Prodigal → KOfam + eggNOG) | 3646472 (Bakta batch 2 runner 3772073 finished 11:56, all OK) | `pipeline.log`, `logs_bakta_batch2.progress`, `PIPELINE_DONE` |
| follow-up: the 35 read-through GrdH proteins vs all 27,981 profiles (starts after `PIPELINE_DONE`) | 4138936 | `followup.log`, `FOLLOWUP_DONE`, `kofam/sec_joined.domtblout.gz` |
| eggNOG-DB KO vocabulary `grep` (child of 3646472) | 787783, reniced 0 → 10 by the reviewer at 12:06 (it was the only un-niced process of this track) | `eggnog/eggnog_db_ko_vocab.txt` (`VOCAB_DONE`) |

Check status:

```bash
ssh -o BatchMode=yes poplar 'cd /scratch1/afreiburger/salterns_2026/bakta_test; cat pipeline.log logs_bakta_batch*.progress; ls PIPELINE_DONE; tail -n 2 logs/emapper.log'
```

Re-run the analysis locally once `PIPELINE_DONE` exists:

```bash
R=/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026
D=<scratch dir>; rsync -a -e "ssh -o BatchMode=yes" --exclude 'fasta' --exclude 'grdscan' --exclude 'tmp' --exclude 'kofam/out' --exclude 'kofam/chunks' \
   --exclude 'kofam/all_proteins.faa' poplar:/scratch1/afreiburger/salterns_2026/bakta_test/ $D/
rsync -a -e "ssh -o BatchMode=yes" poplar:/scratch1/afreiburger/salterns_2026/kofam_db/ko_list $D/
~/Documents/py_venv/bin/python -I $R/annotation_runs/bakta_test/scripts/analyze_bakta_test.py $D $R/annotation_runs/bakta_test/results \
   $R/annotation_choice/marker_ko_table.tsv
```

## Results

### 1. Bakta recovers a full-length Sec GrdH, and KOfam then calls it grdH (interim, from Bakta's CDS stage)

Bakta 1.12 has a selenocysteine step (`bakta/features/cds.py: revise_translational_exceptions`). It joins two
neighbouring Pyrodigal ORFs into one CDS with `U` only when **all** of these hold: same strand and frame; the upstream ORF
ends in TGA; the gap is <100 nt; an Rfam **SECIS** element (ncRNA-region search, `cmscan --cut_tc`, E ≤1e-4) lies between
them; and the joined translation starts with M and has exactly one internal stop. Otherwise the gene stays split, as in
Prodigal.

- **Halanaerobium_mb206:** Bakta found one SECIS_4 (RF01990, score 57.9, E 2.4e-10) at 39,217–39,254, directly adjacent
  to the TGA (39,214–39,216), and emitted a
  439-aa GrdH with U at position 351 (`exception: selenocysteine`). Against the 2026-10 profiles:

  | protein | length | K21579 grdH (thr 624.03) | K10672 grdB (thr 277.97) | KO called |
  |---|---|---|---|---|
  | Prodigal N-terminal piece | 350 aa (ends at TGA) | 623.6 (cov 0.80) | 471.4 | **grdB** (wrong) |
  | Prodigal C-terminal piece | 79 aa | 135.3 | 100.6 | none |
  | **Bakta Sec CDS** | **439 aa, U351** | **808.9** (cov 1.00) | 623.3 | **grdH** |

  So with Bakta's CDS, KOfam passes the grdH threshold by 185 points; with Prodigal it misses by 0.4 and mislabels.
- It is **the only U-containing protein in the 7 finished Bakta runs**; M55B127_mb089's log shows a second Sec join of its
  GrdH (contig_78 702–2,021, U at the TGA 969–971). *grdA* (Sec near the N-terminus) stays a
  108-aa C-terminal piece in both callers (it still passes *grdA*, 175 vs 55). The Desulfobacterota bin's Sec *grdB*-like
  gene (345 aa) is not joined: no SECIS was found there.
- **Generalisation (the 35 Sec-truncated grdH loci of the pre-scan):** `cmscan` of a 2 kb window around each TGA with
  Bakta's own `ncRNA-regions` DB and settings found a **SECIS_4 hit directly after the TGA in all 35** (score 52.9–65.6,
  E 7.8e-12 to 3.3e-09 at -Z 5; `grdscan/sec_windows.fna`, one serial `cmscan --cpu 0`, 4 s, run next to the 48 Bakta
  threads), so Bakta's SECIS condition is met everywhere; whether Pyrodigal also calls the
  C-terminal piece (the other condition) is checked per MAG only where Bakta ran.
- **A Prodigal-side rescue works as well:** translating each of the 35 Prodigal grdH genes through the TGA (as U) to the
  next in-frame stop gives 394–440-aa proteins that score **716–823 on K21579; all 35 pass** (vs 2 of 35 as Prodigal
  pieces), and K21579 becomes the best KO in all 35 (`scripts/sec_readthrough.py`, `results/interim/sec_readthrough_grdH_scores.tsv`).

### 2. Marker calls, Bakta CDS vs Prodigal (interim; 68 marker profiles)

From the Bakta proteins sent to its PSC search (all CDS except exact UPS matches; complete for the novel MAGs, missing
781 exact-match CDS in Halomonas) vs Prodigal on the same contigs (`results/interim/marker_calls_prodigal_vs_bakta_interim.txt`):
- Every marker count is identical in Albimonas, Roseovarius, Methanolobus and the Desulfobacterota bin, including the
  *phn* locus and the MttB-family over-calls (Roseovarius 6, Albimonas 2, Desulfo 3 *mttB* + 1 *mtgB*).
- Bakta has no pyrrolysine read-through. In Methanolobus both callers split *mttB* into the same 333 + 124 aa pieces
  (only the 333-aa piece is called, 304.6) and *mtbB* into 354/375 + 105 aa (Bakta picks a later start); *mtmB* is one
  457-aa gene in both (568.0).
- The only real difference is Halanaerobium **grdH 0→1, grdB 1→0** (section 1). The Halomonas differences are the
  missing UPS-matched CDS in the interim file, not caller differences.

### 3. Bakta's database hits are absent where the MAGs are novel (interim, Bakta logs)

*(Reviewer: these CDS counts come from Bakta's log at the PSC stage, before overlap filtering. The final CDS counts in
Bakta's `.txt` and JSON are in table 4: Roseovarius 3,224, Halomonas 2,624, Albimonas 3,865.)*

| MAG | GTDB | CDS | exact UPS/IPS | PSC (≥90 % id, gives KO/EC) | PSCC (50–90 %, product only) |
|---|---|---|---|---|---|
| Roseovarius_mb121 | Roseovarius sp003298505 | 3,228 | 215 | 2,599 | 263 |
| Halomonas_mb336 | Halomonas fontilapidosi | 2,633 | 781 | 1,414 | 253 |
| Albimonas_cc058 | Albimonas sp035268645 | 3,867 | 1 | **59** | 3,083 |
| Methanolobus_cc101 | Methanolobus sp020055655 | 2,202 | 0 | **0** | 716 |

### 4. Seven finished Bakta runs: KO yield, runtime and names (`results/interim/wave1_bakta_*.tsv`, `scripts/wave1_names.py`)

| MAG | CDS | Bakta `KEGG:` xref | EC | non-hypothetical | exact UPS/IPS | PSC (≥90 %) | PSCC only | no hit | wall | CPU | peak RSS |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Roseovarius_mb121 | 3,224 | **26.4 %** | 23.7 % | 92.9 % | 6.7 % | 80.6 % | 7.6 % | 4.6 % | 36.1 min | 96.5 min | 8.1 GB |
| Halomonas_mb336 | 2,624 | **26.8 %** | 23.6 % | 90.8 % | 29.8 % | 53.8 % | 8.4 % | 6.8 % | 36.6 min | 90.8 min | 9.0 GB |
| Marivita_cc101 | 3,353 | **24.3 %** | 20.7 % | 88.0 % | 62.3 % | 23.4 % | 5.5 % | 7.4 % | 30.3 min | 96.5 min | 8.0 GB |
| Albimonas_cc058 | 3,865 | **0.6 %** | 0.4 % | 77.1 % | 0 | 1.5 % | 76.3 % | 18.7 % | 40.2 min | 133.2 min | 8.2 GB |
| Halanaerobium_mb206 | 2,094 | **2.3 %** | 1.1 % | 82.0 % | 0.2 % | 5.8 % | 77.2 % | 12.2 % | 32.1 min | 131.7 min | 11.5 GB |
| Desulfo_M55B175_mb289 | 2,872 | **0.3 %** | 0.1 % | 71.9 % | 0 | 1.1 % | 72.4 % | 22.0 % | 33.2 min | 133.2 min | 8.0 GB |
| Methanolobus_cc101 | 2,202 | **0 %** | 0 % | 16.4 % | 0 | 0 | 16.5 % | 67.5 % | 44.7 min | 135.9 min | 8.3 GB |

("no hit" includes the 0.6–16 % with only a sub-50 % PSC alignment.) Haloquadratum_mx640 and M55B127_mb089 were
still running at 11:55.

Even with an 80 % PSC rate (Roseovarius, a GTDB species with a reference), only about a quarter of the CDS get a KO,
because only ~18 % of UniRef90 clusters carry one. In every MAG without a close reference (Albimonas, Halanaerobium,
M55B175, Methanolobus) Bakta's KO layer is essentially empty: 0–2.3 % of CDS.

Bakta's names for the paper's loci:
- **Albimonas *phn* locus** (contig k127_5730352): PhnG (98 aa) "Alpha-D-ribose 1-methylphosphonate 5-triphosphate
  synthase subunit PhnG"; **the 494-aa PhnI is "hypothetical protein" with no database hit at all** (KOfam phnI 300 vs
  170); a 72-aa hypothetical; PhnK; then a 321-aa "Magnesium transport protein CorA" between *phnK* and *phnL*; PhnL;
  PhnN; PhnM (328 aa, "…5-triphosphate diphosphatase"), followed by a 66-aa "amidohydrolase family protein". All named
  ones are PSCC (UniRef50 product names, no KO). **Bakta also finds no gene in the *phnH/phnJ* position**, so the
  Albimonas question is not a gene-caller problem.
- **Roseovarius *phn* locus**: *phnFGHIJKLMN* all PSC-level. PhnJ is "Alpha-D-ribose 1-methylphosphonate 5-phosphate C-P
  lyase" with KEGG K06163, so Bakta and KOfam agree here.
- **MttB-family proteins**: Bakta calls all of them **"Methyltransferase"** (no gene symbol). In Roseovarius 5 of 6 carry
  Bakta's own **KEGG K14083** xref (PSC), so Bakta propagates the same *mttB* over-call that KOfam makes. The two in
  Albimonas are PSCC "Methyltransferase" without a KO.
- **Methanolobus (archaeon): Bakta calls McrA, McrB, McrG, MtbB, MtmB, MttC, MtbC, MtbA, MtaA and PylB/C/D/S all
  "hypothetical protein" with no database hit.** Only the 333-aa MttB piece gets a name ("Trimethylamine
  methyltransferase", PSCC). KOfam calls all of them (mcrA 983 vs 741, pylS 718 vs 652, …).
- **Betaine reductase genes in Halanaerobium:** the Sec GrdH is "Betaine reductase complex component B subunit beta"
  (PSCC, so **no KO** from Bakta even though the protein is complete). GrdI is "Betaine reductase component B alpha
  subunit" (PSC, still no KO). GrdA is "glycine/sarcosine/betaine reductase complex selenoprotein A" (K10670, 108-aa
  C-terminal piece). Bakta also emits the 42-aa N-terminal GrdA piece as a separate CDS ("…component A2"), which
  passes *grdA* (64 vs 55), so Bakta + KOfam double-counts *grdA*. GrdC and GrdD are named, with no KO.
- **Desulfobacterota M55B175:** the four MttB-family genes are "Trimethylamine--corrinoid protein Co-methyltransferase"
  or "Trimethylamine methyltransferase" (PSCC, no KO), so Bakta's names repeat the TMA mislabel. The KOfam-*mtgB* gene is
  among them. The Sec *grdB*-like (345 aa) and *grdI/E*-like (446 aa) genes are "hypothetical protein".
- Other disagreements: Halomonas' 828-aa TMAO reductase is KOfam *torZ* (K07812, 1,246.6), Bakta "TorA" with xref
  K08351 (*bisC*). Roseovarius' three KOfam *dmd-tmd* (K00317) genes are Bakta "NADH:flavin oxidoreductase" with K00317.

### Final results

(pending: `PIPELINE_DONE` on poplar, then `scripts/analyze_bakta_test.py`. It was tested end-to-end on a real Bakta
1.12.0 JSON from EmilyKin, CAN_1_bin.107. It writes `results/per_mag_summary.tsv` (CDS counts, KO/EC/hypothetical
fractions, IPS/UPS vs PSC vs PSCC split, U-containing CDS, Bakta/Prodigal wall time, CPU time and peak RSS),
`marker_calls_PBXE.tsv` (per marker × MAG: Prodigal+KOfam / Bakta+KOfam / Bakta KEGG xref / eggNOG), and
`marker_genes_detail.tsv` (every marker-relevant Bakta CDS with Bakta's product/gene/xrefs, the KOfam and eggNOG calls and
the matching Prodigal gene(s)). Bakta's names for the Albimonas *phn* locus and for the MttB-family proteins come from
that last table. `eggnog/eggnog_db_ko_vocab.txt` answers whether eggNOG 5.0.2 contains K21578/9, K22443/4, K25220,
K23996 and the other markers. Sanity check: K06163 (*phnJ*) must be non-zero there. If it is 0, the sqlite fields are
compressed and the vocabulary has to be read from the `.emapper.annotations` output or an sqlite query. *(Reviewer: the
schema of `eggnog.db` was checked. `prots.kegg_ko` is a plain `VARCHAR(256)` text column, so the grep is valid. Its
numbers count occurrences across the whole file, not distinct proteins, so read them as present/absent.)*)

## Conclusion so far: keep Prodigal for the 1,324 MAGs and add a Sec read-through step for the *grd* genes

The question was whether the 1,324-MAG KO calls should use Bakta's CDS to get the selenocysteine GrdH. The four options,
scored on the 35 Sec-truncated *grdH* loci of the pre-scan:

| option | Sec GrdH called as grdH | cost on poplar | other effects |
|---|---|---|---|
| Prodigal + KOfam (current kofam track) | **2 / 35**; the 350-aa piece is called *grdB* | already running | none |
| Prodigal + KOfam, guild rule *grdH* **or** *grdI* (annotation_choice README) | **32 / 35** bins counted as betaine reducers | none | misses UBA8670 maxbin2.Bin.014 (grdI 416), M55B127 maxbin2.Bin.004 R1_C_D1 (321) and UBA8670 metabat2.bin.079 (no grdI); the *grdB* mislabel stays |
| **Prodigal + KOfam + Sec read-through** of TGA-ending *grd*-family genes | **35 / 35** (716–823 vs 624; K21579 the best of the 68 marker KOs in all 35) | seconds | fixes the *grdB* mislabel too |
| Bakta CDS + KOfam | SECIS condition met 35 / 35; full 439-aa U-GrdH emitted in **2 / 2** Bakta runs with a Sec grdH (Halanaerobium; M55B127, contig_78 702–2,021, U at the TGA 969–971) | whole Bakta run. Measured today at load ~500: **30–45 min wall, 91–136 CPU-min, 8–11.5 GB RSS per MAG** at 8 threads (7 MAGs) (the PSC DIAMOND search alone is ~25 min). With 6 jobs (48 threads), 1,324 MAGs is ≈6 days. EmilyKin's 12.7 slot-min per MAG, at lower load, gives ~22 h at 10 jobs. | Bakta's own KO layer is 0–2.3 % of CDS in the novel MAGs (24–27 % even in named species); in Methanolobus, *mcrABG*, *mtbB*, *mtmB* and *pylBCDS* are all "hypothetical protein"; no Pyl read-through; the Sec GrdH gets a name but no KO, so KOfam is still needed; Albimonas *phn* locus and MttB-family calls identical to Prodigal; gene IDs differ from the kofam track's Prodigal IDs |

Caveat on "35 / 35": the 35 loci come from 5 families (M55B127 21, Acidaminobacteraceae 5, Dethiosulfatibacteraceae 5,
BM714 2, Halanaerobiaceae 2), 6 genera (M55B127 21, DATFJG01 5, UBA8670 5, BM714 2, Halanaerobium 1, QPJE01 1) and 10
skani species clusters (`comparison/mags_old_vs_new.tsv`; the largest, Salt_Pond R1_C_D2 metabat.40 = M55B127.1, holds 12).
They are mostly conspecific bins carrying the same locus. The per-locus result is robust, but it rests on about 6 genera
(10 species), not 35 independent loci. *(Corrected by the reviewer: the earlier counts 21/9/7/2/2 were for the 41 bins
with grdH ≥500, which include the 6 bins with a non-Sec, TAA-ending grdH.)*

**Recommendation.** Do not switch the 1,324-MAG KO layer to Bakta for GrdH. Keep the kofam track's Prodigal genes and
add a targeted Sec read-through rescue:
1. take every Prodigal gene whose best hit at any score is a *grd* Sec family (K21579 grdH, K10672 grdB, K21584 grdF,
   K10670 grdA) and that ends in TGA with the frame open downstream;
2. translate through the TGA as U to the next in-frame stop and re-score;
3. accept the joined protein only if the downstream segment matches the C-terminal part of the same profile (the 79-aa
   Halanaerobium piece covers HMM 359–437) or a SECIS_4 lies within ~50 nt of the TGA (`cmscan` with Bakta's
   `ncRNA-regions` DB; found in 35 / 35 here);
4. keep *grdI* as supporting evidence in the guild rule.

Run Bakta only on MAGs you deposit with annotation. Its Sec join is a real advantage for INSDC files, but it does not
change the guild calls if the rescue above is applied. This conclusion does not depend on the pending eggNOG and full
KOfam numbers. Those will show how sparse Bakta's and eggNOG's KOs are on the novel MAGs, not change the GrdH result.

Follow-ups:
- `FOLLOWUP_DONE`: check that K21579 stays the best KO for the 35 joined proteins against all 27,981 profiles.
- The Desulfobacterota *grdB*-like gene (345 aa, grdB 391) has no SECIS_4 hit, so neither Bakta nor the SECIS-gated
  rescue joins it. Check by hand whether it ends in TGA and whether the HMM-gated rescue would join it.

## Review (2026-11-11 ~12:10 CST, reviewer agent)

**Verdict: the running work is correct. Nothing was killed or relaunched.** One process was reniced and the README
was corrected (marked *Reviewer* above).

Checked on poplar and locally:
- **Inputs:** md5 and contig counts of all 9 `fasta/<label>.fna` on poplar match
  `reassembly_2026/work/checkm/<sample>/fasta/<bin>.fna`. Minimum contig lengths are 1,001–2,614 bp, so every contig is
  ≥1,000 bp, as stated.
- **Scripts:** the md5s of `pipeline.sh`, `run_bakta.sh`, `followup.sh` and `grdscan.sh` on poplar equal the local copies in
  `scripts/`. The follow-up input `grdscan/sec_joined.faa` (35 seqs) is identical to
  `results/interim/sec_readthrough_grdH_joined.faa`.
- **KOfam release:** `ko_list` md5 5dca1ddf… on poplar equals the gunzipped local `kofam_ko_list_2026-10-26.gz`.
  Download size and `gzip -t` passed. 27,981 HMMs were searched. The six thresholds quoted above are in `ko_list` with
  `score_type` full. All 68 marker KOs have a profile.
- **KO rule:** `analyze_bakta_test.py` applies DRAM's `sig_scores` rule exactly as `run_kofam_controls.py` does: full or
  domain score per `score_type`, `-` never called, best significant KO by full E-value. It reads the domtblout columns
  correctly. Protein IDs in `kofam/out/*.domtblout` (`bakta|<label>|<locus>`, `prodigal|<label>|<contig>_<n>`) and in
  `eggnog/bakta_all.faa` (`<label>|<locus>`) match the keys the script builds. The Bakta 1.12.0 JSON keys it uses
  (`sequences[].orig_id`, `feature.sequence/locus/start/stop/strand`, `psc.valid`, `exception`) are present in this run's
  output. `psc.valid` is True for 121 of 2,094 Halanaerobium CDS (5.8 %, as in table 4). The E-values are comparable
  across the 39 chunks, because every chunk searches the same 53,219-protein target file.
- **Bakta:** v1.12.0 with DB v6.0 full. All 9 runs are OK. There are exactly 2 U-containing CDS, both 439 aa and named
  "Betaine reductase complex component B subunit beta": Halanaerobium LECAAC_00957 (Sec at 39,214–39,216, SECIS_4 at
  39,217–39,254) and M55B127 DJNIIK_01159 (contig_78 702–2,021, Sec at 969–971). The Halanaerobium Sec CDS carries only a
  UniRef50 xref (PSCC), so it has no KO, as stated.
- **Pre-scan and read-through numbers** (recomputed from `grdscan_ranking.tsv` and
  `sec_readthrough_grdH_scores.tsv`): 35 Sec loci; 2 pass grdH as Prodigal pieces; grdH-or-grdI gives 32; joined proteins
  score 715.7–823.1, all >624.03, with K21579 the best marker KO in all 35. No Desulfobacterota bin passes grdH or grdI
  (best grdI 250.8). Spot checks: Roseovarius has 6 *mttB* calls, 5 of them with a Bakta K14083 xref; Methanolobus
  *mcrABG*/*mtbB*/*mtmB*/*pylBCDS* are all "hypothetical protein" in Bakta.
- **Jobs and threads (12:01–12:08):** orchestrator 3646472 and follow-up 4138936 are alive. There are 39 KOfam
  `hmmsearch --cpu 1` at nice 10 (~68 % CPU each at load ~550) and eggNOG DIAMOND with `--threads 8` at nice 10. With
  the grep that makes 48 threads. KOfam was ~28 % done at 12:08. Host: 3 TB RAM, 2.3 TB available. The kofam track's 40
  `hmmsearch` (`hmm/chunk_NN`, `-T 14.7 -Z 3911146`) are a different process tree and were not touched.

**Fixes made**
1. `renice -n 10 -p 787783` (12:06). The eggNOG-DB vocabulary `grep` was the only un-niced process in this track.
2. README corrections:
   - Profile count: 27,981 HMMs (27,983 directory entries include the two `.hal` files).
   - The "35/35" lineage caveat used the 41-bin family counts. It now gives the 35-locus counts: 6 genera, 10 species
     clusters.
   - Table 3 CDS counts are pre-filter numbers from Bakta's log.
   - The eggnog.db grep is valid: plain-text `kegg_ko` column.
   - The vocab grep is now in the process table.

**Deviations from `annotation_choice/README.md` section 7. These were not fixed; for the user to decide.**
1. **Desulfobacterota pick.** Section 7 chose Desulfobacter.2 for the MttB-family false positive and *dsr* context, not
   for *grd* genes. Its re-assembly counterpart exists: `concoct.bin.049.fasta__R2_restored_C_black` (R2_restored_C_black,
   91.7/0.9, the same skani species cluster as Desulfobacter.2). It has 2 MttB-family genes passing *mttB* (best 551.7)
   and a RAST-truth original. The track's substitute, M55B175 `metabat2.bin.289`, still covers Bakta's MttB-family naming
   (3 *mttB* + 1 *mtgB*), but it is a different genus with no RAST truth. Adding concoct.bin.049 would cost about 35 min of
   Bakta at 8 threads, plus about 1 CPU-h of KOfam on about 6,000 proteins, plus eggNOG. Do it as a separate add-on with
   its own work subfolder and done marker, after `FOLLOWUP_DONE`. Do **not** append it to `labels.txt`:
   `analyze_bakta_test.py` reads that file and would fail on the missing `bakta/<label>/` and KOfam hits.
2. **No Bakta `--meta` arm** (section 7 input 1). Default Bakta trains Pyrodigal in single-genome mode, while the
   Prodigal comparison uses `-p meta`. So any CDS-count or start-site difference between them partly reflects the
   training mode, not Bakta itself. The marker calls are identical in 4 MAGs, so this is low-stakes for the GrdH
   decision.
3. **No original (RAST-truth) MAGs were run.** All 9 are re-assembly bins, so for Bakta, section 7 measure 3
   (found/missed/extra marker calls vs RAST, per tool) is replaced by Bakta-vs-Prodigal. The KOfam arm against RAST is
   already in `annotation_choice/kofam_control_test/`.
4. **Not attempted:** measure 5 (archaeal tRNA/rRNA vs `tRNAscan-SE -A`/barrnap `--kingdom arc`; `--compliant` output
   through `table2asn`/Webin-CLI) and measure 7 (DRAM2 on h100, outside poplar). Measure 1's mapping to KBase gene IDs
   is also not addressed.
5. Still pending as designed: `PIPELINE_DONE` → run `scripts/analyze_bakta_test.py` by hand (command above), then
   `FOLLOWUP_DONE`.

**eggNOG 5.0.2 KO vocabulary (section 7 measure 4; `eggnog/eggnog_db_ko_vocab.txt`, done 12:04).** All 19 KOs
grepped occur in `eggnog.db`. Counts are byte occurrences, not proteins: K21578 88, K21579 43, K10670 73, K10672 68,
K21576 88, K21577 85, K25220 72, K25221 58, K14083 2,536, K22443 97, K22444 180, K20037 104, K20038 73, K23996 87,
K18049 67, K16176 37, K16178 36, K05915 3, K06163 233 (sanity check passed). So eggNOG *can* emit the grd, mtgB, cutC/D
and cntA/B KOs. Whether it does so for these MAGs is in the final `marker_calls_PBXE.tsv` (eggNOG column).
