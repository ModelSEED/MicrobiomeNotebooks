# Alternative annotations for the 57 modelling representatives (KOfam, eggNOG-mapper, Bakta)

*Made 2026-11-20 for the pre-registered RAST-vs-merged comparison in `../ANNOTATION_CHOICE_CRITERIA.md` (arms B–D).
Genomes: the 57 `rep_genome` entries of `../../comparison_rollout/modelling_species.tsv` (31 rollout, 26 co-assembly).
Nothing was written to KBase and nothing was committed.*

**Status: DONE.** eggNOG-mapper ran on 57/57 genomes and Bakta on 57/57, with no failures. KOfam tables for 57/57 were
copied from the existing tracks. Every file collected from poplar passed an md5 check (1,056 files). Run on poplar
2026-11-20 18:51:37 → 21:18:02 UTC (2 h 26 min wall).

## Layout

| path | what |
|---|---|
| `kofam/<mag_id>.tsv.gz` | KOfam per-gene table (Prodigal genes), copied **byte for byte** from the KOfam tracks (source path + md5 in `kofam/sources.tsv`) |
| `kofam/sec_readthrough_57.tsv` | the Sec read-through rescue rows for these 57 genomes, with genomic ORF coordinates (see "KOfam" below) |
| `prodigal/<mag_id>.{faa,gff,genes.tsv}` | the KOfam tracks' Prodigal gene calls, also the eggNOG input (copied from poplar `kofam_rollout/genes/` and `kofam/genes/`) |
| `eggnog/<mag_id>/<mag_id>.eggnog_genes.tsv` | **every Prodigal gene with coordinates** + the eggNOG fields (empty = no annotation; `-` replaced by empty; `ko` = `KEGG_ko` without the `ko:` prefix) |
| `eggnog/<mag_id>/<mag_id>.emapper.{annotations,hits,seed_orthologs}` | that genome's rows of the pooled emapper outputs, original format (`.hits` got a column header line, DIAMOND `--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp`) |
| `eggnog/_pooled/g{1,2}/` | the two raw emapper runs (28 + 29 genomes) |
| `bakta/<mag_id>/` | raw Bakta output (`.json .tsv .gff3 .gbff .embl .faa .ffn .fna .txt .log .inference.tsv .hypotheticals.{tsv,faa}`) |
| `bakta/<mag_id>/<mag_id>.bakta_genes.tsv` | every Bakta feature, parsed from the JSON: locus_tag, type, **contig, start, end, strand**, aa_len, gene, product, hypothetical, pseudogene, truncated, evidence (UPS/IPS, PSC, PSCC, none), ko, ec, go, cog, uniref, pfam, exception (e.g. `selenocysteine@aa351`), db_xrefs |
| `coverage_per_genome.tsv` | the per-genome counts below (all columns) |
| `gene_caller_concordance.tsv` | Bakta CDS vs Prodigal genes per genome: shared stop codons and identical coordinates |
| `runtimes.tsv` | wall, CPU and peak RSS per job (GNU `time -v`) |
| `logs/` | poplar logs: `driver.log`, `monitor.log` (10-min checks of `/tmp`, load, progress), `bakta_<mag_id>.{log,time}`, `bakta_worker_*.progress`, `emapper_g*.{log,time}`, versions, `tmp_baseline.txt`, `fasta.md5` |
| `md5/` | `poplar_md5.txt` (computed on poplar) and `local_paths_md5.txt` (same, mapped to local paths; `md5sum -c` passes) |
| `scripts/` | `driver.sh` (poplar run), `collect.sh` (rsync + md5), `make_sec_table.py`, `summarize.py` (split, coordinate tables, counts), `mags57.tsv` (genome manifest), `bakta_order.txt`, `eggnog_group{1,2}.txt`, `summarize.out` |

poplar work dir: `/scratch1/afreiburger/salterns_2026/annot57/` (inputs `fasta/`, `prodigal/`; outputs `bakta/`, `eggnog/g1`,
`eggnog/g2`; `tmp/` is empty).

## Genomes and names

- `mag_id` is the part of `rep_genome` after `rollout|` or `coassembly|`. Every file here is named by `mag_id`.
- **FASTAs:** rollout from `/scratch1/afreiburger/salterns_2026/rollout/mags_kept/<mag_id>.fa`, co-assembly from
  `/scratch1/afreiburger/salterns_2026/coasm_qc/mags_kept/<mag_id>.fa`. They were copied to `annot57/fasta/` and
  md5-checked against the local `reassembly_2026/work/{rollout,coasm}_mags_kept/` copies.
- **These are the same files the KOfam tracks used.** All 57 md5s equal the md5s in `annotation_runs/kofam_rollout/mags.md5`
  (31 rollout) or `annotation_runs/kofam/mags.md5` (26 co-assembly). So KOfam, eggNOG and Bakta coordinates all refer
  to the same contigs. All contigs are at least 1,000 bp, so the KOfam track's ≥1 kb filter removed nothing.
- **Which KOfam track:** the 31 rollout genomes come from `kofam_rollout`. **All 26 co-assembly genomes come from the
  first-pass `kofam` track; none come from `kofam_coasm38`.** In the first pass, a co-assembly MAG
  `<unit>__<binner>.bin.N` was called `<binner>.bin.N.fasta__<unit>`. **Gene IDs are kept verbatim** (`<kofam_name>|<contig>_<n>`)
  in `kofam/`, `prodigal/` and `eggnog/`, so KOfam and eggNOG join on `gene` directly. The name map:

| mag_id (co-assembly) | KOfam / gene-ID prefix |
|---|---|
| `R1_sediment__concoct.bin.036` | `concoct.bin.036.fasta__R1_sediment` |
| `R1_sediment__concoct.bin.104` | `concoct.bin.104.fasta__R1_sediment` |
| `R1_sediment__maxbin2.Bin.094_sub` | `maxbin2.Bin.094_sub.fasta__R1_sediment` |
| `R1_sediment__metabat2.bin.142_sub` | `metabat2.bin.142_sub.fasta__R1_sediment` |
| `R1_sediment__metabat2.bin.186` | `metabat2.bin.186.fasta__R1_sediment` |
| `R1_sediment__metabat2.bin.206` | `metabat2.bin.206.fasta__R1_sediment` |
| `R1_sediment__metabat2.bin.258` | `metabat2.bin.258.fasta__R1_sediment` |
| `R2A_sediment__concoct.bin.189_sub` | `concoct.bin.189_sub.fasta__R2A_sediment` |
| `R2A_sediment__metabat2.bin.132` | `metabat2.bin.132.fasta__R2A_sediment` |
| `R2_sediment__concoct.bin.112` | `concoct.bin.112.fasta__R2_sediment` |
| `R2_sediment__metabat2.bin.112_sub` | `metabat2.bin.112_sub.fasta__R2_sediment` |
| `R2_sediment__metabat2.bin.128` | `metabat2.bin.128.fasta__R2_sediment` |
| `R2_sediment__metabat2.bin.154` | `metabat2.bin.154.fasta__R2_sediment` |
| `R2_sediment__metabat2.bin.162` | `metabat2.bin.162.fasta__R2_sediment` |
| `R2_sediment__metabat2.bin.405` | `metabat2.bin.405.fasta__R2_sediment` |
| `R2_sediment__metabat2.bin.414_sub` | `metabat2.bin.414_sub.fasta__R2_sediment` |
| `R2_sediment__metabat2.bin.435` | `metabat2.bin.435.fasta__R2_sediment` |
| `R2_sediment__concoct.bin.036_sub` | `concoct.bin.036_sub.fasta__R2_sediment` |
| `R1_sediment__concoct.bin.052_sub` | `concoct.bin.052_sub.fasta__R1_sediment` |
| `R2A_sediment__concoct.bin.151` | `concoct.bin.151.fasta__R2A_sediment` |
| `R2_sediment__concoct.bin.183_sub` | `concoct.bin.183_sub.fasta__R2_sediment` |
| `R1_sediment__metabat2.bin.121_sub` | `metabat2.bin.121_sub.fasta__R1_sediment` |
| `R2_sediment__metabat2.bin.141` | `metabat2.bin.141.fasta__R2_sediment` |
| `R1_sediment__metabat2.bin.356` | `metabat2.bin.356.fasta__R1_sediment` |
| `R2_sediment__metabat2.bin.362_sub` | `metabat2.bin.362_sub.fasta__R2_sediment` |
| `R2_sediment__metabat2.bin.434_sub` | `metabat2.bin.434_sub.fasta__R2_sediment` |

## Tools, versions, databases

| source | tool | database | gene calls |
|---|---|---|---|
| KOfam (existing; not re-run) | HMMER 3.4 `hmmsearch` (all 27,981 profiles), DRAM 1.5 best-KO rule (significant = score above the `ko_list` threshold for its `score_type`; best = lowest full-sequence E-value), Sec read-through rescue | KOfam release 2026-10-26 (`ko_list` md5 5dca1ddf…, see `annotation_runs/kofam_rollout/README.md`) | Prodigal 2.6.3 `-p meta -g 11`, contigs ≥1 kb |
| eggNOG | eggNOG-mapper **2.1.13**, DIAMOND 2.2.0 (`-m diamond`, emapper default `--sensmode sensitive`, `--dmnd_iterate no`), default tax scope / seed-ortholog thresholds | eggNOG DB **5.0.2** (`/scratch1/afreiburger/emilykin/processed/refdata/eggnog`; `eggnog.db`, `eggnog_proteins.dmnd`) | the same Prodigal proteins as KOfam (trailing `*` removed) |
| Bakta | Bakta **1.12.0** (Pyrodigal 3.7.1, tRNAscan-SE 2.0.12, Aragorn 1.2.41, Infernal 1.1.5, PilerCR 1.6.0, AMRFinderPlus 4.2.7, pyhmmer 0.12.1, DIAMOND 2.2.0, BLAST+ 2.17.0) | Bakta DB **v6.0 full**, 2025-02-24 (`/scratch1/afreiburger/emilykin/processed/refdata/bakta/db`; RefSeq 228, UniProtKB/Swiss-Prot 2025_01, Pfam 37.2, Rfam 15, COG 2024, AMRFinderPlus 2024-12-18.1) | Bakta's own (Pyrodigal, **default single-genome mode, not `--meta`**, as in `annotation_runs/bakta_test`) |

eggNOG-mapper and Bakta were both run from the EmilyKin env `544db00dd5c254ecfa7c6335967a046f_`, the same env and commands
as `annotation_runs/bakta_test`. The two changes from bakta_test are `--keep-contig-headers` and `--threads 4`.

**Why eggNOG runs on the KOfam Prodigal genes rather than on Bakta's CDS:** KOfam and eggNOG then annotate the same
genes, with the same IDs and coordinates, so arms B and C differ only in annotation, not in gene calling. It also let
eggNOG run alongside Bakta instead of after it.

## Commands (poplar, `scripts/driver.sh`, detached with `nohup setsid`, everything at `nice -n 10`)

```bash
ANNOT=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_
export PATH=$ANNOT/bin:$PATH
export TMPDIR=/scratch1/afreiburger/salterns_2026/annot57/tmp      # never /tmp (see "Temporary files")
# eggNOG: 2 pooled jobs (28 + 29 genomes, 89,510 + 90,895 proteins), 8 CPUs each
sed 's/\*$//' prodigal/<mag_id>.faa >> eggnog/gN/input.faa         # for each genome of group N
emapper.py --data_dir /scratch1/afreiburger/emilykin/processed/refdata/eggnog -m diamond --itype proteins \
    -i eggnog/gN/input.faa -o annot57_gN --output_dir eggnog/gN --temp_dir $TMPDIR/emapper_gN --cpu 8 \
    --dmnd_iterate no --dbmem --override
# Bakta: one genome per job, 4 threads; 8 workers while eggNOG ran, 12 after it ended (largest genomes first)
bakta --db /scratch1/afreiburger/emilykin/processed/refdata/bakta/db --threads 4 --skip-plot --keep-contig-headers \
    --tmp-dir $TMPDIR/bakta_<mag_id> --output bakta/<mag_id> --prefix <mag_id> fasta/<mag_id>.fa
```

Threads: at most 48 of mine at any time (eggNOG 2×8 + Bakta 8×4, then Bakta 12×4). The 1-min load on the
336 cores was 115–282 before and during the run (`logs/monitor.log`).

Local post-processing (light; `~/Documents/py_venv`):

```bash
cd ~/Documents/MicrobiomeNotebooks/salterns/reassembly_2026
bash modelling_2026/annotations/scripts/collect.sh                  # rsync + md5 check
~/Documents/py_venv/bin/python -I modelling_2026/annotations/scripts/make_sec_table.py $PWD $PWD/modelling_2026/annotations
~/Documents/py_venv/bin/python -I modelling_2026/annotations/scripts/summarize.py modelling_2026/annotations comparison_rollout/modelling_species.tsv
```

## Gene coordinates (for the merge onto RAST features)

All coordinates are 1-based and inclusive, with start < end on both strands, on the **original contig IDs** (`NODE_…`
for rollout, `k127_…` for co-assembly):

| source | table | coordinate columns | annotation columns |
|---|---|---|---|
| KOfam | `kofam/<mag_id>.tsv.gz` | `contig, start, end, strand, partial` (+ `contig_len, aa_len`) | `kofam_id` (best KO), `kofam_score/threshold/score_type/full_evalue/hmm_cov`, `kofam_definition` (KO name, EC in `[EC:…]`), `kofam_all_passing` (`;`-separated) |
| KOfam Sec rescue | `kofam/sec_readthrough_57.tsv` | `contig, strand, orf_start, orf_end` (+ `tga_lo/tga_hi`) | `accepted`, `ko`, `status`, `pieces` (the Prodigal genes the joined ORF replaces), scores, SECIS |
| eggNOG | `eggnog/<mag_id>/<mag_id>.eggnog_genes.tsv` | `contig, start, end, strand, partial` | `ko`, `KEGG_ko`, `EC`, `BiGG_Reaction`, `KEGG_Reaction`, `Preferred_name`, `Description`, `COG_category`, `GOs`, `PFAMs`, … (lists comma-separated, as emapper writes them) |
| Bakta | `bakta/<mag_id>/<mag_id>.bakta_genes.tsv` (or the raw `.tsv`/`.gff3`/`.json`) | `contig, start, end, strand` (`truncated` = contig-edge `5-prime`/`3-prime`) | `ko`, `ec`, `product`, `gene`, `hypothetical`, `evidence`, `go`, `cog`, `uniref`, `pfam` (lists `;`-separated) |

**The KOfam per-gene tables are pre-rescue.** For example, in `R1_C_D1__maxbin2.Bin.003_sub` the TGA-truncated GrdH piece
`NODE_3115…_4` is still `K10672` (grdB). Apply `kofam/sec_readthrough_57.tsv` on top: for rows with
`accepted == True`, the joined ORF (`orf_start..orf_end`) carries `ko` and replaces the genes in `pieces`. These 57
genomes have 73 junctions in 24 genomes: 53 accepted (grdH 17, of which 14 are new calls; grdB 13; grdA 23) and 20
rejected. The ORF coordinates were checked: each ORF contains all its pieces, and where the upstream piece was called,
the ORF shares its start.

**Bakta vs Prodigal genes** (`gene_caller_concordance.tsv`): Bakta calls fewer CDS (169,406 vs 180,405; per genome 0.80–1.00×,
median 0.95×). The biggest differences are in Halanaerobium, M55B127 and other fragmented Bacillota MAGs. Median
per genome: 97.9 % of Bakta CDS share a stop codon with a Prodigal gene, 92.6 % of Prodigal genes share one with a Bakta
CDS, and 87 % of Bakta CDS have identical coordinates. This also confirms that the two coordinate systems agree:
Bakta's `start`/`end` are left/right positions on the same contig IDs as the KOfam/eggNOG tables. Bakta also joins Sec
genes: 28 CDS with a selenocysteine in 18 genomes.

## Coverage

Fractions are of each source's own CDS: Prodigal's 180,405 genes for KOfam and eggNOG, Bakta's 169,406 CDS for Bakta.
"Product" means: KOfam, a KO definition (so the same as KO); eggNOG, a non-empty `Description` or `Preferred_name`;
Bakta, any product other than "hypothetical protein".

| | CDS | with KO | with EC | with BiGG reaction | with product name |
|---|---:|---:|---:|---:|---:|
| KOfam (pre-rescue) | 180,405 | 86,864 (48.1 %) | 46,125 (25.6 %) | – | 86,864 (48.1 %) |
| eggNOG | 180,405 (167,662 = 92.9 % with a seed ortholog) | 100,806 (55.9 %) | 51,936 (28.8 %) | 1,565 (0.9 %) | 156,783 (86.9 %) |
| KOfam ∪ eggNOG (same genes) | 180,405 | 111,004 (61.5 %) | | | |
| Bakta | 169,406 (+26 sORFs, 267 pseudogenes) | 11,685 (6.9 %) | 10,059 (5.9 %) | – | 116,719 (68.9 %) |

Per-genome ranges (min / median / max % of CDS): KOfam KO 29.4 / 49.4 / 64.8; eggNOG KO 40.9 / 57.7 / 68.5; eggNOG EC
21.4 / 29.5 / 35.8; eggNOG BiGG 0.0 / 0.2 / 10.8; eggNOG product 70.2 / 88.2 / 95.0; Bakta KO 0.0 / 1.5 / 26.4; Bakta EC
0.0 / 1.1 / 23.7; Bakta product 7.4 / 73.4 / 93.2.

What this means for the arms:
- **eggNOG BiGG reactions are rare.** 1,565 genes in total, median 9 per genome, none in 9 genomes. The most are in
  WSZY01 `R2A_sediment__concoct.bin.151` (233), Methanolobus (107), DAOWES01 (102) and Desulfotignum (79). eggNOG
  assigns BiGG IDs only through orthologs in BiGG-modelled organisms, so for arm C the KO and EC layers carry almost all
  of eggNOG's added evidence.
- **Bakta's KO/EC layer is sparse on these novel MAGs**, as bakta_test found. Bakta gives a KO or EC only for an exact
  UPS/IPS match or a ≥90 % identity UniRef90 (PSC) hit. Across all CDS: UPS/IPS 3.2 %, PSC 18.7 %, PSCC (UniRef50, product
  name only) 48.5 %, no hit 29.6 %. Bakta KO is ≥20 % only in the 10 genomes with close references (5 Roseovarius,
  Marivita, Roseibacterium, Saliniramus, QPJE01, one Halanaerobium). It is <5 % in 37 genomes and <1 % in 25. Arm D's
  extra evidence is therefore mostly Bakta product names, plus EC for the near-reference genomes.

### Per genome (`coverage_per_genome.tsv` has all columns, including eggNOG `KEGG_Reaction`, Bakta evidence classes, tRNA/rRNA)

"KOfam Sec joins" counts the accepted rescue junctions (`kofam/sec_readthrough_57.tsv`). They are not included in "KOfam KO".

| MAG | set | genus | Prodigal CDS | KOfam KO | KOfam EC | KOfam Sec joins | eggNOG KO | eggNOG EC | eggNOG BiGG | eggNOG product | Bakta CDS | Bakta KO | Bakta EC | Bakta product | Bakta wall (min) |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `R1_sediment__concoct.bin.036` | coassembly | CAZRZN01 | 4,547 | 2,180 | 1,132 | 0 | 2,582 | 1,312 | 1 | 4,065 | 4,229 | 211 | 190 | 3,439 | 54.1 |
| `R1_sediment__concoct.bin.052_sub` | coassembly | CAYRVO01 | 1,439 | 726 | 373 | 1 | 939 | 480 | 29 | 1,308 | 1,186 | 4 | 2 | 748 | 15.0 |
| `R1_sediment__concoct.bin.104` | coassembly | Roseovarius | 3,194 | 1,577 | 863 | 0 | 1,953 | 1,041 | 1 | 2,907 | 2,912 | 113 | 101 | 2,358 | 22.1 |
| `R1_sediment__maxbin2.Bin.094_sub` | coassembly | JAWLGC01 | 2,982 | 1,607 | 789 | 0 | 1,807 | 851 | 65 | 2,661 | 2,836 | 8 | 4 | 2,013 | 20.5 |
| `R1_sediment__metabat2.bin.121_sub` | coassembly | Roseovarius | 3,288 | 1,824 | 1,012 | 0 | 2,019 | 1,083 | 1 | 2,988 | 3,224 | 852 | 765 | 2,996 | 17.5 |
| `R1_sediment__metabat2.bin.142_sub` | coassembly | RI-231 | 3,594 | 1,055 | 591 | 0 | 1,471 | 799 | 19 | 2,950 | 3,174 | 33 | 33 | 1,724 | 21.2 |
| `R1_sediment__metabat2.bin.186` | coassembly | DAOWES01 | 2,210 | 991 | 565 | 4 | 1,279 | 695 | 102 | 1,869 | 2,083 | 5 | 2 | 1,310 | 18.4 |
| `R1_sediment__metabat2.bin.206` | coassembly | Halanaerobium | 2,234 | 1,314 | 684 | 2 | 1,486 | 728 | 41 | 2,067 | 2,094 | 49 | 22 | 1,718 | 20.0 |
| `R1_sediment__metabat2.bin.258` | coassembly | JACKII01 | 4,016 | 2,049 | 1,181 | 0 | 2,096 | 1,141 | 5 | 3,643 | 3,955 | 14 | 11 | 3,084 | 52.1 |
| `R1_sediment__metabat2.bin.356` | coassembly | QPJE01 | 1,813 | 1,174 | 603 | 1 | 1,241 | 594 | 32 | 1,723 | 1,719 | 368 | 271 | 1,594 | 13.7 |
| `R2A_sediment__concoct.bin.151` | coassembly | WSZY01 | 2,165 | 1,136 | 625 | 0 | 1,398 | 712 | 233 | 1,876 | 2,042 | 25 | 24 | 1,389 | 16.5 |
| `R2A_sediment__concoct.bin.189_sub` | coassembly | nan | 3,077 | 1,210 | 620 | 0 | 1,577 | 868 | 9 | 2,685 | 2,763 | 30 | 26 | 1,755 | 20.7 |
| `R2A_sediment__metabat2.bin.132` | coassembly | CAMCWB01 | 4,683 | 2,022 | 1,047 | 0 | 2,351 | 1,238 | 11 | 4,051 | 4,512 | 22 | 14 | 3,117 | 53.9 |
| `R2_sediment__concoct.bin.036_sub` | coassembly | Roseovarius | 2,247 | 981 | 520 | 0 | 1,355 | 720 | 2 | 1,996 | 2,009 | 111 | 106 | 1,475 | 17.9 |
| `R2_sediment__concoct.bin.112` | coassembly | CAJQQK01 | 3,961 | 1,914 | 1,118 | 0 | 2,078 | 1,173 | 9 | 3,445 | 3,968 | 21 | 12 | 2,821 | 51.8 |
| `R2_sediment__concoct.bin.183_sub` | coassembly | JBKDZR01 | 3,836 | 1,947 | 1,026 | 0 | 2,246 | 1,132 | 2 | 3,416 | 3,639 | 56 | 38 | 2,930 | 20.1 |
| `R2_sediment__metabat2.bin.112_sub` | coassembly | Lentibacter | 3,049 | 1,496 | 784 | 0 | 1,801 | 947 | 0 | 2,731 | 2,780 | 513 | 416 | 2,360 | 18.0 |
| `R2_sediment__metabat2.bin.128` | coassembly | DATHSN01 | 3,611 | 1,384 | 805 | 0 | 1,778 | 1,037 | 26 | 3,048 | 3,317 | 17 | 4 | 2,100 | 20.9 |
| `R2_sediment__metabat2.bin.141` | coassembly | PWPR01 | 2,206 | 1,079 | 590 | 1 | 1,221 | 627 | 48 | 1,925 | 2,065 | 1 | 1 | 1,176 | 17.7 |
| `R2_sediment__metabat2.bin.154` | coassembly | JAAAOG01 | 3,256 | 1,405 | 721 | 0 | 1,704 | 843 | 11 | 2,771 | 2,904 | 7 | 8 | 1,669 | 21.5 |
| `R2_sediment__metabat2.bin.162` | coassembly | M55B127 | 2,422 | 1,317 | 677 | 5 | 1,495 | 747 | 45 | 2,191 | 2,263 | 3 | 2 | 1,454 | 19.2 |
| `R2_sediment__metabat2.bin.362_sub` | coassembly | M55B127 | 3,595 | 1,444 | 721 | 3 | 1,868 | 901 | 48 | 2,919 | 3,018 | 4 | 3 | 1,552 | 24.4 |
| `R2_sediment__metabat2.bin.405` | coassembly | Thermohalobaculum | 3,859 | 1,803 | 936 | 0 | 2,160 | 1,092 | 5 | 3,330 | 3,557 | 23 | 10 | 2,552 | 20.7 |
| `R2_sediment__metabat2.bin.414_sub` | coassembly | DAINOF01 | 4,479 | 1,974 | 1,026 | 0 | 2,404 | 1,266 | 3 | 3,847 | 4,080 | 67 | 60 | 3,047 | 55.1 |
| `R2_sediment__metabat2.bin.434_sub` | coassembly | UBA8670 | 3,298 | 1,467 | 788 | 2 | 1,787 | 926 | 55 | 2,761 | 2,893 | 101 | 68 | 2,054 | 24.4 |
| `R2_sediment__metabat2.bin.435` | coassembly | Halonatronomonas | 2,484 | 1,291 | 660 | 4 | 1,575 | 767 | 42 | 2,225 | 2,226 | 228 | 157 | 1,821 | 18.0 |
| `R1_A_D1__concoct.bin.007` | rollout | Methanolobus | 2,375 | 1,130 | 635 | 0 | 1,345 | 692 | 107 | 2,009 | 2,299 | 0 | 0 | 379 | 19.7 |
| `R1_A_D2__maxbin2.Bin.123_sub` | rollout | BM714 | 2,522 | 1,261 | 680 | 5 | 1,455 | 758 | 40 | 2,194 | 2,253 | 245 | 175 | 1,699 | 18.5 |
| `R1_A_D2__metabat2.bin.020` | rollout | DATFJG01 | 3,219 | 1,534 | 825 | 3 | 1,694 | 861 | 57 | 2,766 | 3,064 | 2 | 1 | 1,750 | 22.4 |
| `R1_A_D2__metabat2.bin.053` | rollout | Roseovarius | 3,598 | 1,885 | 1,012 | 0 | 2,121 | 1,102 | 2 | 3,244 | 3,456 | 700 | 652 | 3,123 | 19.5 |
| `R1_B_D1__maxbin2.Bin.073` | rollout | Natrinema | 3,657 | 1,342 | 728 | 0 | 1,660 | 891 | 0 | 2,878 | 3,400 | 1 | 1 | 303 | 24.0 |
| `R1_B_D1__metabat2.bin.017_sub` | rollout | Roseovarius | 3,586 | 1,980 | 1,071 | 0 | 2,100 | 1,107 | 0 | 3,210 | 3,478 | 832 | 732 | 3,232 | 19.1 |
| `R1_B_D1__metabat2.bin.040` | rollout | Roseovarius | 3,349 | 1,770 | 991 | 0 | 1,993 | 1,080 | 0 | 2,995 | 3,259 | 705 | 616 | 2,898 | 18.1 |
| `R1_B_D2__concoct.bin.112_sub` | rollout | Halanaerobium | 2,169 | 1,069 | 549 | 3 | 1,424 | 725 | 40 | 2,005 | 1,735 | 310 | 256 | 1,435 | 16.5 |
| `R1_B_D2__metabat2.bin.065` | rollout | Roseovarius | 3,228 | 1,768 | 944 | 0 | 1,946 | 996 | 1 | 2,929 | 3,116 | 510 | 483 | 2,868 | 18.5 |
| `R1_B_H2O__concoct.bin.053` | rollout | Roseovarius | 3,858 | 1,877 | 995 | 0 | 2,197 | 1,159 | 1 | 3,375 | 3,634 | 618 | 556 | 2,976 | 20.8 |
| `R1_C_D1__concoct.bin.031` | rollout | UBA8670 | 2,406 | 1,274 | 664 | 4 | 1,458 | 753 | 51 | 2,176 | 2,154 | 100 | 62 | 1,689 | 18.4 |
| `R1_C_D1__maxbin2.Bin.003_sub` | rollout | M55B127 | 2,279 | 1,231 | 631 | 2 | 1,336 | 644 | 39 | 2,010 | 2,162 | 4 | 2 | 1,407 | 20.2 |
| `R1_C_D2__metabat2.bin.048` | rollout | Marivita | 3,512 | 1,993 | 1,013 | 0 | 2,133 | 1,040 | 1 | 3,152 | 3,493 | 916 | 785 | 3,254 | 16.2 |
| `R2_A_D1__concoct.bin.120` | rollout | Saliniramus | 3,594 | 1,919 | 1,014 | 0 | 2,068 | 1,056 | 2 | 3,224 | 3,488 | 833 | 742 | 3,061 | 18.0 |
| `R2_A_D1__maxbin2.Bin.018_sub` | rollout | M55B127 | 3,999 | 1,585 | 774 | 4 | 1,766 | 856 | 47 | 2,806 | 3,805 | 2 | 1 | 1,791 | 52.8 |
| `R2_A_D2__metabat2.bin.055` | rollout | Haloplanus | 3,257 | 1,309 | 681 | 0 | 1,448 | 755 | 2 | 2,555 | 3,153 | 0 | 1 | 357 | 23.8 |
| `R2_A_H2O__metabat2.bin.026` | rollout | Halorubrum | 2,336 | 1,047 | 554 | 0 | 1,226 | 624 | 0 | 1,869 | 2,281 | 0 | 1 | 237 | 18.9 |
| `R2_B_D2__maxbin2.Bin.065` | rollout | Roseovarius | 3,012 | 1,705 | 947 | 0 | 1,890 | 1,025 | 1 | 2,760 | 2,930 | 631 | 571 | 2,718 | 19.4 |
| `R2_B_D2__metabat2.bin.006_sub` | rollout | Roseovarius | 3,909 | 1,977 | 1,083 | 0 | 2,246 | 1,210 | 0 | 3,482 | 3,732 | 212 | 197 | 3,259 | 19.7 |
| `R2_B_D2__metabat2.bin.056` | rollout | M55B127 | 3,000 | 1,569 | 769 | 6 | 1,734 | 822 | 44 | 2,662 | 2,888 | 5 | 1 | 1,843 | 21.5 |
| `R2_B_D2__metabat2.bin.063` | rollout | M55B127 | 1,879 | 902 | 444 | 1 | 1,151 | 598 | 32 | 1,651 | 1,557 | 2 | 2 | 954 | 16.3 |
| `R2_B_H2O__metabat2.bin.051` | rollout | Haloquadratum | 2,927 | 1,102 | 586 | 0 | 1,346 | 677 | 0 | 2,111 | 2,770 | 1 | 1 | 215 | 20.0 |
| `R2_C_D1__metabat2.bin.057_sub` | rollout | Desulfotignum | 4,469 | 1,981 | 1,014 | 0 | 2,266 | 1,094 | 77 | 3,880 | 4,318 | 128 | 106 | 3,271 | 55.9 |
| `R2_C_D2__concoct.bin.063` | rollout | Roseibacterium | 3,941 | 2,087 | 1,135 | 0 | 2,268 | 1,169 | 1 | 3,508 | 3,900 | 788 | 740 | 3,540 | 18.4 |
| `R2_C_D2__concoct.bin.072` | rollout | JAVJPJ01 | 3,733 | 1,951 | 1,062 | 0 | 2,169 | 1,152 | 0 | 3,289 | 3,678 | 189 | 155 | 3,216 | 21.2 |
| `R2_C_D2__metabat2.bin.065` | rollout | Halanaerobium | 2,423 | 1,255 | 673 | 2 | 1,641 | 867 | 48 | 2,274 | 2,003 | 484 | 373 | 1,751 | 18.5 |
| `R2_C_H2O__metabat2.bin.005` | rollout | Natronomonas | 3,125 | 1,279 | 694 | 0 | 1,458 | 771 | 1 | 2,511 | 3,083 | 0 | 0 | 309 | 23.6 |
| `R2_C_H2O__metabat2.bin.008` | rollout | Halonotius | 1,970 | 747 | 378 | 0 | 885 | 441 | 0 | 1,514 | 1,906 | 0 | 0 | 141 | 15.7 |
| `R2_restored_C_black__maxbin2.Bin.042` | rollout | Desulfotignum | 3,867 | 1,943 | 1,007 | 0 | 2,092 | 1,039 | 79 | 3,477 | 3,800 | 135 | 102 | 3,173 | 20.4 |
| `R2_restored_C_black__metabat2.bin.011_sub` | rollout | Roseovarius | 3,424 | 1,764 | 925 | 0 | 2,157 | 1,107 | 1 | 3,134 | 3,123 | 418 | 374 | 2,710 | 22.8 |
| `R2_restored_D_shore__concoct.bin.011_sub` | rollout | JASFMS01 | 4,237 | 2,261 | 1,190 | 0 | 2,462 | 1,145 | 46 | 3,735 | 4,040 | 28 | 21 | 2,904 | 53.3 |
| **total** | | | **180,405** | **86,864** | **46,125** | **53** | **100,806** | **51,936** | **1,565** | **156,783** | **169,406** | **11,685** | **10,059** | **116,719** | |

## Run times (`runtimes.tsv`)

| job | wall | CPU | peak RSS |
|---|---|---|---|
| eggNOG g1 (28 genomes, 89,510 proteins, `--cpu 8`) | 66.7 min | 384 CPU-min | 41.1 GB (`--dbmem`) |
| eggNOG g2 (29 genomes, 90,895 proteins, `--cpu 8`) | 67.3 min | 394 CPU-min | 41.2 GB |
| Bakta, 57 jobs (`--threads 4`) | median 20.0 min (13.7–55.9) | 53.8 CPU-h total (38–74 CPU-min each) | 7.6–12.5 GB |
| whole run | 2 h 26 min (18:51:37 → 21:18:02 UTC) | | |

Bakta is I/O-bound, not CPU-bound. The first 8 jobs took 52–56 min each, at 116 % CPU: 42 of those minutes went to the
PSC DIAMOND search reading the 46 GB `psc.dmnd` from NFS (`/scratch1`) at about 20 MB/s, while the jobs waited in D
state. Once the file was in poplar's page cache, the PSC step read at 200–900 MB/s and later jobs took 14–23 min. So
Bakta on more genomes costs about 20 min per genome per worker, plus one cold read of `psc.dmnd`.

## Temporary files

- `TMPDIR` was exported to `annot57/tmp` (on `/scratch1`). Bakta used `--tmp-dir`, and emapper used `--temp_dir`, which
  it passes to DIAMOND `--tmpdir`. Bakta passes the environment to its tools: AMRFinderPlus wrote `amrfinder.*` and
  tRNAscan-SE wrote `tscan*`, both in `annot57/tmp`. tRNAscan-SE's own config says `temp_dir: /tmp`, but its code uses
  `$TMPDIR` when set, and it did.
- `/tmp` on poplar was 9.4–9.5 GB used of 32 GB (30 %) before, during (every 10 min, `logs/monitor.log`) and after the run.
  Files owned by `afreiburger` in `/tmp`: 2,817 before and 2,817 after. These are older and not from this run; they were
  left untouched.
- AMRFinderPlus does not clean up its 57 `amrfinder.*` temp dirs, so they and the tRNAscan-SE leftovers were deleted from
  `annot57/tmp` after the run. `annot57/tmp` is empty.

## Failures and problems

- **None.** Bakta 57/57 OK (all `.json` written, exit 0). eggNOG 2/2 exit 0. KOfam 57/57 tables present, with row counts
  equal to the Prodigal protein counts on poplar.
- Not a failure, but a caveat for the comparison: **Bakta's default single-genome Pyrodigal mode** was used (as in
  bakta_test), while KOfam and eggNOG use Prodigal `-p meta`. So part of the CDS-count difference (above) comes from the
  training mode.
- The **KOfam per-gene tables are pre-Sec-rescue**. Use `kofam/sec_readthrough_57.tsv` with them (see above). eggNOG ran
  on the unjoined Prodigal pieces. Bakta's own Sec joins are in its CDS (`exception` column).
- emapper's `.hits` files have no header line in emapper's output. A header was added to the per-genome split copies only.
  The per-genome `.emapper.annotations` and `.seed_orthologs` copies keep the column header but drop emapper's `##`
  metadata lines (date, command, query count, rate). The `_pooled/` files are unchanged and keep everything.
- Rule deviations during the run. None of them affected any result:
  - One invocation of the system `python3`: a heredoc edit that failed and ran an empty script. Nothing executed. The
    edit was redone with the editor.
  - One read-only `pgrep -f` on poplar with a pattern (`psc.dmnd`) that also matched my own ssh shell, used to read
    DIAMOND `/proc/<pid>/io` counters. No signals were sent. Later checks used `pgrep -x diamond`.
  - Small scratch files written outside `annot57/` on poplar during the two I/O checks (`/dev/shm/x`, `~/.annot57_r1`,
    `~/annot57_rchar.tmp`). Each was removed in the same command, and none exists now.
