# KOfam (DRAM's KO layer) on all 1,324 re-assembly MAGs — poplar run

*Track `kofam`, started 2026-11-11 (US Central 10:31) on poplar. Follows `../../annotation_choice/README.md` §6.1: KO calls for the guilds come from KOfam with KEGG's per-KO thresholds, on contigs ≥1,000 bp, for all 1,324 MAGs (not only the 558 species representatives).*

## Status (2026-11-12): DONE

**ALL_DONE 2026-11-12 10:11Z (40/40 chunks, 1,324 MAGs). Results copied to `results/`; Sec read-through rescue, per-MAG guild calls (`results/guilds_1324.tsv`) and the comparison with the 310 originals are in [`RESULTS.md`](RESULTS.md).** The launch-time status below is kept for the record.

## Status at launch (2026-11-11 11:55 CST)

**ssh / poplar connection.** The ControlMaster connection (`ssh -o BatchMode=yes poplar`) worked from 10:31 CST and stayed up for this whole session. Every poplar job runs under `nohup setsid`, so if the connection drops later the jobs keep running. Only checking on them is blocked: log in again in tmux (`ssh poplar`, Duo), then re-run the status commands below.

**Running: 40 hmmsearch jobs on poplar, started 2026-11-11 11:36:39 CST (17:36:39Z).** Expected finish **2026-11-12, about 05:00–11:00 CST**. Nothing is called for the full profile set until a chunk finishes. Each chunk postprocesses itself when its search ends.

- **Driver:** `run_all.sh`, PID **447804** (nohup + setsid, `poplar:/scratch1/afreiburger/salterns_2026/kofam`). When the last chunk ends it runs `postprocess.py combine` and touches `ALL_DONE`.
- **Per-chunk PIDs:** `run_chunk.sh` PIDs (447810…) are in `logs/chunk_pids.tsv`. The hmmsearch PIDs (447840…448191) are in `poplar_copies/hmmsearch_pids.tsv`.
- **Threads:** 40 single-CPU jobs at nice 10. Each `--cpu 1` hmmsearch has one worker thread plus a mostly idle reader thread, and each job got ≈85 % of a core.
- **Host when launched:** load average 520–640 on 336 cores, 97 % CPU busy, 2.2 TB RAM available. Each job uses ≈30–50 MB RSS.

### Timing and extrapolation

| measurement | value |
|---|---|
| 300 random profiles (seed 1; 0.98 % of total profile length, ΣLENG 161,780 of 16,466,034) × chunk 01 (74,966 searched proteins, 20.2 M residues) | **512 CPU-s**, 8 min 36 s wall, 100 % CPU |
| ⇒ one chunk × all 27,981 profiles, scaled by ΣLENG | **≈14.5 CPU-h** (13.3 h if scaled by profile count) |
| the 40 real chunks, first 16 min (profiles in K-number order, K00001…) | 1.21 % of ΣLENG done on average (1.19–1.40 %); chunk 01 used 14 CPU-min ⇒ **≈19.5 CPU-h/chunk** at this rate |
| **forecast per chunk** | **14.5–19.5 CPU-h ⇒ ≈17–23 h wall at ≈85 % of a core**. The low K numbers are common oxidoreductase families with many hits, so the 16-minute rate is probably pessimistic. `scripts/progress.sh` prints a live ETA. |
| total | ≈580–780 CPU-h. Without de-duplication it would be ≈800–1,080 CPU-h; the README's ±2× pre-estimate was ≈630 CPU-h |
| layout test: 8 of the same profiles × all 2.83 M searched proteins (the "split profiles" layout, same profile×residue work) | **858 CPU-s vs 512**, so the profile-split layout is 1.7× slower. The run uses **protein chunks × all profiles**. |
| postprocessing | 7 s and 145 MB for 368 profiles × chunk 01. A full chunk should take minutes and a few GB, then gzip. |
| disk | ≈1.9 GB of raw hmmsearch tables after 1.2 %, so ≈150 GB uncompressed for all 40 chunks before the per-chunk gzip (NFS has 13 TB free) |

### Inputs and set-up, done

- **MAGs:** 1,324 FASTAs on poplar, and **md5 verified for all 1,324**. Together they hold 752,569 contigs, all ≥1,000 bp, and no contig is shared between two MAGs of the same unit.
- **Genes:** 3,911,146 proteins (1,115,778,826 aa) from Prodigal on all 1,324 MAGs. 8 MAGs failed first on an NFS temp-dir cleanup error (`ENOTEMPTY`) and were re-run after the temp dir moved to local `/tmp`. The pooled gene table is `genes_all.tsv.gz` on poplar (gene → contig → contig_len → MAG → binner, …).
- **Prodigal build difference:** poplar's bioconda build gives 2 extra ~50-aa ORFs on the *Methanolobus* control (2,279 vs 2,277 genes with the local `gtdbtk_venv` Prodigal 2.6.3). All other coordinates are identical.
- **KOfam DB:** 27,981 profiles. All have a ko_list row. 523 ko_list KOs have no profile (e.g. K01977–K01986, K14218–K14239), so they can never be called. 1,870 KOs have threshold `-`. The minimum numeric threshold is 14.70. All 68 markers have a profile. Details are in `poplar_copies/db_check.txt`.
- **De-duplication:** 28.4 % of proteins are exact copies of a protein in another MAG, the same strain binned from a per-sample assembly and a site co-assembly.
  - Chunks are therefore built from **whole species clusters** (558 skani clusters, LPT-balanced). The control clusters are in chunk 01.
  - Each chunk searches **one representative per distinct sequence**: 2,826,115 proteins, 808.5 M aa, **27.5 % less work**. Hits are copied to identical proteins (`chunks/chunk_XX.members.tsv`).
  - This is exact: hmmsearch scores do not depend on the database, and E-values do only through `-Z`, which is fixed at 3,911,146.
  - Every chunk still holds whole MAGs (24–39 MAGs, 12–15 clusters, 20.13–20.29 M searched aa per chunk; `poplar_copies/chunk_stats.txt`, `chunk_assignment.tsv`).

### Postprocessing test on chunk 01 (marker profiles only, `test/markers/`)

This is the same chunk-01 protein file as the production run, searched with the 68 marker profiles only. Like the local control test, it takes the best KO among markers.

- `scripts/compare_controls.py` compares it with `annotation_choice/kofam_control_test/marker_calls_by_mag.tsv`, which used KOfam 2026-05-28, local Prodigal and min_contig 1,000. The comparison covers the 8 control bins: 2 MaxBin2 *Roseovarius*, the best *Methanolobus* and all 5 *Albimonas* bins.
- **544 MAG × marker cells, 3 differ.** All 3 are *phnL* 0 → 1 in Albimonas `concoct.bin.061`, `concoct.bin.018` and `concoct.bin.058_sub`, with scores 223–226.
- The cause is the release: the *phnL* threshold fell from 228.40 (2026-05) to 216.47 (2026-10, domain). There are no other differences, which confirms gene calls, rule and expansion of identical proteins (`test/control_comparison.tsv`).
- **Guild calls on the controls:**
  - both MaxBin2 *Roseovarius*: `mepn_consumer` with all 7 of *phnGHIJKLM*, the locus split over 3 and 4 contigs; `mttB_family_without_mcrA`. R1_A_D2 is also a one-sided TMA candidate (*cntA/yeaW* without *cntB*), as found for the 2020 *Roseovarius*;
  - *Methanolobus* `concoct.bin.101`: `methanogen` and `methylamine_user`;
  - all 5 *Albimonas* bins: **no *phnJ*** call, with 2–5 of *phnGHIJKLM*, as in `annotation_choice/README.md`. The best *phnJ* hits all fall below the 194.6 full-score threshold:
    - `concoct.bin.157_sub`: the 371-aa fragment scores 133.1, with 36 % HMM coverage;
    - `metabat2.bin.062`: 65.7;
    - `concoct.bin.018`: 26.4;
    - `concoct.bin.058_sub` and `concoct.bin.061` have no *phnJ* hit at all (score ≥14.7).
- **Short MaxBin2 contigs:** `maxbin2.Bin.031_sub.fasta__R1_B_D1` gets its only *phnJ* (score 523) from a 1,792 bp MaxBin2 contig, so its `mepn_consumer` call would be lost at DRAM's 2,500 bp default (`guilds_lost_without_short_maxbin2_contigs`).
- The chunk 01 totals at this stage are 3 methanogens, 8 MePn consumers and 0 betaine reducers. They are not final.
- **Full code path** (`test/subset/`, 300 random + 68 marker profiles × chunk 01, tblout + domtblout):
  - 546 hits to threshold-`-` KOs;
  - 1,740 of 34,380 reported (gene, KO) pairs have no domain ≥14.7, and **none of them pass**;
  - marker counts are unchanged by adding non-marker profiles; only `n_genes_with_ko` grows.

## Independent review (2026-11-11, 11:57–12:12 CST)

**Verdict: the run is correct. Nothing was killed or relaunched.** The only correction is a README warning: don't edit the two bash scripts while the jobs run (see "Do not edit" below).

**ssh / your tmux message.** The non-interactive ControlMaster ssh (`ssh -o BatchMode=yes poplar`) worked for every command of this review. The driver is detached from any login: PID 447804 has PPID 1 and is the leader of its own session (SID 447804), and all 40 chunk jobs are in that session. So the "poplar connection ended" in tmux does not affect them. Log in again only to check on them.

**Checked and found correct:**

| what | how | result |
|---|---|---|
| jobs running | `ps` on poplar at 12:00 | driver 447804 alive, 40 `run_chunk.sh` + 40 `hmmsearch` (447840…448191), all nice 10, 86–91 % CPU each, 31–111 MB RSS, no errors in `logs/chunk_*.log`, `run_all.nohup` empty |
| thread / CPU budget | `ps -L` | each `--cpu 1` hmmsearch = 1 worker + 1 master thread. The master re-reads the chunk FASTA for every profile and uses ≈8 % of a core. The track uses ≈36 cores in total, within the 48 limit. |
| scripts | md5 of `scripts/*` and `meta/*` on poplar vs here | identical |
| hmmsearch command | `/proc` cmdline, header of `hmm/chunk_01.out` | HMMER 3.4, `--cpu 1 --noali -T 14.7 --domT 14.7 -Z 3911146`, query `db/kofam_all.hmm`, target `chunks/chunk_XX.faa` |
| KOfam release | md5 on poplar and here | `kofam_db/ko_list.gz` = `697af5b3…` = local `annotation_choice/sources/kegg/kofam_ko_list_2026-10-26.gz`; `profiles.tar.gz` md5 `bc68ad2a…`, gzip CRC OK. The 27,983 files in `profiles/` = 27,981 `K*.hmm` + `eukaryote.hal` + `prokaryote.hal`. |
| thresholds | `db/hmm_info.tsv` vs `kofam_db/ko_list` and `annotation_choice/README.md` | 27,981 rows (17,950 full, 8,684 domain, 1,347 `-`). Spot values: *phnJ* 194.63 full, *phnL* 216.47 domain, *mttB* 218.93, *grdI* 505.90, *grdH* 624.03, *mcrA* 741.47, *cntA* 313.00, *cutC* 988.90; *mpnS*, *phnY\**, *phnZ1* `-`. All match. |
| MAG set | `comparison/mags_old_vs_new.tsv` rows `new|<name>` vs `mag_manifest.tsv` | the same 1,324 names; sample (unit) and species cluster match for all; every local path exists; md5 matches for 25 random local files and 40 random poplar copies |
| genes | the 1,324 `genes/*.done` records | 752,569 contigs in = used (none < 1,000 bp), 3,911,146 genes = `-Z` |
| chunks | `chunk_assignment.tsv`, `members.tsv`, `grep -c '>'` | 1,324 MAGs, each in one chunk; no species cluster split across chunks; for chunks 01, 17 and 40, searched representatives + copied members = the chunk's assigned protein count (100,239 / 104,634 / 85,555) |
| KO rule vs DRAM 1.5 | `test/review/review_check_partial.py` on the real `chunk_01` and `chunk_24` tables so far (667 and 749 complete KOs; 81,077 and 67,100 gene–KO pairs) | tblout "best 1 domain" score = max domtblout domain score for every pair. The significant pairs from postprocess.py's tblout rule equal those from DRAM's per-row domtblout `sig_scores` rule: 3,077 = 3,077 and 2,782 = 2,782, 0 differences. |
| postprocess.py on real full-DB output | `test/review/review_pp_partial.sh`: chunk 01's tables cut at the last complete KO (K00824), in a separate dir `kofam/test/review_partial_pp/` | ran in 6 s with 150 MB peak RSS. 4,530 of 100,239 genes get a KO so far, `passing_pairs_without_domtblout_row` 0. *mcrA* (K00399) is called at 982–983 in the three CONCOCT *Methanolobus*-cluster bins, including the control `concoct.bin.101.fasta__R1_sediment`; a 28-aa fragment (33.3) correctly fails. Outputs are in `test/review/partial_pp_chunk01/`. A full chunk has about 50× more rows, which extrapolates to about 5 min and 1–3 GB per chunk. |
| control comparison | `test/control_comparison.tsv` | 544 cells, 541 equal, 3 differ (*phnL* in 3 *Albimonas* bins), as reported |

**Caveats (documented here; none needs a relaunch):**
- **Score cut on the full-sequence score.** `-T 14.7` reports a gene–KO pair only when its *full* score is ≥14.7. On the real data a domain score is at most 1.0 bit above its full score (8 of 81,077 pairs). So a hit could be missed only for a domain-type KO with a threshold below about 15.7, and only inside that 1-bit window. Three KOs qualify: K02882 (L18Ae, 14.70), K01611 (SAM decarboxylase, 15.37) and K07484 (transposase, 15.43). None is a marker. DRAM's per-MAG default (E ≤ 10, roughly 8 bits) would report such hits, so for these three KOs this run is slightly more conservative at noise level.
- **E-values are computed with Z = 3,911,146**, about 10³ × DRAM's per-MAG Z. The best-KO ranking and the 2-significant-digit tie structure are unchanged; only the printed E-values are larger.
- **The sibling `bakta_test` track is also running hmmsearch.** Its pipeline (`bakta_test/scripts/pipeline.sh`, session 3646472) started 39 single-CPU hmmsearch jobs at about 11:55 CST. They search profile chunks against its own 53,219 test proteins, so they don't duplicate this run and should finish within about an hour. The user then has about 79 single-core jobs on poplar for a while. This track stays within its own 40/48 budget, but expect the ETA from `progress.sh` to wobble.
- `scripts/progress.sh` writes a fixed `/tmp/afr_kofam_cumleng.tsv` and greps the growing `.out` files (≈2.6 GB per chunk at the end) over NFS, so it will get slower near the end. Run it one copy at a time.

**Do not edit `scripts/run_chunk.sh` or `scripts/run_all.sh` on poplar until `ALL_DONE` exists.** bash reads a running script lazily, by byte offset, so editing either file in place would make the 40 waiting `run_chunk.sh` processes, or the driver, run whatever text now sits at their offset when hmmsearch returns. The summary's "later edits to the scripts reach the running jobs" holds only for `postprocess.py`, which each chunk starts as a fresh Python process. `postprocess.py` is safe to edit (keep the poplar and local copies identical). If a bash script must change, write it under a new name.

Review files: `test/review/review_check_partial.py`, `test/review/review_pp_partial.sh`, `test/review/partial_pp_chunk01/` here. Poplar copies: `kofam/test/review_check_partial.py`, `kofam/test/review_pp_partial.sh`, `kofam/test/review_partial_pp/`; that directory is a scratch test and safe to delete.

## What was run, exactly

| step | what | where (poplar) |
|---|---|---|
| inputs | 1,324 MAG FASTAs from `reassembly_2026/work/checkm/<unit>/fasta/<name>.fna` (list = `comparison/mags_old_vs_new.tsv` rows whose `genome` is `new|<name>`), copied with `rsync -aL --partial -e "ssh -o BatchMode=yes"`. **md5 of all 1,324 verified on poplar** (`md5sum -c kofam/mags.md5` → all OK). Local manifest with md5, unit, binner and species cluster: `mag_manifest.tsv`. | `/scratch1/afreiburger/salterns_2026/mags/<name>.fna` |
| KOfam release | **ko_list.gz Last-Modified 2026-10-26 05:26 GMT** (917,028 bytes, md5 `697af5b3b22a919306a1e60f99cec1a6`; identical (same md5) to the local `annotation_choice/sources/kegg/kofam_ko_list_2026-10-26.gz` that the marker table's 2026-10 thresholds come from), **profiles.tar.gz Last-Modified 2026-10-27 01:41 GMT** (1,567,491,980 bytes, md5 `bc68ad2a52ed48812b5496d15c215a7d`, gzip -t OK). Downloaded once (16 parallel HTTP byte ranges, genome.jp serves ~0.4 MB/s per connection) by the sibling `bakta_test` track; this track only reads it. Release name used everywhere: **KOfam 2026-10-26**. The local 2026-05-28 profiles in `~/Documents/EmilyKin/kofam_new/` were **not** used. | `/scratch1/afreiburger/salterns_2026/kofam_db/` (shared, read-only for this track) |
| HMM DB | `scripts/build_db.py`: all profiles concatenated in K-number order into one file, as DRAM does (DRAM searches every profile; KOs with threshold `-` are searched but never called). Also writes `db/hmm_info.tsv` (KO, LENG, threshold, score_type, profile index, definition) and `db/db_check.txt` (profile ↔ ko_list consistency). | `kofam/db/kofam_all.hmm` |
| genes | `scripts/call_genes.py` via `scripts/run_prodigal.sh` (40 × `xargs -P`, nice 10): **Prodigal 2.6.3** (bioconda build `h577a1d6_11`, EmilyKin env) `prodigal -p meta -g 11 -q`, on contigs ≥1,000 bp (**every contig of the 1,324 MAGs is ≥1,000 bp**: 752,569 contigs, minimum 1,000, so nothing is dropped). Protein IDs `<MAG>|<contig>_<n>`. Gene table per MAG: gene, mag, binner (maxbin2/metabat2/concoct from the name), contig, contig_len, start, end, strand, partial, aa_len, rbs_motif. | `kofam/genes/<MAG>.{faa,gff,genes.tsv,done}` |
| chunks | `scripts/make_chunks.py`: 40 chunks of **whole MAGs**, balanced by residues (greedy LPT). The 8 re-assembly control bins of `annotation_choice/kofam_control_test` (`meta/controls_chunk01.txt`) are forced into chunk 01. | `kofam/chunks/chunk_XX.faa`, `chunk_assignment.tsv`, `chunk_stats.txt` |
| search | `scripts/run_chunk.sh XX` for XX = 01..40, launched together by `scripts/run_all.sh` (detached): `nice -n 10 hmmsearch --cpu 1 --noali -T 14.7 --domT 14.7 -Z <total proteins> --tblout … --domtblout … -o … kofam_all.hmm chunk_XX.faa` (**HMMER 3.4**). | `kofam/hmm/chunk_XX.{tblout,domtblout,out}[.gz]`, `.done/.failed` |
| postprocess | `scripts/postprocess.py chunk <W> XX` runs automatically at the end of each chunk; `postprocess.py combine <W>` at the end of `run_all.sh` (stdlib-only Python 3.11 from the EmilyKin env). | `kofam/results/` |

### The KO rule (DRAM 1.5, read from `annotation_choice/sources/dram1_v1.5.0/mag_annotator_utils.py`)

- A (gene, KO) hit is significant when **score > threshold**, with the full-sequence score if ko_list `score_type` is `full`, or the domain score if `domain` (DRAM tests each domain row of the domtblout; the postprocessor uses the hit's best domain score, which is equivalent).
- **KOs whose threshold is `-` (1,870 of 28,504 in 2026-10-26; incl. *mpnS* K18049, *phnY\** K23996, *phnZ1* K05915) are searched but never called** — DRAM 1.5's `sig_scores` skips them. (DRAM2 would call every hit if its preformatted ko_list stores the dash as blank; not followed.) Their hits are kept in `*.no_threshold_hits.tsv.gz` (E ≤ 1e-5) and, for the three marker KOs, in `*.marker_hits.tsv.gz` at any score, for a manual cutoff later.
- **Best KO per gene = the significant KO with the lowest full-sequence E-value** (DRAM 1.5 `generic_hmmscan_formater`, also DRAM2's first sort key); ties (E-values that underflow to 0 for long proteins) broken by the higher full score, then KO id. DRAM itself leaves such ties to pandas' unstable sort.
- **All passing KOs** of a gene are also kept (cf. DRAM2's `kofam_ids`; DRAM2 uses ≥ rather than >, which only differs at exact equality).
- Reporting thresholds: DRAM 1.5 runs `hmmsearch --domtblout` with default reporting (E ≤ 10) per MAG. Pooling 3.9 M proteins would make E ≤ 10 correspond to ~18–20 bits, above the lowest KO thresholds (minimum 14.70; 83 KOs < 25). So hmmsearch reports by score instead: `-T 14.7 --domT 14.7` = the lowest numeric threshold in this ko_list, which guarantees every hit that could pass any KO is written. `-Z` is fixed to the total protein count of all 1,324 MAGs, so E-values do not depend on chunking (and the E-value ranking of KOs within a gene does not depend on Z at all). HMMER's acceleration filters (P-value based, Z-independent) are left at defaults, as in DRAM.
- Cross-check per chunk: number of (gene, KO) pairs in tblout = number aggregated from domtblout (`chunk_XX.qc.txt`, `pairs_match`).

### Guild rules (revised, `annotation_choice/README.md` §6.3), per MAG

Counted on each gene's **best KO** (primary columns) and on **any passing KO** (`*_anypass` columns):

| guild | rule | also reported |
|---|---|---|
| `methanogen` | *mcrA* K00399 | |
| `mepn_consumer` | *phnJ* K06163 | `phn_accessory_counts` (phnG/H/I/K/L/M), `phnGHIJKLM_n_distinct`, `phnGHIJKLM_n_contigs` (>1 = locus split across contigs) |
| `betaine_reducer` | *grdH* K21579 **or** *grdI* K21578 | `grdACD_n_distinct` (support) |
| `methylamine_user` | *mttB* K14083 / *mtbB* K16178 / *mtmB* K16176 **only in an *mcrA* genome** | `mttB_family_without_mcrA` (incl. *mtgB* K25220): the non-pyrrolysine MttB family, not methylamine use |
| `tma_producer` | *cutC*+*cutD* (K20038+K20037) **or** *cntA*+*cntB* (K22443+K22444) | `tma_candidate_one_sided` (e.g. *Roseovarius* CntA without CntB) |

MaxBin2 short-contig flag: every marker gene (any of the 68 in `meta/marker_ko_table.tsv`) on a **MaxBin2 contig < 2,500 bp** is listed in `short_maxbin2_marker_genes`; `guilds_lost_without_short_maxbin2_contigs` names guild calls that would disappear at DRAM's default 2,500 bp minimum. Per-hit flag `short_maxbin2_contig` in the marker-hit table.

## Outputs (poplar `kofam/results/`; copy back when done)

| file | content |
|---|---|
| `per_mag/<MAG>.kofam.tsv.gz` | every gene of the MAG: coordinates, contig length, best KO, score, threshold, score type, E-value, HMM coverage, KO definition, all passing KOs (DRAM-like annotation table) |
| `chunks/chunk_XX.genes_ko.tsv.gz` → `all_genes_ko.tsv.gz` | genes with ≥1 passing KO: best KO + all passing KOs with scores |
| `chunks/chunk_XX.mag_ko_counts.tsv` → `all_mag_ko_counts.tsv` | long MAG × KO table: `n_best`, `n_anypass` |
| `chunks/chunk_XX.mag_markers_guilds.tsv` → `all_mag_markers_guilds.tsv` | one row per MAG: guild calls, flags, `n_<marker>` best-KO counts for all 68 markers |
| `chunks/chunk_XX.marker_hits.tsv.gz` → `all_marker_hits.tsv.gz` | every reported hit (score ≥ 14.7) to a marker profile, passing or not, with the gene's best KO, HMM coverage, contig length and the short-MaxBin2 flag (near misses such as the Sec-truncated *grdH* are visible here) |
| `chunks/chunk_XX.no_threshold_hits.tsv.gz` | hits (E ≤ 1e-5) to the 1,870 KOs without a threshold |
| `chunks/chunk_XX.qc.txt`, `combine_status.txt` | per-chunk counts and the tblout/domtblout cross-check |

## How to check status

```bash
ssh -o BatchMode=yes poplar 'cd /scratch1/afreiburger/salterns_2026/kofam
  ls hmm/*.done 2>/dev/null | wc -l; ls hmm/*.failed 2>/dev/null          # chunks finished / failed (of 40)
  for c in $(seq -w 1 40); do echo -n "$c "; grep -c "^Query:" hmm/chunk_$c.out 2>/dev/null || echo gz; done | column -c 160   # profiles done per chunk (of 27,981)
  tail -3 logs/run_all.log; ls ALL_DONE 2>/dev/null; uptime'
```

Easier: `ssh -o BatchMode=yes poplar 'bash /scratch1/afreiburger/salterns_2026/kofam/scripts/progress.sh'` prints, per chunk, the profiles done, the fraction of total profile length (ΣLENG) done, the elapsed time and a linear ETA (or `done` / `FAILED`). The driver is alive while `kill -0 447804` succeeds; finished chunks have `hmm/chunk_XX.done` and `results/chunks/chunk_XX.postprocessed` (a failed postprocess leaves `results/chunks/chunk_XX.postprocess_failed`; the raw tables stay, re-run it by hand).

Progress of a running chunk = `grep -c '^Query:' hmm/chunk_XX.out` / number of profiles. Profiles are searched in K-number order, which is not uniform in length; `db/hmm_info.tsv` has the cumulative LENG for a better fraction.

While the jobs run, don't edit `scripts/run_chunk.sh` or `scripts/run_all.sh` in place. bash reads them lazily, so an edit changes what the running processes execute (see "Independent review" above). `postprocess.py` is safe to edit.

Restart after a failure or reboot of poplar: `cd /scratch1/afreiburger/salterns_2026/kofam; nohup setsid bash scripts/run_all.sh > logs/run_all.nohup 2>&1 < /dev/null &` — chunks with a `.done` file are skipped; a killed chunk restarts from the beginning (delete its partial `hmm/chunk_XX.*` first). Re-run one chunk's postprocessing with `…/bin/python3 scripts/postprocess.py chunk $PWD XX`.

Detect a chunk whose postprocessing failed (the search is still marked `.done` and `combine` skips it, so `results/combine_status.txt` lists fewer than 40 chunks): `ls results/chunks/*.postprocess_failed`. The full 27,981-profile postprocessing has not run yet; it was tested on chunk 01 with 68 and 368 profiles.

Copy results back (small files only; the pooled gene table `kofam/genes_all.tsv.gz`, 69 MB, is on poplar only):

```bash
rsync -a -e "ssh -o BatchMode=yes" poplar:/scratch1/afreiburger/salterns_2026/kofam/results/ \
  ~/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/annotation_runs/kofam/results/
```

## Files here

| file | what |
|---|---|
| `scripts/` | everything run on poplar (identical copies in `poplar:…/kofam/scripts/`) |
| `mag_manifest.tsv`, `mags.md5` | the 1,324 MAGs: unit, binner, species cluster, CheckM, local path, md5 |
| `meta/marker_ko_table.tsv` | copy of `annotation_choice/marker_ko_table.tsv` (68 markers) used by the postprocessor |
| `meta/controls_chunk01.txt` | the 8 control bins forced into chunk 01 |
| `test/markers/`, `test/subset/`, `test/timing/` | chunk-01 marker-profile test, 368-profile full-code-path test, timing runs (copied back from poplar); `test/control_comparison.tsv` = output of `scripts/compare_controls.py` |
| `poplar_copies/` | `chunk_stats.txt`, `chunk_assignment.tsv`, `chunk_pids.tsv`, `hmmsearch_pids.tsv`, `db_check.txt`, the shared KOfam download script and md5s |
| `scripts/compare_controls.py` | local (py_venv) comparison of a chunk-01 result with `annotation_choice/kofam_control_test` |
| `test/review/` | independent review, 2026-11-11: `review_check_partial.py` (tblout vs domtblout, postprocess rule vs DRAM `sig_scores`, on the partial real output), `review_pp_partial.sh` (postprocess.py on chunk 01's truncated real full-DB tables), `partial_pp_chunk01/` (its outputs and `/usr/bin/time`) |
