# Weisenhorn wetland × KBase-2 MAG integration: handoff review and implementation (2026-11-13)

This review covers Chris Henry's handoff (`HANDOFF.md`, 2026-11-11, run `wetland-mag-integration-20261003T053501Z`). Six
read-only agents checked each claim against primary sources, and a second agent tried to refute each finding. The sources were:
- the KBDL service's source (poplar `KBDLJobRunningPrototype`, HEAD 3568593) and its live API;
- KBase Workspace and EE2;
- KBUtilLib;
- the local GTDB-Tk R232, CheckM and DAS Tool results.

Three reviewers then critiqued the method. The raw per-claim evidence is in the workflow journal (`~/.claude/projects/-home-freiburger-Documents-MicrobiomeNotebooks-wetland/c6b67c5a-c8e2-573a-9763-d5f7f4c5c46d/subagents/workflows/wf_bccd1ab9-687/journal.jsonl`). This
file records the conclusions that change what gets run.

## Corrections to the handoff (verified)
| Handoff says | Actually (2026-11-13) | Consequence |
|---|---|---|
| Build halted: KBDL `build_genome` calls CDS only; RNA only via RAST annotation, which also assigns functions | **Outdated.** Since 2026-10-16, `KBDLBuildGenome` calls RNA genes by default (`call_non_coding=true`). Two callers: `inprocess` (pyrodigal + pyaragorn tRNA/tmRNA + pybarrnap rRNA) and `bakta`, structural only (tRNAscan-SE, rRNA, ncRNA via Rfam, CRISPR). Neither assigns protein function. `KBDLGenomeAnnotation` calls no genes at all. | The halt can be lifted; Phase 1 becomes one choice of `gene_caller` (see Decisions) |
| Pass the SKANI DB names `gtdb_r220_consolidated`, `novel_hq_v2` | **Wrong.** The service resolves objects by id only. Use `ref:cee67d7f…` (GTDB r220) and `ref:1be8415b…` (novel_all_representatives), or `ref:9f4ac5e9…` (novel_hq_v2). | Raw submits with `ref:` ids |
| `novel_hq_v2` = the cumulative novel-genome DB, 7,443 genomes | **Outdated.** It holds 7,442 isolate genomes (AllTheBacteria / BV-BRC) and **no environmental MAGs**. It was replaced on 2026-10-05 by `novel_all_representatives` (78,132 representatives: 64,161 SPIRE MAGs, 9,844 MGnify genomes, 4,116 AllTheBacteria and 11 BV-BRC isolates). | Wrong novel-side reference for wetland MAGs; both are reported below |
| GTDB side = skani vs GTDB r220 | r220 is two releases behind; all 12,859 MAGs already have **GTDB-Tk 2.8.0 / R232** lineages with the proper species rule (ANI ≥ radius **and** AF ≥ 0.5). KBDLSKANI drops AF. | GTDB side taken from R232; r220 is reported for comparison only |
| Census rule: best hit ≥ 95% ANI → species cluster; < 95% to both DBs → its own provisional cluster | There is no AF test, and no MAG-vs-MAG step, so cross-sample duplicates each become a "new species" (8,255 clusters vs ~6,000 with clustering) | Added AF ≥ 50% and de novo 95% clustering of species-novel MAGs |
| "Signal floor": ≥ 95% of MAGs mapped vs both DBs | Undefined, and it would fail for a scientific reason: for 34.8% of MAGs GTDB-Tk reports no R232 reference of any kind, and only ~51% have a species | To be replaced by technical gates plus positive/negative controls (recommended; the controls belong to a KBDL skani run, which has not happened). The census itself has technical checks only: 12,859 MAG files, 12,859 CheckM2 rows, GTDB-Tk summaries covering every MAG, partition checks, and exit code 0 for every skani and CheckM2 step. The share of MAGs with no R232 reference of any kind (34.8%) is reported in `RESULTS.md` as a result |
| No CheckM2 in KBase; ask Andrew for one real CheckM2 report | No released app; a dev-only `kb_CheckM2` 0.0.2 errored on 23 of 30 calls. But **KBDL has `KBDLCheckM2`** with a public CheckM2 DB, and poplar has CheckM2 1.1.0 (EmilyKin bin_qc env). | CheckM2 run on poplar for all MAGs; no narrative report needed |
| MAG collection still forming; DAS Tool / extraction stalled at 147 / 9 sets / 51 MAGs | Nothing stalled. Jobs were queued behind 1,833 binner jobs at the account's ~10-job cap. Today: 709 DAS Tool runs, **12,859 MAGs from 702 of 770 assemblies**, all with R232 lineages; collection is final under the current settings. | Census population frozen at 12,859 |
| Per-MAG FASTAs: ~101 extracted on the pod | All 12,859 FASTAs exist locally (md5-checked against KBase) and on poplar | No re-extraction |
| Landing: `skani_distances` APPEND, idempotent | An Iceberg append is not idempotent (the dedup key is unknown). Landing needs `kbaseincubator` read-write tenant membership, which `freiburgermsu` lacks (403), and runs from a BERDL pod. KBDLStoreLoad is unrelated. | Deferred (see Decisions); census results do not depend on it |
| KBDLBuildGenome attaches species taxonomy | Only for a GTDB r220 hit above 95% ANI, after the companion was wired on 10-31; untested since. Novel DBs have no companion. Content-hash dedup is **permanent and cross-owner**: the first build of a FASTA fixes its gene caller for everyone. | Decide gene_caller **before** any build; pilot first |

Other facts that matter:
- **Assemblies without MAGs:** 68 of 770 have none, almost all low-yield assemblies from deep libraries.
- **Never assembled:** 11 libraries of 31.6–41.8 Gbp. Their single MEGAHIT attempt never finished: 10 went log-silent on 07-30 and were terminated by the user on 08-14, and the 11th (SRR13021855) errored on 07-30 after 6.6 h ("Unexpected Job Error"). `failed_assemblies.csv` lists 41 failed MEGAHIT jobs: 30 earlier failures on other libraries, all recovered when retried, plus these 11 single attempts, which were never retried.
- **Never imported:** 619 runs in `Sample_metadata_wetland_only.csv`.
- **CheckM v1 is complete for all 12,859 MAGs.** The last 3 sets (the 248 newest MAGs) finished on 2026-11-13 19:19 after queuing behind GROW jobs under the account's ~10-job cap. The results were written to KBase and verified 12,859/12,859 (report 272312/20564/1).

## Implemented (no shared writes)
- **Phase 4, quality:** CheckM2 1.1.0 (uniref100.KO.1) on all 12,859 MAGs, on poplar. CheckM2 is authoritative; CheckM v1 and DAS Tool SCG are kept as columns.
- **Phase 3, census:** `census.py`.
  - **Known species** come from R232.
  - **Species-novel MAGs are clustered de novo** (greedy centroid at ANI ≥ 95 / max AF ≥ 50), with stable `WN_<hash>` IDs.
    - They are never moved into a known species through another MAG: the code review found that 109 of 120 such
      "member joins" contradicted GTDB-Tk's own test against the species representative.
    - Instead each one is labelled with the known species it borders and GTDB-Tk's reason for rejecting it: `below_radius`,
      `low_af` (usually a fragment of a known species), or `no_gtdbtk_value`.
    - `qc_known_novel_edges.tsv` lists every species-level and strain-level edge between a known and a new cluster.
  - **Each cluster is labelled** against `novel_all_representatives` and `novel_hq_v2`, and r220 is reported for comparison (skani 0.3.1 on poplar, AF kept).
  - **Strain clusters** at 99%.
  - **Sites** are sampling coordinates merged within 2 km.
  - **Compute:** `run_census_compute.sh`.
  - **Review:** an adversarial code review (12 findings, all fixed) checked the code; an independent recomputation reproduced every headline number of the preliminary census. That check predates the fixes and the CheckM2 re-run, so it does not cover the final numbers in `RESULTS.md`; those were then fact-checked against the outputs (326 claims, 33 corrections applied).
- **Phase 5, representatives:** `select_reps.py`.
  - Eligible: > 70% / < 10%, without GTDB-Tk's multiple-hit flag.
  - Score: completeness − 5×contamination + 0.5·log10 N50.
  - Selection: one MAG per 99% strain, up to 10 per cluster, preferring a new sampling site.
  - No two representatives, in any cluster, share a strain-level edge (≥ 99% ANI with max AF ≥ 50%). 35 representative pairs in different clusters are ≥ 99% ANI only at max AF < 50%. Known species are processed first.
  - Relaxation only for clusters with nothing eligible (one rep at ≥ 50/<10, flagged); otherwise the cluster is census-only.

RESULTS: see `RESULTS.md` (written after the CheckM2 and skani runs finish).

## Decisions needed (not implemented)
1. **`gene_caller` for KBDLBuildGenome (Phase 2).** It is permanent per FASTA content and shared across owners.
   - `bakta`: ncRNA, CRISPR and better CDS on fragmented MAGs; about 3–5 days of queue at the service's cap of 4 builds.
   - `inprocess`: tRNA, tmRNA and rRNA only; about 8 h of gene calling for all MAGs at the cap of 4. Every build also runs a skani search, which brings it to about 14.5–17.5 h, or about 1 day under current poplar load.
   - Never `auto`: it mixes the two callers according to the r220 hit.
2. **Build scope.** All 12,859 MAGs, as the handoff says, or only the representative worklist (~5k) that Phase 6 annotates. The census does not need builds.
3. **Novel-side reference for the record:** `novel_all_representatives` (recommended) or the handoff's `novel_hq_v2`.
4. **Landing to `kbaseincubator.genome_clearhouse.skani_distances`.** It needs:
   - `kbaseincubator` read-write tenant access, requested from a BERDL pod;
   - the `source_catalog` token and the dedup key, agreed with the table owner.
   It also needs a note to Chris before any shared write.
5. **Phase 6 compute route** for the representative subset, e.g. DRAM2 on h100 or KBDLGenomeAnnotation.
6. **Optional, upstream:** re-run MEGAHIT on the 11 never-assembled libraries, and decide on the 619 never-imported runs.
