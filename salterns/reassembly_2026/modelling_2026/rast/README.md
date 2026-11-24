# RAST genomes of the 57 modelling species (salterns 2026 re-assembly)

Written by `rast_driver.py readme` from `ledger.json` (2026-11-20T19:23:00+00:00). KBase narrative ws 275645 (`freiburgermsu:narrative_1793847843029`, "Salterns metagenomes: reassembly and binning (2026)"). Write approved by Andrew on 2026-11-20 for these `m26_*` objects only.

**Status:** 0 of 57 genomes annotated and downloaded; 0 failed.

**In progress** (state counts {'uploading': 2, 'new': 55}). A detached watcher (`watch --scope auto --apply`, output in `watch.log`) polls KBase every 10 min: pilot end to end, automatic pilot check, then the other 55. Check with `python rast_driver.py status`; re-write this file with `python rast_driver.py readme`. To resume after a stop, start the same watch command again (the ledger is resumable; an interrupted submit is reconciled with EE2 before anything is resubmitted). KBase queues were saturated when this started (15-28 h waits on njs and bigmem for this account's jobs, 2026-11-20).

## Inputs

- Genomes: `../../comparison_rollout/modelling_species.tsv`, the species representatives of every MePn-consumer, betaine-reducer, methanogen or methylamine-user species of the final salterns MAG set (31 rollout, 26 co-assembly; 50 Bacteria, 7 Archaea).
- FASTAs: `../../work/rollout_mags_kept/<mag_id>.fa` and `../../work/coasm_mags_kept/<mag_id>.fa`; each md5 equals `../../rollout_final/mags_kept.md5` / `../../coasm_final/mags_kept.md5`.
- Lineage: GTDB-Tk 2.8.0 / GTDB R232, `../../comparison_rollout/mags_three_way.tsv` (`domain`, `classification`).

## What was run

1. **Upload.** Each FASTA was POSTed to the KBase Blobstore (node md5 checked against the local file), then one EE2 job per batch ran `AssemblyUtil.save_assemblies_from_fastas` (AssemblyUtil 3.0.0, tag release, git `59855bcb24014c3030e146b98c8d3e70abac1fa1`) with `workspace_id` 275645 and per input `node` + `assembly_name` `m26_<mag_id>.assembly`; no `min_contig_length`, `type` left at AssemblyUtil's default `Unknown` (as the project's other MAG Assemblies). Batches: the 2 pilot genomes, then the other 55. Each saved Assembly was checked against its FASTA: version 1 of the expected name, contig ids, lengths and md5s, `num_contigs`, `dna_size`, and the AssemblyUtil assembly md5 (md5 of the sorted contig md5s). The Blobstore nodes of step 1 are left in place (AssemblyUtil stores its own copy); they hold only the FASTAs and are not referenced by any object.
2. **RAST.** One EE2 job per genome, app `RAST_SDK/annotate_genome_assembly` ("Annotate Genome/Assembly with RASTtk - v1.073"), method `RAST_SDK.rast_genome_assembly`, RAST_SDK 1.9.5, tag release, git `7171090d87fccc8b7ecf1a1d02398995dcc2dd45` (pinned; the genomes' provenance records service_ver {}). Parameters, as the narrative app sends them except `create_report`:
   - `object_ref`: the verified Assembly UPA; `output_workspace`: the narrative's workspace name; `output_genome_name`: `m26_<mag_id>.RAST`;
   - `genetic_code`: 11 (GTDB-Tk reported table 11 for 23 genomes and none for the 34 placed by ANI screen; no genome uses 4 or 25);
   - `domain`: `A` for the 7 Archaea, `B` for the 50 Bacteria (the app's dropdown values; GTDB R232 domain);
   - `scientific_name`: from GTDB R232: the species name where GTDB names one (34 genomes), else `<genus> sp.` (22), else `<family> bacterium` (1). No `ncbi_taxon_id`, so RAST takes the name as given and no Relation Engine lookup is made;
   - `create_report`: 0. The narrative app sends 1, which also saves a `KBaseReport.Report` under a generated name; 0 keeps the writes to the `m26_*` objects. It changes no annotation step;
   - annotation workflow: RAST_SDK's default for this app (gene calling, rRNA/tRNA/repeat/CRISPR/selenoprotein/pyrrolysoprotein calls, k-mer v2 and similarity-based functions, overlap resolution); this app exposes no workflow options.
   - Concurrency: at most 10 live (queued + running) RAST jobs (client group bigmem, shared with the GROW pipeline). First submitted , last finished .
3. **Download.** `Workspace.get_objects2` of each genome; `genomes/<mag_id>.RAST.json` is the object's `data` (features with locations, translations and `functions`, `cdss`, `mrnas`, `non_coding_features`, `ontology_events`, ...), the same layout as `jplfaria/Cliff MAGs/RASTified_MAG_json/*.RAST.json`. Each file was checked: domain, scientific name, genetic code 11, `assembly_ref` = the uploaded Assembly, contig count and size = the FASTA, protein features present, RAST functions present, and an SSO ontology event.

## Pilot (1 bacterium, 1 archaeon, end to end before the rest)

## Counts per genome

`fn` = protein features with a RAST function; `non-hyp` = of those, not only `hypothetical protein`; `ncRNA/repeat` = non-coding features of type rna / repeat.

| MAG | set | domain | scientific name (GTDB R232) | Assembly | RAST genome | contigs | features | CDS | fn | non-hyp | ncRNA/repeat | state |
|---|---|---|---|---|---|---:|---:|---:|---:|---:|---|---|
| R1_A_D1__concoct.bin.007 | rollout | A | Methanolobus sp020055655 |  |  | 104 |  |  |  |  | / | uploading |
| R1_B_D1__maxbin2.Bin.073 | rollout | A | Natrinema hispanicum_A |  |  | 411 |  |  |  |  | / | new |
| R2_A_D2__metabat2.bin.055 | rollout | A | Haloplanus sp. |  |  | 117 |  |  |  |  | / | new |
| R2_A_H2O__metabat2.bin.026 | rollout | A | Halorubrum sp. |  |  | 92 |  |  |  |  | / | new |
| R2_B_H2O__metabat2.bin.051 | rollout | A | Haloquadratum walsbyi |  |  | 171 |  |  |  |  | / | new |
| R2_C_H2O__metabat2.bin.005 | rollout | A | Natronomonas sp. |  |  | 48 |  |  |  |  | / | new |
| R2_C_H2O__metabat2.bin.008 | rollout | A | Halonotius sp. |  |  | 84 |  |  |  |  | / | new |
| R1_A_D2__maxbin2.Bin.123_sub | rollout | B | BM714 sp035269385 |  |  | 313 |  |  |  |  | / | new |
| R1_A_D2__metabat2.bin.020 | rollout | B | DATFJG01 sp035269305 |  |  | 112 |  |  |  |  | / | new |
| R1_A_D2__metabat2.bin.053 | rollout | B | Roseovarius sp. |  |  | 228 |  |  |  |  | / | new |
| R1_B_D1__metabat2.bin.017_sub | rollout | B | Roseovarius sp964570675 |  |  | 124 |  |  |  |  | / | new |
| R1_B_D1__metabat2.bin.040 | rollout | B | Roseovarius sp035268705 |  |  | 198 |  |  |  |  | / | new |
| R1_B_D2__concoct.bin.112_sub | rollout | B | Halanaerobium sp. |  |  | 443 |  |  |  |  | / | new |
| R1_B_D2__metabat2.bin.065 | rollout | B | Roseovarius sp. |  |  | 150 |  |  |  |  | / | new |
| R1_B_H2O__concoct.bin.053 | rollout | B | Roseovarius sp. |  |  | 418 |  |  |  |  | / | new |
| R1_C_D1__concoct.bin.031 | rollout | B | UBA8670 sp964393495 |  |  | 253 |  |  |  |  | / | new |
| R1_C_D1__maxbin2.Bin.003_sub | rollout | B | M55B127 sp018334565 |  |  | 101 |  |  |  |  | / | new |
| R1_C_D2__metabat2.bin.048 | rollout | B | Marivita sp020055675 |  |  | 58 |  |  |  |  | / | new |
| R1_sediment__concoct.bin.036 | coassembly | B | CAZRZN01 sp964604565 |  |  | 562 |  |  |  |  | / | uploading |
| R1_sediment__concoct.bin.052_sub | coassembly | B | CAYRVO01 sp964427565 |  |  | 280 |  |  |  |  | / | new |
| R1_sediment__concoct.bin.104 | coassembly | B | Roseovarius sp031580875 |  |  | 395 |  |  |  |  | / | new |
| R1_sediment__maxbin2.Bin.094_sub | coassembly | B | JAWLGC01 sp041685265 |  |  | 249 |  |  |  |  | / | new |
| R1_sediment__metabat2.bin.121_sub | coassembly | B | Roseovarius sp003298505 |  |  | 123 |  |  |  |  | / | new |
| R1_sediment__metabat2.bin.142_sub | coassembly | B | RI-231 sp. |  |  | 867 |  |  |  |  | / | new |
| R1_sediment__metabat2.bin.186 | coassembly | B | DAOWES01 sp. |  |  | 137 |  |  |  |  | / | new |
| R1_sediment__metabat2.bin.206 | coassembly | B | Halanaerobium sp. |  |  | 149 |  |  |  |  | / | new |
| R1_sediment__metabat2.bin.258 | coassembly | B | JACKII01 sp964474065 |  |  | 117 |  |  |  |  | / | new |
| R1_sediment__metabat2.bin.356 | coassembly | B | QPJE01 sp964559325 |  |  | 86 |  |  |  |  | / | new |
| R2A_sediment__concoct.bin.151 | coassembly | B | WSZY01 sp. |  |  | 220 |  |  |  |  | / | new |
| R2A_sediment__concoct.bin.189_sub | coassembly | B | UBA8366 bacterium |  |  | 627 |  |  |  |  | / | new |
| R2A_sediment__metabat2.bin.132 | coassembly | B | CAMCWB01 sp. |  |  | 441 |  |  |  |  | / | new |
| R2_A_D1__concoct.bin.120 | rollout | B | Saliniramus sp035311605 |  |  | 191 |  |  |  |  | / | new |
| R2_A_D1__maxbin2.Bin.018_sub | rollout | B | M55B127 sp035312025 |  |  | 203 |  |  |  |  | / | new |
| R2_B_D2__maxbin2.Bin.065 | rollout | B | Roseovarius sp035269125 |  |  | 126 |  |  |  |  | / | new |
| R2_B_D2__metabat2.bin.006_sub | rollout | B | Roseovarius sp024742655 |  |  | 233 |  |  |  |  | / | new |
| R2_B_D2__metabat2.bin.056 | rollout | B | M55B127 sp035269185 |  |  | 105 |  |  |  |  | / | new |
| R2_B_D2__metabat2.bin.063 | rollout | B | M55B127 sp041642015 |  |  | 342 |  |  |  |  | / | new |
| R2_C_D1__metabat2.bin.057_sub | rollout | B | Desulfotignum sp035268905 |  |  | 349 |  |  |  |  | / | new |
| R2_C_D2__concoct.bin.063 | rollout | B | Roseibacterium sp964578115 |  |  | 58 |  |  |  |  | / | new |
| R2_C_D2__concoct.bin.072 | rollout | B | JAVJPJ01 sp040217295 |  |  | 87 |  |  |  |  | / | new |
| R2_C_D2__metabat2.bin.065 | rollout | B | Halanaerobium saccharolyticum |  |  | 417 |  |  |  |  | / | new |
| R2_restored_C_black__maxbin2.Bin.042 | rollout | B | Desulfotignum sp035268565 |  |  | 140 |  |  |  |  | / | new |
| R2_restored_C_black__metabat2.bin.011_sub | rollout | B | Roseovarius sp. |  |  | 431 |  |  |  |  | / | new |
| R2_restored_D_shore__concoct.bin.011_sub | rollout | B | JASFMS01 sp. |  |  | 302 |  |  |  |  | / | new |
| R2_sediment__concoct.bin.036_sub | coassembly | B | Roseovarius sp965656855 |  |  | 434 |  |  |  |  | / | new |
| R2_sediment__concoct.bin.112 | coassembly | B | CAJQQK01 sp964632265 |  |  | 138 |  |  |  |  | / | new |
| R2_sediment__concoct.bin.183_sub | coassembly | B | JBKDZR01 sp. |  |  | 270 |  |  |  |  | / | new |
| R2_sediment__metabat2.bin.112_sub | coassembly | B | Lentibacter sp. |  |  | 428 |  |  |  |  | / | new |
| R2_sediment__metabat2.bin.128 | coassembly | B | DATHSN01 sp. |  |  | 527 |  |  |  |  | / | new |
| R2_sediment__metabat2.bin.141 | coassembly | B | PWPR01 sp964497955 |  |  | 135 |  |  |  |  | / | new |
| R2_sediment__metabat2.bin.154 | coassembly | B | JAAAOG01 sp. |  |  | 590 |  |  |  |  | / | new |
| R2_sediment__metabat2.bin.162 | coassembly | B | M55B127 sp964642575 |  |  | 162 |  |  |  |  | / | new |
| R2_sediment__metabat2.bin.362_sub | coassembly | B | M55B127 sp913055675 |  |  | 565 |  |  |  |  | / | new |
| R2_sediment__metabat2.bin.405 | coassembly | B | Thermohalobaculum sp. |  |  | 408 |  |  |  |  | / | new |
| R2_sediment__metabat2.bin.414_sub | coassembly | B | DAINOF01 sp. |  |  | 568 |  |  |  |  | / | new |
| R2_sediment__metabat2.bin.434_sub | coassembly | B | UBA8670 sp046482875 |  |  | 394 |  |  |  |  | / | new |
| R2_sediment__metabat2.bin.435 | coassembly | B | Halonatronomonas sp. |  |  | 254 |  |  |  |  | / | new |

## Jobs, failures and retries

- Upload job `6b009653b3a2d40340bd94aa` (2 assemblies): queued, submitted 2026-11-20T18:49:53+00:00, finished 
- RAST: no failed or retried jobs.

## Files

- `rast_driver.py`: the driver (plan / advance / pilot-check / watch / status / readme).
- `ledger.json`: per genome the FASTA md5, Blobstore node, upload job, Assembly ref and checks, RAST attempts (job ids, params, timestamps), genome ref, app version, download checks and counts; `events` log.
- `genomes/<mag_id>.RAST.json`: the 57 genome objects (Workspace `data`).
- `logs/`: workspace snapshot before the first write, EE2 logs of the pilot RAST jobs and of any failed job.
- `watch.log`: the detached watcher's output.
