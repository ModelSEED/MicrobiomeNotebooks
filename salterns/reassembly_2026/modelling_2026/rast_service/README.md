# RAST functions for the 57 modelling genomes (RAST annotation service)

Run on 2026-11-20 at Andrew's choice ("RastClient, cancel KBase"), replacing the KBase RASTtk route (`../rast/CANCELLED.md`).

**Genes.** Prodigal 2.6.3 (`-p meta -g 11`) proteins from the KOfam tracks (`sources.tsv`; copies in `faa_src/`), so RAST and KOfam annotate identical genes (same protein ids).
- 31 genomes come from `kofam_rollout`; 26 co-assembly representatives come from the first-pass `kofam` track.

**Service.** ModelSEEDpy 0.4.3 `RastClient` (Andrew's fork, `~/Documents/ModelSEEDpy`) calls `https://tutorial.theseed.org/services/genome_annotation`, `GenomeAnnotation.run_pipeline`, with the stages:
- `annotate_proteins_kmer_v2`;
- `annotate_proteins_kmer_v1`;
- `annotate_proteins_similarity`.

These are RAST's function-assignment stages. KBase's RASTtk app adds gene calling and RNA/repeat stages, which model building doesn't use.

**Command.** `~/Documents/py_venv/bin/python annotate_rast_service.py`, which is resumable. Output:
- `functions/<mag_id>.tsv`: `protein_id`, `function`; an empty function means no RAST call;
- `summary.tsv`: per-genome counts.

**Result.** 57/57 genomes in about 2 s each.

| | Mean | Median | Min | Max |
|---|---|---|---|---|
| Proteins | 3,165 | 3,256 | 1,439 | 4,683 |
| With a RAST function | 2,393 | 2,314 | 1,015 | 3,786 |
| Non-hypothetical | — | 58 % | 34 % | 79 % |

For comparison, the original 310 MAGs' KBase RAST genomes had 61–73 % of features with a function. This set includes 7 archaea and many novel lineages.
