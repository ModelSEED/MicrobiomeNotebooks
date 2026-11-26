# MAG_contigs — original contigs for all 450 salterns MAGs

One FASTA per MAG, named `<KBase Assembly name>.fna`, fetched by `../download_mag_contigs.py`
on 2026-10-19. `index.tsv` maps each MAG to its Assembly ref, FASTA handle, contig count, bp,
md5 and source. The FASTA files are gitignored; this README and `index.tsv` are not.

## Where these come from

| Set | KBase narrative | How the Assemblies were made |
|---|---|---|
| 310 dereplicated | [Metabolism 310 Bins](https://narrative.kbase.us/narrative/77580) (ws 77580, cliff_bueno) | `kb_uploadmethods.batch_import_assemblies_from_staging`, 2020-11-18, FASTA files from Cliff's staging area (`/cliff_bueno`), `min_contig_length=500`, AssemblySet `SSFB_Assembly_Set` |
| 140 unretained | [Unretained MAGs](https://narrative.kbase.us/narrative/163264) (ws 163264, jplfaria) | `kb_uploadmethods.import_fasta_as_assembly_from_staging`, 2023-11-17, FASTA files from a Google Drive download (`drive-download-20231117T205750Z-001/*.contigs.fa`), `min_contig_length=500`; v2 re-saved 2023-11-19 by `kb_gtdbtk` with identical content (md5, bp, contig count) |

Nothing in KBase sits upstream of these Assemblies: they were uploaded as FASTA files, and the
DAS Tool binning was not done in any narrative readable by `freiburgermsu`. Cliff's inventory
(`../data/Files Summary.xlsx`) lists the original bin FASTA files as the "Dereplicated MAGs" and
"Unretained MAGs" folders in the Google Drive folder "JGI Salterns Modeling" (owner
sgtringe@lbl.gov); those folders were not visible to `freiburgermsu` on 2026-10-19.

Every downstream genome — the 2020 RAST genomes in ws 77580, the 2023 copies and the `.dlt`
re-annotations in ws 155805 (Saltern MAG Models), and the RAST genomes in ws 163264 — references
these same Assembly objects.

## Are they unaltered?

- All 145,616 contigs keep their metaSPAdes names (`NODE_<n>_length_<L>_cov_<c>`), and every
  sequence length equals the `<L>` in its name.
- Every MAG's shortest contig is at least 2,500 bp (the binning cutoff), so the 500 bp import
  filter removed nothing.
- 130 of the files are hard links to `../BERDL_MAGs/contigs/` (the same Assembly objects).
