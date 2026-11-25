# Unbinned contigs: identification, merging and prGEMs (design, 2026-11-25)

Goal (Andrew, 2026-11-25): run the unbinned contigs of the 781 assemblies through prFBA. Identify each one through the
existing MAG framework, then either merge it into a binned representative of the same taxon, or, if it is unique,
create a synthetic genome and prGEM for it.

## Inputs
- **Assemblies:** the 781 MEGAHIT assemblies (KBase narrative 261700), with full Bakta annotation (coding and
  non-coding) on poplar. This is the same caller as the MAGs (Andrew's choice).
- **Unbinned contigs:** per assembly, the contigs of at least 2,500 bp (the binning minimum) that are in none of that
  assembly's admitted MAGs. Contig membership comes from the extracted MAG FASTAs (`<bin>__<run>.fna`).
- **Contig floor (checked 2026-11-25):** the KBase MEGAHIT app ran with `min_contig_length = 2000` (provenance), so
  the assemblies hold no contig under 2 kb; every MAG contig is >= 2,500 bp. Contigs of 2.0-2.5 kb are in the
  assemblies but never binnable. Contigs of 1-2 kb exist only in re-assemblies with a 1 kb floor. The threshold
  sensitivity test (`threshold_*.py`, `reassembly/`) decides whether the pipeline should use them.
- **Reference framework:** the 13,164 MAGs, their census species clusters and the 1,206 annotated representatives.

## Identification (Andrew's choice: match to our MAGs first)
1. Match each unbinned contig against MAGs of the **same sample** first, then against **any** MAG. Use skani where the
   contig is long enough, otherwise a protein-level best-hit vote of its Bakta proteins against MAG proteomes.
2. Contigs with no MAG match get GTDB R232-based taxonomy and are candidates for a synthetic genome.

## Abundance rule (Andrew, 2026-11-25): total abundance per sample is preserved
- **Cross-sample matches identify; they never move abundance.** A match to a MAG from another sample may set the
  contig's taxon or species cluster, and nothing more.
- **Merging is same-sample only.** A contig's abundance (read coverage) may be added only to a binned organism, a MAG,
  from the **same sample** whose species cluster it matches.
- **Otherwise the contig stays in its own sample.** With no same-sample binned organism of that taxon, the contig's
  abundance remains in its sample under a synthetic genome / prGEM for that taxon in that sample. It is never pooled
  into another sample's organism or into a cross-sample representative.
- **Check:** for every sample, the sum of abundance over (MAGs + merged contigs + synthetic genomes + still-unassigned
  contigs) must equal that sample's total before the step. This is asserted in code.

## prFBA / prGEM
- **Merge into a same-sample MAG:** the contig's genes join that MAG's genome as extra features, and its model is
  rebuilt as a prGEM. Reactions from contig-only genes get a probability below 1, reflecting assignment confidence.
- **Unique contigs:** contigs of one taxon in one sample are grouped into a synthetic genome, and a prGEM is built
  from it. Built with ModelSEEDpy (msprobability / MSGenome loaded from Bakta calls), reusing the parts of the salterns
  `prfba_augmentation` prototypes that work.

## Open details (to settle on the first samples)
- The match thresholds (contig length for skani vs protein vote; identity and vote share).
- The abundance source: read mapping (CoverM or sylph on the sample's reads) vs MEGAHIT k-mer coverage.
