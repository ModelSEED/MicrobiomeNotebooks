# Wetland MAG census: data scope (what is in, what is out, and why)

Recorded 2026-11-15. Decision (Andrew): **accept that 619 runs of the sample metadata were never imported; do not import them.**

## The sample universe
`wetland/Sample_metadata_wetland_only.csv` lists **1,400 runs** (Peatlands 884, Wetlands 423, Permafrost regions 93).

| Group | Runs | Status |
|---|---|---|
| Imported into KBase narrative 261700 | 781 | Paired-end read libraries. 770 assembled with MEGAHIT; the other 11 are being re-assembled (see below) |
| **Skipped: SPIRE already has MAGs for the BioSample** | **612** | Never imported (deliberate) |
| **Skipped: no INSDC accession** | **7** | Never imported (cannot be fetched from SRA/ENA) |

## Why the 619 were skipped
The import (2026-06-11, 773 `kb_uploadmethods/import_sra_as_reads_from_web` jobs, plus 8 pilot runs imported by Chris
Henry on 2026-06-04) followed one rule:

> **A run was skipped when its BioSample already has MAGs in SPIRE v1.**

The rule comes from the database-membership check in `wetland/db_membership/`: `mag_database_membership.csv` and
`MAG_database_membership_report.md`, committed a2bcf19 on 2026-05-29. Evidence that this rule, and nothing else, decided
the import:
- **It splits all 1,400 rows exactly:**
  - 612 SPIRE-covered runs: 612 skipped, 0 imported;
  - 781 non-covered runs with an SRA accession: 0 skipped, 781 imported;
  - 7 rows without an INSDC accession (DZ1B, DZ2M, H03B, QM01A_25-30, Z2-45-50, LH01, LH02; BioProject "This study";
    Hainan, Guangdong and Liaoning sites): skipped.
- **It holds by BioSample:** of 1,300 BioSamples, the 40 with several runs are never split between imported and skipped.
  It does **not** follow BioProject: 4 projects are split, e.g. PRJNA386568 100 imported / 198 skipped.
- **Data availability was not the reason.** All 612 skipped SRA runs are paired-end WGS metagenomes with FASTQ in ENA
  (608 also on the AWS SRA mirror).
- **Other columns only correlate with SPIRE coverage:** older release years (2016–2019), HiSeq 2000/2500, lower SRR numbers.
- **No written record of the rule was found** (no script, commit message or note). It is reconstructed from the data
  above with very high confidence. The motive is inferred: not re-deriving MAGs that SPIRE has already published for the
  same samples.

What this means for the census:
- **SPIRE holds 10,243 MAGs** (all ≥ 50% complete / ≤ 10% contaminated) from the 578 BioSamples behind the 612 runs.
- **Part of that coverage is already in the census,** through the novel-genome reference `novel_all_representatives`:
  2,140 of those SPIRE MAGs are in it, and 2,420 of their 3,436 SPIRE species clusters have a representative there.
- **Census statements are therefore about the 781 imported runs:** 12,859 MAGs from 702 of 770 assemblies. They are not
  about all 1,400 runs in the metadata.
- **Cost of not importing,** for reference: importing the 612 would add about 8.8 Tbp of reads, about 4 weeks of MEGAHIT
  plus binning, and about 10,000 MAGs. Downloading SPIRE's 10,243 MAGs would be the cheap route if those samples are ever
  wanted.

## The other gaps (in progress or explained)
- **11 imported libraries never assembled** (31.6–41.8 Gbp each). Their only MEGAHIT attempt (2026-07) hung (10) or
  errored (SRR13021855). Re-submitted 2026-11-15 with the identical MEGAHIT release and settings
  (`wetland/reassembly_11/`). Their MAGs will join the census with frozen iterativeIDs.
- **68 of the 770 assemblies yield no MAG.** Almost all are very low-yield assemblies from deep libraries (e.g. 73.9 Gbp
  of reads, 8 kb assembled). They are a re-assembly question, not a binning failure.
- **63 DAS Tool runs produced no bins;** 7 were retried successfully on 2026-11-03 (248 MAGs, now in the census).
