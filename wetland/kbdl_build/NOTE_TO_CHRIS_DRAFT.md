Draft note to Chris Henry: Weisenhorn wetland MAG integration (from Andrew; edit before sending)

Hi Chris. I've taken over the wetland MAG integration from your 2026-11-11 handoff. One request is blocking, and there are
a few questions.

## 1. Blocking: please register GTDB R232 as a SKANIDB with its companion
I want the KBDL builds of all 12,859 wetland MAGs to carry **GTDB R232** taxonomy, not r220: the census is on R232, and
r220 assigns a species to only 3,141 of the MAGs against 6,497 under R232. KBDL has no R232 SKANIDB, and registering a
reference and attaching a companion are operator-only, so this needs you. I've prepared everything on poplar:
- **Sketch DB:** `/scratch1/afreiburger/kbdl_refs/gtdb_r232_skani/`. This is the GTDB-Tk R232 package's skani database
  (index.db, markers.bin, sketches.db, metadata.tsv: skani 0.3.1, consolidated, 199,923 species representatives),
  world-readable. Move it under your own tree first if you prefer.
- **Companion:** `/scratch1/afreiburger/kbdl_refs/gtdb_r232_skani_metadata.json`. 199,923 entries in the same format as
  the r220 companion (951d52…): key `<accession>_genomic.fna.gz` = companion_key(skani Ref_file); value
  {genus (unprefixed), species ("s__…"), lineage, ncbi_taxid: null}, built from GTDB R232 gtdb_taxonomy.tsv with the same
  rules as admin._parse_gtdb_lineage. `build-skani-metadata` can't be used because it needs the separate (*.sketch) layout.

Commands (as chenry on poplar):
```
cd /home/chenry/projects/KBDLJobRunningPrototype
/home/chenry/venvs/kbdl/bin/python -m kbdl_service.admin register-reference --type SKANIDB \
    --path /scratch1/afreiburger/kbdl_refs/gtdb_r232_skani --name gtdb_r232_consolidated \
    --visibility public --owner chenry --skani-version 0.3.1 --expect-genomes 199923
# -> prints the new object_id ref:…; put it in DB below
/home/chenry/venvs/kbdl/bin/python - <<'PY'
import io
from kbdl_service.config import load_settings
from kbdl_service.object_store import references
from kbdl_service.object_store.store import ObjectStore
from kbdl_service.schemas.objects import ObjectType
s = load_settings()
DB = "ref:…"            # from register-reference
payload = open("/scratch1/afreiburger/kbdl_refs/gtdb_r232_skani_metadata.json", "rb").read()
oid = ObjectStore(s).put(io.BytesIO(payload), ObjectType.SKANI_METADATA, "gtdb_r232_consolidated skani metadata", "public", "chenry")
assert references.set_metadata_object_id(s, DB, oid)
print("skani db", DB, "companion", oid)
PY
```
I tested both files with your skani 0.3.1 binary (/scratch1/chenry/opt/skani/0.3.1/skani): companion_key(Ref_file) of
the top hit resolves in the companion, and the lineage equals our GTDB-Tk 2.8.0 R232 call. For example, the canary
archaeon bin.013.fasta__SRR23247753 is Methanoflorens stordalenmirensis in R232, while r220 registered it as
"Bog-38 sp003139855" (bsp1).

Then send me the two ids. My build driver is paused until then.

Also: the 2-MAG canary and the 23-MAG pilot were built on 2026-11-15 against gtdb_r220_consolidated (genome ids vpa1, bsp1, …),
before I switched to R232, so their registry taxonomy is r220. Could you re-classify those 25 registry rows to R232, or tell
me the supported way to do it? Content-hash dedup means a rebuild will not change them.

## 2. Heads-up: a long build run on KBDL (after R232 is registered)
- I'm building all 12,859 MAGs with KBDLBuildGenome:
  - gene_caller="bakta", call_non_coding=true, bakta_db 6.0-full, skani_db = the R232 SKANIDB;
  - private GenomeArchives of 25 MAGs each.
- **The build is running, and since 2026-11-18 it uses all 4 KBDLBuildGenome slots.** At 3 slots it built about
  42 MAGs an hour, so the remaining ~11,900 MAGs need about 9 days at 4.
- **Request: could you raise the KBDLBuildGenome cap?** It is the only real limit on this run. The settings are
  `KBDL_JOB_CONCURRENCY=KBDLBuildGenome=16` (or whatever you're comfortable with) on the workers, which needs a
  worker restart.
  - **Footprint, measured 2026-11-18:** each build runs Bakta as about 4 single-CPU hmmsearch processes, with under
    1 GB RSS per container.
  - **Headroom:** poplar showed about 30–37% CPU idle (100+ cores) and about 1.8 TB of RAM available.
  - **Effect:** at 16 slots the remaining build drops from about 9 days to 2–3. Raising `KBDL_BAKTA_CALL_THREADS`
    instead helps less, because the hmmsearch steps run single-threaded.
- **Searches:** I'm also running up to 16 of the 40 KBDLSKANI slots for a few hours. Every wetland archive is now
  uploaded ahead of the build, so the R232 and novel_all_representatives searches can finish now instead of
  trailing the build for two weeks.
- My driver has been through an adversarial review and failure simulations. It never cancels or blindly resubmits, and
  it verifies every genome (content hash, Bakta provenance, CDS, called RNA).
- Things on your side:
  1. Bakta tempdirs go to ~chenry/.kbutillib/bakta_work on /mnt/homes, which is ~96% full (about 17 GB free).
  2. If a Bakta container hangs, only you can `docker kill` it. May I ping you when my watchdog flags one?
  3. Harvested job results/workdirs are about 100 GB in all in /scratch1/chenry/kbdl/jobs. OK to DELETE /jobs/{id} after
     I've verified a result?
- The archives are private. KBDL has no per-user grant, so I can't share them with you through the API; as operator you
  can read them. Say if you'd rather they were public.
- **Bug report (2026-11-17), KBDLBuildGenome + Bakta:**
  - **What happened:** `bin.001.fasta__SRR11599182` (6.0 Mb, 971 contigs; archive wetland_mags_A0008, job
    93f3c7e0…) failed as a member `tool_failure`:
    `Bakta: bakta 1.11.4 returned an out-of-range interval on contig 'pending.c779' (length 6085): begin=5721 end=6086`.
    The feature ends 1 bp past the contig end, and `_check_called_interval` in genome_normalize.py rejects the whole
    genome by design.
  - **A second case** has the same shape: `bin.001.fasta__SRR11847932`, contig 'pending.c182' (length 11451),
    begin=10879 end=11452. Both times the end is the contig length + 1. That's about 1% of MAGs so far, so expect
    about 100 of the 12,859 to fail this way.
  - **Why it matters:** the failure is deterministic, so a retry won't help.
  - **Suggestion:** clamp a feature that runs over the contig end by 1 bp and mark it partial, or drop that one
    feature and record it in `cds_discarded`. Either way the genome builds.
  - **On my side:** my driver now records these failures (kv `bakta_edge_mags`) and keeps going, unless there are more
    than 2 in a job or they exceed 3% of harvested MAGs. I'll rebuild the listed MAGs once there's a fix.

## 3. Landing in kbaseincubator.genome_clearhouse.skani_distances (planned as a separate final step; my account gets 403)
- C1. Is this table still the target?
- C2. source_catalog: "WEISENHORN" or "Weisenhorn Wetland Study"?
- C3. The dedup key: (source_catalog, genome_id, db_object_id)? Does a MAG with no hits get a row? Is top_ani a fraction
  or a percent? (land.py / SCRIPTS.md would answer all of this.)
  - **Hits come from KBDLSKANI.** I'm running KBDLSKANI over every wetland archive against both
    gtdb_r232_consolidated and novel_all_representatives, the same pattern as your R232 load.
  - **But KBDLSKANI drops the alignment fraction.** Its hits carry only distance = 1 − ANI. For example,
    bin.009.fasta__SRR16547334 has a 99.12% R232 hit at AF 0.484, which GTDB-Tk rejects as a species match.
  - **Heads-up: KBDLSKANI results depend on how many genomes are in one job.**
    - KBUtilLib passes every query of a job to a single `skani search` call. With more than 100 query files, skani
      0.3.1 switches to its marker-hash index (`--no-marker-index` says "default: load index if > 100 query
      files"), and that index loses hits.
    - Measured on novel_all_representatives, with the same 50 MAGs: 131 hits brute-force vs 124 indexed (101
      queries).
    - At 101 queries, all 7 lost hits came from one 31.6 kb MAG. At 125 queries, 375 normal-size MAGs lost nothing
      against either R232 or novel_all_representatives: their pairs were identical to brute force.
    - Larger calls lose more, including normal-size genomes:
      - at 2,000 queries, 2,476 rows (top 5 per MAG) vs 2,372 for the same MAGs inside one 12,859-query call;
      - in that 12,859-query call, a 3.2 Mb MAG lost its only hits (90.75%, AF 42%);
      - re-running our census searches brute-force gave 68 of 12,859 MAGs a novel-DB hit they had lost, all at
        82–91% ANI. No species-level hit changed.
    - Your R232 load batches 5 archives per job. If a job holds several hundred genomes or more, some hits are likely
      missing, mostly for small genomes and at lower ANI. Passing `--no-marker-index`, or capping each call at
      100 queries, avoids it at little cost: my 100-MAG jobs take the same ~10 min as 125-MAG ones.
    - I cap my jobs at 100 MAGs.
  - **Question:** should skani_distances carry AF? If so, could KBDLSKANI return Align_fraction_ref/query (skani
    already reports them)? Until then I'll treat KBDL's ANI as a distance only, and keep species calls on our census
    rule (ANI ≥ 95 and AF ≥ 0.5).
- C4. genome_id: I plan to use the MAG iterativeIDs (e.g. `Nitrosomonas.1_m`, the GTDB-R232-based scheme from our other
  projects), with a mapping to the KBase object names.
- C5. Will you grant me kbaseincubator read-write, or would you rather land a parquet yourself? I'd also like to land
  CheckM2 into genome_quality (rerun with --allmodels for both completeness models) if you agree on the columns.
- May I delete your private weisenhorn_000000 archive (a9b06ab1…)?

Thanks, Andrew
