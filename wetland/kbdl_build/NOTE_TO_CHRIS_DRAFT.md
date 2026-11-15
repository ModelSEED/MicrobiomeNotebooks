Draft note to Chris Henry: Weisenhorn wetland MAG integration (from Andrew; edit before sending)

Hi Chris. I've taken over the wetland MAG integration from your 2026-11-11 handoff. Status and a few questions:

Where things stand
- The MAG set is final: 12,859 MAGs from 702 of the 770 assemblies. All carry local GTDB-Tk 2.8.0 / R232 and CheckM in
  narrative 272312.
- The DAS Tool "stall" was the account's ~10-job concurrency cap, not dead jobs.
- KBDLBuildGenome already calls RNA genes, so the Phase-1 halt is moot.
- I've computed a census: 6,064 species clusters (1,893 R232, 4,171 new) and a 5,356-MAG representative worklist.
  - CheckM2 for all MAGs was run on poplar.
  - The novel side is labelled with novel_all_representatives rather than novel_hq_v2, which is isolates only.
  - Write-up: MicrobiomeNotebooks/wetland/census/ on branch clauded-insights-figure7.

Heads-up: a long build run on KBDL
- I'm building all 12,859 MAGs with KBDLBuildGenome:
  - gene_caller="bakta", call_non_coding=true;
  - skani_db = gtdb_r220_consolidated (ref:cee67d7f…), bakta_db = 6.0-full;
  - private GenomeArchives of 25 MAGs, members named "<MAG>.fna" (e.g. bin.001.fasta__SRR10127268.fna).
- I'm keeping at most 3 of the 4 build slots, so other builds can interleave. Expect roughly 8–14 days.
- A 2-MAG canary and a 23-MAG pilot run first.
- Three things on your side:
  1. Bakta tempdirs go to ~chenry/.kbutillib/bakta_work on /mnt/homes, which is ~96% full (about 17 GB free). If it fills,
     every member fails as "Bakta:".
  2. If a Bakta container hangs, only you (docker group) can `docker kill` it. My driver never cancels and resubmits,
     because the worker keeps building after a cancel. It will flag a job that runs past 3× its predicted time; may I ping
     you then?
  3. Each harvested job's result/workdir is about 200 MB in /scratch1/chenry/kbdl/jobs, about 100 GB in all. OK to
     DELETE /jobs/{id} after I've harvested and verified a result?

Landing in kbaseincubator.genome_clearhouse.skani_distances (not done; my account gets 403 on kbaseincubator)
- C1. Is this table still the target? Or should it be kbaseincubator.clearinghouse, or nowhere?
- C2. source_catalog: "WEISENHORN" (short token, like BVBRC / AllTheBacteria / SPIRE) or "Weisenhorn Wetland Study"?
  What label did the rolled-back wetland rows use?
- C3. The dedup key: (source_catalog, genome_id, db_object_id)? Does a MAG with no hits get a row (n_hits=0) or none?
  Is top_ani a fraction or a percent, and what format is distances? Sharing land.py / SCRIPTS.md would answer all of this.
- C4. genome_id form: the census MAG name ("bin.001.fasta__SRR10127268") or your pod's "SRR…bin.001" form?
- C5. Will you grant me kbaseincubator read-write, or would you rather land a staged parquet yourself?
- Also: may I delete your private weisenhorn_000000 archive (a9b06ab1…)? It has 100 members: 68 are md5-identical to
  current MAGs and 32 are not in the final set.
- DBs I'd land: gtdb_r220_consolidated + novel_all_representatives (not novel_hq_v2).

Thanks, Andrew
