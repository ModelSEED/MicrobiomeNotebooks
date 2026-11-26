"""Drop-in additions for salterns/reassembly_2026/pipeline.py: per-site 'merge' and 'coassemble' stages.

Assumes a Trimmomatic stage named TRIM_STAGE whose ledger output_ref is the '<output_reads_name>_paired'
KBaseFile.PairedEndLibrary (kb_trimmomatic names it output_reads_name + '_paired', Impl:1198), and that the merge
and coassemble tasks use the site (R1, R2, R2A, SF2) as their ledger 'sample' key, e.g. 'merge:R1'.
Only the four sediment sites with six A/B/C x D1/D2 cores are co-assembled (R2_restored_* and water are not).
"""

TRIM_STAGE = "trim"  # rename to the actual Trimmomatic stage key
COASSEMBLY_SITES = ("R1", "R2", "R2A", "SF2")
CORES = [f"{c}_{d}" for c in "ABC" for d in ("D1", "D2")]  # A_D1, A_D2, B_D1, B_D2, C_D1, C_D2 (input_refs[0] donates metadata)

APPS_ADD = {
    "merge": ("kb_ReadsUtilities/KButil_Merge_MultipleReadsLibs_to_OneLibrary", "release"),  # de8955d0a9..., njs
    "coassemble": ("MEGAHIT/run_megahit", "release"),  # e2e6d557da..., MEGAHIT 1.2.9, bigmem
}


def coassembly_names(site):
    return {"merge": f"{site}.sediment_cores6.Trimmomatic_paired.merged.reads",
            "coassemble": f"{site}.sediment_cores6.MEGAHIT.coassembly.Assembly"}


def params_for_coassembly(stage, site, ledger, output_of, ws_name):
    """Service-shape params (one-element list) and output object name, or None if inputs are not ready."""
    n = coassembly_names(site)
    if stage == "merge":
        refs = [output_of(ledger, TRIM_STAGE, f"{site}_{c}") for c in CORES]
        if not all(refs):
            return None
        assert len(set(refs)) == 6, f"{site}: duplicate Trimmomatic refs {refs}"  # the app does not de-duplicate
        return ([{"workspace_name": ws_name,          # string NAME: used as ReadsUtils wsname and KBaseReport workspace_name
                  "input_refs": refs,                 # six 'wsid/objid/ver' PairedEndLibrary refs; >=2 is only logged
                  "output_name": n["merge"],
                  "desc": f"{site} sediment cores A,B,C x D1,D2: 6 Trimmomatic-paired PE libraries merged for a "
                          f"MEGAHIT 1.2.9 co-assembly"}],  # REQUIRED in practice: read unconditionally at Impl:2225
                n["merge"])
    if stage == "coassemble":
        merged = output_of(ledger, "merge", site)
        if not merged:
            return None
        return [{"workspace_name": ws_name, "read_library_ref": merged,
                 "output_contigset_name": n["coassemble"], "min_contig_length": 1000}], n["coassemble"]
    raise ValueError(stage)


def verify_merge(ws, merged_ref, input_refs):
    """End-to-end pairing guard: read_count/total_bases of the merged object must equal the sums of the inputs."""
    infos = ws.get_object_info3({"objects": [{"ref": r} for r in [merged_ref] + list(input_refs)],
                                 "includeMetadata": 1})["infos"]
    meta = [i[10] or {} for i in infos]
    got_reads, got_bases = int(meta[0]["read_count"]), int(meta[0]["total_bases"])
    exp_reads = sum(int(m["read_count"]) for m in meta[1:])
    exp_bases = sum(int(m["total_bases"]) for m in meta[1:])
    ok = got_reads == exp_reads and got_reads % 2 == 0 and got_bases == exp_bases
    return ok, {"merged_reads": got_reads, "expected_reads": exp_reads, "merged_bases": got_bases,
                "expected_bases": exp_bases}
