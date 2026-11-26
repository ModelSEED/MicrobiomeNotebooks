#!/usr/bin/env python
"""Assemble the machine-readable plan TSVs from the raw scans (scan/*.tsv). READ-ONLY w.r.t. the repo.
Outputs in rename_plan/data/:
  files_plan.tsv            one row per in-scope tracked file
  columns_plan.tsv          one row per (file, column / JSON key path / text part); community_modeling provenance
                            subtrees aggregated to one row per subtree
  taxonomy_plan.tsv         per (file, MAG taxonomy column): rank, cells, cells that change in R232, GTDB vs MiDAS
  code_fixes.tsv            hard-coded strings in code (taxon names, IDs, colour keys) with the fix
  binary_for_docs.tsv       xlsx/docx/pptx/pdf/png hits, handed to the documents part
  apply_mag_ids_extension.tsv  (file, delimiter, MAG column, ID column) to add to scripts/apply_mag_ids.py
"""
import collections, csv, json, os, re
csv.field_size_limit(1 << 30)
D = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan/data"
S = os.path.join(D, "scan")
def rd(n): return list(csv.DictReader(open(os.path.join(S, n)), delimiter="\t"))
_cw = csv.DictWriter.writerows
def _clean_rows(self, rows):
    return _cw(self, ({k: (re.sub(r"[\t\r\n]+", " ", v) if isinstance(v, str) else v) for k, v in r.items()} for r in rows))
csv.DictWriter.writerows = _clean_rows
cols = rd("columns_raw.tsv"); files = {r["file"]: r for r in rd("files_raw.tsv")}
gens = {r["file"]: r for r in rd("generators.tsv")}
bk = rd("bin_keyed_columns.tsv"); bkt = rd("bin_keyed_taxonomy.tsv"); taxj = rd("taxonomy_cells_raw.tsv")
names = rd("r220_names_raw.tsv"); code = rd("code_hardcoded_raw.tsv"); binr = rd("binary_raw.tsv")

# ---------------- curated knowledge ----------------------------------------------------------------------------
CURATED_GEN = {   # file -> (generator, how to refresh)
 "gtdbtk.bac120.summary.tsv": ("GTDB-Tk run (poplar processed/mag/gtdbtk_291 + best-of-both hybrid, reprocess_291/bestofboth/compose_static.py)", "replace by results/gtdbtk.r232.summary.tsv (byte copy)"),
 "mag_quality_summary.tsv": ("scripts/build_mag_quality_table.py (poplar /scratch1 paths) + best-of-both compose_static.py", "patch Domain..Species by bin (do NOT re-run the builder: regresses best-of-both)"),
 "mag_abundance_summary.tsv": ("reprocess_291/build_mag_abundance_291.py + bestofboth/compose_abundance_summary.py (poplar)", "patch Domain..Species by bin; mag_iterativeID via apply_mag_ids.py"),
 "mag_abundance_summary.json": ("reprocess_291/build_mag_abundance_291.py + postprocess_291.py (poplar)", "patch taxonomy{Domain..Genus} by bin; mag_iterativeID via apply_mag_ids.py"),
 "mag_iterativeID_old_to_new.csv": ("scripts/recompute_mag_ids.py", "re-run recompute_mag_ids.py after the taxonomy patch (keeps the ancient old columns); keep R220->R232 map separately"),
 "mag_iterativeID_old_to_new.json": ("scripts/recompute_mag_ids.py", "re-run recompute_mag_ids.py after the taxonomy patch"),
 "scripts/_new_ids.json": ("scripts/recompute_mag_ids.py (mag_new / mag_case)", "re-run recompute_mag_ids.py"),
 "mag_abundance_by_day_intersection_with_ids.csv": ("scripts/build_mag_intersection_with_ids.py; IDs by apply_mag_ids.py", "apply_mag_ids.py (MAG, iterativeID); its Kingdom..Species columns are MiDAS (ASV) taxonomy: leave"),
 "mag_percell_abundance.csv": ("scripts/build_percell_abundance.py (local inputs present: counts_291_mag_CAN_*.json, meta/ NanoStats, mag_metagenome_abundance.csv)", "re-run after mag_quality_summary.tsv + map are updated (taxonomy and iterativeID come from those)"),
 "Table_S7_MAG_functional_classification.csv": ("scripts/build_table_s7_classification.py (reads MAG_target_gene_filled.xlsx S1/S3)", "re-run after the xlsx is rebuilt; or apply_mag_ids (bin join) + Genus/Family/Phylum patch"),
 "mag_guild_classification.csv": ("scripts/build_table_s7_classification.py (GUILD_OUT)", "as Table S7 (same run)"),
 "ko_copy_number_matrix_nosZ.csv": ("derived from ko_copy_number_matrix.csv (reprocess_291/build_ko_matrices.py; nosZ split) ; IDs by apply_mag_ids.py", "apply_mag_ids.py (already listed)"),
 "clade_classify/out/allbins_nosz_master.tsv": ("clade_classify pipeline (partE_reconcile / subclade_all_bins); IDs by apply_mag_ids.py", "apply_mag_ids.py (already listed) + genus column patch by bin"),
 "gene_ab_figure/data/taxonomy_labels.tsv": ("gene_ab_figure/01_select_mags_and_genes.py era; IDs by apply_mag_ids.py", "apply_mag_ids.py (already listed) + Genus/Species patch by bin"),
 "gene_ab_figure/offline_bundle/data/taxonomy_labels.tsv": ("gene_ab_figure/05_bundle_offline.py (copy)", "apply_mag_ids.py (already listed) + Genus/Species patch"),
 "gene_ab_figure/data/gene_copy_per_mag.tsv": ("gene_ab_figure/08_gene_copy_per_mag.py", "apply_mag_ids.py (already listed) + genus patch; or re-run 08"),
 "gene_ab_figure/data/gene_copy_per_mag.md": ("gene_ab_figure/08_gene_copy_per_mag.py", "re-run 08 (writes .tsv and .md)"),
 "gene_ab_figure/data/rpkm_long_MAG_gene_sample_phase.tsv": ("scripts/build_rpkm_long_cube.py", "apply by bin join (mag, iterativeID) + genus patch; IDs here are partly stale"),
 "gene_ab_figure/data/rpkm_long_MAG_gene_sample_291_preview.tsv": ("scripts/build_rpkm_cube_291_preview.py", "apply by bin join (mag, iterativeID) + genus patch; preview file, IDs partly stale"),
 "gene_ab_figure/data/rpkm_long_MAG_gene_sample_phase_hq.tsv": ("reprocess_291/build_rpkm_cube_hqthauera.py", "genus patch by bin (no ID column)"),
 "ko_gene_abundance/ko_rpkm_combined.tsv": ("ko_gene_abundance/03_aggregate_ko_rpkm.py; IDs by apply_mag_ids.py", "apply_mag_ids.py (already listed) + genus patch"),
 "ko_gene_abundance/ko_rpkm_per_sample.tsv": ("ko_gene_abundance/03_aggregate_ko_rpkm.py; IDs by apply_mag_ids.py", "apply_mag_ids.py (already listed) + genus patch"),
 "ko_gene_abundance/ko_rpkm_combined_ext.tsv": ("ko_gene_abundance/07_aggregate_ext.py", "ADD to apply_mag_ids.py (MAG, mag_iterativeID); IDs partly stale; + genus patch"),
 "ko_gene_abundance/ko_rpkm_per_sample_ext.tsv": ("ko_gene_abundance/07_aggregate_ext.py", "ADD to apply_mag_ids.py (mag, iterativeID); + genus patch"),
 "ko_gene_abundance/ko_rpkm_root_combined.tsv": ("ko_gene_abundance/04_root_aggregate.py", "re-run 04 after the ko_rpkm tables (roots = MAG ID roots: Rhodobacter -> Pseudogemmobacter etc.)"),
 "ko_gene_abundance/ko_rpkm_root_per_sample.tsv": ("ko_gene_abundance/04_root_aggregate.py", "re-run 04"),
 "Table_S6_skani_novelty.tsv": ("reprocess_291/build_skani_novelty_291.py (poplar) / server_handoff_out/build_skani_ani_class.py", "regenerate mag_gtdb_taxonomy from the R232 summary; skani_best_taxonomy against R232 refs (needs the R232 skani run: skani_preview/ in the prep folder)"),
 "fig1_tree/fig1_annotation.tsv": ("scripts/render_fig1_tree.py --annotation-out (fig1_tree/run_fig1_raxml.sh)", "re-render Fig 1 (no RAxML re-run: tips are bins); MAG taxonomy + IDs from repo files; 61 reference tips still carry R220 names from the gtdbtk_r220_markers resources (21 change in R232, 1 left GTDB)"),
 "fig1_tree/fig1_mags_table.tsv": ("scripts/select_fig1_mags.py", "re-run select_fig1_mags.py (reads the map) or apply by bin join"),
 "network/louvain_module_mag_membership.json": ("scripts/build_module_mag_membership.py", "re-run after the networks (NETWORK_MAG_ONLY decision D2)"),
 "community_modeling/analysis/r13/annotations/taxonomy.tsv": ("community_modeling/analysis/r13/annotations/build_taxonomy.py (bin->iterativeID from kbase_mag_rename_log.json)", "re-run build_taxonomy.py after gtdbtk summary + FLAG/DIDERM fix + KBase rename log (or read the new map)"),
 "community_modeling/analysis/r13/annotations/panel/mag_panel.tsv": ("community_modeling/analysis/r13/annotations/panel/build_panel.py", "re-run build_panel.py (gtdb + kbase_mag_rename_log.json inputs)"),
 "community_modeling/analysis/r13/annotations/calls/summary.tsv": ("community_modeling/analysis/r13/annotations/calls/call_genes.py", "re-run only the summary step, or patch gtdb_phylum/gtdb_order by bin (gene calls do not depend on taxonomy: CODE25 is by bin)"),
 "community_modeling/analysis/r13/annotations/rast/genomes/_manifest.tsv": ("r13 RAST fetch (rast/fetch_rast.py family)", "KBase object-name record: follows the KBase part (genome_name embeds the iterativeID)"),
 "community_modeling/analysis/r13/annotations/kbase_kofam_push/placement_by_bin.tsv": ("kbase_kofam_push/build_payloads.py", "KBase object-name record: follows the KBase part"),
 "kbase_mag_rename_log.json": ("scripts/rename_mags_kbase.py", "KBase part: a new rename run writes the new log; consumers build_taxonomy.py / build_panel.py"),
 "mimag_hq_classification.tsv": ("reprocess_291/build_mimag_hq.py", "patch Phylum/Genus/Species by bin (or re-run if inputs present)"),
 "phase_rpkm_chord/phase_mag_gene_rpkm_long.csv": ("phase_rpkm_chord/01_build_long_table.py", "re-run 01..04 (IDs in this folder are stale pre-408cb39 labels)"),
 "phase_rpkm_chord/mag_category_assignments.csv": ("phase_rpkm_chord/01_build_long_table.py", "re-run 01"),
 "module_pathway_completeness/module_denitrification_genes_phaseI.csv": ("scripts/render_module_denitrification_genes.py", "re-run after networks (MAG link by ID string)"),
 "tablesforpaper/guild_flanking_all.csv": ("spieceasi_modules_flanking.py", "re-run (partner_is_MAG / guild_node_bin come from the ID-string MAG link)"),
}
for ph in ("III", "IV", "V"):
    CURATED_GEN[f"module_pathway_completeness/module_denitrification_genes_phase{ph}.csv"] = CURATED_GEN["module_pathway_completeness/module_denitrification_genes_phaseI.csv"]
for f in ("guild_flanking_nonMAG.csv", "modules_membership_89node.csv", "pao_gao_direct_edges.csv", "pao_gao_guild_series_env.csv", "pao_gao_shared_driver_summary.json"):
    CURATED_GEN["tablesforpaper/" + f] = CURATED_GEN["tablesforpaper/guild_flanking_all.csv"]

DOCS_PART = {"MANUSCRIPT_methods_results_discussion.md"}
ASV_FILES = re.compile(r"^(abundances( \(1\))?\.csv|table_rel_export\.csv|asv_iterativeID_old_to_new\.(csv|json)|iterativeID_[a-z_]+\.json|"
                       r"iterativeIDs\.json|taxonomy\.json|taxonIDs\.json|nonzero_per_day\.json|significantly_connected_organisms\.json|"
                       r"relative_abundance\.json|data_processing\.ipynb|scripts/data_processing\.py|__marimo__/|correlations/|"
                       r"network_method_comparison/|dbRDA/|leftover_reads_analysis/|abundace_correlations/|rel_ab_heatmaps/|qiime/|"
                       r"bvbrc_|synthetic_genomes|network/network_module_membership|network/cooccurrence|network/network_summary|"
                       r"FDR_passing)")
CM_PROV = re.compile(r"^community_modeling/(logs/|runs/|escher_maps/|plan_evidence/|plan_r13/|presentation/|"
                     r"analysis/(r12_final|r12_interim|data_inconsistencies)/|analysis/r13/(?!annotations/))")
CM_ANN = re.compile(r"^community_modeling/analysis/r13/annotations/")
NET_MAGLINK = re.compile(r"^(network/(compositional_networks/|louvain_module_mag_membership)|tablesforpaper/|module_pathway_completeness/|"
                         r"spieceasi_modules_flanking\.py)")
ASV_HDR = re.compile(r"iterativeIDs$|member|asv|node\b|guild_node|^source$|^target$|^name$|partner|^a$|^b$|rep_iterativeID|loop_member|carriers|importers|exporters", re.I)
STALE_OK = re.compile(r"\.bak\.|pre_bb|superseded|/old/|_old\.|^server_handoff_out/|^mag_skani_hits/|^mag_coverage_reanalysis/|test62_", re.I)

def area(f):
    if CM_PROV.match(f): return "community_modeling (r12/r13 run provenance: members = ASV synthetic genomes)"
    if CM_ANN.match(f): return "community_modeling/analysis/r13/annotations (MAG-keyed annotation tables)"
    if f.startswith("community_modeling/"): return "community_modeling (pipeline code/config; members = ASV)"
    if NET_MAGLINK.match(f): return "network MAG-link (ASV nodes linked to MAGs by ID string: NETWORK_MAG_ONLY)"
    if f.startswith("network/"): return "network (ASV co-occurrence)"
    if ASV_FILES.match(f): return "ASV data"
    return f.split("/")[0] if "/" in f else "(repo root)"

bk_idx = {}
for r in bk:
    bk_idx[(r["file"], r["id_col"])] = r

def cm_subtree(f):
    p = f.split("/")
    return "/".join(p[:4]) + "/" if len(p) > 4 else "/".join(p[:-1]) + "/"

def final_ctx(r):
    f, col = r["file"], r["column"]
    hdr = r["header"] or col
    m, i = int(r["old_minted"] or 0), int(r["old_inherit"] or 0) + int(r["old_display"] or 0)
    a = area(f)
    raw = r["inherit_context"]
    # count split of inherit occurrences
    split = collections.Counter()
    if ":" in raw:
        for kv in raw.split(";"):
            k, v = kv.split(":"); split[k] += int(v)
    elif raw:
        split[raw] += i
    basis = r["inherit_context_basis"]
    colname = col.split("=", 1)[1] if "=" in col else col
    bkr = bk_idx.get((f, colname))
    if bkr:
        cur = int(bkr["current_exact"]) + int(bkr["current_embedded"]); n = int(bkr["rows_with_bin"])
        if ASV_HDR.search(colname):
            split = collections.Counter({"ASV": i}); basis = f"bin-keyed table but '{colname}' is an ASV/member column (equal strings only where the MAG inherited the ASV ID)"
        elif cur >= 0.5 * n:
            split = collections.Counter({"MAG": i}); basis = f"bin-keyed MAG ID column: {cur}/{n} rows equal the row bin's current ID ({bkr['stale_or_other_id']} stale)"
    if re.search(r"/kbase_kofam_push/|/rast/genomes/|/rast/fetch_rast|plan_r13/kbase_rast/", f):
        split = collections.Counter({"MAG": i}); basis = "KBase object names (emilykin_MAGs_RAST_code11_b<k>_<iterativeID>) keyed by bin"
        return m, i, split, basis
    if a.startswith("ASV data") or a.startswith("network (ASV"):
        split = collections.Counter({"ASV": i}); basis = "ASV-keyed data (area rule)"
    elif a.startswith("community_modeling (r12/r13") or a.startswith("community_modeling (pipeline"):
        if i and not split.get("MAG"):
            split = collections.Counter({"ASV": i}); basis = "community member names are ASV synthetic genomes (area rule)"
        elif i:
            basis += " | community_modeling: MAG only where a bin / KBase genome name is on the line"
    elif a.startswith("network MAG-link"):
        split = collections.Counter({"ASV": i}); basis = "ASV node IDs whose MAG link is the ID-string identity (NETWORK_MAG_ONLY): labels stay, MAG membership changes"
    return m, i, split, basis

# ---------------- columns_plan ----------------------------------------------------------------------------------
crow = []; cm_agg = collections.defaultdict(lambda: collections.Counter()); cm_ex = {}
file_tot = collections.defaultdict(lambda: collections.Counter())
for r in cols:
    f = r["file"]
    m, i, split, basis = final_ctx(r)
    na = int(r["new_asv_preexisting"] or 0)
    ft = file_tot[f]
    ft["minted"] += m; ft["inherit"] += i; ft["MAG"] += m + split.get("MAG", 0); ft["ASV"] += split.get("ASV", 0)
    ft["ambiguous"] += split.get("ambiguous", 0); ft["new_asv"] += na
    if m + i == 0 and na == 0: continue
    a = area(f)
    if a.startswith("community_modeling (r12/r13"):
        k = cm_subtree(f)
        c = cm_agg[k]; c["files"] += 0; c["minted"] += m; c["inherit"] += i; c["new_asv"] += na
        for kk, v in split.items(): c["inh_" + kk] += v
        c.setdefault("_files", set()).add(f)
        cm_ex.setdefault(k, (f, col := r["column"], r["example"]))
        continue
    if m + i == 0:
        ctx = "ASV (new-ID collision only)"
    else:
        ctx = "MAG" if (m and not (split.get("ASV") or split.get("ambiguous"))) or (split.get("MAG") and not (split.get("ASV") or split.get("ambiguous"))) \
            else "ASV" if not m and split.get("ASV") and not split.get("MAG") and not split.get("ambiguous") \
            else "mixed/ambiguous"
    g = CURATED_GEN.get(f, (gens.get(f, {}).get("writers", ""), ""))[0]
    crow.append(dict(file=f, area=a, column_or_key=r["column"], context=ctx, occ_old_minted=m, occ_old_inherit=i,
                     occ_inherit_MAG=split.get("MAG", 0), occ_inherit_ASV=split.get("ASV", 0),
                     occ_inherit_ambiguous=split.get("ambiguous", 0), occ_new_id_already_ASV=na,
                     old_ids_found=r["old_ids_found"][:200], basis=basis[:160], generated_by=g[:200], example=r["example"][:120]))
for k, c in sorted(cm_agg.items()):
    f, col, ex = cm_ex[k]
    crow.append(dict(file=k + f"  [{len(c['_files'])} files]", area=area(f), column_or_key="(aggregated; e.g. " + f.split("/")[-1] + " :: " + col[:40] + ")",
                     context="ASV (community member names)" if not c["minted"] else "ASV members + MAG KBase names",
                     occ_old_minted=c["minted"], occ_old_inherit=c["inherit"], occ_inherit_MAG=c["inh_MAG"],
                     occ_inherit_ASV=c["inh_ASV"], occ_inherit_ambiguous=c["inh_ambiguous"], occ_new_id_already_ASV=c["new_asv"],
                     old_ids_found="", basis="r12/r13 run provenance; members = ASV synthetic genomes (README line 4)",
                     generated_by="community_modeling pipeline (ek_pipeline.py / analysis scripts)", example=ex[:120]))

def col_action(r):
    f = r["file"]; a = r["area"]; ctx = r["context"]
    if a.startswith("community_modeling (r12/r13"): return "leave (frozen run provenance; members are ASVs; no re-run needed)"
    if STALE_OK.search(f): return "leave (backup / superseded copy)"
    if ctx.startswith("ASV"):
        if a.startswith("network MAG-link"): return "regenerate per decision D2 (labels are ASV IDs and stay; MAG membership changes)"
        return "leave (ASV context)"
    if f in CURATED_GEN: return CURATED_GEN[f][1]
    if a.startswith("community_modeling/analysis/r13/annotations"):
        cn = r["column_or_key"]
        if f.endswith((".log", ".out", ".txt")): return "leave (run log / provenance record)"
        if f.endswith(".md"): return "edit text: MAG sentences only (manual; most IDs here are ASV member names)"
        if f.endswith(".py"): return "fix code strings (see code_fixes.tsv)"
        if re.search(r"genome_name|feature|output_info|output_name|/name|kbase", cn, re.I): return "follow KBase part (object names embed the ID); leave as provenance record"
        return "patch by bin join or re-run its builder"
    if f.endswith((".md", ".txt", ".html")): return "edit text: MAG sentences only (manual)"
    if f.endswith((".py", ".R", ".sh", ".ipynb")): return "fix code strings (see code_fixes.tsv)"
    if ctx == "MAG": return "apply_mag_ids (bin join) / map by simultaneous dict"
    return "manual review"
for r in crow: r["action"] = col_action(r)
with open(os.path.join(D, "columns_plan.tsv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(crow[0]), delimiter="\t", lineterminator="\n"); w.writeheader()
    w.writerows(sorted(crow, key=lambda r: (r["area"], r["file"])))

# ---------------- taxonomy_plan ----------------------------------------------------------------------------------
trow = []
MIDAS_HINT = re.compile(r"^(Kingdom|tax_|asv_|root$)", re.I)
for r in bkt:
    n, nch, nb = int(r["n_cells_r220"]), int(r["n_cells_change_r232"]), int(r["nonblank_cells"])
    frac = float(r["gtdb_match_fraction"])
    f = r["file"]
    kind = "GTDB (MAG)" if frac >= 0.85 or r["rank"] == "lineage" else "MiDAS/ASV or mixed (coincidental matches)"
    if MIDAS_HINT.search(r["tax_col"]) or f.startswith("mag_abundance_by_day"): kind = "MiDAS/ASV (Kingdom..Species of the mapped ASV)"
    if re.search(r"ref_|best_ref|skani_best|closest", r["tax_col"]): kind = "reference genome taxonomy (GTDB release of the reference DB)"
    if n < 2 and nch == 0: continue
    a = area(f)
    if kind.startswith("MiDAS"): act = "leave (ASV taxonomy)"
    elif STALE_OK.search(f) or f.startswith("server_handoff_out/") or f.startswith("mag_skani_hits/"): act = "leave (historical / superseded run record)"
    elif a.startswith("community_modeling (r12/r13"): act = "leave (frozen run provenance)"
    elif kind.startswith("reference"): act = "regenerate against R232 references (needs the R232 skani/GTDB lookup), else leave with a note"
    elif f in CURATED_GEN: act = CURATED_GEN[f][1]
    else: act = "patch by bin join from mag_lineage_r220_r232.tsv"
    trow.append(dict(file=f, area=a, key="bin column: " + r["bin_col"], tax_column=r["tax_col"], rank=r["rank"], kind=kind,
                     cells_r220_match=n, nonblank=nb, match_fraction=frac, cells_change_r232=nch, mags_change=r["n_mags_change"],
                     generated_by=CURATED_GEN.get(f, (gens.get(f, {}).get("writers", ""), ""))[0][:200], action=act))
for r in taxj:     # JSON taxonomy cells (path-level)
    n, nch = int(r["n_mag_cells"]), int(r["n_need_change"])
    f = r["file"]
    if not f.endswith(".json") or n < 2: continue
    a = area(f)
    act = CURATED_GEN.get(f, ("", "patch by bin join"))[1]
    if a.startswith("community_modeling (r12/r13"): act = "leave (frozen run provenance)"
    if f.startswith("bvbrc_") or "skani" in f: act = "leave (reference-hit record)"
    trow.append(dict(file=f, area=a, key="JSON path", tax_column=r["column"], rank=r["rank"], kind="GTDB (MAG) / JSON",
                     cells_r220_match=n, nonblank="", match_fraction="", cells_change_r232=nch, mags_change="",
                     generated_by=CURATED_GEN.get(f, (gens.get(f, {}).get("writers", ""), ""))[0][:200], action=act))
# keep one row per (file, column): the majority rank; minority-rank rows are same-name placeholders across ranks
# (e.g. Order = Family = 'UBA5066'), whose 'change' would be attributed to the wrong rank
best = {}
for r in trow:
    k = (r["file"], r["tax_column"])
    if k not in best or int(r["cells_r220_match"]) > int(best[k]["cells_r220_match"]): best[k] = r
for r in trow:
    k = (r["file"], r["tax_column"])
    if best[k] is not r:
        best[k]["cells_r220_match"] = int(best[k]["cells_r220_match"]) + int(r["cells_r220_match"])
        best[k]["note"] = (best[k].get("note", "") + f" +{r['cells_r220_match']} cells equal the MAG's {r['rank']} name (same placeholder at two ranks; not counted as change)").strip()
trow = list(best.values())
for r in trow:
    r.setdefault("note", "")
    if r["kind"].startswith("MiDAS/ASV or mixed") and r["nonblank"] != "":
        frac = int(r["cells_r220_match"]) / max(1, int(r["nonblank"]))
        r["match_fraction"] = round(frac, 3)
        if frac >= 0.7 and int(r["cells_r220_match"]) >= 5:
            r["kind"] = "GTDB (MAG)"
            f = r["file"]
            r["action"] = (CURATED_GEN[f][1] if f in CURATED_GEN else
                           "leave (historical / superseded run record)" if (STALE_OK.search(f) or f.startswith(("server_handoff_out/", "mag_skani_hits/"))) else
                           "leave (frozen run provenance)" if area(f).startswith("community_modeling (r12/r13") else
                           "patch by bin join from mag_lineage_r220_r232.tsv (genus column; keep the file's own suffix convention)")
with open(os.path.join(D, "taxonomy_plan.tsv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(trow[0]), delimiter="\t", lineterminator="\n"); w.writeheader()
    w.writerows(sorted(trow, key=lambda r: (r["action"].startswith("leave"), r["file"])))

# ---------------- code_fixes ----------------------------------------------------------------------------------------
FIX = {
 ("scripts/render_fig1_tree.py", "Patescibacteria"): "root on r['Phylum'] in {'Patescibacteria','Patescibacteriota'}; add GTDB2PALETTE 'Patescibacteriota':'Patescibacteria', 'Myxococcota_A':'Myxococcota' (palette keys are MiDAS names)",
 ("fig1_tree/compare_trees.py", "Patescibacteria"): "same rooting test as render_fig1_tree.py (accept Patescibacteriota)",
 ("community_modeling/analysis/r13/annotations/build_taxonomy.py", "Patescibacteria"): "FLAG key 'Patescibacteria' -> also 'Patescibacteriota' (else CPR flag silently drops for 21 MAGs)",
 ("community_modeling/analysis/r13/annotations/build_taxonomy.py", "Myxococcota"): "DIDERM add 'Myxococcota_A' (coasm_bin.258) or strip GTDB _X suffixes before lookup",
 ("community_modeling/analysis/r13/annotations/panel/tests/test_build_panel.py", "Patescibacteria"): "test lineage string: update to p__Patescibacteriota (fixture only)",
 ("reprocess_291/run_tree62.sh", "Patescibacteria"): "historical server script (tree62, GTDB-Tk de_novo outgroup p__Patescibacteria); leave or update if ever re-run on R232",
 ("scripts/render_percell_composition.py", "Giesbergeria"): "DOCX_COLORS key Giesbergeria -> Simplicispira (2 MAGs renamed)",
 ("scripts/render_percell_composition.py", "SCMT01"): "DOCX_COLORS key SCMT01 -> Lysobacter (genus_key strips _I of Lysobacter_I)",
 ("scripts/render_percell_composition.py", "JJ008"): "JJ008 keeps 1 of 2 MAGs (other -> CADEDV01): check the top-20 band set after the patch",
 ("scripts/render_percell_composition.py", "Azospira_A"): "Azospira_A -> Proximibacter (2 MAGs): 'Azospira' key no longer covers them; add Proximibacter colour",
 ("scripts/build_fig1_alignment.py", "Azospira_A"): "comment only (reference choice for CAN_4_bin.208); now g__Proximibacter",
 ("fig1_tree/score_ref_candidates.py", "Azospira_A"): "docstring/default genus; now g__Proximibacter",
 ("fig1_tree/run_fig1_raxml.sh", "Azospira_A"): "comment only",
 ("community_modeling/analysis/r13/annotations/calls/call_genes.py", "Patescibacteria"): "docstring only (GTDB wording); logic is CODE25 by bin, unaffected -> optional text update",
 ("community_modeling/analysis/r13/annotations/panel/validate_panel.py", "Patescibacteria"): "message text only (GTDB wording); CODE25 by bin -> optional text update",
 ("network/compositional_networks/render_spieceasi_docx.py", "Patescibacteria"): "leave: module phylum mix of ASV nodes (MiDAS names); numbers/node lists in this script change only under D2a",
 ("network/compositional_networks/render_spieceasi_docx.py", "Myxococcota"): "leave: MiDAS phylum of ASV nodes",
 ("community_modeling/member_reference_phyla.py", "Patescibacteria"): "leave: NCBI/MiDAS phylum names of synthetic-genome references (ASV side)",
 ("community_modeling/ek_pipeline.py", "Patescibacteria"): "leave: NCBI names (ASV synthetic genomes)",
 ("clade_classify/partA_reclassify.py", "Myxococcota"): "comment only; leave",
 ("community_modeling/analysis/r13/annotations/verify_design/verify_design.py", "Ignavibacteria.31"): "leave: Ignavibacteria.31 is the ASV member; CAN_3_bin.221 is its MAG (whose own ID Ignavibacteria.21 is not on this line)",
 ("community_modeling/analysis/r13/annotations/verify_design/v7_storage_gpr_reach.py", "Ignavibacteria.31"): "leave: Ignavibacteria.31 is the ASV member; CAN_3_bin.221 is its MAG",
 ("community_modeling/analysis/r13/annotations/design_answers/compute_answers.py", "Ignavibacteria.31"): "leave: r13 member name (ASV)",
 ("community_modeling/analysis/r13/annotations/verify_design/v8_a52_shares.py", "Ignavibacteria.31"): "leave: r13 member name (ASV)",
 ("community_modeling/curate_templates.py", "Ignavibacteria.13"): "leave: community member name (ASV)",
}
krow = []
for r in code:
    f = r["file"]; a = area(f)
    key = (f, r["string"])
    if key in FIX:
        act = FIX[key]
    elif a.startswith("community_modeling (r12/r13") and r["kind"] != "old_minted":
        act = "leave (provenance script; members are ASVs)"
    elif r["kind"].startswith("r220_taxon_name") and r["context"] == "ASV-name-collision" and key not in FIX:
        act = "check: MiDAS/NCBI (ASV) phylum name in this script -> leave unless it reads GTDB names"
    elif key in FIX: act = FIX[key]
    elif r["kind"] in ("old_inherit", "old_display") and r["context"] in ("ASV", "ambiguous") and ("community_modeling" in f or "spieceasi" in f):
        act = "leave: ASV member/node name (string stays); text describing it as a MAG node changes with D2"
    elif r["kind"] == "new_asv": act = "leave (ASV today); after rename the string is also a MAG ID"
    elif r["kind"] in ("old_inherit", "old_minted"): act = "manual: comment/docstring naming the MAG -> new ID"
    else: act = "manual review"
    krow.append(dict(file=f, line=r["line"], kind=r["kind"], string=r["string"], context=r["context"], area=a, action=act, snippet=r["snippet"][:180]))
with open(os.path.join(D, "code_fixes.tsv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(krow[0]), delimiter="\t", lineterminator="\n"); w.writeheader()
    w.writerows(sorted(krow, key=lambda r: (r["action"].startswith("leave"), r["file"], r["line"])))

# ---------------- binary_for_docs ---------------------------------------------------------------------------------
brow = []
for r in binr:
    hit = int(r["old_minted"] or 0) + int(r["old_inherit"] or 0) + int(r["old_display"] or 0)
    f = r["file"]
    g = gens.get(f, {})
    if r["ext"] in (".png", ".pdf", ".svg") or hit or r["r220_names"]:
        brow.append(dict(file=f, ext=r["ext"], area=area(f), old_minted=r["old_minted"], old_inherit=r["old_inherit"],
                         old_display=r["old_display"], bins=r["bins"], r220_names=r["r220_names"][:200],
                         writers_heuristic=g.get("writers", ""), readers_heuristic=g.get("readers", "")[:150],
                         note=("backup copy" if STALE_OK.search(f) else r["note"])))
with open(os.path.join(D, "binary_for_docs.tsv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(brow[0]), delimiter="\t", lineterminator="\n"); w.writeheader()
    w.writerows(sorted(brow, key=lambda r: (r["ext"], r["file"])))

# ---------------- apply_mag_ids extension ----------------------------------------------------------------------
listed = {("mag_abundance_summary.tsv", "mag_iterativeID"), ("mag_abundance_by_day_intersection_with_ids.csv", "iterativeID"),
          ("gene_ab_figure/data/taxonomy_labels.tsv", "mag_iterativeID"), ("gene_ab_figure/offline_bundle/data/taxonomy_labels.tsv", "mag_iterativeID"),
          ("gene_ab_figure/data/gene_copy_per_mag.tsv", "iterativeID"), ("ko_gene_abundance/ko_rpkm_combined.tsv", "mag_iterativeID"),
          ("ko_gene_abundance/ko_rpkm_per_sample.tsv", "iterativeID"), ("ko_copy_number_matrix_nosZ.csv", "iterativeID"),
          ("clade_classify/out/allbins_nosz_master.tsv", "iterativeID")}
arow = []
for r in bk:
    f, c = r["file"], r["id_col"]
    a = area(f)
    cur = int(r["current_exact"]) + int(r["current_embedded"]); n = int(r["rows_with_bin"])
    if ASV_HDR.search(c) or a.startswith("community_modeling (r12/r13") or a.startswith("network MAG-link"): continue
    if c in ("old_mag_iterativeID", "old_intersection_id", "closest_ref_tip", "candidate", "MAG_j", "mag_id", "replaces_old"): continue
    if cur == 0 and not re.search(r"iter|ID", c): continue
    embedded = int(r["current_embedded"]) > int(r["current_exact"])
    status = "already in apply_mag_ids.py" if (f, c) in listed else ("ADD" if not STALE_OK.search(f) and not f.startswith(("server_handoff_out/", "mag_coverage_reanalysis/")) else "skip (backup/historical)")
    if a.startswith("community_modeling/analysis/r13/annotations") and embedded: status = "follow KBase part (embedded object name)"
    if f in ("mag_iterativeID_old_to_new.csv",): status = "written by recompute_mag_ids.py"
    if f.startswith("fig1_tree/superseded") or "test62" in f: status = "skip (superseded / test)"
    arow.append(dict(file=f, delimiter=r["delim"], mag_col=r["bin_col"], id_col=c, rows_with_bin=n, current_ids=cur,
                     stale_ids=r["stale_or_other_id"], blank=r["blank"], changed_mags_in_file=r["changed_mags_in_file"],
                     value_form="embedded in a longer string" if embedded else "whole cell", status=status))
# JSON MAG-keyed
arow.append(dict(file="mag_abundance_summary.json", delimiter="json", mag_col="<top-level key = bin>", id_col="mag_iterativeID",
                 rows_with_bin=291, current_ids=291, stale_ids=0, blank=0, changed_mags_in_file=51, value_form="whole value", status="already in apply_mag_ids.py"))
arow.append(dict(file="kbase_mag_rename_log.json", delimiter="json", mag_col="<top-level key = bin>", id_col="new_name",
                 rows_with_bin=291, current_ids=291, stale_ids=0, blank=0, changed_mags_in_file=51, value_form="whole value",
                 status="KBase part (rename_mags_kbase.py rewrites it); r13 build_taxonomy.py / build_panel.py read it"))
with open(os.path.join(D, "apply_mag_ids_extension.tsv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(arow[0]), delimiter="\t", lineterminator="\n"); w.writeheader()
    w.writerows(sorted(arow, key=lambda r: (r["status"], r["file"])))

# ---------------- files_plan ----------------------------------------------------------------------------------
tax_by_file = collections.Counter(); tax_kind = collections.defaultdict(set)
for r in trow:
    if not r["action"].startswith("leave"):
        tax_by_file[r["file"]] += int(r["cells_change_r232"] or 0); tax_kind[r["file"]].add(r["tax_column"])
names_by_file = collections.defaultdict(list)
for r in names:
    names_by_file[r["file"]].append(f"{r['r220_name']}:{r['n']}")
def _is_leave(a): return a.startswith("leave") or a.startswith("comment only; leave") or a.startswith("check: MiDAS")
code_by_file = collections.Counter(r["file"] for r in krow if not _is_leave(r["action"]))
code_acts = collections.defaultdict(set)
for r in krow: code_acts[r["file"]].add(r["action"])
col_actions = collections.defaultdict(set)
for r in crow:
    col_actions[r["file"]].add(r["action"])
frow = []
inscope = set(file_tot) | set(tax_by_file) | set(code_by_file) | {r["file"] for r in names if r["file"] in files}
for f in sorted(inscope):
    t = file_tot.get(f, collections.Counter())
    if not (t["minted"] or t["inherit"] or t["new_asv"] or tax_by_file.get(f) or code_by_file.get(f) or names_by_file.get(f)): continue
    a = area(f)
    g = gens.get(f, {})
    cur = CURATED_GEN.get(f)
    shared = int(g.get("basename_shared_by") or 1) > 1
    gen_by = cur[0] if cur else ("" if shared else g.get("writers", ""))
    if cur: generated = "yes (curated)"
    elif shared: generated = f"unknown (basename shared by {g.get('basename_shared_by')} tracked files; not attributed)"
    elif gen_by: generated = "yes (heuristic: script names the file near a write call)"
    else: generated = "no writer found (source/manual, or pattern-named output)"
    if a.startswith("community_modeling (r12/r13"):
        act = "leave (frozen run provenance; members are ASVs)"
    elif STALE_OK.search(f):
        act = "leave (backup / superseded)"
    elif cur:
        act = cur[1]
    elif f in DOCS_PART:
        act = "documents part (manuscript text): MAG sentences only"
    elif re.match(r"^(SERVER_TASKS|SERVER_HANDOFF|paper3/REVIEW|reprocess_291/.*manuscript_stats|mag_skani_hits/REPORT|mag_coverage_reanalysis/)", f):
        act = "leave (historical record; D5)"
    else:
        acts = col_actions.get(f, set())
        nonleave = sorted(x for x in acts if not x.startswith("leave"))
        if nonleave: act = "; ".join(nonleave)
        elif tax_by_file.get(f): act = "patch taxonomy by bin join"
        elif code_by_file.get(f): act = "fix code strings (code_fixes.tsv)"
        elif acts: act = sorted(acts)[0]
        elif names_by_file.get(f): act = "check R220 taxon names in text (MAG sentences only)"
        else: act = "leave"
    if f in code_acts and f.endswith((".py", ".R", ".sh", ".ipynb", ".js")) and not a.startswith("community_modeling (r12/r13"):
        nl = sorted(x for x in code_acts[f] if not _is_leave(x))
        act = ("fix code strings: " + " | ".join(nl)) if nl else "leave (code strings are ASV member/node names or comments; see code_fixes.tsv)"
        if a.startswith("network MAG-link"):
            act = "D2a: update the hard-coded MAG-node counts / edge lists and re-render; D2b: leave (strings are ASV node names)"
    ctxs = []
    if t["MAG"]: ctxs.append(f"MAG {t['MAG']}")
    if t["ASV"]: ctxs.append(f"ASV {t['ASV']}")
    if t["ambiguous"]: ctxs.append(f"ambiguous {t['ambiguous']}")
    frow.append(dict(file=f, area=a, ext=files.get(f, {}).get("ext", os.path.splitext(f)[1]), generated=generated,
                     generated_by=gen_by[:220], readers=g.get("readers", "")[:160],
                     occ_old_minted=t["minted"], occ_old_inherit=t["inherit"], occ_MAG=t["MAG"], occ_ASV=t["ASV"],
                     occ_ambiguous=t["ambiguous"], occ_new_id_already_ASV=t["new_asv"],
                     mag_taxonomy_cells_change=tax_by_file.get(f, 0), r220_name_hits=";".join(names_by_file.get(f, []))[:160],
                     code_strings_to_fix=code_by_file.get(f, 0), context_summary=", ".join(ctxs), action=act))
with open(os.path.join(D, "files_plan.tsv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(frow[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(frow)
print("columns_plan", len(crow), "taxonomy_plan", len(trow), "code_fixes", len(krow), "binary_for_docs", len(brow),
      "apply_ext", len(arow), "files_plan", len(frow))
