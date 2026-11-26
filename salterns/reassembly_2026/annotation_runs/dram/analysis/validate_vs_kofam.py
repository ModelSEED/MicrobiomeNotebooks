"""Validate the DRAM 1.5 run (558 best-per-species MAGs) against the KOfam track (1,324 MAGs).

Run from annotation_runs/dram/ with ~/Documents/py_venv/bin/python -I analysis/validate_vs_kofam.py

Checks
 1. per-MAG gene counts: DRAM annotations vs KOfam Prodigal genes (genes_all.tsv.gz) and run_all.log
 2. gene-level coordinate join (MAG, contig, start, end, strand)
 3. per-gene best KO for the 68 marker KOs: DRAM ko_id vs KOfam best KO (pre-Sec-rescue per_mag tables)
 4. per-MAG marker counts: DRAM vs KOfam n_best (pre-rescue) and vs guilds_1324 n_<marker> (post-rescue)
 5. distill consistency: 'Salterns methane' sheet counts == counts of ko_id in annotations
 6. genome_stats taxonomy / completeness join sanity

Outputs in results/validation/.
"""
import gzip
import os
import re
import sys

import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # annotation_runs/dram
RUNS = os.path.dirname(HERE)
KOF = os.path.join(RUNS, "kofam")
RES = os.path.join(HERE, "results")
OUT = os.path.join(RES, "validation")
os.makedirs(OUT, exist_ok=True)

markers = pd.read_csv(os.path.join(KOF, "meta", "marker_ko_table.tsv"), sep="\t", dtype=str)
assert len(markers) == 68
MK = dict(zip(markers.ko, markers.gene))
NO_THR = set(markers.loc[markers["kofam_threshold_2026-10"] == "-", "ko"])
assert NO_THR == {"K23996", "K05915", "K18049"}, NO_THR

best = pd.read_csv(os.path.join(HERE, "best_per_species_558.tsv"), sep="\t")
MAGS = list(best["name"])
assert len(MAGS) == 558 and len(set(MAGS)) == 558

# ---------------------------------------------------------------- DRAM annotations
print("loading DRAM annotations ...", flush=True)
cols = ["fasta", "scaffold", "gene_position", "start_position", "end_position", "strandedness", "ko_id"]
dram = pd.read_csv(os.path.join(RES, "merged", "annotations.tsv.gz"), sep="\t", index_col=0, usecols=[0] + list(range(1, 9)),
                   dtype=str)
dram = dram[cols].copy()
dram.index.name = "dram_gene"
dram["start_position"] = dram.start_position.astype(int)
dram["end_position"] = dram.end_position.astype(int)
dram["strand"] = dram.strandedness.map({"1": "+", "-1": "-"})
assert dram.strand.notna().all()
assert set(dram.fasta) == set(MAGS), "DRAM MAG set != 558 list"
assert not dram.ko_id.fillna("").str.contains(",").any()

# ---------------------------------------------------------------- 1. gene counts
print("gene counts ...", flush=True)
n_dram = dram.groupby("fasta").size().rename("n_genes_dram")
# DRAM run_all.log 'genes=' (from each MAG's genes.faa)
log_counts = {}
with open(os.path.join(RES, "logs", "run_all.log")) as fh:
    for line in fh:
        m = re.search(r"done (\S+) genes=(\d+)", line)
        if m:
            log_counts[m.group(1)] = int(m.group(2))
n_log = pd.Series(log_counts, name="n_genes_dram_log")

kg = pd.read_csv(os.path.join(KOF, "results", "genes_all.tsv.gz"), sep="\t",
                 usecols=["gene", "mag", "contig", "start", "end", "strand", "partial", "aa_len"], dtype=str)
kg = kg[kg.mag.isin(MAGS)].copy()
assert set(kg.mag) == set(MAGS), "KOfam genes missing for some of the 558"
kg["start"] = kg.start.astype(int)
kg["end"] = kg.end.astype(int)
kg["aa_len"] = kg.aa_len.astype(int)
n_kof = kg.groupby("mag").size().rename("n_genes_kofam")

gc = pd.concat([n_dram, n_log, n_kof], axis=1)
gc.index.name = "mag"
gc["delta_dram_minus_kofam"] = gc.n_genes_dram - gc.n_genes_kofam
gc["log_matches_annotations"] = gc.n_genes_dram == gc.n_genes_dram_log

# ---------------------------------------------------------------- 2. coordinate join
print("coordinate join ...", flush=True)
dkey = dram.reset_index()[["dram_gene", "fasta", "scaffold", "start_position", "end_position", "strand", "ko_id"]]
dkey.columns = ["dram_gene", "mag", "contig", "start", "end", "strand", "dram_ko"]
j = dkey.merge(kg, on=["mag", "contig", "start", "end", "strand"], how="outer", indicator=True)
unmatched = j[j._merge != "both"]
gc["n_genes_coord_matched"] = j[j._merge == "both"].groupby("mag").size().reindex(gc.index).fillna(0).astype(int)
gc["n_dram_only"] = j[j._merge == "left_only"].groupby("mag").size().reindex(gc.index).fillna(0).astype(int)
gc["n_kofam_only"] = j[j._merge == "right_only"].groupby("mag").size().reindex(gc.index).fillna(0).astype(int)
gc.to_csv(os.path.join(OUT, "gene_counts_dram_vs_kofam.tsv"), sep="\t")
unmatched.to_csv(os.path.join(OUT, "genes_unmatched_by_coordinates.tsv"), sep="\t", index=False)

# ---------------------------------------------------------------- 3. per-gene marker KO
print("per-gene best KO (KOfam per_mag tables, pre-rescue) ...", flush=True)
parts = []
for mag in MAGS:
    p = os.path.join(KOF, "results", "per_mag", f"{mag}.kofam.tsv.gz")
    t = pd.read_csv(p, sep="\t", dtype=str,
                    usecols=["gene", "contig", "start", "end", "strand", "kofam_id", "kofam_score", "kofam_threshold",
                             "kofam_full_evalue", "kofam_all_passing"])
    t["mag"] = mag
    parts.append(t)
pm = pd.concat(parts, ignore_index=True)
pm["start"] = pm.start.astype(int)
pm["end"] = pm.end.astype(int)
# per_mag must hold exactly the genes_all genes
assert len(pm) == len(kg), (len(pm), len(kg))

jj = dkey.merge(pm, on=["mag", "contig", "start", "end", "strand"], how="inner")
assert len(jj) == (j._merge == "both").sum()
jj["dram_ko"] = jj.dram_ko.fillna("")
jj["kofam_id"] = jj.kofam_id.fillna("")
jj["kofam_all_passing"] = jj.kofam_all_passing.fillna("")
is_mk = jj.dram_ko.isin(MK) | jj.kofam_id.isin(MK)
# all KOs, not only markers: where do the two best-KO choices differ, and why
allko = jj[(jj.dram_ko != "") | (jj.kofam_id != "")]
allko_dis = allko[allko.dram_ko != allko.kofam_id].copy()
allko_dis["dram_ko_passes_in_kofam"] = [d in {x.split(":")[0] for x in a.split(";") if x}
                                        for d, a in zip(allko_dis.dram_ko, allko_dis.kofam_all_passing)]
allko_dis.to_csv(os.path.join(OUT, "all_ko_best_ko_disagreements.tsv.gz"), sep="\t", index=False)
gm = jj[is_mk].copy()
gm["dram_gene_sym"] = gm.dram_ko.map(MK).fillna(gm.dram_ko)
gm["kofam_gene_sym"] = gm.kofam_id.map(MK).fillna(gm.kofam_id)
gm["agree"] = gm.dram_ko == gm.kofam_id


def reason(r):
    if r.agree:
        return ""
    # per_mag kofam_all_passing lists KO ids only (scores, if present, after ':')
    passing = {x.split(":")[0]: (float(x.split(":")[1]) if ":" in x else float("nan"))
               for x in r.kofam_all_passing.split(";") if x}
    if r.dram_ko == "":
        return "DRAM no KO; KOfam best %s" % r.kofam_id
    if r.dram_ko not in passing:
        return "DRAM KO %s does not pass in KOfam run" % r.dram_ko
    return "both KOs pass in KOfam run (%s %.1f vs %s %.1f); best-KO tie-break" % (
        r.dram_ko, passing[r.dram_ko], r.kofam_id, passing.get(r.kofam_id, float("nan")))


gm["reason"] = gm.apply(reason, axis=1)
gm.to_csv(os.path.join(OUT, "marker_genes_dram_vs_kofam.tsv.gz"), sep="\t", index=False)
dis = gm[~gm.agree]
dis.to_csv(os.path.join(OUT, "marker_gene_disagreements.tsv"), sep="\t", index=False)

# unmatched genes carrying a marker on either side
um_mk = unmatched[unmatched.dram_ko.isin(MK)]
um_mk_k = unmatched[unmatched._merge == "right_only"].merge(pm[["mag", "gene", "kofam_id"]], on=["mag", "gene"], how="left")
um_mk_k = um_mk_k[um_mk_k.kofam_id.isin(MK)]

# ---------------------------------------------------------------- 4. per-MAG marker counts
print("per-MAG marker counts ...", flush=True)
d_cnt = (dram[dram.ko_id.isin(MK)].groupby(["fasta", "ko_id"]).size().unstack(fill_value=0)
         .reindex(index=MAGS, columns=list(MK), fill_value=0))
kc = pd.read_csv(os.path.join(KOF, "results", "all_mag_ko_counts.tsv"), sep="\t", dtype={"mag": str, "ko": str})
kc = kc[kc.mag.isin(MAGS) & kc.ko.isin(MK)]
k_cnt = kc.pivot_table(index="mag", columns="ko", values="n_best", aggfunc="sum", fill_value=0) \
          .reindex(index=MAGS, columns=list(MK), fill_value=0)
diff_pre = (d_cnt - k_cnt)

# second pre-rescue source: the poplar postprocessor's per-MAG table (n_<sym> for all 68 markers)
amg = pd.read_csv(os.path.join(KOF, "results", "all_mag_markers_guilds.tsv"), sep="\t").set_index("mag").loc[MAGS]
k_cnt2 = pd.DataFrame({ko: amg[f"n_{sym}"] for ko, sym in MK.items()})
pre_sources_identical = bool((k_cnt2.values == k_cnt.values).all())

# post-rescue: guilds_1324 n_<sym> where it has the column (the rescue only changes grdH/grdB/grdA/grdF), else pre-rescue
g1324 = pd.read_csv(os.path.join(KOF, "results", "guilds_1324.tsv"), sep="\t").set_index("mag").loc[MAGS]
ALT = {"cntA/yeaW": "cntA", "cntB/yeaX": "cntB"}
post = k_cnt.copy()
post_cols = []
for ko, sym in MK.items():
    c = f"n_{ALT.get(sym, sym)}"
    if c in g1324:
        post[ko] = g1324[c].astype(int)
        post_cols.append(sym)
missing_cols = [sym for sym in MK.values() if sym not in post_cols]
non_sec_changed = [MK[ko] for ko in MK if MK[ko] not in ("grdH", "grdB", "grdA", "grdF") and (post[ko] != k_cnt[ko]).any()]
diff_post = (d_cnt - post)

long_rows = []
for ko, sym in MK.items():
    long_rows.append(dict(ko=ko, gene=sym, threshold=markers.set_index("ko").loc[ko, "kofam_threshold_2026-10"],
                          dram_genes=int(d_cnt[ko].sum()), kofam_genes_pre_rescue=int(k_cnt[ko].sum()),
                          kofam_genes_post_rescue=int(post[ko].sum()),
                          dram_mags=int((d_cnt[ko] > 0).sum()), kofam_mags_pre_rescue=int((k_cnt[ko] > 0).sum()),
                          kofam_mags_post_rescue=int((post[ko] > 0).sum()),
                          mags_differing_pre=int((diff_pre[ko] != 0).sum()),
                          mags_differing_post=int((diff_post[ko] != 0).sum())))
per_ko = pd.DataFrame(long_rows)
per_ko.to_csv(os.path.join(OUT, "marker_ko_totals_dram_vs_kofam.tsv"), sep="\t", index=False)

cells_post = diff_post.stack()
cells_post = cells_post[cells_post != 0].rename("dram_minus_kofam_post_rescue").reset_index()
cells_post.columns = ["mag", "ko", "dram_minus_kofam_post_rescue"]
cells_post["gene"] = cells_post.ko.map(MK)
cells_post["dram"] = [d_cnt.at[m, k] for m, k in zip(cells_post.mag, cells_post.ko)]
cells_post["kofam_pre_rescue"] = [k_cnt.at[m, k] for m, k in zip(cells_post.mag, cells_post.ko)]
cells_post["kofam_post_rescue"] = [post.at[m, k] for m, k in zip(cells_post.mag, cells_post.ko)]
cells_post.to_csv(os.path.join(OUT, "mag_marker_cells_differing_post_rescue.tsv"), sep="\t", index=False)

# ---------------------------------------------------------------- 5. distill sheet consistency
print("distill sheet ...", flush=True)
sheet = pd.read_excel(os.path.join(RES, "distill", "metabolism_summary.xlsx"), sheet_name="Salterns methane")
assert len(sheet) == 65
genome_cols = [c for c in sheet.columns if c not in ("gene_id", "gene_description", "module", "header", "subheader")]
assert set(genome_cols) == set(MAGS), "sheet genome columns != 558"
sh = sheet.set_index("gene_id")[MAGS].astype(int)
ann_cnt = (dram[dram.ko_id.isin(sh.index)].groupby(["ko_id", "fasta"]).size().unstack(fill_value=0)
           .reindex(index=sh.index, columns=MAGS, fill_value=0))
sheet_mismatch = int((sh != ann_cnt).values.sum())
sheet_kos = set(sh.index)
assert sheet_kos | NO_THR == set(MK), "sheet KOs + 3 dash KOs != 68 markers"

# ---------------------------------------------------------------- 6. genome_stats
gs = pd.read_csv(os.path.join(RES, "distill", "genome_stats.tsv"), sep="\t").set_index("genome").loc[MAGS]
tax_in = pd.read_csv(os.path.join(HERE, "gtdb_taxonomy_558.tsv"), sep="\t").set_index("user_genome").loc[MAGS]
q_in = pd.read_csv(os.path.join(HERE, "checkm_quality_558.tsv"), sep="\t").set_index("Bin Id").loc[MAGS]
gs_ok = dict(
    all_tax_start_d=bool(gs.taxonomy.str.startswith("d__").all()),
    tax_equal_input=bool((gs.taxonomy == tax_in.classification).all()),
    completeness_equal_input=bool(((gs["completeness score"] - q_in.Completeness).abs() < 1e-6).all()),
    contamination_equal_input=bool(((gs["contamination score"] - q_in.Contamination).abs() < 1e-6).all()),
    n_completeness_0=int((gs["completeness score"] == 0).sum()),
    n_contamination_100=int((gs["contamination score"] == 100).sum()),
)

# ---------------------------------------------------------------- report
lines = []
P = lines.append
P("# DRAM vs KOfam validation (generated by analysis/validate_vs_kofam.py)")
P("")
P(f"MAGs: DRAM {dram.fasta.nunique()}, KOfam genes for {kg.mag.nunique()} of the 558")
P(f"genes: DRAM {len(dram):,}; KOfam (same 558) {len(kg):,}; run_all.log sum {int(n_log.reindex(MAGS).sum()):,}")
P(f"MAGs with DRAM gene count == KOfam: {(gc.delta_dram_minus_kofam == 0).sum()} / 558; "
  f"delta distribution: {gc.delta_dram_minus_kofam.value_counts().sort_index().to_dict()}")
P(f"MAGs where annotations.tsv rows == run_all.log genes=: {gc.log_matches_annotations.sum()} / 558")
P(f"coordinate join: matched {int((j._merge == 'both').sum()):,}; DRAM-only {int((j._merge == 'left_only').sum()):,}; "
  f"KOfam-only {int((j._merge == 'right_only').sum()):,}")
P(f"unmatched genes carrying a marker KO: DRAM-side {len(um_mk)}, KOfam-side {len(um_mk_k)}")
P("")
P(f"all KOs: genes with a KO on either side {len(allko):,} (DRAM {int((jj.dram_ko != '').sum()):,}, KOfam "
  f"{int((jj.kofam_id != '').sum()):,}); best KO differs for {len(allko_dis):,}; of these DRAM-only KO "
  f"{int((allko_dis.kofam_id == '').sum())}, KOfam-only KO {int((allko_dis.dram_ko == '').sum())}, both called but "
  f"different {int(((allko_dis.dram_ko != '') & (allko_dis.kofam_id != '')).sum())} (DRAM's KO passes in KOfam run: "
  f"{int(allko_dis.dram_ko_passes_in_kofam.sum())}; KOfam full E-value of its best KO == 0: "
  f"{int((allko_dis.kofam_full_evalue.astype(float) == 0).sum())})")
P(f"marker genes (either side has one of the 68 markers as best KO): {len(gm):,}; agree {int(gm.agree.sum()):,}; "
  f"disagree {len(dis)}")
if len(dis):
    P("disagreement reasons: " + str(dis.reason.str.replace(r"\(.*\)", "", regex=True).value_counts().to_dict()))
    P(dis[["mag", "dram_gene", "dram_ko", "kofam_id", "kofam_score", "kofam_full_evalue", "kofam_all_passing", "reason"]]
      .to_string(index=False))
P("")
P(f"per-MAG x marker cells (558 x 68 = {558 * 68}): differing vs KOfam pre-rescue {int((diff_pre != 0).values.sum())}; "
  f"vs KOfam post-rescue {int((diff_post != 0).values.sum())}")
P(f"pre-rescue sources (all_mag_ko_counts n_best vs all_mag_markers_guilds n_<marker>) identical: {pre_sources_identical}")
P(f"markers without an n_ column in guilds_1324 (post-rescue = pre-rescue for these): {len(missing_cols)}; "
  f"non-Sec markers whose guilds_1324 count differs from pre-rescue: {non_sec_changed}")
P(per_ko[(per_ko.mags_differing_pre > 0) | (per_ko.mags_differing_post > 0)].to_string(index=False))
P("")
P(f"'Salterns methane' sheet: 65 KOs x 558 genomes; cells != annotations ko_id counts: {sheet_mismatch}")
P(f"dash-threshold KOs (mpnS/phnY*/phnZ1) called by DRAM: {int(dram.ko_id.isin(NO_THR).sum())}; by KOfam best KO: "
  f"{int(pm.kofam_id.isin(NO_THR).sum())}")
P(f"genome_stats join: {gs_ok}")
rep = "\n".join(lines)
print(rep)
with open(os.path.join(OUT, "validation_report.txt"), "w") as fh:
    fh.write(rep + "\n")
