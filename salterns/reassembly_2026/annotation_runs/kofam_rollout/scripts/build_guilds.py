#!/usr/bin/env python
"""Per-MAG marker counts and guild calls (revised rules, annotation_choice/README.md section 6.3) for the 935 kept
MAGs of the Arm C rollout (KOfam 2026-10-26, full DB, kofam_rollout track), with the Sec read-through rescue
(sec_rescue.py, set name 'rollout').

kofam_rollout copy of kofam/scripts/build_guilds.py (2026-11-19). Changes:
  * MAG list and metadata from kofam_rollout/mag_manifest.tsv (935 MAGs; CheckM2, GTDB R232 from rollout/
    results/mags_kept.tsv), not comparison/mags_old_vs_new.tsv; no 1,324 / 558 asserts;
  * output results/guilds_935.tsv (one row per MAG);
  * species-cluster summary only when the cluster table exists (default
    reassembly_2026/comparison_rollout/mags_three_way.tsv, rows set == 'rollout', genome 'rollout|<name>'):
    results/guilds_935_species_clusters.tsv. Best MAG per cluster = the table's rollout_species_rep flag
    (= highest CheckM2 completeness - 5 x contamination among the cluster's rollout MAGs); asserted to be one per
    cluster and to cover every cluster, whatever the number of clusters;
  * the originals-310 block is dropped (no original-set run in this track);
  * --res / --manifest / --clusters / --sec / --out make it runnable on the validation output.
The rules, per_mag(), gene_calls() and the pre-rescue cross-check against postprocess.py are unchanged.

Rules (per MAG, on each gene's best KO; after the Sec rescue):
  methanogen        mcrA K00399
  mepn_consumer     phnJ K06163                      (+ phnG/H/I/K/L/M counts, n distinct of phnGHIJKLM, n contigs)
  betaine_reducer   grdH K21579 or grdI K21578       (grdH incl. Sec read-through; grdA/C/D reported as support)
  methylamine_user  mttB K14083 / mtbB K16178 / mtmB K16176, only in an mcrA genome
                    (mttB_family_without_mcrA incl. mtgB K25220: the non-Pyl MttB family; pylBCDS count reported)
  tma_producer      cutC K20038 + cutD K20037, or cntA K22443 + cntB K22444; cutC- or cntA-only = tma_candidate
Flag: <guild>_short_maxbin2_only = the call disappears when genes on MaxBin2 contigs < 2,500 bp are dropped.

Run (from annotation_runs/kofam_rollout/):  ~/Documents/py_venv/bin/python -I scripts/build_guilds.py
"""
import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
K = HERE.parent
R26 = K.parents[1]                                     # reassembly_2026/
SHORT = 2500
SYM = {   # guild symbols -> KO
    "mcrA": "K00399", "phnG": "K06166", "phnH": "K06165", "phnI": "K06164", "phnJ": "K06163", "phnK": "K05781",
    "phnL": "K05780", "phnM": "K06162", "grdH": "K21579", "grdI": "K21578", "grdA": "K10670", "grdC": "K21577",
    "grdD": "K21576", "grdB": "K10672", "grdE": "K10671", "grdF": "K21584", "grdG": "K21583",
    "mttB": "K14083", "mtbB": "K16178", "mtmB": "K16176", "mtgB": "K25220", "cutC": "K20038", "cutD": "K20037",
    "cntA": "K22443", "cntB": "K22444", "pylB": "K16180", "pylC": "K16181", "pylD": "K16182", "pylS": "K11627",
}
KO2SYM = {v: k for k, v in SYM.items()}
PHN7 = ["phnG", "phnH", "phnI", "phnJ", "phnK", "phnL", "phnM"]
GUILDS = ["methanogen", "mepn_consumer", "betaine_reducer", "methylamine_user", "tma_producer"]
SEC_COLS = ["set", "status", "pieces", "junction_id", "mag", "binner", "contig", "contig_len", "best_ko_all",
            "all_sig_kos"]
TAX = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def guild_calls(c, con):
    """c: symbol -> n genes; con: symbol -> set of contigs."""
    g = {}
    g["methanogen"] = c["mcrA"] > 0
    g["mepn_consumer"] = c["phnJ"] > 0
    g["phnGHIJKLM_n_distinct"] = sum(c[s] > 0 for s in PHN7)
    g["phnGHIJKLM_n_contigs"] = len(set().union(*(con[s] for s in PHN7)))
    g["betaine_reducer"] = c["grdH"] > 0 or c["grdI"] > 0
    g["grdACD_n_distinct"] = sum(c[s] > 0 for s in ("grdA", "grdC", "grdD"))
    mt = c["mttB"] + c["mtbB"] + c["mtmB"]
    g["methylamine_user"] = mt > 0 and c["mcrA"] > 0
    g["mttB_family_without_mcrA"] = (mt + c["mtgB"]) > 0 and c["mcrA"] == 0
    cut, cnt = c["cutC"] > 0 and c["cutD"] > 0, c["cntA"] > 0 and c["cntB"] > 0
    g["tma_producer"] = cut or cnt
    g["tma_candidate_one_sided"] = (not (cut or cnt)) and (c["cutC"] + c["cntA"]) > 0
    return g


def gene_calls(set_name, genes_best, sec):
    """genes_best: DataFrame gene, mag, binner, contig, contig_len, best_ko, anypass (list of KOs).
    sec: accepted Sec junctions of this set. Returns the gene-level table after the rescue."""
    acc = sec[sec.status.isin(["sec_readthrough", "confirmed"]) & (sec.set == set_name)]
    drop = set(p for ps in acc.pieces.fillna("") for p in ps.split(";") if p)
    gb = genes_best[~genes_best.gene.isin(drop)].assign(source="prodigal")
    add = pd.DataFrame(dict(gene=acc.junction_id, mag=acc.mag, binner=acc.binner, contig=acc.contig,
                            contig_len=acc.contig_len, best_ko=acc.best_ko_all,
                            anypass=acc.all_sig_kos.fillna("").str.split(";"),
                            source=acc.status))
    replaced = genes_best[genes_best.gene.isin(drop)]
    return pd.concat([gb, add], ignore_index=True), replaced


def per_mag(gt, mags):
    gt = gt.copy()
    gt["sym"] = gt.best_ko.map(KO2SYM)
    gt["short"] = (gt.binner == "maxbin2") & (gt.contig_len.astype(int) < SHORT)
    rows = []
    by = {m: x for m, x in gt[gt.sym.notna() | gt.anypass.map(lambda a: any(k in KO2SYM for k in a))].groupby("mag")}
    for m in mags:
        x = by.get(m, gt.iloc[0:0])
        c, con, cl, conl, ca, cona = (defaultdict(int), defaultdict(set), defaultdict(int), defaultdict(set),
                                      defaultdict(int), defaultdict(set))
        sec_c = defaultdict(int)
        short_genes = []
        for r in x.itertuples():
            if isinstance(r.sym, str):
                c[r.sym] += 1; con[r.sym].add(r.contig)
                if not r.short:
                    cl[r.sym] += 1; conl[r.sym].add(r.contig)
                else:
                    short_genes.append(f"{r.sym}@{r.contig}({r.contig_len})")
                if r.source == "sec_readthrough":
                    sec_c[r.sym] += 1
            for k in r.anypass:
                if k in KO2SYM:
                    ca[KO2SYM[k]] += 1; cona[KO2SYM[k]].add(r.contig)
        g, gl, ga = guild_calls(c, con), guild_calls(cl, conl), guild_calls(ca, cona)
        row = {"mag": m}
        row.update({k: int(v) if isinstance(v, bool) else v for k, v in g.items()})
        for k in GUILDS:
            row[f"{k}_short_maxbin2_only"] = int(g[k] and not gl[k])
            row[f"{k}_anypass"] = int(ga[k])
        row["tma_candidate_short_maxbin2_only"] = int(g["tma_candidate_one_sided"] and not gl["tma_candidate_one_sided"]
                                                     and not gl["tma_producer"])
        ev = [s for s in ("grdH", "grdI") if c[s] > 0]
        if sec_c["grdH"] and c["grdH"] == sec_c["grdH"]:
            ev = ["grdH_sec_readthrough" if s == "grdH" else s for s in ev]
        row["betaine_reducer_evidence"] = "+".join(ev)
        row["phn_accessory_counts"] = ";".join(f"{s}:{c[s]}" for s in PHN7 if s != "phnJ")
        row["pylBCDS_n_distinct"] = sum(c[s] > 0 for s in ("pylB", "pylC", "pylD", "pylS"))
        row["marker_genes_on_short_maxbin2_contigs"] = ";".join(short_genes)
        for s in SYM:
            row[f"n_{s}"] = c[s]
        for s in ("grdH", "grdB", "grdA", "grdF"):
            row[f"n_{s}_sec_readthrough"] = sec_c[s]
        rows.append(row)
    return pd.DataFrame(rows)


def load_new(res):
    gk = pd.read_csv(res / "all_genes_ko.tsv.gz", sep="\t",
                     usecols=["gene", "mag", "binner", "contig", "contig_len", "best_ko", "all_passing"])
    gk["anypass"] = gk.all_passing.str.split(";").map(lambda a: [x.split(":")[0] for x in a])
    keep = gk.best_ko.isin(KO2SYM) | gk.anypass.map(lambda a: any(k in KO2SYM for k in a))
    return gk.loc[keep, ["gene", "mag", "binner", "contig", "contig_len", "best_ko", "anypass"]]


def marker_only_best(h):
    """best KO among the marker profiles (DRAM rule restricted to markers) from a marker-hit table."""
    p = h[h.passes == 1].sort_values(["gene", "full_evalue", "full_score", "ko"], ascending=[True, True, False, True])
    p = p.drop_duplicates("gene")
    return p.rename(columns={"ko": "best_ko"}).assign(anypass=lambda d: d.best_ko.map(lambda k: [k]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", default=str(K / "results"), help="results dir copied back from poplar")
    ap.add_argument("--manifest", default=str(K / "mag_manifest.tsv"))
    ap.add_argument("--clusters", default=str(R26 / "comparison_rollout/mags_three_way.tsv"),
                    help="species-cluster table (genome, set, species_cluster, rollout_species_rep); skipped if absent")
    ap.add_argument("--sec", default=None, help="Sec rescue table (default <res>/sec_readthrough/"
                    "sec_readthrough_junctions_scored.tsv); 'none' = no rescue (testing only)")
    ap.add_argument("--set-name", default="rollout")
    ap.add_argument("--out", default="guilds_935.tsv")
    ap.add_argument("--n-mags", type=int, default=935, help="expected MAG count (0 = do not check)")
    a = ap.parse_args()
    res = Path(a.res)
    man = pd.read_csv(a.manifest, sep="\t").set_index("name")
    mags = man.index.tolist()
    if a.n_mags:
        assert len(mags) == a.n_mags and man.index.is_unique, (len(mags), a.n_mags)
    sec_p = Path(a.sec) if a.sec and a.sec != "none" else res / "sec_readthrough/sec_readthrough_junctions_scored.tsv"
    if a.sec == "none":
        sec = pd.DataFrame(columns=SEC_COLS)
        print("WARNING: no Sec rescue table (testing only)")
    else:
        sec = pd.read_csv(sec_p, sep="\t")
        assert set(sec.set) <= {a.set_name}, set(sec.set)

    # ---------------- per MAG
    gn = load_new(res)
    assert set(gn.mag) <= set(mags), sorted(set(gn.mag) - set(mags))[:5]
    pre = per_mag(gn.assign(source="prodigal"), mags).set_index("mag")
    gt, replaced = gene_calls(a.set_name, gn, sec)
    post = per_mag(gt, mags).set_index("mag")
    # cross-check the pre-rescue calls against this track's postprocess.py output
    pp = pd.read_csv(res / "all_mag_markers_guilds.tsv", sep="\t").set_index("mag")
    assert set(pp.index) == set(mags), (len(pp), len(mags))
    chk = ["methanogen", "mepn_consumer", "phnGHIJKLM_n_distinct", "phnGHIJKLM_n_contigs", "betaine_reducer",
           "grdACD_n_distinct", "methylamine_user", "mttB_family_without_mcrA", "tma_producer",
           "tma_candidate_one_sided"] + [f"{g}_anypass" for g in GUILDS]
    diff = {c: int((pre.loc[pp.index, c] != pp[c]).sum()) for c in chk}
    lost_pp = pp.guilds_lost_without_short_maxbin2_contigs.fillna("")
    lost_me = pre.apply(lambda r: ";".join(g for g in GUILDS if r[f"{g}_short_maxbin2_only"]), axis=1)
    diff["guilds_lost_without_short_maxbin2_contigs"] = int((lost_me.loc[pp.index] != lost_pp).sum())
    print("pre-rescue vs postprocess.py, MAGs differing per column:", diff)
    assert sum(diff.values()) == 0

    out = post.copy()
    out["betaine_reducer_without_sec_rescue"] = pre.betaine_reducer
    for s in ("grdH", "grdB", "grdA", "grdF", "grdI"):
        out[f"n_{s}_prodigal_only"] = pre[f"n_{s}"]
    meta_cols = [c for c in ["sample", "binner", "completeness", "contamination", "quality_score", "mimag",
                             "n_contigs", "size_bp", "n50"] + TAX if c in man.columns]
    clu = None
    if a.clusters != "none" and Path(a.clusters).exists():
        t = pd.read_csv(a.clusters, sep="\t")
        t = t[t.set == a.set_name].copy()
        t["name"] = t.genome.str.split("|", n=1).str[1]
        assert t.genome.str.startswith(f"{a.set_name}|").all() and t.name.is_unique
        assert set(t.name) == set(mags), ("cluster table rollout rows != manifest", len(t), len(mags))
        clu = t.set_index("name")
    else:
        print(f"no species-cluster table at {a.clusters}: per-cluster summary skipped (documented next step)")
    for i, c in enumerate(meta_cols):
        out.insert(i, c, man.loc[out.index, c])
    if clu is not None:
        out.insert(2, "species_cluster", clu.loc[out.index, "species_cluster"])
        out.insert(3, "rollout_species_rep", clu.loc[out.index, "rollout_species_rep"].astype(bool).astype(int))
    out.reset_index().to_csv(res / a.out, sep="\t", index=False)
    print(f"{a.out}: {len(out)} MAGs; guild MAGs pre/post rescue:",
          {g: (int(pre[g].sum()), int(post[g].sum())) for g in GUILDS + ["tma_candidate_one_sided", "mttB_family_without_mcrA"]})

    # marker-only best KO vs all-profile best KO (comparability with marker-only runs, e.g. the originals)
    h = pd.read_csv(res / "all_marker_hits.tsv.gz", sep="\t")
    mo = per_mag(marker_only_best(h).assign(source="prodigal")[["gene", "mag", "binner", "contig", "contig_len",
                                                                 "best_ko", "anypass", "source"]], mags).set_index("mag")
    mo_diff = {g: int((mo[g] != pre[g]).sum()) for g in GUILDS + ["tma_candidate_one_sided"]}
    hb = marker_only_best(h).set_index("gene").best_ko
    full = h.drop_duplicates("gene").set_index("gene").gene_best_ko.fillna("")
    gl = hb[hb.isin(KO2SYM)]
    nd = gl[gl != full.reindex(gl.index)]
    print("marker-only vs full-DB best KO, guild-KO genes:", len(gl), "differing:", len(nd),
          nd.map(KO2SYM).value_counts().to_dict(), "MAG guild calls differing:", mo_diff)
    pd.DataFrame({"marker_only_best_ko": nd, "full_db_best_ko": full.reindex(nd.index)}).to_csv(
        res / "marker_only_vs_full_db_best_ko.tsv", sep="\t")

    # Sec rescue: replaced Prodigal calls (bookkeeping)
    if len(sec):
        rep = replaced.assign(sym=replaced.best_ko.map(KO2SYM))
        rep.drop(columns=["anypass"]).to_csv(sec_p.parent / f"replaced_prodigal_piece_calls_{a.set_name}.tsv",
                                             sep="\t", index=False)

    # ---------------- species-cluster summary (only with the cluster table)
    if clu is None:
        return
    out2 = out.reset_index()
    best = out2[out2.rollout_species_rep == 1]
    n_clu = out2.species_cluster.nunique()
    assert best.species_cluster.is_unique and len(best) == n_clu, (len(best), n_clu)
    agg = {"n_mags": ("mag", "size")}
    for g in GUILDS + ["tma_candidate_one_sided", "mttB_family_without_mcrA"]:
        agg[f"{g}_n_mags"] = (g, "sum")
    cs = out2.groupby("species_cluster").agg(**agg)
    b = best.set_index("species_cluster")
    cs["best_mag"] = b.mag
    for c in ["completeness", "contamination", "quality_score", "n50"] + TAX:
        if c in b.columns:
            cs[f"best_{c}" if c in ("completeness", "contamination", "quality_score", "n50") else c] = b[c]
    for g in GUILDS + ["tma_candidate_one_sided", "mttB_family_without_mcrA"]:
        cs[f"{g}_best_mag"] = b[g]
        cs[f"{g}_any_mag"] = (cs[f"{g}_n_mags"] > 0).astype(int)
    for g in GUILDS:
        cs[f"{g}_best_mag_short_maxbin2_only"] = b[f"{g}_short_maxbin2_only"]
    cs["betaine_reducer_evidence_best_mag"] = b.betaine_reducer_evidence
    cs["phnGHIJKLM_n_distinct_best_mag"] = b.phnGHIJKLM_n_distinct
    if "cluster_sets" in clu.columns:
        cs.insert(0, "cluster_sets", clu.drop_duplicates("species_cluster").set_index("species_cluster")
                  .cluster_sets.reindex(cs.index))
    name = a.out.replace(".tsv", "") + "_species_clusters.tsv"
    cs.reset_index().to_csv(res / name, sep="\t", index=False)
    print(f"{name}: {len(cs)} species clusters;",
          {g: (int(cs[f"{g}_best_mag"].sum()), int(cs[f"{g}_any_mag"].sum())) for g in GUILDS + ["tma_candidate_one_sided"]})


if __name__ == "__main__":
    main()
