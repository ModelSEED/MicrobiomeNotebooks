#!/usr/bin/env python
"""Per-MAG marker counts and guild calls (revised rules, annotation_choice/README.md section 6.3) for the 1,324
re-assembly MAGs (KOfam 2026-10-26, full DB, kofam track) and, with the same code, the 310 original MAGs (KOfam
2026-10-26, 68 marker profiles, originals_markers.py), both with the Sec read-through rescue (sec_rescue.py).

Rules (per MAG, on each gene's best KO; after the Sec rescue):
  methanogen        mcrA K00399
  mepn_consumer     phnJ K06163                      (+ phnG/H/I/K/L/M counts, n distinct of phnGHIJKLM, n contigs)
  betaine_reducer   grdH K21579 or grdI K21578       (grdH incl. Sec read-through; grdA/C/D reported as support)
  methylamine_user  mttB K14083 / mtbB K16178 / mtmB K16176, only in an mcrA genome
                    (mttB_family_without_mcrA incl. mtgB K25220: the non-Pyl MttB family; pylBCDS count reported)
  tma_producer      cutC K20038 + cutD K20037, or cntA K22443 + cntB K22444; cutC- or cntA-only = tma_candidate
Flag: <guild>_short_maxbin2_only = the call disappears when genes on MaxBin2 contigs < 2,500 bp are dropped.

Writes (results/): guilds_1324.tsv, guilds_species_clusters_new.tsv, originals_310/guilds_310_kofam.tsv,
and prints the cross-checks (pre-rescue calls vs the kofam track's all_mag_markers_guilds.tsv; best KO among the
68 markers vs over all profiles on the new MAGs).
Run: ~/Documents/py_venv/bin/python -I build_guilds.py
"""
from collections import defaultdict
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
K = HERE.parent
RES = K / "results"
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
DEFINING = {"methanogen": ["mcrA"], "mepn_consumer": ["phnJ"], "betaine_reducer": ["grdH", "grdI"],
            "methylamine_user": ["mttB", "mtbB", "mtmB", "mcrA"], "tma_producer": ["cutC", "cutD", "cntA", "cntB"]}


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


def load_new():
    gk = pd.read_csv(RES / "all_genes_ko.tsv.gz", sep="\t",
                     usecols=["gene", "mag", "binner", "contig", "contig_len", "best_ko", "all_passing"])
    gk["anypass"] = gk.all_passing.str.split(";").map(lambda a: [x.split(":")[0] for x in a])
    keep = gk.best_ko.isin(KO2SYM) | gk.anypass.map(lambda a: any(k in KO2SYM for k in a))
    return gk.loc[keep, ["gene", "mag", "binner", "contig", "contig_len", "best_ko", "anypass"]]


def load_orig():
    h = pd.read_csv(RES / "originals_310/marker_hits.tsv.gz", sep="\t")
    p = h[h.passes == 1]
    ap = p.groupby("gene").ko.apply(list)
    g = p.drop_duplicates("gene").set_index("gene")[["mag", "binner", "contig", "contig_len", "gene_best_ko"]]
    g = g.rename(columns={"gene_best_ko": "best_ko"})
    g["anypass"] = ap
    return g.reset_index()


def marker_only_best(h):
    """best KO among the marker profiles (DRAM rule restricted to markers) from a marker-hit table."""
    p = h[h.passes == 1].sort_values(["gene", "full_evalue", "full_score", "ko"], ascending=[True, True, False, True])
    p = p.drop_duplicates("gene")
    return p.rename(columns={"ko": "best_ko"}).assign(anypass=lambda d: d.best_ko.map(lambda k: [k]))


def main():
    meta = pd.read_csv(R26 / "comparison/mags_old_vs_new.tsv", sep="\t")
    meta["name"] = meta.genome.str.split("|", n=1).str[1]
    meta = meta.set_index("name")
    man = pd.read_csv(K / "mag_manifest.tsv", sep="\t").set_index("name")
    sec = pd.read_csv(RES / "sec_readthrough/sec_readthrough_junctions_scored.tsv", sep="\t")
    info = pd.read_csv(RES / "hmm_info.tsv", sep="\t").set_index("ko")
    tax = ["domain", "phylum", "class", "order", "family", "genus", "species"]

    # ---------------- new 1,324
    gn = load_new()
    mags = man.index.tolist()
    pre = per_mag(gn.assign(source="prodigal"), mags).set_index("mag")
    gt, replaced = gene_calls("new", gn, sec)
    post = per_mag(gt, mags).set_index("mag")
    # cross-check the pre-rescue calls against the kofam track's postprocess.py output
    pp = pd.read_csv(RES / "all_mag_markers_guilds.tsv", sep="\t").set_index("mag")
    chk = ["methanogen", "mepn_consumer", "phnGHIJKLM_n_distinct", "phnGHIJKLM_n_contigs", "betaine_reducer",
           "grdACD_n_distinct", "methylamine_user", "mttB_family_without_mcrA", "tma_producer",
           "tma_candidate_one_sided"] + [f"{g}_anypass" for g in GUILDS]
    diff = {c: int((pre.loc[pp.index, c] != pp[c]).sum()) for c in chk}
    lost_pp = pp.guilds_lost_without_short_maxbin2_contigs.fillna("")
    lost_me = pre.apply(lambda r: ";".join(g for g in GUILDS if r[f"{g}_short_maxbin2_only"]), axis=1)
    diff["guilds_lost_without_short_maxbin2_contigs"] = int((lost_me.loc[pp.index] != lost_pp).sum())
    print("pre-rescue vs kofam-track postprocess, MAGs differing per column:", diff)
    assert sum(diff.values()) == 0

    out = post.copy()
    out["betaine_reducer_without_sec_rescue"] = pre.betaine_reducer
    for s in ("grdH", "grdB", "grdA", "grdF", "grdI"):
        out[f"n_{s}_prodigal_only"] = pre[f"n_{s}"]
    m2 = meta.reindex(out.index)
    out.insert(0, "sample", man.loc[out.index, "sample"])
    out.insert(1, "binner", man.loc[out.index, "binner"])
    out.insert(2, "species_cluster", m2.species_cluster)
    for i, c in enumerate(["completeness", "contamination", "quality", "mimag", "n_contigs", "size_bp", "n50"] + tax):
        out.insert(3 + i, c, m2[c])
    out.reset_index().to_csv(RES / "guilds_1324.tsv", sep="\t", index=False)
    print("new: guild MAGs pre/post rescue:",
          {g: (int(pre[g].sum()), int(post[g].sum())) for g in GUILDS + ["tma_candidate_one_sided", "mttB_family_without_mcrA"]})

    # marker-only best KO vs all-profile best KO on the new MAGs (comparability of the originals' marker-only run)
    h = pd.read_csv(RES / "all_marker_hits.tsv.gz", sep="\t")
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
        RES / "marker_only_vs_full_db_best_ko.tsv", sep="\t")

    # Sec rescue: replaced Prodigal calls (bookkeeping)
    rep = replaced.assign(sym=replaced.best_ko.map(KO2SYM))
    rep.drop(columns=["anypass"]).to_csv(RES / "sec_readthrough/replaced_prodigal_piece_calls_new.tsv", sep="\t",
                                         index=False)

    # ---------------- species-cluster summary (new)
    out2 = out.reset_index()
    # best MAG per cluster = the project's flag (annotation_choice/short_contig_fraction_per_mag.tsv: highest
    # completeness - 5 x contamination; it differs from a max-quality/N50 tie-break only in 3 exact-quality ties)
    bf = pd.read_csv(R26 / "annotation_choice/short_contig_fraction_per_mag.tsv", sep="\t")
    bset = set(bf.loc[bf.best_of_cluster.astype(bool), "mag"])
    best = out2[out2.mag.isin(bset)]
    assert best.species_cluster.is_unique and len(best) == out2.species_cluster.nunique() == 558
    agg = {"n_mags": ("mag", "size")}
    for g in GUILDS + ["tma_candidate_one_sided", "mttB_family_without_mcrA"]:
        agg[f"{g}_n_mags"] = (g, "sum")
    cs = out2.groupby("species_cluster").agg(**agg)
    b = best.set_index("species_cluster")
    cs["best_mag"] = b.mag
    for c in ["completeness", "contamination", "quality", "n50"] + tax:
        cs[f"best_{c}" if c in ("completeness", "contamination", "quality", "n50") else c] = b[c]
    for g in GUILDS + ["tma_candidate_one_sided", "mttB_family_without_mcrA"]:
        cs[f"{g}_best_mag"] = b[g]
        cs[f"{g}_any_mag"] = (cs[f"{g}_n_mags"] > 0).astype(int)
    for g in GUILDS:
        cs[f"{g}_best_mag_short_maxbin2_only"] = b[f"{g}_short_maxbin2_only"]
    cs["betaine_reducer_evidence_best_mag"] = b.betaine_reducer_evidence
    cs["phnGHIJKLM_n_distinct_best_mag"] = b.phnGHIJKLM_n_distinct
    cs.reset_index().to_csv(RES / "guilds_species_clusters_new.tsv", sep="\t", index=False)
    print("new species clusters:", len(cs), {g: (int(cs[f"{g}_best_mag"].sum()), int(cs[f"{g}_any_mag"].sum()))
                                            for g in GUILDS + ["tma_candidate_one_sided"]})

    # ---------------- originals 310 (KOfam 2026-10, 68 markers, same rules + Sec rescue)
    go = load_orig()
    omags = sorted(meta[meta.set.str.startswith("original")].index)
    opre = per_mag(go.assign(source="prodigal"), omags).set_index("mag")
    ogt, orep = gene_calls("orig", go, sec)
    opost = per_mag(ogt, omags).set_index("mag")
    opost["betaine_reducer_without_sec_rescue"] = opre.betaine_reducer
    for s in ("grdH", "grdB", "grdA", "grdF", "grdI"):
        opost[f"n_{s}_prodigal_only"] = opre[f"n_{s}"]
    m3 = meta.loc[opost.index]
    opost.insert(0, "sample", m3["sample"])
    opost.insert(1, "species_cluster", m3.species_cluster)
    for i, c in enumerate(["completeness", "contamination", "quality", "mimag", "n_contigs", "size_bp", "n50"] + tax):
        opost.insert(2 + i, c, m3[c])
    opost.reset_index().to_csv(RES / "originals_310/guilds_310_kofam.tsv", sep="\t", index=False)
    print("orig: guild MAGs pre/post rescue:",
          {g: (int(opre[g].sum()), int(opost[g].sum())) for g in GUILDS + ["tma_candidate_one_sided", "mttB_family_without_mcrA"]})


if __name__ == "__main__":
    main()
