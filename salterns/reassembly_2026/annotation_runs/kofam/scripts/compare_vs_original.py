#!/usr/bin/env python
"""Guild calls, original 310 MAGs vs the 1,324 re-assembly MAGs, at species-cluster level (skani ANI >= 95 %,
AF >= 50 %; comparison/mags_old_vs_new.tsv species_cluster).

Three call sets:
  RAST      correlational/r232/guilds_r232.csv (regex on RAST function strings; mepn = 'c-p lyase', betaine =
            'betaine reductase component b' (GrdH or GrdI), methanogen = 'methyl-coenzyme M reductase',
            tma_producer = 'carnitine monooxygenase' or genes 'clustered with choline TMA-lyase' -- a candidate rule)
  KOfam-orig  results/originals_310/guilds_310_kofam.tsv (KOfam 2026-10, 68 markers, same rules + Sec rescue)
  KOfam-new   results/guilds_1324.tsv (KOfam 2026-10, full DB, + Sec rescue)
For TMA the RAST call is compared with KOfam 'tma_producer OR tma_candidate_one_sided' (the revised rule calls 0
producers on either set).

Per (cluster, guild), when any of RAST / KOfam-orig / KOfam-new calls it:
  recovered_same        RAST member and >= 1 re-assembly MAG of the cluster calls it
  not_recovered         RAST member, cluster has no re-assembly MAG
  lost_rule             RAST member, re-assembly MAGs present but none calls it, and KOfam does not call it in the
                        original genome either (RAST string vs KOfam threshold)
  lost_genome           as above but KOfam calls it in the original (re-assembly MAGs lack the gene / score)
  gained_new_cluster    called in the re-assembly, cluster has no original MAG
  gained_rule           called in the re-assembly and by KOfam in the original, not by RAST
  gained_genome         called in the re-assembly, not in the original by RAST or KOfam
  kofam_orig_only       KOfam calls it in the original only (re-assembly MAGs present or not), RAST does not
Review overrides (2026-11-12, OVERRIDE below; each asserts the computed category first):
  no_callable_phnJ          Albimonas.1 MePn (was lost_rule): no assembly has a callable phnJ (the 371-aa ORF
                            covers 36 % of the PhnJ HMM; fusion/frameshift?), pending a tblastn/frame check
  rast_label_not_supported  BM714.4 betaine (was lost_rule): KOfam's best KO puts the 428-aa "GrdI" in GrdE and
                            the Sec-joined "GrdH" in GrdB (GrdE/GrdB paralogs), so RAST's label is not supported
Near-miss scores: best score of the defining KO in the cluster's re-assembly / original MAGs, any score >= 14.7.

Writes results/comparison_vs_original/: clusters_guilds_old_vs_new.tsv, original_guild_members.tsv,
category_counts.tsv, taxa_of_interest.tsv.
"""
import re
from pathlib import Path

import pandas as pd

K = Path(__file__).resolve().parent.parent
RES = K / "results"
R26 = K.parents[1]
S = R26.parent
OUT = RES / "comparison_vs_original"
PAT = re.compile(r"Salt_Pond_MetaG_?(.+?)_MG_DASTool_bins[._](concoct_out|metabat|maxbin)\.(\d+)")
GUILDS = ["mepn_consumer", "betaine_reducer", "methanogen", "tma"]
DEF_KO = {"mepn_consumer": {"phnJ": "K06163"}, "betaine_reducer": {"grdH": "K21579", "grdI": "K21578"},
          "methanogen": {"mcrA": "K00399"}, "tma": {"cntA": "K22443", "cntB": "K22444", "cutC": "K20038", "cutD": "K20037"}}
# (species_cluster, guild) -> (computed category it must have, label to report); see the docstring
OVERRIDE = {
    ("old|Salt_Pond_MetaG_R2_restored_C_black_MG_DASTool_bins_metabat.8.contigs__", "mepn_consumer"):
        ("lost_rule", "no_callable_phnJ"),                       # Albimonas.1
    ("old|Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.51.contigs__", "betaine_reducer"):
        ("lost_rule", "rast_label_not_supported"),               # BM714.4 (cluster also holds BM714.1)
}
TAXA = [("Rhodobacteraceae", "family", "Rhodobacteraceae"), ("Roseovarius", "genus", "Roseovarius"),
        ("Marivita", "genus", "Marivita"), ("Albimonas", "genus", "Albimonas"), ("Halomonas", "genus", "Halomonas"),
        ("Vreelandella", "genus", "Vreelandella"), ("Methanolobus", "genus", "Methanolobus"),
        ("Methanohalophilus", "genus", "Methanohalophilus"), ("Halanaerobiaceae", "family", "Halanaerobiaceae"),
        ("Desulfobacterota", "phylum", "Desulfobacterota"), ("haloarchaea (c. Halobacteria)", "class", "Halobacteria"),
        ("Bacillota (all)", "phylum", "Bacillota"), ("M55B127 (f.)", "family", "M55B127")]


def key(n):
    m = PAT.search(str(n))
    return f"{m[1]}|{m[2]}|{int(m[3])}" if m else None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(R26 / "comparison/mags_old_vs_new.tsv", sep="\t")
    meta["name"] = meta.genome.str.split("|", n=1).str[1]
    meta["is_old"] = meta.genome.str.startswith("old|")
    # ---- RAST calls, joined by key() (assert 310 <-> 310)
    gr = pd.read_csv(S / "correlational/r232/guilds_r232.csv")
    gr["k"] = gr.mag.map(key)
    old = meta[meta.is_old].copy()
    old["k"] = old.name.map(key)
    assert gr.k.notna().all() and old.k.notna().all() and gr.k.is_unique and old.k.is_unique
    assert set(gr.k) == set(old.k) and len(gr) == len(old) == 310
    gr = gr.set_index("k")
    old = old.set_index("k")
    tot = {c: int(gr[c].sum()) for c in ("mepn_consumer", "betaine_reducer", "tma_producer", "methanogen", "mepn_2020", "betaine_2020")}
    print("RAST totals after join:", tot)
    ko = pd.read_csv(RES / "originals_310/guilds_310_kofam.tsv", sep="\t")
    ko["k"] = ko.mag.map(key)
    ko = ko.set_index("k").loc[old.index]
    o = pd.DataFrame(index=old.index)
    o["mag"] = old.name
    o["iterativeID"] = gr.loc[o.index, "iterativeID_r232"]
    o["species_cluster"] = old.species_cluster
    for c in ["phylum", "class", "family", "genus", "species", "completeness", "contamination"]:
        o[c] = old[c]
    o["rast_mepn_consumer"] = gr.mepn_consumer
    o["rast_betaine_reducer"] = gr.betaine_reducer
    o["rast_methanogen"] = gr.methanogen
    o["rast_tma"] = gr.tma_producer
    o["mepn_2020"] = gr.mepn_2020
    o["betaine_2020"] = gr.betaine_2020
    o["rast_ev_phnJ"] = gr.ev_phnJ
    o["rast_ev_grdH"] = gr.ev_grdH
    o["rast_ev_tma"] = gr.ev_cntA.fillna("") + gr.ev_cutC_cluster.fillna("").map(lambda s: (" | " + s) if s else "")
    o["kofam_mepn_consumer"] = ko.mepn_consumer.astype(bool)
    o["kofam_betaine_reducer"] = ko.betaine_reducer.astype(bool)
    o["kofam_betaine_reducer_evidence"] = ko.betaine_reducer_evidence
    o["kofam_methanogen"] = ko.methanogen.astype(bool)
    o["kofam_tma"] = (ko.tma_producer | ko.tma_candidate_one_sided).astype(bool)
    o["kofam_phnGHIJKLM_n_distinct"] = ko.phnGHIJKLM_n_distinct
    for c in ("rast_mepn_consumer", "rast_betaine_reducer", "rast_methanogen", "rast_tma", "mepn_2020", "betaine_2020"):
        o[c] = o[c].astype(bool)

    n = pd.read_csv(RES / "guilds_1324.tsv", sep="\t")
    n["tma"] = (n.tma_producer | n.tma_candidate_one_sided).astype(int)
    best = pd.read_csv(RES / "guilds_species_clusters_new.tsv", sep="\t").set_index("species_cluster").best_mag

    # ---- near-miss scores (best score of the defining KOs, any reported score)
    hn = pd.read_csv(RES / "all_marker_hits.tsv.gz", sep="\t", usecols=["mag", "ko", "score", "full_score"])
    ho = pd.read_csv(RES / "originals_310/marker_hits.tsv.gz", sep="\t", usecols=["mag", "ko", "score", "full_score"])
    allko = {k for d in DEF_KO.values() for k in d.values()}
    sc = lambda h: h[h.ko.isin(allko)].assign(s=lambda d: d.score.fillna(d.full_score)).groupby(["mag", "ko"]).s.max()
    sn, so = sc(hn), sc(ho)
    j = pd.read_csv(RES / "sec_readthrough/sec_readthrough_junctions_scored.tsv", sep="\t")
    jg = j.groupby("mag").joined_grdH.max()
    thr = pd.read_csv(RES / "hmm_info.tsv", sep="\t").set_index("ko").threshold

    # ---- cluster table
    clusters = sorted(set(meta.species_cluster))
    rows = []
    for c in clusters:
        oc = o[o.species_cluster == c]
        nc = n[n.species_cluster == c]
        r = dict(species_cluster=c, n_orig=len(oc), n_new=len(nc),
                 orig_ids=";".join(oc.iterativeID), best_new_mag=best.get(c, ""))
        src = oc.iloc[0] if len(oc) else n[n.mag == best[c]].iloc[0]
        for t in ("phylum", "class", "family", "genus", "species"):
            r[t] = src[t]
        for g in GUILDS:
            rast = bool(oc[f"rast_{g}"].any()) if len(oc) else False
            kof = bool(oc[f"kofam_{g}"].any()) if len(oc) else False
            nn = int(nc[g].sum()) if len(nc) else 0
            r[f"{g}_rast"] = int(rast)
            r[f"{g}_kofam_orig"] = int(kof)
            r[f"{g}_new_n"] = nn
            r[f"{g}_new_any"] = int(nn > 0)
            r[f"{g}_new_best"] = int(nc.loc[nc.mag == best.get(c, ""), g].sum()) if len(nc) else 0
            if rast and nn > 0:
                cat = "recovered_same"
            elif rast and len(nc) == 0:
                cat = "not_recovered"
            elif rast:
                cat = "lost_genome" if kof else "lost_rule"
            elif nn > 0:
                cat = "gained_new_cluster" if len(oc) == 0 else ("gained_rule" if kof else "gained_genome")
            elif kof:
                cat = "kofam_orig_only"
            else:
                cat = ""
            if (c, g) in OVERRIDE:
                expect, label = OVERRIDE[(c, g)]
                assert cat == expect, (c, g, cat, expect)
                cat = label
            r[f"{g}_category"] = cat
            if cat:
                parts = []
                for sym, k in DEF_KO[g].items():
                    vn = pd.Series([sn.get((m, k), float("nan")) for m in nc.mag], dtype=float).max()
                    vo = pd.Series([so.get((m, k), float("nan")) for m in oc.mag], dtype=float).max()
                    parts.append(f"{sym}(thr {thr[k]}): new {vn:.1f} / orig {vo:.1f}".replace("nan", "-"))
                if g == "betaine_reducer":
                    vj = pd.Series([jg.get(m, float("nan")) for m in nc.mag], dtype=float).max()
                    parts.append(f"grdH Sec-joined: new {vj:.1f}".replace("nan", "-"))
                r[f"{g}_best_scores"] = "; ".join(parts)
            else:
                r[f"{g}_best_scores"] = ""
        if True:                                                # 2020 manuscript sets, context only
            r["mepn_2020"] = int(oc.mepn_2020.any()) if len(oc) else 0
            r["betaine_2020"] = int(oc.betaine_2020.any()) if len(oc) else 0
        rows.append(r)
    ct = pd.DataFrame(rows)
    ct.to_csv(OUT / "clusters_guilds_old_vs_new.tsv", sep="\t", index=False)

    # ---- category counts
    cc = pd.concat({g: ct[f"{g}_category"].value_counts() for g in GUILDS}, axis=1).fillna(0).astype(int)
    cc = cc.drop(index="", errors="ignore")
    cc.to_csv(OUT / "category_counts.tsv", sep="\t")
    print(cc.to_string())

    # ---- original guild members (RAST or KOfam) with their cluster's re-assembly call
    mem = o[o[[f"rast_{g}" for g in GUILDS] + [f"kofam_{g}" for g in GUILDS] + ["mepn_2020", "betaine_2020"]].any(axis=1)].copy()
    cti = ct.set_index("species_cluster")
    for g in GUILDS:
        mem[f"{g}_new_n_of_N"] = [f"{cti.at[c, f'{g}_new_n']}/{cti.at[c, 'n_new']}" for c in mem.species_cluster]
        mem[f"{g}_category"] = [cti.at[c, f"{g}_category"] for c in mem.species_cluster]
    mem.sort_values("iterativeID").to_csv(OUT / "original_guild_members.tsv", sep="\t", index=False)

    # ---- taxa of interest
    trows = []
    for label, rank, val in TAXA:
        oc, nc = o[o[rank] == val], n[n[rank] == val]
        cl = ct[ct[rank] == val]
        r = dict(taxon=label, orig_mags=len(oc), orig_clusters=oc.species_cluster.nunique(),
                 new_mags=len(nc), new_clusters=nc.species_cluster.nunique(),
                 clusters_both=int(((cl.n_orig > 0) & (cl.n_new > 0)).sum()))
        for g in GUILDS:
            r[f"{g}: RAST orig MAGs"] = int(oc[f"rast_{g}"].sum())
            r[f"{g}: KOfam orig MAGs"] = int(oc[f"kofam_{g}"].sum())
            r[f"{g}: new MAGs"] = int(nc[g].sum())
            r[f"{g}: new clusters"] = int((cl[f"{g}_new_any"] == 1).sum())
            cats = cl[f"{g}_category"].value_counts()
            r[f"{g}: categories"] = ";".join(f"{k}:{v}" for k, v in cats.items() if k)
        r["mttB_family_without_mcrA: new MAGs"] = int(nc.mttB_family_without_mcrA.sum())
        r["methylamine_user: new MAGs"] = int(nc.methylamine_user.sum())
        trows.append(r)
    tt = pd.DataFrame(trows)
    tt.to_csv(OUT / "taxa_of_interest.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(tt.set_index("taxon").T.to_string())


if __name__ == "__main__":
    main()
