#!/usr/bin/env python
"""Guild calls of the FINAL salterns MAG set (2026-11-19): the 935 Arm C rollout MAGs + the 520 kept MAGs of the four
first-pass MEGAHIT site co-assemblies (rollout QC), species clusters from compare_rollout.py (skani 95 % ANI / 50 % AF,
joint clustering; comparison_rollout/mags_three_way.tsv).

KOfam calls (KOfam 2026-10-26, DRAM 1.5 best-KO rule, Sec read-through rescue, revised guild rules), per MAG:
  rollout     annotation_runs/kofam_rollout/results/guilds_935.tsv
  coassembly  the 482 co-assembly MAGs that the first pass had KBase-selected: annotation_runs/kofam/results/guilds_1324.tsv
              (first-pass name '<binner>.bin.N.fasta__<unit>' == co-assembly '<unit>__<binner>.bin.N'); the other 38:
              annotation_runs/kofam_coasm38/results/guilds_38.tsv once it exists (else they count as 'missing').
Originals (310): RAST guild calls of correlational/r232/guilds_r232.csv, joined by the key() of
annotation_runs/kofam/scripts/compare_vs_original.py.

Writes comparison_rollout/final_guilds.tsv (per MAG), final_guilds_species.tsv (per species cluster, any-MAG calls),
final_guilds_summary.json (numbers for the manuscript builder). Run: ~/Documents/py_venv/bin/python final_set_guilds.py
"""
import json, os, re

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
S = os.path.dirname(HERE)
OUT = os.path.join(HERE, "comparison_rollout")
G = ["methanogen", "mepn_consumer", "betaine_reducer", "methylamine_user", "tma_producer"]
EXTRA = ["tma_candidate_one_sided", "mttB_family_without_mcrA", "n_grdH", "n_grdH_prodigal_only",
         "n_grdH_sec_readthrough", "n_cutC", "n_cutD", "n_cntA", "n_cntB", "n_mcrA", "n_phnJ", "n_grdI",
         "betaine_reducer_evidence"]
PAT = re.compile(r"Salt_Pond_MetaG_?(.+?)_MG_DASTool_bins[._](concoct_out|metabat|maxbin)\.(\d+)")


def key(n):
    m = PAT.search(str(n))
    return f"{m[1]}|{m[2]}|{int(m[3])}" if m else None


def to_coasm(n):   # 'concoct.bin.002.fasta__R1_sediment' -> 'R1_sediment__concoct.bin.002'
    b, u = n.split("__")
    return f"{u}__{re.sub(r'.fasta$', '', b)}"


def main():
    d = pd.read_csv(os.path.join(OUT, "mags_three_way.tsv"), sep="\t")
    fin = d[d.set.isin(["rollout", "coassembly"])].copy()
    ro = pd.read_csv(os.path.join(HERE, "annotation_runs/kofam_rollout/results/guilds_935.tsv"), sep="\t")
    ro["genome"], ro["kofam_source"] = "rollout|" + ro.mag, "kofam_rollout"
    fp = pd.read_csv(os.path.join(HERE, "annotation_runs/kofam/results/guilds_1324.tsv"), sep="\t")
    co = fp[fp.mag.str.endswith("_sediment")].copy()
    co["genome"], co["kofam_source"] = "coassembly|" + co.mag.map(to_coasm), "kofam (first pass)"
    parts = [ro, co]
    f38 = os.path.join(HERE, "annotation_runs/kofam_coasm38/results/guilds_38.tsv")
    if os.path.exists(f38):
        c38 = pd.read_csv(f38, sep="\t")
        c38["genome"], c38["kofam_source"] = "coassembly|" + c38.mag, "kofam_coasm38"
        parts.append(c38)
    calls = pd.concat([p[["genome", "kofam_source"] + G + EXTRA] for p in parts], ignore_index=True)
    calls = calls[calls.genome.isin(set(fin.genome))]          # first-pass co-assembly MAGs dropped by the rollout QC
    assert calls.genome.is_unique
    m = fin[["genome", "set", "sample", "species_cluster", "final_species_rep", "completeness", "contamination",
             "mimag", "domain", "phylum", "class", "family", "genus", "species"]].merge(calls, on="genome", how="left")
    m["kofam_source"] = m.kofam_source.fillna("missing")
    has = m.kofam_source != "missing"
    for c in G + ["tma_candidate_one_sided", "mttB_family_without_mcrA"]:
        m[c] = m[c].astype("boolean")
    m.to_csv(os.path.join(OUT, "final_guilds.tsv"), sep="\t", index=False)
    x = m[has]
    sp = x.groupby("species_cluster")[G + ["tma_candidate_one_sided", "mttB_family_without_mcrA"]].any()
    tax = m.drop_duplicates("species_cluster").set_index("species_cluster")[["phylum", "class", "family", "genus"]]
    rep = m[m.final_species_rep].set_index("species_cluster")
    for c in ("phylum", "class", "family", "genus"):           # taxonomy of the representative
        tax[c] = rep[c].reindex(tax.index).fillna(tax[c])
    sp = sp.join(tax)
    sp.to_csv(os.path.join(OUT, "final_guilds_species.tsv"), sep="\t")

    def vc(s):
        return {k: int(v) for k, v in s.value_counts().items()}
    mep, bet = sp[sp.mepn_consumer], sp[sp.betaine_reducer]
    # RAST guilds of the originals, in the joint clustering
    gr = pd.read_csv(os.path.join(S, "correlational/r232/guilds_r232.csv"))
    gr["k"] = gr.mag.map(key)
    orig = d[d.set == "original"].copy()
    orig["k"] = orig.genome.str.split("|", n=1).str[1].map(key)
    assert gr.k.is_unique and orig.k.is_unique and set(gr.k) == set(orig.k) and len(gr) == 310
    orig = orig.merge(gr[["k", "mepn_consumer", "betaine_reducer", "methanogen", "iterativeID_r232"]]
                      .rename(columns={"mepn_consumer": "rast_mepn", "betaine_reducer": "rast_betaine",
                                       "methanogen": "rast_methanogen"}), on="k")
    oc = set(orig.species_cluster)

    def rast_vs_final(rcol, gcol):
        r = orig[orig[rcol]]
        cl = set(r.species_cluster)
        rec = {c for c in cl if c in sp.index and sp.loc[c, gcol]}
        return {"rast_mags": int(len(r)), "rast_species": int(len(cl)), "recovered_with_call": int(len(rec)),
                "rast_species_absent_from_final": int(len(cl - set(sp.index))),
                "rast_species_present_without_call": int(len({c for c in cl if c in sp.index} - rec)),
                "further_final_species_without_original": int(sum(1 for c in sp.index[sp[gcol]] if c not in oc)),
                "further_final_species_with_original_not_rast": int(sum(1 for c in sp.index[sp[gcol]]
                                                                      if c in oc and c not in cl))}
    alb = m[m.genus == "Albimonas"]
    bm4 = orig[orig.iterativeID_r232 == "BM714.4"]
    bm4_cl = bm4.species_cluster.iloc[0] if len(bm4) else None
    rh = x[x.family == "Rhodobacteraceae"]
    summ = {
        "n_mags": int(len(m)), "n_mags_with_calls": int(has.sum()), "n_mags_missing_calls": int((~has).sum()),
        "kofam_sources": vc(m.kofam_source), "n_species": int(m.species_cluster.nunique()),
        "guilds": {g: {"mags": int(x[g].sum()), "species": int(sp[g].sum())} for g in G},
        "mepn_species_rhodobacteraceae": int((mep.family == "Rhodobacteraceae").sum()),
        "mepn_species_haloarchaea": int((mep["class"] == "Halobacteria").sum()),
        "mepn_species_by_phylum": vc(mep.phylum),
        "betaine_species_by_phylum": vc(bet.phylum),
        "grdH_genes_prodigal_only": int(x.n_grdH_prodigal_only.sum()), "grdH_genes_after_rescue": int(x.n_grdH.sum()),
        "methanogen_genera": vc(x[x.methanogen].genus),
        "tma_cutC_without_cutD_mags": int(((x.n_cutC > 0) & (x.n_cutD == 0)).sum()),
        "tma_cntA_without_cntB_mags": int(((x.n_cntA > 0) & (x.n_cntB == 0)).sum()),
        "mttB_family_without_mcrA_mags": int(x.mttB_family_without_mcrA.sum()),
        "rhodobacteraceae_mags": int(len(rh)), "rhodobacteraceae_with_mttB_family": int(rh.mttB_family_without_mcrA.sum()),
        "rast_mepn": rast_vs_final("rast_mepn", "mepn_consumer"),
        "rast_betaine": rast_vs_final("rast_betaine", "betaine_reducer"),
        "albimonas": {"mags": int(len(alb)), "mags_with_calls": int((alb.kofam_source != "missing").sum()),
                      "mepn_mags": int(alb.mepn_consumer.fillna(False).sum())},
        "bm714_4_cluster": {"in_final": bool(bm4_cl in sp.index) if bm4_cl else False,
                            "betaine_call": bool(sp.loc[bm4_cl, "betaine_reducer"]) if bm4_cl in sp.index else None},
    }
    json.dump(summ, open(os.path.join(OUT, "final_guilds_summary.json"), "w"), indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
