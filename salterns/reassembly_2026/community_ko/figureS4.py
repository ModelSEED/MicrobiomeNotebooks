#!/usr/bin/env python
"""Figure S4 (2026): Spearman correlations between metabolite concentrations (NMR) and methanogenesis gene abundances.

Samples: the 19 sediment samples with NMR data (all R1, R2 and R2A samples and SF2_A_D1).
Metabolites: the 40 of the 50 NMR metabolites present in at least three of the 19 samples.
Genes (rows, grouped by pathway as in the original figure): acetate activation and CODH/ACS, the H2/CO2 (C1) pathway,
the methyl-H4MPT:CoM methyltransferase, the methylamine, methanol, glycine betaine and methylthiol methyltransferases,
methyl-CoM reductase, and C-P lyase. Abundance: copies per genome equivalent in all reads (read survey named in
../read_survey_current.json). For the methanogenesis pathways only reads whose best hit is a methanogen lineage
(copies_per_GE_methanogen_lineage: Archaea of orders, families or genera named Methano*, ANME excluded) are counted, so
that bacterial and haloarchaeal homologs (acetate kinase, Pta, CODH, Mch, Mer, the MttB superfamily, ...) do not stand in
for methanogens; C-P lyase counts all reads (mostly bacterial).
Correlation: Spearman rho per gene × metabolite; Benjamini-Hochberg FDR over all tests; cells with q < 0.05 are marked.
Columns are clustered as pheatmap does by default (complete linkage, Euclidean distance between metabolite columns).
Writes figures/figureS4_2026.png, results/figureS4_rho.tsv (long: gene, metabolite, rho, p, q), results/figureS4_meta.json.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import heatmap_lib as H  # noqa: E402
import metabolome_taxonomy as MT  # noqa: E402
import profiles as PR  # noqa: E402
import stats_lib as SL  # noqa: E402

# (pathway label, genes) in row order; process and domain per pathway
PATHWAYS = [
    ("Acetate", ["ackA", "pta", "pta_K13788", "pta (K15024)", "eutD", "acs", "cdhA", "cdhB", "cdhC", "cdhD", "cdhE"],
     "Methyl-CoM formation", "Archaeal"),
    ("Hydrogen/CO2", ["fwdA", "fwdB", "fwdC", "fwdD", "fwdE", "fwdF", "fwdG", "fwdH", "ftr", "mch", "mer",
                      "mtrA", "mtrB", "mtrC", "mtrD", "mtrE", "mtrF", "mtrG", "mtrH"], "Methyl-CoM formation", "Archaeal"),
    ("Trimethylamine", ["mttB", "mttC"], "Methyl-CoM formation", "Archaeal"),
    ("Dimethylamine", ["mtbB", "mtbC"], "Methyl-CoM formation", "Archaeal"),
    ("Methylamine", ["mtmB", "mtmC"], "Methyl-CoM formation", "Archaeal"),
    ("Methylamines", ["mtbA"], "Methyl-CoM formation", "Archaeal"),
    ("Methanol", ["mtaB", "mtaC", "mtaA"], "Methyl-CoM formation", "Archaeal"),
    ("Glycine betaine", ["mtgB"], "Methyl-CoM formation", "Archaeal"),
    ("Methylthiols", ["mtsA"], "Methyl-CoM formation", "Archaeal"),
    ("Methyl-CoM reduction", ["mcrA", "mcrB", "mcrG"], "Methanogenesis", "Archaeal"),
    ("C-P lyase", ["phnG", "phnH", "phnI", "phnJ", "phnK", "phnL", "phnM"], "Methylphosphonate", "Bacterial"),
]
PATH_COL = {"Acetate": "#a6cee3", "Hydrogen/CO2": "#1f78b4", "Trimethylamine": "#b2df8a", "Dimethylamine": "#33a02c",
            "Methylamine": "#fb9a99", "Methylamines": "#e31a1c", "Methanol": "#fdbf6f", "Glycine betaine": "#b15928",
            "Methylthiols": "#ffff99", "Methyl-CoM reduction": "#ff7f00", "C-P lyase": "#cab2d6"}
PROC_COL = {"Methyl-CoM formation": "#440154", "Methanogenesis": "#21918c", "Methylphosphonate": "#fde725"}
DOM_COL = {"Archaeal": "#bebebe", "Bacterial": "#000000"}
DOM_LABEL = {"Archaeal": "Methanogen lineages", "Bacterial": "All reads"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=None)
    ap.add_argument("--out", default=os.path.join(HERE, "figures", "figureS4_2026.png"))
    a = ap.parse_args()
    X = MT.metabolome(min_samples=3)
    samples = list(X.index)
    P = PR.panel().set_index("label")
    t = PR.read_profile(a.results)
    have_arch = "copies_per_GE_methanogen_lineage" in t.columns
    W = PR.wide(t, "copies_per_GE", samples)
    WA = PR.wide(t, "copies_per_GE_methanogen_lineage", samples) if have_arch else None
    rows, ann, missing = [], [], []
    for pw, genes, proc, dom in PATHWAYS:
        for g in genes:
            ko = P.loc[g, "ko"]
            src = WA if (dom == "Archaeal" and have_arch) else W
            if ko not in src.index or src.loc[ko].isna().any():
                missing.append(g); continue
            if float(src.loc[ko].astype(float).sum()) == 0:           # e.g. no methanogen-lineage reads at all
                missing.append(g + " (no counted reads)"); continue
            rows.append((g, ko, src.loc[ko].astype(float)))
            ann.append((pw, proc, dom))
    recs = []
    R = np.full((len(rows), X.shape[1]), np.nan)
    for i, (g, ko, v) in enumerate(rows):
        for j, m in enumerate(X.columns):
            if v.nunique() < 2:
                continue
            rho, p = spearmanr(v.values, X[m].values)
            R[i, j] = rho
            recs.append(dict(gene=g, ko=ko, metabolite=m, rho=rho, p=p, pathway=ann[i][0], domain=ann[i][2],
                             methanogen_lineage_only=bool(ann[i][2] == "Archaeal" and have_arch)))
    L = pd.DataFrame(recs)
    L["q"] = SL.bh(L.p.values)
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    L.to_csv(os.path.join(HERE, "results", "figureS4_rho.tsv"), sep="\t", index=False)
    Rf = np.nan_to_num(R, nan=0.0)
    Zc = linkage(Rf.T, method="complete", metric="euclidean")
    order = leaves_list(Zc)
    mets = [X.columns[j] for j in order]
    Q = L.pivot(index="gene", columns="metabolite", values="q")
    mc = L[L.gene == "mcrA"].sort_values("rho")
    meta = dict(survey=PR.survey_dir(a.results), methanogen_lineage_only=have_arch, n_samples=len(samples), samples=samples,
                n_metabolites=X.shape[1], n_genes=len(rows), missing_genes=missing, n_tests=int(len(L)),
                n_q05=int((L.q < 0.05).sum()), metabolite_order=mets,
                mcrA_pos_q05=mc[(mc.q < 0.05) & (mc.rho > 0)].sort_values("rho", ascending=False)[["metabolite", "rho", "q"]].values.tolist(),
                mcrA_neg_q05=mc[(mc.q < 0.05) & (mc.rho < 0)][["metabolite", "rho", "q"]].values.tolist(),
                mcrA_top_pos=mc.tail(6)[::-1][["metabolite", "rho", "q"]].values.tolist(),
                mcrA_top_neg=mc.head(6)[["metabolite", "rho", "q"]].values.tolist())
    json.dump(meta, open(os.path.join(HERE, "results", "figureS4_meta.json"), "w"), indent=1, default=float)
    print(f"{len(rows)} genes x {X.shape[1]} metabolites, {len(samples)} samples; methanogen lineages: {have_arch}; "
          f"q < 0.05: {meta['n_q05']} of {meta['n_tests']}; missing: {missing}")
    print("mcrA +:", [(m, round(r, 2)) for m, r, q in meta["mcrA_top_pos"]])
    print("mcrA -:", [(m, round(r, 2)) for m, r, q in meta["mcrA_top_neg"]])

    # ---- draw (pheatmap layout: dendrogram over the columns, annotation bars at the left, labels at the right)
    n, m = len(rows), len(mets)
    rh, cw = 0.165, 0.2
    ann_w = 0.17
    x0 = 0.3
    hx = x0 + 3 * (ann_w + 0.03) + 0.06
    dend_h = 1.3
    top = 0.15 + dend_h + 0.05
    W_in = hx + cw * m + 1.45 + 2.6
    H_in = top + rh * n + 2.0
    c = H.Canvas(W_in, H_in)
    ad = c.ax(hx, 0.15, cw * m, dend_h)
    H.draw_dendro(ad, Zc, m)
    ah = c.ax(hx, top, cw * m, rh * n)
    M = R[:, order]
    ah.imshow(M, aspect="auto", cmap=H.BWR, vmin=-0.8, vmax=0.8, interpolation="nearest")
    Qm = np.array([[Q.loc[g, mt] if (g in Q.index and mt in Q.columns) else np.nan for mt in mets] for g, _, _ in rows])
    yy, xx = np.where(Qm < 0.05)
    ah.scatter(xx, yy, s=3.5, color="#000000", marker="o", lw=0)
    H.blank(ah)
    for j, (nm, cols) in enumerate([("Reads", [DOM_COL[x[2]] for x in ann]), ("Process", [PROC_COL[x[1]] for x in ann]),
                                    ("Pathway", [PATH_COL[x[0]] for x in ann])]):
        ax_ = c.ax(x0 + j * (ann_w + 0.03), top, ann_w, rh * n)
        H.color_bar(ax_, cols, vertical=True)
        c.fig.text((x0 + j * (ann_w + 0.03) + ann_w / 2) / W_in, 1 - (top + rh * n + 0.06) / H_in, nm, rotation=-50,
                   ha="left", va="top", fontsize=9.5, weight="bold", rotation_mode="anchor")
    for i, (g, ko, _) in enumerate(rows):
        gl = g.split(" ")[0].split("_")[0]
        c.fig.text((hx + cw * m + 0.06) / W_in, 1 - (top + rh * (i + 0.5)) / H_in, f"{ko} {gl}", fontsize=8.3,
                   va="center")
    for j, mt in enumerate(mets):
        c.fig.text((hx + cw * (j + 0.5)) / W_in, 1 - (top + rh * n + 0.06) / H_in, mt, rotation=-50, ha="left",
                   va="top", fontsize=9, rotation_mode="anchor")
    lx = hx + cw * m + 1.45
    y = H.colorbar(c, lx, 0.25, 1.6, H.BWR, -0.8, 0.8, [-0.5, 0, 0.5], title="Spearman ρ")
    c.fig.text(lx / W_in, 1 - (y - 0.05) / H_in, "•  FDR < 0.05", fontsize=8.5, va="top")
    y += 0.3
    y = H.legend_block(c, lx, y, "Reads counted", [(DOM_LABEL[k], v) for k, v in DOM_COL.items()])
    y = H.legend_block(c, lx, y, "Process", list(PROC_COL.items()))
    y = H.legend_block(c, lx, y, "Pathway", [(k, v) for k, v in PATH_COL.items() if k in {x[0] for x in ann}])
    c.fig.savefig(a.out, dpi=300, facecolor="white")
    print("wrote", a.out, f"{W_in:.1f} x {H_in:.1f} in")


if __name__ == "__main__":
    main()
