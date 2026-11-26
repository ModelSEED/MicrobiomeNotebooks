#!/usr/bin/env python
"""Figure 1 (2026): relative abundance of the methylated-compound genes in the 24 sediment metagenomes.

Genes: the 108 KOs of ko_panel_methylated.tsv outside the Figure S4 set (methanogenesis core and acetate genes), kept
when detected (>= 1 read) in at least one of the 24 sediment metagenomes, as in the original Figure 1.
Abundance: gene copies per genome equivalent (GE) in all reads (read survey named in ../read_survey_current.json, or
--results DIR), z-scored per gene across the 24 samples (clipped at ±4).
Columns: Ward clustering (Ward.D2 criterion, Euclidean distance on the z-scores). Rows: grouped by compound.
Statistics (panel_stats.py): Kruskal-Wallis across the three restoration types (unrestored R1 + R2, n = 12; restored
SF2, n = 6; reference R2A, n = 6), Benjamini-Hochberg FDR over all detected panel KOs (shared with Table S1), Dunn's
test with Holm adjustment.
Writes figures/figure1_2026.png, results/figure1_stats.tsv, results/figure1_meta.json.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage, leaves_list

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import heatmap_lib as H  # noqa: E402
import profiles as PR  # noqa: E402
import stats_lib as SL  # noqa: E402

TYPES = ["unrestored", "restored", "reference"]
REST_COL = {"reference": "#440154", "restored": "#21918c", "unrestored": "#fde725"}
CLASS_COL = {"TMA and TMAO": "#E69F00", "DMA": "#56B4E9", "MMA": "#009E73", "Choline and carnitine": "#F0E442",
             "GB, DMG and sarcosine": "#0072B2", "Methanol": "#D55E00", "Formaldehyde and C1": "#CC79A7",
             "DMS and sulfonates": "#9e9e9e", "DMSP and MMPA": "#8C564B", "Phosphonates": "#17BECF",
             "Methane": "#1a1a1a", "Acetate": "#a6cee3", "H2/CO2": "#1f78b4"}
TYPE_COL = {"degradation": "#d9d9d9", "TMA production": "#e7298a", "biosynthesis / production": "#7570b3",
            "transport": "#66a61e", "methyl transfer (corrinoid)": "#1b9e77", "methanogenesis (Mcr)": "#1a1a1a",
            "acetate / CODH": "#a6cee3",
            "methanogenesis (C1 core)": "#1f78b4"}
# the two side-by-side blocks (by compound class) used when there are too many rows for one column
BLOCK_A = ["TMA and TMAO", "DMA", "MMA", "Choline and carnitine", "GB, DMG and sarcosine", "Methane"]
VAL_COL = {"validated": "#ffffff", "flagged": "#e08214", "not testable": "#bababa"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=None)
    ap.add_argument("--source", choices=["read", "assembly"], default="read",
                    help="read: copies per GE in all reads (Figure 1); assembly: KOfam on the assembled genes, depth-"
                         "weighted (concordance check, Table S1)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    tag = "" if a.source == "read" else "_assembly"
    a.out = a.out or os.path.join(HERE, "figures", f"figure1_2026{tag}.png")
    P = PR.panel()
    P1 = P[P.status != "FigS4"].set_index("ko")
    t = PR.profile(a.source, a.results)
    W = PR.wide(t, "copies_per_GE", PR.SED, kos=list(P1.index))
    R = PR.wide(t, "reads", PR.SED, kos=list(P1.index))
    surveyed = W.notna().all(axis=1)
    detected = surveyed & (R.fillna(0).sum(axis=1) > 0)
    kos = [k for k in P1.index if detected[k]]
    Wd = W.loc[kos].astype(float)
    sd = Wd.std(axis=1, ddof=1)
    assert (sd > 0).all()
    Zs = Wd.sub(Wd.mean(axis=1), axis=0).div(sd, axis=0)
    Zl = linkage(Zs.T.values, method="ward")
    order = leaves_list(Zl)
    cols = [Zs.columns[i] for i in order]
    k3 = fcluster(Zl, 3, criterion="maxclust")
    lab = pd.Series(k3, index=Zs.columns)
    typ = pd.Series([SL.site_type(s) for s in Zs.columns], index=Zs.columns)
    ct = pd.crosstab(lab, typ)
    purity = float(ct.max(axis=1).sum() / ct.values.sum())
    exact = bool((ct > 0).sum(axis=1).eq(1).all() and (ct > 0).sum(axis=0).eq(1).all())
    import panel_stats as PS                   # one FDR family over all detected panel KOs (shared with Table S1)
    st = PS.compute(a.source, a.results).loc[kos]
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    out = P1.loc[kos, ["compound_class", "ptype"]].join(st)
    if a.source == "read":
        v_, m_ = PR.read_validation()
        out["read_validation"] = [PR.validation_status(P1.loc[k, "label"], v_, m_) for k in kos]
    out["reads_24"] = R.loc[kos].sum(axis=1).astype(int)
    out["n_samples_detected"] = (R.loc[kos] > 0).sum(axis=1).astype(int)
    out.to_csv(os.path.join(HERE, "results", f"figure1{tag}_stats.tsv"), sep="\t")
    meta = dict(source=a.source, survey=PR.survey_dir(a.results) if a.source == "read" else "assembly_based", n_panel=len(P1), n_surveyed=int(surveyed.sum()),
                n_detected=len(kos), not_detected=[k for k in P1.index if surveyed[k] and not detected[k]],
                not_surveyed=[k for k in P1.index if not surveyed[k]],
                n_q05=int((st.q < 0.05).sum()), ward_k3_crosstab=ct.to_dict(), ward_k3_purity=purity,
                ward_k3_matches_restoration=exact, column_order=cols)
    json.dump(meta, open(os.path.join(HERE, "results", f"figure1{tag}_meta.json"), "w"), indent=1, default=str)
    print(f"{len(kos)} genes plotted ({meta['n_surveyed']} surveyed); KW q < 0.05: {meta['n_q05']}; "
          f"Ward k=3 purity {purity:.2f}, exact {exact}")
    print(ct)

    # ---- draw
    classes = [c for c, _ in PR.CLASS if c in set(P1.loc[kos, "compound_class"])]
    def rows_of(cls):
        return [k for c in cls for k in kos if P1.loc[k, "compound_class"] == c]
    if len(kos) > 60:
        blocks = [rows_of([c for c in classes if c in BLOCK_A]), rows_of([c for c in classes if c not in BLOCK_A])]
    else:
        blocks = [rows_of(classes)]
    rh = 0.118                                     # row height (in)
    cw = 0.155                                     # column width
    nmax = max(len(b) for b in blocks)
    ann_w, gap_ann, lab_w, blk_gap = 0.13, 0.05, 1.05, 0.25
    vdat = PR.read_validation() if a.source == "read" else (None, {})
    vstat = {k: PR.validation_status(P1.loc[k, "label"], *vdat) for k in kos}
    n_ann = 5 if a.source == "read" else 4
    left0 = 0.45
    blk_w = n_ann * (ann_w + 0.02) + gap_ann + cw * 24 + lab_w
    dend_h, bar_h = 0.9, 0.16
    top = 0.15 + dend_h + 0.05 + bar_h + 0.08
    hm_h = rh * nmax
    W_in = left0 + len(blocks) * blk_w + (len(blocks) - 1) * blk_gap + 2.3
    n_types = len(set(P1.loc[kos, "ptype"]))
    leg_h = 0.25 + 0.24 + 1.3 + 0.3 + sum(0.24 + n * 0.19 + 0.12 for n in (3, 2, len(classes), n_types, 3))
    H_in = max(top + hm_h + 1.0, leg_h + 0.3)
    c = H.Canvas(W_in, H_in)
    vmin, vmax = -4, 4
    for bi, rows in enumerate(blocks):
        x0 = left0 + bi * (blk_w + blk_gap)
        hx = x0 + n_ann * (ann_w + 0.02) + gap_ann
        ad = c.ax(hx, 0.15, cw * 24, dend_h)
        H.draw_dendro(ad, Zl, 24)
        ab = c.ax(hx, 0.15 + dend_h + 0.05, cw * 24, bar_h)
        H.color_bar(ab, [REST_COL[SL.site_type(s)] for s in cols], vertical=False)
        if bi == len(blocks) - 1:
            c.fig.text((hx + cw * 24 + 0.08) / W_in, 1 - (0.15 + dend_h + 0.05 + bar_h / 2) / H_in, "Restoration",
                       fontsize=9.5, weight="bold", va="center")
        M = Zs.loc[rows, cols].clip(vmin, vmax).values
        ah = c.ax(hx, top, cw * 24, rh * len(rows))
        ah.imshow(M, aspect="auto", cmap=H.BWR, vmin=vmin, vmax=vmax, interpolation="nearest")
        H.blank(ah)
        anns = [("KW", lambda k: "#000000" if st.loc[k, "q"] < 0.05 else "#bebebe"),
                ("Highest", lambda k: REST_COL[st.loc[k, "highest"]] if st.loc[k, "q"] < 0.05 else "#ffffff"),
                ("Compound", lambda k: CLASS_COL[P1.loc[k, "compound_class"]]),
                ("Type", lambda k: TYPE_COL[P1.loc[k, "ptype"]])]
        if n_ann == 5:
            anns.append(("Validation", lambda k: VAL_COL[vstat[k]]))
        for j, (nm, colf) in enumerate(anns):
            ax_ = c.ax(x0 + j * (ann_w + 0.02), top, ann_w, rh * len(rows))
            H.color_bar(ax_, [colf(k) for k in rows], vertical=True)
            if nm in ("Highest", "Validation"):
                ax_.add_patch(H.Rectangle((0, 0), 1, len(rows), facecolor="none", edgecolor="#bebebe", lw=0.4))
            c.fig.text((x0 + j * (ann_w + 0.02) + ann_w / 2) / W_in, 1 - (top + rh * len(rows) + 0.06) / H_in, nm,
                       rotation=55, ha="right", va="top", fontsize=8.5, weight="bold", rotation_mode="anchor")
        for i, k in enumerate(rows):
            c.fig.text((hx + cw * 24 + 0.05) / W_in, 1 - (top + rh * (i + 0.5)) / H_in,
                       f"{k} {P1.loc[k, 'gene']}", fontsize=6.6, va="center")
        for j, s in enumerate(cols):
            c.fig.text((hx + cw * (j + 0.5)) / W_in, 1 - (top + rh * len(rows) + 0.06) / H_in, s, rotation=55,
                       ha="right", va="top", fontsize=7.2, rotation_mode="anchor")
    lx = left0 + len(blocks) * blk_w + (len(blocks) - 1) * blk_gap + 0.15
    y = H.colorbar(c, lx, 0.25, 1.3, H.BWR, vmin, vmax, [-4, -2, 0, 2, 4], title="z score")
    y = H.legend_block(c, lx, y, "Restoration", [("Reference", REST_COL["reference"]),
                                                 ("Restored", REST_COL["restored"]),
                                                 ("Unrestored", REST_COL["unrestored"])])
    y = H.legend_block(c, lx, y, "KW", [("FDR < 0.05", "#000000"), ("FDR ≥ 0.05", "#bebebe")])
    y = H.legend_block(c, lx, y, "Compound", [(k, CLASS_COL[k]) for k in classes])
    y = H.legend_block(c, lx, y, "Type", [(k, TYPE_COL[k]) for k in TYPE_COL
                                          if k in set(P1.loc[kos, "ptype"])])
    if n_ann == 5:
        y = H.legend_block(c, lx, y, "Validation", [("validated", "#ffffff"), ("flagged or merged", VAL_COL["flagged"]),
                                                     ("not testable", VAL_COL["not testable"])])
    c.fig.savefig(a.out, dpi=300, facecolor="white")
    print("wrote", a.out, f"{W_in:.1f} x {H_in:.1f} in")


if __name__ == "__main__":
    main()
