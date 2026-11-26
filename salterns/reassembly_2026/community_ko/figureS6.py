#!/usr/bin/env python
"""Figure S6 (2026): relative abundance of the cultured genera in the 24 sediment samples.

Rows
  16S            the top row of the original Figure S6 (iTagger / Greengenes), cropped from the original image
                 (figures/figureS6_16S_row_from_original.png); the 16S data are not in this repository.
  Metagenome     % of the reads of 8 single-copy ribosomal-protein marker families (rpsB, rpsC, rpsE, rpsG, rplB, rplC,
  (reads)        rplE, rplF) whose best DIAMOND hit is a reference of the genus (read survey v2,
                 read_markers_v2/results/marker_genus_profile.tsv).
  MAGs           % of all quality-controlled reads that mapped (CoverM) to the genus's species in the final 585-species
                 MAG set (abundance_final/final585_relabund.tsv).
Groups: R1, R2 (unrestored), SF2 (restored), R2A (reference), six samples each (three cores × two depths).
Letters: Kruskal-Wallis across the four sites, then Dunn's test (Holm); sites sharing a letter do not differ (p >= 0.05).
Genus names follow the isolates (Table S3): Halomonas sensu lato includes GTDB Vreelandella (Halomonas species
reassigned in GTDB) and Clostridium includes the GTDB Clostridium_* genera.
Writes figures/figureS6_2026.png and results/figureS6_values.tsv, results/figureS6_stats.tsv.
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
R26 = os.path.dirname(HERE)
SAL = os.path.dirname(R26)
sys.path.insert(0, HERE)
import stats_lib as SL  # noqa: E402

GENERA = ["Halomonas", "Idiomarina", "Salegentibacter", "Clostridium", "Marinilabilia"]
SITES = ["R1", "R2", "SF2", "R2A"]
COL = {"R1": "#fde725", "R2": "#fde725", "SF2": "#21918c", "R2A": "#440154"}
SED = [f"{s}_{c}_{d}" for s in SITES for c in "ABC" for d in ("D1", "D2")]
OUT_FIG = os.path.join(HERE, "figures", "figureS6_2026.png")
OUT_DIR = os.path.join(HERE, "results")


def genus_group(g):
    """Map a GTDB / NCBI genus name to the cultured-genus label used in the figure (None if not one of them)."""
    if not isinstance(g, str):
        return None
    base = g.split("_")[0]
    if base in ("Halomonas", "Vreelandella"):
        return "Halomonas"
    if base in ("Saligentibacter", "Salegentibacter"):
        return "Salegentibacter"
    return base if base in GENERA else None


def mag_row():
    ids = pd.read_csv(os.path.join(SAL, "correlational", "final585", "iterativeIDs_final585.csv"))
    ids["genus"] = ids.gtdb_r232.str.extract(r"g__([^;]*)")[0]
    ids["fig"] = ids.genus.map(genus_group)
    ra = pd.read_csv(os.path.join(R26, "abundance_final", "final585_relabund.tsv"), sep="\t", index_col=0)
    rows = []
    for g in GENERA:
        mags = ids.loc[ids.fig == g, "mag"].tolist()
        v = ra.loc[mags, SED].sum() * 100 if mags else pd.Series(0.0, index=SED)
        for s in SED:
            rows.append(dict(row="MAGs", genus=g, sample=s, site=SL.site_of(s), pct=float(v[s]),
                             n_species=len(mags), members=";".join(ids.loc[ids.fig == g, "iterativeID"])))
    return pd.DataFrame(rows)


def read_row(path):
    t = pd.read_csv(path, sep="\t")
    gcol = next(c for c in ("genus_harmonized", "genus", "genus_gtdb") if c in t.columns)
    scol = next(c for c in ("share_of_marker_reads", "marker_share", "share") if c in t.columns)
    t["fig"] = t[gcol].map(genus_group)
    t = t[t["sample"].isin(SED)]
    piv = t.dropna(subset=["fig"]).groupby(["fig", "sample"])[scol].sum().unstack(fill_value=0.0)
    piv = piv.reindex(index=GENERA, columns=SED, fill_value=0.0) * 100
    rows = []
    for g in GENERA:
        for s in SED:
            rows.append(dict(row="Metagenome (reads)", genus=g, sample=s, site=SL.site_of(s), pct=float(piv.loc[g, s]),
                             n_species=np.nan, members=""))
    return pd.DataFrame(rows)


def read_sensitivity(path):
    """Per genus: share of its marker reads whose best hit is >= 95% identical, and the KW/Dunn result when only those
    reads are counted (denominator: all marker reads of the sample)."""
    t = pd.read_csv(path, sep="\t")
    t = t[t["sample"].isin(SED)]
    t["fig"] = t.genus_harmonized.map(genus_group)
    tot = t.groupby("sample").marker_reads.sum().reindex(SED)
    rows = []
    for g in GENERA:
        d = t[t.fig == g]
        v = (d.groupby("sample").reads_best_hit_ge95.sum().reindex(SED, fill_value=0) / tot * 100)
        r = SL.kw_dunn(v.values, [SL.site_of(s) for s in SED])
        rows.append(dict(genus=g, marker_reads=int(d.marker_reads.sum()),
                         frac_best_hit_ge95=float(d.reads_best_hit_ge95.sum() / max(d.marker_reads.sum(), 1)),
                         gtdb_genera=";".join(sorted(d.genus_harmonized.unique())), H_ge95=r["H"], p_ge95=r["p"],
                         **{f"letter_ge95_{s}": r["letters"][s] for s in SITES}))
    return pd.DataFrame(rows)


def stats_for(df):
    out = []
    for (r, g), d in df.groupby(["row", "genus"], sort=False):
        d = d.set_index("sample").loc[SED]
        res = SL.kw_dunn(d.pct.values, d.site.values)
        out.append(dict(row=r, genus=g, H=res["H"], p=res["p"],
                        **{f"median_{s}": res["medians"][s] for s in SITES},
                        **{f"mean_{s}": float(d.pct[d.site == s].mean()) for s in SITES},
                        **{f"letter_{s}": res["letters"][s] for s in SITES}))
    return pd.DataFrame(out)


ORIG = os.path.join(HERE, "figures", "figureS6_original_2021.png")   # the original Figure S6 (2048 x 1138 px)
# pixel geometry of the original (ggplot facet_grid): panel x ranges, right-strip x range, panel height, row gap
PX = [(104, 398), (417, 710), (729, 1022), (1076, 1369), (1388, 1682)]
STRIP_X = (1683, 1751)
TOP_STRIP = (18, 80)                 # the genus strips above the first row
ROW0 = (80, 543)                     # the 16S panel row in the original
ROW_H, GAP = 463, 19
FONT = "Liberation Sans"


def nice_top(v):
    from matplotlib.ticker import MaxNLocator
    t = MaxNLocator(nbins=3, steps=[1, 2, 2.5, 3, 5, 10]).tick_values(0, max(v, 1e-3))
    return float(t[t >= v * 1.04][0]) if (t >= v * 1.04).any() else float(t[-1])


def draw_row(fig, W, H, y0, df, st, strip_label, last):
    T = nice_top(df.pct.max())
    rng = np.random.default_rng(7)
    for k, g in enumerate(GENERA):
        x0, x1 = PX[k]
        ax = fig.add_axes([x0 / W, 1 - (y0 + ROW_H) / H, (x1 - x0) / W, ROW_H / H])
        d = df[df.genus == g]
        s_ = st[st.genus == g].iloc[0]
        for i, s in enumerate(SITES):
            v = d.loc[d.site == s, "pct"].values
            ax.boxplot(v, positions=[i], widths=0.78, showfliers=False, patch_artist=True,
                       boxprops=dict(facecolor="none", edgecolor=COL[s], lw=2.2),
                       medianprops=dict(color=COL[s], lw=3.0), whiskerprops=dict(color=COL[s], lw=2.2),
                       capprops=dict(lw=0))
            ax.scatter(i + rng.uniform(-0.25, 0.25, len(v)), v, s=250, color=COL[s], alpha=0.72, lw=0, zorder=3)
            if s_.p < 0.05:
                ax.text(i, T, s_[f"letter_{s}"], ha="center", va="center", fontsize=23, family=FONT)
        if d.pct.max() == 0:
            ax.text(1.5, T * 0.5, "no MAG", ha="center", va="center", fontsize=20, color="#7f7f7f",
                    family=FONT, style="italic")
        ax.set_xlim(-40 / 70, 3 + (x1 - x0 - 250) / 70)
        ax.set_ylim(-0.047 * T, 1.047 * T)
        ax.set_xticks(range(4))
        ax.set_xticklabels(SITES if last else [], fontsize=22.5, color="#4d4d4d", family=FONT)
        from matplotlib.ticker import MaxNLocator
        yt = MaxNLocator(nbins=3, steps=[1, 2, 2.5, 3, 5, 10]).tick_values(0, T)
        yt = [t for t in yt if -1e-9 <= t <= T + 1e-9]
        ax.set_yticks(yt)
        ax.set_yticks([(a + b) / 2 for a, b in zip(yt[:-1], yt[1:])], minor=True)
        lab = [("%g" % t) for t in yt]
        ax.set_yticklabels(lab if k == 0 else [], fontsize=22.5, color="#4d4d4d", family=FONT)
        ax.tick_params(axis="both", which="major", length=6, width=1.1, color="#333333", pad=4)
        ax.tick_params(axis="both", which="minor", length=0)
        ax.tick_params(axis="x", length=6 if last else 0)
        if k > 0:
            ax.tick_params(axis="y", which="both", length=0)
        ax.grid(True, which="major", color="#ebebeb", lw=2.0)
        ax.grid(True, which="minor", axis="y", color="#ebebeb", lw=1.0)
        ax.set_axisbelow(True)
        for sp in ax.spines.values():
            sp.set_color("#333333"); sp.set_linewidth(1.4)
    sx0, sx1 = STRIP_X
    sax = fig.add_axes([sx0 / W, 1 - (y0 + ROW_H) / H, (sx1 - sx0) / W, ROW_H / H])
    sax.set_xticks([]); sax.set_yticks([]); sax.set_facecolor("#d9d9d9")
    for sp in sax.spines.values():
        sp.set_color("#333333"); sp.set_linewidth(1.4)
    sax.text(0.5, 0.5, strip_label, rotation=270, ha="center", va="center", fontsize=25, family=FONT,
             color="#1a1a1a", transform=sax.transAxes)
    return T


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    parts = [mag_row()]
    rp = os.path.join(R26, "read_markers_v2", "results", "marker_genus_profile.tsv")
    have_reads = os.path.exists(rp)
    if have_reads:
        parts.insert(0, read_row(rp))
        sens = read_sensitivity(rp)
        sens.to_csv(os.path.join(OUT_DIR, "figureS6_read_sensitivity.tsv"), sep="\t", index=False)
        print(sens.round(4).to_string())
    df = pd.concat(parts, ignore_index=True)
    st = stats_for(df)
    df.to_csv(os.path.join(OUT_DIR, "figureS6_values.tsv"), sep="\t", index=False)
    st.to_csv(os.path.join(OUT_DIR, "figureS6_stats.tsv"), sep="\t", index=False)
    print(st[["row", "genus", "p"] + [f"mean_{s}" for s in SITES] + [f"letter_{s}" for s in SITES]].round(4).to_string())

    rows = ([("Metagenome (reads)", "Metagenome")] if have_reads else []) + [("MAGs", "MAGs")]
    W = 2048
    first_y = ROW0[1] + GAP
    last_bottom = first_y + len(rows) * ROW_H + (len(rows) - 1) * GAP
    H = last_bottom + 112
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100)
    fig.patch.set_facecolor("white")
    tops = {}
    for r, (name, lab) in enumerate(rows):
        y0 = first_y + r * (ROW_H + GAP)
        tops[name] = draw_row(fig, W, H, y0, df[df.row == name], st[st.row == name], lab, r == len(rows) - 1)
    fig.text(((PX[0][0] + PX[-1][1]) / 2) / W, 1 - (last_bottom + 82) / H, "Site", ha="center", va="center",
             fontsize=27, weight="bold", family=FONT)
    fig.text(30 / W, 1 - ((ROW0[0] + last_bottom) / 2) / H, "Relative Abundance (%)", rotation=90, ha="center",
             va="center", fontsize=27, weight="bold", family=FONT)
    tmp = os.path.join(OUT_DIR, "_figS6_rows.png")
    fig.savefig(tmp, dpi=100, facecolor="white")
    plt.close(fig)

    orig = Image.open(ORIG).convert("RGB")
    canvas = Image.open(tmp).convert("RGB")
    canvas.paste(orig.crop((58, 0, 1760, ROW0[1] + 2)), (58, 0))               # genus strips + 16S row + tick labels
    leg = orig.crop((1790, 420, 2040, 685))
    canvas.paste(leg, (1790, int((ROW0[0] + last_bottom) / 2 - (685 - 420) / 2)))
    canvas.save(OUT_FIG, dpi=(200, 200))
    os.remove(tmp)
    json.dump(dict(rows=[r[0] for r in rows], have_reads=have_reads, y_top=tops),
              open(os.path.join(OUT_DIR, "figureS6_meta.json"), "w"), indent=1)
    print("wrote", OUT_FIG, canvas.size)


if __name__ == "__main__":
    main()
