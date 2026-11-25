#!/usr/bin/env python
"""Rebuild 'Table S5 - MAG taxonomy, quality and guilds (GTDB R232).xlsx' (the workbook was first made without a script;
this rebuild, moved here 2026-11-25 from prfba_augmentation/scouting_2026-11-05/inventory/, matches it cell for cell except
that empty taxonomy ranks are written as empty strings)
from the R232 outputs, so that it can be regenerated from an alternative-annotation run.

Inputs (all read-only): <run>/guilds_r232.csv, <run>/mag_methane_correlations_r232.csv  (default <run> = the
committed correlational/r232), data/Saltern_phylogeny_r232.json (Class), MAG_QC/gtdbtk_local_r232/*.summary.tsv
(classification string), data/Cliff_310MAG_relabund.tsv (site-type means: plain mean of the samples of each site type).
--compare checks every cell against the committed workbook.
Run: ~/Documents/py_venv/bin/python -B make_table_s5.py --run DIR --out FILE.xlsx [--compare]
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

S = Path(__file__).resolve().parents[2]          # salterns/
COMMITTED = S / "Table S5 - MAG taxonomy, quality and guilds (GTDB R232).xlsx"


def site_type(sample):
    return "restored" if sample.startswith("SF2") else "reference" if sample.startswith("R2A") else "unrestored"


def build(run):
    g = pd.read_csv(run / "guilds_r232.csv")
    m = pd.read_csv(run / "mag_methane_correlations_r232.csv").set_index("mag")
    ph = json.load(open(S / "data/Saltern_phylogeny_r232.json"))
    gt = pd.concat([pd.read_csv(S / f"MAG_QC/gtdbtk_local_r232/gtdbtk.{d}.summary.tsv", sep="\t") for d in ("ar53", "bac120")])
    gt["mag"] = gt.user_genome.str.replace(".contigs__", "", regex=False)
    cls = gt.set_index("mag")["classification"]
    rel = pd.read_csv(S / "data/Cliff_310MAG_relabund.tsv", sep="\t", index_col=0)
    site = rel.T.groupby(rel.columns.map(site_type)).mean().T * 100
    rows = []
    for r in g.itertuples():
        t = ph[r.mag + ".contigs__"]
        rows.append({
            "MAG (DASTool bin)": r.mag, "iterativeID (GTDB R232)": r.iterativeID_r232,
            "iterativeID (2020)": r.iterativeID_2020 if isinstance(r.iterativeID_2020, str) else np.nan,
            "GTDB R232 classification": cls[r.mag],
            "Domain": t["Kingdom"], "Phylum": t["Phylum"], "Class": t["Class"], "Order": t["Order"],
            "Family": t["Family"], "Genus": t["Genus"], "Species": t["Species"],
            "CheckM completeness (%)": r.completeness, "CheckM contamination (%)": r.contamination,
            "MIMAG medium quality or better": bool(r.completeness >= 50 and r.contamination < 10),
            "MePn consumer (phnJ)": r.mepn_consumer, "MePn consumer (2020 set)": r.mepn_2020,
            "Betaine reducer (grdH)": r.betaine_reducer, "Betaine reducer (2020 set)": r.betaine_2020,
            "Methanogen (mcrA)": r.methanogen, "TMA-route candidate (cntA / cutC-cluster genes)": r.tma_producer,
            "phnJ annotation": r.ev_phnJ, "grdH annotation": r.ev_grdH, "mcrA annotation": r.ev_mcrA,
            "cntA annotation": r.ev_cntA, "cutC-cluster annotation": r.ev_cutC_cluster,
            "Mean rel. abundance, unrestored (%)": site.loc[r.mag, "unrestored"],
            "Mean rel. abundance, restored (%)": site.loc[r.mag, "restored"],
            "Mean rel. abundance, reference (%)": site.loc[r.mag, "reference"],
            "Spearman rho vs CH4, areal (12 cores)": m.loc[r.mag, "rho_raw"], "BH q (areal)": m.loc[r.mag, "q_raw"],
            "Spearman rho vs CH4, DNA-normalised (12 cores)": m.loc[r.mag, "rho_dna"],
            "BH q (DNA-normalised)": m.loc[r.mag, "q_dna"],
        })
    out = pd.DataFrame(rows)
    for c in out.columns:
        if out[c].dtype == object:
            out[c] = out[c].replace("", np.nan)
    # committed row order (verified): archaea first, then phylum, then iterativeID_r232 (plain string sort)
    out = (out.assign(_d=out.Domain.ne("Archaea"))
              .sort_values(["_d", "Phylum", "iterativeID (GTDB R232)"]).drop(columns="_d").reset_index(drop=True))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", type=Path, default=S / "correlational/r232")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--compare", action="store_true")
    a = ap.parse_args()
    t = build(a.run)
    if a.compare:
        ref = pd.read_excel(COMMITTED, sheet_name="Table S5")
        print("columns identical:", list(ref.columns) == list(t.columns))
        t2 = t.set_index("MAG (DASTool bin)").loc[ref["MAG (DASTool bin)"]].reset_index()
        rep = {}
        for c in ref.columns:
            x, y = ref[c], t2[c]
            if pd.api.types.is_float_dtype(x) and pd.api.types.is_numeric_dtype(y):
                d = np.abs(x.values - y.values.astype(float))
                rep[c] = f"max |diff| {np.nanmax(d) if np.isfinite(d).any() else 0:.2e}; NaN pattern same: {bool((x.isna() == y.isna()).all())}"
            else:
                same = ((x == y) | (x.isna() & y.isna())).all()
                rep[c] = "identical" if same else f"{int((~((x == y) | (x.isna() & y.isna()))).sum())} cells differ"
        for c, v in rep.items():
            print(f"  {c:50s} {v}")
        # row order of the committed file
        print("committed order == guilds_r232.csv order:", list(ref["MAG (DASTool bin)"]) == list(t["MAG (DASTool bin)"]))
        t = t2                                       # write in the committed row order
    a.out.parent.mkdir(parents=True, exist_ok=True)
    t.to_excel(a.out, sheet_name="Table S5", index=False)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
