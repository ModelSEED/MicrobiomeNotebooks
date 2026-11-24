#!/usr/bin/env python
"""Summarise genes_vs_rast.tsv (from run_kofam_controls.py) into control_summary.tsv and min_contig_effect.tsv.

control_summary.tsv: per control MAG, for the guild markers, the number of genes RAST names (regex on the RAST
function, as guilds_r232.csv did) next to the number KOfam calls under DRAM's threshold rule, plus the best
score / threshold of the RAST-named genes that KOfam did not call.
min_contig_effect.tsv: KOfam marker calls in the re-assembly control bins with contigs >= 2,500 bp (DRAM default)
versus >= 1,000 bp (the assembly's own minimum).
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
t = pd.read_csv(HERE / "genes_vs_rast.tsv", sep="\t")
t["rast_function"] = t.rast_function.fillna("")
t["kofam_gene"] = t.kofam_gene.fillna("")
RAST = {   # marker -> (RAST regex as in correlational/r232/rerun_correlational_r232.py where one exists, KOfam genes)
    "phnJ (C-P lyase)": (r"5-phosphate C-P lyase|EC 4\.7\.1\.1\b", ["phnJ"]),
    "grdH (betaine reductase B beta, Sec)": (r"Betaine reductase component B beta", ["grdH"]),
    "grdI (betaine reductase B alpha)": (r"Betaine reductase component B alpha", ["grdI"]),
    "mcrA": (r"Methyl coenzyme M reductase alpha", ["mcrA"]),
    "mtbB/mtmB (Pyl)": (r"\[(Dimethylamine|Methylamine)--corrinoid protein\] Co-methyltransferase", ["mtbB", "mtmB"]),
    "mttB, Pyl-type per RAST": (r"\[Trimethylamine--corrinoid protein\] Co-methyltransferase.*pyrrolysine", ["mttB"]),
    "MttB family (RAST 'TMA methyltransferase family')": (r"Trimethylamine methyltransferase family|Mll1230|"
                                                          r"Trimethylamine:corrinoid methyltransferase-like", []),
    "cntA/yeaW": (r"Carnitine monooxygenase, oxygenase", ["cntA/yeaW"]),
}
orig = t[(t.min_contig == 2500) & ~t.mag.str.startswith("new_")]
rows = []
for mag, g in orig.groupby("mag", sort=False):
    for marker, (pat, kos) in RAST.items():
        rast = g[g.rast_function.str.contains(pat, case=False, regex=True)]
        kof = g[g.kofam_gene.isin(kos)] if kos else g.iloc[0:0]
        missed = rast[~rast.kofam_gene.isin(kos)] if kos else rast
        if rast.empty and kof.empty:
            continue
        rows.append({
            "mag": mag, "marker": marker, "rast_genes": len(rast), "kofam_calls": len(kof),
            "rast_genes_kofam_called_as": "; ".join(sorted(set(rast.kofam_gene.replace("", "none")))),
            "missed_best_score/threshold": "; ".join(
                f"{s:.0f}/{float(th):.0f} ({b})" for s, th, b in zip(missed.best_score, missed.best_threshold,
                                                                      missed.best_marker_gene_any_score)
                if pd.notna(s) and str(th) not in ("nan", "-")),
            "kofam_called_genes_rast_says": "; ".join(sorted(set(kof.rast_function.replace("", "(no RAST CDS)"))))[:300],
        })
s = pd.DataFrame(rows)
s.to_csv(HERE / "control_summary.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 60)
print(s.drop(columns="kofam_called_genes_rast_says").to_string(index=False))

# guild retention: every RAST guild member of guilds_r232.csv, re-called with DRAM's KOfam rule
guilds = pd.read_csv(HERE.parents[2] / "correlational/r232/guilds_r232.csv").set_index("iterativeID_r232")
calls = orig.groupby("mag").kofam_gene.apply(set)
ret = []
for lab, g in guilds.iterrows():
    if lab not in calls.index:
        continue
    k = calls[lab]
    ret.append({"mag": lab, "genus": g.genus, "domain": g.domain,
                "rast_mepn_consumer(phnJ)": g.mepn_consumer, "kofam_phnJ": "phnJ" in k,
                "kofam_phnGHIJKLM": "".join(x[-1] for x in ["phnG", "phnH", "phnI", "phnJ", "phnK", "phnL", "phnM"]
                                            if x in k),
                "rast_betaine_reducer(grdH)": g.betaine_reducer, "kofam_grdH": "grdH" in k, "kofam_grdI": "grdI" in k,
                "kofam_grdB_or_E": bool({"grdB", "grdE"} & k),
                "rast_methanogen(mcrA)": g.methanogen, "kofam_mcrA": "mcrA" in k,
                "rast_tma_producer": g.tma_producer, "kofam_cutC": "cutC" in k, "kofam_cntA": "cntA/yeaW" in k,
                "kofam_cntB": "cntB/yeaX" in k, "kofam_mttB_K14083": "mttB" in k, "kofam_mtgB_K25220": "mtgB" in k})
ret = pd.DataFrame(ret)
ret.to_csv(HERE / "guild_retention.tsv", sep="\t", index=False)
mp, br = ret[ret["rast_mepn_consumer(phnJ)"]], ret[ret["rast_betaine_reducer(grdH)"]]
print(f"\nMePn consumers (RAST phnJ): {len(mp)}; KOfam phnJ: {int(mp.kofam_phnJ.sum())}; "
      f"missed: {', '.join(mp.loc[~mp.kofam_phnJ, 'mag'])}")
print(f"KOfam phnJ in MAGs RAST did not call: {', '.join(ret.loc[ret.kofam_phnJ & ~ret['rast_mepn_consumer(phnJ)'], 'mag'])}")
print(f"betaine reducers (RAST grdH|grdI string): {len(br)}; KOfam grdH: {int(br.kofam_grdH.sum())}; grdI: "
      f"{int(br.kofam_grdI.sum())}; grdH|grdI: {int((br.kofam_grdH | br.kofam_grdI).sum())}; "
      f"missed by grdH|grdI: {', '.join(br.loc[~(br.kofam_grdH | br.kofam_grdI), 'mag'])}")
print(f"KOfam mttB (K14083) in non-methanogens: {int((ret.kofam_mttB_K14083 & ~ret.kofam_mcrA).sum())} of "
      f"{int((~ret.kofam_mcrA).sum())} controls")
print(ret.drop(columns=["domain"]).to_string(index=False))

new = t[t.mag.str.startswith("new_") & (t.kofam_gene != "")]
m = new.pivot_table(index="kofam_gene", columns=["mag", "min_contig"], values="gene", aggfunc="count", fill_value=0)
m.to_csv(HERE / "min_contig_effect.tsv", sep="\t")
print()
print(m.to_string())
