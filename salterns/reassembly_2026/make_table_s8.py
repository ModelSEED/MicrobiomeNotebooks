#!/usr/bin/env python
"""Table S8 (separate file, like Table S5): the final re-assembled MAG set, one row per MAG — 935 second-pass (Arm C)
MAGs + 520 kept MAGs of the four first-pass site co-assemblies, all scored with CheckM2 1.1.0 + GUNC 1.1.1, classified
with GTDB-Tk 2.8.0 / R232, clustered into species (skani 95 % ANI / 50 % AF), with KOfam guild calls.
Inputs: comparison_rollout/mags_three_way.tsv, comparison_rollout/final_guilds.tsv (final_set_guilds.py),
rollout_final/results/mags_kept.tsv, coasm_final/coasm_bins_all.tsv, and for the species representatives the labels,
abundance and methane correlations of ../correlational/final585/ (rerun_correlational_final585.py).
Writes ../Table S8 - re-assembled MAG set taxonomy, quality and guilds (GTDB R232).xlsx
Run: ~/Documents/py_venv/bin/python make_table_s8.py
"""
import os

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font

HERE = os.path.dirname(os.path.abspath(__file__))
OUTF = os.path.join(os.path.dirname(HERE), "Table S8 - re-assembled MAG set taxonomy, quality and guilds (GTDB R232).xlsx")


def yn(v):
    return "" if pd.isna(v) else ("yes" if bool(v) else "no")


def main():
    d = pd.read_csv(os.path.join(HERE, "comparison_rollout/mags_three_way.tsv"), sep="\t")
    g = pd.read_csv(os.path.join(HERE, "comparison_rollout/final_guilds.tsv"), sep="\t").set_index("genome")
    k = pd.read_csv(os.path.join(HERE, "rollout_final/results/mags_kept.tsv"), sep="\t").set_index("mag_id")
    c = pd.read_csv(os.path.join(HERE, "coasm_final/coasm_bins_all.tsv"), sep="\t").set_index("mag_id")
    f = d[d.set.isin(["rollout", "coassembly"])].copy().rename(columns={"class": "cls"})
    assert len(f) == len(g) == 1455 and set(f.genome) == set(g.index)
    f["name"] = f.genome.str.split("|", n=1).str[1]
    reps = f[f.final_species_rep].set_index("species_cluster").name
    RF = os.path.join(os.path.dirname(HERE), "correlational", "final585")
    gl = pd.read_csv(os.path.join(RF, "guilds_final585.csv")).set_index("species_cluster")
    ab12 = pd.read_csv(os.path.join(RF, "averaged_normalized_MAG_abundances_final585.csv"), index_col=0)
    mc = pd.read_csv(os.path.join(RF, "mag_methane_correlations_final585.csv")).set_index("mag")
    types = {"unrestored": ("R1_", "R2_"), "restored": ("SF2",), "reference": ("R2A",)}
    mean_t = {t: ab12[[c for c in ab12.columns if c.startswith(pre)]].mean(axis=1) * 100 for t, pre in types.items()}
    assert set(gl.mag) == set(ab12.index) == set(reps.values)
    rows = []
    for r in f.sort_values(["set", "sample", "name"], ascending=[False, True, True]).itertuples():
        q = (k if r.set == "rollout" else c).loc[r.name]
        gg = g.loc[r.genome]
        rep_mag = r.final_species_rep
        sp = gl.loc[r.species_cluster]
        ab = lambda t: round(float(mean_t[t][r.name]), 5) if rep_mag else None
        rows.append({
            "MAG": r.name, "Species label": sp.iterativeID_final585,
            "Source": "second-pass per-sample assembly" if r.set == "rollout" else "first-pass site co-assembly (MEGAHIT)",
            "Sample / co-assembly": r.sample, "Site": r.site,
            "GTDB R232 classification": r.classification, "Domain": r.domain, "Phylum": r.phylum, "Class": r.cls,
            "Order": r.order, "Family": r.family, "Genus": r.genus, "Species": r.species,
            "CheckM2 completeness (%)": q.checkm2_completeness, "CheckM2 contamination (%)": q.checkm2_contamination,
            "GUNC pass": q.gunc_pass, "CheckM (v1.1.3) completeness (%)": q.checkm1_completeness,
            "CheckM (v1.1.3) contamination (%)": q.checkm1_contamination,
            "Quality (MIMAG, CheckM2)": "high" if q.quality == "HQ" else "medium",
            "Size (bp)": int(q.size_bp), "Contigs": int(q.n_contigs), "N50 (bp)": int(q.n50),
            "Species cluster (representative MAG)": reps.get(r.species_cluster, ""),
            "Species representative": yn(r.final_species_rep),
            "In a species with an original (2020) MAG": yn(any(s == "original" for s in str(r.cluster_sets).split("+"))),
            "KOfam calls": {"kofam_rollout": "this study (second pass)", "kofam (first pass)": "this study (first pass)",
                            "kofam_coasm38": "this study (co-assembly bins)", "missing": "pending"}[gg.kofam_source],
            "MePn consumer (phnJ)": yn(gg.mepn_consumer), "Betaine reducer (grdH or grdI)": yn(gg.betaine_reducer),
            "Betaine reducer evidence": "" if pd.isna(gg.betaine_reducer_evidence) else gg.betaine_reducer_evidence,
            "Methanogen (mcrA)": yn(gg.methanogen), "Methylamine user (mttB/mtbB/mtmB with mcrA)": yn(gg.methylamine_user),
            "TMA producer (cutC+cutD or cntA+cntB)": yn(gg.tma_producer),
            "TMA candidate, one gene only (cutC or cntA)": yn(gg.tma_candidate_one_sided),
            "Species: mean rel. abundance, unrestored R1/R2 (% of reads)": ab("unrestored"),
            "Species: mean rel. abundance, restored SF2 (%)": ab("restored"),
            "Species: mean rel. abundance, reference R2A (%)": ab("reference"),
            "Species: Spearman rho with CH4 flux (12 cores, areal)": round(float(mc.loc[r.name, "rho_raw"]), 3) if rep_mag and pd.notna(mc.loc[r.name, "rho_raw"]) else None,
            "Species: BH q": round(float(mc.loc[r.name, "q_raw"]), 4) if rep_mag and pd.notna(mc.loc[r.name, "q_raw"]) else None,
        })
    t = pd.DataFrame(rows)
    t.to_excel(OUTF, sheet_name="Table S8", index=False)
    wb = load_workbook(OUTF); ws = wb["Table S8"]
    for cell in ws[1]:
        cell.font = Font(bold=True); cell.alignment = Alignment(wrap_text=True, vertical="top")
    ws.freeze_panes = "C2"
    for col in ws.columns:
        n = max(len(str(x.value)) if x.value is not None else 0 for x in col[:200])
        ws.column_dimensions[col[0].column_letter].width = min(max(10, n + 2), 60)
    ws.row_dimensions[1].height = 45
    wb.save(OUTF)
    print(f"wrote {OUTF}: {len(t)} MAGs, {t['Species representative'].eq('yes').sum()} species representatives, "
          f"KOfam pending for {t['KOfam calls'].eq('pending').sum()}")


if __name__ == "__main__":
    main()
